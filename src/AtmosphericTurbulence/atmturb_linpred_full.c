// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_linpred_full.c
 * @brief   Full-aperture linear predictor construction and evaluation for predictive control
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef HAVE_OPENBLAS_LAPACKE
#include <cblas.h>
#endif

#include "AtmosphericTurbulence.h"
#include "atmturb_linpred_solve.h"
#include "atmturb_types.h"

/**
 * struct atmturb_linpred_dims_t - Dimension bundle for full-aperture linear predictor build
 * @nxy:      Pixels per 2D wavefront frame (nx * ny).
 * @nbpix:    Number of active pupil pixels.
 * @mvecsize: Measurement history vector length (nbpix * pf_order).
 * @nbmvec:   Number of training samples.
 * @pf_order: Autoregressive history order.
 * @pf_lag:   Integer prediction horizon lag.
 */
typedef struct
{
    long nxy;
    long nbpix;
    long mvecsize;
    long nbmvec;
    int  pf_order;
    long pf_lag;
} atmturb_linpred_dims_t;

/**
 * atmturb_compute_wavefront_piston - Compute per-frame mean piston non-destructively
 * @id_wf:    Input wavefront cube image ID.
 * @mask_idx: Array of active mask pixel indices.
 * @nbpix:    Number of active mask pixels.
 * @nz:       Number of time slices.
 * @nxy:      Frame pixel count (nx * ny).
 *
 * Return: Allocated piston array (length nz), or NULL on failure.
 */
static float *atmturb_compute_wavefront_piston(
    imageID     id_wf,
    const long *mask_idx,
    long        nbpix,
    long        nz,
    long        nxy)
{
    if (nbpix <= 0 || nz <= 0)
    {
        return NULL;
    }

    float *piston = (float *) malloc(sizeof(float) * (size_t) nz);
    if (piston == NULL)
    {
        return NULL;
    }

    const float *wf_data = dcimg[id_wf].array.F;
    double       inv_np  = 1.0 / (double) nbpix;

    #pragma omp parallel for schedule(static)
    for (long kk = 0; kk < nz; kk++)
    {
        double       tot   = 0.0;
        const float *slice = &wf_data[kk * nxy];
        for (long p = 0; p < nbpix; p++)
        {
            tot += (double) slice[mask_idx[p]];
        }
        piston[kk] = (float) (tot * inv_np);
    }
    return piston;
}

/**
 * atmturb_extract_mask_indices - Collect active pupil pixel linear indices
 * @id_mask:   Pupil mask image ID.
 * @nx:        Mask width.
 * @ny:        Mask height.
 * @nbpix_out: Output count of active pixels.
 *
 * Return: Allocated array of pixel indices, or NULL on allocation failure.
 */
static long *atmturb_extract_mask_indices(
    imageID  id_mask,
    long     nx,
    long     ny,
    long    *nbpix_out)
{
    long         count     = 0;
    long         nxy       = nx * ny;
    const float *mask_data = dcimg[id_mask].array.F;

    for (long i = 0; i < nxy; i++)
    {
        if (mask_data[i] > 0.5f)
        {
            count++;
        }
    }

    long *indices = (long *) malloc(sizeof(long) * (size_t) (count > 0 ? count : 1));
    if (indices == NULL)
    {
        *nbpix_out = 0;
        return NULL;
    }

    long idx = 0;
    for (long i = 0; i < nxy; i++)
    {
        if (mask_data[i] > 0.5f)
        {
            indices[idx++] = i;
        }
    }

    *nbpix_out = count;
    return indices;
}

/**
 * atmturb_linpred_populate_matrices - Populate piston-removed training matrices matA and matC
 * @wf_data:  Input 3D wavefront cube data.
 * @mat_a:    Destination history matrix [nbmvec x mvecsize].
 * @mat_c:    Destination target matrix [nbmvec x nbpix].
 * @mask_idx: Active pupil pixel index map [nbpix].
 * @piston:   Per-frame piston array [nz].
 * @d:        Predictor dimension bundle.
 */
static void atmturb_linpred_populate_matrices(
    const float                  *restrict wf_data,
    float                        *restrict mat_a,
    float                        *restrict mat_c,
    const long                   *restrict mask_idx,
    const float                  *restrict piston,
    const atmturb_linpred_dims_t          *d)
{
    #pragma omp parallel for schedule(static)
    for (long m = 0; m < d->nbmvec; m++)
    {
        float *a_row = &mat_a[m * d->mvecsize];
        for (int dt = 0; dt < d->pf_order; dt++)
        {
            long         t     = m + (long) dt;
            float        p_t   = piston[t];
            const float *slice = &wf_data[t * d->nxy];
            float       *a_dst = &a_row[(long) dt * d->nbpix];
            for (long p = 0; p < d->nbpix; p++)
            {
                a_dst[p] = slice[mask_idx[p]] - p_t;
            }
        }

        long         t_targ = m + (long) d->pf_order - 1 + d->pf_lag;
        float        p_targ = piston[t_targ];
        const float *tslice = &wf_data[t_targ * d->nxy];
        float       *c_row  = &mat_c[m * d->nbpix];
        for (long p = 0; p < d->nbpix; p++)
        {
            c_row[p] = tslice[mask_idx[p]] - p_targ;
        }
    }
}

/**
 * AtmosphericTurbulence_Build_LinPredictor_Full - Build full-aperture AR linear prediction matrix
 * @WFin_name:   Input wavefront cube name.
 * @WFmask_name: Pupil mask image name.
 * @PForder:     Autoregressive predictor order (history steps).
 * @PFlag:       Prediction horizon lag.
 * @SVDeps:      Singular value cutoff tolerance.
 * @RegLambda:   Tikhonov regularization parameter.
 *
 * Return: 0 on success, -1 on failure.
 */
int AtmosphericTurbulence_Build_LinPredictor_Full(
    const char *WFin_name,
    const char *WFmask_name,
    int         PForder,
    float       PFlag,
    double      SVDeps,
    double      RegLambda)
{
    imageID id_wfin   = image_ID(WFin_name);
    imageID id_wfmask = image_ID(WFmask_name);
    if (id_wfin < 0 || id_wfmask < 0)
    {
        return -1;
    }

    long  nx       = dcimg[id_wfin].md[0].size[0];
    long  ny       = dcimg[id_wfin].md[0].size[1];
    long  nz       = dcimg[id_wfin].md[0].size[2];
    long  nbpix    = 0;
    long *mask_idx = atmturb_extract_mask_indices(id_wfmask, nx, ny, &nbpix);
    if (mask_idx == NULL || nbpix <= 0)
    {
        free(mask_idx);
        return -1;
    }

    atmturb_linpred_dims_t d = {
        .nxy      = nx * ny,
        .nbpix    = nbpix,
        .mvecsize = nbpix * (long) PForder,
        .nbmvec   = nz - (long) PForder - (long) PFlag,
        .pf_order = PForder,
        .pf_lag   = (long) PFlag
    };
    if (d.nbmvec <= 0 || d.mvecsize <= 0)
    {
        free(mask_idx);
        return -1;
    }

    float *piston = atmturb_compute_wavefront_piston(id_wfin, mask_idx, nbpix, nz, d.nxy);
    float *mat_a  = (float *) malloc(sizeof(float) * (size_t) d.nbmvec * (size_t) d.mvecsize);
    float *mat_c  = (float *) malloc(sizeof(float) * (size_t) d.nbmvec * (size_t) nbpix);
    if (piston == NULL || mat_a == NULL || mat_c == NULL)
    {
        free(piston);
        free(mat_a);
        free(mat_c);
        free(mask_idx);
        return -1;
    }

    atmturb_linpred_populate_matrices(
        dcimg[id_wfin].array.F, mat_a, mat_c, mask_idx, piston, &d);
    free(piston);
    free(mask_idx);

    char filtname[200];
    snprintf(filtname, sizeof(filtname), "Pfilt_%04d_%04ld", PForder, (long) (PFlag + 0.5f));
    imageID id_old = image_ID(filtname);
    if (id_old >= 0)
    {
        delete_image_ID(filtname);
    }
    imageID id_filt = create_2Dimage_ID(filtname, d.mvecsize, nbpix);
    if (id_filt < 0)
    {
        free(mat_a);
        free(mat_c);
        return -1;
    }

    double                    eps  = (RegLambda > 0.0 && RegLambda > SVDeps) ? RegLambda : SVDeps;
    atmturb_linpred_problem_t prob = {
        .mat_a    = mat_a,
        .mat_c    = mat_c,
        .filt_out = dcimg[id_filt].array.F,
        .nbmvec   = d.nbmvec,
        .mvecsize = d.mvecsize,
        .nbpix    = nbpix,
        .reg_eps  = (eps > 0.0) ? eps : 1e-3
    };
    int ret = atmturb_linpred_solve_normal(&prob);

    free(mat_a);
    free(mat_c);
    return ret;
}

/**
 * atmturb_get_or_create_3d_float - Reuse matching 3D float image or create a fresh one
 */
static imageID atmturb_get_or_create_3d_float(
    const char *name,
    long        nx,
    long        ny,
    long        nz)
{
    imageID id = image_ID(name);
    if (id >= 0 && dcimg[id].used == 1 &&
        dcimg[id].md[0].naxis == 3 &&
        dcimg[id].md[0].datatype == _DATATYPE_FLOAT &&
        (long) dcimg[id].md[0].size[0] == nx &&
        (long) dcimg[id].md[0].size[1] == ny &&
        (long) dcimg[id].md[0].size[2] == nz)
    {
        memset(dcimg[id].array.F, 0, sizeof(float) * (size_t) nx * (size_t) ny * (size_t) nz);
        return id;
    }
    if (id >= 0)
    {
        delete_image_ID(name);
    }
    return create_3Dimage_ID(name, nx, ny, nz);
}

typedef struct
{
    const float *wf_in;
    const float *filt;
    const long  *mask_idx;
    float       *wf_outp;
    float       *wf_outf;
    float       *hist_pack;
    long         nxy;
    long         nbpix;
    long         mvecsize;
    int          pf_order;
    long         pf_lag;
} atmturb_linpred_apply_ctx_t;

/**
 * atmturb_linpred_eval_tile16 - Evaluate 4x16 register-blocked SGEMM for up to 16 frames
 */
static void atmturb_linpred_eval_tile16(
    const atmturb_linpred_apply_ctx_t *c,
    long                               t0,
    int                                blk)
{
    memset(c->hist_pack, 0, sizeof(float) * (size_t) c->mvecsize * 16U);
    for (int dt = 0; dt < c->pf_order; dt++)
    {
        for (int b = 0; b < blk; b++)
        {
            long         t      = t0 + (long) b - (long) c->pf_order + (long) dt;
            const float *hslice = &c->wf_in[t * c->nxy];
            long         k_base = (long) dt * c->nbpix;
            for (long q = 0; q < c->nbpix; q++)
            {
                c->hist_pack[(k_base + q) * 16 + b] = hslice[c->mask_idx[q]];
            }
        }
    }

    #pragma omp parallel for schedule(static)
    for (long p0 = 0; p0 < c->nbpix; p0 += 4)
    {
        long p_end = (p0 + 4 <= c->nbpix) ? (p0 + 4) : c->nbpix;
        if (p_end - p0 == 4)
        {
            float acc0[16] = {0.0f};
            float acc1[16] = {0.0f};
            float acc2[16] = {0.0f};
            float acc3[16] = {0.0f};
            const float *restrict r0 = &c->filt[(p0 + 0) * c->mvecsize];
            const float *restrict r1 = &c->filt[(p0 + 1) * c->mvecsize];
            const float *restrict r2 = &c->filt[(p0 + 2) * c->mvecsize];
            const float *restrict r3 = &c->filt[(p0 + 3) * c->mvecsize];

            for (long k = 0; k < c->mvecsize; k++)
            {
                const float *restrict hp = &c->hist_pack[k * 16];
                float                 f0 = r0[k];
                float                 f1 = r1[k];
                float                 f2 = r2[k];
                float                 f3 = r3[k];
                #pragma omp simd
                for (int b = 0; b < 16; b++)
                {
                    float hv = hp[b];
                    acc0[b] += f0 * hv;
                    acc1[b] += f1 * hv;
                    acc2[b] += f2 * hv;
                    acc3[b] += f3 * hv;
                }
            }

            const float *accs[4] = {acc0, acc1, acc2, acc3};
            for (int r = 0; r < 4; r++)
            {
                long mp = c->mask_idx[p0 + r];
                for (int b = 0; b < blk; b++)
                {
                    long  toff = (t0 + (long) b - 1 + c->pf_lag) * c->nxy + mp;
                    float pf   = accs[r][b];
                    c->wf_outp[toff] = pf;
                    c->wf_outf[toff] = c->wf_in[toff] - pf;
                }
            }
        }
        else
        {
            for (long p = p0; p < p_end; p++)
            {
                float                 acc[16] = {0.0f};
                const float *restrict r0      = &c->filt[p * c->mvecsize];
                for (long k = 0; k < c->mvecsize; k++)
                {
                    const float *restrict hp = &c->hist_pack[k * 16];
                    float                 f0 = r0[k];
                    #pragma omp simd
                    for (int b = 0; b < 16; b++)
                    {
                        acc[b] += f0 * hp[b];
                    }
                }
                long mp = c->mask_idx[p];
                for (int b = 0; b < blk; b++)
                {
                    long  toff = (t0 + (long) b - 1 + c->pf_lag) * c->nxy + mp;
                    float pf   = acc[b];
                    c->wf_outp[toff] = pf;
                    c->wf_outf[toff] = c->wf_in[toff] - pf;
                }
            }
        }
    }
}

#ifdef HAVE_OPENBLAS_LAPACKE
/**
 * atmturb_linpred_eval_blas - Evaluate full-aperture predictor across all frames via SGEMM
 */
static int atmturb_linpred_eval_blas(
    const atmturb_linpred_apply_ctx_t *c,
    long                               t_max)
{
    long n_eval = t_max - (long) c->pf_order + 1;
    if (n_eval <= 0)
    {
        return 0;
    }

    float *hist = (float *) malloc(sizeof(float) * (size_t) n_eval * (size_t) c->mvecsize);
    float *pred = (float *) malloc(sizeof(float) * (size_t) n_eval * (size_t) c->nbpix);
    if (hist == NULL || pred == NULL)
    {
        free(hist);
        free(pred);
        return -1;
    }

    #pragma omp parallel for schedule(static)
    for (long e = 0; e < n_eval; e++)
    {
        float *h_row = &hist[e * c->mvecsize];
        for (int dt = 0; dt < c->pf_order; dt++)
        {
            const float *hslice = &c->wf_in[(e + (long) dt) * c->nxy];
            float       *h_dst  = &h_row[(long) dt * c->nbpix];
            for (long q = 0; q < c->nbpix; q++)
            {
                h_dst[q] = hslice[c->mask_idx[q]];
            }
        }
    }

    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                (int) n_eval, (int) c->nbpix, (int) c->mvecsize,
                1.0f, hist, (int) c->mvecsize, c->filt, (int) c->mvecsize,
                0.0f, pred, (int) c->nbpix);

    #pragma omp parallel for schedule(static)
    for (long e = 0; e < n_eval; e++)
    {
        long         t_base = ((long) c->pf_order + e - 1 + c->pf_lag) * c->nxy;
        const float *p_row  = &pred[e * c->nbpix];
        for (long p = 0; p < c->nbpix; p++)
        {
            long  toff       = t_base + c->mask_idx[p];
            float pf         = p_row[p];
            c->wf_outp[toff] = pf;
            c->wf_outf[toff] = c->wf_in[toff] - pf;
        }
    }

    free(hist);
    free(pred);
    return 0;
}
#endif

/**
 * AtmosphericTurbulence_Apply_LinPredictor_Full - Evaluate linear prediction on series
 * @MODE:        0 for direct evaluation.
 * @WFin_name:   Input wavefront cube.
 * @WFmask_name: Active pupil mask.
 * @PForder:     Filter AR order.
 * @PFlag:       Prediction horizon lag.
 * @WFoutp_name: Output predicted wavefront cube name.
 * @WFoutf_name: Output residual error cube name.
 *
 * Return: 0 on success, -1 on failure.
 */
/* size-exempt: public API signature in AtmosphericTurbulence.h requires 7 arguments */
int AtmosphericTurbulence_Apply_LinPredictor_Full(
    int         MODE,
    const char *WFin_name,
    const char *WFmask_name,
    int         PForder,
    float       PFlag,
    const char *WFoutp_name,
    const char *WFoutf_name)
{
    (void) MODE;
    imageID id_wfin   = image_ID(WFin_name);
    imageID id_wfmask = image_ID(WFmask_name);
    char    filtname[200];
    snprintf(filtname, sizeof(filtname), "Pfilt_%04d_%04ld", PForder, (long) (PFlag + 0.5f));
    imageID id_filt = image_ID(filtname);
    if (id_wfin < 0 || id_wfmask < 0 || id_filt < 0)
    {
        return -1;
    }

    long  nx       = dcimg[id_wfin].md[0].size[0];
    long  ny       = dcimg[id_wfin].md[0].size[1];
    long  nz       = dcimg[id_wfin].md[0].size[2];
    long  nbpix    = 0;
    long *mask_idx = atmturb_extract_mask_indices(id_wfmask, nx, ny, &nbpix);
    if (mask_idx == NULL || nbpix <= 0)
    {
        free(mask_idx);
        return -1;
    }

    long    mvecsize  = nbpix * (long) PForder;
    long    t_max     = nz - (long) PFlag;
    imageID id_outp   = atmturb_get_or_create_3d_float(WFoutp_name, nx, ny, nz);
    imageID id_outf   = atmturb_get_or_create_3d_float(WFoutf_name, nx, ny, nz);
    float  *hist_pack = NULL;
    if (id_outp < 0 || id_outf < 0 ||
        posix_memalign((void **) &hist_pack, 64, sizeof(float) * (size_t) mvecsize * 16U) != 0)
    {
        free(mask_idx);
        return -1;
    }

    atmturb_linpred_apply_ctx_t ctx = {
        .wf_in     = dcimg[id_wfin].array.F,
        .filt      = dcimg[id_filt].array.F,
        .mask_idx  = mask_idx,
        .wf_outp   = dcimg[id_outp].array.F,
        .wf_outf   = dcimg[id_outf].array.F,
        .hist_pack = hist_pack,
        .nxy       = nx * ny,
        .nbpix     = nbpix,
        .mvecsize  = mvecsize,
        .pf_order  = PForder,
        .pf_lag    = (long) PFlag
    };

#ifdef HAVE_OPENBLAS_LAPACKE
    if (atmturb_linpred_eval_blas(&ctx, t_max) == 0)
    {
        free(hist_pack);
        free(mask_idx);
        return 0;
    }
#endif

    for (long t0 = (long) PForder; t0 <= t_max; t0 += 16)
    {
        int blk = (t0 + 15 <= t_max) ? 16 : (int) (t_max - t0 + 1);
        atmturb_linpred_eval_tile16(&ctx, t0, blk);
    }

    free(hist_pack);
    free(mask_idx);
    return 0;
}
