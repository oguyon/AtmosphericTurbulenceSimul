// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_linpred_full.c
 * @brief   Full-aperture linear predictor construction and evaluation for predictive control
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "AtmosphericTurbulence.h"
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
 * atmturb_linpred_multiply_ainv_c - Compute filter matrix Pfilt = matC^T * matA_inv^T
 * @a_inv:    Pseudo-inverse matrix [mvecsize x nbmvec].
 * @mat_c:    Target matrix [nbmvec x nbpix].
 * @filt_out: Destination predictor matrix [nbpix x mvecsize].
 * @nbpix:    Number of active pupil pixels.
 * @mvecsize: Measurement history vector length.
 * @nbmvec:   Number of training samples.
 */
static void atmturb_linpred_multiply_ainv_c(
    const float *restrict a_inv,
    const float *restrict mat_c,
    float       *restrict filt_out,
    long                  nbpix,
    long                  mvecsize,
    long                  nbmvec)
{
    float *mat_c_t = (float *) malloc(sizeof(float) * (size_t) (nbpix * nbmvec));
    if (mat_c_t != NULL)
    {
        #pragma omp parallel for schedule(static)
        for (long p = 0; p < nbpix; p++)
        {
            float *c_dst = &mat_c_t[p * nbmvec];
            for (long m = 0; m < nbmvec; m++)
            {
                c_dst[m] = mat_c[m * nbpix + p];
            }
        }

        #pragma omp parallel for schedule(static)
        for (long p0 = 0; p0 < nbpix; p0 += 4)
        {
            if (p0 + 4 <= nbpix)
            {
                const float *restrict c0 = &mat_c_t[(p0 + 0) * nbmvec];
                const float *restrict c1 = &mat_c_t[(p0 + 1) * nbmvec];
                const float *restrict c2 = &mat_c_t[(p0 + 2) * nbmvec];
                const float *restrict c3 = &mat_c_t[(p0 + 3) * nbmvec];
                float       *restrict f0 = &filt_out[(p0 + 0) * mvecsize];
                float       *restrict f1 = &filt_out[(p0 + 1) * mvecsize];
                float       *restrict f2 = &filt_out[(p0 + 2) * mvecsize];
                float       *restrict f3 = &filt_out[(p0 + 3) * mvecsize];

                for (long j = 0; j < mvecsize; j++)
                {
                    const float *restrict a_row = &a_inv[j * nbmvec];
                    float                 s0    = 0.0f;
                    float                 s1    = 0.0f;
                    float                 s2    = 0.0f;
                    float                 s3    = 0.0f;
                    #pragma omp simd reduction(+ : s0, s1, s2, s3)
                    for (long m = 0; m < nbmvec; m++)
                    {
                        float av = a_row[m];
                        s0 += av * c0[m];
                        s1 += av * c1[m];
                        s2 += av * c2[m];
                        s3 += av * c3[m];
                    }
                    f0[j] = s0;
                    f1[j] = s1;
                    f2[j] = s2;
                    f3[j] = s3;
                }
            }
            else
            {
                for (long p = p0; p < nbpix; p++)
                {
                    const float *restrict c_row = &mat_c_t[p * nbmvec];
                    float       *restrict f_row = &filt_out[p * mvecsize];
                    for (long j = 0; j < mvecsize; j++)
                    {
                        const float *restrict a_row = &a_inv[j * nbmvec];
                        float                 sum   = 0.0f;
                        #pragma omp simd reduction(+ : sum)
                        for (long m = 0; m < nbmvec; m++)
                        {
                            sum += a_row[m] * c_row[m];
                        }
                        f_row[j] = sum;
                    }
                }
            }
        }
        free(mat_c_t);
        return;
    }

    #pragma omp parallel for schedule(static)
    for (long p = 0; p < nbpix; p++)
    {
        for (long j = 0; j < mvecsize; j++)
        {
            const float *a_row = &a_inv[j * nbmvec];
            double       sum   = 0.0;
            for (long m = 0; m < nbmvec; m++)
            {
                sum += (double) a_row[m] * (double) mat_c[m * nbpix + p];
            }
            filt_out[p * mvecsize + j] = (float) sum;
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
    if (piston == NULL)
    {
        free(mask_idx);
        return -1;
    }

    imageID id_mata = create_2Dimage_ID("matA", d.mvecsize, d.nbmvec);
    imageID id_matc = create_2Dimage_ID("matC", nbpix, d.nbmvec);
    atmturb_linpred_populate_matrices(
        dcimg[id_wfin].array.F, dcimg[id_mata].array.F, dcimg[id_matc].array.F,
        mask_idx, piston, &d);
    free(piston);

    double eps = (RegLambda > 0.0 && RegLambda > SVDeps) ? RegLambda : SVDeps;
    linopt_compute_SVDpseudoInverse("matA", "matA_inv", eps, 10000, "matA_vt");

    char filtname[200];
    snprintf(filtname, sizeof(filtname), "Pfilt_%04d_%04ld", PForder, (long) (PFlag + 0.5f));
    create_2Dimage_ID(filtname, d.mvecsize, nbpix);
    imageID id_filt = image_ID(filtname);
    imageID id_ainv = image_ID("matA_inv");

    atmturb_linpred_multiply_ainv_c(
        dcimg[id_ainv].array.F, dcimg[id_matc].array.F, dcimg[id_filt].array.F,
        nbpix, d.mvecsize, d.nbmvec);

    delete_image_ID("matA");
    delete_image_ID("matC");
    delete_image_ID("matA_inv");
    delete_image_ID("matA_vt");
    free(mask_idx);
    return 0;
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

    for (long t0 = (long) PForder; t0 <= t_max; t0 += 16)
    {
        int blk = (t0 + 15 <= t_max) ? 16 : (int) (t_max - t0 + 1);
        atmturb_linpred_eval_tile16(&ctx, t0, blk);
    }

    free(hist_pack);
    free(mask_idx);
    return 0;
}
