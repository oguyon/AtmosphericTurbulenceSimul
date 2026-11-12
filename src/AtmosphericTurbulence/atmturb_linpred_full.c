// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_linpred_full.c
 * @brief   Full-aperture linear predictor construction and evaluation for predictive control
 */

#include <stdio.h>
#include <stdlib.h>

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
        for (long p = 0; p < nbpix; p++)
        {
            const float *restrict c_row = &mat_c_t[p * nbmvec];
            float       *restrict f_row = &filt_out[p * mvecsize];
            for (long j = 0; j < mvecsize; j++)
            {
                const float *restrict a_row = &a_inv[j * nbmvec];
                double                sum   = 0.0;
                #pragma omp simd reduction(+ : sum)
                for (long m = 0; m < nbmvec; m++)
                {
                    sum += (double) a_row[m] * (double) c_row[m];
                }
                f_row[j] = (float) sum;
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
    long  nxy      = nx * ny;
    long  nbpix    = 0;
    long *mask_idx = atmturb_extract_mask_indices(id_wfmask, nx, ny, &nbpix);
    if (mask_idx == NULL || nbpix <= 0)
    {
        free(mask_idx);
        return -1;
    }

    long         mvecsize = nbpix * (long) PForder;
    long         t_max    = nz - (long) PFlag;
    imageID      id_outp  = create_3Dimage_ID(WFoutp_name, nx, ny, nz);
    imageID      id_outf  = create_3Dimage_ID(WFoutf_name, nx, ny, nz);
    const float *wf_in    = dcimg[id_wfin].array.F;
    const float *filt     = dcimg[id_filt].array.F;
    float       *wf_outp  = dcimg[id_outp].array.F;
    float       *wf_outf  = dcimg[id_outf].array.F;

    #pragma omp parallel
    {
        float *hist_vec = (float *) malloc(sizeof(float) * (size_t) mvecsize);
        if (hist_vec != NULL)
        {
            #pragma omp for schedule(static)
            for (long t = (long) PForder; t <= t_max; t++)
            {
                for (int dt = 0; dt < PForder; dt++)
                {
                    const float *hslice = &wf_in[(t - (long) PForder + (long) dt) * nxy];
                    float       *hdst   = &hist_vec[(long) dt * nbpix];
                    for (long q = 0; q < nbpix; q++)
                    {
                        hdst[q] = hslice[mask_idx[q]];
                    }
                }

                long         target_off = (t - 1 + (long) PFlag) * nxy;
                const float *tslice     = &wf_in[target_off];
                float       *pslice     = &wf_outp[target_off];
                float       *fslice     = &wf_outf[target_off];

                for (long p = 0; p < nbpix; p++)
                {
                    const float *restrict f_row = &filt[p * mvecsize];
                    double                pred  = 0.0;
                    #pragma omp simd reduction(+ : pred)
                    for (long k = 0; k < mvecsize; k++)
                    {
                        pred += (double) f_row[k] * (double) hist_vec[k];
                    }
                    long  mp   = mask_idx[p];
                    float pf   = (float) pred;
                    pslice[mp] = pf;
                    fslice[mp] = tslice[mp] - pf;
                }
            }
            free(hist_vec);
        }
    }

    free(mask_idx);
    return 0;
}
