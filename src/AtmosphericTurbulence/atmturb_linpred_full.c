// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_linpred_full.c
 * @brief   Full-aperture linear predictor construction and evaluation for predictive control
 */

#include "AtmosphericTurbulence.h"
#include "atmturb_types.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * atmturb_compute_wavefront_piston - Compute per-frame mean piston non-destructively
 * @id_wf: Input wavefront cube image ID.
 * @mask_idx: Array of active mask pixel indices.
 * @nbpix: Number of active mask pixels.
 * @nz: Number of time slices.
 * @nxy: Frame pixel count (nx * ny).
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

    for (long kk = 0; kk < nz; kk++)
    {
        double tot = 0.0;
        const float *slice = &dcimg[id_wf].array.F[kk * nxy];
        for (long p = 0; p < nbpix; p++)
        {
            tot += (double) slice[mask_idx[p]];
        }
        piston[kk] = (float) (tot / (double) nbpix);
    }
    return piston;
}

/**
 * atmturb_extract_mask_indices - Collect active pupil pixel linear indices
 * @id_mask: Pupil mask image ID.
 * @nx: Mask width.
 * @ny: Mask height.
 * @nbpix_out: Output count of active pixels.
 *
 * Return: Allocated array of pixel indices.
 */
static long *atmturb_extract_mask_indices(
    imageID  id_mask,
    long     nx,
    long     ny,
    long    *nbpix_out)
{
    long count = 0;
    long nxy = nx * ny;
    for (long i = 0; i < nxy; i++)
    {
        if (dcimg[id_mask].array.F[i] > 0.5f)
        {
            count++;
        }
    }

    long *indices = malloc(sizeof(long) * (count > 0 ? count : 1));
    long idx = 0;
    for (long i = 0; i < nxy; i++)
    {
        if (dcimg[id_mask].array.F[i] > 0.5f)
        {
            indices[idx++] = i;
        }
    }

    *nbpix_out = count;
    return indices;
}

/**
 * AtmosphericTurbulence_Build_LinPredictor_Full - Build full-aperture AR linear prediction matrix
 * @WFin_name: Input wavefront cube name.
 * @WFmask_name: Pupil mask image name.
 * @PForder: Autoregressive predictor order (history steps).
 * @PFlag: Prediction horizon lag.
 * @SVDeps: Singular value cutoff tolerance.
 * @RegLambda: Tikhonov regularization parameter.
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
    imageID ID_WFin = image_ID(WFin_name);
    imageID ID_WFmask = image_ID(WFmask_name);
    if (ID_WFin < 0 || ID_WFmask < 0)
    {
        return -1;
    }

    long nx = dcimg[ID_WFin].md[0].size[0];
    long ny = dcimg[ID_WFin].md[0].size[1];
    long nz = dcimg[ID_WFin].md[0].size[2];

    long nbpix = 0;
    long *mask_idx = atmturb_extract_mask_indices(ID_WFmask, nx, ny, &nbpix);
    if (mask_idx == NULL || nbpix <= 0)
    {
        free(mask_idx);
        return -1;
    }

    float *piston = atmturb_compute_wavefront_piston(ID_WFin, mask_idx, nbpix, nz, nx * ny);
    if (piston == NULL)
    {
        free(mask_idx);
        return -1;
    }

    long mvecsize = nbpix * PForder;
    long nbmvec = nz - PForder - (long)PFlag;

    if (nbmvec <= 0 || mvecsize <= 0)
    {
        free(piston);
        free(mask_idx);
        return -1;
    }

    imageID IDmatA = create_2Dimage_ID("matA", mvecsize, nbmvec);
    imageID IDmatC = create_2Dimage_ID("matC", nbpix, nbmvec);

    for (long m = 0; m < nbmvec; m++)
    {
        for (long dt = 0; dt < PForder; dt++)
        {
            long t = m + dt;
            float p_t = piston[t];
            for (long p = 0; p < nbpix; p++)
            {
                dcimg[IDmatA].array.F[m * mvecsize + dt * nbpix + p] =
                    dcimg[ID_WFin].array.F[t * nx * ny + mask_idx[p]] - p_t;
            }
        }
        long t_target = m + PForder - 1 + (long)PFlag;
        float p_targ = piston[t_target];
        for (long p = 0; p < nbpix; p++)
        {
            dcimg[IDmatC].array.F[m * nbpix + p] =
                dcimg[ID_WFin].array.F[t_target * nx * ny + mask_idx[p]] - p_targ;
        }
    }
    free(piston);

    double eps = (RegLambda > 0.0 && RegLambda > SVDeps) ? RegLambda : SVDeps;
    linopt_compute_SVDpseudoInverse("matA", "matA_inv", eps, 10000, "matA_vt");

    char filtname[200];
    snprintf(filtname, sizeof(filtname), "Pfilt_%04d_%04ld", PForder, (long)(PFlag + 0.5));
    create_2Dimage_ID(filtname, mvecsize, nbpix);
    imageID ID_filt = image_ID(filtname);

    imageID ID_Ainv = image_ID("matA_inv");
    #pragma omp parallel for schedule(dynamic)
    for (long p = 0; p < nbpix; p++)
    {
        for (long j = 0; j < mvecsize; j++)
        {
            double sum = 0.0;
            const float *a_row = &dcimg[ID_Ainv].array.F[j * nbmvec];
            for (long m = 0; m < nbmvec; m++)
            {
                sum += (double) a_row[m] *
                       (double) dcimg[IDmatC].array.F[m * nbpix + p];
            }
            dcimg[ID_filt].array.F[p * mvecsize + j] = (float) sum;
        }
    }

    delete_image_ID("matA");
    delete_image_ID("matC");
    delete_image_ID("matA_inv");
    delete_image_ID("matA_vt");
    free(mask_idx);

    return 0;
}

/**
 * AtmosphericTurbulence_Apply_LinPredictor_Full - Evaluate linear prediction on series
 * @MODE: 0 for direct evaluation.
 * @WFin_name: Input wavefront cube.
 * @WFmask_name: Active pupil mask.
 * @PForder: Filter AR order.
 * @PFlag: Prediction horizon lag.
 * @WFoutp_name: Output predicted wavefront cube name.
 * @WFoutf_name: Output residual error cube name.
 *
 * Return: 0 on success.
 */
int AtmosphericTurbulence_Apply_LinPredictor_Full(
    int         MODE,
    const char *WFin_name,
    const char *WFmask_name,
    int         PForder,
    float       PFlag,
    const char *WFoutp_name,
    const char *WFoutf_name)
{
    (void)MODE;
    imageID ID_WFin = image_ID(WFin_name);
    imageID ID_WFmask = image_ID(WFmask_name);
    char filtname[200];
    snprintf(filtname, sizeof(filtname), "Pfilt_%04d_%04ld", PForder, (long)(PFlag + 0.5));
    imageID ID_filt = image_ID(filtname);

    if (ID_WFin < 0 || ID_WFmask < 0 || ID_filt < 0)
    {
        return -1;
    }

    long nx = dcimg[ID_WFin].md[0].size[0];
    long ny = dcimg[ID_WFin].md[0].size[1];
    long nz = dcimg[ID_WFin].md[0].size[2];

    long nbpix = 0;
    long *mask_idx = atmturb_extract_mask_indices(ID_WFmask, nx, ny, &nbpix);
    long mvecsize = nbpix * PForder;

    imageID ID_outp = create_3Dimage_ID(WFoutp_name, nx, ny, nz);
    imageID ID_outf = create_3Dimage_ID(WFoutf_name, nx, ny, nz);

    #pragma omp parallel for schedule(dynamic)
    for (long t = PForder; t <= nz - (long)PFlag; t++)
    {
        long target_t = t - 1 + (long)PFlag;
        for (long p = 0; p < nbpix; p++)
        {
            double pred = 0.0;
            for (long dt = 0; dt < PForder; dt++)
            {
                long hist_t = t - PForder + dt;
                for (long q = 0; q < nbpix; q++)
                {
                    pred += dcimg[ID_filt].array.F[p * mvecsize + dt * nbpix + q] *
                            dcimg[ID_WFin].array.F[hist_t * nx * ny + mask_idx[q]];
                }
            }
            dcimg[ID_outp].array.F[target_t * nx * ny + mask_idx[p]] = (float)pred;
            dcimg[ID_outf].array.F[target_t * nx * ny + mask_idx[p]] =
                dcimg[ID_WFin].array.F[target_t * nx * ny + mask_idx[p]] - (float)pred;
        }
    }

    free(mask_idx);
    return 0;
}
