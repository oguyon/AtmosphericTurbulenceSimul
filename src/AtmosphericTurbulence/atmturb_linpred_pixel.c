// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_linpred_pixel.c
 * @brief   Pixel-level and shift-invariant kernel linear prediction filters
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "AtmosphericTurbulence.h"
#include "atmturb_linpred_solve.h"
#include "atmturb_types.h"

/**
 * struct linpred_kern_grid_t - Spatial index layout for 2D shift-invariant kernel ops
 */
typedef struct
{
    const long  *pix_x;
    const long  *pix_y;
    const long  *pix_map;
    const float *src;
    long         nx;
    long         ny;
    long         nbpix;
    long         krad;
    long         ksize;
    long         pforder;
    int          is_2d;
} linpred_kern_grid_t;

/**
 * atmturb_linpred_build_pix_map - Build 2D coordinate lookup map for active pixels
 * @id_mask: Pupil mask image ID.
 * @nx: Pupil mask width.
 * @ny: Pupil mask height.
 * @nbpix: Maximum number of active pixels.
 * @pix_x: Output array for X coordinates.
 * @pix_y: Output array for Y coordinates.
 * @count_out: Output count of populated active pixels.
 *
 * Return: Allocated 2D pixel index map (size nx * ny), or NULL on failure.
 */
static long *atmturb_linpred_build_pix_map(
    imageID  id_mask,
    long     nx,
    long     ny,
    long     nbpix,
    long    *pix_x,
    long    *pix_y,
    long    *count_out)
{
    long *pix_map = (long *) malloc(sizeof(long) * (size_t) (nx * ny));
    if (pix_map == NULL)
    {
        return NULL;
    }
    for (long i = 0; i < nx * ny; i++)
    {
        pix_map[i] = -1;
    }

    long count = 0;
    for (long ii = 0; ii < nx; ii++)
    {
        for (long jj = 0; jj < ny; jj++)
        {
            if (dcimg[id_mask].array.F[jj * nx + ii] > 0.5f && count < nbpix)
            {
                pix_x[count]          = ii;
                pix_y[count]          = jj;
                pix_map[jj * nx + ii] = count;
                count++;
            }
        }
    }
    *count_out = count;
    return pix_map;
}

/**
 * atmturb_linpred_accum_pix_kern - Accumulate one target pixel's neighborhood into local kernel
 */
static void atmturb_linpred_accum_pix_kern(
    const linpred_kern_grid_t *g,
    long                       p,
    float                     *k_acc,
    float                     *c_acc)
{
    long i0  = g->pix_x[p];
    long j0  = g->pix_y[p];
    long ksq = g->ksize * g->ksize;

    for (long dj = -g->krad; dj <= g->krad; dj++)
    {
        long y = j0 + dj;
        if (y < 0 || y >= g->ny)
        {
            continue;
        }
        long kj      = dj + g->krad;
        long row_off = y * g->nx;

        for (long di = -g->krad; di <= g->krad; di++)
        {
            long x = i0 + di;
            long q = (x >= 0 && x < g->nx) ? g->pix_map[row_off + x] : -1;
            if (q < 0)
            {
                continue;
            }
            long ki = di + g->krad;
            for (long dt = 0; dt < g->pforder; dt++)
            {
                long kidx = dt * ksq + kj * g->ksize + ki;
                long fidx = g->is_2d ? (p * (g->pforder * g->nbpix) + dt * g->nbpix + q)
                                     : (dt * g->nbpix * g->nbpix + p * g->nbpix + q);
                k_acc[kidx] += g->src[fidx];
                c_acc[kidx] += 1.0f;
            }
        }
    }
}

/**
 * AtmosphericTurbulence_LinPredictor_filt_2DKernelExtract - Extract shift-invariant 2D kernels
 * @IDfilt_name: Input filter matrix image name.
 * @IDmask_name: Active pupil mask image name.
 * @krad: Extraction neighborhood radius in pixels.
 * @IDkern_name: Output 3D kernel image name.
 *
 * Return: Output image ID on success.
 */
long AtmosphericTurbulence_LinPredictor_filt_2DKernelExtract(
    const char *IDfilt_name,
    const char *IDmask_name,
    long        krad,
    const char *IDkern_name)
{
    imageID IDfilt = image_ID(IDfilt_name);
    imageID IDmask = image_ID(IDmask_name);
    if (IDfilt < 0 || IDmask < 0)
    {
        return -1;
    }

    long  nx      = dcimg[IDmask].md[0].size[0];
    long  ny      = dcimg[IDmask].md[0].size[1];
    int   is_2d   = (dcimg[IDfilt].md[0].naxis <= 2);
    long  nbpix   = is_2d ? dcimg[IDfilt].md[0].size[1] : dcimg[IDfilt].md[0].size[0];
    long  pforder = is_2d ? (dcimg[IDfilt].md[0].size[0] / (nbpix > 0 ? nbpix : 1))
                          : ((dcimg[IDfilt].md[0].naxis > 2) ? dcimg[IDfilt].md[0].size[2] : 1);
    long *pix_x   = (long *) malloc(sizeof(long) * (nbpix > 0 ? (size_t) nbpix : 1));
    long *pix_y   = (long *) malloc(sizeof(long) * (nbpix > 0 ? (size_t) nbpix : 1));
    long  count   = 0;
    long *pix_map = (pix_x != NULL && pix_y != NULL && nbpix > 0)
                        ? atmturb_linpred_build_pix_map(
                              IDmask, nx, ny, nbpix, pix_x, pix_y, &count)
                        : NULL;
    if (pix_map == NULL)
    {
        free(pix_x);
        free(pix_y);
        return -1;
    }

    long    ksize  = 2 * krad + 1;
    long    ktot   = ksize * ksize * pforder;
    imageID IDkern = create_3Dimage_ID(IDkern_name, ksize, ksize, pforder);
    float  *kern_f = dcimg[IDkern].array.F;
    float  *cnt_f  = (float *) calloc((size_t) ktot, sizeof(float));
    linpred_kern_grid_t g = {
        .pix_x = pix_x, .pix_y = pix_y, .pix_map = pix_map, .src = dcimg[IDfilt].array.F,
        .nx = nx, .ny = ny, .nbpix = nbpix, .krad = krad, .ksize = ksize,
        .pforder = pforder, .is_2d = is_2d
    };

    if (cnt_f != NULL)
    {
        #pragma omp parallel
        {
            float *kloc = (float *) calloc((size_t) ktot, sizeof(float));
            float *cloc = (float *) calloc((size_t) ktot, sizeof(float));
            if (kloc != NULL && cloc != NULL)
            {
                #pragma omp for schedule(static)
                for (long p = 0; p < count; p++)
                {
                    atmturb_linpred_accum_pix_kern(&g, p, kloc, cloc);
                }
                #pragma omp critical
                {
                    for (long i = 0; i < ktot; i++)
                    {
                        kern_f[i] += kloc[i];
                        cnt_f[i]  += cloc[i];
                    }
                }
            }
            free(kloc);
            free(cloc);
        }
        for (long i = 0; i < ktot; i++)
        {
            if (cnt_f[i] > 0.5f)
            {
                kern_f[i] /= cnt_f[i];
            }
        }
        free(cnt_f);
    }

    free(pix_map);
    free(pix_x);
    free(pix_y);
    return IDkern;
}

/**
 * atmturb_linpred_expand_pix - Scatter shift-invariant 2D kernel weights for target pixel p
 */
static void atmturb_linpred_expand_pix(
    const linpred_kern_grid_t *g,
    long                       p,
    float                     *dst)
{
    long i0  = g->pix_x[p];
    long j0  = g->pix_y[p];
    long ksq = g->ksize * g->ksize;
    long nsq = g->nbpix * g->nbpix;

    for (long dj = -g->krad; dj <= g->krad; dj++)
    {
        long y = j0 + dj;
        if (y < 0 || y >= g->ny)
        {
            continue;
        }
        long kj      = dj + g->krad;
        long row_off = y * g->nx;

        for (long di = -g->krad; di <= g->krad; di++)
        {
            long x = i0 + di;
            long q = (x >= 0 && x < g->nx) ? g->pix_map[row_off + x] : -1;
            if (q < 0)
            {
                continue;
            }
            long ki = di + g->krad;
            for (long dt = 0; dt < g->pforder; dt++)
            {
                dst[dt * nsq + p * g->nbpix + q] = g->src[dt * ksq + kj * g->ksize + ki];
            }
        }
    }
}

/**
 * AtmosphericTurbulence_LinPredictor_filt_Expand - Expand 2D kernel across full aperture
 * @IDfilt_name: Input 2D/3D shift-invariant kernel.
 * @IDmask_name: Active pupil mask.
 *
 * Return: Output expanded matrix image ID.
 */
long AtmosphericTurbulence_LinPredictor_filt_Expand(
    const char *IDfilt_name,
    const char *IDmask_name)
{
    imageID IDfilt = image_ID(IDfilt_name);
    imageID IDmask = image_ID(IDmask_name);
    if (IDfilt < 0 || IDmask < 0)
    {
        return -1;
    }

    long nx      = dcimg[IDmask].md[0].size[0];
    long ny      = dcimg[IDmask].md[0].size[1];
    long ksize   = dcimg[IDfilt].md[0].size[0];
    long krad    = ksize / 2;
    long pforder = (dcimg[IDfilt].md[0].naxis > 2) ? dcimg[IDfilt].md[0].size[2] : 1;
    long nbpix   = 0;
    for (long i = 0; i < nx * ny; i++)
    {
        if (dcimg[IDmask].array.F[i] > 0.5f)
        {
            nbpix++;
        }
    }

    long *pix_x   = (long *) malloc(sizeof(long) * (nbpix > 0 ? (size_t) nbpix : 1));
    long *pix_y   = (long *) malloc(sizeof(long) * (nbpix > 0 ? (size_t) nbpix : 1));
    long  count   = 0;
    long *pix_map = (pix_x != NULL && pix_y != NULL && nbpix > 0)
                        ? atmturb_linpred_build_pix_map(
                              IDmask, nx, ny, nbpix, pix_x, pix_y, &count)
                        : NULL;
    if (pix_map == NULL)
    {
        free(pix_x);
        free(pix_y);
        return -1;
    }

    imageID             IDout = create_3Dimage_ID("Pfilt_exp", nbpix, nbpix, pforder);
    float              *dst   = dcimg[IDout].array.F;
    linpred_kern_grid_t g     = {
        .pix_x = pix_x, .pix_y = pix_y, .pix_map = pix_map, .src = dcimg[IDfilt].array.F,
        .nx = nx, .ny = ny, .nbpix = nbpix, .krad = krad, .ksize = ksize,
        .pforder = pforder, .is_2d = 0
    };

    #pragma omp parallel for schedule(static)
    for (long p = 0; p < count; p++)
    {
        atmturb_linpred_expand_pix(&g, p, dst);
    }

    free(pix_map);
    free(pix_x);
    free(pix_y);
    return IDout;
}

/**
 * struct linpred_pix_cfg_t - Stencil geometry for single-pixel AR predictor assembly
 */
typedef struct
{
    const float *wf;
    long         nx;
    long         ny;
    long         nz;
    long         lag;
    long         nbstep;
    long         xyrad;
    long         iipix;
    long         jjpix;
    long         kdiam;
    long         n_inputs;
} linpred_pix_cfg_t;

/**
 * atmturb_linpred_fill_pix_sample - Populate one training row for single-pixel AR predictor
 */
static void atmturb_linpred_fill_pix_sample(
    const linpred_pix_cfg_t *cfg,
    long                     s,
    float                   *a_row,
    float                   *b_out)
{
    long nxy = cfg->nx * cfg->ny;
    long ksq = cfg->kdiam * cfg->kdiam;

    for (long dt = 0; dt < cfg->nbstep; dt++)
    {
        long t = s + dt;
        for (long dj = -cfg->xyrad; dj <= cfg->xyrad; dj++)
        {
            long py = cfg->jjpix + dj;
            for (long di = -cfg->xyrad; di <= cfg->xyrad; di++)
            {
                long  px  = cfg->iipix + di;
                long  col = dt * ksq + (dj + cfg->xyrad) * cfg->kdiam + (di + cfg->xyrad);
                float val = 0.0f;
                if (px >= 0 && px < cfg->nx && py >= 0 && py < cfg->ny && t < cfg->nz)
                {
                    val = cfg->wf[t * nxy + py * cfg->nx + px];
                }
                a_row[col] = val;
            }
        }
    }

    long  t_targ = s + cfg->nbstep - 1 + cfg->lag;
    float targ   = 0.0f;
    if (cfg->iipix >= 0 && cfg->iipix < cfg->nx &&
        cfg->jjpix >= 0 && cfg->jjpix < cfg->ny && t_targ < cfg->nz)
    {
        targ = cfg->wf[t_targ * nxy + cfg->jjpix * cfg->nx + cfg->iipix];
    }
    *b_out = targ;
}

/**
 * AtmosphericTurbulence_Build_LinPredictor - Train linear AR predictor for a single pixel
 * @NB_WFstep: Number of wavefront training steps.
 * @WFphaNoise: Measurement phase noise standard deviation.
 * @WFPlag: Prediction horizon lag.
 * @WFP_NBstep: History steps (AR order).
 * @WFP_xyrad: Spatial footprint radius around target pixel.
 * @WFPiipix: Target pixel X coordinate.
 * @WFPjjpix: Target pixel Y coordinate.
 * @slambdaum: Observing wavelength in um.
 *
 * Return: 0 on success, -1 on failure.
 */
int AtmosphericTurbulence_Build_LinPredictor(
    long   NB_WFstep,
    double WFphaNoise,
    long   WFPlag,
    long   WFP_NBstep,
    long   WFP_xyrad,
    long   WFPiipix,
    long   WFPjjpix,
    float  slambdaum)
{
    (void) WFphaNoise;
    (void) slambdaum;

    long    kdiam     = 2 * WFP_xyrad + 1;
    long    n_inputs  = kdiam * kdiam * WFP_NBstep;
    long    n_samples = NB_WFstep - WFP_NBstep - WFPlag;
    imageID ID_wf     = image_ID("outarraypha");
    long    nx        = (ID_wf >= 0) ? dcimg[ID_wf].md[0].size[0] : 0;
    long    ny        = (ID_wf >= 0) ? dcimg[ID_wf].md[0].size[1] : 0;
    long    nz        = (ID_wf >= 0) ? dcimg[ID_wf].md[0].size[2] : 0;
    if (ID_wf < 0 || nz < NB_WFstep || n_samples <= 0 || n_inputs <= 0)
    {
        return -1;
    }

    float *mat_a = (float *) malloc(sizeof(float) * (size_t) n_samples * (size_t) n_inputs);
    float *mat_b = (float *) malloc(sizeof(float) * (size_t) n_samples);
    if (mat_a == NULL || mat_b == NULL)
    {
        free(mat_a);
        free(mat_b);
        return -1;
    }

    linpred_pix_cfg_t cfg = {
        .wf = dcimg[ID_wf].array.F, .nx = nx, .ny = ny, .nz = nz,
        .lag = WFPlag, .nbstep = WFP_NBstep, .xyrad = WFP_xyrad,
        .iipix = WFPiipix, .jjpix = WFPjjpix, .kdiam = kdiam, .n_inputs = n_inputs
    };
    #pragma omp parallel for schedule(static)
    for (long s = 0; s < n_samples; s++)
    {
        atmturb_linpred_fill_pix_sample(&cfg, s, &mat_a[s * n_inputs], &mat_b[s]);
    }

    char out_name[200];
    snprintf(out_name, sizeof(out_name), "WFPfilt_%ld_%ld", WFPiipix, WFPjjpix);
    imageID id_old = image_ID(out_name);
    if (id_old >= 0)
    {
        delete_image_ID(out_name);
    }
    imageID ID_filt = create_2Dimage_ID(out_name, n_inputs, 1);
    int     ret     = -1;
    if (ID_filt >= 0)
    {
        atmturb_linpred_problem_t prob = {
            .mat_a    = mat_a,
            .mat_c    = mat_b,
            .filt_out = dcimg[ID_filt].array.F,
            .nbmvec   = n_samples,
            .mvecsize = n_inputs,
            .nbpix    = 1,
            .reg_eps  = 1e-3
        };
        ret = atmturb_linpred_solve_normal(&prob);
    }

    free(mat_a);
    free(mat_b);
    return ret;
}
