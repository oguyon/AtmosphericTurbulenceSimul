// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_wfs_render_cuda.c
 * @brief   CUDA GPU dispatch wrappers for wavefront series batch rendering
 */

#ifdef HAVE_CUDA

#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "CommandLineInterface/CLIcore.h"
#include "atmturb_cuda.h"
#include "atmturb_cuda_rytov.h"
#include "atmturb_geometry.h"
#include "atmturb_lowfreq.h"
#include "atmturb_profile.h"
#include "atmturb_rolling.h"
#include "atmturb_rytov.h"
#include "atmturb_types.h"
#include "atmturb_wfs_cuda_common.h"
#include "atmturb_wfs_render.h"
#include "atmturb_wfs_render_cuda.h"

/**
 * atmturb_wfs_render_cuda - Render geometric simulation frames on GPU
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @master_size: Master screen dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 * @nbframes: Number of simulation frames.
 * @imgs: Container of output 3D image handles.
 *
 * Flattens atmospheric layers into constant-memory parameters and launches geometric
 * extrusion kernel across all frames simultaneously.
 *
 * Return: 0 on success, -1 on invalid configuration or CUDA error.
 */
int atmturb_wfs_render_cuda(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    long                        master_size,
    long                        pup_size,
    long                        nbframes,
    const atmturb_wfs_images_t *imgs)
{
    if (r == NULL || geom == NULL || r->lowfreq != 0 || CONF_ROLLING != 0 ||
        geom->oversample != 1 || geom->interp != 0)
    {
        return -1;
    }

    int nl = geom->nlayers;
    size_t m_sz = (size_t) nl * sizeof(float *);
    size_t d_sz = (size_t) nl * sizeof(double);
    void *buf = malloc(m_sz + 5 * d_sz);
    if (!buf)
    {
        return -1;
    }
    const float **masters = (const float **) buf;
    double *vx  = (double *) ((char *) buf + m_sz);
    double *vy  = vx + nl;
    double *cn2 = vy + nl;
    double *x0  = cn2 + nl;
    double *y0  = x0 + nl;

    for (int k = 0; k < nl; k++)
    {
        masters[k] = r->layers[k].screens[0].data;
        vx[k]      = geom->layers[k].vx_pix;
        vy[k]      = geom->layers[k].vy_pix;
        cn2[k]     = geom->layers[k].weight * geom->layers[k].weight;
        x0[k]      = geom->layers[k].x0;
        y0[k]      = geom->layers[k].y0;
    }

    double w0 = geom->layers[0].weight;
    double scoeff = (w0 > 0.0) ? (geom->layers[0].weight_s / w0) : 1.0;
    atmturb_cuda_sim_params_t cparams = {
        .nblayers  = nl,
        .msize     = master_size,
        .pup_size  = pup_size,
        .nbframes  = nbframes,
        .Scoeff    = scoeff,
        .h_masters = (const float *const *) masters,
        .vxpix     = vx,
        .vypix     = vy,
        .cn2       = cn2,
        .x0        = x0,
        .y0        = y0
    };
    long total_pixels = nbframes * pup_size * pup_size;
    float *amp_ptr  = &dcimg[imgs->ID_amp].array.F[0];
    float *samp_ptr = &dcimg[imgs->ID_samp].array.F[0];

    #pragma omp parallel for
    for (long i = 0; i < total_pixels; i++)
    {
        amp_ptr[i]  = 1.0f;
        samp_ptr[i] = 1.0f;
    }

    atmturb_cuda_sim_outputs_t outputs = {
        .pha  = &dcimg[imgs->ID_pha].array.F[0],
        .amp  = NULL,
        .spha = &dcimg[imgs->ID_spha].array.F[0],
        .samp = NULL
    };

    int ret = atmturb_wfs_render_frames_cuda(&cparams, &outputs);
    free(buf);
    return ret;
}

/**
 * atmturb_wfs_cuda_build_sublayers - Populate multi-frame sublayer arrays in parallel
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @plan: Active Rytov perturbation plan.
 * @screen_offsets: Screen index offsets per layer.
 * @nbframes: Number of simulation frames.
 * @time_step_s: Time step in seconds.
 * @max_sublayers: Maximum sublayers per superlayer.
 * @sublayers: Output contiguous sublayer array across all frames.
 */
static void atmturb_wfs_cuda_build_sublayers(
    const atmturb_rolling_t       *r,
    const atmturb_geom_t          *geom,
    const atmturb_rytov_plan_t    *plan,
    const int                     *screen_offsets,
    long                           nbframes,
    double                         time_step_s,
    int                            max_sublayers,
    atmturb_cuda_rytov_sublayer_t *sublayers)
{
    double offset_os = (double) plan->guard_pix * (double) geom->oversample;
    size_t frame_stride = (size_t) (plan->nsuper * max_sublayers);

    #pragma omp parallel for
    for (long t = 0; t < nbframes; t++)
    {
        atmturb_wfs_cuda_build_frame_sublayers(
            r, geom, plan->nsuper, plan->supers, screen_offsets,
            t, time_step_s, offset_os, max_sublayers,
            sublayers + (size_t) t * frame_stride);
    }
}

/**
 * atmturb_wfs_cuda_build_lf_modes - Populate multi-frame low-frequency mode arrays in parallel
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @plan: Active Rytov perturbation plan.
 * @nbframes: Number of simulation frames.
 * @time_step_s: Time step in seconds.
 * @max_lf_modes: Maximum low-frequency modes per superlayer.
 * @lf_modes: Output contiguous mode array across all frames.
 */
static void atmturb_wfs_cuda_build_lf_modes(
    const atmturb_rolling_t     *r,
    const atmturb_geom_t        *geom,
    const atmturb_rytov_plan_t  *plan,
    long                         nbframes,
    double                       time_step_s,
    int                          max_lf_modes,
    atmturb_cuda_lowfreq_mode_t *lf_modes)
{
    if (r->lowfreq == 0 || lf_modes == NULL)
    {
        return;
    }

    double offset_os = (double) plan->guard_pix * (double) geom->oversample;
    size_t frame_stride = (size_t) (plan->nsuper * max_lf_modes);

    #pragma omp parallel for
    for (long t = 0; t < nbframes; t++)
    {
        atmturb_wfs_cuda_build_frame_lf_modes(
            r, geom, plan->nsuper, plan->supers,
            t, time_step_s, offset_os, max_lf_modes,
            lf_modes + (size_t) t * frame_stride);
    }
}

/**
 * atmturb_wfs_cuda_init_rytov_lf - Allocate and initialize low-frequency mode buffers
 * @r: Rolling simulation context.
 * @plan: Active Rytov perturbation plan.
 * @nbframes: Number of simulation frames.
 * @out_max_lf: Output pointer for max mode count.
 * @out_super_nlf: Output pointer for mode counts per superlayer.
 * @out_lf_modes: Output pointer for pinned host buffer.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
static int atmturb_wfs_cuda_init_rytov_lf(
    const atmturb_rolling_t      *r,
    const atmturb_rytov_plan_t   *plan,
    long                          nbframes,
    int                          *out_max_lf,
    int                         **out_super_nlf,
    atmturb_cuda_lowfreq_mode_t **out_lf_modes)
{
    *out_max_lf    = 0;
    *out_super_nlf = NULL;
    *out_lf_modes  = NULL;

    if (r->lowfreq == 0)
    {
        return 0;
    }

    int max_lf = 0;
    for (int m = 0; m < plan->nsuper; m++)
    {
        int n_lf = plan->supers[m].nlayers * 24;
        if (n_lf > max_lf)
        {
            max_lf = n_lf;
        }
    }
    *out_max_lf = max_lf;

    if (max_lf > 0)
    {
        int *super_nlf = (int *) malloc(sizeof(int) * (size_t) plan->nsuper);
        size_t n_tot = (size_t) (nbframes * plan->nsuper * max_lf);
        size_t lf_bytes = n_tot * sizeof(atmturb_cuda_lowfreq_mode_t);
        atmturb_cuda_lowfreq_mode_t *lf = NULL;
        if (!super_nlf || atmturb_cuda_host_alloc((void **) &lf, lf_bytes) != 0)
        {
            free(super_nlf);
            return -1;
        }
        memset(lf, 0, lf_bytes);
        for (int m = 0; m < plan->nsuper; m++)
        {
            super_nlf[m] = plan->supers[m].nlayers * 24;
        }
        *out_super_nlf = super_nlf;
        *out_lf_modes  = lf;
    }
    return 0;
}

/**
 * struct atmturb_wfs_cuda_rytov_bufs_t - Memory container for batch Rytov rendering
 * @sublayers: Pre-allocated sublayer coordinate buffer.
 * @super_nlayers: Array of layer counts per superlayer.
 * @super_dist_m: Array of propagation distances per superlayer.
 * @masters: Host pointers to master screen data.
 * @screen_offsets: Screen index offsets per layer.
 * @lf_modes: Low-frequency Fourier mode buffer.
 * @super_nlf_modes: Mode counts per superlayer.
 */
typedef struct
{
    atmturb_cuda_rytov_sublayer_t *sublayers;
    int                           *super_nlayers;
    double                        *super_dist_m;
    const float                  **masters;
    int                           *screen_offsets;
    atmturb_cuda_lowfreq_mode_t   *lf_modes;
    int                           *super_nlf_modes;
} atmturb_wfs_cuda_rytov_bufs_t;

/**
 * atmturb_wfs_rytov_cuda_free_buffers - Free temporary buffers for batch Rytov rendering
 * @b: Pointer to buffer container.
 */
static void atmturb_wfs_rytov_cuda_free_buffers(
    atmturb_wfs_cuda_rytov_bufs_t *b)
{
    if (b->sublayers)
    {
        free(b->sublayers);
        b->sublayers = NULL;
    }
    free(b->super_nlayers);
    b->super_nlayers = NULL;
    free(b->super_dist_m);
    b->super_dist_m = NULL;
    free(b->masters);
    b->masters = NULL;
    free(b->screen_offsets);
    b->screen_offsets = NULL;
    if (b->lf_modes)
    {
        atmturb_cuda_host_free(b->lf_modes);
        b->lf_modes = NULL;
    }
    free(b->super_nlf_modes);
    b->super_nlf_modes = NULL;
}

/**
 * atmturb_wfs_rytov_cuda_alloc_buffers - Allocate temporary working buffers for batch Rytov
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @plan: Active Rytov perturbation plan.
 * @nbframes: Number of simulation frames.
 * @max_sublayers: Maximum sublayers per superlayer.
 * @out_max_lf: Output pointer for max low-frequency modes.
 * @out_total_screens: Output pointer for total master screen count.
 * @b: Destination buffer container.
 *
 * Return: 0 on success, -1 on failure.
 */
static int atmturb_wfs_rytov_cuda_alloc_buffers(
    const atmturb_rolling_t       *r,
    const atmturb_geom_t          *geom,
    const atmturb_rytov_plan_t    *plan,
    long                           nbframes,
    int                            max_sublayers,
    int                           *out_max_lf,
    int                           *out_total_screens,
    atmturb_wfs_cuda_rytov_bufs_t *b)
{
    memset(b, 0, sizeof(*b));

    size_t n_sub_tot = (size_t) (nbframes * plan->nsuper * max_sublayers);
    b->sublayers = (atmturb_cuda_rytov_sublayer_t *) calloc(
        n_sub_tot, sizeof(atmturb_cuda_rytov_sublayer_t));
    b->super_nlayers = (int *) malloc(sizeof(int) * (size_t) plan->nsuper);
    b->super_dist_m  = (double *) malloc(sizeof(double) * (size_t) plan->nsuper);

    if (atmturb_wfs_cuda_init_screen_pool(r, geom->nlayers, &b->masters,
                                         &b->screen_offsets, out_total_screens) != 0)
    {
        atmturb_wfs_rytov_cuda_free_buffers(b);
        return -1;
    }

    if (atmturb_wfs_cuda_init_rytov_lf(r, plan, nbframes, out_max_lf,
                                       &b->super_nlf_modes, &b->lf_modes) != 0)
    {
        atmturb_wfs_rytov_cuda_free_buffers(b);
        return -1;
    }

    if (!b->sublayers || !b->super_nlayers || !b->super_dist_m)
    {
        atmturb_wfs_rytov_cuda_free_buffers(b);
        return -1;
    }

    for (int m = 0; m < plan->nsuper; m++)
    {
        b->super_nlayers[m] = 2 * plan->supers[m].nlayers;
        b->super_dist_m[m]  = plan->supers[m].dist_m;
    }

    return 0;
}

/**
 * atmturb_wfs_render_rytov_cuda - Render Rytov diffractive simulation frames on GPU
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @prof: Active turbulence profile.
 * @params: Observation parameters container.
 * @master_size: Master screen dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 * @nbframes: Number of simulation frames.
 * @imgs: Container of output 3D image handles.
 *
 * Precomputes multi-frame sublayer trajectories and subharmonics on host threads with OpenMP,
 * then dispatches multi-layer GPU Rytov Fourier filtering and pupil assembly.
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_wfs_render_rytov_cuda(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    const atmturb_profile_t    *prof,
    const atmturb_obs_params_t *params,
    long                        master_size,
    long                        pup_size,
    long                        nbframes,
    const atmturb_wfs_images_t *imgs)
{
    if (r == NULL || geom == NULL || prof == NULL || params == NULL || imgs == NULL ||
        nbframes <= 0 || master_size <= 0 || pup_size <= 0 || !atmturb_cuda_device_available())
    {
        return -1;
    }

    atmturb_rytov_plan_t plan;
    double z_bin = (double) CONF_FRESNEL_PROPAGATION_BIN;
    if (atmturb_rytov_plan_init(&plan, prof, geom, pup_size, params->pupil_scale_m,
                                params->lambda_ref_m, params->lambda_s_m, z_bin,
                                CONF_FRESNEL_RYTOV_SEC_EXACT) != 0)
    {
        return -1;
    }

    int max_sublayers = 1;
    for (int m = 0; m < plan.nsuper; m++)
    {
        int n_sub = 2 * plan.supers[m].nlayers;
        if (n_sub > max_sublayers)
        {
            max_sublayers = n_sub;
        }
    }
    if (max_sublayers > 128)
    {
        atmturb_rytov_plan_free(&plan);
        return -1;
    }

    atmturb_wfs_cuda_rytov_bufs_t bufs;
    int max_lf_modes = 0;
    int total_screens = 0;
    if (atmturb_wfs_rytov_cuda_alloc_buffers(
            r, geom, &plan, nbframes, max_sublayers, &max_lf_modes,
            &total_screens, &bufs) != 0)
    {
        atmturb_rytov_plan_free(&plan);
        return -1;
    }

    atmturb_wfs_cuda_build_sublayers(
        r, geom, &plan, bufs.screen_offsets, nbframes, params->time_step_s,
        max_sublayers, bufs.sublayers);
    atmturb_wfs_cuda_build_lf_modes(
        r, geom, &plan, nbframes, params->time_step_s,
        max_lf_modes, bufs.lf_modes);

    atmturb_cuda_rytov_params_t cparams = {
        .msize           = master_size,
        .pup_size        = pup_size,
        .pad_size        = plan.pad_size,
        .guard_pix       = plan.guard_pix,
        .nbframes        = nbframes,
        .nblayers        = total_screens,
        .nsuper          = plan.nsuper,
        .has_sec         = (plan.lambda_s_m > 0.0) ? 1 : 0,
        .sec_shared      = plan.sec_shared,
        .os              = (int) geom->oversample,
        .interp          = geom->interp,
        .use_moisan      = plan.use_moisan,
        .lowfreq         = r->lowfreq,
        .h_masters       = bufs.masters,
        .h_laplace_inv   = plan.laplace_inv,
        .h_exp_x         = (const void *) plan.exp_x,
        .h_exp_y         = (const void *) plan.exp_y,
        .h_filt_a_pri    = (const float *const *) plan.filter_a_pri,
        .h_filt_b_pri    = (const float *const *) plan.filter_b_pri,
        .h_filt_a_sec    = (const float *const *) plan.filter_a_sec,
        .h_filt_b_sec    = (const float *const *) plan.filter_b_sec,
        .h_chrom_ramp    = (const void *const *) plan.chrom_ramp,
        .frame_sublayers = bufs.sublayers,
        .super_nlayers   = bufs.super_nlayers,
        .super_dist_m    = bufs.super_dist_m,
        .max_sublayers   = max_sublayers,
        .frame_lf_modes  = bufs.lf_modes,
        .super_nlf_modes = bufs.super_nlf_modes,
        .max_lf_modes    = max_lf_modes,
        .pha             = (imgs->ID_pha >= 0) ? &dcimg[imgs->ID_pha].array.F[0] : NULL,
        .amp             = (imgs->ID_amp >= 0) ? &dcimg[imgs->ID_amp].array.F[0] : NULL,
        .spha            = (imgs->ID_spha >= 0) ? &dcimg[imgs->ID_spha].array.F[0] : NULL,
        .samp            = (imgs->ID_samp >= 0) ? &dcimg[imgs->ID_samp].array.F[0] : NULL
    };

    int ret = atmturb_cuda_rytov_render(&cparams);

    atmturb_wfs_rytov_cuda_free_buffers(&bufs);
    atmturb_rytov_plan_free(&plan);
    return ret;
}

#endif // HAVE_CUDA
