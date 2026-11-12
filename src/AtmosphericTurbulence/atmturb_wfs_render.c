// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_wfs_render.c
 * @brief   Wavefront series rendering engines (geometric, diffractive, CUDA)
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "AtmosphereModel/AtmosphereModel.h"
#include "AtmosphericTurbulence.h"
#include "atmturb_fresnel.h"
#include "atmturb_geometry.h"
#include "atmturb_lowfreq.h"
#include "atmturb_profile.h"
#include "atmturb_rolling.h"
#include "atmturb_rytov.h"
#include "atmturb_simd.h"
#include "atmturb_types.h"
#include "atmturb_wfs_render.h"

#ifdef HAVE_CUDA
#include "atmturb_wfs_render_cuda.h"
#endif

/**
 * atmturb_wfs_extrude_channel - Extrude one or two rolling epochs into pupil
 * @rev: Evaluated rolling frame state (scrA, scrB, wA, wB).
 * @master_size: Master screen linear dimension.
 * @x: Unwrapped X position in master pixels.
 * @y: Unwrapped Y position in master pixels.
 * @geom_pup: Packed [pup_size, oversample, interp] integers.
 * @w_layer: Base layer amplitude weight.
 * @out_pha: Destination phase array to accumulate into.
 */
static inline void atmturb_wfs_extrude_channel(
    const atmturb_rolling_eval_t *rev,
    long                          master_size,
    double                        x,
    double                        y,
    const atmturb_geom_t         *geom,
    long                          pup_size,
    double                        w_layer,
    float                        *out_pha)
{
    double cur_x = fmod(x, (double) master_size);
    if (cur_x < 0.0)
    {
        cur_x += (double) master_size;
    }
    double cur_y = fmod(y, (double) master_size);
    if (cur_y < 0.0)
    {
        cur_y += (double) master_size;
    }

    int has_b = (rev->wB > 0.0f && rev->scrB != NULL);
    atmturb_extrude_params_t ep = {
        .master   = rev->scrA,
        .msize    = master_size,
        .x0       = cur_x,
        .y0       = cur_y,
        .pup_size = pup_size,
        .os       = (long) geom->oversample,
        .interp   = geom->interp,
        .weight   = (float) (w_layer * (double) rev->wA),
        .out_pha  = out_pha,
        .master_b = has_b ? rev->scrB : NULL,
        .weight_b = has_b ? (float) (w_layer * (double) rev->wB) : 0.0f
    };
    atmturb_extrude_accumulate(&ep);
}

/**
 * atmturb_wfs_render_layer_step - Extrude one layer and optionally fuse subharmonics
 */
static void atmturb_wfs_render_layer_step(
    const atmturb_rolling_t           *r,
    const atmturb_geom_t              *geom,
    int                                k,
    long                               t,
    double                             time_step_s,
    long                               master_size,
    const atmturb_wfs_render_target_t *target,
    atmturb_lf_fused_t                *lf_pri,
    atmturb_lf_fused_t                *lf_sec)
{
    long   pad_size  = target->pup_size + 2 * target->guard_pix;
    double offset_os = (double) target->guard_pix * (double) geom->oversample;

    const atmturb_layer_geom_t *lg = &geom->layers[k];
    atmturb_rolling_eval_t      rev;
    atmturb_rolling_get_frame(r, k, t, time_step_s, &rev);

    double dx     = atmturb_geom_get_layer_dx(lg, t);
    double dy     = atmturb_geom_get_layer_dy(lg, t);
    double wscale = (target->weight_scale > 0.0f) ? (double) target->weight_scale : 1.0;
    double w_pri  = lg->weight * wscale;
    double w_sec  = lg->weight_s * wscale;

    double x = lg->x0 + dx - offset_os;
    double y = lg->y0 + dy - offset_os;
    atmturb_wfs_extrude_channel(&rev, master_size, x, y, geom, pad_size, w_pri, target->pha);

    if (target->spha != NULL)
    {
        double xs = lg->xs0 + dx - offset_os;
        double ys = lg->ys0 + dy - offset_os;
        atmturb_wfs_extrude_channel(&rev, master_size, xs, ys, geom, pad_size,
                                    w_sec, target->spha);
        if (r->lowfreq)
        {
            if (lf_sec != NULL)
            {
                atmturb_lowfreq_params_t lp = {
                    .lf         = &r->layers[k].lf_base,
                    .custom_are = rev.are_eff,
                    .custom_aim = rev.aim_eff,
                    .x0         = xs,
                    .y0         = ys,
                    .weight     = (float) w_sec
                };
                atmturb_lowfreq_fuse(lf_sec, &lp);
            }
            else
            {
                atmturb_lowfreq_accumulate_custom(&r->layers[k].lf_base, rev.are_eff, rev.aim_eff,
                                                  xs, ys, pad_size, (long) geom->oversample,
                                                  (float) w_sec, target->spha);
            }
        }
    }

    if (r->lowfreq)
    {
        if (lf_pri != NULL)
        {
            atmturb_lowfreq_params_t lp = {
                .lf         = &r->layers[k].lf_base,
                .custom_are = rev.are_eff,
                .custom_aim = rev.aim_eff,
                .x0         = x,
                .y0         = y,
                .weight     = (float) w_pri
            };
            atmturb_lowfreq_fuse(lf_pri, &lp);
        }
        else
        {
            atmturb_lowfreq_accumulate_custom(&r->layers[k].lf_base, rev.are_eff, rev.aim_eff,
                                              x, y, pad_size, (long) geom->oversample,
                                              (float) w_pri, target->pha);
        }
    }
}

/**
 * atmturb_wfs_render_layer_target - Render one turbulence layer into target buffer
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @k: Layer index.
 * @t: Frame index.
 * @time_step_s: Time step in seconds.
 * @master_size: Master screen dimension.
 * @target: Render target specifications and buffers.
 */
void atmturb_wfs_render_layer_target(
    const atmturb_rolling_t           *r,
    const atmturb_geom_t              *geom,
    int                                k,
    long                               t,
    double                             time_step_s,
    long                               master_size,
    const atmturb_wfs_render_target_t *target)
{
    atmturb_wfs_render_layer_step(r, geom, k, t, time_step_s, master_size, target, NULL, NULL);
}

/**
 * atmturb_wfs_render_all_layers - Render all geometric layers with fused epochs & subharmonics
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @t: Frame index.
 * @time_step_s: Time step in seconds.
 * @master_size: Master screen dimension.
 * @target: Render target specifications and buffers.
 */
void atmturb_wfs_render_all_layers(
    const atmturb_rolling_t           *r,
    const atmturb_geom_t              *geom,
    long                               t,
    double                             time_step_s,
    long                               master_size,
    const atmturb_wfs_render_target_t *target)
{
    atmturb_lf_fused_t lf_pri = {0};
    atmturb_lf_fused_t lf_sec = {0};

    for (int k = 0; k < geom->nlayers; k++)
    {
        atmturb_wfs_render_layer_step(r, geom, k, t, time_step_s, master_size,
                                      target, &lf_pri, &lf_sec);
    }

    if (r->lowfreq && geom->nlayers > 0)
    {
        long pad_size = target->pup_size + 2 * target->guard_pix;
        long os       = (long) geom->oversample;
        atmturb_lowfreq_flush(&lf_pri, &r->layers[0].lf_base, pad_size, os, target->pha);
        if (target->spha != NULL)
        {
            atmturb_lowfreq_flush(&lf_sec, &r->layers[0].lf_base, pad_size, os, target->spha);
        }
    }
}

/**
 * atmturb_wfs_render_superlayer - Render all layers in a superlayer with fused subharmonics
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @sl: Superlayer descriptor (layer indices and fractional weights).
 * @t: Frame index.
 * @time_step_s: Time step in seconds.
 * @target: Render target specifications and buffers.
 */
void atmturb_wfs_render_superlayer(
    const atmturb_rolling_t           *r,
    const atmturb_geom_t              *geom,
    const atmturb_superlayer_t        *sl,
    long                               t,
    double                             time_step_s,
    const atmturb_wfs_render_target_t *target)
{
    atmturb_lf_fused_t          lf_pri = {0};
    atmturb_lf_fused_t          lf_sec = {0};
    atmturb_wfs_render_target_t sub_t  = *target;

    for (int j = 0; j < sl->nlayers; j++)
    {
        int k              = sl->layer_indices[j];
        sub_t.weight_scale = (sl->layer_weights != NULL) ? sl->layer_weights[j] : 1.0f;
        atmturb_wfs_render_layer_step(r, geom, k, t, time_step_s, r->msize,
                                      &sub_t, &lf_pri, &lf_sec);
    }

    if (r->lowfreq && sl->nlayers > 0)
    {
        int  k0       = sl->layer_indices[0];
        long pad_size = target->pup_size + 2 * target->guard_pix;
        long os       = (long) geom->oversample;
        atmturb_lowfreq_flush(&lf_pri, &r->layers[k0].lf_base, pad_size, os, target->pha);
        if (target->spha != NULL)
        {
            atmturb_lowfreq_flush(&lf_sec, &r->layers[k0].lf_base, pad_size, os, target->spha);
        }
    }
}

/**
 * atmturb_wfs_render_layer - Render one turbulence layer into phase slices for frame t
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @k: Layer index.
 * @t: Frame index.
 * @time_step_s: Time step in seconds.
 * @master_size: Master screen dimension.
 * @pup_size: Pupil dimension.
 * @pha_slice: Primary phase frame slice.
 * @spha_slice: Secondary phase frame slice.
 */
void atmturb_wfs_render_layer(
    const atmturb_rolling_t *r,
    const atmturb_geom_t    *geom,
    int                      k,
    long                     t,
    double                   time_step_s,
    long                     master_size,
    long                     pup_size,
    float                   *pha_slice,
    float                   *spha_slice)
{
    atmturb_wfs_render_target_t target;
    target.pup_size     = pup_size;
    target.guard_pix    = 0;
    target.pha          = pha_slice;
    target.spha         = spha_slice;
    target.weight_scale = 1.0f;

    atmturb_wfs_render_layer_target(r, geom, k, t, time_step_s, master_size, &target);
}

/**
 * atmturb_wfs_render_geometric - Multi-threaded geometric rendering of simulation time steps
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @master_size: Master screen dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 * @nbframes: Number of simulation frames.
 * @time_step_s: Time step between frames in seconds.
 * @imgs: Container of output 3D image handles.
 */
static void atmturb_wfs_render_geometric(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    long                        master_size,
    long                        pup_size,
    long                        nbframes,
    double                      time_step_s,
    const atmturb_wfs_images_t *imgs)
{
    long frame_pixels = pup_size * pup_size;

    #pragma omp parallel for schedule(dynamic, 2)
    for (long t = 0; t < nbframes; t++)
    {
        long slice = t * frame_pixels;
        float *pha_slice  = &dcimg[imgs->ID_pha].array.F[slice];
        float *amp_slice  = (imgs->ID_amp >= 0) ? &dcimg[imgs->ID_amp].array.F[slice] : NULL;
        float *spha_slice = (imgs->ID_spha >= 0) ? &dcimg[imgs->ID_spha].array.F[slice] : NULL;
        float *samp_slice = (imgs->ID_samp >= 0) ? &dcimg[imgs->ID_samp].array.F[slice] : NULL;

        if (amp_slice != NULL)
        {
            atmturb_init_phase_amp(pha_slice, amp_slice, frame_pixels);
        }
        else
        {
            memset(pha_slice, 0, sizeof(float) * (size_t) frame_pixels);
        }

        if (spha_slice != NULL)
        {
            if (samp_slice != NULL)
            {
                atmturb_init_phase_amp(spha_slice, samp_slice, frame_pixels);
            }
            else
            {
                memset(spha_slice, 0, sizeof(float) * (size_t) frame_pixels);
            }
        }

        atmturb_wfs_render_target_t target = {
            .pup_size     = pup_size,
            .guard_pix    = 0,
            .pha          = pha_slice,
            .spha         = spha_slice,
            .weight_scale = 1.0f
        };
        atmturb_wfs_render_all_layers(r, geom, t, time_step_s, master_size, &target);
    }
}

/**
 * atmturb_wfs_pick_threads - Choose outer frame threads and inner FFTW threads
 */
static void atmturb_wfs_pick_threads(
    long grid_size,
    long nbframes,
    int *out_nth,
    int *out_fftw_nth)
{
    int max_t = 1;
#ifdef _OPENMP
    max_t = omp_get_max_threads();
#endif
    int cap   = (grid_size >= 256 && max_t > 16) ? 16 : max_t;
    int nth   = (nbframes > 0 && nbframes < (long) cap) ? (int) nbframes : cap;
    nth       = (nth < 1) ? 1 : nth;
    int f_nth = (nbframes < 4) ? (max_t / nth) : 1;
    f_nth     = (f_nth < 1) ? 1 : ((f_nth > 8) ? 8 : f_nth);

    *out_nth      = nth;
    *out_fftw_nth = f_nth;
}

/**
 * atmturb_wfs_render_diffractive - Multi-threaded diffractive rendering of simulation time steps
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @prof: Active turbulence profile.
 * @params: Observation parameters container.
 * @master_size: Master screen dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 * @nbframes: Number of simulation frames.
 * @imgs: Container of output 3D image handles.
 *
 * Return: 0 on success, -1 on plan or context allocation failure.
 */
static int atmturb_wfs_render_diffractive(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    const atmturb_profile_t    *prof,
    const atmturb_obs_params_t *params,
    long                        master_size,
    long                        pup_size,
    long                        nbframes,
    const atmturb_wfs_images_t *imgs)
{
    atmturb_fresnel_plan_t plan;
    double z_bin = (double) CONF_FRESNEL_PROPAGATION_BIN;
    if (atmturb_fresnel_plan_init(&plan, prof, geom, pup_size, params->pupil_scale_m,
                                  params->lambda_ref_m, params->lambda_s_m, z_bin) != 0)
    {
        return -1;
    }

    int nthreads = 1;
    int fftw_nth = 1;
    atmturb_wfs_pick_threads(pup_size, nbframes, &nthreads, &fftw_nth);

    atmturb_fresnel_ctx_t *ctxs = (atmturb_fresnel_ctx_t *) calloc((size_t) nthreads,
                                                                   sizeof(atmturb_fresnel_ctx_t));
    if (ctxs == NULL)
    {
        atmturb_fresnel_plan_free(&plan);
        return -1;
    }

    for (int tid = 0; tid < nthreads; tid++)
    {
        const atmturb_fresnel_ctx_t *src = (tid > 0) ? &ctxs[0] : NULL;
        if (atmturb_fresnel_ctx_init_mt(&ctxs[tid], pup_size, fftw_nth, src) != 0)
        {
            for (int j = 0; j < tid; j++)
            {
                atmturb_fresnel_ctx_free(&ctxs[j]);
            }
            free(ctxs);
            atmturb_fresnel_plan_free(&plan);
            return -1;
        }
    }

    long frame_pixels = pup_size * pup_size;

    #pragma omp parallel for num_threads(nthreads) schedule(dynamic, 2)
    for (long t = 0; t < nbframes; t++)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        long slice = t * frame_pixels;
        float *pha_slice  = &dcimg[imgs->ID_pha].array.F[slice];
        float *amp_slice  = (imgs->ID_amp >= 0) ? &dcimg[imgs->ID_amp].array.F[slice] : NULL;
        float *spha_slice = (imgs->ID_spha >= 0) ? &dcimg[imgs->ID_spha].array.F[slice] : NULL;
        float *samp_slice = (imgs->ID_samp >= 0) ? &dcimg[imgs->ID_samp].array.F[slice] : NULL;

        atmturb_fresnel_render_step(&ctxs[tid], &plan, r, geom, t, params->time_step_s,
                                    master_size, pup_size, pha_slice, amp_slice,
                                    spha_slice, samp_slice);
    }

    for (int tid = 0; tid < nthreads; tid++)
    {
        atmturb_fresnel_ctx_free(&ctxs[tid]);
    }
    free(ctxs);
    atmturb_fresnel_plan_free(&plan);

    return 0;
}

/**
 * atmturb_wfs_render_rytov - Multi-threaded Rytov diffractive rendering of simulation time steps
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @prof: Active turbulence profile.
 * @params: Observation parameters container.
 * @master_size: Master screen dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 * @nbframes: Number of simulation frames.
 * @imgs: Container of output 3D image handles.
 *
 * Return: 0 on success, -1 on plan or context allocation failure.
 */
static int atmturb_wfs_render_rytov(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    const atmturb_profile_t    *prof,
    const atmturb_obs_params_t *params,
    long                        master_size,
    long                        pup_size,
    long                        nbframes,
    const atmturb_wfs_images_t *imgs)
{
    atmturb_rytov_plan_t plan;
    double z_bin = (double) CONF_FRESNEL_PROPAGATION_BIN;
    if (atmturb_rytov_plan_init(&plan, prof, geom, pup_size, params->pupil_scale_m,
                                params->lambda_ref_m, params->lambda_s_m, z_bin,
                                CONF_FRESNEL_RYTOV_SEC_EXACT) != 0)
    {
        return -1;
    }

    int nthreads = 1;
    int fftw_nth = 1;
    atmturb_wfs_pick_threads(plan.pad_size, nbframes, &nthreads, &fftw_nth);

    atmturb_rytov_ctx_t *ctxs = (atmturb_rytov_ctx_t *) calloc((size_t) nthreads,
                                                               sizeof(atmturb_rytov_ctx_t));
    if (ctxs == NULL)
    {
        atmturb_rytov_plan_free(&plan);
        return -1;
    }

    for (int tid = 0; tid < nthreads; tid++)
    {
        const atmturb_rytov_ctx_t *src = (tid > 0) ? &ctxs[0] : NULL;
        if (atmturb_rytov_ctx_init_mt(&ctxs[tid], plan.pad_size, fftw_nth, src) != 0)
        {
            for (int j = 0; j < tid; j++)
            {
                atmturb_rytov_ctx_free(&ctxs[j]);
            }
            free(ctxs);
            atmturb_rytov_plan_free(&plan);
            return -1;
        }
    }

    long frame_pixels = pup_size * pup_size;

    #pragma omp parallel for num_threads(nthreads) schedule(dynamic, 2)
    for (long t = 0; t < nbframes; t++)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        long   slice      = t * frame_pixels;
        float *pha_slice  = &dcimg[imgs->ID_pha].array.F[slice];
        float *amp_slice  = (imgs->ID_amp >= 0) ? &dcimg[imgs->ID_amp].array.F[slice] : NULL;
        float *spha_slice = (imgs->ID_spha >= 0) ? &dcimg[imgs->ID_spha].array.F[slice] : NULL;
        float *samp_slice = (imgs->ID_samp >= 0) ? &dcimg[imgs->ID_samp].array.F[slice] : NULL;

        atmturb_rytov_render_step(&ctxs[tid], &plan, r, geom, t, params->time_step_s,
                                  master_size, pup_size, pha_slice, amp_slice,
                                  spha_slice, samp_slice);
    }

    for (int tid = 0; tid < nthreads; tid++)
    {
        atmturb_rytov_ctx_free(&ctxs[tid]);
    }
    free(ctxs);
    atmturb_rytov_plan_free(&plan);

    return 0;
}

/**
 * atmturb_wfs_render_frames - Dispatch rendering to diffractive or geometric engine
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @prof: Active turbulence profile.
 * @params: Observation parameters container.
 * @master_size: Master screen dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 * @nbframes: Number of simulation frames.
 * @imgs: Container of output 3D image handles.
 */
void atmturb_wfs_render_frames(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    const atmturb_profile_t    *prof,
    const atmturb_obs_params_t *params,
    long                        master_size,
    long                        pup_size,
    long                        nbframes,
    const atmturb_wfs_images_t *imgs)
{
#ifdef HAVE_CUDA
    int gpu_ok      = atmturb_simd_is_gpu();
    int gpu_forced  = atmturb_simd_is_gpu_forced();
    int use_gpu_dif = gpu_ok && (gpu_forced || pup_size >= 256 ||
                                 (pup_size >= 128 && nbframes >= 32));
    int use_gpu_geo = gpu_ok && (gpu_forced || pup_size >= 256 ||
                                 (pup_size >= 128 && nbframes >= 64));
#endif

    if (CONF_FRESNEL_PROPAGATION == 1)
    {
#ifdef HAVE_CUDA
        if (use_gpu_dif &&
            atmturb_wfs_render_fresnel_cuda(r, geom, prof, params, master_size,
                                            pup_size, nbframes, imgs) == 0)
        {
            return;
        }
#endif
        if (atmturb_wfs_render_diffractive(r, geom, prof, params, master_size,
                                            pup_size, nbframes, imgs) == 0)
        {
            return;
        }
    }
    else if (CONF_FRESNEL_PROPAGATION == 2)
    {
#ifdef HAVE_CUDA
        if (use_gpu_dif &&
            atmturb_wfs_render_rytov_cuda(r, geom, prof, params, master_size,
                                          pup_size, nbframes, imgs) == 0)
        {
            return;
        }
#endif
        if (atmturb_wfs_render_rytov(r, geom, prof, params, master_size,
                                      pup_size, nbframes, imgs) == 0)
        {
            return;
        }
    }

#ifdef HAVE_CUDA
    if (use_gpu_geo &&
        atmturb_wfs_render_cuda(r, geom, master_size, pup_size, nbframes,
                                params->time_step_s, imgs) == 0)
    {
        return;
    }
#endif

    atmturb_wfs_render_geometric(r, geom, master_size, pup_size, nbframes,
                                  params->time_step_s, imgs);
}
