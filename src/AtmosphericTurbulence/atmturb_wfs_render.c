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
 * atmturb_wfs_extrude_channel - Extrude a single turbulence layer phase into pupil
 * @scr: Master screen pointer.
 * @master_size: Master screen linear dimension.
 * @x: Unwrapped X position in master pixels.
 * @y: Unwrapped Y position in master pixels.
 * @pup_size: Linear dimension of pupil.
 * @geom: Observing geometry context.
 * @weight: Scaled layer weight.
 * @out_pha: Destination phase array to accumulate into.
 */
static inline void atmturb_wfs_extrude_channel(
    const float          *scr,
    long                  master_size,
    double                x,
    double                y,
    long                  pup_size,
    const atmturb_geom_t *geom,
    float                 weight,
    float                *out_pha)
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

    atmturb_extrude_params_t ep;
    ep.master   = scr;
    ep.msize    = master_size;
    ep.x0       = cur_x;
    ep.y0       = cur_y;
    ep.pup_size = pup_size;
    ep.os       = (long) geom->oversample;
    ep.interp   = geom->interp;
    ep.weight   = weight;
    ep.out_pha  = out_pha;
    atmturb_extrude_accumulate(&ep);
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
    long pup_size = target->pup_size;
    long guard = target->guard_pix;
    long pad_size = pup_size + 2 * guard;
    double offset_os = (double) guard * (double) geom->oversample;

    const atmturb_layer_geom_t *lg = &geom->layers[k];
    atmturb_rolling_eval_t rev;
    atmturb_rolling_get_frame(r, k, t, time_step_s, &rev);

    double dx = (lg->traj_x != NULL) ? lg->traj_x[t] : ((double) t * lg->vx_pix);
    double dy = (lg->traj_y != NULL) ? lg->traj_y[t] : ((double) t * lg->vy_pix);

    double wscale = (target->weight_scale > 0.0f) ? (double) target->weight_scale : 1.0;
    double w_pri  = lg->weight * wscale;
    double w_sec  = lg->weight_s * wscale;

    double x = lg->x0 + dx - offset_os;
    double y = lg->y0 + dy - offset_os;
    atmturb_wfs_extrude_channel(rev.scrA, master_size, x, y, pad_size, geom,
                                (float) (w_pri * (double) rev.wA), target->pha);
    if (rev.wB > 0.0f && rev.scrB != NULL)
    {
        atmturb_wfs_extrude_channel(rev.scrB, master_size, x, y, pad_size, geom,
                                    (float) (w_pri * (double) rev.wB), target->pha);
    }

    double xs = lg->xs0 + dx - offset_os;
    double ys = lg->ys0 + dy - offset_os;
    atmturb_wfs_extrude_channel(rev.scrA, master_size, xs, ys, pad_size, geom,
                                (float) (w_sec * (double) rev.wA), target->spha);
    if (rev.wB > 0.0f && rev.scrB != NULL)
    {
        atmturb_wfs_extrude_channel(rev.scrB, master_size, xs, ys, pad_size, geom,
                                    (float) (w_sec * (double) rev.wB), target->spha);
    }

    if (r->lowfreq)
    {
        atmturb_lowfreq_accumulate_custom(&r->layers[k].lf_base, rev.are_eff, rev.aim_eff,
                                          x, y, pad_size, (long) geom->oversample,
                                          (float) w_pri, target->pha);

        atmturb_lowfreq_accumulate_custom(&r->layers[k].lf_base, rev.are_eff, rev.aim_eff,
                                          xs, ys, pad_size, (long) geom->oversample,
                                          (float) w_sec, target->spha);
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

    #pragma omp parallel for schedule(dynamic)
    for (long t = 0; t < nbframes; t++)
    {
        long slice = t * frame_pixels;
        float *pha_slice  = &dcimg[imgs->ID_pha].array.F[slice];
        float *amp_slice  = &dcimg[imgs->ID_amp].array.F[slice];
        float *spha_slice = &dcimg[imgs->ID_spha].array.F[slice];
        float *samp_slice = &dcimg[imgs->ID_samp].array.F[slice];

        atmturb_init_phase_amp(pha_slice, amp_slice, frame_pixels);
        atmturb_init_phase_amp(spha_slice, samp_slice, frame_pixels);

        for (int k = 0; k < geom->nlayers; k++)
        {
            atmturb_wfs_render_layer(r, geom, k, t, time_step_s, master_size,
                                     pup_size, pha_slice, spha_slice);
        }
    }
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
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    atmturb_fresnel_ctx_t *ctxs = (atmturb_fresnel_ctx_t *) calloc((size_t) nthreads,
                                                                   sizeof(atmturb_fresnel_ctx_t));
    if (ctxs == NULL)
    {
        atmturb_fresnel_plan_free(&plan);
        return -1;
    }

    for (int tid = 0; tid < nthreads; tid++)
    {
        if (atmturb_fresnel_ctx_init(&ctxs[tid], pup_size) != 0)
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

    #pragma omp parallel for schedule(dynamic)
    for (long t = 0; t < nbframes; t++)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        long slice = t * frame_pixels;
        float *pha_slice  = &dcimg[imgs->ID_pha].array.F[slice];
        float *amp_slice  = &dcimg[imgs->ID_amp].array.F[slice];
        float *spha_slice = &dcimg[imgs->ID_spha].array.F[slice];
        float *samp_slice = &dcimg[imgs->ID_samp].array.F[slice];

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
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    atmturb_rytov_ctx_t *ctxs = (atmturb_rytov_ctx_t *) calloc((size_t) nthreads,
                                                               sizeof(atmturb_rytov_ctx_t));
    if (ctxs == NULL)
    {
        atmturb_rytov_plan_free(&plan);
        return -1;
    }

    for (int tid = 0; tid < nthreads; tid++)
    {
        if (atmturb_rytov_ctx_init(&ctxs[tid], plan.pad_size) != 0)
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

    #pragma omp parallel for schedule(dynamic)
    for (long t = 0; t < nbframes; t++)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        long slice = t * frame_pixels;
        float *pha_slice  = &dcimg[imgs->ID_pha].array.F[slice];
        float *amp_slice  = &dcimg[imgs->ID_amp].array.F[slice];
        float *spha_slice = &dcimg[imgs->ID_spha].array.F[slice];
        float *samp_slice = &dcimg[imgs->ID_samp].array.F[slice];

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
    if (CONF_FRESNEL_PROPAGATION == 1)
    {
        if (atmturb_wfs_render_diffractive(r, geom, prof, params, master_size,
                                           pup_size, nbframes, imgs) == 0)
        {
            return;
        }
    }
    else if (CONF_FRESNEL_PROPAGATION == 2)
    {
#ifdef HAVE_CUDA
        if (atmturb_simd_is_gpu())
        {
            if (atmturb_wfs_render_rytov_cuda(r, geom, prof, params, master_size,
                                              pup_size, nbframes, imgs) == 0)
            {
                return;
            }
        }
#endif
        if (atmturb_wfs_render_rytov(r, geom, prof, params, master_size,
                                     pup_size, nbframes, imgs) == 0)
        {
            return;
        }
    }

#ifdef HAVE_CUDA
    if (atmturb_simd_is_gpu())
    {
        if (atmturb_wfs_render_cuda(r, geom, master_size, pup_size, nbframes, imgs) == 0)
        {
            return;
        }
    }
#endif

    atmturb_wfs_render_geometric(r, geom, master_size, pup_size, nbframes,
                                 params->time_step_s, imgs);
}
