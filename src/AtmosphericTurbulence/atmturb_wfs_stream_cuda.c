// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_wfs_stream_cuda.c
 * @brief   GPU execution contexts for real-time wavefront streaming
 */

#ifdef HAVE_CUDA

#include <stdlib.h>
#include <string.h>

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
 * struct atmturb_cuda_rytov_stream - Persistent context for Rytov streaming
 * @plan: Rytov perturbation propagation plan.
 * @cparams: Configuration parameters for GPU Rytov synthesis engine.
 * @sublayers: Double-buffered host sublayer parameter arrays.
 * @super_nlayers: Number of physical layers per superlayer.
 * @super_dist_m: Propagation distances per superlayer in meters.
 * @lf_modes: Double-buffered low-frequency mode arrays.
 * @super_nlf_modes: Number of low-frequency modes per superlayer.
 * @masters: Host pointers to master screen arrays.
 * @screen_offsets: Screen index offsets per layer.
 * @max_sublayers: Maximum sublayers across all superlayers.
 * @max_lf_modes: Maximum low-frequency modes across all superlayers.
 * @offset_os: Offset in oversampled master pixels for guard band alignment.
 * @ev_done: CUDA completion events for double-buffer slots.
 */
struct atmturb_cuda_rytov_stream
{
    atmturb_rytov_plan_t          plan;
    atmturb_cuda_rytov_params_t   cparams;
    atmturb_cuda_rytov_sublayer_t *sublayers[2];
    int                          *super_nlayers;
    double                       *super_dist_m;
    atmturb_cuda_lowfreq_mode_t  *lf_modes[2];
    int                          *super_nlf_modes;
    const float                 **masters;
    int                          *screen_offsets;
    int                           max_sublayers;
    int                           max_lf_modes;
    double                        offset_os;
    void                         *ev_done[2];
};

/**
 * atmturb_wfs_cuda_init_rytov_lf - Allocate and initialize subharmonic mode arrays
 * @r: Rolling simulation context.
 * @plan: Rytov perturbation plan.
 * @nbframes: Number of frames to allocate.
 * @out_max_lf: Output pointer for maximum low-frequency mode count.
 * @out_super_nlf: Output pointer for mode counts per superlayer.
 * @out_lf_modes: Output pointer for pinned host mode buffer.
 *
 * Return: 0 on success, -1 on failure.
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
 * atmturb_cuda_rytov_stream_free - Free persistent GPU Rytov streaming context
 * @ctx: Stream context to release.
 */
void atmturb_cuda_rytov_stream_free(
    atmturb_cuda_rytov_stream_t *ctx)
{
    if (!ctx)
    {
        return;
    }

    for (int s = 0; s < 2; s++)
    {
        if (ctx->ev_done[s] != NULL)
        {
            atmturb_cuda_rytov_destroy_event(ctx->ev_done[s]);
            ctx->ev_done[s] = NULL;
        }
        free(ctx->sublayers[s]);
        ctx->sublayers[s] = NULL;
        if (ctx->lf_modes[s] != NULL)
        {
            atmturb_cuda_host_free(ctx->lf_modes[s]);
            ctx->lf_modes[s] = NULL;
        }
    }
    free(ctx->super_nlayers);
    free(ctx->super_dist_m);
    free(ctx->super_nlf_modes);
    free(ctx->masters);
    free(ctx->screen_offsets);

    atmturb_rytov_plan_free(&ctx->plan);
    free(ctx);
}

/**
 * atmturb_cuda_rytov_stream_init - Initialize persistent GPU Rytov streaming context
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @prof: Active turbulence profile.
 * @params: Observation parameters container.
 * @master_size: Master screen linear dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 *
 * Return: Allocated streaming context, or NULL on failure.
 */
atmturb_cuda_rytov_stream_t *atmturb_cuda_rytov_stream_init(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    const atmturb_profile_t    *prof,
    const atmturb_obs_params_t *params,
    long                        master_size,
    long                        pup_size)
{
    if (!atmturb_cuda_device_available())
    {
        return NULL;
    }

    atmturb_cuda_rytov_stream_t *ctx =
        (atmturb_cuda_rytov_stream_t *) calloc(1, sizeof(atmturb_cuda_rytov_stream_t));
    if (!ctx)
    {
        return NULL;
    }

    double z_bin = (double) CONF_FRESNEL_PROPAGATION_BIN;
    if (atmturb_rytov_plan_init(&ctx->plan, prof, geom, pup_size, params->pupil_scale_m,
                                params->lambda_ref_m, params->lambda_s_m, z_bin,
                                CONF_FRESNEL_RYTOV_SEC_EXACT) != 0)
    {
        free(ctx);
        return NULL;
    }

    int max_sublayers = 1;
    for (int m = 0; m < ctx->plan.nsuper; m++)
    {
        int n_sub = 2 * ctx->plan.supers[m].nlayers;
        if (n_sub > max_sublayers)
        {
            max_sublayers = n_sub;
        }
    }
    if (max_sublayers > 128)
    {
        atmturb_rytov_plan_free(&ctx->plan);
        free(ctx);
        return NULL;
    }

    int total_screens = 0;
    if (atmturb_wfs_cuda_init_screen_pool(r, geom->nlayers, &ctx->masters,
                                         &ctx->screen_offsets, &total_screens) != 0)
    {
        atmturb_cuda_rytov_stream_free(ctx);
        return NULL;
    }

    ctx->max_sublayers = max_sublayers;
    ctx->offset_os     = (double) ctx->plan.guard_pix * (double) geom->oversample;
    size_t n_sub_tot   = (size_t) (ctx->plan.nsuper * max_sublayers);

    for (int s = 0; s < 2; s++)
    {
        ctx->sublayers[s] = (atmturb_cuda_rytov_sublayer_t *) calloc(
            n_sub_tot, sizeof(atmturb_cuda_rytov_sublayer_t));
        if (!ctx->sublayers[s])
        {
            atmturb_cuda_rytov_stream_free(ctx);
            return NULL;
        }
    }
    ctx->super_nlayers = (int *) malloc(sizeof(int) * (size_t) ctx->plan.nsuper);
    ctx->super_dist_m  = (double *) malloc(sizeof(double) * (size_t) ctx->plan.nsuper);

    if (atmturb_wfs_cuda_init_rytov_lf(r, &ctx->plan, 1, &ctx->max_lf_modes,
                                       &ctx->super_nlf_modes, &ctx->lf_modes[0]) != 0)
    {
        atmturb_cuda_rytov_stream_free(ctx);
        return NULL;
    }
    if (ctx->lf_modes[0] != NULL)
    {
        size_t lf_bytes = (size_t) (ctx->plan.nsuper * ctx->max_lf_modes) *
                          sizeof(atmturb_cuda_lowfreq_mode_t);
        atmturb_cuda_host_alloc((void **) &ctx->lf_modes[1], lf_bytes);
    }

    if (!ctx->super_nlayers || !ctx->super_dist_m)
    {
        atmturb_cuda_rytov_stream_free(ctx);
        return NULL;
    }

    for (int m = 0; m < ctx->plan.nsuper; m++)
    {
        ctx->super_nlayers[m] = 2 * ctx->plan.supers[m].nlayers;
        ctx->super_dist_m[m]  = ctx->plan.supers[m].dist_m;
    }

    ctx->cparams.msize           = master_size;
    ctx->cparams.pup_size        = pup_size;
    ctx->cparams.pad_size        = ctx->plan.pad_size;
    ctx->cparams.guard_pix       = ctx->plan.guard_pix;
    ctx->cparams.nbframes        = 1;
    ctx->cparams.nblayers        = total_screens;
    ctx->cparams.nsuper          = ctx->plan.nsuper;
    ctx->cparams.has_sec         = (ctx->plan.lambda_s_m > 0.0) ? 1 : 0;
    ctx->cparams.sec_shared      = ctx->plan.sec_shared;
    ctx->cparams.os              = (int) geom->oversample;
    ctx->cparams.interp          = geom->interp;
    ctx->cparams.use_moisan      = ctx->plan.use_moisan;
    ctx->cparams.lowfreq         = r->lowfreq;
    ctx->cparams.h_masters       = ctx->masters;
    ctx->cparams.h_laplace_inv   = ctx->plan.laplace_inv;
    ctx->cparams.h_exp_x         = (const void *) ctx->plan.exp_x;
    ctx->cparams.h_exp_y         = (const void *) ctx->plan.exp_y;
    ctx->cparams.h_filt_a_pri    = (const float *const *) ctx->plan.filter_a_pri;
    ctx->cparams.h_filt_b_pri    = (const float *const *) ctx->plan.filter_b_pri;
    ctx->cparams.h_filt_a_sec    = (const float *const *) ctx->plan.filter_a_sec;
    ctx->cparams.h_filt_b_sec    = (const float *const *) ctx->plan.filter_b_sec;
    ctx->cparams.h_chrom_ramp    = (const void *const *) ctx->plan.chrom_ramp;
    ctx->cparams.frame_sublayers = ctx->sublayers[0];
    ctx->cparams.super_nlayers   = ctx->super_nlayers;
    ctx->cparams.super_dist_m    = ctx->super_dist_m;
    ctx->cparams.max_sublayers   = max_sublayers;
    ctx->cparams.frame_lf_modes  = ctx->lf_modes[0];
    ctx->cparams.super_nlf_modes = ctx->super_nlf_modes;
    ctx->cparams.max_lf_modes    = ctx->max_lf_modes;
    ctx->cparams.remove_piston   = 1;

    for (int s = 0; s < 2; s++)
    {
        if (atmturb_cuda_rytov_create_event(&ctx->ev_done[s]) != 0)
        {
            atmturb_cuda_rytov_stream_free(ctx);
            return NULL;
        }
    }

    return ctx;
}

/**
 * atmturb_cuda_rytov_stream_render_launch - Enqueue GPU ops for Rytov frame asynchronously
 * @ctx: Persistent stream context.
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @t: Frame index.
 * @time_step_s: Time step in seconds.
 * @pha: Destination primary phase frame buffer.
 * @amp: Destination primary amplitude frame buffer.
 * @spha: Destination secondary phase frame buffer (or NULL).
 * @samp: Destination secondary amplitude frame buffer (or NULL).
 * @slot: Active double-buffer slot index (0 or 1).
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_cuda_rytov_stream_render_launch(
    atmturb_cuda_rytov_stream_t *ctx,
    const atmturb_rolling_t     *r,
    const atmturb_geom_t        *geom,
    long                         t,
    double                       time_step_s,
    float                       *pha,
    float                       *amp,
    float                       *spha,
    float                       *samp,
    int                          slot)
{
    if (!ctx)
    {
        return -1;
    }
    slot &= 1;

    atmturb_wfs_cuda_build_frame_sublayers(
        r, geom, ctx->plan.nsuper, ctx->plan.supers, ctx->screen_offsets,
        t, time_step_s, ctx->offset_os, ctx->max_sublayers, ctx->sublayers[slot]);

    if (r->lowfreq && ctx->lf_modes[slot] != NULL)
    {
        atmturb_wfs_cuda_build_frame_lf_modes(
            r, geom, ctx->plan.nsuper, ctx->plan.supers,
            t, time_step_s, ctx->offset_os, ctx->max_lf_modes, ctx->lf_modes[slot]);
    }

    ctx->cparams.frame_sublayers = ctx->sublayers[slot];
    ctx->cparams.frame_lf_modes  = ctx->lf_modes[slot];
    ctx->cparams.pha             = pha;
    ctx->cparams.amp             = amp;
    ctx->cparams.spha            = spha;
    ctx->cparams.samp            = samp;
    ctx->cparams.async_stream    = 1;
    ctx->cparams.slot            = slot;

    int ret = atmturb_cuda_rytov_render(&ctx->cparams);
    if (ret != 0)
    {
        return ret;
    }

    return atmturb_cuda_rytov_record_event(ctx->ev_done[slot]);
}

/**
 * atmturb_cuda_rytov_stream_render_wait - Wait for enqueued Rytov frame to finish on GPU
 * @ctx: Persistent stream context.
 * @slot: Buffer slot to synchronize (0 or 1).
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_cuda_rytov_stream_render_wait(
    atmturb_cuda_rytov_stream_t *ctx,
    int                          slot)
{
    if (!ctx)
    {
        return -1;
    }
    return atmturb_cuda_rytov_wait_event(ctx->ev_done[slot & 1]);
}

/**
 * atmturb_cuda_rytov_stream_render_step - Render one stream frame synchronously on GPU
 * @ctx: Persistent stream context.
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @t: Frame index.
 * @time_step_s: Time step in seconds.
 * @pha: Destination primary phase frame buffer.
 * @amp: Destination primary amplitude frame buffer.
 * @spha: Destination secondary phase frame buffer (or NULL).
 * @samp: Destination secondary amplitude frame buffer (or NULL).
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_cuda_rytov_stream_render_step(
    atmturb_cuda_rytov_stream_t *ctx,
    const atmturb_rolling_t     *r,
    const atmturb_geom_t        *geom,
    long                         t,
    double                       time_step_s,
    float                       *pha,
    float                       *amp,
    float                       *spha,
    float                       *samp)
{
    int ret = atmturb_cuda_rytov_stream_render_launch(
        ctx, r, geom, t, time_step_s, pha, amp, spha, samp, 0);
    if (ret != 0)
    {
        return ret;
    }
    return atmturb_cuda_rytov_stream_render_wait(ctx, 0);
}

/**
 * atmturb_wfs_cuda_geom_stream_init - Initialize persistent GPU geometric streaming context
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @master_size: Master screen linear dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 *
 * Return: Allocated stream context, or NULL on failure.
 */
atmturb_cuda_geom_stream_t *atmturb_wfs_cuda_geom_stream_init(
    const atmturb_rolling_t *r,
    const atmturb_geom_t    *geom,
    long                     master_size,
    long                     pup_size)
{
    if (r == NULL || geom == NULL || !atmturb_cuda_device_available() ||
        r->lowfreq != 0 || r->rolling != 0 ||
        geom->oversample != 1 || geom->interp != 0)
    {
        return NULL;
    }

    double w0 = geom->layers[0].weight;
    double scoeff = (w0 > 0.0) ? (geom->layers[0].weight_s / w0) : 1.0;

    for (int k = 0; k < geom->nlayers; k++)
    {
        if (geom->layers[k].traj_x != NULL || geom->layers[k].traj_y != NULL ||
            geom->layers[k].xs0 != geom->layers[k].x0 ||
            geom->layers[k].ys0 != geom->layers[k].y0)
        {
            return NULL;
        }
        double wk = geom->layers[k].weight;
        if (wk > 0.0 && fabs(geom->layers[k].weight_s / wk - scoeff) > 1e-4)
        {
            return NULL;
        }
    }

    int nl = geom->nlayers;
    size_t m_sz = (size_t) nl * sizeof(float *);
    size_t d_sz = (size_t) nl * sizeof(double);
    void *buf = malloc(m_sz + 5 * d_sz);
    if (!buf)
    {
        return NULL;
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

    atmturb_cuda_sim_params_t cparams = {
        .nblayers  = nl,
        .msize     = master_size,
        .pup_size  = pup_size,
        .nbframes  = 1,
        .Scoeff    = scoeff,
        .h_masters = (const float *const *) masters,
        .vxpix     = vx,
        .vypix     = vy,
        .cn2       = cn2,
        .x0        = x0,
        .y0        = y0
    };

    atmturb_cuda_geom_stream_t *ctx = atmturb_cuda_geom_stream_init(&cparams);
    free(buf);
    return ctx;
}

/**
 * atmturb_wfs_cuda_geom_stream_render_step - Render one geometric stream frame on GPU
 * @ctx: Persistent stream context.
 * @t: Frame index.
 * @pha: Destination primary phase frame buffer.
 * @amp: Destination primary amplitude frame buffer (or NULL).
 * @spha: Destination secondary phase frame buffer (or NULL).
 * @samp: Destination secondary amplitude frame buffer (or NULL).
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_wfs_cuda_geom_stream_render_step(
    atmturb_cuda_geom_stream_t *ctx,
    long                        t,
    float                      *pha,
    float                      *amp,
    float                      *spha,
    float                      *samp)
{
    return atmturb_cuda_geom_stream_render_step(ctx, t, pha, amp, spha, samp);
}

/**
 * atmturb_wfs_cuda_geom_stream_free - Free persistent GPU geometric streaming context
 * @ctx: Stream context to release.
 */
void atmturb_wfs_cuda_geom_stream_free(
    atmturb_cuda_geom_stream_t *ctx)
{
    atmturb_cuda_geom_stream_free(ctx);
}

#endif // HAVE_CUDA
