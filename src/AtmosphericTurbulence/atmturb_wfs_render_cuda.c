// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_wfs_render_cuda.c
 * @brief   CUDA GPU dispatch wrappers for wavefront series rendering
 */

#ifdef HAVE_CUDA

#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "CommandLineInterface/CLIcore.h"
#include "atmturb_cuda.h"
#include "atmturb_cuda_rytov.h"
#include "atmturb_geometry.h"
#include "atmturb_profile.h"
#include "atmturb_rolling.h"
#include "atmturb_rytov.h"
#include "atmturb_types.h"
#include "atmturb_wfs_render.h"
#include "atmturb_wfs_render_cuda.h"

int atmturb_wfs_render_cuda(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    long                        master_size,
    long                        pup_size,
    long                        nbframes,
    const atmturb_wfs_images_t *imgs)
{
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

static inline void atmturb_wfs_cuda_build_frame_sublayers(
    const atmturb_rolling_t       *r,
    const atmturb_geom_t          *geom,
    const atmturb_rytov_plan_t    *plan,
    long                           t,
    double                         time_step_s,
    double                         offset_os,
    int                            max_sublayers,
    atmturb_cuda_rytov_sublayer_t *sublayers)
{
    for (int m = 0; m < plan->nsuper; m++)
    {
        const atmturb_superlayer_t *sl = &plan->supers[m];
        size_t base = (size_t) (m * max_sublayers);

        for (int j = 0; j < sl->nlayers; j++)
        {
            int k = sl->layer_indices[j];
            float ws = (sl->layer_weights != NULL) ? sl->layer_weights[j] : 1.0f;
            atmturb_rolling_eval_t rev;
            atmturb_rolling_get_frame(r, k, t, time_step_s, &rev);

            double dx = (geom->layers[k].traj_x != NULL)
                        ? geom->layers[k].traj_x[t]
                        : ((double) t * geom->layers[k].vx_pix);
            double dy = (geom->layers[k].traj_y != NULL)
                        ? geom->layers[k].traj_y[t]
                        : ((double) t * geom->layers[k].vy_pix);

            sublayers[base + j].k     = k;
            sublayers[base + j].w_pri = (float) (geom->layers[k].weight * (double) ws *
                                                 (double) rev.wA);
            sublayers[base + j].w_sec = (float) (geom->layers[k].weight_s * (double) ws *
                                                 (double) rev.wA);
            sublayers[base + j].x     = (float) (geom->layers[k].x0 + dx - offset_os);
            sublayers[base + j].y     = (float) (geom->layers[k].y0 + dy - offset_os);
            sublayers[base + j].xs    = (float) (geom->layers[k].xs0 + dx - offset_os);
            sublayers[base + j].ys    = (float) (geom->layers[k].ys0 + dy - offset_os);
        }
    }
}

static void atmturb_wfs_cuda_build_sublayers(
    const atmturb_rolling_t       *r,
    const atmturb_geom_t          *geom,
    const atmturb_rytov_plan_t    *plan,
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
            r, geom, plan, t, time_step_s, offset_os, max_sublayers,
            sublayers + (size_t) t * frame_stride);
    }
}

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
    if (r->lowfreq != 0 || !atmturb_cuda_device_available())
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
        if (plan.supers[m].nlayers > max_sublayers)
        {
            max_sublayers = plan.supers[m].nlayers;
        }
    }
    if (max_sublayers > 64)
    {
        atmturb_rytov_plan_free(&plan);
        return -1;
    }

    size_t n_sub_tot = (size_t) (nbframes * plan.nsuper * max_sublayers);
    atmturb_cuda_rytov_sublayer_t *sublayers =
        (atmturb_cuda_rytov_sublayer_t *) calloc(n_sub_tot, sizeof(atmturb_cuda_rytov_sublayer_t));
    int *super_nlayers   = (int *) malloc(sizeof(int) * (size_t) plan.nsuper);
    double *super_dist_m = (double *) malloc(sizeof(double) * (size_t) plan.nsuper);
    const float **masters = (const float **) malloc(sizeof(float *) * (size_t) geom->nlayers);

    if (!sublayers || !super_nlayers || !super_dist_m || !masters)
    {
        free(sublayers);
        free(super_nlayers);
        free(super_dist_m);
        free(masters);
        atmturb_rytov_plan_free(&plan);
        return -1;
    }

    for (int k = 0; k < geom->nlayers; k++)
    {
        masters[k] = r->layers[k].screens[0].data;
    }
    for (int m = 0; m < plan.nsuper; m++)
    {
        super_nlayers[m] = plan.supers[m].nlayers;
        super_dist_m[m]  = plan.supers[m].dist_m;
    }

    atmturb_wfs_cuda_build_sublayers(r, geom, &plan, nbframes, params->time_step_s,
                                     max_sublayers, sublayers);

    atmturb_cuda_rytov_params_t cparams = {
        .msize           = master_size,
        .pup_size        = pup_size,
        .pad_size        = plan.pad_size,
        .guard_pix       = plan.guard_pix,
        .nbframes        = nbframes,
        .nblayers        = geom->nlayers,
        .nsuper          = plan.nsuper,
        .has_sec         = (plan.lambda_s_m > 0.0) ? 1 : 0,
        .sec_shared      = plan.sec_shared,
        .os              = (int) geom->oversample,
        .interp          = geom->interp,
        .use_moisan      = plan.use_moisan,
        .h_masters       = masters,
        .h_laplace_inv   = plan.laplace_inv,
        .h_exp_x         = (const void *) plan.exp_x,
        .h_exp_y         = (const void *) plan.exp_y,
        .h_filt_a_pri    = (const float *const *) plan.filter_a_pri,
        .h_filt_b_pri    = (const float *const *) plan.filter_b_pri,
        .h_filt_a_sec    = (const float *const *) plan.filter_a_sec,
        .h_filt_b_sec    = (const float *const *) plan.filter_b_sec,
        .h_chrom_ramp    = (const void *const *) plan.chrom_ramp,
        .frame_sublayers = sublayers,
        .super_nlayers   = super_nlayers,
        .super_dist_m    = super_dist_m,
        .max_sublayers   = max_sublayers,
        .pha             = &dcimg[imgs->ID_pha].array.F[0],
        .amp             = &dcimg[imgs->ID_amp].array.F[0],
        .spha            = (plan.lambda_s_m > 0.0) ? &dcimg[imgs->ID_spha].array.F[0] : NULL,
        .samp            = (plan.lambda_s_m > 0.0) ? &dcimg[imgs->ID_samp].array.F[0] : NULL
    };

    int ret = atmturb_cuda_rytov_render(&cparams);

    free(sublayers);
    free(super_nlayers);
    free(super_dist_m);
    free(masters);
    atmturb_rytov_plan_free(&plan);

    return ret;
}

struct atmturb_cuda_rytov_stream
{
    atmturb_rytov_plan_t          plan;
    atmturb_cuda_rytov_params_t   cparams;
    atmturb_cuda_rytov_sublayer_t *sublayers;
    int                          *super_nlayers;
    double                       *super_dist_m;
    const float                 **masters;
    int                           max_sublayers;
    double                        offset_os;
};

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
    if (r->lowfreq != 0 || !atmturb_cuda_device_available())
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
        if (ctx->plan.supers[m].nlayers > max_sublayers)
        {
            max_sublayers = ctx->plan.supers[m].nlayers;
        }
    }
    if (max_sublayers > 64)
    {
        atmturb_rytov_plan_free(&ctx->plan);
        free(ctx);
        return NULL;
    }

    ctx->max_sublayers = max_sublayers;
    ctx->offset_os     = (double) ctx->plan.guard_pix * (double) geom->oversample;
    size_t n_sub_tot   = (size_t) (ctx->plan.nsuper * max_sublayers);

    ctx->sublayers     = (atmturb_cuda_rytov_sublayer_t *) calloc(
        n_sub_tot, sizeof(atmturb_cuda_rytov_sublayer_t));
    ctx->super_nlayers = (int *) malloc(sizeof(int) * (size_t) ctx->plan.nsuper);
    ctx->super_dist_m  = (double *) malloc(sizeof(double) * (size_t) ctx->plan.nsuper);
    ctx->masters       = (const float **) malloc(sizeof(float *) * (size_t) geom->nlayers);

    if (!ctx->sublayers || !ctx->super_nlayers || !ctx->super_dist_m || !ctx->masters)
    {
        atmturb_cuda_rytov_stream_free(ctx);
        return NULL;
    }

    for (int k = 0; k < geom->nlayers; k++)
    {
        ctx->masters[k] = r->layers[k].screens[0].data;
    }
    for (int m = 0; m < ctx->plan.nsuper; m++)
    {
        ctx->super_nlayers[m] = ctx->plan.supers[m].nlayers;
        ctx->super_dist_m[m]  = ctx->plan.supers[m].dist_m;
    }

    ctx->cparams.msize           = master_size;
    ctx->cparams.pup_size        = pup_size;
    ctx->cparams.pad_size        = ctx->plan.pad_size;
    ctx->cparams.guard_pix       = ctx->plan.guard_pix;
    ctx->cparams.nbframes        = 1;
    ctx->cparams.nblayers        = geom->nlayers;
    ctx->cparams.nsuper          = ctx->plan.nsuper;
    ctx->cparams.has_sec         = (ctx->plan.lambda_s_m > 0.0) ? 1 : 0;
    ctx->cparams.sec_shared      = ctx->plan.sec_shared;
    ctx->cparams.os              = (int) geom->oversample;
    ctx->cparams.interp          = geom->interp;
    ctx->cparams.use_moisan      = ctx->plan.use_moisan;
    ctx->cparams.h_masters       = ctx->masters;
    ctx->cparams.h_laplace_inv   = ctx->plan.laplace_inv;
    ctx->cparams.h_exp_x         = (const void *) ctx->plan.exp_x;
    ctx->cparams.h_exp_y         = (const void *) ctx->plan.exp_y;
    ctx->cparams.h_filt_a_pri    = (const float *const *) ctx->plan.filter_a_pri;
    ctx->cparams.h_filt_b_pri    = (const float *const *) ctx->plan.filter_b_pri;
    ctx->cparams.h_filt_a_sec    = (const float *const *) ctx->plan.filter_a_sec;
    ctx->cparams.h_filt_b_sec    = (const float *const *) ctx->plan.filter_b_sec;
    ctx->cparams.h_chrom_ramp    = (const void *const *) ctx->plan.chrom_ramp;
    ctx->cparams.frame_sublayers = ctx->sublayers;
    ctx->cparams.super_nlayers   = ctx->super_nlayers;
    ctx->cparams.super_dist_m    = ctx->super_dist_m;
    ctx->cparams.max_sublayers   = max_sublayers;

    return ctx;
}

/**
 * atmturb_cuda_rytov_stream_render_step - Render one stream frame on GPU
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
    if (!ctx)
    {
        return -1;
    }

    atmturb_wfs_cuda_build_frame_sublayers(
        r, geom, &ctx->plan, t, time_step_s, ctx->offset_os,
        ctx->max_sublayers, ctx->sublayers);

    ctx->cparams.pha  = pha;
    ctx->cparams.amp  = amp;
    ctx->cparams.spha = spha;
    ctx->cparams.samp = samp;

    return atmturb_cuda_rytov_render(&ctx->cparams);
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
    if (ctx->sublayers)
    {
        free(ctx->sublayers);
    }
    if (ctx->super_nlayers)
    {
        free(ctx->super_nlayers);
    }
    if (ctx->super_dist_m)
    {
        free(ctx->super_dist_m);
    }
    if (ctx->masters)
    {
        free(ctx->masters);
    }
    atmturb_rytov_plan_free(&ctx->plan);
    atmturb_cuda_rytov_cleanup();
    free(ctx);
}

#endif // HAVE_CUDA
