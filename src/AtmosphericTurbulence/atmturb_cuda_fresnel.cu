// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_fresnel.cu
 * @brief   CUDA GPU and cuFFT accelerated multi-layer Fresnel diffractive propagation
 */

#include <cuda_runtime.h>
#include <cufft.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include "atmturb_cuda.h"
#include "atmturb_cuda_fresnel.h"
#include "atmturb_cuda_rytov_internal.h"
#include "atmturb_fresnel.h"
#include "atmturb_geometry.h"
#include "atmturb_lowfreq.h"
#include "atmturb_profile.h"
#include "atmturb_rolling.h"
#include "atmturb_superlayer.h"

extern float CONF_FRESNEL_PROPAGATION_BIN;
}

struct atmturb_cuda_fresnel_stream
{
    atmturb_fresnel_plan_t         plan;
    long                           msize;
    long                           pup_size;
    long                           npix;
    int                            has_sec;
    int                            interp;
    int                            os;
    int                            total_screens;
    int                            max_sublayers;
    int                            lowfreq;
    int                            max_lf_modes;
    int                           *super_nlf_modes;
    const float                  **masters;
    int                           *screen_offsets;
    atmturb_cuda_rytov_sublayer_t *h_sublayers;
    atmturb_cuda_rytov_sublayer_t *d_sublayers;
    atmturb_cuda_lowfreq_mode_t   *h_lf_modes;
    atmturb_cuda_lowfreq_mode_t   *d_lf_modes;
    const float                   *d_masters;

    cudaStream_t                   stream;
    cufftHandle                    plan_c2c;

    cufftComplex                 **d_tf_pri;
    cufftComplex                 **d_tf_sec;

    cufftComplex                  *d_field_pri;
    cufftComplex                  *d_field_sec;
    float                         *d_super_pha;
    float                         *d_super_spha;
    float                         *d_frame_pha;
    float                         *d_frame_amp;
    float                         *d_frame_spha;
    float                         *d_frame_samp;
};

__global__ static void atmturb_cuda_fresnel_init_modulate_kernel(
    cufftComplex *d_field,
    const float  *d_super_pha,
    int           npix)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= npix)
    {
        return;
    }
    float p = d_super_pha[i];
    float s, c;
    __sincosf(p, &s, &c);
    d_field[i] = make_cuFloatComplex(c, s);
}

__global__ static void atmturb_cuda_fresnel_modulate_kernel(
    cufftComplex *d_field,
    const float  *d_super_pha,
    int           npix)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= npix)
    {
        return;
    }
    float p = d_super_pha[i];
    float s, c;
    __sincosf(p, &s, &c);
    cufftComplex f = d_field[i];
    d_field[i] = make_cuFloatComplex(f.x * c - f.y * s, f.y * c + f.x * s);
}

__global__ static void atmturb_cuda_fresnel_tf_mul_kernel(
    cufftComplex       *d_field,
    const cufftComplex *d_tf,
    int                 npix)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= npix)
    {
        return;
    }
    cufftComplex f = d_field[i];
    cufftComplex t = d_tf[i];
    d_field[i] = make_cuFloatComplex(f.x * t.x - f.y * t.y, f.y * t.x + f.x * t.y);
}

__global__ static void atmturb_cuda_fresnel_extract_kernel(
    const cufftComplex *d_field,
    float              *d_pha,
    float              *d_amp,
    int                 npix)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= npix)
    {
        return;
    }
    cufftComplex f = d_field[i];
    float re = f.x;
    float im = f.y;
    if (d_amp != NULL)
    {
        d_amp[i] = hypotf(re, im);
    }
    float phi_geo = d_pha[i];
    float s_geo, c_geo;
    __sincosf(phi_geo, &s_geo, &c_geo);
    float d_re = re * c_geo + im * s_geo;
    float d_im = im * c_geo - re * s_geo;
    d_pha[i]   = phi_geo + atan2f(d_im, d_re);
}

static inline void atmturb_fresnel_cuda_build_sublayers(
    const atmturb_rolling_t       *r,
    const atmturb_geom_t          *geom,
    const atmturb_fresnel_plan_t  *plan,
    const int                     *screen_offsets,
    long                           t,
    double                         time_step_s,
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

            double dx = atmturb_geom_get_layer_dx(&geom->layers[k], t);
            double dy = atmturb_geom_get_layer_dy(&geom->layers[k], t);

            float cur_x  = (float) (geom->layers[k].x0 + dx);
            float cur_y  = (float) (geom->layers[k].y0 + dy);
            float cur_xs = (float) (geom->layers[k].xs0 + dx);
            float cur_ys = (float) (geom->layers[k].ys0 + dy);

            double w_base   = geom->layers[k].weight * (double) ws;
            double ws_base  = geom->layers[k].weight_s * (double) ws;
            int    k_offset = screen_offsets ? screen_offsets[k] : k;

            size_t idx_a = base + (size_t) (2 * j);
            sublayers[idx_a].k     = k_offset + rev.idxA;
            sublayers[idx_a].w_pri = (float) (w_base * (double) rev.wA);
            sublayers[idx_a].w_sec = (float) (ws_base * (double) rev.wA);
            sublayers[idx_a].x     = cur_x;
            sublayers[idx_a].y     = cur_y;
            sublayers[idx_a].xs    = cur_xs;
            sublayers[idx_a].ys    = cur_ys;

            size_t idx_b = base + (size_t) (2 * j + 1);
            if (rev.wB > 0.0f && rev.idxB >= 0)
            {
                sublayers[idx_b].k     = k_offset + rev.idxB;
                sublayers[idx_b].w_pri = (float) (w_base * (double) rev.wB);
                sublayers[idx_b].w_sec = (float) (ws_base * (double) rev.wB);
                sublayers[idx_b].x     = cur_x;
                sublayers[idx_b].y     = cur_y;
                sublayers[idx_b].xs    = cur_xs;
                sublayers[idx_b].ys    = cur_ys;
            }
            else
            {
                memset(&sublayers[idx_b], 0, sizeof(atmturb_cuda_rytov_sublayer_t));
            }
        }
    }
}

static inline void atmturb_fresnel_cuda_build_lf_modes(
    const atmturb_rolling_t      *r,
    const atmturb_geom_t         *geom,
    const atmturb_fresnel_plan_t *plan,
    long                          t,
    double                        time_step_s,
    int                           max_lf_modes,
    atmturb_cuda_lowfreq_mode_t  *lf_modes)
{
    if (r->lowfreq == 0 || lf_modes == NULL)
    {
        return;
    }

    for (int m = 0; m < plan->nsuper; m++)
    {
        const atmturb_superlayer_t *sl = &plan->supers[m];
        size_t base = (size_t) (m * max_lf_modes);
        int mode_idx = 0;

        for (int j = 0; j < sl->nlayers; j++)
        {
            int k = sl->layer_indices[j];
            float ws = (sl->layer_weights != NULL) ? sl->layer_weights[j] : 1.0f;
            atmturb_rolling_eval_t rev;
            atmturb_rolling_get_frame(r, k, t, time_step_s, &rev);

            const atmturb_layer_geom_t *lg = &geom->layers[k];
            double dx = atmturb_geom_get_layer_dx(lg, t);
            double dy = atmturb_geom_get_layer_dy(lg, t);

            double x  = lg->x0 + dx;
            double y  = lg->y0 + dy;
            double xs = lg->xs0 + dx;
            double ys = lg->ys0 + dy;

            double w_pri = lg->weight * (double) ws;
            double w_sec = lg->weight_s * (double) ws;
            const atmturb_lowfreq_t *lf_base = &r->layers[k].lf_base;

            for (int lm = 0; lm < lf_base->nmodes; lm++)
            {
                float kx = lf_base->kx[lm];
                float ky = lf_base->ky[lm];
                float are = rev.are_eff[lm];
                float aim = rev.aim_eff[lm];

                atmturb_cuda_lowfreq_mode_t *out = &lf_modes[base + mode_idx++];
                out->are_pri = (float) (w_pri * (double) are);
                out->aim_pri = (float) (w_pri * (double) aim);
                out->are_sec = (float) (w_sec * (double) are);
                out->aim_sec = (float) (w_sec * (double) aim);
                out->th0_pri = (float) ((double) kx * x + (double) ky * y);
                out->th0_sec = (float) ((double) kx * xs + (double) ky * ys);
                out->dkx     = (float) ((double) kx * (double) geom->oversample);
                out->dky     = (float) ((double) ky * (double) geom->oversample);
            }
        }
    }
}

static int atmturb_fresnel_cuda_init_screen_pool(
    const atmturb_rolling_t  *r,
    int                       nlayers,
    const float            ***out_masters,
    int                     **out_offsets,
    int                      *out_total)
{
    int total = 0;
    for (int k = 0; k < nlayers; k++)
    {
        total += r->layers[k].nscreens;
    }
    const float **masters = (const float **) malloc(sizeof(float *) * (size_t) total);
    int *offsets = (int *) malloc(sizeof(int) * (size_t) nlayers);
    if (!masters || !offsets)
    {
        free(masters);
        free(offsets);
        return -1;
    }
    int off = 0;
    for (int k = 0; k < nlayers; k++)
    {
        offsets[k] = off;
        for (int s = 0; s < r->layers[k].nscreens; s++)
        {
            masters[off + s] = r->layers[k].screens[s].data;
        }
        off += r->layers[k].nscreens;
    }
    *out_masters = masters;
    *out_offsets = offsets;
    *out_total   = total;
    return 0;
}

/**
 * atmturb_cuda_fresnel_stream_free - Free persistent GPU Fresnel streaming context
 * @ctx: Stream context to release.
 */
void atmturb_cuda_fresnel_stream_free(
    atmturb_cuda_fresnel_stream_t *ctx)
{
    if (ctx == NULL)
    {
        return;
    }

    if (ctx->stream != 0)
    {
        cudaStreamSynchronize(ctx->stream);
        cudaStreamDestroy(ctx->stream);
        ctx->stream = 0;
    }
    if (ctx->plan_c2c != 0)
    {
        cufftDestroy(ctx->plan_c2c);
        ctx->plan_c2c = 0;
    }

    if (ctx->d_tf_pri != NULL)
    {
        for (int m = 0; m < ctx->plan.nsuper; m++)
        {
            if (ctx->d_tf_pri[m] != NULL)
            {
                cudaFree(ctx->d_tf_pri[m]);
            }
        }
        free(ctx->d_tf_pri);
        ctx->d_tf_pri = NULL;
    }

    if (ctx->d_tf_sec != NULL)
    {
        for (int m = 0; m < ctx->plan.nsuper; m++)
        {
            if (ctx->d_tf_sec[m] != NULL)
            {
                cudaFree(ctx->d_tf_sec[m]);
            }
        }
        free(ctx->d_tf_sec);
        ctx->d_tf_sec = NULL;
    }

    cudaFree(ctx->d_field_pri);
    cudaFree(ctx->d_field_sec);
    cudaFree(ctx->d_super_pha);
    cudaFree(ctx->d_super_spha);
    cudaFree(ctx->d_frame_pha);
    cudaFree(ctx->d_frame_amp);
    cudaFree(ctx->d_frame_spha);
    cudaFree(ctx->d_frame_samp);
    cudaFree(ctx->d_sublayers);

    if (ctx->d_lf_modes != NULL)
    {
        cudaFree(ctx->d_lf_modes);
        ctx->d_lf_modes = NULL;
    }
    if (ctx->h_lf_modes != NULL)
    {
        free(ctx->h_lf_modes);
        ctx->h_lf_modes = NULL;
    }
    if (ctx->super_nlf_modes != NULL)
    {
        free(ctx->super_nlf_modes);
        ctx->super_nlf_modes = NULL;
    }

    free(ctx->h_sublayers);
    free(ctx->masters);
    free(ctx->screen_offsets);

    atmturb_fresnel_plan_free(&ctx->plan);
    free(ctx);
}

static int atmturb_cuda_fresnel_init_tf(
    atmturb_cuda_fresnel_stream_t *ctx,
    size_t                         cpix_bytes)
{
    ctx->d_tf_pri = (cufftComplex **) calloc((size_t) ctx->plan.nsuper, sizeof(cufftComplex *));
    if (ctx->has_sec)
    {
        ctx->d_tf_sec = (cufftComplex **) calloc((size_t) ctx->plan.nsuper,
                                                 sizeof(cufftComplex *));
    }

    for (int m = 0; m < ctx->plan.nsuper; m++)
    {
        if (cudaMalloc((void **) &ctx->d_tf_pri[m], cpix_bytes) != cudaSuccess)
        {
            return -1;
        }
        cudaMemcpyAsync(ctx->d_tf_pri[m], ctx->plan.tf_pri[m], cpix_bytes,
                        cudaMemcpyHostToDevice, ctx->stream);

        if (ctx->has_sec && ctx->plan.tf_sec != NULL && ctx->plan.tf_sec[m] != NULL)
        {
            if (cudaMalloc((void **) &ctx->d_tf_sec[m], cpix_bytes) != cudaSuccess)
            {
                return -1;
            }
            cudaMemcpyAsync(ctx->d_tf_sec[m], ctx->plan.tf_sec[m], cpix_bytes,
                            cudaMemcpyHostToDevice, ctx->stream);
        }
    }
    return 0;
}

static int atmturb_cuda_fresnel_init_lf(
    atmturb_cuda_fresnel_stream_t *ctx,
    const atmturb_rolling_t       *r)
{
    ctx->lowfreq = r->lowfreq;
    if (r->lowfreq == 0)
    {
        return 0;
    }

    int max_lf_modes = 0;
    for (int m = 0; m < ctx->plan.nsuper; m++)
    {
        int n_lf = ctx->plan.supers[m].nlayers * 24;
        if (n_lf > max_lf_modes)
        {
            max_lf_modes = n_lf;
        }
    }
    ctx->max_lf_modes = max_lf_modes;

    if (max_lf_modes > 0)
    {
        size_t n_lf_tot = (size_t) (ctx->plan.nsuper * max_lf_modes);
        ctx->h_lf_modes = (atmturb_cuda_lowfreq_mode_t *) calloc(
            n_lf_tot, sizeof(atmturb_cuda_lowfreq_mode_t));
        ctx->super_nlf_modes = (int *) malloc(sizeof(int) * (size_t) ctx->plan.nsuper);
        if (!ctx->h_lf_modes || !ctx->super_nlf_modes ||
            cudaMalloc((void **) &ctx->d_lf_modes,
                       n_lf_tot * sizeof(atmturb_cuda_lowfreq_mode_t)) != cudaSuccess)
        {
            return -1;
        }
        for (int m = 0; m < ctx->plan.nsuper; m++)
        {
            ctx->super_nlf_modes[m] = ctx->plan.supers[m].nlayers * 24;
        }
    }
    return 0;
}

/**
 * atmturb_cuda_fresnel_stream_init - Initialize persistent GPU Fresnel streaming context
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @prof: Active turbulence profile.
 * @params: Observation parameters container.
 * @master_size: Master screen linear dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 *
 * Return: Allocated streaming context, or NULL on failure.
 */
atmturb_cuda_fresnel_stream_t *atmturb_cuda_fresnel_stream_init(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    const atmturb_profile_t    *prof,
    const atmturb_obs_params_t *params,
    long                        master_size,
    long                        pup_size)
{
    if (r == NULL || geom == NULL || prof == NULL || params == NULL ||
        master_size <= 0 || pup_size <= 0 || !atmturb_cuda_device_available())
    {
        return NULL;
    }

    atmturb_cuda_fresnel_stream_t *ctx =
        (atmturb_cuda_fresnel_stream_t *) calloc(1, sizeof(atmturb_cuda_fresnel_stream_t));
    if (ctx == NULL)
    {
        return NULL;
    }

    double z_bin = (double) CONF_FRESNEL_PROPAGATION_BIN;
    if (atmturb_fresnel_plan_init(&ctx->plan, prof, geom, pup_size, params->pupil_scale_m,
                                  params->lambda_ref_m, params->lambda_s_m, z_bin) != 0)
    {
        free(ctx);
        return NULL;
    }

    ctx->msize    = master_size;
    ctx->pup_size = pup_size;
    ctx->npix     = pup_size * pup_size;
    ctx->has_sec  = (params->lambda_s_m > 0.0) ? 1 : 0;
    ctx->interp   = geom->interp;
    ctx->os       = (int) geom->oversample;

    int max_sublayers = 1;
    for (int m = 0; m < ctx->plan.nsuper; m++)
    {
        int n_sub = 2 * ctx->plan.supers[m].nlayers;
        if (n_sub > max_sublayers)
        {
            max_sublayers = n_sub;
        }
    }
    ctx->max_sublayers = max_sublayers;

    if (atmturb_fresnel_cuda_init_screen_pool(r, geom->nlayers, &ctx->masters,
                                              &ctx->screen_offsets, &ctx->total_screens) != 0)
    {
        atmturb_cuda_fresnel_stream_free(ctx);
        return NULL;
    }

    ctx->d_masters = atmturb_cuda_sync_device_masters(
        master_size, ctx->total_screens, ctx->masters);
    if (ctx->d_masters == NULL)
    {
        atmturb_cuda_fresnel_stream_free(ctx);
        return NULL;
    }

    if (cudaStreamCreate(&ctx->stream) != cudaSuccess)
    {
        atmturb_cuda_fresnel_stream_free(ctx);
        return NULL;
    }

    if (cufftPlan2d(&ctx->plan_c2c, (int) pup_size, (int) pup_size, CUFFT_C2C) != CUFFT_SUCCESS)
    {
        atmturb_cuda_fresnel_stream_free(ctx);
        return NULL;
    }
    cufftSetStream(ctx->plan_c2c, ctx->stream);

    size_t npix_bytes = (size_t) ctx->npix * sizeof(float);
    size_t cpix_bytes = (size_t) ctx->npix * sizeof(cufftComplex);

    if (atmturb_cuda_fresnel_init_tf(ctx, cpix_bytes) != 0 ||
        atmturb_cuda_fresnel_init_lf(ctx, r) != 0)
    {
        atmturb_cuda_fresnel_stream_free(ctx);
        return NULL;
    }

    size_t n_sub_tot = (size_t) (ctx->plan.nsuper * max_sublayers);
    ctx->h_sublayers = (atmturb_cuda_rytov_sublayer_t *) calloc(
        n_sub_tot, sizeof(atmturb_cuda_rytov_sublayer_t));

    if (cudaMalloc((void **) &ctx->d_field_pri, cpix_bytes) != cudaSuccess ||
        cudaMalloc((void **) &ctx->d_super_pha, npix_bytes) != cudaSuccess ||
        cudaMalloc((void **) &ctx->d_frame_pha, npix_bytes) != cudaSuccess ||
        cudaMalloc((void **) &ctx->d_frame_amp, npix_bytes) != cudaSuccess ||
        cudaMalloc((void **) &ctx->d_sublayers,
                   n_sub_tot * sizeof(atmturb_cuda_rytov_sublayer_t)) != cudaSuccess)
    {
        atmturb_cuda_fresnel_stream_free(ctx);
        return NULL;
    }

    if (ctx->has_sec)
    {
        if (cudaMalloc((void **) &ctx->d_field_sec, cpix_bytes) != cudaSuccess ||
            cudaMalloc((void **) &ctx->d_super_spha, npix_bytes) != cudaSuccess ||
            cudaMalloc((void **) &ctx->d_frame_spha, npix_bytes) != cudaSuccess ||
            cudaMalloc((void **) &ctx->d_frame_samp, npix_bytes) != cudaSuccess)
        {
            atmturb_cuda_fresnel_stream_free(ctx);
            return NULL;
        }
    }

    cudaStreamSynchronize(ctx->stream);
    return ctx;
}

/**
 * atmturb_cuda_fresnel_stream_render_step - Render one Fresnel stream frame on GPU
 * @ctx: Persistent stream context.
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @t: Frame index.
 * @time_step_s: Time step in seconds.
 * @pha: Destination primary phase frame buffer.
 * @amp: Destination primary amplitude frame buffer (or NULL).
 * @spha: Destination secondary phase frame buffer (or NULL).
 * @samp: Destination secondary amplitude frame buffer (or NULL).
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_cuda_fresnel_stream_render_step(
    atmturb_cuda_fresnel_stream_t *ctx,
    const atmturb_rolling_t       *r,
    const atmturb_geom_t          *geom,
    long                           t,
    double                         time_step_s,
    float                         *pha,
    float                         *amp,
    float                         *spha,
    float                         *samp)
{
    if (ctx == NULL || pha == NULL)
    {
        return -1;
    }

    atmturb_fresnel_cuda_build_sublayers(
        r, geom, &ctx->plan, ctx->screen_offsets, t, time_step_s,
        ctx->max_sublayers, ctx->h_sublayers);

    size_t sub_bytes = (size_t) (ctx->plan.nsuper * ctx->max_sublayers) *
                       sizeof(atmturb_cuda_rytov_sublayer_t);
    cudaMemcpyAsync(ctx->d_sublayers, ctx->h_sublayers, sub_bytes,
                    cudaMemcpyHostToDevice, ctx->stream);

    if (ctx->lowfreq && ctx->h_lf_modes != NULL)
    {
        atmturb_fresnel_cuda_build_lf_modes(
            r, geom, &ctx->plan, t, time_step_s,
            ctx->max_lf_modes, ctx->h_lf_modes);
        size_t lf_bytes = (size_t) (ctx->plan.nsuper * ctx->max_lf_modes) *
                          sizeof(atmturb_cuda_lowfreq_mode_t);
        cudaMemcpyAsync(ctx->d_lf_modes, ctx->h_lf_modes, lf_bytes,
                        cudaMemcpyHostToDevice, ctx->stream);
    }

    size_t npix_bytes = (size_t) ctx->npix * sizeof(float);
    cudaMemsetAsync(ctx->d_frame_pha, 0, npix_bytes, ctx->stream);
    if (ctx->has_sec && spha != NULL)
    {
        cudaMemsetAsync(ctx->d_frame_spha, 0, npix_bytes, ctx->stream);
    }

    dim3 block(256);
    dim3 grid(((int) ctx->npix + block.x - 1) / block.x);

    for (int m = 0; m < ctx->plan.nsuper; m++)
    {
        const atmturb_superlayer_t *sl = &ctx->plan.supers[m];
        size_t base = (size_t) (m * ctx->max_sublayers);
        const atmturb_cuda_rytov_sublayer_t *d_sl = ctx->d_sublayers + base;
        int n_sub = 2 * sl->nlayers;

        atmturb_cuda_extrude_sl_launch(
            ctx->interp, ctx->d_masters, (int) ctx->msize, (int) ctx->pup_size,
            (int) ctx->pup_size, 0, n_sub, ctx->has_sec, ctx->os, d_sl,
            ctx->d_super_pha, ctx->d_super_spha, ctx->d_frame_pha, ctx->d_frame_spha,
            ctx->stream);

        if (ctx->lowfreq && ctx->super_nlf_modes != NULL &&
            ctx->super_nlf_modes[m] > 0 && ctx->d_lf_modes != NULL)
        {
            const atmturb_cuda_lowfreq_mode_t *d_lf =
                ctx->d_lf_modes + (size_t) (m * ctx->max_lf_modes);
            atmturb_cuda_lowfreq_accumulate_launch(
                (int) ctx->pup_size, (int) ctx->pup_size, 0,
                ctx->super_nlf_modes[m], ctx->has_sec, d_lf,
                ctx->d_super_pha, ctx->d_super_spha,
                ctx->d_frame_pha, ctx->d_frame_spha,
                ctx->stream);
        }

        if (m == 0)
        {
            atmturb_cuda_fresnel_init_modulate_kernel<<<grid, block, 0, ctx->stream>>>(
                ctx->d_field_pri, ctx->d_super_pha, (int) ctx->npix);
        }
        else
        {
            atmturb_cuda_fresnel_modulate_kernel<<<grid, block, 0, ctx->stream>>>(
                ctx->d_field_pri, ctx->d_super_pha, (int) ctx->npix);
        }

        cufftExecC2C(ctx->plan_c2c, ctx->d_field_pri, ctx->d_field_pri, CUFFT_FORWARD);
        atmturb_cuda_fresnel_tf_mul_kernel<<<grid, block, 0, ctx->stream>>>(
            ctx->d_field_pri, ctx->d_tf_pri[m], (int) ctx->npix);
        cufftExecC2C(ctx->plan_c2c, ctx->d_field_pri, ctx->d_field_pri, CUFFT_INVERSE);

        if (ctx->has_sec && (spha != NULL || samp != NULL) && ctx->d_tf_sec[m] != NULL)
        {
            if (m == 0)
            {
                atmturb_cuda_fresnel_init_modulate_kernel<<<grid, block, 0, ctx->stream>>>(
                    ctx->d_field_sec, ctx->d_super_spha, (int) ctx->npix);
            }
            else
            {
                atmturb_cuda_fresnel_modulate_kernel<<<grid, block, 0, ctx->stream>>>(
                    ctx->d_field_sec, ctx->d_super_spha, (int) ctx->npix);
            }

            cufftExecC2C(ctx->plan_c2c, ctx->d_field_sec, ctx->d_field_sec, CUFFT_FORWARD);
            atmturb_cuda_fresnel_tf_mul_kernel<<<grid, block, 0, ctx->stream>>>(
                ctx->d_field_sec, ctx->d_tf_sec[m], (int) ctx->npix);
            cufftExecC2C(ctx->plan_c2c, ctx->d_field_sec, ctx->d_field_sec, CUFFT_INVERSE);
        }
    }

    atmturb_cuda_fresnel_extract_kernel<<<grid, block, 0, ctx->stream>>>(
        ctx->d_field_pri, ctx->d_frame_pha, (amp != NULL) ? ctx->d_frame_amp : NULL,
        (int) ctx->npix);

    if (ctx->has_sec && (spha != NULL || samp != NULL))
    {
        atmturb_cuda_fresnel_extract_kernel<<<grid, block, 0, ctx->stream>>>(
            ctx->d_field_sec, ctx->d_frame_spha, (samp != NULL) ? ctx->d_frame_samp : NULL,
            (int) ctx->npix);
    }

    cudaMemcpyAsync(pha, ctx->d_frame_pha, npix_bytes, cudaMemcpyDeviceToHost, ctx->stream);
    if (amp != NULL)
    {
        cudaMemcpyAsync(amp, ctx->d_frame_amp, npix_bytes, cudaMemcpyDeviceToHost, ctx->stream);
    }
    if (ctx->has_sec && spha != NULL)
    {
        cudaMemcpyAsync(spha, ctx->d_frame_spha, npix_bytes, cudaMemcpyDeviceToHost, ctx->stream);
    }
    if (ctx->has_sec && samp != NULL)
    {
        cudaMemcpyAsync(samp, ctx->d_frame_samp, npix_bytes, cudaMemcpyDeviceToHost, ctx->stream);
    }

    cudaStreamSynchronize(ctx->stream);
    return 0;
}
