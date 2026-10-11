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
    atmturb_cuda_rytov_sublayer_t *h_sublayers[2];
    atmturb_cuda_rytov_sublayer_t *d_sublayers[2];
    atmturb_cuda_lowfreq_mode_t   *h_lf_modes[2];
    atmturb_cuda_lowfreq_mode_t   *d_lf_modes[2];
    const float                   *d_masters;

    cudaStream_t                   stream;
    cudaEvent_t                    ev_done[2];
    cufftHandle                    plan_c2c;
    cufftHandle                    plan_batch2;

    cufftComplex                 **d_tf_pri;
    cufftComplex                 **d_tf_sec;

    cufftComplex                  *d_field_pri;
    cufftComplex                  *d_field_sec;
    float                         *d_super_pha;
    float                         *d_super_spha;
    float                         *d_frame_pha[2];
    float                         *d_frame_amp[2];
    float                         *d_frame_spha[2];
    float                         *d_frame_samp[2];
    float                         *d_sum_buf[2];
    cudaGraphExec_t                graph_exec[2];
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

__global__ static void atmturb_cuda_fresnel_tf_mul_vec2_kernel(
    float4       *d_field,
    const float4 *d_tf,
    int           npairs)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= npairs)
    {
        return;
    }
    float4 f = d_field[i];
    float4 t = d_tf[i];
    float4 out;
    out.x = f.x * t.x - f.y * t.y;
    out.y = f.y * t.x + f.x * t.y;
    out.z = f.z * t.z - f.w * t.w;
    out.w = f.w * t.z + f.z * t.w;
    d_field[i] = out;
}

__global__ static void atmturb_cuda_fresnel_extract_kernel(
    const cufftComplex *d_field,
    float              *d_pha,
    float              *d_amp,
    float              *d_sum,
    int                 npix)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float p = 0.0f;
    if (i < npix)
    {
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
        p = phi_geo + atan2f(d_im, d_re);
        d_pha[i] = p;
    }

    if (d_sum != NULL)
    {
        for (int offset = 16; offset > 0; offset /= 2)
        {
            p += __shfl_down_sync(0xffffffff, p, offset);
        }
        __shared__ float s_warp_sums[32];
        int lane = threadIdx.x % 32;
        int warp_id = threadIdx.x / 32;
        if (lane == 0)
        {
            s_warp_sums[warp_id] = p;
        }
        __syncthreads();

        if (threadIdx.x < 32)
        {
            int nwarps = blockDim.x / 32;
            float bval = (threadIdx.x < nwarps) ? s_warp_sums[threadIdx.x] : 0.0f;
            for (int offset = 16; offset > 0; offset /= 2)
            {
                bval += __shfl_down_sync(0xffffffff, bval, offset);
            }
            if (threadIdx.x == 0)
            {
                atomicAdd(d_sum, bval);
            }
        }
    }
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

    for (int s = 0; s < 2; s++)
    {
        if (ctx->ev_done[s] != 0)
        {
            cudaEventDestroy(ctx->ev_done[s]);
            ctx->ev_done[s] = 0;
        }
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
    if (ctx->plan_batch2 != 0)
    {
        cufftDestroy(ctx->plan_batch2);
        ctx->plan_batch2 = 0;
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
    ctx->d_field_pri = NULL;
    ctx->d_field_sec = NULL;
    cudaFree(ctx->d_super_pha);
    cudaFree(ctx->d_super_spha);

    for (int s = 0; s < 2; s++)
    {
        if (ctx->graph_exec[s] != NULL)
        {
            cudaGraphExecDestroy(ctx->graph_exec[s]);
            ctx->graph_exec[s] = NULL;
        }
        cudaFree(ctx->d_frame_pha[s]);
        cudaFree(ctx->d_frame_amp[s]);
        cudaFree(ctx->d_frame_spha[s]);
        cudaFree(ctx->d_frame_samp[s]);
        cudaFree(ctx->d_sum_buf[s]);
        cudaFree(ctx->d_sublayers[s]);
        if (ctx->h_sublayers[s] != NULL)
        {
            cudaFreeHost(ctx->h_sublayers[s]);
            ctx->h_sublayers[s] = NULL;
        }
        if (ctx->d_lf_modes[s] != NULL)
        {
            cudaFree(ctx->d_lf_modes[s]);
            ctx->d_lf_modes[s] = NULL;
        }
        if (ctx->h_lf_modes[s] != NULL)
        {
            cudaFreeHost(ctx->h_lf_modes[s]);
            ctx->h_lf_modes[s] = NULL;
        }
    }

    if (ctx->super_nlf_modes != NULL)
    {
        free(ctx->super_nlf_modes);
        ctx->super_nlf_modes = NULL;
    }

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
        size_t lf_bytes = n_lf_tot * sizeof(atmturb_cuda_lowfreq_mode_t);
        ctx->super_nlf_modes = (int *) malloc(sizeof(int) * (size_t) ctx->plan.nsuper);
        if (!ctx->super_nlf_modes)
        {
            return -1;
        }
        for (int s = 0; s < 2; s++)
        {
            if (cudaMallocHost((void **) &ctx->h_lf_modes[s], lf_bytes) != cudaSuccess ||
                cudaMalloc((void **) &ctx->d_lf_modes[s], lf_bytes) != cudaSuccess)
            {
                return -1;
            }
            memset(ctx->h_lf_modes[s], 0, lf_bytes);
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

    for (int s = 0; s < 2; s++)
    {
        if (cudaEventCreateWithFlags(&ctx->ev_done[s], cudaEventDisableTiming) != cudaSuccess)
        {
            atmturb_cuda_fresnel_stream_free(ctx);
            return NULL;
        }
    }

    if (cufftPlan2d(&ctx->plan_c2c, (int) pup_size, (int) pup_size, CUFFT_C2C) != CUFFT_SUCCESS)
    {
        atmturb_cuda_fresnel_stream_free(ctx);
        return NULL;
    }
    cufftSetStream(ctx->plan_c2c, ctx->stream);

    if (ctx->has_sec)
    {
        int n[2] = {(int) pup_size, (int) pup_size};
        if (cufftPlanMany(&ctx->plan_batch2, 2, n, NULL, 1, (int) ctx->npix,
                          NULL, 1, (int) ctx->npix, CUFFT_C2C, 2) != CUFFT_SUCCESS)
        {
            atmturb_cuda_fresnel_stream_free(ctx);
            return NULL;
        }
        cufftSetStream(ctx->plan_batch2, ctx->stream);
    }

    size_t npix_bytes = (size_t) ctx->npix * sizeof(float);
    size_t cpix_bytes = (size_t) ctx->npix * sizeof(cufftComplex);
    size_t field_alloc_bytes = (ctx->has_sec ? 2 : 1) * cpix_bytes;

    if (atmturb_cuda_fresnel_init_tf(ctx, cpix_bytes) != 0 ||
        atmturb_cuda_fresnel_init_lf(ctx, r) != 0)
    {
        atmturb_cuda_fresnel_stream_free(ctx);
        return NULL;
    }

    size_t n_sub_tot = (size_t) (ctx->plan.nsuper * max_sublayers);
    size_t sub_bytes = n_sub_tot * sizeof(atmturb_cuda_rytov_sublayer_t);

    for (int s = 0; s < 2; s++)
    {
        if (cudaMallocHost((void **) &ctx->h_sublayers[s], sub_bytes) != cudaSuccess ||
            cudaMalloc((void **) &ctx->d_sublayers[s], sub_bytes) != cudaSuccess ||
            cudaMalloc((void **) &ctx->d_frame_pha[s], npix_bytes) != cudaSuccess ||
            cudaMalloc((void **) &ctx->d_frame_amp[s], npix_bytes) != cudaSuccess ||
            cudaMalloc((void **) &ctx->d_sum_buf[s], sizeof(float)) != cudaSuccess)
        {
            atmturb_cuda_fresnel_stream_free(ctx);
            return NULL;
        }
        memset(ctx->h_sublayers[s], 0, sub_bytes);

        if (ctx->has_sec)
        {
            if (cudaMalloc((void **) &ctx->d_frame_spha[s], npix_bytes) != cudaSuccess ||
                cudaMalloc((void **) &ctx->d_frame_samp[s], npix_bytes) != cudaSuccess)
            {
                atmturb_cuda_fresnel_stream_free(ctx);
                return NULL;
            }
        }
    }

    if (cudaMalloc((void **) &ctx->d_field_pri, field_alloc_bytes) != cudaSuccess ||
        cudaMalloc((void **) &ctx->d_super_pha, npix_bytes) != cudaSuccess)
    {
        atmturb_cuda_fresnel_stream_free(ctx);
        return NULL;
    }

    if (ctx->has_sec)
    {
        ctx->d_field_sec = ctx->d_field_pri + ctx->npix;
        if (cudaMalloc((void **) &ctx->d_super_spha, npix_bytes) != cudaSuccess)
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
/**
 * atmturb_cuda_fresnel_propagate_superlayer - Modulate, propagate via cuFFT, and filter
 * @ctx: Persistent stream context.
 * @m: Superlayer index.
 * @do_batch2: 1 to execute batched 2-channel cuFFT, 0 for single channel.
 * @grid: Grid dimensions for element-wise kernels.
 * @block: Block dimensions for element-wise kernels.
 */
static void atmturb_cuda_fresnel_propagate_superlayer(
    atmturb_cuda_fresnel_stream_t *ctx,
    int                            m,
    int                            do_batch2,
    dim3                           grid,
    dim3                           block)
{
    if (do_batch2)
    {
        if (m == 0)
        {
            atmturb_cuda_fresnel_init_modulate_kernel<<<grid, block, 0, ctx->stream>>>(
                ctx->d_field_pri, ctx->d_super_pha, (int) ctx->npix);
            atmturb_cuda_fresnel_init_modulate_kernel<<<grid, block, 0, ctx->stream>>>(
                ctx->d_field_sec, ctx->d_super_spha, (int) ctx->npix);
        }
        else
        {
            atmturb_cuda_fresnel_modulate_kernel<<<grid, block, 0, ctx->stream>>>(
                ctx->d_field_pri, ctx->d_super_pha, (int) ctx->npix);
            atmturb_cuda_fresnel_modulate_kernel<<<grid, block, 0, ctx->stream>>>(
                ctx->d_field_sec, ctx->d_super_spha, (int) ctx->npix);
        }

        cufftExecC2C(ctx->plan_batch2, ctx->d_field_pri, ctx->d_field_pri, CUFFT_FORWARD);

        int npairs = (int) ctx->npix / 2;
        dim3 b_pairs(256);
        dim3 g_pairs((npairs + 256 - 1) / 256);

        atmturb_cuda_fresnel_tf_mul_vec2_kernel<<<g_pairs, b_pairs, 0, ctx->stream>>>(
            (float4 *) ctx->d_field_pri, (const float4 *) ctx->d_tf_pri[m], npairs);
        atmturb_cuda_fresnel_tf_mul_vec2_kernel<<<g_pairs, b_pairs, 0, ctx->stream>>>(
            (float4 *) ctx->d_field_sec, (const float4 *) ctx->d_tf_sec[m], npairs);

        cufftExecC2C(ctx->plan_batch2, ctx->d_field_pri, ctx->d_field_pri, CUFFT_INVERSE);
    }
    else
    {
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
        int npairs = (int) ctx->npix / 2;
        dim3 b_pairs(256);
        dim3 g_pairs((npairs + 256 - 1) / 256);
        atmturb_cuda_fresnel_tf_mul_vec2_kernel<<<g_pairs, b_pairs, 0, ctx->stream>>>(
            (float4 *) ctx->d_field_pri, (const float4 *) ctx->d_tf_pri[m], npairs);
        cufftExecC2C(ctx->plan_c2c, ctx->d_field_pri, ctx->d_field_pri, CUFFT_INVERSE);
    }
}

/**
 * atmturb_cuda_fresnel_step_ops - Enqueue GPU operations for one Fresnel propagation frame
 * @ctx: Persistent stream context.
 * @slot: Active double-buffer slot index (0 or 1).
 * @has_amp: 1 if amplitude is requested, 0 otherwise.
 * @has_sec_out: 1 if secondary wavelength output is active, 0 otherwise.
 */
static void atmturb_cuda_fresnel_step_ops(
    atmturb_cuda_fresnel_stream_t *ctx,
    int                            slot,
    int                            has_amp,
    int                            has_sec_out)
{
    size_t npix_bytes = (size_t) ctx->npix * sizeof(float);
    cudaMemsetAsync(ctx->d_frame_pha[slot], 0, npix_bytes, ctx->stream);
    if (ctx->has_sec && has_sec_out)
    {
        cudaMemsetAsync(ctx->d_frame_spha[slot], 0, npix_bytes, ctx->stream);
    }

    dim3 block(256);
    dim3 grid(((int) ctx->npix + block.x - 1) / block.x);

    for (int m = 0; m < ctx->plan.nsuper; m++)
    {
        const atmturb_superlayer_t *sl = &ctx->plan.supers[m];
        size_t base = (size_t) (m * ctx->max_sublayers);
        const atmturb_cuda_rytov_sublayer_t *d_sl = ctx->d_sublayers[slot] + base;
        int n_sub = 2 * sl->nlayers;

        int n_lf = (ctx->lowfreq && ctx->super_nlf_modes != NULL && ctx->d_lf_modes[slot] != NULL)
                       ? ctx->super_nlf_modes[m] : 0;
        const atmturb_cuda_lowfreq_mode_t *d_lf = (n_lf > 0)
                       ? ctx->d_lf_modes[slot] + (size_t) (m * ctx->max_lf_modes) : NULL;

        atmturb_cuda_extrude_sl_launch(
            ctx->interp, ctx->d_masters, (int) ctx->msize, (int) ctx->pup_size,
            (int) ctx->pup_size, 0, n_sub, ctx->has_sec, ctx->os, d_sl,
            n_lf, d_lf,
            ctx->d_super_pha, ctx->d_super_spha, ctx->d_frame_pha[slot], ctx->d_frame_spha[slot],
            ctx->stream);

        int do_batch2 = (ctx->has_sec && has_sec_out && ctx->d_tf_sec[m] != NULL);
        atmturb_cuda_fresnel_propagate_superlayer(ctx, m, do_batch2, grid, block);
    }

    if (ctx->d_sum_buf[slot] != NULL)
    {
        cudaMemsetAsync(ctx->d_sum_buf[slot], 0, sizeof(float), ctx->stream);
    }

    atmturb_cuda_fresnel_extract_kernel<<<grid, block, 0, ctx->stream>>>(
        ctx->d_field_pri, ctx->d_frame_pha[slot], (has_amp) ? ctx->d_frame_amp[slot] : NULL,
        ctx->d_sum_buf[slot], (int) ctx->npix);

    if (ctx->has_sec && has_sec_out)
    {
        atmturb_cuda_fresnel_extract_kernel<<<grid, block, 0, ctx->stream>>>(
            ctx->d_field_sec, ctx->d_frame_spha[slot], (has_amp) ? ctx->d_frame_samp[slot] : NULL,
            NULL, (int) ctx->npix);
    }

    if (ctx->d_sum_buf[slot] != NULL)
    {
        atmturb_cuda_piston_sub_launch(
            ctx->d_frame_pha[slot], ctx->d_sum_buf[slot], (int) ctx->npix, ctx->stream);
        if (ctx->has_sec && has_sec_out)
        {
            atmturb_cuda_remove_piston_launch(
                ctx->d_frame_spha[slot], ctx->d_sum_buf[slot], (int) ctx->npix, ctx->stream);
        }
    }
}

int atmturb_cuda_fresnel_stream_render_launch(
    atmturb_cuda_fresnel_stream_t *ctx,
    const atmturb_rolling_t       *r,
    const atmturb_geom_t          *geom,
    long                           t,
    double                         time_step_s,
    float                         *pha,
    float                         *amp,
    float                         *spha,
    float                         *samp,
    int                            slot)
{
    if (ctx == NULL || pha == NULL)
    {
        return -1;
    }
    int s = slot & 1;

    atmturb_fresnel_cuda_build_sublayers(
        r, geom, &ctx->plan, ctx->screen_offsets, t, time_step_s,
        ctx->max_sublayers, ctx->h_sublayers[s]);

    size_t sub_bytes = (size_t) (ctx->plan.nsuper * ctx->max_sublayers) *
                       sizeof(atmturb_cuda_rytov_sublayer_t);
    cudaMemcpyAsync(ctx->d_sublayers[s], ctx->h_sublayers[s], sub_bytes,
                    cudaMemcpyHostToDevice, ctx->stream);

    if (ctx->lowfreq && ctx->h_lf_modes[s] != NULL)
    {
        atmturb_fresnel_cuda_build_lf_modes(
            r, geom, &ctx->plan, t, time_step_s,
            ctx->max_lf_modes, ctx->h_lf_modes[s]);
        size_t lf_bytes = (size_t) (ctx->plan.nsuper * ctx->max_lf_modes) *
                          sizeof(atmturb_cuda_lowfreq_mode_t);
        cudaMemcpyAsync(ctx->d_lf_modes[s], ctx->h_lf_modes[s], lf_bytes,
                        cudaMemcpyHostToDevice, ctx->stream);
    }

    int has_amp     = (amp != NULL) ? 1 : 0;
    int has_sec_out = (ctx->has_sec && (spha != NULL || samp != NULL)) ? 1 : 0;

    if (ctx->graph_exec[s] == NULL)
    {
        cudaGraph_t graph;
        if (cudaStreamBeginCapture(ctx->stream, cudaStreamCaptureModeGlobal) == cudaSuccess)
        {
            atmturb_cuda_fresnel_step_ops(ctx, s, has_amp, has_sec_out);
            if (cudaStreamEndCapture(ctx->stream, &graph) == cudaSuccess)
            {
                if (cudaGraphInstantiate(&ctx->graph_exec[s], graph, NULL, NULL, 0) != cudaSuccess)
                {
                    ctx->graph_exec[s] = NULL;
                }
                cudaGraphDestroy(graph);
            }
        }
    }

    if (ctx->graph_exec[s] != NULL)
    {
        cudaGraphLaunch(ctx->graph_exec[s], ctx->stream);
    }
    else
    {
        atmturb_cuda_fresnel_step_ops(ctx, s, has_amp, has_sec_out);
    }

    size_t npix_bytes = (size_t) ctx->npix * sizeof(float);
    cudaMemcpyAsync(pha, ctx->d_frame_pha[s], npix_bytes, cudaMemcpyDeviceToHost, ctx->stream);
    if (amp != NULL)
    {
        cudaMemcpyAsync(amp, ctx->d_frame_amp[s], npix_bytes, cudaMemcpyDeviceToHost, ctx->stream);
    }
    if (ctx->has_sec && spha != NULL)
    {
        cudaMemcpyAsync(spha, ctx->d_frame_spha[s], npix_bytes, cudaMemcpyDeviceToHost, ctx->stream);
    }
    if (ctx->has_sec && samp != NULL)
    {
        cudaMemcpyAsync(samp, ctx->d_frame_samp[s], npix_bytes, cudaMemcpyDeviceToHost, ctx->stream);
    }

    cudaEventRecord(ctx->ev_done[s], ctx->stream);
    return 0;
}

int atmturb_cuda_fresnel_stream_render_wait(
    atmturb_cuda_fresnel_stream_t *ctx,
    int                            slot)
{
    if (ctx == NULL)
    {
        return -1;
    }
    int s = slot & 1;
    return (cudaEventSynchronize(ctx->ev_done[s]) == cudaSuccess) ? 0 : -1;
}

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
    int ret = atmturb_cuda_fresnel_stream_render_launch(
        ctx, r, geom, t, time_step_s, pha, amp, spha, samp, 0);
    if (ret != 0)
    {
        return ret;
    }
    return atmturb_cuda_fresnel_stream_render_wait(ctx, 0);
}
