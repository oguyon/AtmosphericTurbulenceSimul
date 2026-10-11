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
#include "atmturb_cuda_fresnel_kernels.h"
#include "atmturb_cuda_piston.h"
#include "atmturb_cuda_rytov_internal.h"
#include "atmturb_fresnel.h"
#include "atmturb_geometry.h"
#include "atmturb_lowfreq.h"
#include "atmturb_profile.h"
#include "atmturb_rolling.h"
#include "atmturb_superlayer.h"
#include "atmturb_wfs_cuda_common.h"

extern float CONF_FRESNEL_PROPAGATION_BIN;
}

/**
 * struct atmturb_cuda_fresnel_stream - Persistent GPU context for Fresnel propagation
 * @plan: Superlayer Fresnel diffraction plan.
 * @msize: Linear dimension of master phase screen.
 * @pup_size: Linear dimension of synthesized output pupil.
 * @npix: Total number of pixels on pupil grid (pup_size * pup_size).
 * @has_sec: 1 if secondary chromatic wavelength is enabled, 0 otherwise.
 * @interp: Interpolation scheme (0=bilinear, 1=Keys bicubic).
 * @os: Oversampling factor from master grid to simulation grid.
 * @total_screens: Total number of allocated master screens across all layers.
 * @max_sublayers: Maximum number of sublayers per superlayer.
 * @lowfreq: Flag indicating whether subharmonic low-frequency modes are enabled.
 * @max_lf_modes: Maximum low-frequency modes per superlayer.
 * @super_nlf_modes: Array of mode counts per superlayer.
 * @masters: Host pointers to master screen arrays.
 * @screen_offsets: Screen index offsets per layer.
 * @h_sublayers: Pinned host sublayer coordinate buffers (double buffered).
 * @d_sublayers: Device sublayer coordinate buffers (double buffered).
 * @h_lf_modes: Pinned host low-frequency mode buffers (double buffered).
 * @d_lf_modes: Device low-frequency mode buffers (double buffered).
 * @d_masters: Device master screens array.
 * @stream: CUDA execution stream.
 * @ev_done: CUDA completion events for double-buffer slots.
 * @plan_c2c: cuFFT 2D forward/inverse plan for single-channel propagation.
 * @plan_batch2: cuFFT 2D batched plan for simultaneous dual-wavelength propagation.
 * @d_tf_pri: Array of device pointers to primary wavelength transfer functions.
 * @d_tf_sec: Array of device pointers to secondary wavelength transfer functions.
 * @d_field_pri: Device complex electric field buffer for primary wavelength.
 * @d_field_sec: Device complex electric field buffer for secondary wavelength.
 * @d_super_pha: Device scratch buffer for primary superlayer phase screen.
 * @d_super_spha: Device scratch buffer for secondary superlayer phase screen.
 * @d_frame_pha: Double-buffered device destination buffers for primary phase.
 * @d_frame_amp: Double-buffered device destination buffers for primary amplitude.
 * @d_frame_spha: Double-buffered device destination buffers for secondary phase.
 * @d_frame_samp: Double-buffered device destination buffers for secondary amplitude.
 * @d_sum_buf: Double-buffered device atomic reduction accumulators for piston.
 * @graph_exec: Captured CUDA execution graphs for zero-overhead kernel replay.
 */
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

/**
 * atmturb_cuda_fresnel_init_tf - Upload Fresnel optical transfer functions to GPU
 * @ctx: Persistent stream context.
 * @cpix_bytes: Size in bytes of complex 2D frequency grid.
 *
 * Return: 0 on success, -1 on CUDA allocation failure.
 */
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

/**
 * atmturb_cuda_fresnel_init_lf - Allocate host and device subharmonic low-frequency mode buffers
 * @ctx: Persistent stream context.
 * @r: Rolling simulation context.
 *
 * Return: 0 on success, -1 on failure.
 */
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

    if (atmturb_wfs_cuda_init_screen_pool(r, geom->nlayers, &ctx->masters,
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
 * atmturb_cuda_fresnel_propagate_superlayer - Modulate field and apply Fresnel diffraction
 * @ctx: Persistent stream context.
 * @m: Superlayer index.
 * @do_batch2: 1 to execute batched 2-channel cuFFT, 0 for single channel.
 */
static void atmturb_cuda_fresnel_propagate_superlayer(
    atmturb_cuda_fresnel_stream_t *ctx,
    int                            m,
    int                            do_batch2)
{
    if (do_batch2)
    {
        if (m == 0)
        {
            atmturb_cuda_fresnel_init_modulate_launch(
                ctx->d_field_pri, ctx->d_super_pha, (int) ctx->npix, ctx->stream);
            atmturb_cuda_fresnel_init_modulate_launch(
                ctx->d_field_sec, ctx->d_super_spha, (int) ctx->npix, ctx->stream);
        }
        else
        {
            atmturb_cuda_fresnel_modulate_launch(
                ctx->d_field_pri, ctx->d_super_pha, (int) ctx->npix, ctx->stream);
            atmturb_cuda_fresnel_modulate_launch(
                ctx->d_field_sec, ctx->d_super_spha, (int) ctx->npix, ctx->stream);
        }

        cufftExecC2C(ctx->plan_batch2, ctx->d_field_pri, ctx->d_field_pri, CUFFT_FORWARD);

        atmturb_cuda_fresnel_tf_mul_vec2_launch(
            ctx->d_field_pri, ctx->d_tf_pri[m], (int) ctx->npix, ctx->stream);
        atmturb_cuda_fresnel_tf_mul_vec2_launch(
            ctx->d_field_sec, ctx->d_tf_sec[m], (int) ctx->npix, ctx->stream);

        cufftExecC2C(ctx->plan_batch2, ctx->d_field_pri, ctx->d_field_pri, CUFFT_INVERSE);
    }
    else
    {
        if (m == 0)
        {
            atmturb_cuda_fresnel_init_modulate_launch(
                ctx->d_field_pri, ctx->d_super_pha, (int) ctx->npix, ctx->stream);
        }
        else
        {
            atmturb_cuda_fresnel_modulate_launch(
                ctx->d_field_pri, ctx->d_super_pha, (int) ctx->npix, ctx->stream);
        }

        cufftExecC2C(ctx->plan_c2c, ctx->d_field_pri, ctx->d_field_pri, CUFFT_FORWARD);
        atmturb_cuda_fresnel_tf_mul_vec2_launch(
            ctx->d_field_pri, ctx->d_tf_pri[m], (int) ctx->npix, ctx->stream);
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
        atmturb_cuda_fresnel_propagate_superlayer(ctx, m, do_batch2);
    }

    if (ctx->d_sum_buf[slot] != NULL)
    {
        cudaMemsetAsync(ctx->d_sum_buf[slot], 0, sizeof(float), ctx->stream);
    }

    atmturb_cuda_fresnel_extract_launch(
        ctx->d_field_pri, ctx->d_frame_pha[slot], (has_amp) ? ctx->d_frame_amp[slot] : NULL,
        ctx->d_sum_buf[slot], (int) ctx->npix, ctx->stream);

    if (ctx->has_sec && has_sec_out)
    {
        atmturb_cuda_fresnel_extract_launch(
            ctx->d_field_sec, ctx->d_frame_spha[slot], (has_amp) ? ctx->d_frame_samp[slot] : NULL,
            NULL, (int) ctx->npix, ctx->stream);
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

/**
 * atmturb_cuda_fresnel_stream_render_launch - Asynchronously dispatch one Fresnel frame on GPU
 * @ctx: Persistent stream context.
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @t: Simulation frame index.
 * @time_step_s: Time step in seconds.
 * @pha: Destination host buffer for primary phase.
 * @amp: Destination host buffer for primary amplitude (or NULL).
 * @spha: Destination host buffer for secondary phase (or NULL).
 * @samp: Destination host buffer for secondary amplitude (or NULL).
 * @slot: Active double-buffer slot index (0 or 1).
 *
 * Return: 0 on success, -1 on failure.
 */
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

    atmturb_wfs_cuda_build_frame_sublayers(
        r, geom, ctx->plan.nsuper, ctx->plan.supers, ctx->screen_offsets,
        t, time_step_s, 0.0, ctx->max_sublayers, ctx->h_sublayers[s]);

    size_t sub_bytes = (size_t) (ctx->plan.nsuper * ctx->max_sublayers) *
                       sizeof(atmturb_cuda_rytov_sublayer_t);
    cudaMemcpyAsync(ctx->d_sublayers[s], ctx->h_sublayers[s], sub_bytes,
                    cudaMemcpyHostToDevice, ctx->stream);

    if (ctx->lowfreq && ctx->h_lf_modes[s] != NULL)
    {
        atmturb_wfs_cuda_build_frame_lf_modes(
            r, geom, ctx->plan.nsuper, ctx->plan.supers,
            t, time_step_s, 0.0, ctx->max_lf_modes, ctx->h_lf_modes[s]);
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
    cudaMemcpyAsync(
        pha, ctx->d_frame_pha[s], npix_bytes, cudaMemcpyDeviceToHost, ctx->stream);
    if (amp != NULL)
    {
        cudaMemcpyAsync(
            amp, ctx->d_frame_amp[s], npix_bytes, cudaMemcpyDeviceToHost, ctx->stream);
    }
    if (ctx->has_sec && spha != NULL)
    {
        cudaMemcpyAsync(
            spha, ctx->d_frame_spha[s], npix_bytes, cudaMemcpyDeviceToHost, ctx->stream);
    }
    if (ctx->has_sec && samp != NULL)
    {
        cudaMemcpyAsync(
            samp, ctx->d_frame_samp[s], npix_bytes, cudaMemcpyDeviceToHost, ctx->stream);
    }

    cudaEventRecord(ctx->ev_done[s], ctx->stream);
    return 0;
}

/**
 * atmturb_cuda_fresnel_stream_render_wait - Wait for in-flight GPU frame to finish
 * @ctx: Persistent stream context.
 * @slot: Double-buffer slot index (0 or 1).
 *
 * Return: 0 on success, -1 on failure.
 */
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

/**
 * atmturb_cuda_fresnel_stream_render_step - Synchronously render one Fresnel frame on GPU
 * @ctx: Persistent stream context.
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @t: Simulation frame index.
 * @time_step_s: Time step in seconds.
 * @pha: Destination host buffer for primary phase.
 * @amp: Destination host buffer for primary amplitude (or NULL).
 * @spha: Destination host buffer for secondary phase (or NULL).
 * @samp: Destination host buffer for secondary amplitude (or NULL).
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
    int ret = atmturb_cuda_fresnel_stream_render_launch(
        ctx, r, geom, t, time_step_s, pha, amp, spha, samp, 0);
    if (ret != 0)
    {
        return ret;
    }
    return atmturb_cuda_fresnel_stream_render_wait(ctx, 0);
}
