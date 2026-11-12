// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_fresnel_batch.cu
 * @brief   Multi-frame batched cuFFT 3D wavefront cube renderer for Fresnel propagation
 */

#include <cuda_runtime.h>
#include <cufft.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

extern "C" {
#include "atmturb_cuda.h"
#include "atmturb_cuda_fresnel.h"
#include "atmturb_cuda_fresnel_kernels.h"
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

typedef struct
{
    atmturb_fresnel_plan_t         plan;
    long                           msize;
    long                           pup_size;
    long                           npix;
    int                            batch_b;
    int                            batch_rem;
    int                            has_sec;
    int                            has_amp;
    int                            n_ch;
    int                            interp;
    int                            os;
    int                            total_screens;
    int                            max_sublayers;
    int                            lowfreq;
    int                            max_lf_modes;
    const float                  **masters;
    int                           *screen_offsets;
    const float                   *d_masters;
    cudaStream_t                   stream;
    cufftHandle                    plan_many;
    cufftHandle                    plan_rem;
    cufftComplex                 **d_tf_pri;
    cufftComplex                 **d_tf_sec;
    atmturb_cuda_rytov_sublayer_t *h_sub;
    atmturb_cuda_rytov_sublayer_t *d_sub;
    atmturb_cuda_lowfreq_mode_t   *h_lf;
    atmturb_cuda_lowfreq_mode_t   *d_lf;
    cufftComplex                  *d_fields;
    float                         *d_super_pha;
    float                         *d_super_spha;
    float                         *d_frame_pha;
    float                         *d_frame_amp;
    float                         *d_frame_spha;
    float                         *d_frame_samp;
} atmturb_fresnel_batch_ctx_t;

/**
 * atmturb_fresnel_tf_mul_batch_vec2_kernel - Multiply even-sized field batch by OTF via float4
 */
__global__ static void atmturb_fresnel_tf_mul_batch_vec2_kernel(
    float4       *d_field,
    const float4 *d_tf,
    int           npairs,
    int           count)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int b = (int) blockIdx.y;
    if (i >= npairs || b >= count)
    {
        return;
    }

    size_t idx = (size_t) b * (size_t) npairs + (size_t) i;
    float4 f   = d_field[idx];
    float4 t   = __ldg(&d_tf[i]);
    float4 out;
    out.x = f.x * t.x - f.y * t.y;
    out.y = f.y * t.x + f.x * t.y;
    out.z = f.z * t.z - f.w * t.w;
    out.w = f.w * t.z + f.z * t.w;
    d_field[idx] = out;
}

/**
 * atmturb_fresnel_tf_mul_batch_c2_kernel - Multiply odd-sized field batch by OTF via cufftComplex
 */
__global__ static void atmturb_fresnel_tf_mul_batch_c2_kernel(
    cufftComplex       *d_field,
    const cufftComplex *d_tf,
    int                 npix,
    int                 count)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int b = (int) blockIdx.y;
    if (i >= npix || b >= count)
    {
        return;
    }

    size_t       idx = (size_t) b * (size_t) npix + (size_t) i;
    cufftComplex f   = d_field[idx];
    cufftComplex t   = d_tf[i];
    d_field[idx]     = make_cuFloatComplex(f.x * t.x - f.y * t.y, f.y * t.x + f.x * t.y);
}

/**
 * atmturb_fresnel_tf_mul_batch_launch - Dispatch batched OTF complex product across frames
 */
static void atmturb_fresnel_tf_mul_batch_launch(
    cufftComplex       *d_field,
    const cufftComplex *d_tf,
    int                 npix,
    int                 count,
    cudaStream_t        stream)
{
    if ((npix & 1) == 0)
    {
        int  npairs = npix >> 1;
        dim3 b_2d(256, 1, 1);
        dim3 g_2d((unsigned int) ((npairs + 255) / 256), (unsigned int) count, 1);
        atmturb_fresnel_tf_mul_batch_vec2_kernel<<<g_2d, b_2d, 0, stream>>>(
            (float4 *) d_field, (const float4 *) d_tf, npairs, count);
    }
    else
    {
        dim3 b_2d(256, 1, 1);
        dim3 g_2d((unsigned int) ((npix + 255) / 256), (unsigned int) count, 1);
        atmturb_fresnel_tf_mul_batch_c2_kernel<<<g_2d, b_2d, 0, stream>>>(
            d_field, d_tf, npix, count);
    }
}

/**
 * atmturb_fresnel_batch_ctx_free - Release all GPU plans, transfer functions, and buffers
 */
static void atmturb_fresnel_batch_ctx_free(
    atmturb_fresnel_batch_ctx_t *c)
{
    if (c == NULL)
    {
        return;
    }
    if (c->stream != 0)
    {
        cudaStreamSynchronize(c->stream);
        cudaStreamDestroy(c->stream);
    }
    if (c->plan_many != 0)
    {
        cufftDestroy(c->plan_many);
    }
    if (c->plan_rem != 0)
    {
        cufftDestroy(c->plan_rem);
    }
    for (int m = 0; m < c->plan.nsuper; m++)
    {
        if (c->d_tf_pri != NULL && c->d_tf_pri[m] != NULL)
        {
            cudaFree(c->d_tf_pri[m]);
        }
        if (c->d_tf_sec != NULL && c->d_tf_sec[m] != NULL)
        {
            cudaFree(c->d_tf_sec[m]);
        }
    }
    free(c->d_tf_pri);
    free(c->d_tf_sec);
    cudaFreeHost(c->h_sub);
    cudaFree(c->d_sub);
    cudaFreeHost(c->h_lf);
    cudaFree(c->d_lf);
    cudaFree(c->d_fields);
    cudaFree(c->d_super_pha);
    cudaFree(c->d_super_spha);
    cudaFree(c->d_frame_pha);
    cudaFree(c->d_frame_amp);
    cudaFree(c->d_frame_spha);
    cudaFree(c->d_frame_samp);
    free(c->masters);
    free(c->screen_offsets);
    atmturb_fresnel_plan_free(&c->plan);
}

/**
 * atmturb_fresnel_batch_alloc_buffers - Allocate batched device arrays and cuFFT plans
 */
static int atmturb_fresnel_batch_alloc_buffers(
    atmturb_fresnel_batch_ctx_t *c)
{
    size_t npix       = (size_t) c->npix;
    size_t B          = (size_t) c->batch_b;
    size_t fbytes     = sizeof(float) * npix * B;
    size_t cpix_bytes = sizeof(cufftComplex) * npix;
    size_t sub_bytes  = sizeof(atmturb_cuda_rytov_sublayer_t) *
                        (size_t) c->plan.nsuper * B * (size_t) c->max_sublayers;

    if (cudaMallocHost((void **) &c->h_sub, sub_bytes) != cudaSuccess ||
        cudaMalloc((void **) &c->d_sub, sub_bytes) != cudaSuccess ||
        cudaMalloc((void **) &c->d_fields, cpix_bytes * B * (size_t) c->n_ch) != cudaSuccess ||
        cudaMalloc((void **) &c->d_super_pha, fbytes) != cudaSuccess ||
        cudaMalloc((void **) &c->d_frame_pha, fbytes) != cudaSuccess ||
        (c->has_amp && cudaMalloc((void **) &c->d_frame_amp, fbytes) != cudaSuccess))
    {
        return -1;
    }

    if (c->has_sec &&
        (cudaMalloc((void **) &c->d_super_spha, fbytes) != cudaSuccess ||
         cudaMalloc((void **) &c->d_frame_spha, fbytes) != cudaSuccess ||
         (c->has_amp && cudaMalloc((void **) &c->d_frame_samp, fbytes) != cudaSuccess)))
    {
        return -1;
    }

    if (c->lowfreq && c->max_lf_modes > 0)
    {
        size_t lf_bytes = sizeof(atmturb_cuda_lowfreq_mode_t) *
                          (size_t) c->plan.nsuper * B * (size_t) c->max_lf_modes;
        if (cudaMallocHost((void **) &c->h_lf, lf_bytes) != cudaSuccess ||
            cudaMalloc((void **) &c->d_lf, lf_bytes) != cudaSuccess)
        {
            return -1;
        }
    }

    int dims[2] = {(int) c->pup_size, (int) c->pup_size};
    int dist    = (int) c->npix;
    if (cufftPlanMany(&c->plan_many, 2, dims, dims, 1, dist, dims, 1, dist,
                      CUFFT_C2C, c->batch_b * c->n_ch) != CUFFT_SUCCESS ||
        cufftSetStream(c->plan_many, c->stream) != CUFFT_SUCCESS)
    {
        return -1;
    }
    if (c->batch_rem > 0 &&
        (cufftPlanMany(&c->plan_rem, 2, dims, dims, 1, dist, dims, 1, dist,
                       CUFFT_C2C, c->batch_rem * c->n_ch) != CUFFT_SUCCESS ||
         cufftSetStream(c->plan_rem, c->stream) != CUFFT_SUCCESS))
    {
        return -1;
    }
    return 0;
}

/**
 * atmturb_fresnel_batch_ctx_init - Initialize superlayer plan, OTFs, master cache, and buffers
 */
static int atmturb_fresnel_batch_ctx_init(
    atmturb_fresnel_batch_ctx_t            *c,
    const atmturb_cuda_fresnel_batch_req_t *req)
{
    memset(c, 0, sizeof(*c));
    c->msize    = req->master_size;
    c->pup_size = req->pup_size;
    c->npix     = req->pup_size * req->pup_size;
    c->has_sec  = (req->params->lambda_s_m > 0.0 && (req->spha != NULL || req->samp != NULL));
    c->has_amp  = (req->amp != NULL || req->samp != NULL) ? 1 : 0;
    c->n_ch     = c->has_sec ? 2 : 1;
    c->interp   = req->geom->interp;
    c->os       = (int) req->geom->oversample;
    c->lowfreq  = req->r->lowfreq;

    int max_b = (req->pup_size <= 256) ? 16 : ((req->pup_size <= 512) ? 8 : 4);
    c->batch_b   = (req->nbframes < (long) max_b) ? (int) req->nbframes : max_b;
    c->batch_rem = (int) (req->nbframes % (long) c->batch_b);

    double lam_s = c->has_sec ? req->params->lambda_s_m : 0.0;
    if (atmturb_fresnel_plan_init(&c->plan, req->prof, req->geom, req->pup_size,
                                  req->params->pupil_scale_m, req->params->lambda_ref_m,
                                  lam_s, (double) CONF_FRESNEL_PROPAGATION_BIN) != 0)
    {
        return -1;
    }

    c->max_sublayers = 1;
    for (int m = 0; m < c->plan.nsuper; m++)
    {
        int n_sub = 2 * c->plan.supers[m].nlayers;
        if (n_sub > c->max_sublayers)
        {
            c->max_sublayers = n_sub;
        }
        int n_lf = c->plan.supers[m].nlayers * 24;
        if (n_lf > c->max_lf_modes)
        {
            c->max_lf_modes = n_lf;
        }
    }

    if (atmturb_wfs_cuda_init_screen_pool(req->r, req->geom->nlayers, &c->masters,
                                          &c->screen_offsets, &c->total_screens) != 0)
    {
        return -1;
    }
    c->d_masters = atmturb_cuda_sync_device_masters(c->msize, c->total_screens, c->masters);
    if (c->d_masters == NULL || cudaStreamCreate(&c->stream) != cudaSuccess)
    {
        return -1;
    }

    size_t cpix_bytes = sizeof(cufftComplex) * (size_t) c->npix;
    c->d_tf_pri = (cufftComplex **) calloc((size_t) c->plan.nsuper, sizeof(cufftComplex *));
    c->d_tf_sec = c->has_sec
                      ? (cufftComplex **) calloc((size_t) c->plan.nsuper, sizeof(cufftComplex *))
                      : NULL;
    if (c->d_tf_pri == NULL || (c->has_sec && c->d_tf_sec == NULL))
    {
        return -1;
    }
    for (int m = 0; m < c->plan.nsuper; m++)
    {
        if (cudaMalloc((void **) &c->d_tf_pri[m], cpix_bytes) != cudaSuccess ||
            cudaMemcpyAsync(c->d_tf_pri[m], c->plan.tf_pri[m], cpix_bytes,
                            cudaMemcpyHostToDevice, c->stream) != cudaSuccess)
        {
            return -1;
        }
        if (c->has_sec &&
            (cudaMalloc((void **) &c->d_tf_sec[m], cpix_bytes) != cudaSuccess ||
             cudaMemcpyAsync(c->d_tf_sec[m], c->plan.tf_sec[m], cpix_bytes,
                             cudaMemcpyHostToDevice, c->stream) != cudaSuccess))
        {
            return -1;
        }
    }
    return atmturb_fresnel_batch_alloc_buffers(c);
}

/**
 * atmturb_fresnel_batch_pack_chunk - Build host sublayer & lowfreq tables for [t0 .. t0+count-1]
 */
static void atmturb_fresnel_batch_pack_chunk(
    atmturb_fresnel_batch_ctx_t            *c,
    const atmturb_cuda_fresnel_batch_req_t *req,
    long                                    t0,
    int                                     count)
{
    int nsuper  = c->plan.nsuper;
    int max_sub = c->max_sublayers;
    int max_lf  = c->max_lf_modes;

#ifdef _OPENMP
    #pragma omp parallel for if (count >= 4) schedule(static)
#endif
    for (int b = 0; b < count; b++)
    {
        atmturb_cuda_rytov_sublayer_t tmp_sub[256];
        atmturb_cuda_rytov_sublayer_t *sub_buf = tmp_sub;
        if (nsuper * max_sub > 256)
        {
            sub_buf = (atmturb_cuda_rytov_sublayer_t *) malloc(
                sizeof(*sub_buf) * (size_t) (nsuper * max_sub));
        }
        atmturb_wfs_cuda_build_frame_sublayers(
            req->r, req->geom, nsuper, c->plan.supers, c->screen_offsets,
            t0 + (long) b, req->params->time_step_s, 0.0, max_sub, sub_buf);
        for (int m = 0; m < nsuper; m++)
        {
            memcpy(&c->h_sub[(size_t) (m * c->batch_b + b) * (size_t) max_sub],
                   &sub_buf[(size_t) (m * max_sub)],
                   sizeof(*sub_buf) * (size_t) max_sub);
        }
        if (sub_buf != tmp_sub)
        {
            free(sub_buf);
        }

        if (c->lowfreq && max_lf > 0 && c->h_lf != NULL)
        {
            atmturb_cuda_lowfreq_mode_t *lf_buf = (atmturb_cuda_lowfreq_mode_t *) malloc(
                sizeof(*lf_buf) * (size_t) (nsuper * max_lf));
            if (lf_buf != NULL)
            {
                atmturb_wfs_cuda_build_frame_lf_modes(
                    req->r, req->geom, nsuper, c->plan.supers,
                    t0 + (long) b, req->params->time_step_s, 0.0, max_lf, lf_buf);
                for (int m = 0; m < nsuper; m++)
                {
                    memcpy(&c->h_lf[(size_t) (m * c->batch_b + b) * (size_t) max_lf],
                           &lf_buf[(size_t) (m * max_lf)],
                           sizeof(*lf_buf) * (size_t) max_lf);
                }
                free(lf_buf);
            }
        }
    }
}

/**
 * atmturb_fresnel_batch_prop_supers - Extrude, modulate, and propagate superlayers for a chunk
 */
static void atmturb_fresnel_batch_prop_supers(
    atmturb_fresnel_batch_ctx_t *c,
    int                          count,
    int                          tot_pix,
    cufftComplex                *d_pri,
    cufftComplex                *d_sec,
    cufftHandle                  fft_plan)
{
    for (int m = 0; m < c->plan.nsuper; m++)
    {
        int n_sub = 2 * c->plan.supers[m].nlayers;
        int n_lf  = (c->lowfreq && c->d_lf != NULL) ? (c->plan.supers[m].nlayers * 24) : 0;
        atmturb_cuda_extrude_batch_args_t ext;
        ext.interp_mode  = c->interp;
        ext.d_masters    = c->d_masters;
        ext.msize        = (int) c->msize;
        ext.pad_size     = (int) c->pup_size;
        ext.pup_size     = (int) c->pup_size;
        ext.guard        = 0;
        ext.n_sub        = n_sub;
        ext.sub_stride   = c->max_sublayers;
        ext.has_sec      = c->has_sec;
        ext.os           = c->os;
        ext.batch_count  = count;
        ext.d_sublayers  = c->d_sub + (size_t) (m * c->batch_b * c->max_sublayers);
        ext.n_lf_modes   = n_lf;
        ext.lf_stride    = c->max_lf_modes;
        ext.d_lf_modes   = (n_lf > 0) ? (c->d_lf + (size_t) (m * c->batch_b * c->max_lf_modes))
                                      : NULL;
        ext.d_super_pha  = c->d_super_pha;
        ext.d_super_spha = c->d_super_spha;
        ext.d_frame_pha  = c->d_frame_pha;
        ext.d_frame_spha = c->d_frame_spha;
        atmturb_cuda_extrude_sl_batch_launch(&ext, c->stream);

        if (m == 0)
        {
            atmturb_cuda_fresnel_init_modulate_launch(d_pri, c->d_super_pha, tot_pix, c->stream);
            if (c->has_sec)
            {
                atmturb_cuda_fresnel_init_modulate_launch(d_sec, c->d_super_spha, tot_pix,
                                                          c->stream);
            }
        }
        else
        {
            atmturb_cuda_fresnel_modulate_launch(d_pri, c->d_super_pha, tot_pix, c->stream);
            if (c->has_sec)
            {
                atmturb_cuda_fresnel_modulate_launch(d_sec, c->d_super_spha, tot_pix, c->stream);
            }
        }

        if (fabs(c->plan.supers[m].step_dist_m) > 0.0)
        {
            cufftExecC2C(fft_plan, d_pri, d_pri, CUFFT_FORWARD);
            atmturb_fresnel_tf_mul_batch_launch(d_pri, c->d_tf_pri[m], (int) c->npix, count,
                                                c->stream);
            if (c->has_sec)
            {
                atmturb_fresnel_tf_mul_batch_launch(d_sec, c->d_tf_sec[m], (int) c->npix, count,
                                                    c->stream);
            }
            cufftExecC2C(fft_plan, d_pri, d_pri, CUFFT_INVERSE);
        }
    }
}

/**
 * atmturb_fresnel_batch_run_chunk - Enqueue and copy one multi-frame Fresnel batch chunk
 */
static int atmturb_fresnel_batch_run_chunk(
    atmturb_fresnel_batch_ctx_t            *c,
    const atmturb_cuda_fresnel_batch_req_t *req,
    long                                    t0,
    int                                     count)
{
    atmturb_fresnel_batch_pack_chunk(c, req, t0, count);

    size_t sub_bytes = sizeof(atmturb_cuda_rytov_sublayer_t) *
                       (size_t) c->plan.nsuper * (size_t) c->batch_b * (size_t) c->max_sublayers;
    cudaMemcpyAsync(c->d_sub, c->h_sub, sub_bytes, cudaMemcpyHostToDevice, c->stream);
    if (c->lowfreq && c->max_lf_modes > 0 && c->d_lf != NULL)
    {
        size_t lf_bytes = sizeof(atmturb_cuda_lowfreq_mode_t) *
                          (size_t) c->plan.nsuper * (size_t) c->batch_b * (size_t) c->max_lf_modes;
        cudaMemcpyAsync(c->d_lf, c->h_lf, lf_bytes, cudaMemcpyHostToDevice, c->stream);
    }

    int           tot_pix    = count * (int) c->npix;
    size_t        chunk_fbyt = sizeof(float) * (size_t) tot_pix;
    cufftComplex *d_pri      = c->d_fields;
    cufftComplex *d_sec      = c->has_sec ? (c->d_fields + (size_t) tot_pix) : NULL;
    cufftHandle   fft_plan   = (count == c->batch_b) ? c->plan_many : c->plan_rem;

    cudaMemsetAsync(c->d_frame_pha, 0, chunk_fbyt, c->stream);
    if (c->has_sec)
    {
        cudaMemsetAsync(c->d_frame_spha, 0, chunk_fbyt, c->stream);
    }

    atmturb_fresnel_batch_prop_supers(c, count, tot_pix, d_pri, d_sec, fft_plan);

    atmturb_cuda_fresnel_extract_launch(
        d_pri, c->d_frame_pha, c->has_amp ? c->d_frame_amp : NULL, NULL, tot_pix, c->stream);
    if (c->has_sec)
    {
        atmturb_cuda_fresnel_extract_launch(
            d_sec, c->d_frame_spha, c->has_amp ? c->d_frame_samp : NULL, NULL, tot_pix, c->stream);
    }

    size_t dst_off = (size_t) t0 * (size_t) c->npix;
    if (req->pha != NULL &&
        cudaMemcpyAsync(req->pha + dst_off, c->d_frame_pha, chunk_fbyt,
                        cudaMemcpyDeviceToHost, c->stream) != cudaSuccess)
    {
        return -1;
    }
    if (req->amp != NULL &&
        cudaMemcpyAsync(req->amp + dst_off, c->d_frame_amp, chunk_fbyt,
                        cudaMemcpyDeviceToHost, c->stream) != cudaSuccess)
    {
        return -1;
    }
    if (c->has_sec && req->spha != NULL &&
        cudaMemcpyAsync(req->spha + dst_off, c->d_frame_spha, chunk_fbyt,
                        cudaMemcpyDeviceToHost, c->stream) != cudaSuccess)
    {
        return -1;
    }
    if (c->has_sec && req->samp != NULL &&
        cudaMemcpyAsync(req->samp + dst_off, c->d_frame_samp, chunk_fbyt,
                        cudaMemcpyDeviceToHost, c->stream) != cudaSuccess)
    {
        return -1;
    }
    return (cudaStreamSynchronize(c->stream) == cudaSuccess) ? 0 : -1;
}

int atmturb_cuda_fresnel_render_batch(
    const atmturb_cuda_fresnel_batch_req_t *req)
{
    if (req == NULL || req->r == NULL || req->geom == NULL || req->prof == NULL ||
        req->params == NULL || req->nbframes <= 0 || req->master_size <= 0 ||
        req->pup_size <= 0 || !atmturb_cuda_device_available())
    {
        return -1;
    }

    atmturb_fresnel_batch_ctx_t ctx;
    if (atmturb_fresnel_batch_ctx_init(&ctx, req) != 0)
    {
        atmturb_fresnel_batch_ctx_free(&ctx);
        return -1;
    }

    int ret = 0;
    for (long t0 = 0; t0 < req->nbframes; t0 += (long) ctx.batch_b)
    {
        long rem   = req->nbframes - t0;
        int  count = (rem < (long) ctx.batch_b) ? (int) rem : ctx.batch_b;
        if (atmturb_fresnel_batch_run_chunk(&ctx, req, t0, count) != 0)
        {
            ret = -1;
            break;
        }
    }

    atmturb_fresnel_batch_ctx_free(&ctx);
    return ret;
}
