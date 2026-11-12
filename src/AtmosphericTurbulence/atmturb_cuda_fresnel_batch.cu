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
#include "atmturb_cuda_dma_ring.h"
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
    cudaStream_t stream;
    cufftHandle  plan_many;
    int          many_pup;
    int          many_cnt;
    cufftHandle  plan_rem;
    int          rem_pup;
    int          rem_cnt;
    char        *h_buf;
    size_t       h_cap;
    char        *d_buf;
    size_t       d_cap;
} fresnel_batch_cache_t;

static fresnel_batch_cache_t s_fbat = {0, 0, 0, 0, 0, 0, 0, NULL, 0, NULL, 0};

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
    atmturb_cuda_rytov_sublayer_t *h_sub[2];
    atmturb_cuda_rytov_sublayer_t *d_sub[2];
    atmturb_cuda_lowfreq_mode_t   *h_lf[2];
    atmturb_cuda_lowfreq_mode_t   *d_lf[2];
    cufftComplex                  *d_fields;
    float                         *d_super_pha;
    float                         *d_super_spha;
    float                         *d_frame_pha[2];
    float                         *d_frame_amp[2];
    float                         *d_frame_spha[2];
    float                         *d_frame_samp[2];
} atmturb_fresnel_batch_ctx_t;

static inline size_t fb_align256(
    size_t n)
{
    return (n + 255) & ~(size_t) 255;
}

void atmturb_cuda_fresnel_batch_cleanup(void)
{
    if (s_fbat.stream != 0)
    {
        cudaStreamSynchronize(s_fbat.stream);
        cudaStreamDestroy(s_fbat.stream);
        s_fbat.stream = 0;
    }
    if (s_fbat.plan_many != 0)
    {
        cufftDestroy(s_fbat.plan_many);
        s_fbat.plan_many = 0;
    }
    if (s_fbat.plan_rem != 0)
    {
        cufftDestroy(s_fbat.plan_rem);
        s_fbat.plan_rem = 0;
    }
    cudaFreeHost(s_fbat.h_buf);
    cudaFree(s_fbat.d_buf);
    memset(&s_fbat, 0, sizeof(s_fbat));
}

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
    out.x        = f.x * t.x - f.y * t.y;
    out.y        = f.y * t.x + f.x * t.y;
    out.z        = f.z * t.z - f.w * t.w;
    out.w        = f.w * t.z + f.z * t.w;
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
 * atmturb_fresnel_batch_ctx_free - Release per-call host descriptors while preserving GPU cache
 */
static void atmturb_fresnel_batch_ctx_free(
    atmturb_fresnel_batch_ctx_t *c)
{
    if (c == NULL)
    {
        return;
    }
    free(c->d_tf_pri);
    free(c->d_tf_sec);
    free(c->masters);
    free(c->screen_offsets);
    atmturb_fresnel_plan_free(&c->plan);
}

/**
 * atmturb_fresnel_batch_ensure_plans - Reuse or create cached batched cuFFT plans
 */
static int atmturb_fresnel_batch_ensure_plans(
    atmturb_fresnel_batch_ctx_t *c)
{
    int dims[2]  = {(int) c->pup_size, (int) c->pup_size};
    int dist     = (int) c->npix;
    int req_many = c->batch_b * c->n_ch;
    int req_rem  = c->batch_rem * c->n_ch;

    if (s_fbat.plan_many == 0 || s_fbat.many_pup != (int) c->pup_size ||
        s_fbat.many_cnt != req_many)
    {
        if (s_fbat.plan_many != 0)
        {
            cufftDestroy(s_fbat.plan_many);
            s_fbat.plan_many = 0;
        }
        if (cufftPlanMany(&s_fbat.plan_many, 2, dims, dims, 1, dist, dims, 1, dist,
                          CUFFT_C2C, req_many) != CUFFT_SUCCESS ||
            cufftSetStream(s_fbat.plan_many, c->stream) != CUFFT_SUCCESS)
        {
            return -1;
        }
        s_fbat.many_pup = (int) c->pup_size;
        s_fbat.many_cnt = req_many;
    }
    c->plan_many = s_fbat.plan_many;

    if (req_rem > 0)
    {
        if (s_fbat.plan_rem == 0 || s_fbat.rem_pup != (int) c->pup_size ||
            s_fbat.rem_cnt != req_rem)
        {
            if (s_fbat.plan_rem != 0)
            {
                cufftDestroy(s_fbat.plan_rem);
                s_fbat.plan_rem = 0;
            }
            if (cufftPlanMany(&s_fbat.plan_rem, 2, dims, dims, 1, dist, dims, 1, dist,
                              CUFFT_C2C, req_rem) != CUFFT_SUCCESS ||
                cufftSetStream(s_fbat.plan_rem, c->stream) != CUFFT_SUCCESS)
            {
                return -1;
            }
            s_fbat.rem_pup = (int) c->pup_size;
            s_fbat.rem_cnt = req_rem;
        }
        c->plan_rem = s_fbat.plan_rem;
    }
    return 0;
}

/**
 * atmturb_fresnel_batch_slice_workspace - Partition cached host and device arenas into slots
 */
static void atmturb_fresnel_batch_slice_workspace(
    atmturb_fresnel_batch_ctx_t *c,
    size_t                       fld_bytes,
    size_t                       fbytes,
    size_t                       sub_bytes,
    size_t                       lf_bytes,
    size_t                       tf_bytes)
{
    char *hp = s_fbat.h_buf;
    char *dp = s_fbat.d_buf;

    c->d_fields = (cufftComplex *) dp;
    dp         += fld_bytes;
    c->d_super_pha = (float *) dp;
    dp            += fbytes;
    c->d_super_spha = c->has_sec ? (float *) dp : NULL;
    if (c->has_sec)
    {
        dp += fbytes;
    }

    for (int s = 0; s < 2; s++)
    {
        c->h_sub[s] = (atmturb_cuda_rytov_sublayer_t *) hp;
        hp         += sub_bytes;
        c->h_lf[s]  = (lf_bytes > 0) ? (atmturb_cuda_lowfreq_mode_t *) hp : NULL;
        hp         += lf_bytes;

        c->d_sub[s] = (atmturb_cuda_rytov_sublayer_t *) dp;
        dp         += sub_bytes;
        c->d_lf[s]  = (lf_bytes > 0) ? (atmturb_cuda_lowfreq_mode_t *) dp : NULL;
        dp         += lf_bytes;

        c->d_frame_pha[s] = (float *) dp;
        dp               += fbytes;
        c->d_frame_amp[s] = c->has_amp ? (float *) dp : NULL;
        if (c->has_amp) { dp += fbytes; }
        c->d_frame_spha[s] = c->has_sec ? (float *) dp : NULL;
        if (c->has_sec) { dp += fbytes; }
        c->d_frame_samp[s] = (c->has_sec && c->has_amp) ? (float *) dp : NULL;
        if (c->has_sec && c->has_amp) { dp += fbytes; }
    }

    for (int m = 0; m < c->plan.nsuper; m++)
    {
        c->d_tf_pri[m] = (cufftComplex *) dp;
        dp            += tf_bytes;
        if (c->has_sec)
        {
            c->d_tf_sec[m] = (cufftComplex *) dp;
            dp            += tf_bytes;
        }
    }
}

/**
 * atmturb_fresnel_batch_alloc_buffers - Ensure cached device/pinned workspace and upload OTFs
 */
static int atmturb_fresnel_batch_alloc_buffers(
    atmturb_fresnel_batch_ctx_t *c)
{
    size_t npix      = (size_t) c->npix;
    size_t B         = (size_t) c->batch_b;
    size_t raw_fbyt  = sizeof(float) * npix * B;
    size_t fbytes    = fb_align256(raw_fbyt);
    size_t fld_bytes = fb_align256(sizeof(cufftComplex) * npix * B * (size_t) c->n_ch);
    size_t tf_raw    = sizeof(cufftComplex) * npix;
    size_t tf_bytes  = fb_align256(tf_raw);
    size_t sub_bytes = fb_align256(sizeof(atmturb_cuda_rytov_sublayer_t) *
                                   (size_t) c->plan.nsuper * B * (size_t) c->max_sublayers);
    size_t lf_bytes  = (c->lowfreq && c->max_lf_modes > 0)
                           ? fb_align256(sizeof(atmturb_cuda_lowfreq_mode_t) *
                                         (size_t) c->plan.nsuper * B * (size_t) c->max_lf_modes)
                           : 0;
    int    n_out_ch  = (1 + c->has_amp) * c->n_ch;
    size_t h_need    = 2 * (sub_bytes + lf_bytes);
    size_t d_need    = fld_bytes + (size_t) c->n_ch * fbytes +
                       2 * (sub_bytes + lf_bytes + (size_t) n_out_ch * fbytes) +
                       (size_t) (c->plan.nsuper * c->n_ch) * tf_bytes;

    if (atmturb_cuda_dma_ring_prepare(raw_fbyt) != 0)
    {
        return -1;
    }
    if (s_fbat.h_cap < h_need)
    {
        cudaFreeHost(s_fbat.h_buf);
        s_fbat.h_buf = NULL;
        if (cudaMallocHost((void **) &s_fbat.h_buf, h_need) != cudaSuccess)
        {
            s_fbat.h_cap = 0;
            return -1;
        }
        s_fbat.h_cap = h_need;
    }
    if (s_fbat.d_cap < d_need)
    {
        cudaFree(s_fbat.d_buf);
        s_fbat.d_buf = NULL;
        if (cudaMalloc((void **) &s_fbat.d_buf, d_need) != cudaSuccess)
        {
            s_fbat.d_cap = 0;
            return -1;
        }
        s_fbat.d_cap = d_need;
    }

    atmturb_fresnel_batch_slice_workspace(c, fld_bytes, fbytes, sub_bytes, lf_bytes, tf_bytes);
    for (int m = 0; m < c->plan.nsuper; m++)
    {
        if (fabs(c->plan.supers[m].step_dist_m) > 0.0)
        {
            cudaMemcpyAsync(c->d_tf_pri[m], c->plan.tf_pri[m], tf_raw,
                            cudaMemcpyHostToDevice, c->stream);
            if (c->has_sec)
            {
                cudaMemcpyAsync(c->d_tf_sec[m], c->plan.tf_sec[m], tf_raw,
                                cudaMemcpyHostToDevice, c->stream);
            }
        }
    }
    return atmturb_fresnel_batch_ensure_plans(c);
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

    int max_b    = (req->pup_size <= 256) ? 32 : ((req->pup_size <= 512) ? 16 : 8);
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
    c->max_lf_modes  = c->lowfreq ? ATMTURB_LOWFREQ_NMODES : 0;
    for (int m = 0; m < c->plan.nsuper; m++)
    {
        int n_sub = 2 * c->plan.supers[m].nlayers;
        if (n_sub > c->max_sublayers)
        {
            c->max_sublayers = n_sub;
        }
    }

    if (atmturb_wfs_cuda_init_screen_pool(req->r, req->geom->nlayers, &c->masters,
                                          &c->screen_offsets, &c->total_screens) != 0)
    {
        return -1;
    }
    c->d_masters = atmturb_cuda_sync_device_masters(c->msize, c->total_screens, c->masters);
    if (c->d_masters == NULL)
    {
        return -1;
    }
    if (s_fbat.stream == 0 &&
        cudaStreamCreateWithFlags(&s_fbat.stream, cudaStreamNonBlocking) != cudaSuccess)
    {
        return -1;
    }
    c->stream = s_fbat.stream;

    c->d_tf_pri = (cufftComplex **) calloc((size_t) c->plan.nsuper, sizeof(cufftComplex *));
    c->d_tf_sec = c->has_sec
                      ? (cufftComplex **) calloc((size_t) c->plan.nsuper, sizeof(cufftComplex *))
                      : NULL;
    if (c->d_tf_pri == NULL || (c->has_sec && c->d_tf_sec == NULL))
    {
        return -1;
    }
    return atmturb_fresnel_batch_alloc_buffers(c);
}

/**
 * atmturb_fresnel_batch_pack_chunk - Build host sublayer & lowfreq tables for [t0 .. t0+count-1]
 */
static void atmturb_fresnel_batch_pack_chunk(
    atmturb_fresnel_batch_ctx_t            *c,
    const atmturb_cuda_fresnel_batch_req_t *req,
    int                                     slot,
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
        atmturb_cuda_rytov_sublayer_t  tmp_sub[256];
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
            memcpy(&c->h_sub[slot][(size_t) (m * c->batch_b + b) * (size_t) max_sub],
                   &sub_buf[(size_t) (m * max_sub)], sizeof(*sub_buf) * (size_t) max_sub);
        }
        if (sub_buf != tmp_sub)
        {
            free(sub_buf);
        }

        if (c->lowfreq && max_lf > 0 && c->h_lf[slot] != NULL)
        {
            atmturb_cuda_lowfreq_mode_t  tmp_lf[256];
            atmturb_cuda_lowfreq_mode_t *lf_buf = tmp_lf;
            if (nsuper * max_lf > 256)
            {
                lf_buf = (atmturb_cuda_lowfreq_mode_t *) malloc(
                    sizeof(*lf_buf) * (size_t) (nsuper * max_lf));
            }
            atmturb_wfs_cuda_build_frame_lf_modes(
                req->r, req->geom, nsuper, c->plan.supers,
                t0 + (long) b, req->params->time_step_s, 0.0, max_lf, lf_buf);
            for (int m = 0; m < nsuper; m++)
            {
                memcpy(&c->h_lf[slot][(size_t) (m * c->batch_b + b) * (size_t) max_lf],
                       &lf_buf[(size_t) (m * max_lf)], sizeof(*lf_buf) * (size_t) max_lf);
            }
            if (lf_buf != tmp_lf)
            {
                free(lf_buf);
            }
        }
    }
}

/**
 * atmturb_fresnel_batch_extrude_one - Extrude one superlayer @m for @count frames in @slot
 */
static void atmturb_fresnel_batch_extrude_one(
    atmturb_fresnel_batch_ctx_t *c,
    int                          slot,
    int                          m,
    int                          count,
    int                          direct_overwrite,
    int                          accum_only)
{
    int n_sub = 2 * c->plan.supers[m].nlayers;
    int n_lf  = (c->lowfreq && c->d_lf[slot] != NULL && c->plan.supers[m].nlayers > 0)
                    ? c->max_lf_modes : 0;

    atmturb_cuda_extrude_batch_args_t ext = {
        .interp_mode  = c->interp,
        .d_masters    = c->d_masters,
        .msize        = (int) c->msize,
        .pad_size     = (int) c->pup_size,
        .pup_size     = (int) c->pup_size,
        .guard        = 0,
        .n_sub        = n_sub,
        .sub_stride   = c->max_sublayers,
        .has_sec      = c->has_sec,
        .os           = c->os,
        .batch_count  = count,
        .d_sublayers  = c->d_sub[slot] + (size_t) (m * c->batch_b * c->max_sublayers),
        .n_lf_modes   = n_lf,
        .lf_stride    = c->max_lf_modes,
        .d_lf_modes   = (n_lf > 0) ? (c->d_lf[slot] + (size_t) (m * c->batch_b * c->max_lf_modes))
                                   : NULL,
        .d_super_pha  = direct_overwrite ? c->d_frame_pha[slot]
                                         : (accum_only ? NULL : c->d_super_pha),
        .d_super_spha = direct_overwrite ? c->d_frame_spha[slot]
                                         : (accum_only ? NULL : c->d_super_spha),
        .d_frame_pha  = direct_overwrite ? NULL : c->d_frame_pha[slot],
        .d_frame_spha = direct_overwrite ? NULL : c->d_frame_spha[slot]
    };
    atmturb_cuda_extrude_sl_batch_launch(&ext, c->stream);
}

/**
 * atmturb_fresnel_batch_prop_supers - Extrude, modulate, and propagate superlayers for a chunk
 */
static void atmturb_fresnel_batch_prop_supers(
    atmturb_fresnel_batch_ctx_t *c,
    int                          slot,
    int                          n_prop,
    int                          count,
    int                          tot_pix,
    cufftComplex                *d_pri,
    cufftComplex                *d_sec,
    cufftHandle                  fft_plan)
{
    for (int m = 0; m < n_prop; m++)
    {
        atmturb_fresnel_batch_extrude_one(c, slot, m, count, (m == 0), 0);

        if (m == 0)
        {
            atmturb_cuda_fresnel_init_modulate_launch(d_pri, c->d_frame_pha[slot], tot_pix,
                                                      c->stream);
            if (c->has_sec)
            {
                atmturb_cuda_fresnel_init_modulate_launch(d_sec, c->d_frame_spha[slot], tot_pix,
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
 * atmturb_fresnel_batch_run_chunk - Enqueue GPU Fresnel chunk and hand off to async D2H ring
 */
static int atmturb_fresnel_batch_run_chunk(
    atmturb_fresnel_batch_ctx_t            *c,
    const atmturb_cuda_fresnel_batch_req_t *req,
    int                                     slot,
    long                                    t0,
    int                                     count)
{
    if (atmturb_cuda_dma_ring_drain(slot) != 0)
    {
        return -1;
    }

    atmturb_fresnel_batch_pack_chunk(c, req, slot, t0, count);

    size_t sub_bytes = sizeof(atmturb_cuda_rytov_sublayer_t) *
                       (size_t) c->plan.nsuper * (size_t) c->batch_b * (size_t) c->max_sublayers;
    cudaMemcpyAsync(c->d_sub[slot], c->h_sub[slot], sub_bytes, cudaMemcpyHostToDevice, c->stream);
    if (c->lowfreq && c->max_lf_modes > 0 && c->d_lf[slot] != NULL)
    {
        size_t lf_bytes = sizeof(atmturb_cuda_lowfreq_mode_t) *
                          (size_t) c->plan.nsuper * (size_t) c->batch_b * (size_t) c->max_lf_modes;
        cudaMemcpyAsync(c->d_lf[slot], c->h_lf[slot], lf_bytes, cudaMemcpyHostToDevice, c->stream);
    }

    int           tot_pix        = count * (int) c->npix;
    size_t        chunk_fbyt     = sizeof(float) * (size_t) tot_pix;
    cufftComplex *d_pri          = c->d_fields;
    cufftComplex *d_sec          = c->has_sec ? (c->d_fields + (size_t) tot_pix) : NULL;
    cufftHandle   fft_plan       = (count == c->batch_b) ? c->plan_many : c->plan_rem;
    int           has_trail_zero = (c->plan.nsuper > 0 &&
                                    fabs(c->plan.supers[c->plan.nsuper - 1].step_dist_m) <= 0.0);
    int           n_prop         = has_trail_zero ? (c->plan.nsuper - 1) : c->plan.nsuper;

    if (n_prop <= 0)
    {
        if (has_trail_zero)
        {
            atmturb_fresnel_batch_extrude_one(c, slot, 0, count, 1, 0);
        }
        else
        {
            cudaMemsetAsync(c->d_frame_pha[slot], 0, chunk_fbyt, c->stream);
            if (c->has_sec)
            {
                cudaMemsetAsync(c->d_frame_spha[slot], 0, chunk_fbyt, c->stream);
            }
        }
        if (c->has_amp)
        {
            atmturb_cuda_fill_ones_launch(c->d_frame_amp[slot], tot_pix, c->stream);
            if (c->has_sec)
            {
                atmturb_cuda_fill_ones_launch(c->d_frame_samp[slot], tot_pix, c->stream);
            }
        }
    }
    else
    {
        atmturb_fresnel_batch_prop_supers(c, slot, n_prop, count, tot_pix, d_pri, d_sec, fft_plan);
        atmturb_cuda_fresnel_extract_launch(
            d_pri, c->d_frame_pha[slot], c->has_amp ? c->d_frame_amp[slot] : NULL,
            NULL, tot_pix, c->stream);
        if (c->has_sec)
        {
            atmturb_cuda_fresnel_extract_launch(
                d_sec, c->d_frame_spha[slot], c->has_amp ? c->d_frame_samp[slot] : NULL,
                NULL, tot_pix, c->stream);
        }
        if (has_trail_zero)
        {
            atmturb_fresnel_batch_extrude_one(c, slot, c->plan.nsuper - 1, count, 0, 1);
        }
    }

    size_t             dst_off  = (size_t) t0 * (size_t) c->npix;
    const float *const d_src[4] = {
        c->d_frame_pha[slot], c->d_frame_amp[slot],
        c->d_frame_spha[slot], c->d_frame_samp[slot]
    };
    float *const h_dst[4] = {
        (req->pha != NULL) ? (req->pha + dst_off) : NULL,
        (req->amp != NULL) ? (req->amp + dst_off) : NULL,
        (c->has_sec && req->spha != NULL) ? (req->spha + dst_off) : NULL,
        (c->has_sec && req->samp != NULL) ? (req->samp + dst_off) : NULL
    };
    return atmturb_cuda_dma_ring_enqueue(slot, c->stream, d_src, h_dst, chunk_fbyt);
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

    int ret   = 0;
    int chunk = 0;
    for (long t0 = 0; t0 < req->nbframes; t0 += (long) ctx.batch_b, chunk++)
    {
        long rem   = req->nbframes - t0;
        int  count = (rem < (long) ctx.batch_b) ? (int) rem : ctx.batch_b;
        if (atmturb_fresnel_batch_run_chunk(&ctx, req, chunk & 1, t0, count) != 0)
        {
            ret = -1;
            break;
        }
    }

    if (atmturb_cuda_dma_ring_flush() != 0)
    {
        ret = -1;
    }
    atmturb_fresnel_batch_ctx_free(&ctx);
    return ret;
}
