// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_rytov_batch.cu
 * @brief   Multi-frame batched cuFFT 3D wavefront cube renderer for Rytov propagation
 */

#include <cuda_runtime.h>
#include <cufft.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include "atmturb_cuda.h"
#include "atmturb_cuda_dma_ring.h"
#include "atmturb_cuda_rytov.h"
#include "atmturb_cuda_rytov_internal.h"
}

typedef struct
{
    cufftComplex       *d_acc_batch;
    const cufftComplex *d_spec;
    const cufftComplex *d_hat_bounds;
    const cufftComplex *d_ramp;
    const cufftComplex *d_exp_x;
    const cufftComplex *d_exp_y;
    const float        *d_laplace_inv;
    const float        *d_filt_a_pri;
    const float        *d_filt_b_pri;
    const float        *d_filt_a_sec;
    const float        *d_filt_b_sec;
    int                 pad;
    int                 n_half;
    int                 n_inv;
    int                 ch_base;
    int                 count;
} rytov_batch_mac_args_t;

typedef struct
{
    long                           pad;
    long                           pup;
    int                            has_sec;
    int                            nsuper;
    int                            batch_b;
    int                            batch_rem;
    size_t                         sub_cap;
    size_t                         lf_cap;
    cudaStream_t                   stream;
    cufftHandle                    plan_1d_many;
    cufftHandle                    plan_r2c_many;
    cufftHandle                    plan_c2r_many;
    cufftHandle                    plan_1d_rem;
    cufftHandle                    plan_r2c_rem;
    cufftHandle                    plan_c2r_rem;
    float                         *d_laplace_inv;
    cufftComplex                  *d_exp_x;
    cufftComplex                  *d_exp_y;
    float                         *d_filt_a_pri;
    float                         *d_filt_b_pri;
    float                         *d_filt_a_sec;
    float                         *d_filt_b_sec;
    cufftComplex                  *d_chrom_ramp;
    float                         *d_super_pha;
    float                         *d_super_spha;
    float                         *d_bounds;
    cufftComplex                  *d_hat_bounds;
    cufftComplex                  *d_spec;
    cufftComplex                  *d_acc_batch;
    float                         *d_inv_out;
    float                         *d_sum_I;
    atmturb_cuda_rytov_sublayer_t *h_sub[2];
    atmturb_cuda_rytov_sublayer_t *d_sub[2];
    atmturb_cuda_lowfreq_mode_t   *h_lf[2];
    atmturb_cuda_lowfreq_mode_t   *d_lf[2];
    float                         *d_frame_pha[2];
    float                         *d_frame_amp[2];
    float                         *d_frame_spha[2];
    float                         *d_frame_samp[2];
} rytov_batch_state_t;

static rytov_batch_state_t s_bat;

/**
 * rytov_batch_bounds_kernel - Extract top-bottom and left-right boundary jumps across batch
 */
__global__ static void rytov_batch_bounds_kernel(
    const float *d_super_pha,
    float       *d_bounds,
    int          pad,
    int          count)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int b   = (int) blockIdx.y;
    if (idx >= pad || b >= count)
    {
        return;
    }

    const float *src = d_super_pha + (size_t) b * (size_t) pad * (size_t) pad;
    float       *dst = d_bounds + (size_t) b * (size_t) (2 * pad);
    dst[idx]         = src[idx] - src[(pad - 1) * pad + idx];
    dst[pad + idx]   = src[idx * pad] - src[idx * pad + (pad - 1)];
}

/**
 * rytov_batch_moisan_val - Evaluate Poisson periodic boundary adjustment at (i, j)
 */
__device__ static inline cufftComplex rytov_batch_moisan_val(
    cufftComplex        s,
    const cufftComplex *hat_bounds,
    const cufftComplex *exp_x,
    const cufftComplex *exp_y,
    const float        *laplace_inv,
    int                 i,
    int                 j,
    int                 pad,
    int                 n_half)
{
    const cufftComplex *hat_a = hat_bounds;
    const cufftComplex *hat_b = hat_bounds + n_half;
    float               b_re  = (j <= pad / 2) ? hat_b[j].x : hat_b[pad - j].x;
    float               b_im  = (j <= pad / 2) ? hat_b[j].y : -hat_b[pad - j].y;
    cufftComplex        ey    = exp_y[j];
    cufftComplex        ex    = exp_x[i];
    cufftComplex        ha    = hat_a[i];

    float v_re = (ha.x * ey.x - ha.y * ey.y) + (b_re * ex.x - b_im * ex.y);
    float v_im = (ha.x * ey.y + ha.y * ey.x) + (b_re * ex.y + b_im * ex.x);
    float linv = laplace_inv[j * n_half + i];
    s.x += v_re * linv;
    s.y += v_im * linv;
    return s;
}

/**
 * rytov_batch_spec_mac_kernel - Fused batched Moisan subtraction and Rytov filter MAC
 */
template <int USE_MOISAN, int IS_DUAL>
__global__ static void rytov_batch_spec_mac_kernel(
    rytov_batch_mac_args_t args)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int b = (int) blockIdx.z;
    if (i >= args.n_half || j >= args.pad || b >= args.count)
    {
        return;
    }

    int          idx    = j * args.n_half + i;
    size_t       n_spec = (size_t) args.pad * (size_t) args.n_half;
    cufftComplex s      = args.d_spec[(size_t) b * n_spec + (size_t) idx];

    if (USE_MOISAN)
    {
        s = rytov_batch_moisan_val(
            s, args.d_hat_bounds + (size_t) b * (size_t) (2 * args.n_half),
            args.d_exp_x, args.d_exp_y, args.d_laplace_inv, i, j, args.pad, args.n_half);
    }

    cufftComplex *acc = args.d_acc_batch + (size_t) (b * args.n_inv + args.ch_base) * n_spec;
    float         fa  = __ldg(&args.d_filt_a_pri[idx]);
    float         fb  = __ldg(&args.d_filt_b_pri[idx]);

    acc[idx].x          += s.x * fa;
    acc[idx].y          += s.y * fa;
    acc[n_spec + idx].x += s.x * fb;
    acc[n_spec + idx].y += s.y * fb;

    if (IS_DUAL)
    {
        cufftComplex r    = __ldg(&args.d_ramp[idx]);
        float        ps_x = s.x * r.x - s.y * r.y;
        float        ps_y = s.x * r.y + s.y * r.x;
        float        fa_s = __ldg(&args.d_filt_a_sec[idx]);
        float        fb_s = __ldg(&args.d_filt_b_sec[idx]);

        acc[2 * n_spec + idx].x += ps_x * fa_s;
        acc[2 * n_spec + idx].y += ps_y * fa_s;
        acc[3 * n_spec + idx].x += ps_x * fb_s;
        acc[3 * n_spec + idx].y += ps_y * fb_s;
    }
}

/**
 * rytov_batch_assemble_kernel - Assemble pupil phase & amplitude and warp-reduce sum(A^2)
 */
__global__ static void rytov_batch_assemble_kernel(
    int          pup_size,
    int          guard,
    int          pad_size,
    int          n_inv,
    int          ch_base,
    float       *d_frame_pha,
    float       *d_frame_amp,
    const float *d_inv_out,
    float       *d_sum_I)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int b = (int) blockIdx.z;

    float a_sq = 0.0f;
    if (x < pup_size && y < pup_size)
    {
        int max_clamp = (pup_size > 2) ? (pup_size - 2) : 0;
        int sx        = (guard == 0 && (x == 0 || x == pup_size - 1))
                            ? ((x == 0) ? ((pup_size > 1) ? 1 : 0) : max_clamp) : x;
        int sy        = (guard == 0 && (y == 0 || y == pup_size - 1))
                            ? ((y == 0) ? ((pup_size > 1) ? 1 : 0) : max_clamp) : y;

        size_t       pad_sq = (size_t) pad_size * (size_t) pad_size;
        size_t       pup_sq = (size_t) pup_size * (size_t) pup_size;
        const float *d_dphi = d_inv_out + (size_t) (b * n_inv + ch_base) * pad_sq;
        const float *d_chi  = d_inv_out + (size_t) (b * n_inv + ch_base + 1) * pad_sq;
        size_t       src_i  = (size_t) (sy + guard) * (size_t) pad_size + (size_t) (sx + guard);
        size_t       dst_i  = (size_t) b * pup_sq + (size_t) y * (size_t) pup_size + (size_t) x;

        d_frame_pha[dst_i] += d_dphi[src_i];
        float c             = fminf(fmaxf(d_chi[src_i], -3.0f), 3.0f);
        float a             = expf(c);
        d_frame_amp[dst_i]  = a;
        a_sq                = a * a;
    }

    for (int offset = 16; offset > 0; offset /= 2)
    {
        a_sq += __shfl_down_sync(0xffffffff, a_sq, offset);
    }

    __shared__ float s_warp[32];
    int              lane = threadIdx.x + threadIdx.y * blockDim.x;
    if ((lane & 31) == 0)
    {
        s_warp[lane >> 5] = a_sq;
    }
    __syncthreads();

    if (lane < 32)
    {
        int   nwarps = (blockDim.x * blockDim.y + 31) / 32;
        float val    = (lane < nwarps) ? s_warp[lane] : 0.0f;
        for (int offset = 16; offset > 0; offset /= 2)
        {
            val += __shfl_down_sync(0xffffffff, val, offset);
        }
        if (lane == 0)
        {
            atomicAdd(&d_sum_I[b], val);
        }
    }
}

/**
 * rytov_batch_normalize_kernel - Normalize each frame's amplitude by its RMS flux
 */
__global__ static void rytov_batch_normalize_kernel(
    float       *d_frame_amp,
    int          npix,
    int          count,
    const float *d_sum_I)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int b   = (int) blockIdx.y;
    if (idx >= npix || b >= count)
    {
        return;
    }

    float sum_i = d_sum_I[b];
    float norm  = (sum_i > 0.0f) ? rsqrtf(sum_i / (float) npix) : 1.0f;
    d_frame_amp[(size_t) b * (size_t) npix + (size_t) idx] *= norm;
}

void atmturb_cuda_rytov_batch_cleanup(void)
{
    if (s_bat.stream != 0)
    {
        cudaStreamSynchronize(s_bat.stream);
        cudaStreamDestroy(s_bat.stream);
    }
    if (s_bat.plan_1d_many != 0)  { cufftDestroy(s_bat.plan_1d_many); }
    if (s_bat.plan_r2c_many != 0) { cufftDestroy(s_bat.plan_r2c_many); }
    if (s_bat.plan_c2r_many != 0) { cufftDestroy(s_bat.plan_c2r_many); }
    if (s_bat.plan_1d_rem != 0)   { cufftDestroy(s_bat.plan_1d_rem); }
    if (s_bat.plan_r2c_rem != 0)  { cufftDestroy(s_bat.plan_r2c_rem); }
    if (s_bat.plan_c2r_rem != 0)  { cufftDestroy(s_bat.plan_c2r_rem); }

    cudaFree(s_bat.d_laplace_inv);
    cudaFree(s_bat.d_exp_x);
    cudaFree(s_bat.d_exp_y);
    cudaFree(s_bat.d_filt_a_pri);
    cudaFree(s_bat.d_filt_b_pri);
    cudaFree(s_bat.d_filt_a_sec);
    cudaFree(s_bat.d_filt_b_sec);
    cudaFree(s_bat.d_chrom_ramp);
    cudaFree(s_bat.d_super_pha);
    cudaFree(s_bat.d_super_spha);
    cudaFree(s_bat.d_bounds);
    cudaFree(s_bat.d_hat_bounds);
    cudaFree(s_bat.d_spec);
    cudaFree(s_bat.d_acc_batch);
    cudaFree(s_bat.d_inv_out);
    cudaFree(s_bat.d_sum_I);
    for (int s = 0; s < 2; s++)
    {
        cudaFreeHost(s_bat.h_sub[s]);
        cudaFree(s_bat.d_sub[s]);
        cudaFreeHost(s_bat.h_lf[s]);
        cudaFree(s_bat.d_lf[s]);
        cudaFree(s_bat.d_frame_pha[s]);
        cudaFree(s_bat.d_frame_amp[s]);
        cudaFree(s_bat.d_frame_spha[s]);
        cudaFree(s_bat.d_frame_samp[s]);
    }
    memset(&s_bat, 0, sizeof(s_bat));
}

/**
 * rytov_batch_make_plans - Build 1D R2C, 2D R2C, and 2D C2R cuFFT plans for @cnt frames
 */
static int rytov_batch_make_plans(
    int          pad,
    int          n_inv,
    int          cnt,
    cudaStream_t stream,
    cufftHandle *p_1d,
    cufftHandle *p_r2c,
    cufftHandle *p_c2r)
{
    int n1[1]   = {pad};
    int n2[2]   = {pad, pad};
    int n_half  = pad / 2 + 1;
    int n_spec  = pad * n_half;
    int pad_pix = pad * pad;

    if (cufftPlanMany(p_1d, 1, n1, NULL, 1, pad, NULL, 1, n_half, CUFFT_R2C, 2 * cnt) != 0 ||
        cufftPlanMany(p_r2c, 2, n2, NULL, 1, pad_pix, NULL, 1, n_spec, CUFFT_R2C, cnt) != 0 ||
        cufftPlanMany(p_c2r, 2, n2, NULL, 1, n_spec, NULL, 1, pad_pix,
                      CUFFT_C2R, cnt * n_inv) != 0 ||
        cufftSetStream(*p_1d, stream) != 0 ||
        cufftSetStream(*p_r2c, stream) != 0 ||
        cufftSetStream(*p_c2r, stream) != 0)
    {
        return -1;
    }
    return 0;
}

/**
 * rytov_batch_ensure_dynamic - Ensure H2D sublayer and lowfreq slot buffers meet capacity
 */
static int rytov_batch_ensure_dynamic(
    size_t sub_bytes,
    size_t lf_bytes)
{
    if (sub_bytes > s_bat.sub_cap)
    {
        for (int s = 0; s < 2; s++)
        {
            cudaFreeHost(s_bat.h_sub[s]); s_bat.h_sub[s] = NULL;
            cudaFree(s_bat.d_sub[s]);     s_bat.d_sub[s] = NULL;
            if (cudaMallocHost((void **) &s_bat.h_sub[s], sub_bytes) != cudaSuccess ||
                cudaMalloc((void **) &s_bat.d_sub[s], sub_bytes) != cudaSuccess)
            {
                return -1;
            }
        }
        s_bat.sub_cap = sub_bytes;
    }
    if (lf_bytes > s_bat.lf_cap)
    {
        for (int s = 0; s < 2; s++)
        {
            cudaFreeHost(s_bat.h_lf[s]); s_bat.h_lf[s] = NULL;
            cudaFree(s_bat.d_lf[s]);     s_bat.d_lf[s] = NULL;
            if (cudaMallocHost((void **) &s_bat.h_lf[s], lf_bytes) != cudaSuccess ||
                cudaMalloc((void **) &s_bat.d_lf[s], lf_bytes) != cudaSuccess)
            {
                return -1;
            }
        }
        s_bat.lf_cap = lf_bytes;
    }
    return 0;
}

/**
 * rytov_batch_ensure_state - Allocate or reuse persistent batched Rytov GPU workspace
 */
static int rytov_batch_ensure_state(
    const atmturb_cuda_rytov_params_t *p,
    int                                batch_b,
    int                                batch_rem)
{
    size_t sub_byt = sizeof(atmturb_cuda_rytov_sublayer_t) *
                     (size_t) p->nsuper * (size_t) batch_b * (size_t) p->max_sublayers;
    size_t lf_byt  = (p->lowfreq && p->max_lf_modes > 0)
                         ? (sizeof(atmturb_cuda_lowfreq_mode_t) *
                            (size_t) p->nsuper * (size_t) batch_b * (size_t) p->max_lf_modes)
                         : 0;
    size_t fbytes  = sizeof(float) * (size_t) (p->pup_size * p->pup_size) * (size_t) batch_b;

    if (atmturb_cuda_dma_ring_prepare(fbytes) != 0)
    {
        return -1;
    }

    if (s_bat.stream != 0 && s_bat.pad == p->pad_size && s_bat.pup == p->pup_size &&
        s_bat.has_sec == p->has_sec && s_bat.nsuper >= p->nsuper &&
        s_bat.batch_b == batch_b && s_bat.batch_rem == batch_rem)
    {
        return rytov_batch_ensure_dynamic(sub_byt, lf_byt);
    }

    atmturb_cuda_rytov_batch_cleanup();

    long   pad    = p->pad_size;
    long   n_half = pad / 2 + 1;
    size_t n_spec = (size_t) pad * (size_t) n_half;
    size_t pad_sq = (size_t) pad * (size_t) pad;
    size_t B      = (size_t) batch_b;
    size_t ns     = (size_t) p->nsuper;
    int    n_inv  = p->has_sec ? 4 : 2;

    if (cudaStreamCreateWithFlags(&s_bat.stream, cudaStreamNonBlocking) != cudaSuccess ||
        rytov_batch_make_plans((int) pad, n_inv, batch_b, s_bat.stream,
                               &s_bat.plan_1d_many, &s_bat.plan_r2c_many,
                               &s_bat.plan_c2r_many) != 0 ||
        (batch_rem > 0 &&
         rytov_batch_make_plans((int) pad, n_inv, batch_rem, s_bat.stream,
                                &s_bat.plan_1d_rem, &s_bat.plan_r2c_rem,
                                &s_bat.plan_c2r_rem) != 0))
    {
        return -1;
    }

    if (cudaMalloc((void **) &s_bat.d_laplace_inv, sizeof(float) * n_spec) != cudaSuccess ||
        cudaMalloc((void **) &s_bat.d_exp_x, sizeof(cufftComplex) * (size_t) n_half) != 0 ||
        cudaMalloc((void **) &s_bat.d_exp_y, sizeof(cufftComplex) * (size_t) pad) != 0 ||
        cudaMalloc((void **) &s_bat.d_filt_a_pri, sizeof(float) * ns * n_spec) != 0 ||
        cudaMalloc((void **) &s_bat.d_filt_b_pri, sizeof(float) * ns * n_spec) != 0 ||
        cudaMalloc((void **) &s_bat.d_filt_a_sec, sizeof(float) * ns * n_spec) != 0 ||
        cudaMalloc((void **) &s_bat.d_filt_b_sec, sizeof(float) * ns * n_spec) != 0 ||
        cudaMalloc((void **) &s_bat.d_chrom_ramp, sizeof(cufftComplex) * ns * n_spec) != 0 ||
        cudaMalloc((void **) &s_bat.d_super_pha, sizeof(float) * B * pad_sq) != 0 ||
        cudaMalloc((void **) &s_bat.d_super_spha, sizeof(float) * B * pad_sq) != 0 ||
        cudaMalloc((void **) &s_bat.d_bounds, sizeof(float) * B * (size_t) (2 * pad)) != 0 ||
        cudaMalloc((void **) &s_bat.d_hat_bounds,
                   sizeof(cufftComplex) * B * (size_t) (2 * n_half)) != 0 ||
        cudaMalloc((void **) &s_bat.d_spec, sizeof(cufftComplex) * B * n_spec) != 0 ||
        cudaMalloc((void **) &s_bat.d_acc_batch,
                   sizeof(cufftComplex) * B * (size_t) n_inv * n_spec) != 0 ||
        cudaMalloc((void **) &s_bat.d_inv_out, sizeof(float) * B * (size_t) n_inv * pad_sq) != 0 ||
        cudaMalloc((void **) &s_bat.d_sum_I, sizeof(float) * B) != 0)
    {
        return -1;
    }

    for (int s = 0; s < 2; s++)
    {
        if (cudaMalloc((void **) &s_bat.d_frame_pha[s], fbytes) != cudaSuccess ||
            cudaMalloc((void **) &s_bat.d_frame_amp[s], fbytes) != cudaSuccess ||
            (p->has_sec &&
             (cudaMalloc((void **) &s_bat.d_frame_spha[s], fbytes) != cudaSuccess ||
              cudaMalloc((void **) &s_bat.d_frame_samp[s], fbytes) != cudaSuccess)))
        {
            return -1;
        }
    }

    s_bat.pad       = p->pad_size;
    s_bat.pup       = p->pup_size;
    s_bat.has_sec   = p->has_sec;
    s_bat.nsuper    = p->nsuper;
    s_bat.batch_b   = batch_b;
    s_bat.batch_rem = batch_rem;
    return rytov_batch_ensure_dynamic(sub_byt, lf_byt);
}

/**
 * rytov_batch_upload_filters - Asynchronously upload static Rytov filters to device
 */
static void rytov_batch_upload_filters(
    const atmturb_cuda_rytov_params_t *p)
{
    long   pad    = p->pad_size;
    long   n_half = pad / 2 + 1;
    size_t n_spec = (size_t) pad * (size_t) n_half;

    if (p->use_moisan)
    {
        cudaMemcpyAsync(s_bat.d_laplace_inv, p->h_laplace_inv, sizeof(float) * n_spec,
                        cudaMemcpyHostToDevice, s_bat.stream);
        cudaMemcpyAsync(s_bat.d_exp_x, p->h_exp_x, sizeof(cufftComplex) * (size_t) n_half,
                        cudaMemcpyHostToDevice, s_bat.stream);
        cudaMemcpyAsync(s_bat.d_exp_y, p->h_exp_y, sizeof(cufftComplex) * (size_t) pad,
                        cudaMemcpyHostToDevice, s_bat.stream);
    }

    for (int m = 0; m < p->nsuper; m++)
    {
        if (p->super_dist_m[m] <= 0.0)
        {
            continue;
        }
        size_t off = (size_t) m * n_spec;
        cudaMemcpyAsync(s_bat.d_filt_a_pri + off, p->h_filt_a_pri[m], sizeof(float) * n_spec,
                        cudaMemcpyHostToDevice, s_bat.stream);
        cudaMemcpyAsync(s_bat.d_filt_b_pri + off, p->h_filt_b_pri[m], sizeof(float) * n_spec,
                        cudaMemcpyHostToDevice, s_bat.stream);
        if (p->has_sec)
        {
            cudaMemcpyAsync(s_bat.d_filt_a_sec + off, p->h_filt_a_sec[m], sizeof(float) * n_spec,
                            cudaMemcpyHostToDevice, s_bat.stream);
            cudaMemcpyAsync(s_bat.d_filt_b_sec + off, p->h_filt_b_sec[m], sizeof(float) * n_spec,
                            cudaMemcpyHostToDevice, s_bat.stream);
            if (p->sec_shared && p->h_chrom_ramp != NULL && p->h_chrom_ramp[m] != NULL)
            {
                cudaMemcpyAsync(s_bat.d_chrom_ramp + off, p->h_chrom_ramp[m],
                                sizeof(cufftComplex) * n_spec, cudaMemcpyHostToDevice,
                                s_bat.stream);
            }
        }
    }
}

/**
 * rytov_batch_accum_super - Transform and accumulate spectrum for superlayer @m across @count
 */
static void rytov_batch_accum_super(
    const atmturb_cuda_rytov_params_t *p,
    int                                m,
    int                                count,
    cufftHandle                        p_1d,
    cufftHandle                        p_r2c)
{
    int    pad     = (int) p->pad_size;
    int    n_half  = pad / 2 + 1;
    size_t off     = (size_t) m * (size_t) pad * (size_t) n_half;
    int    is_dual = (p->has_sec && p->sec_shared && p->h_chrom_ramp != NULL &&
                      p->h_chrom_ramp[m] != NULL);
    dim3   b_1d(256, 1, 1);
    dim3   g_1d((unsigned int) ((pad + 255) / 256), (unsigned int) count, 1);
    dim3   b_sp(16, 16, 1);
    dim3   g_sp((unsigned int) ((n_half + 15) / 16), (unsigned int) ((pad + 15) / 16),
                (unsigned int) count);

    if (p->use_moisan)
    {
        rytov_batch_bounds_kernel<<<g_1d, b_1d, 0, s_bat.stream>>>(
            s_bat.d_super_pha, s_bat.d_bounds, pad, count);
        cufftExecR2C(p_1d, (cufftReal *) s_bat.d_bounds, (cufftComplex *) s_bat.d_hat_bounds);
    }
    cufftExecR2C(p_r2c, (cufftReal *) s_bat.d_super_pha, (cufftComplex *) s_bat.d_spec);

    rytov_batch_mac_args_t ma = {
        .d_acc_batch   = s_bat.d_acc_batch,
        .d_spec        = s_bat.d_spec,
        .d_hat_bounds  = s_bat.d_hat_bounds,
        .d_ramp        = s_bat.d_chrom_ramp + off,
        .d_exp_x       = s_bat.d_exp_x,
        .d_exp_y       = s_bat.d_exp_y,
        .d_laplace_inv = s_bat.d_laplace_inv,
        .d_filt_a_pri  = s_bat.d_filt_a_pri + off,
        .d_filt_b_pri  = s_bat.d_filt_b_pri + off,
        .d_filt_a_sec  = s_bat.d_filt_a_sec + off,
        .d_filt_b_sec  = s_bat.d_filt_b_sec + off,
        .pad           = pad,
        .n_half        = n_half,
        .n_inv         = p->has_sec ? 4 : 2,
        .ch_base       = 0,
        .count         = count
    };

    if (p->use_moisan)
    {
        if (is_dual)
        {
            rytov_batch_spec_mac_kernel<1, 1><<<g_sp, b_sp, 0, s_bat.stream>>>(ma);
            return;
        }
        rytov_batch_spec_mac_kernel<1, 0><<<g_sp, b_sp, 0, s_bat.stream>>>(ma);
        if (p->has_sec)
        {
            rytov_batch_bounds_kernel<<<g_1d, b_1d, 0, s_bat.stream>>>(
                s_bat.d_super_spha, s_bat.d_bounds, pad, count);
            cufftExecR2C(p_1d, (cufftReal *) s_bat.d_bounds, (cufftComplex *) s_bat.d_hat_bounds);
            cufftExecR2C(p_r2c, (cufftReal *) s_bat.d_super_spha, (cufftComplex *) s_bat.d_spec);
            ma.ch_base      = 2;
            ma.d_filt_a_pri = s_bat.d_filt_a_sec + off;
            ma.d_filt_b_pri = s_bat.d_filt_b_sec + off;
            rytov_batch_spec_mac_kernel<1, 0><<<g_sp, b_sp, 0, s_bat.stream>>>(ma);
        }
    }
    else
    {
        if (is_dual)
        {
            rytov_batch_spec_mac_kernel<0, 1><<<g_sp, b_sp, 0, s_bat.stream>>>(ma);
            return;
        }
        rytov_batch_spec_mac_kernel<0, 0><<<g_sp, b_sp, 0, s_bat.stream>>>(ma);
        if (p->has_sec)
        {
            cufftExecR2C(p_r2c, (cufftReal *) s_bat.d_super_spha, (cufftComplex *) s_bat.d_spec);
            ma.ch_base      = 2;
            ma.d_filt_a_pri = s_bat.d_filt_a_sec + off;
            ma.d_filt_b_pri = s_bat.d_filt_b_sec + off;
            rytov_batch_spec_mac_kernel<0, 0><<<g_sp, b_sp, 0, s_bat.stream>>>(ma);
        }
    }
}

/**
 * rytov_batch_finalize_chunk - Inverse C2R, pupil assembly, flux norm, and optional piston sub
 */
static void rytov_batch_finalize_chunk(
    const atmturb_cuda_rytov_params_t *p,
    int                                slot,
    int                                count,
    int                                any_diff,
    cufftHandle                        p_c2r)
{
    int pup  = (int) p->pup_size;
    int npix = pup * pup;
    int n_inv = p->has_sec ? 4 : 2;

    if (any_diff)
    {
        cufftExecC2R(p_c2r, (cufftComplex *) s_bat.d_acc_batch, (cufftReal *) s_bat.d_inv_out);

        dim3 b_ass(16, 16, 1);
        dim3 g_ass((unsigned int) ((pup + 15) / 16), (unsigned int) ((pup + 15) / 16),
                   (unsigned int) count);
        dim3 b_nrm(256, 1, 1);
        dim3 g_nrm((unsigned int) ((npix + 255) / 256), (unsigned int) count, 1);

        cudaMemsetAsync(s_bat.d_sum_I, 0, sizeof(float) * (size_t) count, s_bat.stream);
        rytov_batch_assemble_kernel<<<g_ass, b_ass, 0, s_bat.stream>>>(
            pup, (int) p->guard_pix, (int) p->pad_size, n_inv, 0,
            s_bat.d_frame_pha[slot], s_bat.d_frame_amp[slot], s_bat.d_inv_out, s_bat.d_sum_I);
        rytov_batch_normalize_kernel<<<g_nrm, b_nrm, 0, s_bat.stream>>>(
            s_bat.d_frame_amp[slot], npix, count, s_bat.d_sum_I);

        if (p->has_sec)
        {
            cudaMemsetAsync(s_bat.d_sum_I, 0, sizeof(float) * (size_t) count, s_bat.stream);
            rytov_batch_assemble_kernel<<<g_ass, b_ass, 0, s_bat.stream>>>(
                pup, (int) p->guard_pix, (int) p->pad_size, n_inv, 2,
                s_bat.d_frame_spha[slot], s_bat.d_frame_samp[slot], s_bat.d_inv_out,
                s_bat.d_sum_I);
            rytov_batch_normalize_kernel<<<g_nrm, b_nrm, 0, s_bat.stream>>>(
                s_bat.d_frame_samp[slot], npix, count, s_bat.d_sum_I);
        }
    }
    else
    {
        atmturb_cuda_fill_ones_launch(s_bat.d_frame_amp[slot], count * npix, s_bat.stream);
        if (p->has_sec)
        {
            atmturb_cuda_fill_ones_launch(s_bat.d_frame_samp[slot], count * npix, s_bat.stream);
        }
    }

    if (p->remove_piston)
    {
        for (int b = 0; b < count; b++)
        {
            atmturb_cuda_remove_piston_launch(
                s_bat.d_frame_pha[slot] + (size_t) b * (size_t) npix,
                &s_bat.d_sum_I[b], npix, s_bat.stream);
            if (p->has_sec)
            {
                atmturb_cuda_remove_piston_launch(
                    s_bat.d_frame_spha[slot] + (size_t) b * (size_t) npix,
                    &s_bat.d_sum_I[b], npix, s_bat.stream);
            }
        }
    }
}

/**
 * rytov_batch_run_chunk - Enqueue one chunk of @count frames and hand off to async D2H ring
 */
static int rytov_batch_run_chunk(
    const atmturb_cuda_rytov_params_t *p,
    const float                       *d_masters,
    int                                slot,
    long                               t0,
    int                                count,
    int                                any_diff)
{
    if (atmturb_cuda_dma_ring_drain(slot) != 0)
    {
        return -1;
    }

    size_t sub_frame = (size_t) (p->nsuper * p->max_sublayers);
    size_t sub_bytes = sizeof(atmturb_cuda_rytov_sublayer_t) * (size_t) count * sub_frame;
    memcpy(s_bat.h_sub[slot], p->frame_sublayers + (size_t) t0 * sub_frame, sub_bytes);
    cudaMemcpyAsync(s_bat.d_sub[slot], s_bat.h_sub[slot], sub_bytes,
                    cudaMemcpyHostToDevice, s_bat.stream);

    size_t lf_frame = (size_t) (p->nsuper * p->max_lf_modes);
    if (p->lowfreq && p->max_lf_modes > 0 && p->frame_lf_modes != NULL)
    {
        size_t lf_bytes = sizeof(atmturb_cuda_lowfreq_mode_t) * (size_t) count * lf_frame;
        memcpy(s_bat.h_lf[slot], p->frame_lf_modes + (size_t) t0 * lf_frame, lf_bytes);
        cudaMemcpyAsync(s_bat.d_lf[slot], s_bat.h_lf[slot], lf_bytes,
                        cudaMemcpyHostToDevice, s_bat.stream);
    }

    size_t pup_sq     = (size_t) p->pup_size * (size_t) p->pup_size;
    size_t chunk_fbyt = sizeof(float) * (size_t) count * pup_sq;
    int    n_inv      = p->has_sec ? 4 : 2;
    size_t n_spec     = (size_t) p->pad_size * (size_t) (p->pad_size / 2 + 1);

    cudaMemsetAsync(s_bat.d_frame_pha[slot], 0, chunk_fbyt, s_bat.stream);
    if (p->has_sec)
    {
        cudaMemsetAsync(s_bat.d_frame_spha[slot], 0, chunk_fbyt, s_bat.stream);
    }
    if (any_diff)
    {
        cudaMemsetAsync(s_bat.d_acc_batch, 0,
                        sizeof(cufftComplex) * (size_t) (count * n_inv) * n_spec, s_bat.stream);
    }

    int         is_full = (count == s_bat.batch_b);
    cufftHandle p_1d    = is_full ? s_bat.plan_1d_many : s_bat.plan_1d_rem;
    cufftHandle p_r2c   = is_full ? s_bat.plan_r2c_many : s_bat.plan_r2c_rem;
    cufftHandle p_c2r   = is_full ? s_bat.plan_c2r_many : s_bat.plan_c2r_rem;

    for (int m = 0; m < p->nsuper; m++)
    {
        int is_diff = (p->super_dist_m[m] > 0.0) ? 1 : 0;
        int is_dual = (p->has_sec && p->sec_shared && p->h_chrom_ramp != NULL &&
                       p->h_chrom_ramp[m] != NULL);
        int n_lf    = (p->lowfreq && p->super_nlf_modes != NULL && p->frame_lf_modes != NULL)
                          ? p->super_nlf_modes[m] : 0;

        atmturb_cuda_extrude_batch_args_t ext = {
            .interp_mode  = p->interp,
            .d_masters    = d_masters,
            .msize        = (int) p->msize,
            .pad_size     = (int) p->pad_size,
            .pup_size     = (int) p->pup_size,
            .guard        = (int) p->guard_pix,
            .n_sub        = p->super_nlayers[m],
            .sub_stride   = (int) sub_frame,
            .has_sec      = p->has_sec,
            .os           = (p->os > 1) ? p->os : 1,
            .batch_count  = count,
            .d_sublayers  = s_bat.d_sub[slot] + (size_t) (m * p->max_sublayers),
            .n_lf_modes   = n_lf,
            .lf_stride    = (int) lf_frame,
            .d_lf_modes   = (n_lf > 0) ? (s_bat.d_lf[slot] + (size_t) (m * p->max_lf_modes)) : NULL,
            .d_super_pha  = is_diff ? s_bat.d_super_pha : NULL,
            .d_super_spha = (is_diff && p->has_sec && !is_dual) ? s_bat.d_super_spha : NULL,
            .d_frame_pha  = s_bat.d_frame_pha[slot],
            .d_frame_spha = p->has_sec ? s_bat.d_frame_spha[slot] : NULL
        };
        atmturb_cuda_extrude_sl_batch_launch(&ext, s_bat.stream);

        if (is_diff)
        {
            rytov_batch_accum_super(p, m, count, p_1d, p_r2c);
        }
    }

    rytov_batch_finalize_chunk(p, slot, count, any_diff, p_c2r);

    size_t             dst_off  = (size_t) t0 * pup_sq;
    const float *const d_src[4] = {
        s_bat.d_frame_pha[slot], s_bat.d_frame_amp[slot],
        s_bat.d_frame_spha[slot], s_bat.d_frame_samp[slot]
    };
    float *const h_dst[4] = {
        (p->pha != NULL) ? (p->pha + dst_off) : NULL,
        (p->amp != NULL) ? (p->amp + dst_off) : NULL,
        (p->has_sec && p->spha != NULL) ? (p->spha + dst_off) : NULL,
        (p->has_sec && p->samp != NULL) ? (p->samp + dst_off) : NULL
    };
    return atmturb_cuda_dma_ring_enqueue(slot, s_bat.stream, d_src, h_dst, chunk_fbyt);
}

int atmturb_cuda_rytov_render_batch(
    const atmturb_cuda_rytov_params_t *params)
{
    if (params == NULL || params->nblayers <= 0 || params->pup_size <= 0 || params->nbframes < 2)
    {
        return -1;
    }

    const float *d_masters = atmturb_cuda_sync_device_masters(
        params->msize, params->nblayers, params->h_masters);
    if (d_masters == NULL)
    {
        return -1;
    }

    int max_b     = (params->pad_size <= 384) ? 16 : ((params->pad_size <= 640) ? 8 : 4);
    int batch_b   = (params->nbframes < (long) max_b) ? (int) params->nbframes : max_b;
    int batch_rem = (int) (params->nbframes % (long) batch_b);

    if (rytov_batch_ensure_state(params, batch_b, batch_rem) != 0)
    {
        atmturb_cuda_rytov_batch_cleanup();
        return -1;
    }
    rytov_batch_upload_filters(params);

    int any_diff = 0;
    for (int m = 0; m < params->nsuper; m++)
    {
        if (params->super_dist_m[m] > 0.0)
        {
            any_diff = 1;
            break;
        }
    }

    int ret   = 0;
    int chunk = 0;
    for (long t0 = 0; t0 < params->nbframes; t0 += (long) batch_b, chunk++)
    {
        long rem   = params->nbframes - t0;
        int  count = (rem < (long) batch_b) ? (int) rem : batch_b;
        if (rytov_batch_run_chunk(params, d_masters, chunk & 1, t0, count, any_diff) != 0)
        {
            ret = -1;
            break;
        }
    }

    if (atmturb_cuda_dma_ring_flush() != 0)
    {
        ret = -1;
    }
    return ret;
}
