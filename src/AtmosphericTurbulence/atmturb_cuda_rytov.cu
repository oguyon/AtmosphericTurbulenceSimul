// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_rytov.cu
 * @brief   CUDA and cuFFT accelerated multi-layer Rytov propagation
 */

#include <cuda_runtime.h>
#include <cufft.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "atmturb_cuda.h"
#include "atmturb_cuda_rytov.h"

#define ATMTURB_CUDA_MAX_SUB_LAYERS 64

__constant__ atmturb_cuda_rytov_sublayer_t c_sublayers[ATMTURB_CUDA_MAX_SUB_LAYERS];

// Persistent cuFFT plans and configuration cache
static cufftHandle s_plan_1d     = 0;
static cufftHandle s_plan_2d_r2c = 0;
static cufftHandle s_plan_2d_c2r = 0;
static long        s_cached_pad  = 0;

// Persistent device precomputed filter buffers
static float        *s_d_laplace_inv  = NULL;
static cufftComplex *s_d_exp_x        = NULL;
static cufftComplex *s_d_exp_y        = NULL;
static float        *s_d_filt_a_pri   = NULL;
static float        *s_d_filt_b_pri   = NULL;
static float        *s_d_filt_a_sec   = NULL;
static float        *s_d_filt_b_sec   = NULL;
static cufftComplex *s_d_chrom_ramp   = NULL;
static size_t        s_cached_filters = 0;

// Persistent device frame working buffers
static float        *s_d_super_pha    = NULL;
static float        *s_d_super_spha   = NULL;
static float        *s_d_bounds       = NULL;
static cufftComplex *s_d_hat_bounds   = NULL;
static cufftComplex *s_d_spec         = NULL;
static cufftComplex *s_d_acc_dphi_pri = NULL;
static cufftComplex *s_d_acc_chi_pri  = NULL;
static cufftComplex *s_d_acc_dphi_sec = NULL;
static cufftComplex *s_d_acc_chi_sec  = NULL;
static float        *s_d_dphi_out     = NULL;
static float        *s_d_chi_out      = NULL;
static float        *s_d_frame_pha    = NULL;
static float        *s_d_frame_amp    = NULL;
static float        *s_d_frame_spha   = NULL;
static float        *s_d_frame_samp   = NULL;
static float        *s_d_sum_I        = NULL;
static size_t        s_cached_work    = 0;

__device__ static inline float atmturb_cuda_sample_screen(
    const float *d_masters,
    int          msize,
    int          is_pow2,
    int          mask,
    int          k,
    float        cur_x,
    float        cur_y)
{
    int ix = __float2int_rd(cur_x);
    int iy = __float2int_rd(cur_y);
    float fx = cur_x - (float) ix;
    float fy = cur_y - (float) iy;

    int ix0, iy0, ix1, iy1;
    if (is_pow2)
    {
        ix0 = ix & mask;
        iy0 = iy & mask;
        ix1 = (ix0 + 1) & mask;
        iy1 = (iy0 + 1) & mask;
    }
    else
    {
        ix0 = ix % msize;
        if (ix0 < 0) ix0 += msize;
        iy0 = iy % msize;
        if (iy0 < 0) iy0 += msize;
        ix1 = (ix0 + 1 == msize) ? 0 : (ix0 + 1);
        iy1 = (iy0 + 1 == msize) ? 0 : (iy0 + 1);
    }

    size_t msize_sq = (size_t) msize * (size_t) msize;
    const float *screen = d_masters + (size_t) k * msize_sq;
    float v00 = __ldg(&screen[iy0 * msize + ix0]);
    float v10 = __ldg(&screen[iy0 * msize + ix1]);
    float v01 = __ldg(&screen[iy1 * msize + ix0]);
    float v11 = __ldg(&screen[iy1 * msize + ix1]);

    float top = v00 + fx * (v10 - v00);
    float bot = v01 + fx * (v11 - v01);
    return top + fy * (bot - top);
}

__device__ static inline float atmturb_cuda_sample_screen_bicubic(
    const float *d_masters,
    int          msize,
    int          is_pow2,
    int          mask,
    int          k,
    float        cur_x,
    float        cur_y)
{
    int ix = __float2int_rd(cur_x);
    int iy = __float2int_rd(cur_y);
    float fx = cur_x - (float) ix;
    float fy = cur_y - (float) iy;

    float fx2 = fx * fx;
    float fx3 = fx2 * fx;
    float wx[4];
    wx[0] = -0.5f * fx + fx2 - 0.5f * fx3;
    wx[1] = 1.0f - 2.5f * fx2 + 1.5f * fx3;
    wx[2] = 0.5f * fx + 2.0f * fx2 - 1.5f * fx3;
    wx[3] = -0.5f * fx2 + 0.5f * fx3;

    float fy2 = fy * fy;
    float fy3 = fy2 * fy;
    float wy[4];
    wy[0] = -0.5f * fy + fy2 - 0.5f * fy3;
    wy[1] = 1.0f - 2.5f * fy2 + 1.5f * fy3;
    wy[2] = 0.5f * fy + 2.0f * fy2 - 1.5f * fy3;
    wy[3] = -0.5f * fy2 + 0.5f * fy3;

    int ix_arr[4], iy_arr[4];
    if (is_pow2)
    {
        for (int m = 0; m < 4; m++)
        {
            ix_arr[m] = (ix + m - 1) & mask;
            iy_arr[m] = (iy + m - 1) & mask;
        }
    }
    else
    {
        for (int m = 0; m < 4; m++)
        {
            int xm = (ix + m - 1) % msize;
            if (xm < 0)
            {
                xm += msize;
            }
            ix_arr[m] = xm;

            int ym = (iy + m - 1) % msize;
            if (ym < 0)
            {
                ym += msize;
            }
            iy_arr[m] = ym;
        }
    }

    size_t msize_sq = (size_t) msize * (size_t) msize;
    const float *screen = d_masters + (size_t) k * msize_sq;

    float val = 0.0f;
    for (int n = 0; n < 4; n++)
    {
        const float *row = screen + (size_t) iy_arr[n] * msize;
        float h = wx[0] * __ldg(&row[ix_arr[0]]) +
                  wx[1] * __ldg(&row[ix_arr[1]]) +
                  wx[2] * __ldg(&row[ix_arr[2]]) +
                  wx[3] * __ldg(&row[ix_arr[3]]);
        val += wy[n] * h;
    }
    return val;
}

template <int interp_mode>
__global__ static void atmturb_cuda_extrude_sl_kernel(
    const float *d_masters,
    int          msize,
    int          pad_size,
    int          pup_size,
    int          guard,
    int          n_sub,
    int          has_sec,
    int          os,
    float       *d_super_pha,
    float       *d_super_spha,
    float       *d_frame_pha,
    float       *d_frame_spha)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= pad_size || j >= pad_size)
    {
        return;
    }

    int is_pow2 = ((msize & (msize - 1)) == 0);
    int mask    = msize - 1;

    float p_val  = 0.0f;
    float sp_val = 0.0f;

    for (int l = 0; l < n_sub; l++)
    {
        atmturb_cuda_rytov_sublayer_t info = c_sublayers[l];
        float cur_x  = info.x + (float) (i * os);
        float cur_y  = info.y + (float) (j * os);
        if (interp_mode == 1)
        {
            p_val += info.w_pri * atmturb_cuda_sample_screen_bicubic(d_masters, msize, is_pow2,
                                                                     mask, info.k, cur_x, cur_y);
        }
        else
        {
            p_val += info.w_pri * atmturb_cuda_sample_screen(d_masters, msize, is_pow2,
                                                             mask, info.k, cur_x, cur_y);
        }
        if (has_sec)
        {
            float cur_xs = info.xs + (float) (i * os);
            float cur_ys = info.ys + (float) (j * os);
            if (interp_mode == 1)
            {
                sp_val += info.w_sec * atmturb_cuda_sample_screen_bicubic(
                    d_masters, msize, is_pow2, mask, info.k, cur_xs, cur_ys);
            }
            else
            {
                sp_val += info.w_sec * atmturb_cuda_sample_screen(
                    d_masters, msize, is_pow2, mask, info.k, cur_xs, cur_ys);
            }
        }
    }

    int pad_idx = j * pad_size + i;
    d_super_pha[pad_idx] = p_val;
    if (has_sec)
    {
        d_super_spha[pad_idx] = sp_val;
    }

    if (i >= guard && i < guard + pup_size && j >= guard && j < guard + pup_size)
    {
        int pup_idx = (j - guard) * pup_size + (i - guard);
        d_frame_pha[pup_idx] += p_val;
        if (has_sec)
        {
            d_frame_spha[pup_idx] += sp_val;
        }
    }
}

__global__ static void atmturb_cuda_extract_bounds_kernel(
    const float *d_super_pha,
    float       *d_bounds,
    int          pad_size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= pad_size)
    {
        return;
    }

    d_bounds[0 * pad_size + idx] = d_super_pha[0 * pad_size + idx] -
                                   d_super_pha[(pad_size - 1) * pad_size + idx];
    d_bounds[1 * pad_size + idx] = d_super_pha[idx * pad_size + 0] -
                                   d_super_pha[idx * pad_size + (pad_size - 1)];
}

__global__ static void atmturb_cuda_moisan_adjust_kernel(
    cufftComplex       *d_spec,
    const cufftComplex *d_hat_bounds,
    const cufftComplex *d_exp_x,
    const cufftComplex *d_exp_y,
    const float        *d_laplace_inv,
    int                 pad_size,
    int                 n_half)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= n_half || j >= pad_size)
    {
        return;
    }

    const cufftComplex *hat_a = d_hat_bounds;
    const cufftComplex *hat_b = d_hat_bounds + n_half;

    float b_re, b_im;
    if (j <= pad_size / 2)
    {
        b_re = hat_b[j].x;
        b_im = hat_b[j].y;
    }
    else
    {
        int j_sym = pad_size - j;
        b_re =  hat_b[j_sym].x;
        b_im = -hat_b[j_sym].y;
    }

    cufftComplex ey = d_exp_y[j];
    cufftComplex ex = d_exp_x[i];
    cufftComplex a  = hat_a[i];

    float v_re = (a.x * ey.x - a.y * ey.y) + (b_re * ex.x - b_im * ex.y);
    float v_im = (a.x * ey.y + a.y * ey.x) + (b_re * ex.y + b_im * ex.x);

    int spec_idx = j * n_half + i;
    float linv   = d_laplace_inv[spec_idx];
    d_spec[spec_idx].x += v_re * linv;
    d_spec[spec_idx].y += v_im * linv;
}

__global__ static void atmturb_cuda_filter_mac_kernel(
    cufftComplex       *d_acc_dphi,
    cufftComplex       *d_acc_chi,
    const cufftComplex *d_spec,
    const float        *d_filt_a,
    const float        *d_filt_b,
    int                 ntot)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= ntot)
    {
        return;
    }

    cufftComplex s = d_spec[k];
    float fa = d_filt_a[k];
    float fb = d_filt_b[k];

    d_acc_dphi[k].x += s.x * fa;
    d_acc_dphi[k].y += s.y * fa;
    d_acc_chi[k].x  += s.x * fb;
    d_acc_chi[k].y  += s.y * fb;
}

__global__ static void atmturb_cuda_filter_mac_rotated_kernel(
    cufftComplex       *d_acc_dphi_sec,
    cufftComplex       *d_acc_chi_sec,
    const cufftComplex *d_spec,
    const cufftComplex *d_ramp,
    const float        *d_filt_a_sec,
    const float        *d_filt_b_sec,
    int                 ntot)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= ntot)
    {
        return;
    }

    cufftComplex s = d_spec[k];
    cufftComplex r = d_ramp[k];

    float ps_x = s.x * r.x - s.y * r.y;
    float ps_y = s.x * r.y + s.y * r.x;

    float fa = d_filt_a_sec[k];
    float fb = d_filt_b_sec[k];

    d_acc_dphi_sec[k].x += ps_x * fa;
    d_acc_dphi_sec[k].y += ps_y * fa;
    d_acc_chi_sec[k].x  += ps_x * fb;
    d_acc_chi_sec[k].y  += ps_y * fb;
}

__global__ static void atmturb_cuda_assemble_kernel(
    int          pup_size,
    int          guard,
    int          pad_size,
    float       *d_frame_pha,
    float       *d_frame_amp,
    const float *d_dphi_out,
    const float *d_chi_out,
    float       *d_sum_I)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    float a_sq = 0.0f;
    if (x < pup_size && y < pup_size)
    {
        int src_idx = (y + guard) * pad_size + (x + guard);
        int dst_idx = y * pup_size + x;

        d_frame_pha[dst_idx] += d_dphi_out[src_idx];
        float a = expf(d_chi_out[src_idx]);
        d_frame_amp[dst_idx] = a;
        a_sq = a * a;
    }

    for (int offset = 16; offset > 0; offset /= 2)
    {
        a_sq += __shfl_down_sync(0xffffffff, a_sq, offset);
    }

    __shared__ float s_warp_sums[32];
    int lane = threadIdx.x + threadIdx.y * blockDim.x;
    int warp_id = lane / 32;
    int lane_id = lane % 32;

    if (lane_id == 0)
    {
        s_warp_sums[warp_id] = a_sq;
    }
    __syncthreads();

    if (lane < 32)
    {
        int nwarps = (blockDim.x * blockDim.y + 31) / 32;
        float val = (lane < nwarps) ? s_warp_sums[lane] : 0.0f;
        for (int offset = 16; offset > 0; offset /= 2)
        {
            val += __shfl_down_sync(0xffffffff, val, offset);
        }
        if (lane == 0)
        {
            atomicAdd(d_sum_I, val);
        }
    }
}

__global__ static void atmturb_cuda_normalize_kernel(
    float       *d_frame_amp,
    int          npix,
    const float *d_sum_I)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= npix)
    {
        return;
    }

    float sum_i = *d_sum_I;
    float norm = (sum_i > 0.0f) ? (1.0f / sqrtf(sum_i / (float) npix)) : 1.0f;
    d_frame_amp[idx] *= norm;
}

void atmturb_cuda_rytov_cleanup(void)
{
    if (s_plan_1d != 0)     { cufftDestroy(s_plan_1d);     s_plan_1d = 0; }
    if (s_plan_2d_r2c != 0) { cufftDestroy(s_plan_2d_r2c); s_plan_2d_r2c = 0; }
    if (s_plan_2d_c2r != 0) { cufftDestroy(s_plan_2d_c2r); s_plan_2d_c2r = 0; }
    s_cached_pad = 0;

    if (s_d_laplace_inv != NULL) { cudaFree(s_d_laplace_inv); s_d_laplace_inv = NULL; }
    if (s_d_exp_x != NULL)       { cudaFree(s_d_exp_x);       s_d_exp_x = NULL; }
    if (s_d_exp_y != NULL)       { cudaFree(s_d_exp_y);       s_d_exp_y = NULL; }
    if (s_d_filt_a_pri != NULL)  { cudaFree(s_d_filt_a_pri);  s_d_filt_a_pri = NULL; }
    if (s_d_filt_b_pri != NULL)  { cudaFree(s_d_filt_b_pri);  s_d_filt_b_pri = NULL; }
    if (s_d_filt_a_sec != NULL)  { cudaFree(s_d_filt_a_sec);  s_d_filt_a_sec = NULL; }
    if (s_d_filt_b_sec != NULL)  { cudaFree(s_d_filt_b_sec);  s_d_filt_b_sec = NULL; }
    if (s_d_chrom_ramp != NULL)  { cudaFree(s_d_chrom_ramp);  s_d_chrom_ramp = NULL; }
    s_cached_filters = 0;

    if (s_d_super_pha != NULL)    { cudaFree(s_d_super_pha);    s_d_super_pha = NULL; }
    if (s_d_super_spha != NULL)   { cudaFree(s_d_super_spha);   s_d_super_spha = NULL; }
    if (s_d_bounds != NULL)       { cudaFree(s_d_bounds);       s_d_bounds = NULL; }
    if (s_d_hat_bounds != NULL)   { cudaFree(s_d_hat_bounds);   s_d_hat_bounds = NULL; }
    if (s_d_spec != NULL)         { cudaFree(s_d_spec);         s_d_spec = NULL; }
    if (s_d_acc_dphi_pri != NULL) { cudaFree(s_d_acc_dphi_pri); s_d_acc_dphi_pri = NULL; }
    if (s_d_acc_chi_pri != NULL)  { cudaFree(s_d_acc_chi_pri);  s_d_acc_chi_pri = NULL; }
    if (s_d_acc_dphi_sec != NULL) { cudaFree(s_d_acc_dphi_sec); s_d_acc_dphi_sec = NULL; }
    if (s_d_acc_chi_sec != NULL)  { cudaFree(s_d_acc_chi_sec);  s_d_acc_chi_sec = NULL; }
    if (s_d_dphi_out != NULL)     { cudaFree(s_d_dphi_out);     s_d_dphi_out = NULL; }
    if (s_d_chi_out != NULL)      { cudaFree(s_d_chi_out);      s_d_chi_out = NULL; }
    if (s_d_frame_pha != NULL)    { cudaFree(s_d_frame_pha);    s_d_frame_pha = NULL; }
    if (s_d_frame_amp != NULL)    { cudaFree(s_d_frame_amp);    s_d_frame_amp = NULL; }
    if (s_d_frame_spha != NULL)   { cudaFree(s_d_frame_spha);   s_d_frame_spha = NULL; }
    if (s_d_frame_samp != NULL)   { cudaFree(s_d_frame_samp);   s_d_frame_samp = NULL; }
    if (s_d_sum_I != NULL)        { cudaFree(s_d_sum_I);        s_d_sum_I = NULL; }
    s_cached_work = 0;
}

static int atmturb_cuda_rytov_init_plans(long pad_size)
{
    if (s_plan_1d != 0 && s_cached_pad == pad_size)
    {
        return 0;
    }
    atmturb_cuda_rytov_cleanup();

    if (cufftPlan1d(&s_plan_1d, (int) pad_size, CUFFT_R2C, 2) != CUFFT_SUCCESS ||
        cufftPlan2d(&s_plan_2d_r2c, (int) pad_size, (int) pad_size, CUFFT_R2C) != CUFFT_SUCCESS ||
        cufftPlan2d(&s_plan_2d_c2r, (int) pad_size, (int) pad_size, CUFFT_C2R) != CUFFT_SUCCESS)
    {
        atmturb_cuda_rytov_cleanup();
        return -1;
    }
    s_cached_pad = pad_size;
    return 0;
}

static int atmturb_cuda_rytov_init_filters(
    const atmturb_cuda_rytov_params_t *params)
{
    long pad = params->pad_size;
    long n_half = pad / 2 + 1;
    long n_spec = pad * n_half;
    int ns = params->nsuper;
    size_t f_bytes = (size_t) (ns * n_spec) * sizeof(float);

    if (s_d_filt_a_pri != NULL && s_cached_filters == f_bytes)
    {
        return 0;
    }

    cudaMalloc((void **) &s_d_laplace_inv, sizeof(float) * (size_t) n_spec);
    cudaMalloc((void **) &s_d_exp_x, sizeof(cufftComplex) * (size_t) n_half);
    cudaMalloc((void **) &s_d_exp_y, sizeof(cufftComplex) * (size_t) pad);
    cudaMalloc((void **) &s_d_filt_a_pri, f_bytes);
    cudaMalloc((void **) &s_d_filt_b_pri, f_bytes);
    cudaMalloc((void **) &s_d_filt_a_sec, f_bytes);
    cudaMalloc((void **) &s_d_filt_b_sec, f_bytes);
    cudaMalloc((void **) &s_d_chrom_ramp, (size_t) (ns * n_spec) * sizeof(cufftComplex));

    cudaMemcpy(s_d_laplace_inv, params->h_laplace_inv, sizeof(float) * (size_t) n_spec,
               cudaMemcpyHostToDevice);
    cudaMemcpy(s_d_exp_x, params->h_exp_x, sizeof(cufftComplex) * (size_t) n_half,
               cudaMemcpyHostToDevice);
    cudaMemcpy(s_d_exp_y, params->h_exp_y, sizeof(cufftComplex) * (size_t) pad,
               cudaMemcpyHostToDevice);

    for (int m = 0; m < ns; m++)
    {
        size_t off = (size_t) (m * n_spec);
        cudaMemcpy(s_d_filt_a_pri + off, params->h_filt_a_pri[m],
                   sizeof(float) * (size_t) n_spec, cudaMemcpyHostToDevice);
        cudaMemcpy(s_d_filt_b_pri + off, params->h_filt_b_pri[m],
                   sizeof(float) * (size_t) n_spec, cudaMemcpyHostToDevice);
        if (params->has_sec)
        {
            cudaMemcpy(s_d_filt_a_sec + off, params->h_filt_a_sec[m],
                       sizeof(float) * (size_t) n_spec, cudaMemcpyHostToDevice);
            cudaMemcpy(s_d_filt_b_sec + off, params->h_filt_b_sec[m],
                       sizeof(float) * (size_t) n_spec, cudaMemcpyHostToDevice);
            if (params->sec_shared && params->h_chrom_ramp != NULL &&
                params->h_chrom_ramp[m] != NULL)
            {
                cudaMemcpy(s_d_chrom_ramp + off, params->h_chrom_ramp[m],
                           sizeof(cufftComplex) * (size_t) n_spec, cudaMemcpyHostToDevice);
            }
        }
    }

    s_cached_filters = f_bytes;
    return 0;
}

static int atmturb_cuda_rytov_init_workspace(
    long pad_size,
    long pup_size)
{
    long pad_pixels = pad_size * pad_size;
    long pup_pixels = pup_size * pup_size;
    long n_half     = pad_size / 2 + 1;
    long n_spec     = pad_size * n_half;

    if (s_d_super_pha != NULL && s_cached_work == (size_t) pad_pixels)
    {
        return 0;
    }

    cudaMalloc((void **) &s_d_super_pha, sizeof(float) * (size_t) pad_pixels);
    cudaMalloc((void **) &s_d_super_spha, sizeof(float) * (size_t) pad_pixels);
    cudaMalloc((void **) &s_d_bounds, sizeof(float) * (size_t) (2 * pad_size));
    cudaMalloc((void **) &s_d_hat_bounds, sizeof(cufftComplex) * (size_t) (2 * n_half));
    cudaMalloc((void **) &s_d_spec, sizeof(cufftComplex) * (size_t) n_spec);
    cudaMalloc((void **) &s_d_acc_dphi_pri, sizeof(cufftComplex) * (size_t) n_spec);
    cudaMalloc((void **) &s_d_acc_chi_pri, sizeof(cufftComplex) * (size_t) n_spec);
    cudaMalloc((void **) &s_d_acc_dphi_sec, sizeof(cufftComplex) * (size_t) n_spec);
    cudaMalloc((void **) &s_d_acc_chi_sec, sizeof(cufftComplex) * (size_t) n_spec);
    cudaMalloc((void **) &s_d_dphi_out, sizeof(float) * (size_t) pad_pixels);
    cudaMalloc((void **) &s_d_chi_out, sizeof(float) * (size_t) pad_pixels);
    cudaMalloc((void **) &s_d_frame_pha, sizeof(float) * (size_t) pup_pixels);
    cudaMalloc((void **) &s_d_frame_amp, sizeof(float) * (size_t) pup_pixels);
    cudaMalloc((void **) &s_d_frame_spha, sizeof(float) * (size_t) pup_pixels);
    cudaMalloc((void **) &s_d_frame_samp, sizeof(float) * (size_t) pup_pixels);
    cudaMalloc((void **) &s_d_sum_I, sizeof(float));

    s_cached_work = (size_t) pad_pixels;
    return 0;
}

static void atmturb_cuda_rytov_render_superlayer(
    const atmturb_cuda_rytov_params_t *params,
    int                                m,
    int                                pad,
    int                                n_half,
    int                                n_spec)
{
    dim3 b_1d(256);
    dim3 g_1d(((int) pad + b_1d.x - 1) / b_1d.x);
    atmturb_cuda_extract_bounds_kernel<<<g_1d, b_1d>>>(s_d_super_pha, s_d_bounds, pad);

    cufftExecR2C(s_plan_1d, (cufftReal *) s_d_bounds, (cufftComplex *) s_d_hat_bounds);
    cufftExecR2C(s_plan_2d_r2c, (cufftReal *) s_d_super_pha, (cufftComplex *) s_d_spec);

    dim3 b_spec(16, 16);
    dim3 g_spec(((int) n_half + b_spec.x - 1) / b_spec.x,
                ((int) pad + b_spec.y - 1) / b_spec.y);
    atmturb_cuda_moisan_adjust_kernel<<<g_spec, b_spec>>>(
        s_d_spec, s_d_hat_bounds, s_d_exp_x, s_d_exp_y, s_d_laplace_inv, pad, n_half);

    dim3 b_mac(256);
    dim3 g_mac(((int) n_spec + b_mac.x - 1) / b_mac.x);
    size_t off = (size_t) (m * n_spec);
    atmturb_cuda_filter_mac_kernel<<<g_mac, b_mac>>>(
        s_d_acc_dphi_pri, s_d_acc_chi_pri, s_d_spec,
        s_d_filt_a_pri + off, s_d_filt_b_pri + off, n_spec);

    if (params->has_sec)
    {
        if (params->sec_shared && params->h_chrom_ramp != NULL &&
            params->h_chrom_ramp[m] != NULL)
        {
            atmturb_cuda_filter_mac_rotated_kernel<<<g_mac, b_mac>>>(
                s_d_acc_dphi_sec, s_d_acc_chi_sec, s_d_spec, s_d_chrom_ramp + off,
                s_d_filt_a_sec + off, s_d_filt_b_sec + off, n_spec);
        }
        else
        {
            atmturb_cuda_extract_bounds_kernel<<<g_1d, b_1d>>>(s_d_super_spha, s_d_bounds, pad);
            cufftExecR2C(s_plan_1d, (cufftReal *) s_d_bounds, (cufftComplex *) s_d_hat_bounds);
            cufftExecR2C(s_plan_2d_r2c, (cufftReal *) s_d_super_spha, (cufftComplex *) s_d_spec);
            atmturb_cuda_moisan_adjust_kernel<<<g_spec, b_spec>>>(
                s_d_spec, s_d_hat_bounds, s_d_exp_x, s_d_exp_y, s_d_laplace_inv, pad, n_half);
            atmturb_cuda_filter_mac_kernel<<<g_mac, b_mac>>>(
                s_d_acc_dphi_sec, s_d_acc_chi_sec, s_d_spec,
                s_d_filt_a_sec + off, s_d_filt_b_sec + off, n_spec);
        }
    }
}

int atmturb_cuda_rytov_render(
    const atmturb_cuda_rytov_params_t *params)
{
    if (params->nblayers <= 0 || params->pup_size <= 0 || params->nbframes <= 0)
    {
        return -1;
    }

    float *d_masters = atmturb_cuda_sync_device_masters(params->msize, params->nblayers,
                                                        params->h_masters);
    if (!d_masters)
    {
        return -1;
    }

    long pad        = params->pad_size;
    long guard      = params->guard_pix;
    long pup_size   = params->pup_size;
    long pup_pixels = pup_size * pup_size;
    long n_half     = pad / 2 + 1;
    long n_spec     = pad * n_half;
    int  has_sec    = params->has_sec;
    int  os         = (params->os > 1) ? params->os : 1;

    if (atmturb_cuda_rytov_init_plans(pad) != 0 ||
        atmturb_cuda_rytov_init_filters(params) != 0 ||
        atmturb_cuda_rytov_init_workspace(pad, pup_size) != 0)
    {
        return -1;
    }

    dim3 b_ext(16, 16);
    dim3 g_ext(((int) pad + b_ext.x - 1) / b_ext.x, ((int) pad + b_ext.y - 1) / b_ext.y);
    dim3 b_ass(16, 16);
    dim3 g_ass(((int) pup_size + b_ass.x - 1) / b_ass.x,
               ((int) pup_size + b_ass.y - 1) / b_ass.y);
    dim3 b_norm(256);
    dim3 g_norm(((int) pup_pixels + b_norm.x - 1) / b_norm.x);

    for (long t = 0; t < params->nbframes; t++)
    {
        cudaMemset(s_d_frame_pha, 0, sizeof(float) * (size_t) pup_pixels);
        cudaMemset(s_d_frame_spha, 0, sizeof(float) * (size_t) pup_pixels);
        cudaMemset(s_d_acc_dphi_pri, 0, sizeof(cufftComplex) * (size_t) n_spec);
        cudaMemset(s_d_acc_chi_pri, 0, sizeof(cufftComplex) * (size_t) n_spec);
        if (has_sec)
        {
            cudaMemset(s_d_acc_dphi_sec, 0, sizeof(cufftComplex) * (size_t) n_spec);
            cudaMemset(s_d_acc_chi_sec, 0, sizeof(cufftComplex) * (size_t) n_spec);
        }

        for (int m = 0; m < params->nsuper; m++)
        {
            int n_sub = params->super_nlayers[m];
            size_t sub_idx = (size_t) ((t * params->nsuper + m) * params->max_sublayers);
            cudaMemcpyToSymbol(c_sublayers, params->frame_sublayers + sub_idx,
                               sizeof(atmturb_cuda_rytov_sublayer_t) * (size_t) n_sub);

            if (params->interp == 1)
            {
                atmturb_cuda_extrude_sl_kernel<1><<<g_ext, b_ext>>>(
                    d_masters, (int) params->msize, (int) pad, (int) pup_size, (int) guard,
                    n_sub, has_sec, os, s_d_super_pha, s_d_super_spha,
                    s_d_frame_pha, s_d_frame_spha);
            }
            else
            {
                atmturb_cuda_extrude_sl_kernel<0><<<g_ext, b_ext>>>(
                    d_masters, (int) params->msize, (int) pad, (int) pup_size, (int) guard,
                    n_sub, has_sec, os, s_d_super_pha, s_d_super_spha,
                    s_d_frame_pha, s_d_frame_spha);
            }

            if (params->super_dist_m[m] > 0.0)
            {
                atmturb_cuda_rytov_render_superlayer(params, m, (int) pad, (int) n_half,
                                                     (int) n_spec);
            }
        }

        cufftExecC2R(s_plan_2d_c2r, (cufftComplex *) s_d_acc_dphi_pri, (cufftReal *) s_d_dphi_out);
        cufftExecC2R(s_plan_2d_c2r, (cufftComplex *) s_d_acc_chi_pri, (cufftReal *) s_d_chi_out);

        cudaMemset(s_d_sum_I, 0, sizeof(float));
        atmturb_cuda_assemble_kernel<<<g_ass, b_ass>>>(
            (int) pup_size, (int) guard, (int) pad, s_d_frame_pha, s_d_frame_amp,
            s_d_dphi_out, s_d_chi_out, s_d_sum_I);
        atmturb_cuda_normalize_kernel<<<g_norm, b_norm>>>(
            s_d_frame_amp, (int) pup_pixels, s_d_sum_I);

        size_t slice_off = (size_t) (t * pup_pixels);
        cudaMemcpy(params->pha + slice_off, s_d_frame_pha,
                   sizeof(float) * (size_t) pup_pixels, cudaMemcpyDeviceToHost);
        cudaMemcpy(params->amp + slice_off, s_d_frame_amp,
                   sizeof(float) * (size_t) pup_pixels, cudaMemcpyDeviceToHost);

        if (has_sec)
        {
            cufftExecC2R(s_plan_2d_c2r, (cufftComplex *) s_d_acc_dphi_sec,
                         (cufftReal *) s_d_dphi_out);
            cufftExecC2R(s_plan_2d_c2r, (cufftComplex *) s_d_acc_chi_sec,
                         (cufftReal *) s_d_chi_out);

            cudaMemset(s_d_sum_I, 0, sizeof(float));
            atmturb_cuda_assemble_kernel<<<g_ass, b_ass>>>(
                (int) pup_size, (int) guard, (int) pad, s_d_frame_spha, s_d_frame_samp,
                s_d_dphi_out, s_d_chi_out, s_d_sum_I);
            atmturb_cuda_normalize_kernel<<<g_norm, b_norm>>>(
                s_d_frame_samp, (int) pup_pixels, s_d_sum_I);

            cudaMemcpy(params->spha + slice_off, s_d_frame_spha,
                       sizeof(float) * (size_t) pup_pixels, cudaMemcpyDeviceToHost);
            cudaMemcpy(params->samp + slice_off, s_d_frame_samp,
                       sizeof(float) * (size_t) pup_pixels, cudaMemcpyDeviceToHost);
        }
    }

    return 0;
}
