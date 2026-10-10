// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_rytov_kernels.cu
 * @brief   GPU device kernels for multi-layer Rytov wavefront synthesis
 */

#include <cuda_runtime.h>
#include <cufft.h>
#include <math.h>

#include "atmturb_cuda_rytov_internal.h"

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
            if (xm < 0) xm += msize;
            ix_arr[m] = xm;

            int ym = (iy + m - 1) % msize;
            if (ym < 0) ym += msize;
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
    const float                         *d_masters,
    int                                  msize,
    int                                  pad_size,
    int                                  pup_size,
    int                                  guard,
    int                                  n_sub,
    int                                  has_sec,
    int                                  os,
    const atmturb_cuda_rytov_sublayer_t *d_sublayers,
    float                               *d_super_pha,
    float                               *d_super_spha,
    float                               *d_frame_pha,
    float                               *d_frame_spha)
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
        atmturb_cuda_rytov_sublayer_t info = d_sublayers[l];
        float cur_x  = info.x + (float) (i * os);
        float cur_y  = info.y + (float) (j * os);
        if (info.w_pri != 0.0f)
        {
            if (interp_mode == 1)
            {
                p_val += info.w_pri * atmturb_cuda_sample_screen_bicubic(
                    d_masters, msize, is_pow2, mask, info.k, cur_x, cur_y);
            }
            else
            {
                p_val += info.w_pri * atmturb_cuda_sample_screen(
                    d_masters, msize, is_pow2, mask, info.k, cur_x, cur_y);
            }
        }
        if (has_sec && info.w_sec != 0.0f)
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

void atmturb_cuda_extrude_sl_launch(
    int                                  interp_mode,
    const float                         *d_masters,
    int                                  msize,
    int                                  pad_size,
    int                                  pup_size,
    int                                  guard,
    int                                  n_sub,
    int                                  has_sec,
    int                                  os,
    const atmturb_cuda_rytov_sublayer_t *d_sublayers,
    float                               *d_super_pha,
    float                               *d_super_spha,
    float                               *d_frame_pha,
    float                               *d_frame_spha,
    cudaStream_t                         stream)
{
    dim3 b_ext(16, 16);
    dim3 g_ext(((int) pad_size + b_ext.x - 1) / b_ext.x,
               ((int) pad_size + b_ext.y - 1) / b_ext.y);
    if (interp_mode == 1)
    {
        atmturb_cuda_extrude_sl_kernel<1><<<g_ext, b_ext, 0, stream>>>(
            d_masters, msize, pad_size, pup_size, guard, n_sub, has_sec, os,
            d_sublayers, d_super_pha, d_super_spha, d_frame_pha, d_frame_spha);
    }
    else
    {
        atmturb_cuda_extrude_sl_kernel<0><<<g_ext, b_ext, 0, stream>>>(
            d_masters, msize, pad_size, pup_size, guard, n_sub, has_sec, os,
            d_sublayers, d_super_pha, d_super_spha, d_frame_pha, d_frame_spha);
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

void atmturb_cuda_extract_bounds_launch(
    const float *d_super_pha,
    float       *d_bounds,
    int          pad_size,
    cudaStream_t stream)
{
    dim3 b_1d(256);
    dim3 g_1d(((int) pad_size + b_1d.x - 1) / b_1d.x);
    atmturb_cuda_extract_bounds_kernel<<<g_1d, b_1d, 0, stream>>>(
        d_super_pha, d_bounds, pad_size);
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
    cufftComplex ha = hat_a[i];

    float v_re = (ha.x * ey.x - ha.y * ey.y) + (b_re * ex.x - b_im * ex.y);
    float v_im = (ha.x * ey.y + ha.y * ey.x) + (b_re * ex.y + b_im * ex.x);

    int idx = j * n_half + i;
    float linv = d_laplace_inv[idx];
    d_spec[idx].x += v_re * linv;
    d_spec[idx].y += v_im * linv;
}

void atmturb_cuda_moisan_adjust_launch(
    cufftComplex       *d_spec,
    const cufftComplex *d_hat_bounds,
    const cufftComplex *d_exp_x,
    const cufftComplex *d_exp_y,
    const float        *d_laplace_inv,
    int                 pad_size,
    int                 n_half,
    cudaStream_t        stream)
{
    dim3 b_spec(16, 16);
    dim3 g_spec(((int) n_half + b_spec.x - 1) / b_spec.x,
                ((int) pad_size + b_spec.y - 1) / b_spec.y);
    atmturb_cuda_moisan_adjust_kernel<<<g_spec, b_spec, 0, stream>>>(
        d_spec, d_hat_bounds, d_exp_x, d_exp_y, d_laplace_inv, pad_size, n_half);
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

void atmturb_cuda_filter_mac_launch(
    cufftComplex       *d_acc_dphi,
    cufftComplex       *d_acc_chi,
    const cufftComplex *d_spec,
    const float        *d_filt_a,
    const float        *d_filt_b,
    int                 ntot,
    cudaStream_t        stream)
{
    dim3 b_mac(256);
    dim3 g_mac(((int) ntot + b_mac.x - 1) / b_mac.x);
    atmturb_cuda_filter_mac_kernel<<<g_mac, b_mac, 0, stream>>>(
        d_acc_dphi, d_acc_chi, d_spec, d_filt_a, d_filt_b, ntot);
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

void atmturb_cuda_filter_mac_rotated_launch(
    cufftComplex       *d_acc_dphi_sec,
    cufftComplex       *d_acc_chi_sec,
    const cufftComplex *d_spec,
    const cufftComplex *d_ramp,
    const float        *d_filt_a_sec,
    const float        *d_filt_b_sec,
    int                 ntot,
    cudaStream_t        stream)
{
    dim3 b_mac(256);
    dim3 g_mac(((int) ntot + b_mac.x - 1) / b_mac.x);
    atmturb_cuda_filter_mac_rotated_kernel<<<g_mac, b_mac, 0, stream>>>(
        d_acc_dphi_sec, d_acc_chi_sec, d_spec, d_ramp,
        d_filt_a_sec, d_filt_b_sec, ntot);
}

__global__ static void atmturb_cuda_filter_mac_dual_kernel(
    cufftComplex       *d_acc_dphi_pri,
    cufftComplex       *d_acc_chi_pri,
    cufftComplex       *d_acc_dphi_sec,
    cufftComplex       *d_acc_chi_sec,
    const cufftComplex *d_spec,
    const cufftComplex *d_ramp,
    const float        *d_filt_a_pri,
    const float        *d_filt_b_pri,
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

    float fa_pri = d_filt_a_pri[k];
    float fb_pri = d_filt_b_pri[k];
    d_acc_dphi_pri[k].x += s.x * fa_pri;
    d_acc_dphi_pri[k].y += s.y * fa_pri;
    d_acc_chi_pri[k].x  += s.x * fb_pri;
    d_acc_chi_pri[k].y  += s.y * fb_pri;

    cufftComplex r = d_ramp[k];
    float ps_x = s.x * r.x - s.y * r.y;
    float ps_y = s.x * r.y + s.y * r.x;

    float fa_sec = d_filt_a_sec[k];
    float fb_sec = d_filt_b_sec[k];
    d_acc_dphi_sec[k].x += ps_x * fa_sec;
    d_acc_dphi_sec[k].y += ps_y * fa_sec;
    d_acc_chi_sec[k].x  += ps_x * fb_sec;
    d_acc_chi_sec[k].y  += ps_y * fb_sec;
}

void atmturb_cuda_filter_mac_dual_launch(
    cufftComplex       *d_acc_dphi_pri,
    cufftComplex       *d_acc_chi_pri,
    cufftComplex       *d_acc_dphi_sec,
    cufftComplex       *d_acc_chi_sec,
    const cufftComplex *d_spec,
    const cufftComplex *d_ramp,
    const float        *d_filt_a_pri,
    const float        *d_filt_b_pri,
    const float        *d_filt_a_sec,
    const float        *d_filt_b_sec,
    int                 ntot,
    cudaStream_t        stream)
{
    dim3 b_mac(256);
    dim3 g_mac(((int) ntot + b_mac.x - 1) / b_mac.x);
    atmturb_cuda_filter_mac_dual_kernel<<<g_mac, b_mac, 0, stream>>>(
        d_acc_dphi_pri, d_acc_chi_pri, d_acc_dphi_sec, d_acc_chi_sec,
        d_spec, d_ramp, d_filt_a_pri, d_filt_b_pri, d_filt_a_sec, d_filt_b_sec, ntot);
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

void atmturb_cuda_assemble_launch(
    int          pup_size,
    int          guard,
    int          pad_size,
    float       *d_frame_pha,
    float       *d_frame_amp,
    const float *d_dphi_out,
    const float *d_chi_out,
    float       *d_sum_I,
    cudaStream_t stream)
{
    dim3 b_ass(16, 16);
    dim3 g_ass(((int) pup_size + b_ass.x - 1) / b_ass.x,
               ((int) pup_size + b_ass.y - 1) / b_ass.y);
    atmturb_cuda_assemble_kernel<<<g_ass, b_ass, 0, stream>>>(
        pup_size, guard, pad_size, d_frame_pha, d_frame_amp,
        d_dphi_out, d_chi_out, d_sum_I);
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
    float norm = (sum_i > 0.0f) ? rsqrtf(sum_i / (float) npix) : 1.0f;
    d_frame_amp[idx] *= norm;
}

void atmturb_cuda_normalize_launch(
    float       *d_frame_amp,
    int          npix,
    const float *d_sum_I,
    cudaStream_t stream)
{
    dim3 b_norm(256);
    dim3 g_norm(((int) npix + b_norm.x - 1) / b_norm.x);
    atmturb_cuda_normalize_kernel<<<g_norm, b_norm, 0, stream>>>(
        d_frame_amp, npix, d_sum_I);
}
