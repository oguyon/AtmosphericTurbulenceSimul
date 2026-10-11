// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_rytov_extrude.cu
 * @brief   GPU device kernels for superlayer extrusion and subharmonic phase synthesis
 */

#include <cuda_runtime.h>
#include <math.h>

#include "atmturb_cuda_rytov_internal.h"

/**
 * atmturb_cuda_sample_screen - Bilinear interpolation sample from toroidal master screen
 * @d_masters: Array of device pointers to master phase screens.
 * @msize: Master screen dimension in pixels.
 * @is_pow2: 1 if msize is a power of 2, 0 otherwise.
 * @mask: Bitmask for power-of-2 wrapping (msize - 1).
 * @k: Master screen layer index.
 * @cur_x: Fractional sampling coordinate along X in master pixels.
 * @cur_y: Fractional sampling coordinate along Y in master pixels.
 *
 * Evaluates bilinear interpolation with periodic toroidal boundary conditions.
 *
 * Return: Interpolated phase value in radians.
 */
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

/**
 * atmturb_cuda_sample_screen_bicubic - Keys bicubic interpolation from toroidal master screen
 * @d_masters: Array of device pointers to master phase screens.
 * @msize: Master screen dimension in pixels.
 * @is_pow2: 1 if msize is a power of 2, 0 otherwise.
 * @mask: Bitmask for power-of-2 wrapping (msize - 1).
 * @k: Master screen layer index.
 * @cur_x: Fractional sampling coordinate along X in master pixels.
 * @cur_y: Fractional sampling coordinate along Y in master pixels.
 *
 * Evaluates Catmull-Rom / Keys bicubic convolution (alpha = -0.5) over a 4x4 pixel support.
 *
 * Return: Interpolated phase value in radians.
 */
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
    float wx0 = -0.5f * fx + fx2 - 0.5f * fx3;
    float wx1 = 1.0f - 2.5f * fx2 + 1.5f * fx3;
    float wx2 = 0.5f * fx + 2.0f * fx2 - 1.5f * fx3;
    float wx3 = -0.5f * fx2 + 0.5f * fx3;

    float fy2 = fy * fy;
    float fy3 = fy2 * fy;
    float wy0 = -0.5f * fy + fy2 - 0.5f * fy3;
    float wy1 = 1.0f - 2.5f * fy2 + 1.5f * fy3;
    float wy2 = 0.5f * fy + 2.0f * fy2 - 1.5f * fy3;
    float wy3 = -0.5f * fy2 + 0.5f * fy3;

    int ix0, ix1, ix2, ix3;
    int iy0, iy1, iy2, iy3;
    if (is_pow2)
    {
        ix0 = (ix - 1) & mask;
        ix1 = ix       & mask;
        ix2 = (ix + 1) & mask;
        ix3 = (ix + 2) & mask;
        iy0 = (iy - 1) & mask;
        iy1 = iy       & mask;
        iy2 = (iy + 1) & mask;
        iy3 = (iy + 2) & mask;
    }
    else
    {
        ix0 = (ix - 1) % msize; if (ix0 < 0) ix0 += msize;
        ix1 = ix       % msize; if (ix1 < 0) ix1 += msize;
        ix2 = (ix + 1) % msize; if (ix2 < 0) ix2 += msize;
        ix3 = (ix + 2) % msize; if (ix3 < 0) ix3 += msize;
        iy0 = (iy - 1) % msize; if (iy0 < 0) iy0 += msize;
        iy1 = iy       % msize; if (iy1 < 0) iy1 += msize;
        iy2 = (iy + 1) % msize; if (iy2 < 0) iy2 += msize;
        iy3 = (iy + 2) % msize; if (iy3 < 0) iy3 += msize;
    }

    size_t msize_sq = (size_t) msize * (size_t) msize;
    const float *screen = d_masters + (size_t) k * msize_sq;
    const float *row0 = screen + (size_t) iy0 * (size_t) msize;
    const float *row1 = screen + (size_t) iy1 * (size_t) msize;
    const float *row2 = screen + (size_t) iy2 * (size_t) msize;
    const float *row3 = screen + (size_t) iy3 * (size_t) msize;

    float h0 = wx0 * __ldg(&row0[ix0]) + wx1 * __ldg(&row0[ix1]) +
               wx2 * __ldg(&row0[ix2]) + wx3 * __ldg(&row0[ix3]);
    float h1 = wx0 * __ldg(&row1[ix0]) + wx1 * __ldg(&row1[ix1]) +
               wx2 * __ldg(&row1[ix2]) + wx3 * __ldg(&row1[ix3]);
    float h2 = wx0 * __ldg(&row2[ix0]) + wx1 * __ldg(&row2[ix1]) +
               wx2 * __ldg(&row2[ix2]) + wx3 * __ldg(&row2[ix3]);
    float h3 = wx0 * __ldg(&row3[ix0]) + wx1 * __ldg(&row3[ix1]) +
               wx2 * __ldg(&row3[ix2]) + wx3 * __ldg(&row3[ix3]);

    return wy0 * h0 + wy1 * h1 + wy2 * h2 + wy3 * h3;
}

/**
 * atmturb_cuda_accumulate_lowfreq_shared - Shared-memory batched accumulation of subharmonics
 * @n_modes: Total number of low-frequency Fourier modes.
 * @has_sec: 1 if secondary chromatic wavelength is enabled, 0 otherwise.
 * @d_modes: Array of device low-frequency mode definitions.
 * @i: Grid pixel column index.
 * @j: Grid pixel row index.
 * @pad_size: Linear dimension of padded simulation grid.
 * @acc_pri_out: Destination accumulator for primary wavelength phase.
 * @acc_sec_out: Destination accumulator for secondary wavelength phase.
 *
 * Utilizes 2D threadblock shared memory tiles to evaluate complex Fourier basis exponentials
 * across thread ranks with minimal redundant trig evaluations.
 */
__device__ static inline void atmturb_cuda_accumulate_lowfreq_shared(
    int                                n_modes,
    int                                has_sec,
    const atmturb_cuda_lowfreq_mode_t *d_modes,
    int                                i,
    int                                j,
    int                                pad_size,
    float                             *acc_pri_out,
    float                             *acc_sec_out)
{
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    int valid = (i < pad_size && j < pad_size);

    __shared__ float s_X0[4][16];
    __shared__ float s_X1[4][16];
    __shared__ float s_Yc[4][16];
    __shared__ float s_Ys[4][16];
    __shared__ float s_X0_sec[4][16];
    __shared__ float s_X1_sec[4][16];

    float acc_pri = 0.0f;
    float acc_sec = 0.0f;

    for (int base_m = 0; base_m < n_modes; base_m += 4)
    {
        int count = n_modes - base_m;
        if (count > 4)
        {
            count = 4;
        }

        if (ty < count)
        {
            atmturb_cuda_lowfreq_mode_t mode = d_modes[base_m + ty];
            float ca, sa;
            __sincosf(mode.th0_pri + (float) i * mode.dkx, &sa, &ca);
            s_X0[ty][tx] = mode.are_pri * ca - mode.aim_pri * sa;
            s_X1[ty][tx] = mode.are_pri * sa + mode.aim_pri * ca;

            if (has_sec)
            {
                if (mode.th0_sec == mode.th0_pri)
                {
                    s_X0_sec[ty][tx] = mode.are_sec * ca - mode.aim_sec * sa;
                    s_X1_sec[ty][tx] = mode.are_sec * sa + mode.aim_sec * ca;
                }
                else
                {
                    float ca_s, sa_s;
                    __sincosf(mode.th0_sec + (float) i * mode.dkx, &sa_s, &ca_s);
                    s_X0_sec[ty][tx] = mode.are_sec * ca_s - mode.aim_sec * sa_s;
                    s_X1_sec[ty][tx] = mode.are_sec * sa_s + mode.aim_sec * ca_s;
                }
            }
        }

        if (tx < count)
        {
            atmturb_cuda_lowfreq_mode_t mode = d_modes[base_m + tx];
            float cb, sb;
            __sincosf((float) j * mode.dky, &sb, &cb);
            s_Yc[tx][ty] = cb;
            s_Ys[tx][ty] = sb;
        }

        __syncthreads();

        if (valid)
        {
            #pragma unroll
            for (int mb = 0; mb < count; mb++)
            {
                float yc = s_Yc[mb][ty];
                float ys = s_Ys[mb][ty];
                acc_pri += s_X0[mb][tx] * yc - s_X1[mb][tx] * ys;
                if (has_sec)
                {
                    acc_sec += s_X0_sec[mb][tx] * yc - s_X1_sec[mb][tx] * ys;
                }
            }
        }

        __syncthreads();
    }

    *acc_pri_out += acc_pri;
    *acc_sec_out += acc_sec;
}

/**
 * atmturb_cuda_extrude_sl_kernel - Device kernel extruding sublayers onto padded grid
 * @d_masters: Array of device pointers to master phase screens.
 * @msize: Master phase screen linear dimension in pixels.
 * @pad_size: Linear dimension of padded simulation grid in pixels.
 * @pup_size: Linear dimension of synthesized pupil in pixels.
 * @guard: Guard band margin in pixels surrounding pupil.
 * @n_sub: Number of active sublayers in current superlayer.
 * @has_sec: 1 if secondary wavelength is enabled, 0 otherwise.
 * @os: Master screen oversampling factor (1 or 2).
 * @d_sublayers: Array of device sublayer coordinates and interpolation weights.
 * @n_lf_modes: Number of active low-frequency subharmonic modes.
 * @d_lf_modes: Array of device subharmonic mode parameters.
 * @d_super_pha: Destination superlayer primary phase array.
 * @d_super_spha: Destination superlayer secondary phase array.
 * @d_frame_pha: Cumulative pupil primary phase array.
 * @d_frame_spha: Cumulative pupil secondary phase array.
 */
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
    int                                  n_lf_modes,
    const atmturb_cuda_lowfreq_mode_t   *d_lf_modes,
    float                               *d_super_pha,
    float                               *d_super_spha,
    float                               *d_frame_pha,
    float                               *d_frame_spha)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int valid = (i < pad_size && j < pad_size);

    int is_pow2 = ((msize & (msize - 1)) == 0);
    int mask    = msize - 1;

    float p_val  = 0.0f;
    float sp_val = 0.0f;

    if (valid)
    {
        for (int l = 0; l < n_sub; l++)
        {
            atmturb_cuda_rytov_sublayer_t info = d_sublayers[l];
            float cur_x = info.x + (float) (i * os);
            float cur_y = info.y + (float) (j * os);

            if (has_sec && info.w_sec != 0.0f && info.xs == info.x && info.ys == info.y)
            {
                float samp = (interp_mode == 1)
                    ? atmturb_cuda_sample_screen_bicubic(
                        d_masters, msize, is_pow2, mask, info.k, cur_x, cur_y)
                    : atmturb_cuda_sample_screen(
                        d_masters, msize, is_pow2, mask, info.k, cur_x, cur_y);
                if (info.w_pri != 0.0f)
                {
                    p_val += info.w_pri * samp;
                }
                sp_val += info.w_sec * samp;
            }
            else
            {
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
        }
    }

    if (n_lf_modes > 0 && d_lf_modes != NULL)
    {
        atmturb_cuda_accumulate_lowfreq_shared(
            n_lf_modes, has_sec, d_lf_modes, i, j, pad_size, &p_val, &sp_val);
    }

    if (!valid)
    {
        return;
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

/**
 * atmturb_cuda_extrude_sl_launch - Launch superlayer phase screen extrusion and subharmonics
 * @interp_mode: Interpolation scheme (0=bilinear, 1=Keys bicubic).
 * @d_masters: Array of device pointers to master phase screens.
 * @msize: Master phase screen linear dimension in pixels.
 * @pad_size: Linear dimension of padded simulation grid in pixels.
 * @pup_size: Linear dimension of synthesized pupil in pixels.
 * @guard: Guard band margin in pixels surrounding pupil.
 * @n_sub: Number of active sublayers in current superlayer.
 * @has_sec: 1 if secondary wavelength is enabled, 0 otherwise.
 * @os: Master screen oversampling factor (1 or 2).
 * @d_sublayers: Array of device sublayer coordinates and interpolation weights.
 * @n_lf_modes: Number of active low-frequency subharmonic modes.
 * @d_lf_modes: Array of device subharmonic mode parameters.
 * @d_super_pha: Destination superlayer primary phase array.
 * @d_super_spha: Destination superlayer secondary phase array.
 * @d_frame_pha: Cumulative pupil primary phase array.
 * @d_frame_spha: Cumulative pupil secondary phase array.
 * @stream: CUDA execution stream.
 */
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
    int                                  n_lf_modes,
    const atmturb_cuda_lowfreq_mode_t   *d_lf_modes,
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
            d_sublayers, n_lf_modes, d_lf_modes, d_super_pha, d_super_spha,
            d_frame_pha, d_frame_spha);
    }
    else
    {
        atmturb_cuda_extrude_sl_kernel<0><<<g_ext, b_ext, 0, stream>>>(
            d_masters, msize, pad_size, pup_size, guard, n_sub, has_sec, os,
            d_sublayers, n_lf_modes, d_lf_modes, d_super_pha, d_super_spha,
            d_frame_pha, d_frame_spha);
    }
}

/**
 * atmturb_cuda_lowfreq_accumulate_kernel - Device kernel for subharmonic low-order modes
 * @pad_size: Linear dimension of padded simulation grid.
 * @pup_size: Linear dimension of synthesized pupil.
 * @guard: Guard margin in pixels.
 * @n_modes: Number of active low-frequency modes to evaluate.
 * @has_sec: 1 if secondary wavelength enabled, 0 otherwise.
 * @d_modes: Array of device low-frequency mode parameters.
 * @d_super_pha: Destination superlayer primary phase array.
 * @d_super_spha: Destination superlayer secondary phase array.
 * @d_frame_pha: Cumulative pupil primary phase array.
 * @d_frame_spha: Cumulative pupil secondary phase array.
 */
__global__ static void atmturb_cuda_lowfreq_accumulate_kernel(
    int                                pad_size,
    int                                pup_size,
    int                                guard,
    int                                n_modes,
    int                                has_sec,
    const atmturb_cuda_lowfreq_mode_t *d_modes,
    float                             *d_super_pha,
    float                             *d_super_spha,
    float                             *d_frame_pha,
    float                             *d_frame_spha)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    float acc_pri = 0.0f;
    float acc_sec = 0.0f;

    atmturb_cuda_accumulate_lowfreq_shared(
        n_modes, has_sec, d_modes, i, j, pad_size, &acc_pri, &acc_sec);

    if (i >= pad_size || j >= pad_size)
    {
        return;
    }

    int pad_idx = j * pad_size + i;
    d_super_pha[pad_idx] += acc_pri;
    if (has_sec && d_super_spha != NULL)
    {
        d_super_spha[pad_idx] += acc_sec;
    }

    if (i >= guard && i < guard + pup_size && j >= guard && j < guard + pup_size)
    {
        int pup_idx = (j - guard) * pup_size + (i - guard);
        if (d_frame_pha != NULL)
        {
            d_frame_pha[pup_idx] += acc_pri;
        }
        if (has_sec && d_frame_spha != NULL)
        {
            d_frame_spha[pup_idx] += acc_sec;
        }
    }
}

/**
 * atmturb_cuda_lowfreq_accumulate_launch - Launch subharmonic low-order accumulation
 * @pad_size: Linear dimension of padded simulation grid.
 * @pup_size: Linear dimension of synthesized pupil.
 * @guard: Guard margin in pixels.
 * @n_modes: Number of active low-frequency modes to evaluate.
 * @has_sec: 1 if secondary wavelength enabled, 0 otherwise.
 * @d_modes: Array of device low-frequency mode parameters.
 * @d_super_pha: Destination superlayer primary phase array.
 * @d_super_spha: Destination superlayer secondary phase array.
 * @d_frame_pha: Cumulative pupil primary phase array.
 * @d_frame_spha: Cumulative pupil secondary phase array.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_lowfreq_accumulate_launch(
    int                                pad_size,
    int                                pup_size,
    int                                guard,
    int                                n_modes,
    int                                has_sec,
    const atmturb_cuda_lowfreq_mode_t *d_modes,
    float                             *d_super_pha,
    float                             *d_super_spha,
    float                             *d_frame_pha,
    float                             *d_frame_spha,
    cudaStream_t                       stream)
{
    if (n_modes <= 0 || d_modes == NULL)
    {
        return;
    }

    dim3 b_lf(16, 16);
    dim3 g_lf(((int) pad_size + b_lf.x - 1) / b_lf.x,
              ((int) pad_size + b_lf.y - 1) / b_lf.y);
    atmturb_cuda_lowfreq_accumulate_kernel<<<g_lf, b_lf, 0, stream>>>(
        pad_size, pup_size, guard, n_modes, has_sec, d_modes,
        d_super_pha, d_super_spha, d_frame_pha, d_frame_spha);
}
