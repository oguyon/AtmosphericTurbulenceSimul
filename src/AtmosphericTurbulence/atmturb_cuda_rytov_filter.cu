// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_rytov_filter.cu
 * @brief   GPU device kernels for Moisan decomposition, Fourier filtering, and pupil assembly
 */

#include <cuda_runtime.h>
#include <cufft.h>
#include <math.h>

#include "atmturb_cuda_rytov_internal.h"

/**
 * atmturb_cuda_extract_bounds_kernel - Extract boundary differences for Moisan decomposition
 * @d_super_pha: Input superlayer phase screen on padded grid.
 * @d_bounds: Output boundary gradient array (row 0: top-bottom diff, row 1: left-right diff).
 * @pad_size: Linear dimension of padded grid in pixels.
 *
 * Evaluates boundary discontinuities between opposing edges of the extruded phase screen:
 * delta_top_bottom(x) = phi(x, 0) - phi(x, N-1)
 * delta_left_right(y) = phi(0, y) - phi(N-1, y)
 * These form the source term for Poisson's equation to decompose the screen into periodic
 * and smooth components without windowing artifacts.
 */
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

/**
 * atmturb_cuda_extract_bounds_launch - Dispatch 1D boundary difference extraction on GPU
 * @d_super_pha: Input superlayer phase screen on padded grid.
 * @d_bounds: Output boundary gradient array.
 * @pad_size: Linear dimension of padded grid in pixels.
 * @stream: CUDA execution stream.
 */
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

/**
 * atmturb_cuda_moisan_calc - Device helper computing Poisson boundary correction in frequency space
 * @s: Input complex Fourier spectrum coefficient at frequency (i, j).
 * @hat_bounds: 1D Fourier transforms of boundary differences.
 * @exp_x: Precomputed complex exponential roots along X.
 * @exp_y: Precomputed complex exponential roots along Y.
 * @laplace_inv: Inverse discrete Laplacian grid multiplier array.
 * @i: Spatial frequency index along X [0, n_half - 1].
 * @j: Spatial frequency index along Y [0, pad_size - 1].
 * @pad_size: Padded grid linear dimension.
 * @n_half: Half-complex linear dimension along X (pad_size/2 + 1).
 *
 * Implements Lionel Moisan's periodic plus smooth decomposition in Fourier space:
 * P(u, v) = S(u, v) - V(u, v) / (2 * cos(2*pi*u/N) + 2 * cos(2*pi*v/N) - 4).
 *
 * Return: Corrected periodic component Fourier coefficient.
 */
__device__ static inline cufftComplex atmturb_cuda_moisan_calc(
    cufftComplex       s,
    const cufftComplex *hat_bounds,
    const cufftComplex *exp_x,
    const cufftComplex *exp_y,
    const float        *laplace_inv,
    int                 i,
    int                 j,
    int                 pad_size,
    int                 n_half)
{
    const cufftComplex *hat_a = hat_bounds;
    const cufftComplex *hat_b = hat_bounds + n_half;

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

    cufftComplex ey = exp_y[j];
    cufftComplex ex = exp_x[i];
    cufftComplex ha = hat_a[i];

    float v_re = (ha.x * ey.x - ha.y * ey.y) + (b_re * ex.x - b_im * ex.y);
    float v_im = (ha.x * ey.y + ha.y * ey.x) + (b_re * ex.y + b_im * ex.x);

    int idx = j * n_half + i;
    float linv = laplace_inv[idx];
    s.x += v_re * linv;
    s.y += v_im * linv;
    return s;
}

/**
 * atmturb_cuda_moisan_adjust_kernel - Device kernel for in-place Moisan boundary adjustment
 * @d_spec: In-place 2D Fourier spectrum of phase screen.
 * @d_hat_bounds: 1D Fourier transforms of boundary row/col differences.
 * @d_exp_x: Precomputed complex roots along X.
 * @d_exp_y: Precomputed complex roots along Y.
 * @d_laplace_inv: Inverse discrete Laplacian grid values.
 * @pad_size: Linear dimension of padded grid.
 * @n_half: Half-complex dimension along X.
 */
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

    int idx = j * n_half + i;
    d_spec[idx] = atmturb_cuda_moisan_calc(
        d_spec[idx], d_hat_bounds, d_exp_x, d_exp_y, d_laplace_inv, i, j, pad_size, n_half);
}

/**
 * atmturb_cuda_moisan_adjust_launch - Launch in-place Moisan boundary correction on GPU
 * @d_spec: In-place 2D Fourier spectrum of phase screen.
 * @d_hat_bounds: 1D Fourier transforms of boundary row/col differences.
 * @d_exp_x: Precomputed complex roots along X.
 * @d_exp_y: Precomputed complex roots along Y.
 * @d_laplace_inv: Inverse discrete Laplacian grid values.
 * @pad_size: Linear dimension of padded grid.
 * @n_half: Half-complex dimension along X.
 * @stream: CUDA execution stream.
 */
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

/**
 * atmturb_cuda_moisan_filter_mac_kernel - Fused Moisan correction and transfer function MAC
 * @args: Struct containing parameter arrays and destination accumulators.
 *
 * Fuses the Fourier-domain boundary correction and the Fresnel transfer function application
 * (filters A and B for phase perturbation dphi and log-amplitude chi) into a single kernel pass,
 * eliminating intermediate spectrum DRAM round-trips.
 */
__global__ static void atmturb_cuda_moisan_filter_mac_kernel(
    atmturb_cuda_moisan_mac_t args)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= args.n_half || j >= args.pad_size)
    {
        return;
    }

    int idx = j * args.n_half + i;
    cufftComplex s = atmturb_cuda_moisan_calc(
        args.d_spec[idx], args.d_hat_bounds, args.d_exp_x, args.d_exp_y,
        args.d_laplace_inv, i, j, args.pad_size, args.n_half);

    float fa = args.d_filt_a[idx];
    float fb = args.d_filt_b[idx];

    args.d_acc_dphi[idx].x += s.x * fa;
    args.d_acc_dphi[idx].y += s.y * fa;
    args.d_acc_chi[idx].x  += s.x * fb;
    args.d_acc_chi[idx].y  += s.y * fb;
}

/**
 * atmturb_cuda_moisan_filter_mac_launch - Launch fused Moisan filter MAC on GPU
 * @args: Pointer to launcher arguments struct.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_moisan_filter_mac_launch(
    const atmturb_cuda_moisan_mac_t *args,
    cudaStream_t                     stream)
{
    dim3 b_spec(16, 16);
    dim3 g_spec(((int) args->n_half + b_spec.x - 1) / b_spec.x,
                ((int) args->pad_size + b_spec.y - 1) / b_spec.y);
    atmturb_cuda_moisan_filter_mac_kernel<<<g_spec, b_spec, 0, stream>>>(*args);
}

/**
 * atmturb_cuda_moisan_dual_mac_kernel - Fused Moisan correction and dual-wavelength MAC
 * @args: Struct containing parameter arrays and dual destination accumulators.
 *
 * Simultaneously evaluates primary and secondary wavelength Fresnel transfer functions
 * with chromatic dispersion phase tilt ramps directly from the boundary-corrected spectrum.
 */
__global__ static void atmturb_cuda_moisan_dual_mac_kernel(
    atmturb_cuda_moisan_dual_mac_t args)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= args.n_half || j >= args.pad_size)
    {
        return;
    }

    int idx = j * args.n_half + i;
    cufftComplex s = atmturb_cuda_moisan_calc(
        args.d_spec[idx], args.d_hat_bounds, args.d_exp_x, args.d_exp_y,
        args.d_laplace_inv, i, j, args.pad_size, args.n_half);

    float fa_pri = args.d_filt_a_pri[idx];
    float fb_pri = args.d_filt_b_pri[idx];

    args.d_acc_dphi_pri[idx].x += s.x * fa_pri;
    args.d_acc_dphi_pri[idx].y += s.y * fa_pri;
    args.d_acc_chi_pri[idx].x  += s.x * fb_pri;
    args.d_acc_chi_pri[idx].y  += s.y * fb_pri;

    cufftComplex r = args.d_ramp[idx];
    float ps_x = s.x * r.x - s.y * r.y;
    float ps_y = s.x * r.y + s.y * r.x;
    float fa_sec = args.d_filt_a_sec[idx];
    float fb_sec = args.d_filt_b_sec[idx];

    args.d_acc_dphi_sec[idx].x += ps_x * fa_sec;
    args.d_acc_dphi_sec[idx].y += ps_y * fa_sec;
    args.d_acc_chi_sec[idx].x  += ps_x * fb_sec;
    args.d_acc_chi_sec[idx].y  += ps_y * fb_sec;
}

/**
 * atmturb_cuda_moisan_dual_mac_launch - Launch fused Moisan dual MAC on GPU
 * @args: Pointer to dual launcher arguments struct.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_moisan_dual_mac_launch(
    const atmturb_cuda_moisan_dual_mac_t *args,
    cudaStream_t                          stream)
{
    dim3 b_spec(16, 16);
    dim3 g_spec(((int) args->n_half + b_spec.x - 1) / b_spec.x,
                ((int) args->pad_size + b_spec.y - 1) / b_spec.y);
    atmturb_cuda_moisan_dual_mac_kernel<<<g_spec, b_spec, 0, stream>>>(*args);
}

/**
 * atmturb_cuda_filter_mac_kernel - Vectorized single-wavelength filter multiply-accumulate
 * @d_acc_dphi: Phase perturbation spectrum accumulator.
 * @d_acc_chi: Log-amplitude spectrum accumulator.
 * @d_spec: Layer Fourier spectrum.
 * @d_filt_a: Fresnel filter A multiplier array.
 * @d_filt_b: Fresnel filter B multiplier array.
 * @ntot: Total number of half-complex frequency elements.
 */
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

/**
 * atmturb_cuda_filter_mac_launch - Launch single-wavelength filter MAC on GPU
 * @d_acc_dphi: Phase perturbation spectrum accumulator.
 * @d_acc_chi: Log-amplitude spectrum accumulator.
 * @d_spec: Layer Fourier spectrum.
 * @d_filt_a: Fresnel filter A multiplier array.
 * @d_filt_b: Fresnel filter B multiplier array.
 * @ntot: Total number of half-complex frequency elements.
 * @stream: CUDA execution stream.
 */
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

/**
 * atmturb_cuda_filter_mac_dual_kernel - Vectorized dual-wavelength filter multiply-accumulate
 * @d_acc_dphi_pri: Primary phase perturbation spectrum accumulator.
 * @d_acc_chi_pri: Primary log-amplitude spectrum accumulator.
 * @d_acc_dphi_sec: Secondary phase perturbation spectrum accumulator.
 * @d_acc_chi_sec: Secondary log-amplitude spectrum accumulator.
 * @d_spec: Layer Fourier spectrum.
 * @d_ramp: Chromatic dispersion Fourier phase ramp.
 * @d_filt_a_pri: Primary filter A array.
 * @d_filt_b_pri: Primary filter B array.
 * @d_filt_a_sec: Secondary filter A array.
 * @d_filt_b_sec: Secondary filter B array.
 * @ntot: Total number of half-complex frequency elements.
 */
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

/**
 * atmturb_cuda_filter_mac_dual_launch - Launch dual-wavelength filter MAC on GPU
 * @d_acc_dphi_pri: Primary phase perturbation spectrum accumulator.
 * @d_acc_chi_pri: Primary log-amplitude spectrum accumulator.
 * @d_acc_dphi_sec: Secondary phase perturbation spectrum accumulator.
 * @d_acc_chi_sec: Secondary log-amplitude spectrum accumulator.
 * @d_spec: Layer Fourier spectrum.
 * @d_ramp: Chromatic dispersion Fourier phase ramp.
 * @d_filt_a_pri: Primary filter A array.
 * @d_filt_b_pri: Primary filter B array.
 * @d_filt_a_sec: Secondary filter A array.
 * @d_filt_b_sec: Secondary filter B array.
 * @ntot: Total number of half-complex frequency elements.
 * @stream: CUDA execution stream.
 */
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

/**
 * atmturb_cuda_assemble_kernel - Assemble pupil phase and amplitude with parallel flux sum
 * @pup_size: Linear dimension of synthesized pupil.
 * @guard: Guard margin in pixels.
 * @pad_size: Linear dimension of padded simulation grid.
 * @d_frame_pha: Cumulative pupil phase array (modified in place).
 * @d_frame_amp: Destination pupil amplitude array.
 * @d_dphi_out: Synthesized diffractive phase perturbation array.
 * @d_chi_out: Synthesized diffractive log-amplitude perturbation array.
 * @d_sum_I: Pointer to device atomic intensity reduction accumulator.
 *
 * Adds Fresnel diffractive phase perturbation dphi to the geometric phase screen, clamps
 * log-amplitude chi to [-3, 3] to prevent numerical overflow, computes amplitude A = exp(chi),
 * and performs block-level warp reduction to accumulate sum(A^2) for flux normalization.
 */
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
        int max_clamp = (pup_size > 2) ? (pup_size - 2) : 0;
        int sx = (guard == 0 && (x == 0 || x == pup_size - 1))
                     ? ((x == 0) ? ((pup_size > 1) ? 1 : 0) : max_clamp)
                     : x;
        int sy = (guard == 0 && (y == 0 || y == pup_size - 1))
                     ? ((y == 0) ? ((pup_size > 1) ? 1 : 0) : max_clamp)
                     : y;
        int src_idx = (sy + guard) * pad_size + (sx + guard);
        int dst_idx = y * pup_size + x;

        d_frame_pha[dst_idx] += d_dphi_out[src_idx];
        float c = d_chi_out[src_idx];
        if (c > 3.0f)
        {
            c = 3.0f;
        }
        else if (c < -3.0f)
        {
            c = -3.0f;
        }
        float a = expf(c);
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

/**
 * atmturb_cuda_assemble_launch - Launch pupil phase and amplitude assembly on GPU
 * @pup_size: Linear dimension of synthesized pupil.
 * @guard: Guard margin in pixels.
 * @pad_size: Linear dimension of padded simulation grid.
 * @d_frame_pha: Cumulative pupil phase array.
 * @d_frame_amp: Destination pupil amplitude array.
 * @d_dphi_out: Synthesized diffractive phase perturbation array.
 * @d_chi_out: Synthesized diffractive log-amplitude perturbation array.
 * @d_sum_I: Pointer to device intensity accumulator.
 * @stream: CUDA execution stream.
 */
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

/**
 * atmturb_cuda_normalize_kernel - Normalize pupil amplitude array by root-mean-square intensity
 * @d_frame_amp: Pupil amplitude array (modified in place).
 * @npix: Total number of pupil pixels.
 * @d_sum_I: Pointer to device total intensity reduction result.
 */
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

/**
 * atmturb_cuda_normalize_launch - Launch pupil amplitude normalization on GPU
 * @d_frame_amp: Pupil amplitude array (modified in place).
 * @npix: Total number of pupil pixels.
 * @d_sum_I: Pointer to device total intensity reduction result.
 * @stream: CUDA execution stream.
 */
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
