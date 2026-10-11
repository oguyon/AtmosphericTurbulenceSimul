// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_fresnel_kernels.cu
 * @brief   GPU device kernels for multi-layer Fresnel optical diffraction propagation
 */

#include <cuda_runtime.h>
#include <cufft.h>
#include <math.h>

#include "atmturb_cuda_fresnel_kernels.h"

/**
 * atmturb_cuda_fresnel_init_modulate_kernel - Initialize optical field from first phase screen
 * @d_field: Destination complex electric field array.
 * @d_super_pha: Source phase screen array on padded grid.
 * @npix: Total number of grid pixels.
 *
 * Converts the first atmospheric superlayer phase screen phi(x, y) into a complex optical
 * field of unit amplitude: E(x, y) = cos(phi) + i * sin(phi).
 */
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

/**
 * atmturb_cuda_fresnel_init_modulate_launch - Launch complex field initialization on GPU
 * @d_field: Destination complex electric field array.
 * @d_super_pha: Source phase screen array on padded grid.
 * @npix: Total number of grid pixels.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_fresnel_init_modulate_launch(
    cufftComplex *d_field,
    const float  *d_super_pha,
    int           npix,
    cudaStream_t  stream)
{
    dim3 b_1d(256);
    dim3 g_1d(((int) npix + b_1d.x - 1) / b_1d.x);
    atmturb_cuda_fresnel_init_modulate_kernel<<<g_1d, b_1d, 0, stream>>>(
        d_field, d_super_pha, npix);
}

/**
 * atmturb_cuda_fresnel_modulate_kernel - Multiply propagating complex field by thin phase screen
 * @d_field: Complex electric field array (modified in place).
 * @d_super_pha: Superlayer phase screen array.
 * @npix: Total number of grid pixels.
 *
 * Implements the thin phase screen approximation for an intermediate atmospheric layer:
 * E_out(x, y) = E_in(x, y) * exp(i * phi(x, y))
 *             = (re * cos(phi) - im * sin(phi)) + i * (im * cos(phi) + re * sin(phi)).
 */
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

/**
 * atmturb_cuda_fresnel_modulate_launch - Launch thin phase screen modulation on GPU
 * @d_field: Complex electric field array (modified in place).
 * @d_super_pha: Superlayer phase screen array.
 * @npix: Total number of grid pixels.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_fresnel_modulate_launch(
    cufftComplex *d_field,
    const float  *d_super_pha,
    int           npix,
    cudaStream_t  stream)
{
    dim3 b_1d(256);
    dim3 g_1d(((int) npix + b_1d.x - 1) / b_1d.x);
    atmturb_cuda_fresnel_modulate_kernel<<<g_1d, b_1d, 0, stream>>>(
        d_field, d_super_pha, npix);
}

/**
 * atmturb_cuda_fresnel_tf_mul_vec2_kernel - Vectorized optical transfer function multiply
 * @d_field: Spatial frequency complex field array (cast to float4).
 * @d_tf: Optical transfer function complex array (cast to float4).
 * @npairs: Number of complex pairs (npix / 2).
 *
 * Multiplies two independent complex numbers per thread using 128-bit vector float4 memory
 * transactions to maximize GPU memory bus saturation during diffractive Fresnel propagation:
 * E(u, v) <- E(u, v) * exp(-i * pi * lambda * z * (u^2 + v^2)).
 */
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

__global__ static void atmturb_cuda_fresnel_tf_mul_tail_kernel(
    cufftComplex       *d_field,
    const cufftComplex *d_tf,
    int                 idx)
{
    cufftComplex f = d_field[idx];
    cufftComplex t = d_tf[idx];
    cufftComplex out;
    out.x = f.x * t.x - f.y * t.y;
    out.y = f.y * t.x + f.x * t.y;
    d_field[idx] = out;
}

/**
 * atmturb_cuda_fresnel_tf_mul_vec2_launch - Launch vectorized optical transfer function multiply
 * @d_field: Complex spatial frequency field array.
 * @d_tf: Complex optical transfer function array.
 * @npix: Total number of complex elements.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_fresnel_tf_mul_vec2_launch(
    cufftComplex       *d_field,
    const cufftComplex *d_tf,
    int                 npix,
    cudaStream_t        stream)
{
    int npairs = npix / 2;
    if (npairs > 0)
    {
        dim3 b_1d(256);
        dim3 g_1d(((int) npairs + b_1d.x - 1) / b_1d.x);
        atmturb_cuda_fresnel_tf_mul_vec2_kernel<<<g_1d, b_1d, 0, stream>>>(
            (float4 *) d_field, (const float4 *) d_tf, npairs);
    }
    if ((npix & 1) != 0)
    {
        atmturb_cuda_fresnel_tf_mul_tail_kernel<<<1, 1, 0, stream>>>(
            d_field, d_tf, npix - 1);
    }
}

/**
 * atmturb_cuda_fresnel_extract_kernel - Extract amplitude, unwrapped phase, and parallel piston sum
 * @d_field: Propagated complex electric field array.
 * @d_pha: Input geometric reference phase and output total unwrapped phase.
 * @d_amp: Destination amplitude array (or NULL).
 * @d_sum: Optional device memory pointer for atomic piston reduction (or NULL).
 * @npix: Total number of grid pixels.
 *
 * Extracts amplitude |E| = hypot(re, im), demodulates carrier phase using geometric reference
 * phase phi_geo to eliminate 2*pi wrapping ambiguity:
 *   delta_phi = atan2(im * cos(phi_geo) - re * sin(phi_geo), re * cos(phi_geo) + im * sin(phi_geo))
 *   phi_out = phi_geo + delta_phi
 * and performs block-level warp reduction to accumulate the phase sum for zero-mean piston removal.
 */
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

/**
 * atmturb_cuda_fresnel_extract_launch - Launch amplitude and phase extraction on GPU
 * @d_field: Input propagated complex electric field array.
 * @d_pha: Input geometric reference and output total unwrapped phase.
 * @d_amp: Destination amplitude array (or NULL).
 * @d_sum: Optional device memory pointer for atomic piston reduction (or NULL).
 * @npix: Total number of grid pixels.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_fresnel_extract_launch(
    const cufftComplex *d_field,
    float              *d_pha,
    float              *d_amp,
    float              *d_sum,
    int                 npix,
    cudaStream_t        stream)
{
    dim3 b_1d(256);
    dim3 g_1d(((int) npix + b_1d.x - 1) / b_1d.x);
    atmturb_cuda_fresnel_extract_kernel<<<g_1d, b_1d, 0, stream>>>(
        d_field, d_pha, d_amp, d_sum, npix);
}
