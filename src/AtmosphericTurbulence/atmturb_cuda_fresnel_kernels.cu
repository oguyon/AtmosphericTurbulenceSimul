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
 * atmturb_cuda_fresnel_init_modulate_kernel - Initialize unit-amplitude field without sponge
 */
__global__ static void atmturb_cuda_fresnel_init_modulate_kernel(
    cufftComplex *d_field,
    const float  *d_super_pha,
    int           tot_pix)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= tot_pix)
    {
        return;
    }

    float p = d_super_pha[i];
    float s, c;
    __sincosf(p, &s, &c);
    d_field[i] = make_cuFloatComplex(c, s);
}

/**
 * atmturb_cuda_fresnel_init_sponge_kernel - Initialize complex field with 1D rim absorber
 */
__global__ static void atmturb_cuda_fresnel_init_sponge_kernel(
    cufftComplex *d_field,
    const float  *d_super_pha,
    const float  *d_sponge_1d,
    int           pad_size,
    int           count)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int b = (int) blockIdx.z;
    if (x >= pad_size || y >= pad_size || b >= count)
    {
        return;
    }

    size_t idx = ((size_t) b * (size_t) pad_size + (size_t) y) * (size_t) pad_size + (size_t) x;
    float  w   = __ldg(&d_sponge_1d[x]) * __ldg(&d_sponge_1d[y]);
    float  p   = d_super_pha[idx];
    float  s, c;
    __sincosf(p, &s, &c);
    d_field[idx] = make_cuFloatComplex(w * c, w * s);
}

/**
 * atmturb_cuda_fresnel_init_modulate_launch - Launch complex field initialization on GPU
 */
void atmturb_cuda_fresnel_init_modulate_launch(
    cufftComplex *d_field,
    const float  *d_super_pha,
    const float  *d_sponge_1d,
    int           pad_size,
    int           count,
    cudaStream_t  stream)
{
    if (pad_size <= 0 || count <= 0)
    {
        return;
    }
    if (d_sponge_1d == NULL)
    {
        int  tot  = count * pad_size * pad_size;
        dim3 b_1d(256);
        dim3 g_1d((unsigned int) ((tot + 255) / 256));
        atmturb_cuda_fresnel_init_modulate_kernel<<<g_1d, b_1d, 0, stream>>>(
            d_field, d_super_pha, tot);
        return;
    }

    dim3 b_3d(32, 8, 1);
    dim3 g_3d(((unsigned int) pad_size + 31u) / 32u,
              ((unsigned int) pad_size + 7u) / 8u,
              (unsigned int) count);
    atmturb_cuda_fresnel_init_sponge_kernel<<<g_3d, b_3d, 0, stream>>>(
        d_field, d_super_pha, d_sponge_1d, pad_size, count);
}

/**
 * atmturb_cuda_fresnel_modulate_kernel - Multiply complex field by thin phase screen
 */
__global__ static void atmturb_cuda_fresnel_modulate_kernel(
    cufftComplex *d_field,
    const float  *d_super_pha,
    int           tot_pix)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= tot_pix)
    {
        return;
    }

    float p = d_super_pha[i];
    float s, c;
    __sincosf(p, &s, &c);
    cufftComplex f = d_field[i];
    d_field[i]     = make_cuFloatComplex(f.x * c - f.y * s, f.y * c + f.x * s);
}

/**
 * atmturb_cuda_fresnel_modulate_sponge_kernel - Modulate complex field and apply rim absorber
 */
__global__ static void atmturb_cuda_fresnel_modulate_sponge_kernel(
    cufftComplex *d_field,
    const float  *d_super_pha,
    const float  *d_sponge_1d,
    int           pad_size,
    int           count)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int b = (int) blockIdx.z;
    if (x >= pad_size || y >= pad_size || b >= count)
    {
        return;
    }

    size_t idx = ((size_t) b * (size_t) pad_size + (size_t) y) * (size_t) pad_size + (size_t) x;
    float  w   = __ldg(&d_sponge_1d[x]) * __ldg(&d_sponge_1d[y]);
    float  p   = d_super_pha[idx];
    float  s, c;
    __sincosf(p, &s, &c);
    cufftComplex f = d_field[idx];
    d_field[idx]   = make_cuFloatComplex(w * (f.x * c - f.y * s), w * (f.y * c + f.x * s));
}

/**
 * atmturb_cuda_fresnel_modulate_launch - Launch thin phase screen modulation on GPU
 */
void atmturb_cuda_fresnel_modulate_launch(
    cufftComplex *d_field,
    const float  *d_super_pha,
    const float  *d_sponge_1d,
    int           pad_size,
    int           count,
    cudaStream_t  stream)
{
    if (pad_size <= 0 || count <= 0)
    {
        return;
    }
    if (d_sponge_1d == NULL)
    {
        int  tot  = count * pad_size * pad_size;
        dim3 b_1d(256);
        dim3 g_1d((unsigned int) ((tot + 255) / 256));
        atmturb_cuda_fresnel_modulate_kernel<<<g_1d, b_1d, 0, stream>>>(
            d_field, d_super_pha, tot);
        return;
    }

    dim3 b_3d(32, 8, 1);
    dim3 g_3d(((unsigned int) pad_size + 31u) / 32u,
              ((unsigned int) pad_size + 7u) / 8u,
              (unsigned int) count);
    atmturb_cuda_fresnel_modulate_sponge_kernel<<<g_3d, b_3d, 0, stream>>>(
        d_field, d_super_pha, d_sponge_1d, pad_size, count);
}

/**
 * atmturb_cuda_fresnel_tf1d_mul_kernel - Multiply 2D Fourier spectrum by separable 1D chirp
 */
__global__ static void atmturb_cuda_fresnel_tf1d_mul_kernel(
    cufftComplex       *d_field,
    const cufftComplex *d_tf1d,
    int                 pad_size,
    int                 count)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int b = (int) blockIdx.z;
    if (x >= pad_size || y >= pad_size || b >= count)
    {
        return;
    }

    cufftComplex hx  = __ldg(&d_tf1d[x]);
    cufftComplex hy  = __ldg(&d_tf1d[y]);
    float        tre = hy.x * hx.x - hy.y * hx.y;
    float        tim = hy.x * hx.y + hy.y * hx.x;
    size_t       idx = ((size_t) b * (size_t) pad_size + (size_t) y) * (size_t) pad_size +
                       (size_t) x;
    cufftComplex f   = d_field[idx];
    d_field[idx]     = make_cuFloatComplex(f.x * tre - f.y * tim, f.y * tre + f.x * tim);
}

/**
 * atmturb_cuda_fresnel_tf1d_mul_launch - Launch separable 1D Fresnel chirp multiplication
 */
void atmturb_cuda_fresnel_tf1d_mul_launch(
    cufftComplex       *d_field,
    const cufftComplex *d_tf1d,
    int                 pad_size,
    int                 count,
    cudaStream_t        stream)
{
    if (pad_size <= 0 || count <= 0)
    {
        return;
    }
    dim3 b_3d(32, 8, 1);
    dim3 g_3d(((unsigned int) pad_size + 31u) / 32u,
              ((unsigned int) pad_size + 7u) / 8u,
              (unsigned int) count);
    atmturb_cuda_fresnel_tf1d_mul_kernel<<<g_3d, b_3d, 0, stream>>>(
        d_field, d_tf1d, pad_size, count);
}

/**
 * atmturb_cuda_fresnel_extract_kernel - Extract cropped pupil amplitude, phase, and piston sum
 */
__global__ static void atmturb_cuda_fresnel_extract_kernel(
    atmturb_cuda_fresnel_extract_args_t a)
{
    int   pup_sq  = a.pup_size * a.pup_size;
    int   tot_pup = a.count * pup_sq;
    int   i       = blockIdx.x * blockDim.x + threadIdx.x;
    float p       = 0.0f;

    if (i < tot_pup)
    {
        size_t pad_idx = (size_t) i;
        if (a.guard > 0)
        {
            int b   = i / pup_sq;
            int rem = i - b * pup_sq;
            int py  = rem / a.pup_size;
            int px  = rem - py * a.pup_size;
            pad_idx = ((size_t) b * (size_t) a.pad_size + (size_t) (py + a.guard)) *
                      (size_t) a.pad_size + (size_t) (px + a.guard);
        }

        cufftComplex f  = a.d_field[pad_idx];
        float        re = f.x;
        float        im = f.y;
        if (a.d_amp != NULL)
        {
            a.d_amp[i] = hypotf(re, im);
        }
        float phi_geo = a.d_pha[i];
        float s_geo, c_geo;
        __sincosf(phi_geo, &s_geo, &c_geo);
        p          = phi_geo + atan2f(im * c_geo - re * s_geo, re * c_geo + im * s_geo);
        a.d_pha[i] = p;
    }

    if (a.d_sum != NULL)
    {
        for (int offset = 16; offset > 0; offset /= 2)
        {
            p += __shfl_down_sync(0xffffffff, p, offset);
        }
        __shared__ float s_warp_sums[32];
        int              lane    = threadIdx.x % 32;
        int              warp_id = threadIdx.x / 32;
        if (lane == 0)
        {
            s_warp_sums[warp_id] = p;
        }
        __syncthreads();

        if (threadIdx.x < 32)
        {
            int   nwarps = blockDim.x / 32;
            float bval   = (threadIdx.x < nwarps) ? s_warp_sums[threadIdx.x] : 0.0f;
            for (int offset = 16; offset > 0; offset /= 2)
            {
                bval += __shfl_down_sync(0xffffffff, bval, offset);
            }
            if (threadIdx.x == 0)
            {
                atomicAdd(a.d_sum, bval);
            }
        }
    } // if (a.d_sum != NULL)
}

/**
 * atmturb_cuda_fresnel_extract_launch - Launch pupil amplitude and unwrapped phase extraction
 */
void atmturb_cuda_fresnel_extract_launch(
    const atmturb_cuda_fresnel_extract_args_t *args,
    cudaStream_t                               stream)
{
    if (args == NULL || args->pup_size <= 0 || args->count <= 0)
    {
        return;
    }
    int  tot  = args->count * args->pup_size * args->pup_size;
    dim3 b_1d(256);
    dim3 g_1d((unsigned int) ((tot + 255) / 256));
    atmturb_cuda_fresnel_extract_kernel<<<g_1d, b_1d, 0, stream>>>(*args);
}
