// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    wfprop_fresnel_cuda.cu
 * @brief   CUDA cuFFT diffractive propagation implementation
 */

#include <cuda_runtime.h>
#include <cufft.h>
#include <math.h>
#include <stdio.h>

#include "wfprop_fresnel_cuda.h"

#ifndef PI
#    define PI 3.14159265358979323846264338328
#endif

/**
 * fresnel_c2c_kernel - Apply centered quadratic phase transfer in Fourier domain
 * @d_data: Device complex float array.
 * @nx: Grid width.
 * @ny: Grid height.
 * @coeff: Quadratic phase coefficient.
 * @scale: Normalization scale factor (1.0 / (nx * ny)).
 */
__global__ static void fresnel_c2c_kernel(cufftComplex *d_data, long nx, long ny,
                                          float coeff, float scale)
{
    long ii = blockIdx.x * blockDim.x + threadIdx.x;
    long jj = blockIdx.y * blockDim.y + threadIdx.y;
    if (ii >= nx || jj >= ny)
    {
        return;
    }

    long u = (ii < nx / 2) ? ii : (ii - nx);
    long v = (jj < ny / 2) ? jj : (jj - ny);
    float sqdist = (float)(u * u + v * v);
    float angle  = -coeff * sqdist;

    float s, c;
    __sincosf(angle, &s, &c);

    long idx = jj * nx + ii;
    float re = d_data[idx].x * scale;
    float im = d_data[idx].y * scale;

    d_data[idx].x = re * c - im * s;
    d_data[idx].y = re * s + im * c;
}

/**
 * fresnel_z2z_kernel - Apply centered quadratic phase transfer (double precision)
 * @d_data: Device complex double array.
 * @nx: Grid width.
 * @ny: Grid height.
 * @coeff: Quadratic phase coefficient.
 * @scale: Normalization scale factor (1.0 / (nx * ny)).
 */
__global__ static void fresnel_z2z_kernel(cufftDoubleComplex *d_data, long nx, long ny,
                                          double coeff, double scale)
{
    long ii = blockIdx.x * blockDim.x + threadIdx.x;
    long jj = blockIdx.y * blockDim.y + threadIdx.y;
    if (ii >= nx || jj >= ny)
    {
        return;
    }

    long u = (ii < nx / 2) ? ii : (ii - nx);
    long v = (jj < ny / 2) ? jj : (jj - ny);
    double sqdist = (double)(u * u + v * v);
    double angle  = -coeff * sqdist;

    double s, c;
    sincos(angle, &s, &c);

    long idx = jj * nx + ii;
    double re = d_data[idx].x * scale;
    double im = d_data[idx].y * scale;

    d_data[idx].x = re * c - im * s;
    d_data[idx].y = re * s + im * c;
}

/**
 * wfprop_fresnel_device_available - Test CUDA GPU availability
 *
 * Return: 1 if available, 0 otherwise.
 */
int wfprop_fresnel_device_available(void)
{
    int count = 0;
    cudaError_t err = cudaGetDeviceCount(&count);
    return (err == cudaSuccess && count > 0) ? 1 : 0;
}

/**
 * wfprop_fresnel_propagate_float - Single-precision cuFFT Fresnel propagation
 * @h_in: Input host array.
 * @h_out: Output host array.
 * @nx: Horizontal size.
 * @ny: Vertical size.
 * @coeff: Quadratic phase coefficient.
 *
 * Return: 0 on success, -1 on CUDA error.
 */
static int wfprop_fresnel_propagate_float(const void *h_in, void *h_out,
                                          long nx, long ny, float coeff)
{
    size_t nbytes = sizeof(cufftComplex) * nx * ny;
    cufftComplex *d_buf = NULL;
    if (cudaMalloc((void **)&d_buf, nbytes) != cudaSuccess)
    {
        return -1;
    }

    cudaMemcpy(d_buf, h_in, nbytes, cudaMemcpyHostToDevice);

    cufftHandle plan;
    cufftPlan2d(&plan, (int)ny, (int)nx, CUFFT_C2C);
    cufftExecC2C(plan, d_buf, d_buf, CUFFT_FORWARD);

    dim3 block(16, 16);
    dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y);
    float scale = 1.0f / (float)(nx * ny);
    fresnel_c2c_kernel<<<grid, block>>>(d_buf, nx, ny, coeff, scale);

    cufftExecC2C(plan, d_buf, d_buf, CUFFT_INVERSE);
    cudaMemcpy(h_out, d_buf, nbytes, cudaMemcpyDeviceToHost);

    cufftDestroy(plan);
    cudaFree(d_buf);
    return 0;
}

/**
 * wfprop_fresnel_propagate_double - Double-precision cuFFT Fresnel propagation
 * @h_in: Input host array.
 * @h_out: Output host array.
 * @nx: Horizontal size.
 * @ny: Vertical size.
 * @coeff: Quadratic phase coefficient.
 *
 * Return: 0 on success, -1 on CUDA error.
 */
static int wfprop_fresnel_propagate_double(const void *h_in, void *h_out,
                                           long nx, long ny, double coeff)
{
    size_t nbytes = sizeof(cufftDoubleComplex) * nx * ny;
    cufftDoubleComplex *d_buf = NULL;
    if (cudaMalloc((void **)&d_buf, nbytes) != cudaSuccess)
    {
        return -1;
    }

    cudaMemcpy(d_buf, h_in, nbytes, cudaMemcpyHostToDevice);

    cufftHandle plan;
    cufftPlan2d(&plan, (int)ny, (int)nx, CUFFT_Z2Z);
    cufftExecZ2Z(plan, d_buf, d_buf, CUFFT_FORWARD);

    dim3 block(16, 16);
    dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y);
    double scale = 1.0 / (double)(nx * ny);
    fresnel_z2z_kernel<<<grid, block>>>(d_buf, nx, ny, coeff, scale);

    cufftExecZ2Z(plan, d_buf, d_buf, CUFFT_INVERSE);
    cudaMemcpy(h_out, d_buf, nbytes, cudaMemcpyDeviceToHost);

    cufftDestroy(plan);
    cudaFree(d_buf);
    return 0;
}

/**
 * wfprop_fresnel_propagate_cuda - Top-level CUDA Fresnel propagation dispatcher
 * @h_in: Input host array.
 * @h_out: Output host array.
 * @nx: Horizontal grid dimension.
 * @ny: Vertical grid dimension.
 * @pupil_scale: Physical pixel sampling scale [m/pixel].
 * @z: Propagation distance [m].
 * @lambda: Optical wavelength [m].
 * @is_double: 1 for double precision, 0 for single precision.
 *
 * Return: 0 on success, -1 on CUDA error.
 */
int wfprop_fresnel_propagate_cuda(const void *h_in, void *h_out, long nx, long ny,
                                  double pupil_scale, double z, double lambda,
                                  int is_double)
{
    double coeff = PI * z * lambda / (pupil_scale * nx) / (pupil_scale * nx);
    if (is_double)
    {
        return wfprop_fresnel_propagate_double(h_in, h_out, nx, ny, coeff);
    }
    return wfprop_fresnel_propagate_float(h_in, h_out, nx, ny, (float)coeff);
}
