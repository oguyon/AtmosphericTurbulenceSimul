// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_piston.cu
 * @brief   GPU-accelerated aperture piston reduction and subtraction kernels
 *
 * Implements high-throughput parallel reduction using warp-shuffle instructions
 * (__shfl_down_sync) across 32 threads, followed by inter-warp accumulation in
 * shared memory and a single atomicAdd to global memory. This avoids host roundtrips
 * and CPU synchronization stalls during simulation runs.
 */

#include "atmturb_cuda_piston.h"

#include <cuda_runtime.h>

/**
 * atmturb_cuda_piston_sum_kernel - Parallel reduction kernel computing total aperture phase
 * @d_pha: Input device phase array.
 * @d_sum: Output device float accumulator receiving total sum.
 * @npix: Total number of pixels in pupil aperture.
 *
 * Each thread accumulates a grid-stride strided sum of elements, followed by a warp-level
 * shuffle reduction down to lane 0. Lane 0 of each warp writes to shared memory, which the
 * first warp then reduces across all warps before atomically adding to @d_sum.
 */
__global__ static void atmturb_cuda_piston_sum_kernel(
    const float *d_pha,
    float       *d_sum,
    int          npix)
{
    float val = 0.0f;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = idx; i < npix; i += stride)
    {
        val += d_pha[i];
    }

    /* Warp-level parallel reduction */
    for (int offset = 16; offset > 0; offset /= 2)
    {
        val += __shfl_down_sync(0xffffffff, val, offset);
    }

    __shared__ float s_warp_sums[32];
    int lane = threadIdx.x % 32;
    int warp_id = threadIdx.x / 32;
    if (lane == 0)
    {
        s_warp_sums[warp_id] = val;
    }
    __syncthreads();

    /* Inter-warp reduction across warp leaders */
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

/**
 * atmturb_cuda_piston_sub_kernel - Subtract aperture mean piston from each phase pixel
 * @d_pha: In-place device phase array to zero-center.
 * @d_sum: Pointer to device float holding precomputed aperture sum.
 * @npix: Total number of pixels in pupil aperture.
 */
__global__ static void atmturb_cuda_piston_sub_kernel(
    float       *d_pha,
    const float *d_sum,
    int          npix)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= npix)
    {
        return;
    }
    float mean = (*d_sum) / (float) npix;
    d_pha[idx] -= mean;
}

/**
 * atmturb_cuda_remove_piston_launch - Calculate and subtract aperture-mean piston on GPU
 * @d_pha: In-place phase screen array on device.
 * @d_sum_buf: Single-float working buffer on device for reduction sum.
 * @npix: Total number of active pupil pixels in phase screen.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_remove_piston_launch(
    float       *d_pha,
    float       *d_sum_buf,
    int          npix,
    cudaStream_t stream)
{
    if (d_pha == NULL || d_sum_buf == NULL || npix <= 0)
    {
        return;
    }

    cudaMemsetAsync(d_sum_buf, 0, sizeof(float), stream);

    dim3 block(256);
    int num_blocks = (npix + 255) / 256;
    if (num_blocks > 64)
    {
        num_blocks = 64;
    }
    dim3 grid(num_blocks);

    atmturb_cuda_piston_sum_kernel<<<grid, block, 0, stream>>>(d_pha, d_sum_buf, npix);

    dim3 full_grid((npix + 255) / 256);
    atmturb_cuda_piston_sub_kernel<<<full_grid, block, 0, stream>>>(d_pha, d_sum_buf, npix);
}

/**
 * atmturb_cuda_piston_sub_launch - Subtract precomputed mean piston from phase screen on GPU
 * @d_pha: In-place phase screen array on device.
 * @d_sum_buf: Pointer to precomputed reduction sum of phase values on device.
 * @npix: Total number of active pupil pixels in phase screen.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_piston_sub_launch(
    float       *d_pha,
    const float *d_sum_buf,
    int          npix,
    cudaStream_t stream)
{
    if (d_pha == NULL || d_sum_buf == NULL || npix <= 0)
    {
        return;
    }

    dim3 block(256);
    dim3 grid((npix + 255) / 256);
    atmturb_cuda_piston_sub_kernel<<<grid, block, 0, stream>>>(d_pha, d_sum_buf, npix);
}
