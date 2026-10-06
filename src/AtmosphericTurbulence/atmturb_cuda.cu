// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda.cu
 * @brief   CUDA-accelerated multi-layer wavefront extrusion
 */

#include <cuda_runtime.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "atmturb_cuda.h"

/**
 * atmturb_render_kernel - GPU kernel synthesizing multi-layer extruded wavefront frames
 * @d_masters: Device master screen data concatenated.
 * @msize: Master screen dimension.
 * @pup_size: Synthesized pupil dimension.
 * @nbframes: Number of frames.
 * @nblayers: Number of turbulence layers.
 * @d_vx: X velocity per layer.
 * @d_vy: Y velocity per layer.
 * @d_weights: Layer Cn2 amplitude weights.
 * @Scoeff: Secondary wavelength scale factor.
 * @d_pha: Primary phase output array.
 * @d_amp: Primary amplitude output array.
 * @d_spha: Secondary phase output array.
 * @d_samp: Secondary amplitude output array.
 */
__global__ static void atmturb_render_kernel(const float *d_masters, long msize,
                                             long pup_size, long nbframes,
                                             long nblayers, const float *d_vx,
                                             const float *d_vy, const float *d_weights,
                                             float Scoeff, float *d_pha,
                                             float *d_amp, float *d_spha,
                                             float *d_samp)
{
    long ii = blockIdx.x * blockDim.x + threadIdx.x;
    long jj = blockIdx.y * blockDim.y + threadIdx.y;
    long t  = blockIdx.z;

    if (ii >= pup_size || jj >= pup_size || t >= nbframes)
    {
        return;
    }

    float total_pha = 0.0f;
    for (long k = 0; k < nblayers; k++)
    {
        float cur_x = (float)(t + 1) * d_vx[k] + (float)ii;
        float cur_y = (float)(t + 1) * d_vy[k] + (float)jj;

        float floor_x = floorf(cur_x);
        float floor_y = floorf(cur_y);
        float fx      = cur_x - floor_x;
        float fy      = cur_y - floor_y;

        long ix0 = ((long)floor_x) % msize;
        if (ix0 < 0) ix0 += msize;
        long iy0 = ((long)floor_y) % msize;
        if (iy0 < 0) iy0 += msize;

        long ix1 = (ix0 + 1) % msize;
        long iy1 = (iy0 + 1) % msize;

        const float *screen = d_masters + k * (msize * msize);
        float v00 = screen[iy0 * msize + ix0];
        float v10 = screen[iy0 * msize + ix1];
        float v01 = screen[iy1 * msize + ix0];
        float v11 = screen[iy1 * msize + ix1];

        float val = (1.0f - fx) * (1.0f - fy) * v00 +
                    fx * (1.0f - fy) * v10 +
                    (1.0f - fx) * fy * v01 +
                    fx * fy * v11;
        total_pha += d_weights[k] * val;
    }

    long out_idx     = t * (pup_size * pup_size) + jj * pup_size + ii;
    d_pha[out_idx]   = total_pha;
    d_amp[out_idx]   = 1.0f;
    d_spha[out_idx]  = total_pha * Scoeff;
    d_samp[out_idx]  = 1.0f;
}

/**
 * atmturb_cuda_device_available - Check if CUDA GPU is present and ready
 *
 * Return: 1 if CUDA device available, 0 otherwise.
 */
int atmturb_cuda_device_available(void)
{
    int count = 0;
    cudaError_t err = cudaGetDeviceCount(&count);
    return (err == cudaSuccess && count > 0) ? 1 : 0;
}

/**
 * atmturb_cuda_alloc_upload_inputs - Upload master screens and kinematics to GPU
 * @p: Simulation configuration.
 * @d_m: Pointer to allocated device master screens.
 * @d_vx: Pointer to allocated device X velocity array.
 * @d_vy: Pointer to allocated device Y velocity array.
 * @d_w: Pointer to allocated device weight array.
 *
 * Return: 0 on success, -1 on allocation/copy failure.
 */
static int atmturb_cuda_alloc_upload_inputs(const atmturb_cuda_sim_params_t *p,
                                           float **d_m, float **d_vx,
                                           float **d_vy, float **d_w)
{
    size_t screen_bytes = sizeof(float) * p->msize * p->msize;
    size_t total_mbytes = screen_bytes * p->nblayers;
    if (cudaMalloc((void **)d_m, total_mbytes) != cudaSuccess) return -1;

    for (long k = 0; k < p->nblayers; k++)
    {
        cudaMemcpy(*d_m + k * p->msize * p->msize, p->h_masters[k],
                   screen_bytes, cudaMemcpyHostToDevice);
    }

    float *h_vx = (float *)malloc(sizeof(float) * p->nblayers);
    float *h_vy = (float *)malloc(sizeof(float) * p->nblayers);
    float *h_w  = (float *)malloc(sizeof(float) * p->nblayers);
    for (long k = 0; k < p->nblayers; k++)
    {
        h_vx[k] = (float)p->vxpix[k];
        h_vy[k] = (float)p->vypix[k];
        h_w[k]  = (float)sqrt(p->cn2[k]);
    }

    size_t vec_bytes = sizeof(float) * p->nblayers;
    cudaMalloc((void **)d_vx, vec_bytes);
    cudaMalloc((void **)d_vy, vec_bytes);
    cudaMalloc((void **)d_w, vec_bytes);
    cudaMemcpy(*d_vx, h_vx, vec_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(*d_vy, h_vy, vec_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(*d_w, h_w, vec_bytes, cudaMemcpyHostToDevice);

    free(h_vx);
    free(h_vy);
    free(h_w);
    return 0;
}

/**
 * atmturb_wfs_render_frames_cuda - CUDA GPU accelerated wavefront time-series rendering
 * @params: Simulation parameters and input screen pointers.
 * @outputs: Destination host arrays for synthesized frames.
 *
 * Return: 0 on success, -1 on CUDA runtime error.
 */
int atmturb_wfs_render_frames_cuda(const atmturb_cuda_sim_params_t *params,
                                   atmturb_cuda_sim_outputs_t *outputs)
{
    float *d_m = NULL, *d_vx = NULL, *d_vy = NULL, *d_w = NULL;
    if (atmturb_cuda_alloc_upload_inputs(params, &d_m, &d_vx, &d_vy, &d_w) != 0)
    {
        return -1;
    }

    size_t out_bytes = sizeof(float) * params->nbframes * params->pup_size * params->pup_size;
    float *d_pha = NULL, *d_amp = NULL, *d_spha = NULL, *d_samp = NULL;
    cudaMalloc((void **)&d_pha, out_bytes);
    cudaMalloc((void **)&d_amp, out_bytes);
    cudaMalloc((void **)&d_spha, out_bytes);
    cudaMalloc((void **)&d_samp, out_bytes);

    dim3 block(16, 16);
    dim3 grid((params->pup_size + block.x - 1) / block.x,
              (params->pup_size + block.y - 1) / block.y,
              params->nbframes);

    atmturb_render_kernel<<<grid, block>>>(d_m, params->msize, params->pup_size,
                                           params->nbframes, params->nblayers,
                                           d_vx, d_vy, d_w, (float)params->Scoeff,
                                           d_pha, d_amp, d_spha, d_samp);

    cudaMemcpy(outputs->pha, d_pha, out_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(outputs->amp, d_amp, out_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(outputs->spha, d_spha, out_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(outputs->samp, d_samp, out_bytes, cudaMemcpyDeviceToHost);

    cudaFree(d_m);
    cudaFree(d_vx);
    cudaFree(d_vy);
    cudaFree(d_w);
    cudaFree(d_pha);
    cudaFree(d_amp);
    cudaFree(d_spha);
    cudaFree(d_samp);
    return 0;
}
