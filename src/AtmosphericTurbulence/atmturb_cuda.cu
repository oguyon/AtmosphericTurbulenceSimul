// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda.cu
 * @brief   CUDA-accelerated multi-layer wavefront extrusion with persistent caching
 */

#include <cuda_runtime.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "atmturb_cuda.h"

// Maximum supported layers in constant memory cache
#define ATMTURB_CUDA_MAX_LAYERS 64
__constant__ float c_vx[ATMTURB_CUDA_MAX_LAYERS];
__constant__ float c_vy[ATMTURB_CUDA_MAX_LAYERS];
__constant__ float c_x0[ATMTURB_CUDA_MAX_LAYERS];
__constant__ float c_y0[ATMTURB_CUDA_MAX_LAYERS];
__constant__ float c_weights[ATMTURB_CUDA_MAX_LAYERS];

// Persistent GPU memory caches to avoid per-render cudaMalloc/free overhead
static float              *s_d_m = NULL;
static size_t              s_cached_mbytes = 0;
static const float *const *s_cached_h_masters = NULL;

static float *s_d_pha  = NULL;
static float *s_d_amp  = NULL;
static float *s_d_spha = NULL;
static float *s_d_samp = NULL;
static size_t s_cached_out_bytes = 0;

/**
 * atmturb_render_kernel - GPU kernel synthesizing multi-layer extruded wavefront frames
 * @d_masters: Device master screen data concatenated.
 * @msize: Master screen dimension.
 * @pup_size: Synthesized pupil dimension.
 * @nbframes: Number of frames.
 * @nblayers: Number of turbulence layers.
 * @Scoeff: Secondary wavelength scale factor.
 * @d_pha: Primary phase output array.
 * @d_amp: Primary amplitude output array.
 * @d_spha: Secondary phase output array.
 * @d_samp: Secondary amplitude output array.
 */
__global__ static void atmturb_render_kernel(
    const float *d_masters,
    int          msize,
    int          pup_size,
    int          nbframes,
    int          nblayers,
    float        Scoeff,
    float       *d_pha,
    float       *d_amp,
    float       *d_spha,
    float       *d_samp)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;
    int jj = blockIdx.y * blockDim.y + threadIdx.y;
    int t  = blockIdx.z;

    if (ii >= pup_size || jj >= pup_size || t >= nbframes)
    {
        return;
    }

    int is_pow2 = ((msize & (msize - 1)) == 0);
    int mask = msize - 1;
    int msize_sq = msize * msize;
    float total_pha = 0.0f;

    for (int k = 0; k < nblayers; k++)
    {
        float cur_x = c_x0[k] + (float) t * c_vx[k] + (float) ii;
        float cur_y = c_y0[k] + (float) t * c_vy[k] + (float) jj;

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

        const float *screen = d_masters + k * msize_sq;
        float v00 = __ldg(&screen[iy0 * msize + ix0]);
        float v10 = __ldg(&screen[iy0 * msize + ix1]);
        float v01 = __ldg(&screen[iy1 * msize + ix0]);
        float v11 = __ldg(&screen[iy1 * msize + ix1]);

        float top = v00 + fx * (v10 - v00);
        float bot = v01 + fx * (v11 - v01);
        float val = top + fy * (bot - top);
        total_pha += c_weights[k] * val;
    }

    int out_idx     = t * (pup_size * pup_size) + jj * pup_size + ii;
    d_pha[out_idx]  = total_pha;
    if (d_amp != NULL)
    {
        d_amp[out_idx]  = 1.0f;
    }
    d_spha[out_idx] = total_pha * Scoeff;
    if (d_samp != NULL)
    {
        d_samp[out_idx] = 1.0f;
    }
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
 * atmturb_cuda_cleanup - Release persistent GPU buffers and context
 */
void atmturb_cuda_cleanup(void)
{
    if (s_d_m != NULL)    { cudaFree(s_d_m);    s_d_m = NULL; }
    if (s_d_pha != NULL)  { cudaFree(s_d_pha);  s_d_pha = NULL; }
    if (s_d_amp != NULL)  { cudaFree(s_d_amp);  s_d_amp = NULL; }
    if (s_d_spha != NULL) { cudaFree(s_d_spha); s_d_spha = NULL; }
    if (s_d_samp != NULL) { cudaFree(s_d_samp); s_d_samp = NULL; }
    s_cached_mbytes = 0;
    s_cached_h_masters = NULL;
    s_cached_out_bytes = 0;
}

/**
 * atmturb_cuda_sync_masters - Upload master screens if changed or not yet cached
 * @p: Simulation parameters.
 *
 * Return: 0 on success, -1 on allocation or copy failure.
 */
static int atmturb_cuda_sync_masters(
    const atmturb_cuda_sim_params_t *p)
{
    size_t screen_bytes = sizeof(float) * (size_t) (p->msize * p->msize);
    size_t total_mbytes = screen_bytes * (size_t) p->nblayers;

    if (s_d_m == NULL || s_cached_mbytes != total_mbytes ||
        s_cached_h_masters != p->h_masters)
    {
        if (s_d_m != NULL)
        {
            cudaFree(s_d_m);
        }
        if (cudaMalloc((void **) &s_d_m, total_mbytes) != cudaSuccess)
        {
            return -1;
        }
        for (long k = 0; k < p->nblayers; k++)
        {
            cudaMemcpy(s_d_m + k * p->msize * p->msize, p->h_masters[k],
                       screen_bytes, cudaMemcpyHostToDevice);
        }
        s_cached_mbytes = total_mbytes;
        s_cached_h_masters = p->h_masters;
    }
    return 0;
}

/**
 * atmturb_cuda_sync_kinematics - Upload layer velocity and position offsets to constant memory
 * @p: Simulation parameters.
 *
 * Return: 0 on success, -1 if layers exceed constant buffer capacity.
 */
static int atmturb_cuda_sync_kinematics(
    const atmturb_cuda_sim_params_t *p)
{
    if (p->nblayers > ATMTURB_CUDA_MAX_LAYERS)
    {
        return -1;
    }

    float h_vx[ATMTURB_CUDA_MAX_LAYERS];
    float h_vy[ATMTURB_CUDA_MAX_LAYERS];
    float h_w[ATMTURB_CUDA_MAX_LAYERS];
    float h_x0[ATMTURB_CUDA_MAX_LAYERS];
    float h_y0[ATMTURB_CUDA_MAX_LAYERS];

    for (long k = 0; k < p->nblayers; k++)
    {
        h_vx[k] = (float) p->vxpix[k];
        h_vy[k] = (float) p->vypix[k];
        h_w[k]  = (float) sqrt(p->cn2[k]);
        h_x0[k] = (p->x0 != NULL) ? (float) p->x0[k] : (float) p->vxpix[k];
        h_y0[k] = (p->y0 != NULL) ? (float) p->y0[k] : (float) p->vypix[k];
    }

    size_t vec_bytes = sizeof(float) * (size_t) p->nblayers;
    cudaMemcpyToSymbol(c_vx, h_vx, vec_bytes);
    cudaMemcpyToSymbol(c_vy, h_vy, vec_bytes);
    cudaMemcpyToSymbol(c_weights, h_w, vec_bytes);
    cudaMemcpyToSymbol(c_x0, h_x0, vec_bytes);
    cudaMemcpyToSymbol(c_y0, h_y0, vec_bytes);
    return 0;
}

/**
 * atmturb_cuda_sync_outputs - Allocate device output buffers if capacity insufficient
 * @out_bytes: Required output buffer size in bytes.
 * @need_amp: Flag indicating if amplitude buffers are required on device.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
static int atmturb_cuda_sync_outputs(
    size_t out_bytes,
    int    need_amp)
{
    if (s_d_pha == NULL || s_cached_out_bytes < out_bytes)
    {
        if (s_d_pha != NULL)
        {
            cudaFree(s_d_pha);
            cudaFree(s_d_amp);
            cudaFree(s_d_spha);
            cudaFree(s_d_samp);
            s_d_amp = NULL;
            s_d_samp = NULL;
        }
        if (cudaMalloc((void **) &s_d_pha, out_bytes) != cudaSuccess ||
            cudaMalloc((void **) &s_d_spha, out_bytes) != cudaSuccess)
        {
            return -1;
        }
        if (need_amp)
        {
            if (cudaMalloc((void **) &s_d_amp, out_bytes) != cudaSuccess ||
                cudaMalloc((void **) &s_d_samp, out_bytes) != cudaSuccess)
            {
                return -1;
            }
        }
        s_cached_out_bytes = out_bytes;
    }
    return 0;
}

/**
 * atmturb_wfs_render_frames_cuda - CUDA GPU accelerated wavefront time-series rendering
 * @params: Simulation parameters and input screen pointers.
 * @outputs: Destination host arrays for synthesized frames.
 *
 * Return: 0 on success, -1 on CUDA runtime error.
 */
int atmturb_wfs_render_frames_cuda(
    const atmturb_cuda_sim_params_t *params,
    atmturb_cuda_sim_outputs_t      *outputs)
{
    if (atmturb_cuda_sync_masters(params) != 0 ||
        atmturb_cuda_sync_kinematics(params) != 0)
    {
        return -1;
    }

    size_t npix = (size_t) (params->nbframes * params->pup_size * params->pup_size);
    size_t out_bytes = sizeof(float) * npix;
    int need_amp = (outputs->amp != NULL || outputs->samp != NULL) ? 1 : 0;
    if (atmturb_cuda_sync_outputs(out_bytes, need_amp) != 0)
    {
        return -1;
    }

    dim3 block(32, 8);
    dim3 grid(((int) params->pup_size + block.x - 1) / block.x,
              ((int) params->pup_size + block.y - 1) / block.y,
              (int) params->nbframes);

    atmturb_render_kernel<<<grid, block>>>(
        s_d_m, (int) params->msize, (int) params->pup_size,
        (int) params->nbframes, (int) params->nblayers,
        (float) params->Scoeff,
        s_d_pha, s_d_amp, s_d_spha, s_d_samp);

    cudaMemcpy(outputs->pha, s_d_pha, out_bytes, cudaMemcpyDeviceToHost);
    if (outputs->amp != NULL && s_d_amp != NULL)
    {
        cudaMemcpy(outputs->amp, s_d_amp, out_bytes, cudaMemcpyDeviceToHost);
    }
    if (outputs->spha != NULL)
    {
        cudaMemcpy(outputs->spha, s_d_spha, out_bytes, cudaMemcpyDeviceToHost);
    }
    if (outputs->samp != NULL && s_d_samp != NULL)
    {
        cudaMemcpy(outputs->samp, s_d_samp, out_bytes, cudaMemcpyDeviceToHost);
    }

    return 0;
}
