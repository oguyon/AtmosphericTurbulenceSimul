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
 * @t_offset: Starting frame offset for streaming synthesis.
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
    int          t_offset,
    float       *d_pha,
    float       *d_amp,
    float       *d_spha,
    float       *d_samp)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;
    int jj = blockIdx.y * blockDim.y + threadIdx.y;
    int f  = blockIdx.z;

    if (ii >= pup_size || jj >= pup_size || f >= nbframes)
    {
        return;
    }

    int t = t_offset + f;

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

    int out_idx     = f * (pup_size * pup_size) + jj * pup_size + ii;
    d_pha[out_idx]  = total_pha;
    if (d_amp != NULL)
    {
        d_amp[out_idx]  = 1.0f;
    }
    if (d_spha != NULL)
    {
        d_spha[out_idx] = total_pha * Scoeff;
    }
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
 * atmturb_cuda_sync_device_masters - Upload and cache master screens in GPU device memory
 * @msize: Master screen dimension in pixels.
 * @nblayers: Number of simulation layers.
 * @h_masters: Array of host screen pointers.
 *
 * Return: Pointer to device master screens array, or NULL on failure.
 */
float *atmturb_cuda_sync_device_masters(
    long                msize,
    long                nblayers,
    const float *const *h_masters)
{
    size_t screen_bytes = sizeof(float) * (size_t) (msize * msize);
    size_t total_mbytes = screen_bytes * (size_t) nblayers;

    if (s_d_m == NULL || s_cached_mbytes != total_mbytes ||
        s_cached_h_masters != h_masters)
    {
        if (s_d_m != NULL)
        {
            cudaFree(s_d_m);
        }
        if (cudaMalloc((void **) &s_d_m, total_mbytes) != cudaSuccess)
        {
            return NULL;
        }
        for (long k = 0; k < nblayers; k++)
        {
            cudaMemcpy(s_d_m + k * msize * msize, h_masters[k],
                       screen_bytes, cudaMemcpyHostToDevice);
        }
        s_cached_mbytes = total_mbytes;
        s_cached_h_masters = h_masters;
    }
    return s_d_m;
}

static int atmturb_cuda_sync_masters(
    const atmturb_cuda_sim_params_t *p)
{
    return (atmturb_cuda_sync_device_masters(p->msize, p->nblayers,
                                            p->h_masters) != NULL) ? 0 : -1;
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
        (float) params->Scoeff, 0,
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

__global__ static void atmturb_cuda_geom_reduce_sum(
    const float *d_in,
    int          npix,
    float       *d_sum)
{
    float val = 0.0f;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = idx; i < npix; i += stride)
    {
        val += d_in[i];
    }

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

__global__ static void atmturb_cuda_geom_sub_mean(
    float       *d_out,
    int          npix,
    const float *d_sum)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= npix)
    {
        return;
    }
    float mean = (*d_sum) / (float) npix;
    d_out[idx] -= mean;
}

struct atmturb_cuda_geom_stream
{
    long          nblayers;
    long          msize;
    long          pup_size;
    float         Scoeff;
    cudaStream_t  stream;
    float        *d_pha;
    float        *d_amp;
    float        *d_spha;
    float        *d_samp;
    float        *d_sum;
    size_t        npix;
};

/**
 * atmturb_cuda_geom_stream_init - Initialize GPU context for 2D geometric streaming
 * @params: Simulation parameters and input screen pointers.
 *
 * Return: Allocated stream context, or NULL on failure.
 */
atmturb_cuda_geom_stream_t *atmturb_cuda_geom_stream_init(
    const atmturb_cuda_sim_params_t *params)
{
    if (params == NULL || params->pup_size <= 0 || params->msize <= 0)
    {
        return NULL;
    }
    if (atmturb_cuda_sync_masters(params) != 0 ||
        atmturb_cuda_sync_kinematics(params) != 0)
    {
        return NULL;
    }

    atmturb_cuda_geom_stream_t *ctx =
        (atmturb_cuda_geom_stream_t *) calloc(1, sizeof(atmturb_cuda_geom_stream_t));
    if (ctx == NULL)
    {
        return NULL;
    }

    ctx->nblayers = params->nblayers;
    ctx->msize    = params->msize;
    ctx->pup_size = params->pup_size;
    ctx->Scoeff   = (float) params->Scoeff;
    ctx->npix     = (size_t) (params->pup_size * params->pup_size);
    size_t nbytes = sizeof(float) * ctx->npix;

    if (cudaStreamCreate(&ctx->stream) != cudaSuccess)
    {
        free(ctx);
        return NULL;
    }

    if (cudaMalloc((void **) &ctx->d_pha, nbytes) != cudaSuccess ||
        cudaMalloc((void **) &ctx->d_amp, nbytes) != cudaSuccess ||
        cudaMalloc((void **) &ctx->d_spha, nbytes) != cudaSuccess ||
        cudaMalloc((void **) &ctx->d_samp, nbytes) != cudaSuccess ||
        cudaMalloc((void **) &ctx->d_sum, sizeof(float)) != cudaSuccess)
    {
        atmturb_cuda_geom_stream_free(ctx);
        return NULL;
    }

    return ctx;
}

/**
 * atmturb_cuda_geom_stream_render_step - Render a single geometric frame on GPU
 * @ctx: Persistent stream context.
 * @t: Simulation frame index.
 * @pha: Destination host buffer for primary phase.
 * @amp: Destination host buffer for primary amplitude (or NULL).
 * @spha: Destination host buffer for secondary phase (or NULL).
 * @samp: Destination host buffer for secondary amplitude (or NULL).
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_cuda_geom_stream_render_step(
    atmturb_cuda_geom_stream_t *ctx,
    long                        t,
    float                      *pha,
    float                      *amp,
    float                      *spha,
    float                      *samp)
{
    if (ctx == NULL || pha == NULL)
    {
        return -1;
    }

    dim3 block(32, 8);
    dim3 grid(((int) ctx->pup_size + block.x - 1) / block.x,
              ((int) ctx->pup_size + block.y - 1) / block.y,
              1);

    atmturb_render_kernel<<<grid, block, 0, ctx->stream>>>(
        s_d_m, (int) ctx->msize, (int) ctx->pup_size,
        1, (int) ctx->nblayers, ctx->Scoeff, (int) t,
        ctx->d_pha, (amp != NULL) ? ctx->d_amp : NULL,
        (spha != NULL) ? ctx->d_spha : NULL, (samp != NULL) ? ctx->d_samp : NULL);

    cudaMemsetAsync(ctx->d_sum, 0, sizeof(float), ctx->stream);
    dim3 r_block(256);
    dim3 r_grid(((int) ctx->npix + 255) / 256);
    atmturb_cuda_geom_reduce_sum<<<r_grid, r_block, 0, ctx->stream>>>(
        ctx->d_pha, (int) ctx->npix, ctx->d_sum);
    atmturb_cuda_geom_sub_mean<<<r_grid, r_block, 0, ctx->stream>>>(
        ctx->d_pha, (int) ctx->npix, ctx->d_sum);

    if (spha != NULL)
    {
        cudaMemsetAsync(ctx->d_sum, 0, sizeof(float), ctx->stream);
        atmturb_cuda_geom_reduce_sum<<<r_grid, r_block, 0, ctx->stream>>>(
            ctx->d_spha, (int) ctx->npix, ctx->d_sum);
        atmturb_cuda_geom_sub_mean<<<r_grid, r_block, 0, ctx->stream>>>(
            ctx->d_spha, (int) ctx->npix, ctx->d_sum);
    }

    size_t nbytes = sizeof(float) * ctx->npix;
    cudaMemcpyAsync(pha, ctx->d_pha, nbytes, cudaMemcpyDeviceToHost, ctx->stream);
    if (amp != NULL)
    {
        cudaMemcpyAsync(amp, ctx->d_amp, nbytes, cudaMemcpyDeviceToHost, ctx->stream);
    }
    if (spha != NULL)
    {
        cudaMemcpyAsync(spha, ctx->d_spha, nbytes, cudaMemcpyDeviceToHost, ctx->stream);
    }
    if (samp != NULL)
    {
        cudaMemcpyAsync(samp, ctx->d_samp, nbytes, cudaMemcpyDeviceToHost, ctx->stream);
    }

    cudaStreamSynchronize(ctx->stream);
    return 0;
}

/**
 * atmturb_cuda_geom_stream_free - Release GPU geometric streaming context
 * @ctx: Stream context to release.
 */
void atmturb_cuda_geom_stream_free(
    atmturb_cuda_geom_stream_t *ctx)
{
    if (ctx == NULL)
    {
        return;
    }
    if (ctx->d_pha != NULL)  cudaFree(ctx->d_pha);
    if (ctx->d_amp != NULL)  cudaFree(ctx->d_amp);
    if (ctx->d_spha != NULL) cudaFree(ctx->d_spha);
    if (ctx->d_samp != NULL) cudaFree(ctx->d_samp);
    if (ctx->d_sum != NULL)  cudaFree(ctx->d_sum);
    if (ctx->stream != 0)    cudaStreamDestroy(ctx->stream);
    free(ctx);
}

