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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "atmturb_cuda.h"
#include "atmturb_cuda_dma_ring.h"
#include "atmturb_cuda_fresnel.h"
#include "atmturb_cuda_geom_ext.h"
#include "atmturb_cuda_piston.h"
#include "atmturb_cuda_psf.h"
#include "atmturb_cuda_screens.h"
#include "atmturb_linpred_solve.h"
extern "C" {
#include "atmturb_simd.h"
}

// Maximum supported layers in constant memory cache and cached rolling screens
#define ATMTURB_CUDA_MAX_LAYERS  64
#define ATMTURB_CUDA_MAX_SCREENS 512

typedef struct
{
    float vx[ATMTURB_CUDA_MAX_LAYERS];
    float vy[ATMTURB_CUDA_MAX_LAYERS];
    float x0[ATMTURB_CUDA_MAX_LAYERS];
    float y0[ATMTURB_CUDA_MAX_LAYERS];
    float xs0[ATMTURB_CUDA_MAX_LAYERS];
    float ys0[ATMTURB_CUDA_MAX_LAYERS];
    float weights[ATMTURB_CUDA_MAX_LAYERS];
    float weights_s[ATMTURB_CUDA_MAX_LAYERS];
} atmturb_cuda_const_kin_t;

__constant__ atmturb_cuda_const_kin_t c_kin;

// Persistent GPU memory and 2D layered texture caches
static float              *s_d_m              = NULL;
static size_t              s_cached_mbytes    = 0;
static long                s_cached_nblayers  = 0;
static const float        *s_cached_screen_ptrs[ATMTURB_CUDA_MAX_SCREENS];
static cudaArray_t         s_master_array     = NULL;
static cudaTextureObject_t s_master_tex       = 0;
static int                 s_has_chrom_shift  = 0;

static float       *s_d_pha[2]          = {NULL, NULL};
static float       *s_d_spha[2]         = {NULL, NULL};
static size_t       s_cached_pha_bytes  = 0;
static size_t       s_cached_spha_bytes = 0;
static cudaStream_t s_batch_stream      = 0;

/**
 * atmturb_cuda_sample_bili - Bilinear lookup on master screen layer @k via texture or __ldg
 */
__device__ static inline float atmturb_cuda_sample_bili(
    const float        *d_masters,
    cudaTextureObject_t master_tex,
    int                 msize,
    int                 is_pow2,
    int                 mask,
    size_t              msize_sq,
    int                 k,
    float               cur_x,
    float               cur_y)
{
    int   ix = __float2int_rd(cur_x);
    int   iy = __float2int_rd(cur_y);
    float fx = cur_x - (float) ix;
    float fy = cur_y - (float) iy;
    int   ix0, iy0, ix1, iy1;

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

    float v00, v10, v01, v11;
    if (master_tex != 0)
    {
        float tx0 = (float) ix0 + 0.5f;
        float tx1 = (float) ix1 + 0.5f;
        float ty0 = (float) iy0 + 0.5f;
        float ty1 = (float) iy1 + 0.5f;
        v00 = tex2DLayered<float>(master_tex, tx0, ty0, k);
        v10 = tex2DLayered<float>(master_tex, tx1, ty0, k);
        v01 = tex2DLayered<float>(master_tex, tx0, ty1, k);
        v11 = tex2DLayered<float>(master_tex, tx1, ty1, k);
    }
    else
    {
        const float *screen = d_masters + (size_t) k * msize_sq;
        v00 = __ldg(&screen[(size_t) iy0 * (size_t) msize + (size_t) ix0]);
        v10 = __ldg(&screen[(size_t) iy0 * (size_t) msize + (size_t) ix1]);
        v01 = __ldg(&screen[(size_t) iy1 * (size_t) msize + (size_t) ix0]);
        v11 = __ldg(&screen[(size_t) iy1 * (size_t) msize + (size_t) ix1]);
    }

    float top = v00 + fx * (v10 - v00);
    float bot = v01 + fx * (v11 - v01);
    return top + fy * (bot - top);
}

/**
 * atmturb_render_kernel - GPU kernel synthesizing multi-layer extruded wavefront frames
 */
__global__ static void atmturb_render_kernel(
    const float        *d_masters,
    cudaTextureObject_t master_tex,
    int                 msize,
    int                 pup_size,
    int                 nbframes,
    int                 nblayers,
    int                 has_chrom_shift,
    int                 t_offset,
    float              *d_pha,
    float              *d_spha)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;
    int jj = blockIdx.y * blockDim.y + threadIdx.y;
    if (ii >= pup_size || jj >= pup_size)
    {
        return;
    }

    int    is_pow2  = ((msize & (msize - 1)) == 0);
    int    mask     = msize - 1;
    size_t msize_sq = (size_t) msize * (size_t) msize;
    size_t pup_sq   = (size_t) pup_size * (size_t) pup_size;
    int    chunk    = (nbframes + (int) gridDim.z - 1) / (int) gridDim.z;
    int    f_beg    = (int) blockIdx.z * chunk;
    int    f_end    = (f_beg + chunk > nbframes) ? nbframes : (f_beg + chunk);

    for (int f = f_beg; f < f_end; f++)
    {
        int   t          = t_offset + f;
        float total_pha  = 0.0f;
        float total_spha = 0.0f;

        for (int k = 0; k < nblayers; k++)
        {
            float dx  = (float) t * c_kin.vx[k] + (float) ii;
            float dy  = (float) t * c_kin.vy[k] + (float) jj;
            float val = atmturb_cuda_sample_bili(
                d_masters, master_tex, msize, is_pow2, mask, msize_sq,
                k, c_kin.x0[k] + dx, c_kin.y0[k] + dy);
            total_pha += c_kin.weights[k] * val;

            if (d_spha != NULL)
            {
                if (has_chrom_shift &&
                    (c_kin.xs0[k] != c_kin.x0[k] || c_kin.ys0[k] != c_kin.y0[k]))
                {
                    float val_s = atmturb_cuda_sample_bili(
                        d_masters, master_tex, msize, is_pow2, mask, msize_sq,
                        k, c_kin.xs0[k] + dx, c_kin.ys0[k] + dy);
                    total_spha += c_kin.weights_s[k] * val_s;
                }
                else
                {
                    total_spha += c_kin.weights_s[k] * val;
                }
            }
        }

        size_t out_idx = (size_t) f * pup_sq + (size_t) jj * (size_t) pup_size + (size_t) ii;
        d_pha[out_idx] = total_pha;
        if (d_spha != NULL)
        {
            d_spha[out_idx] = total_spha;
        }
    }
}

/**
 * atmturb_cuda_free_master_tex - Destroy cached 2D layered texture object and array
 */
static void atmturb_cuda_free_master_tex(void)
{
    if (s_master_tex != 0)
    {
        cudaDestroyTextureObject(s_master_tex);
        s_master_tex = 0;
    }
    if (s_master_array != NULL)
    {
        cudaFreeArray(s_master_array);
        s_master_array = NULL;
    }
}

/**
 * atmturb_cuda_rebuild_master_tex - Populate 2D layered Morton-tiled texture from s_d_m
 * @msize: Master screen dimension in pixels.
 * @nblayers: Number of screens stored in s_d_m.
 */
static void atmturb_cuda_rebuild_master_tex(
    long msize,
    long nblayers)
{
    atmturb_cuda_free_master_tex();
    if (s_d_m == NULL || msize <= 0 || msize > 16384 || nblayers <= 0 || nblayers > 2048)
    {
        return;
    }

    cudaChannelFormatDesc desc = cudaCreateChannelDesc<float>();
    cudaExtent            ext  = make_cudaExtent((size_t) msize, (size_t) msize, (size_t) nblayers);
    if (cudaMalloc3DArray(&s_master_array, &desc, ext, cudaArrayLayered) != cudaSuccess)
    {
        s_master_array = NULL;
        return;
    }

    cudaMemcpy3DParms cp;
    memset(&cp, 0, sizeof(cp));
    cp.srcPtr = make_cudaPitchedPtr(s_d_m, sizeof(float) * (size_t) msize,
                                    (size_t) msize, (size_t) msize);
    cp.dstArray = s_master_array;
    cp.extent   = ext;
    cp.kind     = cudaMemcpyDeviceToDevice;
    if (cudaMemcpy3D(&cp) != cudaSuccess)
    {
        atmturb_cuda_free_master_tex();
        return;
    }

    cudaResourceDesc res_desc;
    memset(&res_desc, 0, sizeof(res_desc));
    res_desc.resType         = cudaResourceTypeArray;
    res_desc.res.array.array = s_master_array;

    cudaTextureDesc tex_desc;
    memset(&tex_desc, 0, sizeof(tex_desc));
    tex_desc.addressMode[0]   = cudaAddressModeClamp;
    tex_desc.addressMode[1]   = cudaAddressModeClamp;
    tex_desc.filterMode       = cudaFilterModePoint;
    tex_desc.readMode         = cudaReadModeElementType;
    tex_desc.normalizedCoords = 0;

    if (cudaCreateTextureObject(&s_master_tex, &res_desc, &tex_desc, NULL) != cudaSuccess)
    {
        s_master_tex = 0;
        atmturb_cuda_free_master_tex();
    }
}

unsigned long long atmturb_cuda_get_master_tex(void)
{
    return (unsigned long long) s_master_tex;
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
    if (s_batch_stream != 0)
    {
        cudaStreamSynchronize(s_batch_stream);
        cudaStreamDestroy(s_batch_stream);
        s_batch_stream = 0;
    }
    atmturb_cuda_dma_ring_cleanup();
    atmturb_cuda_geom_ext_cleanup();
    atmturb_cuda_fresnel_batch_cleanup();
    atmturb_cuda_linpred_cleanup();
    atmturb_cuda_screens_cleanup();
    atmturb_cuda_psf_cleanup();
    atmturb_cuda_free_master_tex();
    if (s_d_m != NULL)
    {
        cudaFree(s_d_m);
        s_d_m = NULL;
    }
    for (int s = 0; s < 2; s++)
    {
        cudaFree(s_d_pha[s]);
        s_d_pha[s] = NULL;
        cudaFree(s_d_spha[s]);
        s_d_spha[s] = NULL;
    }
    s_cached_mbytes = 0;
    s_cached_nblayers = 0;
    memset(s_cached_screen_ptrs, 0, sizeof(s_cached_screen_ptrs));
    s_cached_pha_bytes  = 0;
    s_cached_spha_bytes = 0;
}

/**
 * atmturb_cuda_alloc_master_cache - Allocate linear device master screen cache for GPU synthesis
 * @msize: Master screen dimension in pixels.
 * @nblayers: Total number of master screens across layers.
 *
 * Return: Device buffer pointer, or NULL on allocation failure.
 */
float *atmturb_cuda_alloc_master_cache(
    long msize,
    long nblayers)
{
    if (msize <= 0 || nblayers <= 0 || nblayers > ATMTURB_CUDA_MAX_SCREENS)
    {
        return NULL;
    }

    size_t total_mbytes = sizeof(float) * (size_t) msize * (size_t) msize * (size_t) nblayers;
    if (s_d_m == NULL || s_cached_mbytes != total_mbytes)
    {
        atmturb_cuda_free_master_tex();
        if (s_d_m != NULL)
        {
            cudaFree(s_d_m);
            s_d_m = NULL;
        }
        if (cudaMalloc((void **) &s_d_m, total_mbytes) != cudaSuccess)
        {
            s_cached_mbytes = 0;
            s_cached_nblayers = 0;
            return NULL;
        }
        s_cached_mbytes = total_mbytes;
    }
    s_cached_nblayers = 0;
    memset(s_cached_screen_ptrs, 0, sizeof(s_cached_screen_ptrs));
    return s_d_m;
}

/**
 * atmturb_cuda_commit_master_cache - Bind GPU-synthesized master screens and build 2D texture
 * @msize: Master screen dimension in pixels.
 * @nblayers: Total number of master screens across layers.
 * @h_masters: Host pointer tags matching dcimg screen buffers for cache-hit tracking.
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_cuda_commit_master_cache(
    long                msize,
    long                nblayers,
    const float *const *h_masters)
{
    if (s_d_m == NULL || h_masters == NULL || msize <= 0 ||
        nblayers <= 0 || nblayers > ATMTURB_CUDA_MAX_SCREENS)
    {
        return -1;
    }

    for (long k = 0; k < nblayers; k++)
    {
        s_cached_screen_ptrs[k] = h_masters[k];
    }
    s_cached_nblayers = nblayers;
    atmturb_cuda_rebuild_master_tex(msize, nblayers);
    return 0;
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
    if (msize <= 0 || nblayers <= 0 || nblayers > ATMTURB_CUDA_MAX_SCREENS || h_masters == NULL)
    {
        return NULL;
    }

    size_t screen_bytes = sizeof(float) * (size_t) msize * (size_t) msize;
    size_t total_mbytes = screen_bytes * (size_t) nblayers;
    int    need_upload  = 0;

    if (s_d_m == NULL || s_cached_mbytes != total_mbytes || s_cached_nblayers != nblayers)
    {
        need_upload = 1;
    }
    else
    {
        for (long k = 0; k < nblayers; k++)
        {
            if (s_cached_screen_ptrs[k] != h_masters[k])
            {
                need_upload = 1;
                break;
            }
        }
    }

    if (need_upload)
    {
        atmturb_cuda_free_master_tex();
        if (s_d_m != NULL)
        {
            cudaFree(s_d_m);
            s_d_m = NULL;
        }
        if (cudaMalloc((void **) &s_d_m, total_mbytes) != cudaSuccess)
        {
            s_cached_mbytes = 0;
            s_cached_nblayers = 0;
            return NULL;
        }
        for (long k = 0; k < nblayers; k++)
        {
            if (cudaMemcpy(s_d_m + (size_t) k * (size_t) msize * (size_t) msize, h_masters[k],
                           screen_bytes, cudaMemcpyHostToDevice) != cudaSuccess)
            {
                cudaFree(s_d_m);
                s_d_m = NULL;
                s_cached_mbytes = 0;
                s_cached_nblayers = 0;
                return NULL;
            }
            s_cached_screen_ptrs[k] = h_masters[k];
        }
        s_cached_mbytes = total_mbytes;
        s_cached_nblayers = nblayers;
        atmturb_cuda_rebuild_master_tex(msize, nblayers);
    }
    return s_d_m;
}

/**
 * atmturb_cuda_sync_masters - Upload simulation master phase screens to cached device memory
 * @p: Simulation parameters containing master screen dimensions and host pointers.
 *
 * Return: 0 on success, -1 on device allocation or transfer failure.
 */
static int atmturb_cuda_sync_masters(
    const atmturb_cuda_sim_params_t *p)
{
    return (atmturb_cuda_sync_device_masters(p->msize, p->nblayers,
                                             p->h_masters) != NULL) ? 0 : -1;
}

/**
 * atmturb_cuda_sync_kinematics - Upload layer velocity and position offsets to constant memory
 * @p: Simulation parameters.
 * @stream: CUDA stream for asynchronous constant-memory upload (or 0 for synchronous).
 *
 * Return: 0 on success, -1 if layers exceed constant buffer capacity.
 */
static int atmturb_cuda_sync_kinematics(
    const atmturb_cuda_sim_params_t *p,
    cudaStream_t                     stream)
{
    if (p->nblayers > ATMTURB_CUDA_MAX_LAYERS)
    {
        return -1;
    }

    atmturb_cuda_const_kin_t h_kin;
    memset(&h_kin, 0, sizeof(h_kin));
    int has_shift = 0;

    for (long k = 0; k < p->nblayers; k++)
    {
        h_kin.vx[k]        = (float) p->vxpix[k];
        h_kin.vy[k]        = (float) p->vypix[k];
        h_kin.weights[k]   = (float) sqrt(p->cn2[k]);
        h_kin.weights_s[k] = (p->cn2_s != NULL) ? (float) sqrt(p->cn2_s[k])
                                                : (h_kin.weights[k] * (float) p->Scoeff);
        h_kin.x0[k]        = (p->x0 != NULL) ? (float) p->x0[k] : (float) p->vxpix[k];
        h_kin.y0[k]        = (p->y0 != NULL) ? (float) p->y0[k] : (float) p->vypix[k];
        h_kin.xs0[k]       = (p->xs0 != NULL) ? (float) p->xs0[k] : h_kin.x0[k];
        h_kin.ys0[k]       = (p->ys0 != NULL) ? (float) p->ys0[k] : h_kin.y0[k];
        if (h_kin.xs0[k] != h_kin.x0[k] || h_kin.ys0[k] != h_kin.y0[k])
        {
            has_shift = 1;
        }
    }

    s_has_chrom_shift = has_shift;
    cudaError_t err = (stream != 0)
                          ? cudaMemcpyToSymbolAsync(c_kin, &h_kin, sizeof(h_kin), 0,
                                                    cudaMemcpyHostToDevice, stream)
                          : cudaMemcpyToSymbol(c_kin, &h_kin, sizeof(h_kin));
    return (err == cudaSuccess) ? 0 : -1;
}

/**
 * atmturb_cuda_fill_amp_pair - Multi-threaded host 1.0f fill overlapped with PCIe DMA
 */
static void atmturb_cuda_fill_amp_pair(
    float *amp,
    float *samp,
    size_t npix)
{
    if (amp == NULL && samp == NULL)
    {
        return;
    }
#ifdef _OPENMP
    if (npix >= 65536 && !omp_in_parallel())
    {
        int max_t = omp_get_max_threads();
        #pragma omp parallel num_threads((max_t > 8) ? 8 : max_t)
        {
            int  tid = omp_get_thread_num();
            int  nt  = omp_get_num_threads();
            long i0  = (long) ((npix * (size_t) tid) / (size_t) nt);
            long i1  = (long) ((npix * (size_t) (tid + 1)) / (size_t) nt);
            if (i1 > i0)
            {
                if (amp != NULL)
                {
                    atmturb_init_phase_amp(NULL, amp + i0, i1 - i0);
                }
                if (samp != NULL)
                {
                    atmturb_init_phase_amp(NULL, samp + i0, i1 - i0);
                }
            }
        }
        return;
    }
#endif
    if (amp != NULL)
    {
        atmturb_init_phase_amp(NULL, amp, (long) npix);
    }
    if (samp != NULL)
    {
        atmturb_init_phase_amp(NULL, samp, (long) npix);
    }
}

/**
 * atmturb_cuda_sync_outputs - Allocate double-buffered device phase output slots lazily
 */
static int atmturb_cuda_sync_outputs(
    size_t slot_bytes,
    int    need_spha)
{
    if (s_d_pha[0] == NULL || s_cached_pha_bytes < slot_bytes)
    {
        for (int s = 0; s < 2; s++)
        {
            cudaFree(s_d_pha[s]);
            s_d_pha[s] = NULL;
            if (cudaMalloc((void **) &s_d_pha[s], slot_bytes) != cudaSuccess)
            {
                s_cached_pha_bytes = 0;
                return -1;
            }
        }
        s_cached_pha_bytes = slot_bytes;
    }

    if (need_spha && (s_d_spha[0] == NULL || s_cached_spha_bytes < slot_bytes))
    {
        for (int s = 0; s < 2; s++)
        {
            cudaFree(s_d_spha[s]);
            s_d_spha[s] = NULL;
            if (cudaMalloc((void **) &s_d_spha[s], slot_bytes) != cudaSuccess)
            {
                s_cached_spha_bytes = 0;
                return -1;
            }
        }
        s_cached_spha_bytes = slot_bytes;
    }
    return 0;
}

/**
 * atmturb_cuda_make_batch_grid - Compute 3D launch grid dimensions for batch extrusion
 */
static dim3 atmturb_cuda_make_batch_grid(
    long pup_size,
    long nbframes,
    dim3 block)
{
    int gx        = ((int) pup_size + (int) block.x - 1) / (int) block.x;
    int gy        = ((int) pup_size + (int) block.y - 1) / (int) block.y;
    int xy_blocks = gx * gy;
    int min_gz    = (512 + xy_blocks - 1) / xy_blocks;
    int chunk_gz  = ((int) nbframes + 15) / 16;
    int grid_z    = (min_gz > chunk_gz) ? min_gz : chunk_gz;
    if (grid_z > (int) nbframes)
    {
        grid_z = (int) nbframes;
    }
    if (grid_z > 65535)
    {
        grid_z = 65535;
    }
    return dim3(gx, gy, grid_z);
}

/**
 * atmturb_wfs_render_frames_cuda - CUDA GPU accelerated wavefront time-series rendering
 */
int atmturb_wfs_render_frames_cuda(
    const atmturb_cuda_sim_params_t *params,
    atmturb_cuda_sim_outputs_t      *outputs)
{
    if (s_batch_stream == 0 &&
        cudaStreamCreateWithFlags(&s_batch_stream, cudaStreamNonBlocking) != cudaSuccess)
    {
        return -1;
    }
    if (atmturb_cuda_sync_masters(params) != 0 ||
        atmturb_cuda_sync_kinematics(params, s_batch_stream) != 0)
    {
        return -1;
    }

    size_t pup_sq    = (size_t) params->pup_size * (size_t) params->pup_size;
    size_t npix      = (size_t) params->nbframes * pup_sq;
    size_t out_bytes = sizeof(float) * npix;
    int    need_spha = (outputs->spha != NULL) ? 1 : 0;

    if (out_bytes <= (size_t) (512 * 1024))
    {
        if (atmturb_cuda_sync_outputs(out_bytes, need_spha) != 0)
        {
            return -1;
        }
        dim3 block(32, 8);
        dim3 grid = atmturb_cuda_make_batch_grid(params->pup_size, params->nbframes, block);
        atmturb_render_kernel<<<grid, block, 0, s_batch_stream>>>(
            s_d_m, s_master_tex, (int) params->msize, (int) params->pup_size,
            (int) params->nbframes, (int) params->nblayers, s_has_chrom_shift, 0,
            s_d_pha[0], need_spha ? s_d_spha[0] : NULL);
        if (outputs->pha != NULL)
        {
            cudaMemcpyAsync(outputs->pha, s_d_pha[0], out_bytes,
                            cudaMemcpyDeviceToHost, s_batch_stream);
        }
        if (need_spha)
        {
            cudaMemcpyAsync(outputs->spha, s_d_spha[0], out_bytes,
                            cudaMemcpyDeviceToHost, s_batch_stream);
        }
        atmturb_cuda_fill_amp_pair(outputs->amp, outputs->samp, npix);
        return (cudaStreamSynchronize(s_batch_stream) == cudaSuccess) ? 0 : -1;
    }

    long   batch_b    = (params->nbframes < 32) ? params->nbframes : 32;
    size_t slot_bytes = sizeof(float) * (size_t) batch_b * pup_sq;
    if (atmturb_cuda_sync_outputs(slot_bytes, need_spha) != 0 ||
        atmturb_cuda_dma_ring_prepare(slot_bytes) != 0)
    {
        return -1;
    }

    dim3 block(32, 8);
    int  ret        = 0;
    int  amp_filled = 0;
    for (long t0 = 0, chunk = 0; t0 < params->nbframes; t0 += batch_b, chunk++)
    {
        int  slot  = (int) (chunk & 1);
        long count = (params->nbframes - t0 < batch_b) ? (params->nbframes - t0) : batch_b;
        if (atmturb_cuda_dma_ring_drain(slot) != 0)
        {
            ret = -1;
            break;
        }

        dim3 grid = atmturb_cuda_make_batch_grid(params->pup_size, count, block);
        atmturb_render_kernel<<<grid, block, 0, s_batch_stream>>>(
            s_d_m, s_master_tex, (int) params->msize, (int) params->pup_size,
            (int) count, (int) params->nblayers, s_has_chrom_shift, (int) t0,
            s_d_pha[slot], need_spha ? s_d_spha[slot] : NULL);

        size_t             dst_off  = (size_t) t0 * pup_sq;
        size_t             cbytes   = sizeof(float) * (size_t) count * pup_sq;
        const float *const d_src[4] = {
            s_d_pha[slot], NULL, need_spha ? s_d_spha[slot] : NULL, NULL
        };
        float *const h_dst[4] = {
            (outputs->pha != NULL) ? (outputs->pha + dst_off) : NULL, NULL,
            need_spha ? (outputs->spha + dst_off) : NULL, NULL
        };
        if (atmturb_cuda_dma_ring_enqueue(slot, s_batch_stream, d_src, h_dst, cbytes) != 0)
        {
            ret = -1;
            break;
        }
        if (chunk == 1)
        {
            atmturb_cuda_fill_amp_pair(outputs->amp, outputs->samp, npix);
            amp_filled = 1;
        }
    }
    if (!amp_filled && ret == 0)
    {
        atmturb_cuda_fill_amp_pair(outputs->amp, outputs->samp, npix);
    }
    if (atmturb_cuda_dma_ring_flush() != 0)
    {
        ret = -1;
    }
    return ret;
}

struct atmturb_cuda_geom_stream
{
    long          nblayers;
    long          msize;
    long          pup_size;
    int           has_chrom_shift;
    cudaStream_t  streams[2];
    cudaEvent_t   ev_done[2];
    float        *d_buf;
    float        *d_pha[2];
    float        *d_spha[2];
    float        *d_sum[2];
    size_t        npix;
};

/**
 * atmturb_cuda_geom_stream_init - Initialize GPU context for 2D geometric streaming
 */
atmturb_cuda_geom_stream_t *atmturb_cuda_geom_stream_init(
    const atmturb_cuda_sim_params_t *params)
{
    if (params == NULL || params->pup_size <= 0 || params->msize <= 0)
    {
        return NULL;
    }
    if (atmturb_cuda_sync_masters(params) != 0 ||
        atmturb_cuda_sync_kinematics(params, 0) != 0)
    {
        return NULL;
    }

    atmturb_cuda_geom_stream_t *ctx =
        (atmturb_cuda_geom_stream_t *) calloc(1, sizeof(atmturb_cuda_geom_stream_t));
    if (ctx == NULL)
    {
        return NULL;
    }

    ctx->nblayers        = params->nblayers;
    ctx->msize           = params->msize;
    ctx->pup_size        = params->pup_size;
    ctx->has_chrom_shift = s_has_chrom_shift;
    ctx->npix            = (size_t) params->pup_size * (size_t) params->pup_size;
    size_t stride        = (ctx->npix + 63) & ~(size_t) 63;
    size_t tot_floats    = 4 * stride + 128;

    if (cudaMalloc((void **) &ctx->d_buf, sizeof(float) * tot_floats) != cudaSuccess)
    {
        atmturb_cuda_geom_stream_free(ctx);
        return NULL;
    }

    for (int s = 0; s < 2; s++)
    {
        if (cudaStreamCreateWithFlags(&ctx->streams[s], cudaStreamNonBlocking) != cudaSuccess ||
            cudaEventCreateWithFlags(&ctx->ev_done[s], cudaEventDisableTiming) != cudaSuccess)
        {
            atmturb_cuda_geom_stream_free(ctx);
            return NULL;
        }
        ctx->d_pha[s]  = ctx->d_buf + (size_t) (2 * s) * stride;
        ctx->d_spha[s] = ctx->d_buf + (size_t) (2 * s + 1) * stride;
        ctx->d_sum[s]  = ctx->d_buf + 4 * stride + (size_t) (s * 64);
    }

    return ctx;
}

int atmturb_cuda_geom_stream_render_launch(
    atmturb_cuda_geom_stream_t *ctx,
    long                        t,
    float                      *pha,
    float                      *spha,
    int                         slot)
{
    if (ctx == NULL || pha == NULL)
    {
        return -1;
    }
    int          s  = slot & 1;
    cudaStream_t st = ctx->streams[s];
    dim3         block(32, 8);
    dim3         grid(((int) ctx->pup_size + block.x - 1) / block.x,
                      ((int) ctx->pup_size + block.y - 1) / block.y, 1);

    atmturb_render_kernel<<<grid, block, 0, st>>>(
        s_d_m, s_master_tex, (int) ctx->msize, (int) ctx->pup_size,
        1, (int) ctx->nblayers, ctx->has_chrom_shift, (int) t,
        ctx->d_pha[s], (spha != NULL) ? ctx->d_spha[s] : NULL);

    atmturb_cuda_remove_piston_launch(ctx->d_pha[s], ctx->d_sum[s], (int) ctx->npix, st);
    if (spha != NULL)
    {
        atmturb_cuda_remove_piston_launch(ctx->d_spha[s], ctx->d_sum[s], (int) ctx->npix, st);
    }

    size_t nbytes = sizeof(float) * ctx->npix;
    cudaMemcpyAsync(pha, ctx->d_pha[s], nbytes, cudaMemcpyDeviceToHost, st);
    if (spha != NULL)
    {
        cudaMemcpyAsync(spha, ctx->d_spha[s], nbytes, cudaMemcpyDeviceToHost, st);
    }
    cudaEventRecord(ctx->ev_done[s], st);
    return 0;
}

int atmturb_cuda_geom_stream_render_wait(
    atmturb_cuda_geom_stream_t *ctx,
    int                         slot)
{
    if (ctx == NULL)
    {
        return -1;
    }
    return (cudaEventSynchronize(ctx->ev_done[slot & 1]) == cudaSuccess) ? 0 : -1;
}

int atmturb_cuda_geom_stream_render_step(
    atmturb_cuda_geom_stream_t *ctx,
    long                        t,
    float                      *pha,
    float                      *amp,
    float                      *spha,
    float                      *samp)
{
    if (atmturb_cuda_geom_stream_render_launch(ctx, t, pha, spha, 0) != 0)
    {
        return -1;
    }
    atmturb_cuda_fill_amp_pair(amp, samp, ctx->npix);
    return atmturb_cuda_geom_stream_render_wait(ctx, 0);
}

void atmturb_cuda_geom_stream_free(
    atmturb_cuda_geom_stream_t *ctx)
{
    if (ctx == NULL)
    {
        return;
    }
    for (int s = 0; s < 2; s++)
    {
        if (ctx->streams[s] != 0)
        {
            cudaStreamSynchronize(ctx->streams[s]);
            cudaStreamDestroy(ctx->streams[s]);
        }
        if (ctx->ev_done[s] != 0)
        {
            cudaEventDestroy(ctx->ev_done[s]);
        }
    }
    cudaFree(ctx->d_buf);
    free(ctx);
}

/**
 * atmturb_cuda_host_register - Page-lock an existing host buffer for asynchronous CUDA DMA
 * @ptr: Host memory pointer (e.g. ImageStreamIO shared memory region).
 * @size: Size of region in bytes.
 *
 * Return: 0 on success, -1 on invalid arguments or CUDA registration failure.
 */
int atmturb_cuda_host_register(
    void   *ptr,
    size_t  size)
{
    if (ptr == NULL || size == 0)
    {
        return -1;
    }
    cudaError_t err = cudaHostRegister(ptr, size, cudaHostRegisterDefault);
    return (err == cudaSuccess) ? 0 : -1;
}

/**
 * atmturb_cuda_host_unregister - Unpin a previously registered host memory buffer
 * @ptr: Host memory pointer (safe to pass NULL).
 */
void atmturb_cuda_host_unregister(
    void *ptr)
{
    if (ptr != NULL)
    {
        cudaHostUnregister(ptr);
    }
}

/**
 * atmturb_cuda_host_alloc - Allocate page-locked pinned host memory for zero-copy/DMA transfers
 * @ptr: Output pointer receiving allocated pinned memory address.
 * @size: Requested allocation size in bytes.
 *
 * Return: 0 on success, -1 on invalid arguments or CUDA allocation failure.
 */
int atmturb_cuda_host_alloc(
    void  **ptr,
    size_t  size)
{
    if (ptr == NULL || size == 0)
    {
        return -1;
    }
    cudaError_t err = cudaMallocHost(ptr, size);
    if (err != cudaSuccess)
    {
        *ptr = NULL;
        return -1;
    }
    return 0;
}

/**
 * atmturb_cuda_host_free - Free page-locked host memory allocated by atmturb_cuda_host_alloc()
 * @ptr: Pinned host memory pointer (safe to pass NULL).
 */
void atmturb_cuda_host_free(
    void *ptr)
{
    if (ptr != NULL)
    {
        cudaFreeHost(ptr);
    }
}

