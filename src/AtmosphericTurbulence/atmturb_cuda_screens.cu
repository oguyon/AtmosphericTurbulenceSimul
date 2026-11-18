// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_screens.cu
 * @brief   CUDA/cuFFT batch synthesis of von Karman master turbulence phase screens
 *
 * Uses counter-addressed SplitMix64 sub-streams to reproduce the exact realization sequence
 * of `atmturb_fill_spectrum_row_float()` directly in VRAM, executes batched 2D forward cuFFTs,
 * and populates the persistent GPU master screen + 2D layered texture cache without H2D copies.
 */

#include <cuda_runtime.h>
#include <cufft.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#include "atmturb_cuda.h"
#include "atmturb_cuda_screens.h"

#define ATMTURB_SCR_CUDA_MAX_BATCH 4

typedef struct
{
    uint64_t base_seed;
    float    k0_sq;
    float    sqrt_pref;
    float    inv_two_km2;
    int      slot_a;
    int      slot_b;
} atmturb_cuda_scr_item_t;

typedef struct
{
    atmturb_cuda_scr_item_t items[ATMTURB_SCR_CUDA_MAX_BATCH];
} atmturb_cuda_scr_pack_t;

static long          s_scr_msize     = 0;
static int           s_scr_batch_max = 0;
static cudaStream_t  s_scr_stream    = 0;
static cufftHandle   s_scr_plan_b    = 0;
static cufftHandle   s_scr_plan_1    = 0;
static cufftComplex *s_d_scr_complex = NULL;
static float        *s_d_scr_tmp     = NULL;

/**
 * atmturb_cuda_mix64 - Stateless SplitMix64 output mixer
 * @z: Pre-incremented 64-bit state word.
 *
 * Return: Mixed 64-bit pseudorandom word.
 */
__host__ __device__ static inline uint64_t atmturb_cuda_mix64(
    uint64_t z)
{
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

/**
 * atmturb_cuda_resolve_seed - Resolve 0 user seed via real-time clock matching CPU helper
 * @seed: User seed (0 = clock-derived).
 *
 * Return: Non-zero 64-bit seed.
 */
static uint64_t atmturb_cuda_resolve_seed(
    uint64_t seed)
{
    if (seed != 0)
    {
        return seed;
    }
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    uint64_t s = (uint64_t) ts.tv_sec * 1000000000ULL + (uint64_t) ts.tv_nsec;
    return atmturb_cuda_mix64(s + 0x9e3779b97f4a7c15ULL) | 1ULL;
}

/**
 * atmturb_cuda_fill_spectrum_kernel - Synthesize batched 2D von Karman complex spectra in VRAM
 * @d_spec: Output batched complex spectrum array (`count * size * size` elements).
 * @size: Grid linear dimension in pixels (multiple of 16).
 * @pack: Batch item spectral parameters and PRNG base seeds.
 */
__global__ static void atmturb_cuda_fill_spectrum_kernel(
    cufftComplex            *d_spec,
    int                      size,
    atmturb_cuda_scr_pack_t  pack)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;
    int jj = blockIdx.y * blockDim.y + threadIdx.y;
    int b  = (int) blockIdx.z;
    if (ii >= size || jj >= size)
    {
        return;
    }

    const atmturb_cuda_scr_item_t it = pack.items[b];
    size_t idx = ((size_t) b * (size_t) size + (size_t) jj) * (size_t) size + (size_t) ii;
    if (ii == 0 && jj == 0)
    {
        d_spec[idx] = make_cuFloatComplex(0.0f, 0.0f);
        return;
    }

    const uint64_t gamma = 0x9e3779b97f4a7c15ULL;
    uint64_t s_row = it.base_seed ^ (((uint64_t) jj + 1ULL) * 0xd1b54a32d192ed03ULL);
    uint64_t rng0  = atmturb_cuda_mix64(s_row + 2ULL * gamma);
    uint64_t r1    = atmturb_cuda_mix64(rng0 + (uint64_t) (2 * ii + 1) * gamma);
    uint64_t r2_u  = atmturb_cuda_mix64(rng0 + (uint64_t) (2 * ii + 2) * gamma);

    const double inv_2p53 = 1.0 / 9007199254740992.0;
    double d1 = (double) ((r1 >> 11) + 1ULL) * inv_2p53;
    double d2 = (double) (r2_u >> 11) * inv_2p53;

    int   half = size >> 1;
    int   fx   = (ii < half) ? ii : (ii - size);
    int   fy   = (jj < half) ? jj : (jj - size);
    float fy2  = (float) (fy * fy);
    float r2   = (float) (fx * fx) + fy2;
    float amp  = it.sqrt_pref * powf(r2 + it.k0_sq, -11.0f / 12.0f);
    if (it.inv_two_km2 > 0.0f)
    {
        amp *= expf(-r2 * it.inv_two_km2);
    }

    if ((fy2 >= 256.0f) || (ii >= 16 && ii < size - 16))
    {
        float u1 = (float) d1;
        if (u1 > 0.99999994f)
        {
            u1 = 0.99999994f;
        }
        float r_bm = sqrtf(-2.0f * logf(u1));
        float s, c;
        __sincosf(6.2831853071795864769f * (float) d2, &s, &c);
        d_spec[idx] = make_cuFloatComplex(amp * r_bm * c, amp * r_bm * s);
    }
    else
    {
        double r_bm = sqrt(-2.0 * log(d1));
        double s, c;
        sincos(6.283185307179586476925286766559 * d2, &s, &c);
        d_spec[idx] = make_cuFloatComplex(amp * (float) (r_bm * c), amp * (float) (r_bm * s));
    }
}

/**
 * atmturb_cuda_deinterleave_screens_kernel - Split complex DFT outputs into real/imag screens
 * @d_spec: Input batched complex spatial screens (`count * npix`).
 * @d_cache: Optional persistent GPU master cache (`total_vram_screens * npix`, or NULL).
 * @d_tmp: Contiguous device staging buffer (`2 * count * npix`) when `d_cache == NULL`.
 * @npix: Pixels per screen (`size * size`).
 * @pack: Batch item destination VRAM slot indices.
 */
__global__ static void atmturb_cuda_deinterleave_screens_kernel(
    const cufftComplex      *d_spec,
    float                   *d_cache,
    float                   *d_tmp,
    size_t                   npix,
    atmturb_cuda_scr_pack_t  pack)
{
    size_t idx = (size_t) blockIdx.x * (size_t) blockDim.x + (size_t) threadIdx.x;
    int    b   = (int) blockIdx.y;
    if (idx >= npix)
    {
        return;
    }

    cufftComplex val = d_spec[(size_t) b * npix + idx];
    if (d_cache != NULL)
    {
        int sa = pack.items[b].slot_a;
        int sb = pack.items[b].slot_b;
        if (sa >= 0)
        {
            d_cache[(size_t) sa * npix + idx] = val.x;
        }
        if (sb >= 0)
        {
            d_cache[(size_t) sb * npix + idx] = val.y;
        }
    }
    else if (d_tmp != NULL)
    {
        d_tmp[(size_t) (2 * b) * npix + idx]     = val.x;
        d_tmp[(size_t) (2 * b + 1) * npix + idx] = val.y;
    }
}

void atmturb_cuda_screens_cleanup(void)
{
    if (s_scr_stream != 0)
    {
        cudaStreamSynchronize(s_scr_stream);
        cudaStreamDestroy(s_scr_stream);
        s_scr_stream = 0;
    }
    if (s_scr_plan_b != 0)
    {
        cufftDestroy(s_scr_plan_b);
        s_scr_plan_b = 0;
    }
    if (s_scr_plan_1 != 0)
    {
        cufftDestroy(s_scr_plan_1);
        s_scr_plan_1 = 0;
    }
    if (s_d_scr_complex != NULL)
    {
        cudaFree(s_d_scr_complex);
        s_d_scr_complex = NULL;
    }
    if (s_d_scr_tmp != NULL)
    {
        cudaFree(s_d_scr_tmp);
        s_d_scr_tmp = NULL;
    }
    s_scr_msize     = 0;
    s_scr_batch_max = 0;
}

/**
 * atmturb_cuda_screens_ensure_ctx - Lazily allocate cuFFT plans and scratch buffers for msize
 * @msize: Master screen dimension in pixels.
 * @need_tmp: 1 if temporary real/imag staging buffer is required (`d_cache == NULL`), 0 otherwise.
 *
 * Return: 0 on success, -1 on CUDA/cuFFT failure.
 */
static int atmturb_cuda_screens_ensure_ctx(
    long msize,
    int  need_tmp)
{
    int batch_b = (msize >= 4096) ? 2 : ATMTURB_SCR_CUDA_MAX_BATCH;
    size_t npix = (size_t) msize * (size_t) msize;

    if (s_scr_msize != msize || s_scr_batch_max != batch_b || s_d_scr_complex == NULL)
    {
        atmturb_cuda_screens_cleanup();
        if (cudaStreamCreate(&s_scr_stream) != cudaSuccess ||
            cudaMalloc((void **) &s_d_scr_complex,
                       sizeof(cufftComplex) * npix * (size_t) batch_b) != cudaSuccess)
        {
            atmturb_cuda_screens_cleanup();
            return -1;
        }

        int n_dims[2] = {(int) msize, (int) msize};
        int inembed[2] = {(int) msize, (int) msize};
        int dist = (int) npix;
        if (cufftPlanMany(&s_scr_plan_b, 2, n_dims, inembed, 1, dist,
                          inembed, 1, dist, CUFFT_C2C, batch_b) != CUFFT_SUCCESS ||
            cufftPlan2d(&s_scr_plan_1, (int) msize, (int) msize, CUFFT_C2C) != CUFFT_SUCCESS ||
            cufftSetStream(s_scr_plan_b, s_scr_stream) != CUFFT_SUCCESS ||
            cufftSetStream(s_scr_plan_1, s_scr_stream) != CUFFT_SUCCESS)
        {
            atmturb_cuda_screens_cleanup();
            return -1;
        }
        s_scr_msize     = msize;
        s_scr_batch_max = batch_b;
    }

    if (need_tmp && s_d_scr_tmp == NULL)
    {
        if (cudaMalloc((void **) &s_d_scr_tmp,
                       sizeof(float) * npix * (size_t) (2 * batch_b)) != cudaSuccess)
        {
            return -1;
        }
    }
    return 0;
}

/**
 * atmturb_cuda_pack_batch_specs - Populate spectral constants for a chunk of screen pairs
 */
static void atmturb_cuda_pack_batch_specs(
    const atmturb_screen_batch_req_t *req,
    int                               start_p,
    int                               count,
    atmturb_cuda_scr_pack_t          *pack)
{
    memset(pack, 0, sizeof(*pack));
    for (int b = 0; b < count; b++)
    {
        int p = start_p + b;
        const atmturb_screen_spec_t *sp = &req->specs[p];
        float sz       = (float) sp->size;
        float r0_pix   = (sp->r0_pix > 0.0) ? (float) sp->r0_pix : powf(6.88f, 0.6f);
        float k0       = (sp->L0_pix > 0.0) ? (sz / (float) sp->L0_pix) : 0.0f;
        float km_scale = 5.92f / 6.2831853071795864769f;
        float km       = (sp->l0_pix > 0.0) ? (km_scale * sz / (float) sp->l0_pix) : 0.0f;
        float pref     = 0.023f * powf(r0_pix, -5.0f / 3.0f) * powf(sz, 5.0f / 3.0f);

        pack->items[b].base_seed   = atmturb_cuda_resolve_seed(sp->seed);
        pack->items[b].k0_sq       = k0 * k0;
        pack->items[b].sqrt_pref   = sqrtf(pref);
        pack->items[b].inv_two_km2 = (km > 0.0f) ? (0.5f / (km * km)) : 0.0f;
        pack->items[b].slot_a      = (req->vram_slot_a != NULL) ? req->vram_slot_a[p] : -1;
        pack->items[b].slot_b      = (req->vram_slot_b != NULL) ? req->vram_slot_b[p] : -1;
    }
}

/**
 * atmturb_cuda_run_screen_chunk - Execute spectrum synthesis, cuFFT, deinterleave, and D2H copy
 */
static int atmturb_cuda_run_screen_chunk(
    const atmturb_screen_batch_req_t *req,
    float                            *d_cache,
    long                              msize,
    int                               start_p,
    int                               count)
{
    atmturb_cuda_scr_pack_t pack;
    atmturb_cuda_pack_batch_specs(req, start_p, count, &pack);

    size_t npix = (size_t) msize * (size_t) msize;
    dim3 b_spec(32, 8, 1);
    dim3 g_spec(((int) msize + 31) / 32, ((int) msize + 7) / 8, (unsigned int) count);
    atmturb_cuda_fill_spectrum_kernel<<<g_spec, b_spec, 0, s_scr_stream>>>(
        s_d_scr_complex, (int) msize, pack);

    if (count == s_scr_batch_max)
    {
        if (cufftExecC2C(s_scr_plan_b, s_d_scr_complex, s_d_scr_complex,
                         CUFFT_FORWARD) != CUFFT_SUCCESS)
        {
            return -1;
        }
    }
    else
    {
        for (int b = 0; b < count; b++)
        {
            cufftComplex *ptr = s_d_scr_complex + (size_t) b * npix;
            if (cufftExecC2C(s_scr_plan_1, ptr, ptr, CUFFT_FORWARD) != CUFFT_SUCCESS)
            {
                return -1;
            }
        }
    }

    dim3 b_de(256, 1, 1);
    dim3 g_de((unsigned int) ((npix + 255) / 256), (unsigned int) count, 1);
    atmturb_cuda_deinterleave_screens_kernel<<<g_de, b_de, 0, s_scr_stream>>>(
        s_d_scr_complex, d_cache, s_d_scr_tmp, npix, pack);

    size_t scr_bytes = sizeof(float) * npix;
    for (int b = 0; b < count; b++)
    {
        int p = start_p + b;
        const float *src_a = (d_cache != NULL)
                                 ? (d_cache + (size_t) pack.items[b].slot_a * npix)
                                 : (s_d_scr_tmp + (size_t) (2 * b) * npix);
        const float *src_b = (d_cache != NULL)
                                 ? ((pack.items[b].slot_b >= 0)
                                        ? (d_cache + (size_t) pack.items[b].slot_b * npix)
                                        : NULL)
                                 : (s_d_scr_tmp + (size_t) (2 * b + 1) * npix);
        if (req->out_a[p] != NULL &&
            cudaMemcpyAsync(req->out_a[p], src_a, scr_bytes,
                            cudaMemcpyDeviceToHost, s_scr_stream) != cudaSuccess)
        {
            return -1;
        }
        if (req->out_b != NULL && req->out_b[p] != NULL && src_b != NULL &&
            cudaMemcpyAsync(req->out_b[p], src_b, scr_bytes,
                            cudaMemcpyDeviceToHost, s_scr_stream) != cudaSuccess)
        {
            return -1;
        }
    }
    return 0;
}

int atmturb_cuda_generate_screen_pairs(
    const atmturb_screen_batch_req_t *req)
{
    if (req == NULL || req->specs == NULL || req->out_a == NULL || req->npairs <= 0)
    {
        return -1;
    }

    long msize = req->specs[0].size;
    if (msize < 256 || (msize & 15) != 0 || msize > 8192)
    {
        return -1;
    }
    for (int p = 0; p < req->npairs; p++)
    {
        if (req->specs[p].size != msize || req->specs[p].precision != 0)
        {
            return -1;
        }
    }

    float *d_cache = NULL;
    if (req->total_vram_screens > 0 && req->vram_slot_a != NULL && req->vram_host_tags != NULL)
    {
        d_cache = atmturb_cuda_alloc_master_cache(msize, req->total_vram_screens);
    }

    if (atmturb_cuda_screens_ensure_ctx(msize, (d_cache == NULL) ? 1 : 0) != 0)
    {
        return -1;
    }

    for (int p = 0; p < req->npairs; p += s_scr_batch_max)
    {
        int count = req->npairs - p;
        if (count > s_scr_batch_max)
        {
            count = s_scr_batch_max;
        }
        if (atmturb_cuda_run_screen_chunk(req, d_cache, msize, p, count) != 0)
        {
            return -1;
        }
    }

    if (cudaStreamSynchronize(s_scr_stream) != cudaSuccess)
    {
        return -1;
    }
    if (d_cache != NULL)
    {
        atmturb_cuda_commit_master_cache(msize, req->total_vram_screens, req->vram_host_tags);
    }
    return 0;
}
