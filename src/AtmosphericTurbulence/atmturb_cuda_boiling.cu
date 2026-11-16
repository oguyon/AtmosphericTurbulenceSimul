// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_boiling.cu
 * @brief   CUDA/cuFFT Fourier AR(1)/AR(2) boiling synthesis and layered texture slice updates
 */

#include <cuda_runtime.h>
#include <cufft.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "atmturb_cuda.h"
#include "atmturb_cuda_boiling.h"

struct atmturb_cuda_boiling
{
    int                             nlayers;
    long                            msize;
    long                            hcols;
    int                             ar_order;
    int                             total_vram_screens;
    cudaStream_t                    stream;
    cufftHandle                     c2r_plan;
    float                          *d_cache;
    cufftComplex                   *d_c2r_work;
    float                          *d_screen_out;
    atmturb_cuda_boil_layer_spec_t *specs;
    cufftComplex                  **d_spec_cur;
    cufftComplex                  **d_spec_aux;
};

/**
 * atmturb_cuda_boil_mix64 - Stateless SplitMix64 bit-avalanching mixer
 */
__device__ static inline uint64_t atmturb_cuda_boil_mix64(
    uint64_t z)
{
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

/**
 * atmturb_cuda_boil_gauss - Draw standard Gaussian pair at stream index and draw offset
 */
__device__ static inline float2 atmturb_cuda_boil_gauss(
    uint64_t seed,
    uint64_t stream_id,
    int      draw_idx)
{
    const uint64_t gamma    = 0x9e3779b97f4a7c15ULL;
    const double   inv_2p53 = 1.0 / 9007199254740992.0;
    uint64_t       s0       = seed ^ ((stream_id + 1ULL) * 0xd1b54a32d192ed03ULL);
    uint64_t       rng0     = atmturb_cuda_boil_mix64(s0 + 2ULL * gamma);
    uint64_t       r1       = atmturb_cuda_boil_mix64(rng0 + (uint64_t) (2 * draw_idx + 1) * gamma);
    uint64_t       r2       = atmturb_cuda_boil_mix64(rng0 + (uint64_t) (2 * draw_idx + 2) * gamma);

    float u1 = (float) (((double) (r1 >> 11) + 1.0) * inv_2p53);
    if (u1 > 0.99999994f)
    {
        u1 = 0.99999994f;
    }
    float u2 = (float) ((double) (r2 >> 11) * inv_2p53);
    float rad = sqrtf(-2.0f * logf(u1));
    float sn, cs;
    __sincosf(6.2831853071795864769f * u2, &sn, &cs);
    return make_float2(rad * cs, rad * sn);
}

/**
 * atmturb_cuda_boil_omr2 - Cancellation-free evaluation of 1 - exp(-2 * eta)
 */
__device__ static inline float atmturb_cuda_boil_omr2(
    float eta,
    float rp2)
{
    float v = (eta < 0.02f) ? (2.0f * eta * (1.0f - eta * (1.0f - 0.6666667f * eta)))
                            : (1.0f - rp2);
    return (v > 0.0f) ? v : 0.0f;
}

/**
 * atmturb_cuda_boil_mode_update - Apply stationary init or AR(1)/AR(2) step to one mode
 */
__device__ static inline cufftComplex atmturb_cuda_boil_mode_update(
    cufftComplex *d_cur,
    cufftComplex *d_aux,
    size_t        idx,
    float2        z0,
    float2        z1,
    float         eta,
    int           ar_order,
    int           is_init)
{
    float rp   = expf(-eta);
    float rp2  = rp * rp;
    float omr2 = atmturb_cuda_boil_omr2(eta, rp2);
    cufftComplex out;

    if (is_init)
    {
        if (ar_order == 2 && d_aux != NULL)
        {
            float inv_1pr2 = 1.0f / (1.0f + rp2);
            float c0       = 2.0f * rp * inv_1pr2;
            float c1       = omr2 * inv_1pr2;
            out.x          = c0 * z0.x + c1 * z1.x;
            out.y          = c0 * z0.y + c1 * z1.y;
            d_aux[idx]     = make_cuFloatComplex(c1 * (rp * z0.x + z1.x),
                                                 c1 * (rp * z0.y + z1.y));
        }
        else
        {
            out = make_cuFloatComplex(z0.x, z0.y);
        }
        d_cur[idx] = out;
        return out;
    }

    cufftComplex cur = d_cur[idx];
    if (ar_order == 2 && d_aux != NULL)
    {
        cufftComplex aux = d_aux[idx];
        float b   = omr2 * sqrtf(omr2 / (1.0f + rp2));
        float ux  = rp * aux.x + b * z0.x;
        float uy  = rp * aux.y + b * z0.y;
        out.x     = rp * cur.x + ux;
        out.y     = rp * cur.y + uy;
        d_aux[idx] = make_cuFloatComplex(ux, uy);
    }
    else
    {
        float b = sqrtf(omr2);
        out.x   = rp * cur.x + b * z0.x;
        out.y   = rp * cur.y + b * z0.y;
    }
    d_cur[idx] = out;
    return out;
}

/**
 * atmturb_cuda_boil_interior_kernel - Update interior Hermitian columns 1..half_end-1
 */
__global__ static void atmturb_cuda_boil_interior_kernel(
    cufftComplex                  *d_cur,
    cufftComplex                  *d_aux,
    cufftComplex                  *d_work,
    int                            msize,
    int                            ar_order,
    int                            is_init,
    uint64_t                       seed0,
    uint64_t                       seed1,
    atmturb_cuda_boil_layer_spec_t sp)
{
    int ii       = (int) (blockIdx.x * blockDim.x + threadIdx.x) + 1;
    int jj       = (int) (blockIdx.y * blockDim.y + threadIdx.y);
    int hcols    = msize / 2 + 1;
    int half_end = (msize % 2 == 0) ? (msize / 2) : hcols;
    if (ii >= half_end || jj >= msize)
    {
        return;
    }

    int    fy  = (jj <= msize / 2) ? jj : (jj - msize);
    float  r2  = (float) (ii * ii + fy * fy);
    float  amp = sp.sqrt_pref_pair * powf(r2 + sp.k0_sq, -11.0f / 12.0f);
    if (sp.inv_two_km2 > 0.0f)
    {
        amp *= expf(-r2 * sp.inv_two_km2);
    }

    float2 g0 = atmturb_cuda_boil_gauss(seed0, (uint64_t) jj, ii - 1);
    float2 g1 = (is_init && ar_order == 2)
                    ? atmturb_cuda_boil_gauss(seed1, (uint64_t) jj, ii - 1)
                    : make_float2(0.0f, 0.0f);
    float2 z0  = make_float2(amp * g0.x, amp * g0.y);
    float2 z1  = make_float2(amp * g1.x, amp * g1.y);
    size_t idx = (size_t) jj * (size_t) hcols + (size_t) ii;
    d_work[idx] = atmturb_cuda_boil_mode_update(
        d_cur, d_aux, idx, z0, z1, sp.alpha_dt * sqrtf(r2), ar_order, is_init);
}

/**
 * atmturb_cuda_boil_axis_kernel - Update self-conjugate columns {0, msize/2} and mirror rows
 */
__global__ static void atmturb_cuda_boil_axis_kernel(
    cufftComplex                  *d_cur,
    cufftComplex                  *d_aux,
    cufftComplex                  *d_work,
    int                            msize,
    int                            ii,
    int                            ar_order,
    int                            is_init,
    uint64_t                       seed0,
    uint64_t                       seed1,
    atmturb_cuda_boil_layer_spec_t sp)
{
    int jj       = (int) (blockIdx.x * blockDim.x + threadIdx.x);
    int half_lim = msize / 2;
    int hcols    = half_lim + 1;
    if (jj > half_lim)
    {
        return;
    }

    size_t idx = (size_t) jj * (size_t) hcols + (size_t) ii;
    if (ii == 0 && jj == 0)
    {
        d_cur[idx]  = make_cuFloatComplex(0.0f, 0.0f);
        d_work[idx] = make_cuFloatComplex(0.0f, 0.0f);
        return;
    }

    int      col_tag  = (ii == 0) ? 0 : 1;
    int      draw_idx = (ii == 0) ? (jj - 1) : jj;
    uint64_t st_id    = (uint64_t) (msize + col_tag);
    int      is_nyq   = (jj == 0) || ((msize % 2 == 0) && (jj == half_lim));
    float    r2       = (float) (ii * ii + jj * jj);
    float    spref    = is_nyq ? sp.sqrt_pref_real : sp.sqrt_pref_pair;
    float    amp      = spref * powf(r2 + sp.k0_sq, -11.0f / 12.0f);
    if (sp.inv_two_km2 > 0.0f)
    {
        amp *= expf(-r2 * sp.inv_two_km2);
    }

    float2 g0 = atmturb_cuda_boil_gauss(seed0, st_id, draw_idx);
    float2 g1 = (is_init && ar_order == 2)
                    ? atmturb_cuda_boil_gauss(seed1, st_id, draw_idx)
                    : make_float2(0.0f, 0.0f);
    float2 z0 = make_float2(amp * g0.x, is_nyq ? 0.0f : (amp * g0.y));
    float2 z1 = make_float2(amp * g1.x, is_nyq ? 0.0f : (amp * g1.y));

    cufftComplex val = atmturb_cuda_boil_mode_update(
        d_cur, d_aux, idx, z0, z1, sp.alpha_dt * sqrtf(r2), ar_order, is_init);
    d_work[idx] = val;

    if (!is_nyq && jj > 0)
    {
        size_t mir  = (size_t) (msize - jj) * (size_t) hcols + (size_t) ii;
        cufftComplex conj_v = make_cuFloatComplex(val.x, -val.y);
        d_cur[mir]  = conj_v;
        d_work[mir] = conj_v;
    }
}

atmturb_cuda_boiling_t *atmturb_cuda_boiling_init(
    int                                   nlayers,
    long                                  msize,
    int                                   ar_order,
    int                                   total_vram_screens,
    const atmturb_cuda_boil_layer_spec_t *specs)
{
    if (nlayers <= 0 || msize < 64 || msize > 8192 || total_vram_screens <= 0 ||
        specs == NULL || !atmturb_cuda_device_available())
    {
        return NULL;
    }

    atmturb_cuda_boiling_t *ctx = (atmturb_cuda_boiling_t *) calloc(1, sizeof(*ctx));
    if (ctx == NULL)
    {
        return NULL;
    }
    ctx->nlayers            = nlayers;
    ctx->msize              = msize;
    ctx->hcols              = msize / 2 + 1;
    ctx->ar_order           = ar_order;
    ctx->total_vram_screens = total_vram_screens;
    ctx->d_cache            = atmturb_cuda_alloc_master_cache(msize, total_vram_screens);
    ctx->specs              = (atmturb_cuda_boil_layer_spec_t *) calloc(
        (size_t) nlayers, sizeof(*ctx->specs));
    ctx->d_spec_cur = (cufftComplex **) calloc((size_t) nlayers, sizeof(cufftComplex *));
    ctx->d_spec_aux = (cufftComplex **) calloc((size_t) nlayers, sizeof(cufftComplex *));
    if (!ctx->d_cache || !ctx->specs || !ctx->d_spec_cur || !ctx->d_spec_aux)
    {
        atmturb_cuda_boiling_free(ctx);
        return NULL;
    }
    memcpy(ctx->specs, specs, sizeof(*specs) * (size_t) nlayers);

    size_t ncomp = (size_t) (msize * ctx->hcols);
    size_t npix  = (size_t) (msize * msize);
    if (cudaStreamCreate(&ctx->stream) != cudaSuccess ||
        cudaMalloc((void **) &ctx->d_c2r_work, sizeof(cufftComplex) * ncomp) != cudaSuccess ||
        cudaMalloc((void **) &ctx->d_screen_out, sizeof(float) * npix) != cudaSuccess ||
        cufftPlan2d(&ctx->c2r_plan, (int) msize, (int) msize, CUFFT_C2R) != CUFFT_SUCCESS ||
        cufftSetStream(ctx->c2r_plan, ctx->stream) != CUFFT_SUCCESS)
    {
        atmturb_cuda_boiling_free(ctx);
        return NULL;
    }

    for (int k = 0; k < nlayers; k++)
    {
        if (!specs[k].active)
        {
            continue;
        }
        if (cudaMalloc((void **) &ctx->d_spec_cur[k],
                       sizeof(cufftComplex) * ncomp) != cudaSuccess ||
            (ar_order == 2 &&
             cudaMalloc((void **) &ctx->d_spec_aux[k],
                        sizeof(cufftComplex) * ncomp) != cudaSuccess))
        {
            atmturb_cuda_boiling_free(ctx);
            return NULL;
        }
    }
    return ctx;
}

int atmturb_cuda_boiling_step_layer(
    atmturb_cuda_boiling_t *ctx,
    int                     layer_idx,
    int                     slot,
    int                     is_init,
    const uint64_t          seeds[2],
    float                  *h_dst_opt)
{
    if (ctx == NULL || layer_idx < 0 || layer_idx >= ctx->nlayers || slot < 0)
    {
        return -1;
    }
    const atmturb_cuda_boil_layer_spec_t sp = ctx->specs[layer_idx];
    if (!sp.active)
    {
        return 0;
    }

    int  msize = (int) ctx->msize;
    dim3 b_int(32, 8, 1);
    dim3 g_int((unsigned int) ((msize / 2 + 31) / 32), (unsigned int) ((msize + 7) / 8), 1);
    atmturb_cuda_boil_interior_kernel<<<g_int, b_int, 0, ctx->stream>>>(
        ctx->d_spec_cur[layer_idx], ctx->d_spec_aux[layer_idx], ctx->d_c2r_work,
        msize, ctx->ar_order, is_init, seeds[0], seeds[1], sp);

    unsigned int g_ax = (unsigned int) ((msize / 2 + 128) / 128);
    atmturb_cuda_boil_axis_kernel<<<g_ax, 128, 0, ctx->stream>>>(
        ctx->d_spec_cur[layer_idx], ctx->d_spec_aux[layer_idx], ctx->d_c2r_work,
        msize, 0, ctx->ar_order, is_init, seeds[0], seeds[1], sp);
    if (msize % 2 == 0)
    {
        atmturb_cuda_boil_axis_kernel<<<g_ax, 128, 0, ctx->stream>>>(
            ctx->d_spec_cur[layer_idx], ctx->d_spec_aux[layer_idx], ctx->d_c2r_work,
            msize, msize / 2, ctx->ar_order, is_init, seeds[0], seeds[1], sp);
    }

    if (cufftExecC2R(ctx->c2r_plan, ctx->d_c2r_work, ctx->d_screen_out) != CUFFT_SUCCESS ||
        atmturb_cuda_update_master_slice(ctx->msize, sp.vram_offset + slot,
                                         ctx->d_screen_out, ctx->stream) != 0)
    {
        return -1;
    }
    if (h_dst_opt != NULL &&
        cudaMemcpyAsync(h_dst_opt, ctx->d_screen_out,
                        sizeof(float) * (size_t) (ctx->msize * ctx->msize),
                        cudaMemcpyDeviceToHost, ctx->stream) != cudaSuccess)
    {
        return -1;
    }
    return (cudaStreamSynchronize(ctx->stream) == cudaSuccess) ? 0 : -1;
}

int atmturb_cuda_boiling_upload_static(
    atmturb_cuda_boiling_t *ctx,
    int                     vram_slot,
    const float            *h_src)
{
    if (ctx == NULL || ctx->d_cache == NULL || h_src == NULL || vram_slot < 0)
    {
        return -1;
    }
    size_t np = (size_t) ctx->msize * (size_t) ctx->msize;
    return (cudaMemcpy(ctx->d_cache + (size_t) vram_slot * np, h_src,
                       sizeof(float) * np, cudaMemcpyHostToDevice) == cudaSuccess) ? 0 : -1;
}

int atmturb_cuda_boiling_commit_cache(
    atmturb_cuda_boiling_t *ctx,
    const float *const     *vram_host_tags)
{
    if (ctx == NULL || vram_host_tags == NULL)
    {
        return -1;
    }
    return atmturb_cuda_commit_master_cache(ctx->msize, ctx->total_vram_screens, vram_host_tags);
}

void atmturb_cuda_boiling_free(
    atmturb_cuda_boiling_t *ctx)
{
    if (ctx == NULL)
    {
        return;
    }
    if (ctx->stream != 0)
    {
        cudaStreamSynchronize(ctx->stream);
        cudaStreamDestroy(ctx->stream);
    }
    if (ctx->c2r_plan != 0)
    {
        cufftDestroy(ctx->c2r_plan);
    }
    cudaFree(ctx->d_c2r_work);
    cudaFree(ctx->d_screen_out);
    for (int k = 0; k < ctx->nlayers; k++)
    {
        if (ctx->d_spec_cur != NULL)
        {
            cudaFree(ctx->d_spec_cur[k]);
        }
        if (ctx->d_spec_aux != NULL)
        {
            cudaFree(ctx->d_spec_aux[k]);
        }
    }
    free(ctx->d_spec_cur);
    free(ctx->d_spec_aux);
    free(ctx->specs);
    free(ctx);
}
