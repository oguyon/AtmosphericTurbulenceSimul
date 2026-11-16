// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_boiling.h
 * @brief   CUDA/cuFFT Fourier AR(1)/AR(2) boiling synthesis and layered texture slice updates
 */

#ifndef ATMTURB_CUDA_BOILING_H
#define ATMTURB_CUDA_BOILING_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct atmturb_cuda_boiling atmturb_cuda_boiling_t;

/**
 * struct atmturb_cuda_boil_layer_spec_t - Per-layer spectral parameters for GPU boiling
 * @k0_sq: Squared outer-scale spatial frequency index (msize / L0_pix)^2.
 * @sqrt_pref_pair: Spectral amplitude prefactor for conjugate-paired modes.
 * @sqrt_pref_real: Spectral amplitude prefactor for real Nyquist modes.
 * @inv_two_km2: Inner-scale Gaussian damping factor 0.5 / km^2 (0 if disabled).
 * @alpha_dt: Dimensionless decay scale 2*pi*sigma_boil*dt_kf / (msize*dx_master).
 * @vram_offset: Base slice index of this layer in the GPU master screen cache.
 * @active: Non-zero if this layer undergoes Fourier boiling.
 */
typedef struct
{
    float k0_sq;
    float sqrt_pref_pair;
    float sqrt_pref_real;
    float inv_two_km2;
    float alpha_dt;
    int   vram_offset;
    int   active;
} atmturb_cuda_boil_layer_spec_t;

#ifdef HAVE_CUDA

/**
 * atmturb_cuda_boiling_init - Allocate GPU AR state, cuFFT C2R plan, and VRAM master cache
 * @nlayers: Number of atmospheric layers.
 * @msize: Master screen linear dimension [pixels].
 * @ar_order: Autoregressive order (1 for AR(1), 2 for critically damped AR(2)).
 * @total_vram_screens: Total screen slots across all layers in GPU master cache.
 * @specs: Per-layer spectral parameter array (nlayers elements).
 *
 * Return: Allocated GPU boiling context, or NULL on unsupported config or CUDA failure.
 */
atmturb_cuda_boiling_t *atmturb_cuda_boiling_init(
    int                                   nlayers,
    long                                  msize,
    int                                   ar_order,
    int                                   total_vram_screens,
    const atmturb_cuda_boil_layer_spec_t *specs);

/**
 * atmturb_cuda_boiling_step_layer - Synthesize initial or next AR keyframe on GPU
 * @ctx: Active GPU boiling context.
 * @layer_idx: Layer index (0..nlayers - 1).
 * @slot: Target screen slot within this layer (0..nscreens - 1).
 * @is_init: Non-zero for stationary t=0 initialization, 0 for incremental AR step.
 * @seeds: Two-element PRNG seed array (seeds[0] primary/Z0, seeds[1] Z1 for AR(2) init).
 * @h_dst_opt: Optional host destination buffer (msize * msize floats, NULL to skip D2H).
 *
 * Return: 0 on success, -1 on CUDA/cuFFT execution failure.
 */
int atmturb_cuda_boiling_step_layer(
    atmturb_cuda_boiling_t *ctx,
    int                     layer_idx,
    int                     slot,
    int                     is_init,
    const uint64_t          seeds[2],
    float                  *h_dst_opt);

/**
 * atmturb_cuda_boiling_upload_static - Upload a non-boiling static layer screen to VRAM cache
 * @ctx: Active GPU boiling context.
 * @vram_slot: Global VRAM slot index.
 * @h_src: Host float buffer (msize * msize elements).
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_cuda_boiling_upload_static(
    atmturb_cuda_boiling_t *ctx,
    int                     vram_slot,
    const float            *h_src);

/**
 * atmturb_cuda_boiling_commit_cache - Finalize initial VRAM master screen cache and 2D texture
 * @ctx: Active GPU boiling context.
 * @vram_host_tags: Array of total_vram_screens host pointer tags for cache-hit tracking.
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_cuda_boiling_commit_cache(
    atmturb_cuda_boiling_t *ctx,
    const float *const     *vram_host_tags);

/**
 * atmturb_cuda_boiling_free - Release GPU AR spectra, scratch buffers, and cuFFT plan
 * @ctx: GPU boiling context to free (safe on NULL).
 */
void atmturb_cuda_boiling_free(
    atmturb_cuda_boiling_t *ctx);

#else /* !HAVE_CUDA */

static inline atmturb_cuda_boiling_t *atmturb_cuda_boiling_init(
    int                                   nlayers,
    long                                  msize,
    int                                   ar_order,
    int                                   total_vram_screens,
    const atmturb_cuda_boil_layer_spec_t *specs)
{
    (void) nlayers;
    (void) msize;
    (void) ar_order;
    (void) total_vram_screens;
    (void) specs;
    return NULL;
}

static inline int atmturb_cuda_boiling_step_layer(
    atmturb_cuda_boiling_t *ctx,
    int                     layer_idx,
    int                     slot,
    int                     is_init,
    const uint64_t          seeds[2],
    float                  *h_dst_opt)
{
    (void) ctx;
    (void) layer_idx;
    (void) slot;
    (void) is_init;
    (void) seeds;
    (void) h_dst_opt;
    return -1;
}

static inline int atmturb_cuda_boiling_upload_static(
    atmturb_cuda_boiling_t *ctx,
    int                     vram_slot,
    const float            *h_src)
{
    (void) ctx;
    (void) vram_slot;
    (void) h_src;
    return -1;
}

static inline int atmturb_cuda_boiling_commit_cache(
    atmturb_cuda_boiling_t *ctx,
    const float *const     *vram_host_tags)
{
    (void) ctx;
    (void) vram_host_tags;
    return -1;
}

static inline void atmturb_cuda_boiling_free(
    atmturb_cuda_boiling_t *ctx)
{
    (void) ctx;
}

#endif /* HAVE_CUDA */

#ifdef __cplusplus
}
#endif

#endif /* ATMTURB_CUDA_BOILING_H */
