// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_screens.h
 * @brief   CUDA/cuFFT batch synthesis of von Karman master turbulence phase screens
 */

#ifndef ATMTURB_CUDA_SCREENS_H
#define ATMTURB_CUDA_SCREENS_H

#include "AtmosphericTurbulence.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * struct atmturb_screen_batch_req - Batch master screen synthesis request descriptor
 * @specs: Array of `npairs` screen specifications (all sharing identical size and precision).
 * @out_a: Array of `npairs` host output buffers for real screen channel (non-NULL).
 * @out_b: Array of `npairs` host output buffers for imaginary screen channel (entries may be NULL).
 * @vram_slot_a: Optional `npairs` VRAM cache slot indices for `out_a` (or NULL).
 * @vram_slot_b: Optional `npairs` VRAM cache slot indices for `out_b` (-1 to skip, or NULL).
 * @npairs: Number of complex Fourier screen pairs to synthesize.
 * @total_vram_screens: Total screen count to commit in persistent GPU master cache (0 to skip).
 * @vram_host_tags: Array of `total_vram_screens` host pointer tags for GPU cache-hit tracking.
 */
typedef struct
{
    const atmturb_screen_spec_t *specs;
    float *const                *out_a;
    float *const                *out_b;
    const int                   *vram_slot_a;
    const int                   *vram_slot_b;
    int                          npairs;
    int                          total_vram_screens;
    const float *const          *vram_host_tags;
} atmturb_screen_batch_req_t;

#ifdef HAVE_CUDA

/**
 * atmturb_cuda_generate_screen_pairs - Synthesize master screen pairs on CUDA GPU via cuFFT
 * @req: Batch request descriptor with screen specs, host buffers, and optional VRAM cache slots.
 *
 * Return: 0 on success, -1 on unsupported configuration or CUDA/cuFFT failure.
 */
int atmturb_cuda_generate_screen_pairs(
    const atmturb_screen_batch_req_t *req);

/**
 * atmturb_cuda_screens_cleanup - Free cached cuFFT plans and scratch buffers for screen synthesis
 */
void atmturb_cuda_screens_cleanup(void);

#else /* !HAVE_CUDA */

static inline int atmturb_cuda_generate_screen_pairs(
    const atmturb_screen_batch_req_t *req)
{
    (void) req;
    return -1;
}

static inline void atmturb_cuda_screens_cleanup(void)
{
}

#endif /* HAVE_CUDA */

#ifdef __cplusplus
}
#endif

#endif /* ATMTURB_CUDA_SCREENS_H */
