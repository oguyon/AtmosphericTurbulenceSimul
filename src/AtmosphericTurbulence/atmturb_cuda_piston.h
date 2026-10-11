// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_piston.h
 * @brief   GPU-accelerated aperture piston reduction and subtraction
 *
 * Provides shared parallel reduction kernels using warp-shuffle primitives
 * and block-shared memory to compute and remove aperture-mean piston across
 * geometric, Fresnel, and Rytov simulation engines.
 */

#ifndef ATMTURB_CUDA_PISTON_H
#define ATMTURB_CUDA_PISTON_H

#include <cuda_runtime.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef HAVE_CUDA

static inline void atmturb_cuda_remove_piston_launch(
    float       *d_pha,
    float       *d_sum_buf,
    int          npix,
    cudaStream_t stream)
{
    (void) d_pha;
    (void) d_sum_buf;
    (void) npix;
    (void) stream;
}

static inline void atmturb_cuda_piston_sub_launch(
    float       *d_pha,
    const float *d_sum_buf,
    int          npix,
    cudaStream_t stream)
{
    (void) d_pha;
    (void) d_sum_buf;
    (void) npix;
    (void) stream;
}

#else

/**
 * atmturb_cuda_remove_piston_launch - Calculate and subtract aperture-mean piston on GPU
 * @d_pha: In-place phase screen array on device.
 * @d_sum_buf: Single-float working buffer on device for reduction sum.
 * @npix: Total number of active pupil pixels in phase screen.
 * @stream: CUDA execution stream.
 *
 * Computes the aperture sum using warp-shuffle reduction and atomic addition into
 * @d_sum_buf, followed by an in-place subtraction of the mean (sum / npix).
 */
void atmturb_cuda_remove_piston_launch(
    float       *d_pha,
    float       *d_sum_buf,
    int          npix,
    cudaStream_t stream);

/**
 * atmturb_cuda_piston_sub_launch - Subtract precomputed mean piston from phase screen on GPU
 * @d_pha: In-place phase screen array on device.
 * @d_sum_buf: Pointer to precomputed reduction sum of phase values on device.
 * @npix: Total number of active pupil pixels in phase screen.
 * @stream: CUDA execution stream.
 *
 * Used when a prior compute kernel (such as Fresnel extraction) has already accumulated
 * the aperture sum into @d_sum_buf during its output pass, eliminating the redundant
 * full-grid reduction pass.
 */
void atmturb_cuda_piston_sub_launch(
    float       *d_pha,
    const float *d_sum_buf,
    int          npix,
    cudaStream_t stream);

#endif // HAVE_CUDA

#ifdef __cplusplus
}
#endif

#endif // ATMTURB_CUDA_PISTON_H
