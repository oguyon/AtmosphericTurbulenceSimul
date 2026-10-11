// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_fresnel_kernels.h
 * @brief   GPU device kernel prototypes and launch wrappers for Fresnel propagation
 */

#ifndef ATMTURB_CUDA_FRESNEL_KERNELS_H
#define ATMTURB_CUDA_FRESNEL_KERNELS_H

#include <cuda_runtime.h>
#include <cufft.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * atmturb_cuda_fresnel_init_modulate_launch - Initialize complex field from initial phase screen
 * @d_field: Destination complex electric field array.
 * @d_super_pha: Source phase screen on padded grid.
 * @npix: Total number of grid pixels.
 * @stream: CUDA execution stream.
 *
 * Initializes unit-amplitude electric field E(x, y) = exp(i * phi(x, y)) on the simulation grid.
 */
void atmturb_cuda_fresnel_init_modulate_launch(
    cufftComplex *d_field,
    const float  *d_super_pha,
    int           npix,
    cudaStream_t  stream);

/**
 * atmturb_cuda_fresnel_modulate_launch - Modulate propagating complex field by superlayer phase
 * @d_field: In-place complex electric field array.
 * @d_super_pha: Superlayer phase screen to apply.
 * @npix: Total number of grid pixels.
 * @stream: CUDA execution stream.
 *
 * Applies thin phase screen modulation E(x, y) <- E(x, y) * exp(i * phi(x, y)).
 */
void atmturb_cuda_fresnel_modulate_launch(
    cufftComplex *d_field,
    const float  *d_super_pha,
    int           npix,
    cudaStream_t  stream);

/**
 * atmturb_cuda_fresnel_tf_mul_vec2_launch - Vectorized optical transfer function multiplication
 * @d_field: Complex spatial frequency field array (modified in place).
 * @d_tf: Complex Fresnel optical transfer function array.
 * @npix: Total number of complex elements (must be even).
 * @stream: CUDA execution stream.
 *
 * Multiplies Fourier-transformed optical field by Fresnel diffractive propagation kernel H(u, v)
 * using 128-bit vector float4 transactions.
 */
void atmturb_cuda_fresnel_tf_mul_vec2_launch(
    cufftComplex       *d_field,
    const cufftComplex *d_tf,
    int                 npix,
    cudaStream_t        stream);

/**
 * atmturb_cuda_fresnel_extract_launch - Extract amplitude, unwrapped phase, and piston sum
 * @d_field: Input propagated complex electric field.
 * @d_pha: Input geometric phase and output total unwrapped phase.
 * @d_amp: Destination amplitude array (or NULL if disabled).
 * @d_sum: Optional device memory pointer for parallel piston sum reduction (or NULL).
 * @npix: Total number of grid pixels.
 * @stream: CUDA execution stream.
 *
 * Extracts amplitude |E| and unwraps diffractive phase using geometric reference phase.
 */
void atmturb_cuda_fresnel_extract_launch(
    const cufftComplex *d_field,
    float              *d_pha,
    float              *d_amp,
    float              *d_sum,
    int                 npix,
    cudaStream_t        stream);

#ifdef __cplusplus
}
#endif

#endif /* ATMTURB_CUDA_FRESNEL_KERNELS_H */
