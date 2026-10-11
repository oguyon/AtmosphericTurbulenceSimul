// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_fresnel_simd.h
 * @brief   SIMD-accelerated field modulation and extraction kernels for Fresnel propagation
 */

#ifndef ATMTURB_FRESNEL_SIMD_H
#define ATMTURB_FRESNEL_SIMD_H

#include <fftw3.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * atmturb_fresnel_simd_modulate_field - Modulate complex optical field with phase screen
 * @field: In-out complex field array (npix elements, interleaved re/im).
 * @phase: Real phase screen array in radians.
 * @npix: Total number of pixels.
 */
void atmturb_fresnel_simd_modulate_field(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    npix);

/**
 * atmturb_fresnel_simd_init_modulate_field - Initialize field from unit plane wave and phase
 * @field: Output complex field array (npix elements, interleaved re/im).
 * @phase: Real phase screen array in radians.
 * @npix: Total number of pixels.
 */
void atmturb_fresnel_simd_init_modulate_field(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    npix);

/**
 * atmturb_fresnel_simd_extract_unwrapped - Extract amplitude and unwrapped phase
 * @field: Complex optical field at pupil.
 * @pha: Phase array holding accumulated geometric phase (updated in-place).
 * @amp: Destination amplitude array (can be NULL if not needed).
 * @npix: Total number of pixels.
 */
void atmturb_fresnel_simd_extract_unwrapped(
    const fftwf_complex *restrict field,
    float               *restrict pha,
    float               *restrict amp,
    long                          npix);

/* Scalar reference implementations */
void atmturb_fresnel_modulate_field_scalar(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    npix);

void atmturb_fresnel_init_modulate_field_scalar(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    npix);

void atmturb_fresnel_extract_unwrapped_scalar(
    const fftwf_complex *restrict field,
    float               *restrict pha,
    float               *restrict amp,
    long                          npix);

/* AVX2 implementations */
void atmturb_fresnel_modulate_field_avx2(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    npix);

void atmturb_fresnel_init_modulate_field_avx2(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    npix);

void atmturb_fresnel_extract_unwrapped_avx2(
    const fftwf_complex *restrict field,
    float               *restrict pha,
    float               *restrict amp,
    long                          npix);

#ifdef __cplusplus
}
#endif

#endif // ATMTURB_FRESNEL_SIMD_H
