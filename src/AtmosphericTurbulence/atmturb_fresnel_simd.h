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

/**
 * atmturb_fresnel_simd_phasor_masked - Synthesize masked pupil complex field m * exp(i * scale * x)
 * @field: Output complex field array (npix elements, interleaved re/im).
 * @opd: Input phase or OPD array (NULL for unmodulated flat wavefront).
 * @mask: Real aperture transmission mask array.
 * @scale: Multiplicative scale converting @opd to radians.
 * @npix: Number of contiguous pixels to populate.
 */
void atmturb_fresnel_simd_phasor_masked(
    fftwf_complex *restrict field,
    const float   *restrict opd,
    const float   *restrict mask,
    float                   scale,
    long                    npix);

/**
 * atmturb_fresnel_simd_spectrum_strip - Synthesize 8 von Karman Fourier modes for r2 >= 256
 * @out: Output complex spectrum strip (8 elements).
 * @u1: Uniform (0, 1) Box-Muller radial samples (8 floats).
 * @u2: Uniform [0, 1) Box-Muller angular samples (8 floats).
 * @fx0: Starting X spatial frequency index (lanes use fx0 + 0..7).
 * @spec_p: 4-float parameter tuple {fy2, k0_sq, sqrt_pref, inv_two_km2}.
 */
void atmturb_fresnel_simd_spectrum_strip(
    fftwf_complex *restrict out,
    const float   *restrict u1,
    const float   *restrict u2,
    float                   fx0,
    const float   *restrict spec_p);

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

void atmturb_fresnel_phasor_masked_scalar(
    fftwf_complex *restrict field,
    const float   *restrict opd,
    const float   *restrict mask,
    float                   scale,
    long                    npix);

void atmturb_fresnel_spectrum_strip_scalar(
    fftwf_complex *restrict out,
    const float   *restrict u1,
    const float   *restrict u2,
    float                   fx0,
    const float   *restrict spec_p);

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

void atmturb_fresnel_phasor_masked_avx2(
    fftwf_complex *restrict field,
    const float   *restrict opd,
    const float   *restrict mask,
    float                   scale,
    long                    npix);

void atmturb_fresnel_spectrum_strip_avx2(
    fftwf_complex *restrict out,
    const float   *restrict u1,
    const float   *restrict u2,
    float                   fx0,
    const float   *restrict spec_p);

#ifdef __cplusplus
}
#endif

#endif // ATMTURB_FRESNEL_SIMD_H
