// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_fresnel_simd_avx2.c
 * @brief   Compiler-vectorized AVX2 field modulation and extraction kernels for Fresnel propagation
 */

#include "atmturb_fresnel_simd.h"

#if defined(__AVX2__)
#include "atmturb_fresnel_simd_impl.h"

/**
 * atmturb_fresnel_modulate_field_avx2 - AVX2 accelerated optical field modulation
 * @field: In-out complex field array.
 * @phase: Real phase screen array in radians.
 * @npix: Total number of pixels.
 */
void atmturb_fresnel_modulate_field_avx2(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    npix)
{
    atmturb_fresnel_impl_modulate_field(field, phase, npix);
}

/**
 * atmturb_fresnel_init_modulate_field_avx2 - AVX2 accelerated unit plane wave init
 * @field: Output complex field array.
 * @phase: Real phase screen array in radians.
 * @npix: Total number of pixels.
 */
void atmturb_fresnel_init_modulate_field_avx2(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    npix)
{
    atmturb_fresnel_impl_init_modulate_field(field, phase, npix);
}

/**
 * atmturb_fresnel_extract_unwrapped_avx2 - AVX2 accelerated unwrapped extraction
 * @field: Complex optical field at pupil.
 * @pha: Phase array holding accumulated geometric phase (updated in-place).
 * @amp: Destination amplitude array (can be NULL).
 * @npix: Total number of pixels.
 */
void atmturb_fresnel_extract_unwrapped_avx2(
    const fftwf_complex *restrict field,
    float               *restrict pha,
    float               *restrict amp,
    long                          npix)
{
    atmturb_fresnel_impl_extract_unwrapped(field, pha, amp, npix);
}

/**
 * atmturb_fresnel_phasor_masked_avx2 - AVX2 masked pupil complex field synthesis
 * @field: Output complex field array.
 * @opd: Input phase or OPD array (NULL for flat wavefront).
 * @mask: Real aperture mask array.
 * @scale: Multiplicative scale converting @opd to radians.
 * @npix: Number of contiguous pixels.
 */
void atmturb_fresnel_phasor_masked_avx2(
    fftwf_complex *restrict field,
    const float   *restrict opd,
    const float   *restrict mask,
    float                   scale,
    long                    npix)
{
    atmturb_fresnel_impl_phasor_masked(field, opd, mask, scale, npix);
}

/**
 * atmturb_fresnel_spectrum_strip_avx2 - AVX2 8-lane von Karman Fourier mode synthesis
 * @out: Output complex spectrum strip (8 elements).
 * @u1: Uniform (0, 1) Box-Muller radial samples (8 floats).
 * @u2: Uniform [0, 1) Box-Muller angular samples (8 floats).
 * @fx0: Starting X spatial frequency index.
 * @spec_p: 4-float parameter tuple {fy2, k0_sq, sqrt_pref, inv_two_km2}.
 */
void atmturb_fresnel_spectrum_strip_avx2(
    fftwf_complex *restrict out,
    const float   *restrict u1,
    const float   *restrict u2,
    float                   fx0,
    const float   *restrict spec_p)
{
    atmturb_fresnel_impl_spectrum_strip(out, u1, u2, fx0, spec_p);
}

#else /* !defined(__AVX2__) */

void atmturb_fresnel_modulate_field_avx2(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    npix)
{
    atmturb_fresnel_modulate_field_scalar(field, phase, npix);
}

void atmturb_fresnel_init_modulate_field_avx2(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    npix)
{
    atmturb_fresnel_init_modulate_field_scalar(field, phase, npix);
}

void atmturb_fresnel_extract_unwrapped_avx2(
    const fftwf_complex *restrict field,
    float               *restrict pha,
    float               *restrict amp,
    long                          npix)
{
    atmturb_fresnel_extract_unwrapped_scalar(field, pha, amp, npix);
}

void atmturb_fresnel_phasor_masked_avx2(
    fftwf_complex *restrict field,
    const float   *restrict opd,
    const float   *restrict mask,
    float                   scale,
    long                    npix)
{
    atmturb_fresnel_phasor_masked_scalar(field, opd, mask, scale, npix);
}

void atmturb_fresnel_spectrum_strip_avx2(
    fftwf_complex *restrict out,
    const float   *restrict u1,
    const float   *restrict u2,
    float                   fx0,
    const float   *restrict spec_p)
{
    atmturb_fresnel_spectrum_strip_scalar(out, u1, u2, fx0, spec_p);
}

#endif // defined(__AVX2__)
