// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_simd_avx2_ops.c
 * @brief   AVX2 target instantiation of compiler-vectorized array, lowfreq, and piston kernels
 *
 * Compiled with `-mavx2 -mfma -fno-math-errno` so the compiler lowers `#pragma omp simd`
 * loops in `atmturb_simd_ops_impl.h` to 256-bit AVX2 + FMA instructions, including `vfmaddsubps`
 * complex products, 4-way ILP piston reductions, and L1-tiled 24-mode subharmonic synthesis.
 */

#define _GNU_SOURCE
#include "atmturb_simd.h"

#if defined(__AVX2__)

#include "atmturb_simd_ops_impl.h"

/**
 * atmturb_scale_float_array_avx2 - Vectorized array multiplication by scalar (AVX2)
 * @dest: Output float array.
 * @src: Input float array.
 * @scale: Scalar multiplier.
 * @n: Number of elements.
 */
void atmturb_scale_float_array_avx2(
    float       *dest,
    const float *src,
    float        scale,
    long         n)
{
    atmturb_scale_float_array_impl(dest, src, scale, n);
}

/**
 * atmturb_init_phase_amp_avx2 - Vectorized phase zeroing and amplitude one-filling (AVX2)
 * @pha: Phase array to zero (or NULL to skip).
 * @amp: Amplitude array to fill with 1.0f (or NULL to skip).
 * @n: Number of elements.
 */
void atmturb_init_phase_amp_avx2(
    float *pha,
    float *amp,
    long   n)
{
    atmturb_init_phase_amp_impl(pha, amp, n);
}

/**
 * atmturb_add_float_array_avx2 - Vectorized array accumulation (dest[i] += src[i])
 * @dest: Output/accumulator float array.
 * @src: Input float array.
 * @n: Number of elements.
 */
void atmturb_add_float_array_avx2(
    float       *dest,
    const float *src,
    long         n)
{
    atmturb_add_float_array_impl(dest, src, n);
}

/**
 * atmturb_complex_mul_array_avx2 - Vectorized complex array point-wise product (AVX2)
 * @dest: Output interleaved complex float array (length 2 * n_complex).
 * @src1: First input interleaved complex float array.
 * @src2: Second input interleaved complex float array.
 * @n_complex: Number of complex elements.
 */
void atmturb_complex_mul_array_avx2(
    float       *dest,
    const float *src1,
    const float *src2,
    long         n_complex)
{
    atmturb_complex_mul_array_impl(dest, src1, src2, n_complex);
}

/**
 * atmturb_extrude_lowfreq_avx2 - AVX2 L1-tiled low-order subharmonic mode accumulation
 * @params: Low-frequency configuration and data pointers.
 */
void atmturb_extrude_lowfreq_avx2(
    const atmturb_lowfreq_params_t *params)
{
    atmturb_extrude_lowfreq_impl(params);
}

/**
 * atmturb_remove_piston_stream_avx2 - AVX2 vectorized piston removal with stream copy
 * @dst: Destination phase array (may alias src).
 * @src: Source phase array.
 * @npix: Total number of pixels.
 */
void atmturb_remove_piston_stream_avx2(
    float       *dst,
    const float *src,
    long         npix)
{
    atmturb_remove_piston_stream_impl(dst, src, npix);
}

#else /* !defined(__AVX2__) */

void atmturb_scale_float_array_avx2(
    float       *dest,
    const float *src,
    float        scale,
    long         n)
{
    atmturb_scale_float_array_scalar(dest, src, scale, n);
}

void atmturb_init_phase_amp_avx2(
    float *pha,
    float *amp,
    long   n)
{
    atmturb_init_phase_amp_scalar(pha, amp, n);
}

void atmturb_extrude_lowfreq_avx2(
    const atmturb_lowfreq_params_t *params)
{
    atmturb_extrude_lowfreq_scalar(params);
}

void atmturb_add_float_array_avx2(
    float       *dest,
    const float *src,
    long         n)
{
    atmturb_add_float_array_scalar(dest, src, n);
}

void atmturb_complex_mul_array_avx2(
    float       *dest,
    const float *src1,
    const float *src2,
    long         n_complex)
{
    atmturb_complex_mul_array_scalar(dest, src1, src2, n_complex);
}

void atmturb_remove_piston_stream_avx2(
    float       *dst,
    const float *src,
    long         npix)
{
    atmturb_remove_piston_stream_scalar(dst, src, npix);
}

#endif /* defined(__AVX2__) */
