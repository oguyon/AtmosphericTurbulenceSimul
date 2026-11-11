// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_simd_avx512.c
 * @brief   AVX-512 target instantiation of compiler-vectorized turbulence compute kernels
 *
 * Compiled with `-mavx512f -mavx512dq -mavx512bw -mavx512vl -mfma -fno-math-errno` so the
 * compiler lowers shared `#pragma omp simd` kernels in `atmturb_simd_extrude_impl.h` and
 * `atmturb_simd_ops_impl.h` to 512-bit AVX-512 vector instructions with identical OpenMP
 * threading and toroidal seam-splitting behavior as the AVX2 backend.
 */

#define _GNU_SOURCE
#include "atmturb_simd.h"

#if defined(__AVX512F__)

#include "atmturb_simd_extrude_impl.h"
#include "atmturb_simd_ops_impl.h"

/**
 * atmturb_extrude_accumulate_bilinear_avx512 - AVX-512 accelerated bilinear extrusion
 * @params: Extrusion configuration and data pointers.
 */
void atmturb_extrude_accumulate_bilinear_avx512(
    const atmturb_extrude_params_t *params)
{
    atmturb_extrude_bilinear_impl(params);
}

/**
 * atmturb_extrude_accumulate_bicubic_avx512 - AVX-512 accelerated Keys bicubic extrusion
 * @params: Extrusion configuration and data pointers.
 */
void atmturb_extrude_accumulate_bicubic_avx512(
    const atmturb_extrude_params_t *params)
{
    atmturb_extrude_bicubic_impl(params);
}

/**
 * atmturb_extrude_accumulate_avx512 - AVX-512 extrusion dispatcher (bilinear/bicubic)
 * @params: Extrusion configuration and data pointers.
 */
void atmturb_extrude_accumulate_avx512(
    const atmturb_extrude_params_t *params)
{
    if (params->interp == ATMTURB_INTERP_BICUBIC)
    {
        atmturb_extrude_bicubic_impl(params);
    }
    else
    {
        atmturb_extrude_bilinear_impl(params);
    }
}

/**
 * atmturb_scale_float_array_avx512 - Vectorized array multiplication by scalar (AVX-512)
 * @dest: Output float array.
 * @src: Input float array.
 * @scale: Scalar multiplier.
 * @n: Number of elements.
 */
void atmturb_scale_float_array_avx512(
    float       *dest,
    const float *src,
    float        scale,
    long         n)
{
    atmturb_scale_float_array_impl(dest, src, scale, n);
}

/**
 * atmturb_init_phase_amp_avx512 - Vectorized phase zeroing and amplitude one-filling (AVX-512)
 * @pha: Phase array.
 * @amp: Amplitude array.
 * @n: Number of elements.
 */
void atmturb_init_phase_amp_avx512(
    float *pha,
    float *amp,
    long   n)
{
    atmturb_init_phase_amp_impl(pha, amp, n);
}

/**
 * atmturb_extrude_lowfreq_avx512 - AVX-512 L1-tiled low-order subharmonic mode accumulation
 * @params: Low-frequency configuration and data pointers.
 */
void atmturb_extrude_lowfreq_avx512(
    const atmturb_lowfreq_params_t *params)
{
    atmturb_extrude_lowfreq_impl(params);
}

/**
 * atmturb_add_float_array_avx512 - Vectorized array accumulation (dest[i] += src[i])
 * @dest: Output/accumulator float array.
 * @src: Input float array.
 * @n: Number of elements.
 */
void atmturb_add_float_array_avx512(
    float       *dest,
    const float *src,
    long         n)
{
    atmturb_add_float_array_impl(dest, src, n);
}

/**
 * atmturb_complex_mul_array_avx512 - Vectorized complex array product (AVX-512)
 * @dest: Output complex float array (length 2 * n_complex).
 * @src1: First input complex float array.
 * @src2: Second input complex float array.
 * @n_complex: Number of complex elements.
 */
void atmturb_complex_mul_array_avx512(
    float       *dest,
    const float *src1,
    const float *src2,
    long         n_complex)
{
    atmturb_complex_mul_array_impl(dest, src1, src2, n_complex);
}

/**
 * atmturb_remove_piston_stream_avx512 - AVX-512 vectorized piston removal with stream copy
 * @dst: Destination phase array (may alias src).
 * @src: Source phase array.
 * @npix: Total number of pixels.
 */
void atmturb_remove_piston_stream_avx512(
    float       *dst,
    const float *src,
    long         npix)
{
    atmturb_remove_piston_stream_impl(dst, src, npix);
}

#else /* !defined(__AVX512F__) */

void atmturb_extrude_accumulate_avx512(
    const atmturb_extrude_params_t *params)
{
    atmturb_extrude_accumulate_scalar(params);
}

void atmturb_extrude_accumulate_bilinear_avx512(
    const atmturb_extrude_params_t *params)
{
    atmturb_extrude_accumulate_bilinear_scalar(params);
}

void atmturb_extrude_accumulate_bicubic_avx512(
    const atmturb_extrude_params_t *params)
{
    atmturb_extrude_accumulate_bicubic_scalar(params);
}

void atmturb_scale_float_array_avx512(
    float       *dest,
    const float *src,
    float        scale,
    long         n)
{
    atmturb_scale_float_array_scalar(dest, src, scale, n);
}

void atmturb_init_phase_amp_avx512(
    float *pha,
    float *amp,
    long   n)
{
    atmturb_init_phase_amp_scalar(pha, amp, n);
}

void atmturb_extrude_lowfreq_avx512(
    const atmturb_lowfreq_params_t *params)
{
    atmturb_extrude_lowfreq_scalar(params);
}

void atmturb_add_float_array_avx512(
    float       *dest,
    const float *src,
    long         n)
{
    atmturb_add_float_array_scalar(dest, src, n);
}

void atmturb_complex_mul_array_avx512(
    float       *dest,
    const float *src1,
    const float *src2,
    long         n_complex)
{
    atmturb_complex_mul_array_scalar(dest, src1, src2, n_complex);
}

void atmturb_remove_piston_stream_avx512(
    float       *dst,
    const float *src,
    long         npix)
{
    atmturb_remove_piston_stream_scalar(dst, src, npix);
}

#endif /* defined(__AVX512F__) */
