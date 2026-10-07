// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_simd.h
 * @brief   SIMD-accelerated compute kernels for atmospheric turbulence simulation
 */

#ifndef ATMTURB_SIMD_H
#define ATMTURB_SIMD_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * enum atmturb_interp_mode_t - Sub-pixel interpolation schemes for extrusion
 * @ATMTURB_INTERP_BILINEAR: Legacy 4-point bilinear interpolation.
 * @ATMTURB_INTERP_BICUBIC: Separable Keys bicubic convolution kernel (a = -0.5).
 */
typedef enum
{
    ATMTURB_INTERP_BILINEAR = 0,
    ATMTURB_INTERP_BICUBIC  = 1
} atmturb_interp_mode_t;

/**
 * struct atmturb_extrude_params_t - Parameter bundle for phase screen extrusion
 * @master: Pointer to master phase screen data.
 * @msize: Master screen linear dimension in pixels.
 * @x0: Sub-pixel X coordinate offset in master pixels.
 * @y0: Sub-pixel Y coordinate offset in master pixels.
 * @pup_size: Linear dimension of the extracted pupil in pixels.
 * @os: Master grid oversampling factor (stride = os).
 * @interp: Interpolation scheme (0=bilinear, 1=Keys bicubic).
 * @weight: Layer Cn2 amplitude weight factor.
 * @out_pha: Output pupil phase array to accumulate into.
 */
typedef struct
{
    const float *master;
    long         msize;
    double       x0;
    double       y0;
    long         pup_size;
    long         os;
    int          interp;
    float        weight;
    float       *out_pha;
} atmturb_extrude_params_t;

/**
 * atmturb_extrude_accumulate - Phase screen extrusion with weighted accumulation
 * @params: Extrusion configuration and data pointers.
 */
void atmturb_extrude_accumulate(
    const atmturb_extrude_params_t *params);

/**
 * atmturb_scale_float_array - Multiply float array by scalar constant
 * @dest: Output float array.
 * @src: Input float array.
 * @scale: Multiplicative scale factor.
 * @n: Number of elements.
 */
void atmturb_scale_float_array(
    float       *dest,
    const float *src,
    float        scale,
    long         n);

/**
 * atmturb_init_phase_amp - Vectorized zero-init for phase and one-init for amplitude
 * @pha: Output phase array.
 * @amp: Output amplitude array.
 * @n: Number of elements.
 */
void atmturb_init_phase_amp(
    float *pha,
    float *amp,
    long   n);

/**
 * atmturb_simd_active_isa - Query name of active vectorized ISA implementation
 *
 * Return: String name of active ISA ("AVX-512", "AVX2", or "Scalar").
 */
const char *atmturb_simd_active_isa(void);

/* Scalar reference implementations */
void atmturb_extrude_accumulate_scalar(
    const atmturb_extrude_params_t *params);
void atmturb_extrude_accumulate_bilinear_scalar(
    const atmturb_extrude_params_t *params);
void atmturb_extrude_accumulate_bicubic_scalar(
    const atmturb_extrude_params_t *params);
void atmturb_scale_float_array_scalar(
    float       *dest,
    const float *src,
    float        scale,
    long         n);
void atmturb_init_phase_amp_scalar(
    float *pha,
    float *amp,
    long   n);

/* AVX2 implementations */
void atmturb_extrude_accumulate_avx2(
    const atmturb_extrude_params_t *params);
void atmturb_extrude_accumulate_bilinear_avx2(
    const atmturb_extrude_params_t *params);
void atmturb_extrude_accumulate_bicubic_avx2(
    const atmturb_extrude_params_t *params);
void atmturb_scale_float_array_avx2(
    float       *dest,
    const float *src,
    float        scale,
    long         n);
void atmturb_init_phase_amp_avx2(
    float *pha,
    float *amp,
    long   n);

/* AVX-512 implementations */
void atmturb_extrude_accumulate_avx512(
    const atmturb_extrude_params_t *params);
void atmturb_extrude_accumulate_bilinear_avx512(
    const atmturb_extrude_params_t *params);
void atmturb_extrude_accumulate_bicubic_avx512(
    const atmturb_extrude_params_t *params);
void atmturb_scale_float_array_avx512(
    float       *dest,
    const float *src,
    float        scale,
    long         n);
void atmturb_init_phase_amp_avx512(
    float *pha,
    float *amp,
    long   n);

#ifdef __cplusplus
}
#endif

#endif // ATMTURB_SIMD_H
