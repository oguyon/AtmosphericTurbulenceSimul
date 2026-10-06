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
 * atmturb_extrude_accumulate - Bilinear phase screen extrusion with weighted accumulation
 * @master: Pointer to master phase screen data.
 * @msize: Master screen dimension (width = height = msize).
 * @x0: Sub-pixel X coordinate offset.
 * @y0: Sub-pixel Y coordinate offset.
 * @pup_size: Linear dimension of the extracted pupil.
 * @weight: Layer Cn2 amplitude weight factor.
 * @out_pha: Output pupil phase array to accumulate into.
 */
void atmturb_extrude_accumulate(const float *master, long msize, double x0, double y0,
                                long pup_size, float weight, float *out_pha);

/**
 * atmturb_scale_float_array - Multiply float array by scalar constant
 * @dest: Output float array.
 * @src: Input float array.
 * @scale: Multiplicative scale factor.
 * @n: Number of elements.
 */
void atmturb_scale_float_array(float *dest, const float *src, float scale, long n);

/**
 * atmturb_init_phase_amp - Vectorized zero-init for phase and one-init for amplitude
 * @pha: Output phase array.
 * @amp: Output amplitude array.
 * @n: Number of elements.
 */
void atmturb_init_phase_amp(float *pha, float *amp, long n);

#ifdef __cplusplus
}
#endif

#endif // ATMTURB_SIMD_H
