// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    wfprop_fresnel_engine.h
 * @brief   Thread-safe 2D Fresnel diffractive propagation engine
 */

#ifndef WFPROP_FRESNEL_ENGINE_H
#define WFPROP_FRESNEL_ENGINE_H

#include <fftw3.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * struct wfprop_fresnel_engine_t - Thread-local FFTW Fresnel propagation context
 * @n: Grid dimension in pixels (grid is n x n).
 * @buf: Pre-allocated complex working buffer (n * n elements).
 * @fwd_plan: Pre-computed FFTW forward 2D DFT plan.
 * @inv_plan: Pre-computed FFTW inverse 2D DFT plan.
 */
typedef struct
{
    long           n;
    fftwf_complex *buf;
    fftwf_plan     fwd_plan;
    fftwf_plan     inv_plan;
} wfprop_fresnel_engine_t;

/**
 * wfprop_fresnel_engine_init - Initialize thread-local Fresnel propagation engine
 * @eng: Pointer to engine structure to initialize.
 * @n: Grid dimension in pixels.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int wfprop_fresnel_engine_init(
    wfprop_fresnel_engine_t *eng,
    long                     n);

/**
 * wfprop_fresnel_engine_free - Release resources held by Fresnel propagation engine
 * @eng: Pointer to engine structure to tear down.
 */
void wfprop_fresnel_engine_free(
    wfprop_fresnel_engine_t *eng);

/**
 * wfprop_fresnel_tf_build - Precompute 2D angular spectrum transfer function
 * @tf: Output complex transfer function array (n * n elements).
 * @n: Grid dimension in pixels.
 * @pixscale_m: Grid physical sampling scale [m/pixel].
 * @z_m: Propagation distance along optical path [m].
 * @lambda_m: Optical wavelength [m].
 * @cutoff_rad: Frequency mask cutoff radius in pixels (<= 0 for Nyquist disc).
 *
 * Return: 0 on success.
 */
int wfprop_fresnel_tf_build(
    fftwf_complex *tf,
    long           n,
    double         pixscale_m,
    double         z_m,
    double         lambda_m,
    double         cutoff_rad);

/**
 * wfprop_fresnel_engine_apply - Apply optical diffraction step using precomputed transfer function
 * @eng: Initialized Fresnel engine.
 * @field: In-out 2D complex optical field (n * n elements).
 * @tf: Precomputed 2D complex transfer function (n * n elements).
 */
void wfprop_fresnel_engine_apply(
    wfprop_fresnel_engine_t *eng,
    fftwf_complex           *field,
    const fftwf_complex     *tf);

#ifdef __cplusplus
}
#endif

#endif // WFPROP_FRESNEL_ENGINE_H
