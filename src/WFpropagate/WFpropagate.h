// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    WFpropagate.h
 * @brief   Diffractive Fresnel optical propagation public interface
 */

#ifndef WFPROPAGATE_H
#define WFPROPAGATE_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * init_WFpropagate - Initialize WFpropagate module state
 *
 * Return: 0 on success.
 */
int init_WFpropagate(void);

/**
 * Fresnel_propagate_wavefront - Fresnel propagate complex optical field
 * @in: Name of input complex image.
 * @out: Name of output complex image.
 * @PUPIL_SCALE: Physical pixel scale [m/pixel].
 * @z: Propagation distance [m].
 * @lambda: Optical wavelength [m].
 *
 * Return: 0 on success.
 */
int Fresnel_propagate_wavefront(
    const char *in,
    const char *out,
    double      PUPIL_SCALE,
    double      z,
    double      lambda);

/**
 * Init_Fresnel_propagate_wavefront - Initialize anti-aliased Fresnel transfer kernel
 * @Cim: Output complex kernel image name.
 * @size: Grid dimension in pixels.
 * @PUPIL_SCALE: Physical pixel scale [m/pixel].
 * @z: Propagation distance [m].
 * @lambda: Optical wavelength [m].
 * @FPMASKRAD: Focal plane mask radius cutoff.
 * @Precision: 0 for single precision float, 1 for double precision.
 *
 * Return: 0 on success.
 */
int Init_Fresnel_propagate_wavefront(
    const char *Cim,
    long        size,
    double      PUPIL_SCALE,
    double      z,
    double      lambda,
    double      FPMASKRAD,
    int         Precision);

/**
 * Fresnel_propagate_wavefront1 - Apply precomputed transfer function kernel
 * @in: Name of input complex optical field.
 * @out: Name of output complex optical field.
 * @Cin: Name of precomputed complex transfer kernel.
 *
 * Return: 0 on success.
 */
int Fresnel_propagate_wavefront1(
    const char *in,
    const char *out,
    const char *Cin);

/**
 * Fresnel_propagate_cube - Propagate optical field through a series of distances
 * @IDcin_name: Input complex field name.
 * @IDout_name_amp: Output amplitude 3D cube image name.
 * @IDout_name_pha: Output phase 3D cube image name.
 * @PUPIL_SCALE: Physical pixel scale [m/pixel].
 * @zstart: Starting propagation distance [m].
 * @zend: Ending propagation distance [m].
 * @NBzpts: Number of distance sampling points.
 * @lambda: Optical wavelength [m].
 *
 * Return: 0 on success.
 */
long Fresnel_propagate_cube(
    const char *IDcin_name,
    const char *IDout_name_amp,
    const char *IDout_name_pha,
    double      PUPIL_SCALE,
    double      zstart,
    double      zend,
    long        NBzpts,
    double      lambda);

/**
 * WFpropagate_run - Standalone test harness for Lyot propagation
 *
 * Return: 0 on success.
 */
long WFpropagate_run(void);

#ifdef __cplusplus
}
#endif

#endif // WFPROPAGATE_H
