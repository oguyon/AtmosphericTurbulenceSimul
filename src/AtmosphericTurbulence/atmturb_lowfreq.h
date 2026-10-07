// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_lowfreq.h
 * @brief   Analytic Johansson-Gavel subharmonic low-order modes for turbulence screens
 */

#ifndef ATMTURB_LOWFREQ_H
#define ATMTURB_LOWFREQ_H

#include <stdint.h>

#define ATMTURB_LOWFREQ_NMODES 24

/**
 * struct atmturb_lowfreq_t - Precomputed subharmonic Fourier modes container
 * @nmodes: Number of active subharmonic modes (24 for 3 levels).
 * @kx: Spatial frequency along X in radians per master pixel.
 * @ky: Spatial frequency along Y in radians per master pixel.
 * @are: Real amplitude for screen 0 (p0).
 * @aim: Imaginary amplitude for screen 0 (p0).
 * @bre: Real amplitude for screen 1 (p1).
 * @bim: Imaginary amplitude for screen 1 (p1).
 */
typedef struct atmturb_lowfreq_t
{
    int   nmodes;
    float kx[ATMTURB_LOWFREQ_NMODES];
    float ky[ATMTURB_LOWFREQ_NMODES];
    float are[ATMTURB_LOWFREQ_NMODES];
    float aim[ATMTURB_LOWFREQ_NMODES];
    float bre[ATMTURB_LOWFREQ_NMODES];
    float bim[ATMTURB_LOWFREQ_NMODES];
} atmturb_lowfreq_t;

/**
 * atmturb_lowfreq_init - Initialize 24 subharmonic Fourier modes for a screen pair
 * @lf: Output low-order modes container.
 * @msize: Master screen size in pixels.
 * @r0_pix: Fried parameter in master pixels.
 * @L0_pix: Outer scale in master pixels (<=0 for infinite).
 * @l0_pix: Inner scale in master pixels (<=0 for none).
 * @seed: Random seed for mode amplitudes.
 */
void atmturb_lowfreq_init(
    atmturb_lowfreq_t *lf,
    long               msize,
    double             r0_pix,
    double             L0_pix,
    double             l0_pix,
    uint64_t           seed);

/**
 * atmturb_lowfreq_draw_amplitudes - Draw random complex Gaussian amplitudes for 24 modes
 * @are: Destination real amplitudes array (ATMTURB_LOWFREQ_NMODES floats).
 * @aim: Destination imaginary amplitudes array (ATMTURB_LOWFREQ_NMODES floats).
 * @msize: Master screen size in pixels.
 * @r0_pix: Fried parameter in master pixels.
 * @L0_pix: Outer scale in master pixels (<=0 for infinite).
 * @l0_pix: Inner scale in master pixels (<=0 for none).
 * @seed: Random seed for mode amplitudes.
 */
void atmturb_lowfreq_draw_amplitudes(
    float    *are,
    float    *aim,
    long      msize,
    double    r0_pix,
    double    L0_pix,
    double    l0_pix,
    uint64_t  seed);

/**
 * atmturb_lowfreq_accumulate - Analytically evaluate and accumulate low-order modes
 * @lf: Low-order modes container.
 * @screen_idx: 0 for screen A (p0), 1 for screen B (p1).
 * @x0: Continuous X position in master pixels.
 * @y0: Continuous Y position in master pixels.
 * @pup_size: Linear dimension of pupil.
 * @os: Oversampling factor (stride).
 * @weight: Layer scaling weight.
 * @out_pha: Destination phase array (pup_size * pup_size floats).
 */
void atmturb_lowfreq_accumulate(
    const atmturb_lowfreq_t *lf,
    int                      screen_idx,
    double                   x0,
    double                   y0,
    long                     pup_size,
    long                     os,
    float                    weight,
    float                   *out_pha);

/**
 * atmturb_lowfreq_accumulate_custom - Evaluate and accumulate modes with custom amplitudes
 * @lf: Low-order modes geometry container (kx, ky).
 * @are: Real amplitudes array (ATMTURB_LOWFREQ_NMODES floats).
 * @aim: Imaginary amplitudes array (ATMTURB_LOWFREQ_NMODES floats).
 * @x0: Continuous X position in master pixels.
 * @y0: Continuous Y position in master pixels.
 * @pup_size: Linear dimension of pupil.
 * @os: Oversampling factor (stride).
 * @weight: Layer scaling weight.
 * @out_pha: Destination phase array (pup_size * pup_size floats).
 */
void atmturb_lowfreq_accumulate_custom(
    const atmturb_lowfreq_t *lf,
    const float             *are,
    const float             *aim,
    double                   x0,
    double                   y0,
    long                     pup_size,
    long                     os,
    float                    weight,
    float                   *out_pha);

#endif /* ATMTURB_LOWFREQ_H */
