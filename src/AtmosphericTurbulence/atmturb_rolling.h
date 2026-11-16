// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_rolling.h
 * @brief   Multi-epoch rolling screen management and cross-fading for turbulence screens
 */

#ifndef ATMTURB_ROLLING_H
#define ATMTURB_ROLLING_H

#include <stdint.h>
#include "atmturb_geometry.h"
#include "atmturb_lowfreq.h"
#include "atmturb_profile.h"

/**
 * struct atmturb_rolling_screen_t - Single phase screen with low-frequency mode amplitudes
 * @data: Pointer to 2D master screen floating-point array (msize * msize).
 * @are: Real subharmonic amplitudes for this screen (ATMTURB_LOWFREQ_NMODES floats).
 * @aim: Imaginary subharmonic amplitudes for this screen (ATMTURB_LOWFREQ_NMODES floats).
 * @image_id: Image ID in milk dcimg table (-1 if not registered in dcimg).
 */
typedef struct
{
    float *data;
    float  are[ATMTURB_LOWFREQ_NMODES];
    float  aim[ATMTURB_LOWFREQ_NMODES];
    long   image_id;
} atmturb_rolling_screen_t;

/**
 * struct atmturb_rolling_layer_t - Per-layer rolling screens container
 * @nscreens: Number of active screens allocated for this layer.
 * @screens: Array of rolling screens for this layer.
 * @lf_base: Base subharmonic mode geometry (kx, ky, nmodes).
 * @t_dec_s: Epoch decorrelation time in seconds.
 */
typedef struct
{
    int                       nscreens;
    atmturb_rolling_screen_t *screens;
    atmturb_lowfreq_t         lf_base;
    double                    t_dec_s;
} atmturb_rolling_layer_t;

struct atmturb_boiling_t;

/**
 * struct atmturb_rolling_t - Top-level rolling screen context for all simulation layers
 * @nlayers: Number of simulation layers.
 * @msize: Master screen linear dimension.
 * @rolling: Flag indicating if rolling cross-fading is enabled (0 or 1).
 * @lowfreq: Flag indicating if analytic low-frequency modes are enabled (0 or 1).
 * @layers: Array of per-layer rolling context (nlayers elements).
 * @boil: Optional co-moving Fourier AR(1)/AR(2) boiling context (NULL when disabled).
 */
typedef struct
{
    int                       nlayers;
    long                      msize;
    int                       rolling;
    int                       lowfreq;
    atmturb_rolling_layer_t  *layers;
    struct atmturb_boiling_t *boil;
} atmturb_rolling_t;

/**
 * struct atmturb_rolling_eval_t - Frame evaluation state with screen pointers and blended modes
 * @scrA: Pointer to primary active screen array.
 * @scrB: Pointer to secondary active screen array (NULL if wB == 0).
 * @wA: Cross-fading weight for screen A (cos(theta)).
 * @wB: Cross-fading weight for screen B (sin(theta)).
 * @idxA: Index of primary active screen in layer screens array.
 * @idxB: Index of secondary active screen in layer screens array.
 * @are_eff: Blended lowfreq real amplitudes (ATMTURB_LOWFREQ_NMODES floats).
 * @aim_eff: Blended lowfreq imaginary amplitudes (ATMTURB_LOWFREQ_NMODES floats).
 */
typedef struct
{
    const float *scrA;
    const float *scrB;
    float        wA;
    float        wB;
    int          idxA;
    int          idxB;
    float        are_eff[ATMTURB_LOWFREQ_NMODES];
    float        aim_eff[ATMTURB_LOWFREQ_NMODES];
} atmturb_rolling_eval_t;

/**
 * atmturb_rolling_init - Allocate and generate rolling screens and low-order modes
 * @r: Pointer to rolling context to initialize.
 * @prof: Active atmospheric profile.
 * @geom: Computed observing geometry.
 * @msize: Master screen dimension in pixels.
 * @nbframes: Total number of frames in simulation.
 * @time_step_s: Time step per frame in seconds.
 * @precision: FFT precision (0=single, 1=double).
 * @seed: Master PRNG seed.
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_rolling_init(
    atmturb_rolling_t       *r,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    long                     msize,
    long                     nbframes,
    double                   time_step_s,
    long                     precision,
    uint64_t                 seed);

/**
 * atmturb_rolling_free - Release all resources associated with rolling context
 * @r: Pointer to rolling context.
 */
void atmturb_rolling_free(
    atmturb_rolling_t *r);

/**
 * atmturb_rolling_get_frame - Evaluate screen pointers, weights, and mode blending for frame t
 * @r: Rolling simulation context.
 * @layer_idx: Layer index (0..nlayers-1).
 * @t: Frame index (0..nbframes-1).
 * @time_step_s: Time step in seconds.
 * @out: Evaluated pointers and blended amplitudes.
 */
void atmturb_rolling_get_frame(
    const atmturb_rolling_t *r,
    int                      layer_idx,
    long                     t,
    double                   time_step_s,
    atmturb_rolling_eval_t  *out);

/**
 * atmturb_rolling_advance - Advance active Fourier boiling state up to frame target_t
 * @r: Rolling simulation context.
 * @target_t: Target frame index to ensure materialized.
 * @gpu_only: Non-zero to skip D2H screen copies when streaming exclusively on GPU.
 *
 * Return: 0 on success (or no-op when boiling disabled), -1 on failure.
 */
int atmturb_rolling_advance(
    atmturb_rolling_t *r,
    long               target_t,
    int                gpu_only);

/**
 * atmturb_rolling_next_window - Advance and compute next collision-free 3D cube frame window
 * @r: Rolling simulation context.
 * @t0: Inclusive window start frame index.
 * @nbframes: Total cube frame count.
 * @out_t1: Exclusive window end frame index.
 *
 * Return: 1 when [t0, *out_t1) is ready to render, 0 when t0 >= nbframes.
 */
int atmturb_rolling_next_window(
    atmturb_rolling_t *r,
    long               t0,
    long               nbframes,
    long              *out_t1);

/**
 * atmturb_rolling_is_windowed - Query whether 3D cube rendering requires multi-window stepping
 * @r: Rolling simulation context.
 *
 * Return: Non-zero if boiling pool is smaller than cube keyframes, 0 otherwise.
 */
int atmturb_rolling_is_windowed(
    const atmturb_rolling_t *r);

#endif /* ATMTURB_ROLLING_H */
