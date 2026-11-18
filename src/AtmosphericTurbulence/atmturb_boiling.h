// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_boiling.h
 * @brief   Co-moving Fourier AR(1)/AR(2) eddy boiling for multi-layer turbulence screens
 */

#ifndef ATMTURB_BOILING_H
#define ATMTURB_BOILING_H

#include <stdint.h>
#include <fftw3.h>

#include "atmturb_cuda_boiling.h"
#include "atmturb_geometry.h"
#include "atmturb_lowfreq.h"
#include "atmturb_profile.h"
#include "atmturb_rolling.h"

/**
 * struct atmturb_boil_layer_t - Per-layer co-moving Fourier AR(1)/AR(2) state
 * @active: Non-zero if this layer undergoes Fourier boiling (sigma_boil > 0).
 * @k_sub: Simulation frames per master-screen keyframe step (>= 1).
 * @vram_offset: Base slice index of this layer in the GPU master cache.
 * @max_kf_ready: Highest keyframe index materialized into screens[] (-1 before init).
 * @max_lf_ready: Highest frame index materialized into subharmonic ring (-1 before init).
 * @lf_cap: Capacity of per-frame subharmonic ring buffer (frames).
 * @rho_bar: PSD-weighted 1-step correlation between consecutive keyframes.
 * @k0_sq: Squared outer-scale wavenumber index (msize / L0_pix)^2.
 * @sqrt_pref_pair: Spectral amplitude factor for conjugate-paired complex modes.
 * @sqrt_pref_real: Spectral amplitude factor for self-conjugate Nyquist modes.
 * @inv_two_km2: Inner-scale damping coefficient 0.5 / km^2.
 * @alpha_dt: Dimensionless keyframe decay prefactor 2*pi*sigma_boil*dt_kf / (msize*dx).
 * @spec_cur: Persistent complex half-spectrum state (msize * (msize/2 + 1)).
 * @spec_aux: Auxiliary cascaded state u_m = psi_m - r*psi_{m-1} for AR(2) (or NULL).
 * @lf_amp: Subharmonic mode standard deviations (24 floats).
 * @lf_r: Subharmonic single-frame pole decay exp(-gamma_m * dt_frame) (24 floats).
 * @lf_omr2: Subharmonic 1 - r^2 factor (24 floats).
 * @lf_b: Subharmonic innovation gain factor (24 floats).
 * @lf_cur_re: Current real subharmonic state (24 floats).
 * @lf_cur_im: Current imaginary subharmonic state (24 floats).
 * @lf_aux_re: Auxiliary AR(2) real subharmonic state (24 floats).
 * @lf_aux_im: Auxiliary AR(2) imaginary subharmonic state (24 floats).
 * @lf_ring_re: Ring buffer of frame-rate real subharmonic amplitudes (lf_cap * 24 floats).
 * @lf_ring_im: Ring buffer of frame-rate imaginary subharmonic amplitudes (lf_cap * 24).
 */
typedef struct
{
    int            active;
    int            k_sub;
    int            vram_offset;
    long           max_kf_ready;
    long           max_lf_ready;
    long           lf_cap;
    double         rho_bar;
    float          k0_sq;
    float          sqrt_pref_pair;
    float          sqrt_pref_real;
    float          inv_two_km2;
    float          alpha_dt;
    fftwf_complex *spec_cur;
    fftwf_complex *spec_aux;
    float          lf_amp[ATMTURB_LOWFREQ_NMODES];
    float          lf_r[ATMTURB_LOWFREQ_NMODES];
    float          lf_omr2[ATMTURB_LOWFREQ_NMODES];
    float          lf_b[ATMTURB_LOWFREQ_NMODES];
    float          lf_cur_re[ATMTURB_LOWFREQ_NMODES];
    float          lf_cur_im[ATMTURB_LOWFREQ_NMODES];
    float          lf_aux_re[ATMTURB_LOWFREQ_NMODES];
    float          lf_aux_im[ATMTURB_LOWFREQ_NMODES];
    float         *lf_ring_re;
    float         *lf_ring_im;
} atmturb_boil_layer_t;

/**
 * struct atmturb_boiling_t - Top-level multi-layer Fourier boiling context
 * @ar_order: Autoregressive order (1 for AR(1), 2 for critically damped AR(2)).
 * @windowed: Non-zero when 3D cube keyframes exceed pool capacity and render in windows.
 * @msize: Master screen linear dimension [pixels].
 * @hcols: Half-spectrum column count (msize / 2 + 1).
 * @nlayers: Number of atmospheric layers.
 * @seed: Resolved base PRNG seed.
 * @c2r_plan: Cached single-precision 2D complex-to-real FFTW inverse plan.
 * @c2r_work: Scratch half-spectrum buffer preserved across destructive FFTW C2R calls.
 * @gpu: Optional GPU boiling acceleration context (NULL on CPU).
 * @layers: Array of per-layer AR states (nlayers elements).
 */
struct atmturb_boiling_t
{
    int                     ar_order;
    int                     windowed;
    long                    msize;
    long                    hcols;
    int                     nlayers;
    uint64_t                seed;
    fftwf_plan              c2r_plan;
    fftwf_complex          *c2r_work;
    atmturb_cuda_boiling_t *gpu;
    atmturb_boil_layer_t   *layers;
};

typedef struct atmturb_boiling_t atmturb_boiling_t;

int atmturb_boiling_init(
    atmturb_rolling_t       *r,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    long                     msize,
    long                     nbframes,
    double                   time_step_s,
    long                     precision,
    uint64_t                 seed);

void atmturb_boiling_free(
    atmturb_boiling_t *b);

int atmturb_boiling_advance(
    atmturb_rolling_t *r,
    long               target_t,
    int                gpu_only);

int atmturb_boiling_next_window(
    atmturb_rolling_t *r,
    long               t0,
    long               nbframes,
    long              *out_t1);

void atmturb_boiling_get_frame(
    const atmturb_rolling_t *r,
    int                      layer_idx,
    long                     t,
    atmturb_rolling_eval_t  *out);

typedef struct
{
    long     msize;
    double   r0_pix;
    double   L0_pix;
    double   l0_pix;
    double   gamma_scale;
    int      ar_order;
    uint64_t seed0;
    uint64_t seed1;
} atmturb_boil_lf_spec_t;

/* Internal leaf compute kernels (atmturb_boiling_kernels.c) */
double atmturb_boil_compute_rho1(
    const atmturb_boil_layer_t *bl,
    long                        msize,
    int                         ar_order);

void atmturb_boil_init_spectrum_cpu(
    atmturb_boil_layer_t *bl,
    fftwf_complex        *c2r_work,
    long                  msize,
    int                   ar_order,
    const uint64_t        seeds[2]);

void atmturb_boil_step_spectrum_cpu(
    atmturb_boil_layer_t *bl,
    fftwf_complex        *c2r_work,
    long                  msize,
    int                   ar_order,
    uint64_t              step_seed);

void atmturb_boil_init_lowfreq_layer(
    atmturb_boil_layer_t         *bl,
    const atmturb_lowfreq_t      *lf_base,
    const atmturb_boil_lf_spec_t *spec);

void atmturb_boil_step_lowfreq_layer(
    atmturb_boil_layer_t *bl,
    int                   ar_order,
    long                  frame_t,
    uint64_t              step_seed);

#endif /* ATMTURB_BOILING_H */
