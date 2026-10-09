// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_rytov.h
 * @brief   Fourier-space Rytov diffractive propagation engine
 */

#ifndef ATMTURB_RYTOV_H
#define ATMTURB_RYTOV_H

#include <fftw3.h>

#include "atmturb_geometry.h"
#include "atmturb_profile.h"
#include "atmturb_rolling.h"
#include "atmturb_superlayer.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * struct atmturb_rytov_plan_t - Precomputed Rytov Fourier propagation plan
 * @nsuper: Number of super-layers.
 * @supers: Array of super-layer definitions and statistics.
 * @grid_size: Linear dimension of wavefront grid in pixels.
 * @pixscale_m: Grid physical pixel size [meters/pixel].
 * @lambda_ref_m: Primary reference wavelength [meters].
 * @lambda_s_m: Secondary observing wavelength [meters].
 * @sec_shared: 1 if secondary reuses primary spectrum with chromatic ramp (Option B).
 * @filter_a_pri: Array of nsuper primary phase diffractive filters [half-spectrum].
 * @filter_b_pri: Array of nsuper primary amplitude diffractive filters [half-spectrum].
 * @filter_a_sec: Array of nsuper secondary phase diffractive filters [half-spectrum].
 * @filter_b_sec: Array of nsuper secondary amplitude diffractive filters [half-spectrum].
 * @chrom_ramp: Precomputed complex chromatic phase factor arrays for Option B.
 * @laplace_inv: Precomputed inverse discrete Laplacian lookup table.
 * @exp_y: Precomputed vertical boundary Fourier factor array.
 * @exp_x: Precomputed horizontal boundary Fourier factor array.
 */
typedef struct
{
    int                   nsuper;
    atmturb_superlayer_t *supers;
    long                  grid_size;
    double                pixscale_m;
    double                lambda_ref_m;
    double                lambda_s_m;
    int                   sec_shared;
    float               **filter_a_pri;
    float               **filter_b_pri;
    float               **filter_a_sec;
    float               **filter_b_sec;
    fftwf_complex       **chrom_ramp;
    float                *laplace_inv;
    fftwf_complex        *exp_y;
    fftwf_complex        *exp_x;
} atmturb_rytov_plan_t;

/**
 * struct atmturb_rytov_ctx_t - Thread-local scratchpad and plans for Rytov render
 * @grid_size: Linear dimension of wavefront grid in pixels.
 * @plan_r2c: 2D real-to-complex FFTW plan.
 * @plan_c2r: 2D complex-to-real FFTW plan.
 * @plan_1d_a: 1D real-to-complex FFTW plan along horizontal boundary.
 * @plan_1d_b: 1D real-to-complex FFTW plan along vertical boundary.
 * @real_in: Real working array for 2D forward FFT.
 * @spec: Complex working half-spectrum.
 * @bound_a: Horizontal boundary difference vector.
 * @bound_b: Vertical boundary difference vector.
 * @hat_a: 1D half-spectrum of horizontal boundary.
 * @hat_b: 1D half-spectrum of vertical boundary.
 * @acc_dphi_pri: Primary phase correction accumulator.
 * @acc_chi_pri: Primary log-amplitude accumulator.
 * @acc_dphi_sec: Secondary phase correction accumulator.
 * @acc_chi_sec: Secondary log-amplitude accumulator.
 * @super_pha: Scratchpad for super-layer primary phase screen.
 * @super_spha: Scratchpad for super-layer secondary phase screen.
 * @dphi_out: Scratchpad for reconstructed phase correction.
 * @chi_out: Scratchpad for reconstructed log-amplitude.
 */
typedef struct
{
    long           grid_size;
    fftwf_plan     plan_r2c;
    fftwf_plan     plan_c2r;
    fftwf_plan     plan_1d_a;
    fftwf_plan     plan_1d_b;
    float         *real_in;
    fftwf_complex *spec;
    float         *bound_a;
    float         *bound_b;
    fftwf_complex *hat_a;
    fftwf_complex *hat_b;
    fftwf_complex *acc_dphi_pri;
    fftwf_complex *acc_chi_pri;
    fftwf_complex *acc_dphi_sec;
    fftwf_complex *acc_chi_sec;
    float         *super_pha;
    float         *super_spha;
    float         *dphi_out;
    float         *chi_out;
} atmturb_rytov_ctx_t;

/**
 * atmturb_rytov_plan_init - Initialize Rytov diffractive filters and lookup tables
 * @plan: Propagation plan container to populate.
 * @prof: Active turbulence profile.
 * @geom: Computed observing geometry.
 * @pup_size: Linear dimension of wavefront in pixels.
 * @pixscale_m: Physical pixel scale [m/pixel].
 * @lambda_ref_m: Primary reference wavelength [m].
 * @lambda_s_m: Secondary observing wavelength [m].
 * @z_bin_m: Altitude binning distance threshold [m].
 * @force_exact_sec: Force Option C (separate secondary extrusion and FFT).
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_rytov_plan_init(
    atmturb_rytov_plan_t    *plan,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    long                     pup_size,
    double                   pixscale_m,
    double                   lambda_ref_m,
    double                   lambda_s_m,
    double                   z_bin_m,
    int                      force_exact_sec);

/**
 * atmturb_rytov_plan_free - Release resources held by Rytov propagation plan
 * @plan: Propagation plan container to tear down.
 */
void atmturb_rytov_plan_free(
    atmturb_rytov_plan_t *plan);

/**
 * atmturb_rytov_ctx_init - Initialize thread-local Rytov render context and FFTW plans
 * @ctx: Thread context structure to initialize.
 * @pup_size: Linear dimension of wavefront in pixels.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_rytov_ctx_init(
    atmturb_rytov_ctx_t *ctx,
    long                 pup_size);

/**
 * atmturb_rytov_ctx_free - Release resources held by thread context
 * @ctx: Thread context structure to tear down.
 */
void atmturb_rytov_ctx_free(
    atmturb_rytov_ctx_t *ctx);

/**
 * atmturb_rytov_render_step - Render one frame using Fourier-space Rytov accumulation
 * @ctx: Thread-local execution context.
 * @plan: Precomputed Rytov propagation plan.
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @t: Simulation frame index.
 * @time_step_s: Time step between frames in seconds.
 * @master_size: Master screen dimension in pixels.
 * @pup_size: Linear dimension of wavefront in pixels.
 * @pha_slice: Destination primary phase slice.
 * @amp_slice: Destination primary amplitude slice.
 * @spha_slice: Destination secondary phase slice.
 * @samp_slice: Destination secondary amplitude slice.
 */
void atmturb_rytov_render_step(
    atmturb_rytov_ctx_t        *ctx,
    const atmturb_rytov_plan_t *plan,
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    long                        t,
    double                      time_step_s,
    long                        master_size,
    long                        pup_size,
    float                      *pha_slice,
    float                      *amp_slice,
    float                      *spha_slice,
    float                      *samp_slice);

#ifdef __cplusplus
}
#endif

#endif // ATMTURB_RYTOV_H
