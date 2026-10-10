// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_rytov_internal.h
 * @brief   Internal compute kernels for Rytov diffractive propagation
 */

#ifndef ATMTURB_RYTOV_INTERNAL_H
#define ATMTURB_RYTOV_INTERNAL_H

#include "atmturb_rytov.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * atmturb_rytov_assemble_output - Add diffractive phase and normalize amplitude
 * @pup_size: Linear dimension of pupil.
 * @guard_pix: Guard band margin in pixels.
 * @pad_size: Linear dimension of padded compute grid.
 * @pha: Accumulated geometric phase (updated in-place with diffractive correction).
 * @amp: Destination amplitude array.
 * @dphi: Reconstructed diffractive phase correction.
 * @chi: Reconstructed log-amplitude array.
 */
void atmturb_rytov_assemble_output(
    long         pup_size,
    long         guard_pix,
    long         pad_size,
    float       *pha,
    float       *amp,
    const float *dphi,
    const float *chi);

/**
 * atmturb_rytov_decompose_periodic - Subtract harmonic smooth component from 2D spectrum
 * @ctx: Thread execution context.
 * @plan: Propagation plan.
 * @phi: Input 2D real phase screen.
 */
void atmturb_rytov_decompose_periodic(
    atmturb_rytov_ctx_t        *ctx,
    const atmturb_rytov_plan_t *plan,
    const float                *phi);

/**
 * atmturb_rytov_accumulate_filters - Multiply spectrum by filters and accumulate
 * @acc_dphi: Phase correction accumulator.
 * @acc_chi: Log-amplitude accumulator.
 * @spec: Periodic spectrum input.
 * @filt_a: Phase filter array.
 * @filt_b: Amplitude filter array.
 * @ntot: Elements in half-spectrum.
 */
void atmturb_rytov_accumulate_filters(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          ntot);

/**
 * atmturb_rytov_accumulate_filters_rotated - Apply chromatic ramp and accumulate filters
 * @acc_dphi: Phase correction accumulator.
 * @acc_chi: Log-amplitude accumulator.
 * @spec: Periodic spectrum input.
 * @ramp: Complex chromatic phase factor array.
 * @filt_a: Phase filter array.
 * @filt_b: Amplitude filter array.
 * @ntot: Elements in half-spectrum.
 */
void atmturb_rytov_accumulate_filters_rotated(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const fftwf_complex *restrict ramp,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          ntot);

/**
 * atmturb_rytov_init_tables - Build Poisson solver lookup tables
 * @plan: Propagation plan to populate.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_rytov_init_tables(
    atmturb_rytov_plan_t *plan);

/**
 * atmturb_rytov_build_layer_filters - Compute transfer function arrays for one super-layer
 * @plan: Propagation plan container.
 * @m: Super-layer index.
 * @npix: Total pixels in half-spectrum.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_rytov_build_layer_filters(
    atmturb_rytov_plan_t *plan,
    int                   m,
    size_t                npix);

/**
 * atmturb_rytov_wisdom_load - Import pre-computed FFTW wisdom from disk
 */
void atmturb_rytov_wisdom_load(void);

/**
 * atmturb_rytov_wisdom_save - Export updated FFTW wisdom to disk
 */
void atmturb_rytov_wisdom_save(void);

#ifdef __cplusplus
}
#endif

#endif // ATMTURB_RYTOV_INTERNAL_H
