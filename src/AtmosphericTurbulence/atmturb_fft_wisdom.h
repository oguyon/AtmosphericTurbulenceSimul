// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_fft_wisdom.h
 * @brief   Thread-safe and process-atomic FFTW3 single-precision wisdom caching
 */

#ifndef ATMTURB_FFT_WISDOM_H
#define ATMTURB_FFT_WISDOM_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * atmturb_fftwf_plan_flags - Query FFTW planner flags for persistent 2D/1D transforms
 *
 * Defaults to FFTW_MEASURE unless ATMTURB_FFTW_MEASURE=0 is set in the environment.
 *
 * Return: FFTW planner bitmask (FFTW_MEASURE or FFTW_ESTIMATE).
 */
unsigned int atmturb_fftwf_plan_flags(void);

/**
 * atmturb_fftwf_wisdom_load - Import persistent single-precision FFTW wisdom from disk
 */
void atmturb_fftwf_wisdom_load(void);

/**
 * atmturb_fftwf_wisdom_save - Atomically export updated FFTW wisdom if new plans were measured
 */
void atmturb_fftwf_wisdom_save(void);

#ifdef __cplusplus
}
#endif

#endif // ATMTURB_FFT_WISDOM_H
