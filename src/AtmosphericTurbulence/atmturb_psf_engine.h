// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_psf_engine.h
 * @brief   Pruned 2-stage 1D FFTW3 Fourier optics engine for Nyquist science PSF synthesis
 */

#ifndef ATMTURB_PSF_ENGINE_H
#define ATMTURB_PSF_ENGINE_H

#include <fftw3.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * struct atmturb_ao_psf_engine_t - Pruned 2-stage 1D Fourier optics engine for Nyquist PSF
 * @pup_size: Pupil linear dimension in pixels.
 * @n_fft: Zero-padded 1D FFT dimension (power of 2, >= 2 * pup_size).
 * @n_psf: Output cropped PSF dimension (min(n_fft, 512)).
 * @col_stride: Cache-padded stride (`n_fft + 8`) preventing power-of-2 cache-set aliasing.
 * @norm_factor: Parseval intensity normalization factor (`1 / n_fft^2`).
 * @row_in: Zero-padded horizontal input rows (`pup_size * n_fft`).
 * @row_out: Horizontal 1D DFT output rows (`pup_size * n_fft`).
 * @col_in: Pruned and FFT-shifted vertical input columns (`n_psf * col_stride`).
 * @col_out: Vertical 1D DFT focal output (`n_psf * col_stride`).
 * @plan_rows: Cached FFTW batched 1D plan for `pup_size` active horizontal rows.
 * @plan_cols: Cached FFTW batched 1D plan for `n_psf` cropped vertical columns.
 */
typedef struct
{
    long           pup_size;
    long           n_fft;
    long           n_psf;
    long           col_stride;
    double         norm_factor;
    fftwf_complex *row_in;
    fftwf_complex *row_out;
    fftwf_complex *col_in;
    fftwf_complex *col_out;
    fftwf_plan     plan_rows;
    fftwf_plan     plan_cols;
    fftwf_plan     plan_2d;
} atmturb_ao_psf_engine_t;

/**
 * atmturb_ao_psf_alloc_bufs - Allocate zero-padded pruned 2-stage FFT working arrays
 * @eng: PSF engine descriptor with dimensions populated.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_ao_psf_alloc_bufs(
    atmturb_ao_psf_engine_t *eng);

/**
 * atmturb_ao_psf_free_bufs - Release complex working arrays of PSF engine
 * @eng: PSF engine descriptor.
 */
void atmturb_ao_psf_free_bufs(
    atmturb_ao_psf_engine_t *eng);

/**
 * atmturb_ao_psf_engine_init - Allocate pruned 2-stage buffers and create cached FFTW plans
 * @eng: PSF engine structure to initialize.
 * @pup_size: Pupil dimension in pixels.
 * @n_fft: Zero-padded FFT dimension (power of 2, >= 2 * pup_size).
 * @fftw_nth: Number of FFTW worker threads per transform.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_ao_psf_engine_init(
    atmturb_ao_psf_engine_t *eng,
    long                     pup_size,
    long                     n_fft,
    int                      fftw_nth);

/**
 * atmturb_ao_psf_engine_free - Release FFTW plans and complex arrays
 * @eng: Initialized PSF engine.
 */
void atmturb_ao_psf_engine_free(
    atmturb_ao_psf_engine_t *eng);

/**
 * atmturb_ao_render_psf - Form instantaneous zero-padded centered PSF via pruned 2-stage FFT
 * @eng: Initialized PSF engine (zero-padded borders pre-cleared).
 * @phase_in: Input OPD or phase array [pup_size * pup_size], or NULL for flat wavefront.
 * @pupil_mask: Circular pupil mask [pup_size * pup_size].
 * @pup_size: Pupil grid linear dimension.
 * @phase_scale: Multiplicative factor converting @phase_in to science radians.
 * @psf_buf: Destination or accumulation 2D PSF buffer [n_psf * n_psf].
 * @accum: Non-zero to add into @psf_buf, 0 to overwrite @psf_buf.
 *
 * Return: Peak pixel intensity of the instantaneous PSF slice.
 */
double atmturb_ao_render_psf(
    const atmturb_ao_psf_engine_t *eng,
    const float                   *phase_in,
    const float                   *pupil_mask,
    long                           pup_size,
    float                          phase_scale,
    float                         *psf_buf,
    int                            accum);

/**
 * atmturb_ao_render_psf_amp - Form centered PSF with optional scintillation amplitude map
 * @eng: Initialized PSF engine (`eng->pup_size` defines pupil dimension).
 * @phase_in: Input OPD or phase array [pup_size * pup_size], or NULL for flat wavefront.
 * @amp_in: Optional pupil electric-field amplitude [pup_size * pup_size], or NULL for 1.0.
 * @pupil_mask: Circular pupil mask [pup_size * pup_size].
 * @phase_scale: Multiplicative factor converting @phase_in to science radians.
 * @psf_buf: Destination or accumulation 2D PSF buffer [n_psf * n_psf].
 * @accum: Non-zero to add into @psf_buf, 0 to overwrite @psf_buf.
 *
 * Return: Peak pixel intensity of the instantaneous PSF slice.
 */
double atmturb_ao_render_psf_amp(
    const atmturb_ao_psf_engine_t *eng,
    const float                   *phase_in,
    const float                   *amp_in,
    const float                   *pupil_mask,
    float                          phase_scale,
    float                         *psf_buf,
    int                            accum);

#ifdef __cplusplus
}
#endif

#endif /* ATMTURB_PSF_ENGINE_H */
