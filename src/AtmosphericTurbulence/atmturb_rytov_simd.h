// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_rytov_simd.h
 * @brief   SIMD-accelerated transfer function kernels for Rytov propagation
 */

#ifndef ATMTURB_RYTOV_SIMD_H
#define ATMTURB_RYTOV_SIMD_H

#include <fftw3.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * atmturb_rytov_simd_accumulate_filters - Multiply spectrum by filters and accumulate
 * @acc_dphi: Phase correction accumulator.
 * @acc_chi: Log-amplitude accumulator.
 * @spec: Periodic spectrum input.
 * @filt_a: Phase filter array.
 * @filt_b: Amplitude filter array.
 * @ntot: Elements in half-spectrum.
 */
void atmturb_rytov_simd_accumulate_filters(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          ntot);

/**
 * atmturb_rytov_simd_accumulate_filters_rotated - Apply chromatic ramp and accumulate
 * @acc_dphi: Phase correction accumulator.
 * @acc_chi: Log-amplitude accumulator.
 * @spec: Periodic spectrum input.
 * @ramp: Complex chromatic phase factor array.
 * @filt_a: Phase filter array.
 * @filt_b: Amplitude filter array.
 * @ntot: Elements in half-spectrum.
 */
void atmturb_rytov_simd_accumulate_filters_rotated(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const fftwf_complex *restrict ramp,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          ntot);

/**
 * atmturb_rytov_simd_decompose_row - Decompose boundary harmonics for one spectrum row
 * @spec_row: Target spectrum row complex array.
 * @hat_a: 1D Fourier transform of boundary A.
 * @exp_x: Precomputed complex exponential array along X.
 * @laplace_inv_row: Inverse discrete Laplacian row weights.
 * @ey_re: Real part of exp_y for current row.
 * @ey_im: Imaginary part of exp_y for current row.
 * @b_re: Real part of hat_b for current row.
 * @b_im: Imaginary part of hat_b for current row.
 * @n_half: Number of complex elements in half-row.
 */
void atmturb_rytov_simd_decompose_row(
    fftwf_complex *restrict       spec_row,
    const fftwf_complex *restrict hat_a,
    const fftwf_complex *restrict exp_x,
    const float *restrict         laplace_inv_row,
    float                         ey_re,
    float                         ey_im,
    float                         b_re,
    float                         b_im,
    long                          n_half);

/**
 * atmturb_rytov_simd_assemble_row - Add diffractive phase and evaluate Rytov amplitude row
 * @pha_row: Destination/in-out phase row (or NULL).
 * @amp_row: Destination amplitude row (or NULL).
 * @dphi_row: Source diffractive phase correction row.
 * @chi_row: Source log-amplitude row.
 * @len: Number of contiguous pixels in row span.
 *
 * Return: Sum of squared amplitudes across the row span.
 */
double atmturb_rytov_simd_assemble_row(
    float       *restrict pha_row,
    float       *restrict amp_row,
    const float *restrict dphi_row,
    const float *restrict chi_row,
    long                  len);

/* Scalar reference implementations */
void atmturb_rytov_accumulate_filters_scalar(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          ntot);

void atmturb_rytov_accumulate_filters_rotated_scalar(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const fftwf_complex *restrict ramp,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          ntot);

void atmturb_rytov_decompose_row_scalar(
    fftwf_complex *restrict       spec_row,
    const fftwf_complex *restrict hat_a,
    const fftwf_complex *restrict exp_x,
    const float *restrict         laplace_inv_row,
    float                         ey_re,
    float                         ey_im,
    float                         b_re,
    float                         b_im,
    long                          n_half);

double atmturb_rytov_assemble_row_scalar(
    float       *restrict pha_row,
    float       *restrict amp_row,
    const float *restrict dphi_row,
    const float *restrict chi_row,
    long                  len);

/* AVX2 implementations */
void atmturb_rytov_accumulate_filters_avx2(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          ntot);

void atmturb_rytov_accumulate_filters_rotated_avx2(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const fftwf_complex *restrict ramp,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          ntot);

void atmturb_rytov_decompose_row_avx2(
    fftwf_complex *restrict       spec_row,
    const fftwf_complex *restrict hat_a,
    const fftwf_complex *restrict exp_x,
    const float *restrict         laplace_inv_row,
    float                         ey_re,
    float                         ey_im,
    float                         b_re,
    float                         b_im,
    long                          n_half);

double atmturb_rytov_assemble_row_avx2(
    float       *restrict pha_row,
    float       *restrict amp_row,
    const float *restrict dphi_row,
    const float *restrict chi_row,
    long                  len);

/* AVX-512 implementations */
void atmturb_rytov_accumulate_filters_avx512(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          ntot);

void atmturb_rytov_accumulate_filters_rotated_avx512(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const fftwf_complex *restrict ramp,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          ntot);

void atmturb_rytov_decompose_row_avx512(
    fftwf_complex *restrict       spec_row,
    const fftwf_complex *restrict hat_a,
    const fftwf_complex *restrict exp_x,
    const float *restrict         laplace_inv_row,
    float                         ey_re,
    float                         ey_im,
    float                         b_re,
    float                         b_im,
    long                          n_half);

double atmturb_rytov_assemble_row_avx512(
    float       *restrict pha_row,
    float       *restrict amp_row,
    const float *restrict dphi_row,
    const float *restrict chi_row,
    long                  len);

#ifdef __cplusplus
}
#endif

#endif // ATMTURB_RYTOV_SIMD_H
