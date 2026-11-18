// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_rytov_simd_avx2.c
 * @brief   Compiler-vectorized AVX2 transfer-function and output kernels for Rytov propagation
 */

#include "atmturb_rytov_simd.h"

#if defined(__AVX2__)
#include "atmturb_rytov_simd_impl.h"

/**
 * atmturb_rytov_accumulate_filters_avx2 - AVX2 vectorized filter accumulation
 * @acc_dphi: Phase correction accumulator.
 * @acc_chi: Log-amplitude accumulator.
 * @spec: Periodic spectrum input.
 * @filt_a: Phase filter array.
 * @filt_b: Amplitude filter array.
 * @ntot: Elements in half-spectrum.
 */
void atmturb_rytov_accumulate_filters_avx2(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          ntot)
{
    atmturb_rytov_impl_accumulate_filters(acc_dphi, acc_chi, spec, filt_a, filt_b, ntot);
}

/**
 * atmturb_rytov_accumulate_filters_rotated_avx2 - AVX2 rotated filter accumulation
 * @acc_dphi: Phase correction accumulator.
 * @acc_chi: Log-amplitude accumulator.
 * @spec: Periodic spectrum input.
 * @ramp: Complex chromatic phase factor array.
 * @filt_a: Phase filter array.
 * @filt_b: Amplitude filter array.
 * @ntot: Elements in half-spectrum.
 */
void atmturb_rytov_accumulate_filters_rotated_avx2(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const fftwf_complex *restrict ramp,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          ntot)
{
    atmturb_rytov_impl_accumulate_filters_rotated(acc_dphi, acc_chi, spec, ramp,
                                                  filt_a, filt_b, ntot);
}

/**
 * atmturb_rytov_decompose_row_avx2 - AVX2 vectorized harmonic decomposition
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
void atmturb_rytov_decompose_row_avx2(
    fftwf_complex *restrict       spec_row,
    const fftwf_complex *restrict hat_a,
    const fftwf_complex *restrict exp_x,
    const float *restrict         laplace_inv_row,
    float                         ey_re,
    float                         ey_im,
    float                         b_re,
    float                         b_im,
    long                          n_half)
{
    atmturb_rytov_impl_decompose_row(spec_row, hat_a, exp_x, laplace_inv_row,
                                     ey_re, ey_im, b_re, b_im, n_half);
}

/**
 * atmturb_rytov_assemble_row_avx2 - AVX2 vectorized row phase/amplitude assembly
 * @pha_row: Destination/in-out phase row (or NULL).
 * @amp_row: Destination amplitude row (or NULL).
 * @dphi_row: Source diffractive phase correction row.
 * @chi_row: Source log-amplitude row.
 * @len: Number of contiguous pixels in row span.
 *
 * Return: Sum of squared amplitudes across the row span.
 */
double atmturb_rytov_assemble_row_avx2(
    float       *restrict pha_row,
    float       *restrict amp_row,
    const float *restrict dphi_row,
    const float *restrict chi_row,
    long                  len)
{
    return atmturb_rytov_impl_assemble_row(pha_row, amp_row, dphi_row, chi_row, len);
}

#else /* !defined(__AVX2__) */

void atmturb_rytov_accumulate_filters_avx2(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          ntot)
{
    atmturb_rytov_accumulate_filters_scalar(acc_dphi, acc_chi, spec, filt_a, filt_b, ntot);
}

void atmturb_rytov_accumulate_filters_rotated_avx2(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const fftwf_complex *restrict ramp,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          ntot)
{
    atmturb_rytov_accumulate_filters_rotated_scalar(acc_dphi, acc_chi, spec, ramp,
                                                    filt_a, filt_b, ntot);
}

void atmturb_rytov_decompose_row_avx2(
    fftwf_complex *restrict       spec_row,
    const fftwf_complex *restrict hat_a,
    const fftwf_complex *restrict exp_x,
    const float *restrict         laplace_inv_row,
    float                         ey_re,
    float                         ey_im,
    float                         b_re,
    float                         b_im,
    long                          n_half)
{
    atmturb_rytov_decompose_row_scalar(spec_row, hat_a, exp_x, laplace_inv_row,
                                       ey_re, ey_im, b_re, b_im, n_half);
}

double atmturb_rytov_assemble_row_avx2(
    float       *restrict pha_row,
    float       *restrict amp_row,
    const float *restrict dphi_row,
    const float *restrict chi_row,
    long                  len)
{
    return atmturb_rytov_assemble_row_scalar(pha_row, amp_row, dphi_row, chi_row, len);
}

#endif // defined(__AVX2__)
