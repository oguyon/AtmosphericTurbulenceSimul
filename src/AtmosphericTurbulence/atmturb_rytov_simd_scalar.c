// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_rytov_simd_scalar.c
 * @brief   Scalar fallback implementations for Rytov diffractive kernels
 */

#include <math.h>
#include <stddef.h>

#include "atmturb_rytov_simd.h"

/**
 * atmturb_rytov_accumulate_filters_scalar - Scalar fallback for filter accumulation
 * @acc_dphi: Phase correction accumulator.
 * @acc_chi: Log-amplitude accumulator.
 * @spec: Periodic spectrum input.
 * @filt_a: Phase filter array.
 * @filt_b: Amplitude filter array.
 * @ntot: Elements in half-spectrum.
 */
void atmturb_rytov_accumulate_filters_scalar(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          ntot)
{
    for (long k = 0; k < ntot; k++)
    {
        float p_re = spec[k][0];
        float p_im = spec[k][1];
        float fa   = filt_a[k];
        float fb   = filt_b[k];

        acc_dphi[k][0] += p_re * fa;
        acc_dphi[k][1] += p_im * fa;
        acc_chi[k][0]  += p_re * fb;
        acc_chi[k][1]  += p_im * fb;
    }
}

/**
 * atmturb_rytov_accumulate_filters_rotated_scalar - Scalar fallback for rotated accumulation
 * @acc_dphi: Phase correction accumulator.
 * @acc_chi: Log-amplitude accumulator.
 * @spec: Periodic spectrum input.
 * @ramp: Complex chromatic phase factor array.
 * @filt_a: Phase filter array.
 * @filt_b: Amplitude filter array.
 * @ntot: Elements in half-spectrum.
 */
void atmturb_rytov_accumulate_filters_rotated_scalar(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const fftwf_complex *restrict ramp,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          ntot)
{
    for (long k = 0; k < ntot; k++)
    {
        float p_re = spec[k][0];
        float p_im = spec[k][1];
        float cb   = ramp[k][0];
        float sb   = ramp[k][1];

        float ps_re = p_re * cb - p_im * sb;
        float ps_im = p_re * sb + p_im * cb;

        float fa = filt_a[k];
        float fb = filt_b[k];

        acc_dphi[k][0] += ps_re * fa;
        acc_dphi[k][1] += ps_im * fa;
        acc_chi[k][0]  += ps_re * fb;
        acc_chi[k][1]  += ps_im * fb;
    }
}

/**
 * atmturb_rytov_decompose_row_scalar - Scalar fallback for Moisan harmonic decomposition
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
void atmturb_rytov_decompose_row_scalar(
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
    for (long i = 0; i < n_half; i++)
    {
        float a_re = hat_a[i][0];
        float a_im = hat_a[i][1];
        float ex_re = exp_x[i][0];
        float ex_im = exp_x[i][1];

        float v_re = (a_re * ey_re - a_im * ey_im) + (b_re * ex_re - b_im * ex_im);
        float v_im = (a_re * ey_im + a_im * ey_re) + (b_re * ex_im + b_im * ex_re);

        float linv = laplace_inv_row[i];
        spec_row[i][0] += v_re * linv;
        spec_row[i][1] += v_im * linv;
    }
}

/**
 * atmturb_rytov_assemble_row_scalar - Scalar reference row phase/amplitude assembly
 * @pha_row: Destination/in-out phase row (or NULL).
 * @amp_row: Destination amplitude row (or NULL).
 * @dphi_row: Source diffractive phase correction row.
 * @chi_row: Source log-amplitude row.
 * @len: Number of contiguous pixels in row span.
 *
 * Return: Sum of squared amplitudes across the row span.
 */
double atmturb_rytov_assemble_row_scalar(
    float       *restrict pha_row,
    float       *restrict amp_row,
    const float *restrict dphi_row,
    const float *restrict chi_row,
    long                  len)
{
    double sum_i = 0.0;

    if (pha_row != NULL)
    {
        for (long x = 0; x < len; x++)
        {
            pha_row[x] += dphi_row[x];
        }
    }
    if (amp_row != NULL)
    {
        for (long x = 0; x < len; x++)
        {
            float c = chi_row[x];
            if (c > 3.0f)
            {
                c = 3.0f;
            }
            else if (c < -3.0f)
            {
                c = -3.0f;
            }
            float a    = expf(c);
            amp_row[x] = a;
            sum_i += (double) (a * a);
        }
    }
    return sum_i;
}
