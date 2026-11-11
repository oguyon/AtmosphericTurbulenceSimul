// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_rytov_simd_avx2.c
 * @brief   Compiler-vectorized AVX2 transfer-function and output kernels for Rytov propagation
 *
 * Compiled with `-mavx2 -mfma -fno-math-errno` so `#pragma omp simd` loops lower to 256-bit
 * AVX2 + FMA instructions without architecture-specific intrinsics.
 */

#include <stddef.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "atmturb_rytov_simd.h"

#if defined(__AVX2__)

/**
 * atmturb_rytov_fast_exp_clamped - Pure-FP squaring exp(clamp(x, -3, 3)) (~1.5 ULP)
 * @x: Input log-amplitude value.
 *
 * Uses scaling u = clamp(x, -3, 3) / 8 in [-0.375, 0.375], a degree-7 Horner polynomial, and
 * three squarings ((exp(u)^2)^2)^2 so the compiler lowers the entire evaluation to branchless
 * 256-bit AVX2 vmaxps/vminps/vmulps/vfmaddps instructions with zero integer conversions.
 *
 * Return: Clamped exponential value exp(clamp(x, -3, 3)).
 */
static inline float atmturb_rytov_fast_exp_clamped(
    float x)
{
    float c  = __builtin_fminf(3.0f, __builtin_fmaxf(-3.0f, x));
    float u  = c * 0.125f;
    float p  = 1.9841270e-4f;
    p = p * u + 1.3888889e-3f;
    p = p * u + 8.3333338e-3f;
    p = p * u + 4.1666668e-2f;
    p = p * u + 1.6666667e-1f;
    p = p * u + 5.0000000e-1f;
    p = p * u + 1.0f;
    p = p * u + 1.0f;

    float p2 = p * p;
    float p4 = p2 * p2;
    return p4 * p4;
}

/**
 * atmturb_rytov_acc_filt_span - Vectorized filter accumulation over [k0, k1)
 */
static inline void atmturb_rytov_acc_filt_span(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          k0,
    long                          k1)
{
    #pragma omp simd
    for (long k = k0; k < k1; k++)
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
#ifdef _OPENMP
    if (ntot >= 100000 && !omp_in_parallel())
    {
        int max_t = omp_get_max_threads();
        int nth   = (max_t > 8) ? 8 : max_t;
        #pragma omp parallel num_threads(nth)
        {
            int  tid = omp_get_thread_num();
            int  nt  = omp_get_num_threads();
            long k0  = (ntot * (long) tid) / (long) nt;
            long k1  = (ntot * (long) (tid + 1)) / (long) nt;
            atmturb_rytov_acc_filt_span(acc_dphi, acc_chi, spec, filt_a, filt_b, k0, k1);
        }
        return;
    }
#endif
    atmturb_rytov_acc_filt_span(acc_dphi, acc_chi, spec, filt_a, filt_b, 0, ntot);
}

/**
 * atmturb_rytov_acc_filt_rot_span - Vectorized rotated filter accumulation over [k0, k1)
 */
static inline void atmturb_rytov_acc_filt_rot_span(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const fftwf_complex *restrict ramp,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          k0,
    long                          k1)
{
    #pragma omp simd
    for (long k = k0; k < k1; k++)
    {
        float p_re  = spec[k][0];
        float p_im  = spec[k][1];
        float cb    = ramp[k][0];
        float sb    = ramp[k][1];
        float ps_re = p_re * cb - p_im * sb;
        float ps_im = p_re * sb + p_im * cb;
        float fa    = filt_a[k];
        float fb    = filt_b[k];

        acc_dphi[k][0] += ps_re * fa;
        acc_dphi[k][1] += ps_im * fa;
        acc_chi[k][0]  += ps_re * fb;
        acc_chi[k][1]  += ps_im * fb;
    }
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
#ifdef _OPENMP
    if (ntot >= 100000 && !omp_in_parallel())
    {
        int max_t = omp_get_max_threads();
        int nth   = (max_t > 8) ? 8 : max_t;
        #pragma omp parallel num_threads(nth)
        {
            int  tid = omp_get_thread_num();
            int  nt  = omp_get_num_threads();
            long k0  = (ntot * (long) tid) / (long) nt;
            long k1  = (ntot * (long) (tid + 1)) / (long) nt;
            atmturb_rytov_acc_filt_rot_span(acc_dphi, acc_chi, spec, ramp,
                                            filt_a, filt_b, k0, k1);
        }
        return;
    }
#endif
    atmturb_rytov_acc_filt_rot_span(acc_dphi, acc_chi, spec, ramp, filt_a, filt_b, 0, ntot);
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
    #pragma omp simd
    for (long i = 0; i < n_half; i++)
    {
        float a_re  = hat_a[i][0];
        float a_im  = hat_a[i][1];
        float ex_re = exp_x[i][0];
        float ex_im = exp_x[i][1];

        float v_re = (a_re * ey_re - a_im * ey_im) + (b_re * ex_re - b_im * ex_im);
        float v_im = (a_re * ey_im + a_im * ey_re) + (b_re * ex_im + b_im * ex_re);

        float linv      = laplace_inv_row[i];
        spec_row[i][0] += v_re * linv;
        spec_row[i][1] += v_im * linv;
    }
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
    if (pha_row != NULL)
    {
        #pragma omp simd
        for (long x = 0; x < len; x++)
        {
            pha_row[x] += dphi_row[x];
        }
    }
    if (amp_row == NULL)
    {
        return 0.0;
    }

    float sum_f = 0.0f;
    #pragma omp simd reduction(+ : sum_f)
    for (long x = 0; x < len; x++)
    {
        float a    = atmturb_rytov_fast_exp_clamped(chi_row[x]);
        amp_row[x] = a;
        sum_f += a * a;
    }
    return (double) sum_f;
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

#endif
