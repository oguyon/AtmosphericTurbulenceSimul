// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_rytov_simd_avx2.c
 * @brief   AVX2-accelerated transfer function kernels for Rytov propagation
 */

#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif

#include "atmturb_rytov_simd.h"

#if defined(__AVX2__)

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
    const float *spec_flt = (const float *) spec;
    float *dphi_flt       = (float *) acc_dphi;
    float *chi_flt        = (float *) acc_chi;

    static const int s_lo[8] = {0, 0, 1, 1, 2, 2, 3, 3};
    static const int s_hi[8] = {4, 4, 5, 5, 6, 6, 7, 7};
    __m256i idx_lo = _mm256_loadu_si256((const __m256i *) s_lo);
    __m256i idx_hi = _mm256_loadu_si256((const __m256i *) s_hi);

    long k = 0;
    for (; k <= ntot - 8; k += 8)
    {
        __m256 v_fa = _mm256_loadu_ps(&filt_a[k]);
        __m256 v_fb = _mm256_loadu_ps(&filt_b[k]);
        __m256 fa_lo = _mm256_permutevar8x32_ps(v_fa, idx_lo);
        __m256 fa_hi = _mm256_permutevar8x32_ps(v_fa, idx_hi);
        __m256 fb_lo = _mm256_permutevar8x32_ps(v_fb, idx_lo);
        __m256 fb_hi = _mm256_permutevar8x32_ps(v_fb, idx_hi);

        long idx = 2 * k;
        __m256 sp_lo = _mm256_loadu_ps(&spec_flt[idx]);
        __m256 sp_hi = _mm256_loadu_ps(&spec_flt[idx + 8]);

        __m256 dp_lo = _mm256_loadu_ps(&dphi_flt[idx]);
        __m256 dp_hi = _mm256_loadu_ps(&dphi_flt[idx + 8]);
        _mm256_storeu_ps(&dphi_flt[idx],     _mm256_fmadd_ps(sp_lo, fa_lo, dp_lo));
        _mm256_storeu_ps(&dphi_flt[idx + 8], _mm256_fmadd_ps(sp_hi, fa_hi, dp_hi));

        __m256 ch_lo = _mm256_loadu_ps(&chi_flt[idx]);
        __m256 ch_hi = _mm256_loadu_ps(&chi_flt[idx + 8]);
        _mm256_storeu_ps(&chi_flt[idx],     _mm256_fmadd_ps(sp_lo, fb_lo, ch_lo));
        _mm256_storeu_ps(&chi_flt[idx + 8], _mm256_fmadd_ps(sp_hi, fb_hi, ch_hi));
    }

    for (; k < ntot; k++)
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
    const float *spec_flt = (const float *) spec;
    const float *ramp_flt = (const float *) ramp;
    float *dphi_flt       = (float *) acc_dphi;
    float *chi_flt        = (float *) acc_chi;

    static const int s_lo[8] = {0, 0, 1, 1, 2, 2, 3, 3};
    static const int s_hi[8] = {4, 4, 5, 5, 6, 6, 7, 7};
    __m256i idx_lo = _mm256_loadu_si256((const __m256i *) s_lo);
    __m256i idx_hi = _mm256_loadu_si256((const __m256i *) s_hi);

    long k = 0;
    for (; k <= ntot - 8; k += 8)
    {
        __m256 v_fa = _mm256_loadu_ps(&filt_a[k]);
        __m256 v_fb = _mm256_loadu_ps(&filt_b[k]);
        __m256 fa_lo = _mm256_permutevar8x32_ps(v_fa, idx_lo);
        __m256 fa_hi = _mm256_permutevar8x32_ps(v_fa, idx_hi);
        __m256 fb_lo = _mm256_permutevar8x32_ps(v_fb, idx_lo);
        __m256 fb_hi = _mm256_permutevar8x32_ps(v_fb, idx_hi);

        long idx = 2 * k;
        __m256 sp_lo = _mm256_loadu_ps(&spec_flt[idx]);
        __m256 sp_hi = _mm256_loadu_ps(&spec_flt[idx + 8]);
        __m256 rm_lo = _mm256_loadu_ps(&ramp_flt[idx]);
        __m256 rm_hi = _mm256_loadu_ps(&ramp_flt[idx + 8]);

        __m256 sp_re0 = _mm256_moveldup_ps(sp_lo);
        __m256 sp_im0 = _mm256_movehdup_ps(sp_lo);
        __m256 rm_sw0 = _mm256_permute_ps(rm_lo, _MM_SHUFFLE(2, 3, 0, 1));
        __m256 ps_lo = _mm256_fmaddsub_ps(sp_re0, rm_lo, _mm256_mul_ps(sp_im0, rm_sw0));

        __m256 sp_re1 = _mm256_moveldup_ps(sp_hi);
        __m256 sp_im1 = _mm256_movehdup_ps(sp_hi);
        __m256 rm_sw1 = _mm256_permute_ps(rm_hi, _MM_SHUFFLE(2, 3, 0, 1));
        __m256 ps_hi = _mm256_fmaddsub_ps(sp_re1, rm_hi, _mm256_mul_ps(sp_im1, rm_sw1));

        __m256 dp_lo = _mm256_loadu_ps(&dphi_flt[idx]);
        __m256 dp_hi = _mm256_loadu_ps(&dphi_flt[idx + 8]);
        _mm256_storeu_ps(&dphi_flt[idx],     _mm256_fmadd_ps(ps_lo, fa_lo, dp_lo));
        _mm256_storeu_ps(&dphi_flt[idx + 8], _mm256_fmadd_ps(ps_hi, fa_hi, dp_hi));

        __m256 ch_lo = _mm256_loadu_ps(&chi_flt[idx]);
        __m256 ch_hi = _mm256_loadu_ps(&chi_flt[idx + 8]);
        _mm256_storeu_ps(&chi_flt[idx],     _mm256_fmadd_ps(ps_lo, fb_lo, ch_lo));
        _mm256_storeu_ps(&chi_flt[idx + 8], _mm256_fmadd_ps(ps_hi, fb_hi, ch_hi));
    }

    for (; k < ntot; k++)
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
    __m256 v_ey_comb = _mm256_setr_ps(ey_re, ey_im, ey_re, ey_im, ey_re, ey_im, ey_re, ey_im);
    __m256 v_ey_swap = _mm256_setr_ps(-ey_im, ey_re, -ey_im, ey_re, -ey_im, ey_re, -ey_im, ey_re);
    __m256 v_b_comb  = _mm256_setr_ps(b_re, b_im, b_re, b_im, b_re, b_im, b_re, b_im);
    __m256 v_b_swap  = _mm256_setr_ps(-b_im, b_re, -b_im, b_re, -b_im, b_re, -b_im, b_re);

    long i = 0;
    for (; i <= n_half - 4; i += 4)
    {
        __m256 va = _mm256_loadu_ps((const float *) &hat_a[i]);
        __m256 a_re = _mm256_moveldup_ps(va);
        __m256 a_im = _mm256_movehdup_ps(va);
        __m256 term1 = _mm256_fmadd_ps(a_im, v_ey_swap, _mm256_mul_ps(a_re, v_ey_comb));

        __m256 vex = _mm256_loadu_ps((const float *) &exp_x[i]);
        __m256 ex_re = _mm256_moveldup_ps(vex);
        __m256 ex_im = _mm256_movehdup_ps(vex);
        __m256 term2 = _mm256_fmadd_ps(ex_im, v_b_swap, _mm256_mul_ps(ex_re, v_b_comb));

        __m256 v_sum = _mm256_add_ps(term1, term2);

        __m128 l4 = _mm_loadu_ps(&laplace_inv_row[i]);
        __m256 l_dup = _mm256_castps128_ps256(_mm_unpacklo_ps(l4, l4));
        __m128 l_hi = _mm_unpackhi_ps(l4, l4);
        l_dup = _mm256_insertf128_ps(l_dup, l_hi, 1);

        __m256 sp = _mm256_loadu_ps((const float *) &spec_row[i]);
        _mm256_storeu_ps((float *) &spec_row[i], _mm256_fmadd_ps(v_sum, l_dup, sp));
    }

    for (; i < n_half; i++)
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

#else

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

#endif
