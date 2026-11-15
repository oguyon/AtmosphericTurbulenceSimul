// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_rytov_simd_impl.h
 * @brief   Shared compiler-vectorized OpenMP SIMD kernels for Rytov propagation
 *
 * Included by target-specific translation units (`_avx2.c`, `_avx512.c`) compiled with
 * ISA flags so `#pragma omp simd` lowers to 256-bit or 512-bit FMA instructions without
 * explicit `<immintrin.h>` intrinsics.
 */

#ifndef ATMTURB_RYTOV_SIMD_IMPL_H
#define ATMTURB_RYTOV_SIMD_IMPL_H

#include <stddef.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "atmturb_rytov_simd.h"

/**
 * atmturb_rytov_fast_exp_clamped - Pure-FP squaring exp(clamp(x, -3, 3)) (~1.5 ULP)
 * @x: Input log-amplitude value.
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
 * atmturb_rytov_impl_accumulate_filters - Threaded + SIMD filter accumulation
 */
static inline void atmturb_rytov_impl_accumulate_filters(
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
 * atmturb_rytov_impl_accumulate_filters_rotated - Threaded + SIMD rotated accumulation
 */
static inline void atmturb_rytov_impl_accumulate_filters_rotated(
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
 * atmturb_rytov_impl_decompose_row - Vectorized Moisan harmonic row decomposition
 */
static inline void atmturb_rytov_impl_decompose_row(
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
 * atmturb_rytov_impl_assemble_row - Vectorized row phase and amplitude assembly
 */
static inline double atmturb_rytov_impl_assemble_row(
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

#endif // ATMTURB_RYTOV_SIMD_IMPL_H
