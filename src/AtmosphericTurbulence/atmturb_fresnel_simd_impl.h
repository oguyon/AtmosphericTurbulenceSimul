// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_fresnel_simd_impl.h
 * @brief   Shared compiler-vectorized OpenMP SIMD kernels for Fresnel propagation
 *
 * Included by target-specific translation units (`_avx2.c`, `_avx512.c`) compiled with
 * ISA flags so `#pragma omp simd` lowers branchless minimax polynomials to 256-bit or
 * 512-bit FMA instructions without architecture-specific `<immintrin.h>` intrinsics.
 */

#ifndef ATMTURB_FRESNEL_SIMD_IMPL_H
#define ATMTURB_FRESNEL_SIMD_IMPL_H

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "atmturb_fresnel_simd.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef union
{
    float    f;
    uint32_t u;
    int32_t  i;
} atmturb_fp32_bits_t;

/**
 * atmturb_fresnel_fast_sincos - Branchless single-precision sine and cosine (~1.2 ULP)
 * @v: Input angle in radians.
 * @s_out: Output sine value.
 * @c_out: Output cosine value.
 */
static inline void atmturb_fresnel_fast_sincos(
    float  v,
    float *s_out,
    float *c_out)
{
    float   fk = __builtin_rintf(v * (float) (2.0 / M_PI));
    int32_t k  = (int32_t) fk;

    float r = v - fk * 1.5703125f;
    r      -= fk * 0.0004837512969970703125f;
    r      -= fk * 7.549789948768648e-8f;

    float z      = r * r;
    float poly_s = -1.9515295891e-4f * z + 8.3321608736e-3f;
    poly_s       = poly_s * z - 1.6666654611e-1f;
    poly_s       = (poly_s * z) * r + r;

    float poly_c = -1.38537042643e-3f * z + 4.16635846761e-2f;
    poly_c       = poly_c * z - 4.99999053455e-1f;
    poly_c       = poly_c * z + 1.0f;

    float v_sin = ((k & 1) != 0) ? poly_c : poly_s;
    float v_cos = ((k & 1) != 0) ? poly_s : poly_c;

    *s_out = ((k & 2) != 0) ? -v_sin : v_sin;
    *c_out = (((k + 1) & 2) != 0) ? -v_cos : v_cos;
}

/**
 * atmturb_fresnel_fast_atan2 - Branchless minimax atan2f(y, x) (~2.4e-7 rad max error)
 * @y: Imaginary / numerator value.
 * @x: Real / denominator value.
 *
 * Return: Principal angle in [-pi, pi].
 */
static inline float atmturb_fresnel_fast_atan2(
    float y,
    float x)
{
    float ax = __builtin_fabsf(x);
    float ay = __builtin_fabsf(y);
    float mn = __builtin_fminf(ax, ay);
    float mx = __builtin_fmaxf(__builtin_fmaxf(ax, ay), 1.0e-30f);
    float a  = mn / mx;
    float s  = a * a;

    float r = 0.0028662257f;
    r       = r * s - 0.0161657367f;
    r       = r * s + 0.0429096138f;
    r       = r * s - 0.0752896400f;
    r       = r * s + 0.1065626393f;
    r       = r * s - 0.1420889944f;
    r       = r * s + 0.1999355085f;
    r       = r * s - 0.3333314528f;
    r       = (r * s) * a + a;

    r = (ay > ax) ? (1.5707963267948966f - r) : r;
    r = (x < 0.0f) ? (3.1415926535897932f - r) : r;
    return (y < 0.0f) ? -r : r;
}

/**
 * atmturb_fresnel_fast_log - Branchless natural logarithm for positive floats (~1.2 ULP)
 * @x: Positive normalized float input.
 *
 * Return: Natural logarithm ln(x).
 */
static inline float atmturb_fresnel_fast_log(
    float x)
{
    atmturb_fp32_bits_t bx = {.f = x};
    int32_t             ie = (int32_t) (bx.u >> 23) - 127;
    float               e  = (float) ie;

    atmturb_fp32_bits_t bm = {.u = (bx.u & 0x007fffffu) | 0x3f000000u};
    float               m  = bm.f;
    int                 lt = (m < 0.707106781186547524f);
    float               ad = lt ? m : 0.0f;
    m                      = (m - 1.0f) + ad;
    e                     += lt ? 0.0f : 1.0f;

    float z = m * m;
    float y = 7.0376836292e-2f;
    y       = y * m - 1.1514610310e-1f;
    y       = y * m + 1.1676998740e-1f;
    y       = y * m - 1.2420140846e-1f;
    y       = y * m + 1.4249322787e-1f;
    y       = y * m - 1.6668057665e-1f;
    y       = y * m + 2.0000714765e-1f;
    y       = y * m - 2.4999993993e-1f;
    y       = y * m + 3.3333331174e-1f;
    y       = y * (m * z);
    y       = e * (-2.12194440e-4f) + y;
    y       = y - 0.5f * z;
    y       = m + y;
    return e * 0.693359375f + y;
}

/**
 * atmturb_fresnel_fast_exp - Branchless natural exponential with clean underflow (~1.2 ULP)
 * @x: Input exponent.
 *
 * Return: Exponential exp(x).
 */
static inline float atmturb_fresnel_fast_exp(
    float x)
{
    const float lo_lim = -87.33654475f;
    const float hi_lim = 88.37626266f;
    int         uflow  = (x < lo_lim);
    float       xc     = __builtin_fminf(__builtin_fmaxf(x, lo_lim), hi_lim);

    float fx = __builtin_rintf(xc * 1.44269504088896341f);
    float r  = xc - fx * 0.693359375f;
    r       -= fx * (-2.12194440e-4f);

    float z = r * r;
    float y = 1.9875691500e-4f;
    y       = y * r + 1.3981999507e-3f;
    y       = y * r + 8.3334519073e-3f;
    y       = y * r + 4.1665795894e-2f;
    y       = y * r + 1.6666665459e-1f;
    y       = y * r + 5.0000001201e-1f;
    y       = y * z + (r + 1.0f);

    int32_t             ik = (int32_t) fx + 127;
    atmturb_fp32_bits_t bp = {.u = ((uint32_t) ik) << 23};
    return uflow ? 0.0f : (y * bp.f);
}

/**
 * atmturb_fresnel_modulate_span - Vectorized field modulation over [i0, i1)
 */
static inline void atmturb_fresnel_modulate_span(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    i0,
    long                    i1)
{
    #pragma omp simd
    for (long i = i0; i < i1; i++)
    {
        float s, c;
        atmturb_fresnel_fast_sincos(phase[i], &s, &c);
        float re    = field[i][0];
        float im    = field[i][1];
        field[i][0] = re * c - im * s;
        field[i][1] = re * s + im * c;
    }
}

/**
 * atmturb_fresnel_impl_modulate_field - Threaded + SIMD complex field modulation
 */
static inline void atmturb_fresnel_impl_modulate_field(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    npix)
{
#ifdef _OPENMP
    if (npix >= 262144 && !omp_in_parallel())
    {
        int max_t = omp_get_max_threads();
        int nth   = (max_t > 8) ? 8 : max_t;
        #pragma omp parallel num_threads(nth)
        {
            int  tid = omp_get_thread_num();
            int  nt  = omp_get_num_threads();
            long i0  = (npix * (long) tid) / (long) nt;
            long i1  = (npix * (long) (tid + 1)) / (long) nt;
            atmturb_fresnel_modulate_span(field, phase, i0, i1);
        }
        return;
    }
#endif
    atmturb_fresnel_modulate_span(field, phase, 0, npix);
}

/**
 * atmturb_fresnel_init_span - Vectorized unit-amplitude field initialization over [i0, i1)
 */
static inline void atmturb_fresnel_init_span(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    i0,
    long                    i1)
{
    #pragma omp simd
    for (long i = i0; i < i1; i++)
    {
        float s, c;
        atmturb_fresnel_fast_sincos(phase[i], &s, &c);
        field[i][0] = c;
        field[i][1] = s;
    }
}

/**
 * atmturb_fresnel_impl_init_modulate_field - Threaded + SIMD unit plane wave init
 */
static inline void atmturb_fresnel_impl_init_modulate_field(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    npix)
{
#ifdef _OPENMP
    if (npix >= 262144 && !omp_in_parallel())
    {
        int max_t = omp_get_max_threads();
        int nth   = (max_t > 8) ? 8 : max_t;
        #pragma omp parallel num_threads(nth)
        {
            int  tid = omp_get_thread_num();
            int  nt  = omp_get_num_threads();
            long i0  = (npix * (long) tid) / (long) nt;
            long i1  = (npix * (long) (tid + 1)) / (long) nt;
            atmturb_fresnel_init_span(field, phase, i0, i1);
        }
        return;
    }
#endif
    atmturb_fresnel_init_span(field, phase, 0, npix);
}

/**
 * atmturb_fresnel_extract_span - Vectorized unwrapped phase and amplitude span [i0, i1)
 */
static inline void atmturb_fresnel_extract_span(
    const fftwf_complex *restrict field,
    float               *restrict pha,
    float               *restrict amp,
    long                          i0,
    long                          i1)
{
    if (amp != NULL)
    {
        #pragma omp simd
        for (long i = i0; i < i1; i++)
        {
            float re      = field[i][0];
            float im      = field[i][1];
            amp[i]        = sqrtf(re * re + im * im);
            float phi_geo = pha[i];
            float s_geo, c_geo;
            atmturb_fresnel_fast_sincos(phi_geo, &s_geo, &c_geo);
            float d_re = re * c_geo + im * s_geo;
            float d_im = im * c_geo - re * s_geo;
            pha[i]     = phi_geo + atmturb_fresnel_fast_atan2(d_im, d_re);
        }
    }
    else
    {
        #pragma omp simd
        for (long i = i0; i < i1; i++)
        {
            float re      = field[i][0];
            float im      = field[i][1];
            float phi_geo = pha[i];
            float s_geo, c_geo;
            atmturb_fresnel_fast_sincos(phi_geo, &s_geo, &c_geo);
            float d_re = re * c_geo + im * s_geo;
            float d_im = im * c_geo - re * s_geo;
            pha[i]     = phi_geo + atmturb_fresnel_fast_atan2(d_im, d_re);
        }
    }
}

/**
 * atmturb_fresnel_impl_extract_unwrapped - Threaded + SIMD unwrapped extraction
 */
static inline void atmturb_fresnel_impl_extract_unwrapped(
    const fftwf_complex *restrict field,
    float               *restrict pha,
    float               *restrict amp,
    long                          npix)
{
#ifdef _OPENMP
    if (npix >= 262144 && !omp_in_parallel())
    {
        int max_t = omp_get_max_threads();
        int nth   = (max_t > 8) ? 8 : max_t;
        #pragma omp parallel num_threads(nth)
        {
            int  tid = omp_get_thread_num();
            int  nt  = omp_get_num_threads();
            long i0  = (npix * (long) tid) / (long) nt;
            long i1  = (npix * (long) (tid + 1)) / (long) nt;
            atmturb_fresnel_extract_span(field, pha, amp, i0, i1);
        }
        return;
    }
#endif
    atmturb_fresnel_extract_span(field, pha, amp, 0, npix);
}

/**
 * atmturb_fresnel_impl_phasor_masked - Vectorized masked pupil complex field synthesis
 */
static inline void atmturb_fresnel_impl_phasor_masked(
    fftwf_complex *restrict field,
    const float   *restrict opd,
    const float   *restrict mask,
    float                   scale,
    long                    npix)
{
    if (opd == NULL)
    {
        #pragma omp simd
        for (long i = 0; i < npix; i++)
        {
            field[i][0] = mask[i];
            field[i][1] = 0.0f;
        }
        return;
    }

    #pragma omp simd
    for (long i = 0; i < npix; i++)
    {
        float m = mask[i];
        float s, c;
        atmturb_fresnel_fast_sincos(opd[i] * scale, &s, &c);
        field[i][0] = m * c;
        field[i][1] = m * s;
    }
}

/**
 * atmturb_fresnel_impl_spectrum_strip - Vectorized 8-lane von Karman Fourier mode strip
 */
static inline void atmturb_fresnel_impl_spectrum_strip(
    fftwf_complex *restrict out,
    const float   *restrict u1,
    const float   *restrict u2,
    float                   fx0,
    const float   *restrict spec_p)
{
    const float fy2         = spec_p[0];
    const float k0_sq       = spec_p[1];
    const float sqrt_pref   = spec_p[2];
    const float inv_two_km2 = spec_p[3];
    const float two_pi      = 2.0f * (float) M_PI;

    #pragma omp simd
    for (int lane = 0; lane < 8; lane++)
    {
        float fx    = fx0 + (float) lane;
        float r2    = fx * fx + fy2;
        float dist2 = r2 + k0_sq;
        float r_bm  = sqrtf(-2.0f * atmturb_fresnel_fast_log(u1[lane]));
        float s, c;
        atmturb_fresnel_fast_sincos(two_pi * u2[lane], &s, &c);

        float exp_a = (-11.0f / 12.0f) * atmturb_fresnel_fast_log(dist2) - r2 * inv_two_km2;
        float mag   = sqrt_pref * atmturb_fresnel_fast_exp(exp_a) * r_bm;

        out[lane][0] = mag * c;
        out[lane][1] = mag * s;
    }
}

#endif // ATMTURB_FRESNEL_SIMD_IMPL_H
