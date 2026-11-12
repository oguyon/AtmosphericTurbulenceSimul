// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_fresnel_simd_avx2.c
 * @brief   AVX2-accelerated field modulation and extraction kernels for Fresnel propagation
 */

#define _GNU_SOURCE
#include <math.h>
#include <stdint.h>

#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "atmturb_fresnel_simd.h"

#if defined(__AVX2__)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * sincos256_ps - Vectorized single-precision sine and cosine on 8 floats
 * @v: Input radians vector.
 * @s: Output sine vector.
 * @c: Output cosine vector.
 */
static inline void sincos256_ps(
    __m256  v,
    __m256 *s,
    __m256 *c)
{
    const __m256 inv_pio2 = _mm256_set1_ps((float) (2.0 / M_PI));
    const __m256 p1       = _mm256_set1_ps(1.5703125f);
    const __m256 p2       = _mm256_set1_ps(0.0004837512969970703125f);
    const __m256 p3       = _mm256_set1_ps(7.549789948768648e-8f);

    __m256 v_k = _mm256_round_ps(_mm256_mul_ps(v, inv_pio2),
                                 _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    __m256i vi_k = _mm256_cvtps_epi32(v_k);

    __m256 r = _mm256_fnmadd_ps(v_k, p1, v);
    r        = _mm256_fnmadd_ps(v_k, p2, r);
    r        = _mm256_fnmadd_ps(v_k, p3, r);

    __m256 z = _mm256_mul_ps(r, r);

    const __m256 s3     = _mm256_set1_ps(-1.9515295891e-4f);
    const __m256 s2     = _mm256_set1_ps(+8.3321608736e-3f);
    const __m256 s1     = _mm256_set1_ps(-1.6666654611e-1f);
    __m256       poly_s = _mm256_fmadd_ps(z, s3, s2);
    poly_s              = _mm256_fmadd_ps(z, poly_s, s1);
    poly_s              = _mm256_fmadd_ps(_mm256_mul_ps(poly_s, z), r, r);

    const __m256 c3     = _mm256_set1_ps(-1.38537042643e-3f);
    const __m256 c2     = _mm256_set1_ps(+4.16635846761e-2f);
    const __m256 c1     = _mm256_set1_ps(-4.99999053455e-1f);
    __m256       poly_c = _mm256_fmadd_ps(z, c3, c2);
    poly_c              = _mm256_fmadd_ps(z, poly_c, c1);
    poly_c              = _mm256_fmadd_ps(poly_c, z, _mm256_set1_ps(1.0f));

    __m256i swap_mask = _mm256_slli_epi32(_mm256_and_si256(vi_k, _mm256_set1_epi32(1)), 31);
    __m256  v_sin     = _mm256_blendv_ps(poly_s, poly_c, _mm256_castsi256_ps(swap_mask));
    __m256  v_cos     = _mm256_blendv_ps(poly_c, poly_s, _mm256_castsi256_ps(swap_mask));

    __m256i sign_s = _mm256_slli_epi32(_mm256_and_si256(vi_k, _mm256_set1_epi32(2)), 30);
    __m256i sign_c = _mm256_slli_epi32(_mm256_and_si256(
        _mm256_add_epi32(vi_k, _mm256_set1_epi32(1)), _mm256_set1_epi32(2)), 30);

    *s = _mm256_xor_ps(v_sin, _mm256_castsi256_ps(sign_s));
    *c = _mm256_xor_ps(v_cos, _mm256_castsi256_ps(sign_c));
}

/**
 * atan2_256_ps - Vectorized 8-lane minimax atan2f(y, x) (~2.4e-7 rad max error)
 * @y: Imaginary / numerator vector.
 * @x: Real / denominator vector.
 *
 * Return: 8-lane angle vector in [-pi, pi].
 */
static inline __m256 atan2_256_ps(
    __m256 y,
    __m256 x)
{
    const __m256 sign_mask = _mm256_set1_ps(-0.0f);
    const __m256 abs_mask  = _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff));
    const __m256 pio2      = _mm256_set1_ps(1.5707963267948966f);
    const __m256 pi        = _mm256_set1_ps(3.1415926535897932f);
    const __m256 tiny      = _mm256_set1_ps(1.0e-30f);

    __m256 ax = _mm256_and_ps(x, abs_mask);
    __m256 ay = _mm256_and_ps(y, abs_mask);
    __m256 mn = _mm256_min_ps(ax, ay);
    __m256 mx = _mm256_max_ps(_mm256_max_ps(ax, ay), tiny);
    __m256 a  = _mm256_div_ps(mn, mx);
    __m256 s  = _mm256_mul_ps(a, a);

    __m256 r = _mm256_set1_ps(0.0028662257f);
    r        = _mm256_fmadd_ps(r, s, _mm256_set1_ps(-0.0161657367f));
    r        = _mm256_fmadd_ps(r, s, _mm256_set1_ps(0.0429096138f));
    r        = _mm256_fmadd_ps(r, s, _mm256_set1_ps(-0.0752896400f));
    r        = _mm256_fmadd_ps(r, s, _mm256_set1_ps(0.1065626393f));
    r        = _mm256_fmadd_ps(r, s, _mm256_set1_ps(-0.1420889944f));
    r        = _mm256_fmadd_ps(r, s, _mm256_set1_ps(0.1999355085f));
    r        = _mm256_fmadd_ps(r, s, _mm256_set1_ps(-0.3333314528f));
    r        = _mm256_fmadd_ps(_mm256_mul_ps(r, s), a, a);

    __m256 swap  = _mm256_cmp_ps(ay, ax, _CMP_GT_OQ);
    r            = _mm256_blendv_ps(r, _mm256_sub_ps(pio2, r), swap);
    __m256 x_neg = _mm256_cmp_ps(x, _mm256_setzero_ps(), _CMP_LT_OQ);
    r            = _mm256_blendv_ps(r, _mm256_sub_ps(pi, r), x_neg);

    return _mm256_xor_ps(r, _mm256_and_ps(y, sign_mask));
}

/**
 * atmturb_fresnel_modulate_blocks - Modulate 8-pixel blocks [b0, b1) of complex field
 */
static inline void atmturb_fresnel_modulate_blocks(
    float       *restrict f_flt,
    const float *restrict phase,
    long                  b0,
    long                  b1)
{
    static const int s_idx_lo[8] = {0, 0, 1, 1, 2, 2, 3, 3};
    static const int s_idx_hi[8] = {4, 4, 5, 5, 6, 6, 7, 7};
    __m256i          idx_lo      = _mm256_loadu_si256((const __m256i *) s_idx_lo);
    __m256i          idx_hi      = _mm256_loadu_si256((const __m256i *) s_idx_hi);

    for (long b = b0; b < b1; b++)
    {
        long   i  = b * 8;
        __m256 vp = _mm256_loadu_ps(&phase[i]);
        __m256 vs, vc;
        sincos256_ps(vp, &vs, &vc);

        __m256 c_lo  = _mm256_permutevar8x32_ps(vc, idx_lo);
        __m256 s_lo  = _mm256_permutevar8x32_ps(vs, idx_lo);
        __m256 c_hi  = _mm256_permutevar8x32_ps(vc, idx_hi);
        __m256 s_hi  = _mm256_permutevar8x32_ps(vs, idx_hi);
        __m256 cs_lo = _mm256_blend_ps(c_lo, s_lo, 0xAA);
        __m256 sc_lo = _mm256_blend_ps(s_lo, c_lo, 0xAA);
        __m256 cs_hi = _mm256_blend_ps(c_hi, s_hi, 0xAA);
        __m256 sc_hi = _mm256_blend_ps(s_hi, c_hi, 0xAA);

        long   idx = 2 * i;
        __m256 f0  = _mm256_loadu_ps(&f_flt[idx]);
        __m256 f1  = _mm256_loadu_ps(&f_flt[idx + 8]);
        __m256 re0 = _mm256_moveldup_ps(f0);
        __m256 im0 = _mm256_movehdup_ps(f0);
        __m256 re1 = _mm256_moveldup_ps(f1);
        __m256 im1 = _mm256_movehdup_ps(f1);

        _mm256_storeu_ps(&f_flt[idx], _mm256_fmaddsub_ps(re0, cs_lo, _mm256_mul_ps(im0, sc_lo)));
        _mm256_storeu_ps(&f_flt[idx + 8],
                         _mm256_fmaddsub_ps(re1, cs_hi, _mm256_mul_ps(im1, sc_hi)));
    }
}

/**
 * atmturb_fresnel_modulate_field_avx2 - AVX2 accelerated optical field modulation
 * @field: In-out complex field array.
 * @phase: Real phase screen array in radians.
 * @npix: Total number of pixels.
 */
void atmturb_fresnel_modulate_field_avx2(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    npix)
{
    float *f_flt   = (float *) field;
    long   nblocks = npix / 8;

#ifdef _OPENMP
    if (npix >= 262144 && !omp_in_parallel())
    {
        int max_t = omp_get_max_threads();
        int nth   = (max_t > 8) ? 8 : max_t;
        #pragma omp parallel num_threads(nth)
        {
            int  tid = omp_get_thread_num();
            int  nt  = omp_get_num_threads();
            long b0  = (nblocks * (long) tid) / (long) nt;
            long b1  = (nblocks * (long) (tid + 1)) / (long) nt;
            atmturb_fresnel_modulate_blocks(f_flt, phase, b0, b1);
        }
    }
    else
#endif
    {
        atmturb_fresnel_modulate_blocks(f_flt, phase, 0, nblocks);
    }

    for (long i = nblocks * 8; i < npix; i++)
    {
        float s, c;
        sincosf(phase[i], &s, &c);
        float re    = field[i][0];
        float im    = field[i][1];
        field[i][0] = re * c - im * s;
        field[i][1] = re * s + im * c;
    }
}

/**
 * atmturb_fresnel_init_blocks - Initialize 8-pixel blocks [b0, b1) of unit-amplitude field
 */
static inline void atmturb_fresnel_init_blocks(
    float       *restrict f_flt,
    const float *restrict phase,
    long                  b0,
    long                  b1)
{
    static const int s_idx_lo[8] = {0, 0, 1, 1, 2, 2, 3, 3};
    static const int s_idx_hi[8] = {4, 4, 5, 5, 6, 6, 7, 7};
    __m256i          idx_lo      = _mm256_loadu_si256((const __m256i *) s_idx_lo);
    __m256i          idx_hi      = _mm256_loadu_si256((const __m256i *) s_idx_hi);

    for (long b = b0; b < b1; b++)
    {
        long   i  = b * 8;
        __m256 vp = _mm256_loadu_ps(&phase[i]);
        __m256 vs, vc;
        sincos256_ps(vp, &vs, &vc);

        __m256 c_lo = _mm256_permutevar8x32_ps(vc, idx_lo);
        __m256 s_lo = _mm256_permutevar8x32_ps(vs, idx_lo);
        __m256 c_hi = _mm256_permutevar8x32_ps(vc, idx_hi);
        __m256 s_hi = _mm256_permutevar8x32_ps(vs, idx_hi);

        long idx = 2 * i;
        _mm256_storeu_ps(&f_flt[idx],     _mm256_blend_ps(c_lo, s_lo, 0xAA));
        _mm256_storeu_ps(&f_flt[idx + 8], _mm256_blend_ps(c_hi, s_hi, 0xAA));
    }
}

/**
 * atmturb_fresnel_init_modulate_field_avx2 - AVX2 accelerated unit plane wave init
 * @field: Output complex field array.
 * @phase: Real phase screen array in radians.
 * @npix: Total number of pixels.
 */
void atmturb_fresnel_init_modulate_field_avx2(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    npix)
{
    float *f_flt   = (float *) field;
    long   nblocks = npix / 8;

#ifdef _OPENMP
    if (npix >= 262144 && !omp_in_parallel())
    {
        int max_t = omp_get_max_threads();
        int nth   = (max_t > 8) ? 8 : max_t;
        #pragma omp parallel num_threads(nth)
        {
            int  tid = omp_get_thread_num();
            int  nt  = omp_get_num_threads();
            long b0  = (nblocks * (long) tid) / (long) nt;
            long b1  = (nblocks * (long) (tid + 1)) / (long) nt;
            atmturb_fresnel_init_blocks(f_flt, phase, b0, b1);
        }
    }
    else
#endif
    {
        atmturb_fresnel_init_blocks(f_flt, phase, 0, nblocks);
    }

    for (long i = nblocks * 8; i < npix; i++)
    {
        sincosf(phase[i], &field[i][1], &field[i][0]);
    }
}

/**
 * atmturb_fresnel_extract_blocks - Extract unwrapped phase and amplitude over blocks [b0, b1)
 */
static inline void atmturb_fresnel_extract_blocks(
    const float *restrict f_flt,
    float       *restrict pha,
    float       *restrict amp,
    long                  b0,
    long                  b1)
{
    for (long b = b0; b < b1; b++)
    {
        long   i   = b * 8;
        long   idx = 2 * i;
        __m256 f0  = _mm256_loadu_ps(&f_flt[idx]);
        __m256 f1  = _mm256_loadu_ps(&f_flt[idx + 8]);

        __m256 re = _mm256_shuffle_ps(f0, f1, _MM_SHUFFLE(2, 0, 2, 0));
        __m256 im = _mm256_shuffle_ps(f0, f1, _MM_SHUFFLE(3, 1, 3, 1));
        re = _mm256_castpd_ps(
            _mm256_permute4x64_pd(_mm256_castps_pd(re), _MM_SHUFFLE(3, 1, 2, 0)));
        im = _mm256_castpd_ps(
            _mm256_permute4x64_pd(_mm256_castps_pd(im), _MM_SHUFFLE(3, 1, 2, 0)));

        if (amp != NULL)
        {
            __m256 sq = _mm256_fmadd_ps(re, re, _mm256_mul_ps(im, im));
            _mm256_storeu_ps(&amp[i], _mm256_sqrt_ps(sq));
        }

        __m256 vphi = _mm256_loadu_ps(&pha[i]);
        __m256 vs, vc;
        sincos256_ps(vphi, &vs, &vc);

        __m256 d_re = _mm256_fmadd_ps(re, vc, _mm256_mul_ps(im, vs));
        __m256 d_im = _mm256_fmsub_ps(im, vc, _mm256_mul_ps(re, vs));
        _mm256_storeu_ps(&pha[i], _mm256_add_ps(vphi, atan2_256_ps(d_im, d_re)));
    }
}

/**
 * atmturb_fresnel_extract_unwrapped_avx2 - AVX2 accelerated unwrapped extraction
 * @field: Complex optical field at pupil.
 * @pha: Phase array holding accumulated geometric phase (updated in-place).
 * @amp: Destination amplitude array (can be NULL).
 * @npix: Total number of pixels.
 */
void atmturb_fresnel_extract_unwrapped_avx2(
    const fftwf_complex *restrict field,
    float               *restrict pha,
    float               *restrict amp,
    long                          npix)
{
    const float *f_flt   = (const float *) field;
    long         nblocks = npix / 8;

#ifdef _OPENMP
    if (npix >= 262144 && !omp_in_parallel())
    {
        int max_t = omp_get_max_threads();
        int nth   = (max_t > 8) ? 8 : max_t;
        #pragma omp parallel num_threads(nth)
        {
            int  tid = omp_get_thread_num();
            int  nt  = omp_get_num_threads();
            long b0  = (nblocks * (long) tid) / (long) nt;
            long b1  = (nblocks * (long) (tid + 1)) / (long) nt;
            atmturb_fresnel_extract_blocks(f_flt, pha, amp, b0, b1);
        }
    }
    else
#endif
    {
        atmturb_fresnel_extract_blocks(f_flt, pha, amp, 0, nblocks);
    }

    for (long i = nblocks * 8; i < npix; i++)
    {
        float re = field[i][0];
        float im = field[i][1];
        if (amp != NULL)
        {
            amp[i] = sqrtf(re * re + im * im);
        }
        float phi_geo = pha[i];
        float s_geo, c_geo;
        sincosf(phi_geo, &s_geo, &c_geo);
        float d_re = re * c_geo + im * s_geo;
        float d_im = im * c_geo - re * s_geo;
        pha[i]     = phi_geo + atan2f(d_im, d_re);
    }
}

/**
 * log256_ps - Vectorized 8-lane natural logarithm for positive normalized floats (~1.2 ULP)
 */
static inline __m256 log256_ps(
    __m256 x)
{
    const __m256i mant_mask = _mm256_set1_epi32(0x007fffff);
    const __m256i one_bits  = _mm256_set1_epi32(0x3f800000);
    const __m256  sqrt_half = _mm256_set1_ps(0.707106781186547524f);
    const __m256  one       = _mm256_set1_ps(1.0f);

    __m256i ix = _mm256_castps_si256(x);
    __m256i ie = _mm256_sub_epi32(_mm256_srli_epi32(ix, 23), _mm256_set1_epi32(0x7f));
    __m256  e  = _mm256_cvtepi32_ps(ie);
    __m256  m  = _mm256_castsi256_ps(
        _mm256_or_si256(_mm256_and_si256(ix, mant_mask), _mm256_set1_epi32(0x3f000000)));

    __m256 mask = _mm256_cmp_ps(m, sqrt_half, _CMP_LT_OQ);
    __m256 addm = _mm256_and_ps(m, mask);
    m           = _mm256_sub_ps(m, one);
    e           = _mm256_add_ps(e, _mm256_andnot_ps(mask, one));
    m           = _mm256_add_ps(m, addm);
    (void) one_bits;

    __m256 z = _mm256_mul_ps(m, m);
    __m256 y = _mm256_set1_ps(7.0376836292e-2f);
    y        = _mm256_fmadd_ps(y, m, _mm256_set1_ps(-1.1514610310e-1f));
    y        = _mm256_fmadd_ps(y, m, _mm256_set1_ps(1.1676998740e-1f));
    y        = _mm256_fmadd_ps(y, m, _mm256_set1_ps(-1.2420140846e-1f));
    y        = _mm256_fmadd_ps(y, m, _mm256_set1_ps(1.4249322787e-1f));
    y        = _mm256_fmadd_ps(y, m, _mm256_set1_ps(-1.6668057665e-1f));
    y        = _mm256_fmadd_ps(y, m, _mm256_set1_ps(2.0000714765e-1f));
    y        = _mm256_fmadd_ps(y, m, _mm256_set1_ps(-2.4999993993e-1f));
    y        = _mm256_fmadd_ps(y, m, _mm256_set1_ps(3.3333331174e-1f));
    y        = _mm256_mul_ps(y, _mm256_mul_ps(m, z));
    y        = _mm256_fmadd_ps(e, _mm256_set1_ps(-2.12194440e-4f), y);
    y        = _mm256_fnmadd_ps(z, _mm256_set1_ps(0.5f), y);
    y        = _mm256_add_ps(m, y);
    return _mm256_fmadd_ps(e, _mm256_set1_ps(0.693359375f), y);
}

/**
 * exp256_ps - Vectorized 8-lane natural exponential with clean underflow zeroing (~1.2 ULP)
 */
static inline __m256 exp256_ps(
    __m256 x)
{
    const __m256 lo_lim = _mm256_set1_ps(-87.33654475f);
    const __m256 hi_lim = _mm256_set1_ps(88.37626266f);
    __m256       uflow  = _mm256_cmp_ps(x, lo_lim, _CMP_LT_OQ);
    __m256       xc     = _mm256_min_ps(_mm256_max_ps(x, lo_lim), hi_lim);

    __m256 fx = _mm256_round_ps(_mm256_mul_ps(xc, _mm256_set1_ps(1.44269504088896341f)),
                                _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    __m256 r  = _mm256_fnmadd_ps(fx, _mm256_set1_ps(0.693359375f), xc);
    r         = _mm256_fnmadd_ps(fx, _mm256_set1_ps(-2.12194440e-4f), r);

    __m256 z = _mm256_mul_ps(r, r);
    __m256 y = _mm256_set1_ps(1.9875691500e-4f);
    y        = _mm256_fmadd_ps(y, r, _mm256_set1_ps(1.3981999507e-3f));
    y        = _mm256_fmadd_ps(y, r, _mm256_set1_ps(8.3334519073e-3f));
    y        = _mm256_fmadd_ps(y, r, _mm256_set1_ps(4.1665795894e-2f));
    y        = _mm256_fmadd_ps(y, r, _mm256_set1_ps(1.6666665459e-1f));
    y        = _mm256_fmadd_ps(y, r, _mm256_set1_ps(5.0000001201e-1f));
    y        = _mm256_fmadd_ps(y, z, _mm256_add_ps(r, _mm256_set1_ps(1.0f)));

    __m256i ik    = _mm256_add_epi32(_mm256_cvtps_epi32(fx), _mm256_set1_epi32(0x7f));
    __m256  pow2k = _mm256_castsi256_ps(_mm256_slli_epi32(ik, 23));
    return _mm256_andnot_ps(uflow, _mm256_mul_ps(y, pow2k));
}

/**
 * atmturb_fresnel_phasor_masked_avx2 - AVX2 masked pupil complex field synthesis
 */
void atmturb_fresnel_phasor_masked_avx2(
    fftwf_complex *restrict field,
    const float   *restrict opd,
    const float   *restrict mask,
    float                   scale,
    long                    npix)
{
    static const int s_idx_lo[8] = {0, 0, 1, 1, 2, 2, 3, 3};
    static const int s_idx_hi[8] = {4, 4, 5, 5, 6, 6, 7, 7};
    __m256i          idx_lo      = _mm256_loadu_si256((const __m256i *) s_idx_lo);
    __m256i          idx_hi      = _mm256_loadu_si256((const __m256i *) s_idx_hi);
    float           *f_flt       = (float *) field;
    long             nblocks     = npix / 8;

    if (opd == NULL)
    {
        __m256 vz = _mm256_setzero_ps();
        for (long b = 0; b < nblocks; b++)
        {
            long   i    = b * 8;
            __m256 vm   = _mm256_loadu_ps(&mask[i]);
            __m256 m_lo = _mm256_permutevar8x32_ps(vm, idx_lo);
            __m256 m_hi = _mm256_permutevar8x32_ps(vm, idx_hi);
            _mm256_storeu_ps(&f_flt[2 * i],     _mm256_blend_ps(m_lo, vz, 0xAA));
            _mm256_storeu_ps(&f_flt[2 * i + 8], _mm256_blend_ps(m_hi, vz, 0xAA));
        }
        atmturb_fresnel_phasor_masked_scalar(&field[nblocks * 8], NULL, &mask[nblocks * 8],
                                             scale, npix - nblocks * 8);
        return;
    }

    __m256 vscale = _mm256_set1_ps(scale);
    for (long b = 0; b < nblocks; b++)
    {
        long   i  = b * 8;
        __m256 vm = _mm256_loadu_ps(&mask[i]);
        __m256 vp = _mm256_mul_ps(_mm256_loadu_ps(&opd[i]), vscale);
        __m256 vs, vc;
        sincos256_ps(vp, &vs, &vc);
        vc = _mm256_mul_ps(vc, vm);
        vs = _mm256_mul_ps(vs, vm);

        __m256 c_lo = _mm256_permutevar8x32_ps(vc, idx_lo);
        __m256 s_lo = _mm256_permutevar8x32_ps(vs, idx_lo);
        __m256 c_hi = _mm256_permutevar8x32_ps(vc, idx_hi);
        __m256 s_hi = _mm256_permutevar8x32_ps(vs, idx_hi);
        _mm256_storeu_ps(&f_flt[2 * i],     _mm256_blend_ps(c_lo, s_lo, 0xAA));
        _mm256_storeu_ps(&f_flt[2 * i + 8], _mm256_blend_ps(c_hi, s_hi, 0xAA));
    }
    atmturb_fresnel_phasor_masked_scalar(&field[nblocks * 8], &opd[nblocks * 8],
                                         &mask[nblocks * 8], scale, npix - nblocks * 8);
}

/**
 * atmturb_fresnel_spectrum_strip_avx2 - AVX2 8-lane von Karman Fourier mode synthesis
 */
void atmturb_fresnel_spectrum_strip_avx2(
    fftwf_complex *restrict out,
    const float   *restrict u1,
    const float   *restrict u2,
    float                   fx0,
    const float   *restrict spec_p)
{
    static const int   s_idx_lo[8] = {0, 0, 1, 1, 2, 2, 3, 3};
    static const int   s_idx_hi[8] = {4, 4, 5, 5, 6, 6, 7, 7};
    static const float s_lanes[8]  = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f};

    __m256 vfx    = _mm256_add_ps(_mm256_set1_ps(fx0), _mm256_loadu_ps(s_lanes));
    __m256 vr2    = _mm256_fmadd_ps(vfx, vfx, _mm256_set1_ps(spec_p[0]));
    __m256 vdist2 = _mm256_add_ps(vr2, _mm256_set1_ps(spec_p[1]));

    __m256 vrbm   = _mm256_sqrt_ps(_mm256_mul_ps(_mm256_set1_ps(-2.0f),
                                                 log256_ps(_mm256_loadu_ps(u1))));
    __m256 vtheta = _mm256_mul_ps(_mm256_set1_ps(2.0f * (float) M_PI), _mm256_loadu_ps(u2));
    __m256 vs, vc;
    sincos256_ps(vtheta, &vs, &vc);

    __m256 vlog_d = log256_ps(vdist2);
    __m256 vexp_a = _mm256_mul_ps(_mm256_set1_ps(-11.0f / 12.0f), vlog_d);
    if (spec_p[3] > 0.0f)
    {
        vexp_a = _mm256_fnmadd_ps(vr2, _mm256_set1_ps(spec_p[3]), vexp_a);
    }
    __m256 vamp = _mm256_mul_ps(_mm256_set1_ps(spec_p[2]), exp256_ps(vexp_a));
    __m256 vmag = _mm256_mul_ps(vamp, vrbm);
    __m256 vre  = _mm256_mul_ps(vmag, vc);
    __m256 vim  = _mm256_mul_ps(vmag, vs);

    __m256i idx_lo = _mm256_loadu_si256((const __m256i *) s_idx_lo);
    __m256i idx_hi = _mm256_loadu_si256((const __m256i *) s_idx_hi);
    float  *f_out  = (float *) out;
    _mm256_storeu_ps(&f_out[0], _mm256_blend_ps(_mm256_permutevar8x32_ps(vre, idx_lo),
                                                _mm256_permutevar8x32_ps(vim, idx_lo), 0xAA));
    _mm256_storeu_ps(&f_out[8], _mm256_blend_ps(_mm256_permutevar8x32_ps(vre, idx_hi),
                                                _mm256_permutevar8x32_ps(vim, idx_hi), 0xAA));
}

#else /* !defined(__AVX2__) */

void atmturb_fresnel_modulate_field_avx2(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    npix)
{
    atmturb_fresnel_modulate_field_scalar(field, phase, npix);
}

void atmturb_fresnel_init_modulate_field_avx2(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    npix)
{
    atmturb_fresnel_init_modulate_field_scalar(field, phase, npix);
}

void atmturb_fresnel_extract_unwrapped_avx2(
    const fftwf_complex *restrict field,
    float               *restrict pha,
    float               *restrict amp,
    long                          npix)
{
    atmturb_fresnel_extract_unwrapped_scalar(field, pha, amp, npix);
}

void atmturb_fresnel_phasor_masked_avx2(
    fftwf_complex *restrict field,
    const float   *restrict opd,
    const float   *restrict mask,
    float                   scale,
    long                    npix)
{
    atmturb_fresnel_phasor_masked_scalar(field, opd, mask, scale, npix);
}

void atmturb_fresnel_spectrum_strip_avx2(
    fftwf_complex *restrict out,
    const float   *restrict u1,
    const float   *restrict u2,
    float                   fx0,
    const float   *restrict spec_p)
{
    atmturb_fresnel_spectrum_strip_scalar(out, u1, u2, fx0, spec_p);
}

#endif
