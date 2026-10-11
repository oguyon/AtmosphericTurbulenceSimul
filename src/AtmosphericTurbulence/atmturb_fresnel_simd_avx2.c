// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_fresnel_simd_avx2.c
 * @brief   AVX2-accelerated field modulation and extraction kernels for Fresnel propagation
 */

#define _GNU_SOURCE
#include <math.h>

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
    const __m256 p1 = _mm256_set1_ps(1.5703125f);
    const __m256 p2 = _mm256_set1_ps(0.0004837512969970703125f);
    const __m256 p3 = _mm256_set1_ps(7.549789948768648e-8f);

    __m256 v_k = _mm256_round_ps(_mm256_mul_ps(v, inv_pio2),
                                 _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    __m256i vi_k = _mm256_cvtps_epi32(v_k);

    __m256 r = _mm256_fnmadd_ps(v_k, p1, v);
    r = _mm256_fnmadd_ps(v_k, p2, r);
    r = _mm256_fnmadd_ps(v_k, p3, r);

    __m256 z = _mm256_mul_ps(r, r);

    const __m256 s3 = _mm256_set1_ps(-1.9515295891e-4f);
    const __m256 s2 = _mm256_set1_ps(+8.3321608736e-3f);
    const __m256 s1 = _mm256_set1_ps(-1.6666654611e-1f);
    __m256 poly_s = _mm256_fmadd_ps(z, s3, s2);
    poly_s = _mm256_fmadd_ps(z, poly_s, s1);
    poly_s = _mm256_fmadd_ps(_mm256_mul_ps(poly_s, z), r, r);

    const __m256 c3 = _mm256_set1_ps(-1.38537042643e-3f);
    const __m256 c2 = _mm256_set1_ps(+4.16635846761e-2f);
    const __m256 c1 = _mm256_set1_ps(-4.99999053455e-1f);
    __m256 poly_c = _mm256_fmadd_ps(z, c3, c2);
    poly_c = _mm256_fmadd_ps(z, poly_c, c1);
    poly_c = _mm256_fmadd_ps(poly_c, z, _mm256_set1_ps(1.0f));

    __m256i swap_mask = _mm256_slli_epi32(_mm256_and_si256(vi_k, _mm256_set1_epi32(1)), 31);
    __m256 v_sin = _mm256_blendv_ps(poly_s, poly_c, _mm256_castsi256_ps(swap_mask));
    __m256 v_cos = _mm256_blendv_ps(poly_c, poly_s, _mm256_castsi256_ps(swap_mask));

    __m256i sign_s = _mm256_slli_epi32(_mm256_and_si256(vi_k, _mm256_set1_epi32(2)), 30);
    __m256i sign_c = _mm256_slli_epi32(_mm256_and_si256(
        _mm256_add_epi32(vi_k, _mm256_set1_epi32(1)), _mm256_set1_epi32(2)), 30);

    *s = _mm256_xor_ps(v_sin, _mm256_castsi256_ps(sign_s));
    *c = _mm256_xor_ps(v_cos, _mm256_castsi256_ps(sign_c));
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
    float *f_flt = (float *) field;
    static const int s_idx_lo[8] = {0, 0, 1, 1, 2, 2, 3, 3};
    static const int s_idx_hi[8] = {4, 4, 5, 5, 6, 6, 7, 7};
    __m256i idx_lo = _mm256_loadu_si256((const __m256i *) s_idx_lo);
    __m256i idx_hi = _mm256_loadu_si256((const __m256i *) s_idx_hi);

    long i = 0;
    for (; i <= npix - 8; i += 8)
    {
        __m256 vp = _mm256_loadu_ps(&phase[i]);
        __m256 vs, vc;
        sincos256_ps(vp, &vs, &vc);

        __m256 c_lo = _mm256_permutevar8x32_ps(vc, idx_lo);
        __m256 s_lo = _mm256_permutevar8x32_ps(vs, idx_lo);
        __m256 c_hi = _mm256_permutevar8x32_ps(vc, idx_hi);
        __m256 s_hi = _mm256_permutevar8x32_ps(vs, idx_hi);

        __m256 cs_lo = _mm256_blend_ps(c_lo, s_lo, 0xAA);
        __m256 sc_lo = _mm256_blend_ps(s_lo, c_lo, 0xAA);
        __m256 cs_hi = _mm256_blend_ps(c_hi, s_hi, 0xAA);
        __m256 sc_hi = _mm256_blend_ps(s_hi, c_hi, 0xAA);

        long idx = 2 * i;
        __m256 f0 = _mm256_loadu_ps(&f_flt[idx]);
        __m256 f1 = _mm256_loadu_ps(&f_flt[idx + 8]);

        __m256 re0 = _mm256_moveldup_ps(f0);
        __m256 im0 = _mm256_movehdup_ps(f0);
        __m256 re1 = _mm256_moveldup_ps(f1);
        __m256 im1 = _mm256_movehdup_ps(f1);

        __m256 res0 = _mm256_fmaddsub_ps(re0, cs_lo, _mm256_mul_ps(im0, sc_lo));
        __m256 res1 = _mm256_fmaddsub_ps(re1, cs_hi, _mm256_mul_ps(im1, sc_hi));

        _mm256_storeu_ps(&f_flt[idx],     res0);
        _mm256_storeu_ps(&f_flt[idx + 8], res1);
    }
    for (; i < npix; i++)
    {
        float p = phase[i];
        float s, c;
        sincosf(p, &s, &c);
        float re = field[i][0];
        float im = field[i][1];
        field[i][0] = re * c - im * s;
        field[i][1] = re * s + im * c;
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
    float *f_flt = (float *) field;
    static const int s_idx_lo[8] = {0, 0, 1, 1, 2, 2, 3, 3};
    static const int s_idx_hi[8] = {4, 4, 5, 5, 6, 6, 7, 7};
    __m256i idx_lo = _mm256_loadu_si256((const __m256i *) s_idx_lo);
    __m256i idx_hi = _mm256_loadu_si256((const __m256i *) s_idx_hi);

    long i = 0;
    for (; i <= npix - 8; i += 8)
    {
        __m256 vp = _mm256_loadu_ps(&phase[i]);
        __m256 vs, vc;
        sincos256_ps(vp, &vs, &vc);

        __m256 c_lo = _mm256_permutevar8x32_ps(vc, idx_lo);
        __m256 s_lo = _mm256_permutevar8x32_ps(vs, idx_lo);
        __m256 c_hi = _mm256_permutevar8x32_ps(vc, idx_hi);
        __m256 s_hi = _mm256_permutevar8x32_ps(vs, idx_hi);

        __m256 cs_lo = _mm256_blend_ps(c_lo, s_lo, 0xAA);
        __m256 cs_hi = _mm256_blend_ps(c_hi, s_hi, 0xAA);

        long idx = 2 * i;
        _mm256_storeu_ps(&f_flt[idx],     cs_lo);
        _mm256_storeu_ps(&f_flt[idx + 8], cs_hi);
    }
    for (; i < npix; i++)
    {
        float p = phase[i];
        float s, c;
        sincosf(p, &s, &c);
        field[i][0] = c;
        field[i][1] = s;
    }
}

/**
 * atan2_256_ps - Vectorized 8-wide atan2 using 15th-order minimax polynomial
 * @y: Numerator vector.
 * @x: Denominator vector.
 *
 * Return: 8-wide vector of atan2(y, x) in radians within 4.8e-7 max error.
 */
static inline __m256 atan2_256_ps(
    __m256 y,
    __m256 x)
{
    __m256 sign_mask = _mm256_set1_ps(-0.0f);
    __m256 ax        = _mm256_andnot_ps(sign_mask, x);
    __m256 ay        = _mm256_andnot_ps(sign_mask, y);

    __m256 cmp       = _mm256_cmp_ps(ay, ax, _CMP_GT_OQ);
    __m256 num       = _mm256_blendv_ps(ay, ax, cmp);
    __m256 den       = _mm256_blendv_ps(ax, ay, cmp);
    __m256 den_safe  = _mm256_max_ps(den, _mm256_set1_ps(1e-30f));
    __m256 z         = _mm256_div_ps(num, den_safe);
    __m256 z2        = _mm256_mul_ps(z, z);

    __m256 poly = _mm256_set1_ps(-4.2577745724e-03f);
    poly = _mm256_fmadd_ps(poly, z2, _mm256_set1_ps(2.2568064116e-02f));
    poly = _mm256_fmadd_ps(poly, z2, _mm256_set1_ps(-5.6882557658e-02f));
    poly = _mm256_fmadd_ps(poly, z2, _mm256_set1_ps(9.7094637892e-02f));
    poly = _mm256_fmadd_ps(poly, z2, _mm256_set1_ps(-1.3933270093e-01f));
    poly = _mm256_fmadd_ps(poly, z2, _mm256_set1_ps(1.9951118324e-01f));
    poly = _mm256_fmadd_ps(poly, z2, _mm256_set1_ps(-3.3330222945e-01f));
    poly = _mm256_fmadd_ps(poly, z2, _mm256_set1_ps(9.9999941614e-01f));

    __m256 at = _mm256_mul_ps(poly, z);

    __m256 piover2 = _mm256_set1_ps(1.5707963267948966f);
    __m256 angle   = _mm256_blendv_ps(at, _mm256_sub_ps(piover2, at), cmp);

    __m256 pi    = _mm256_set1_ps(3.1415926535897932f);
    __m256 x_neg = _mm256_cmp_ps(x, _mm256_setzero_ps(), _CMP_LT_OQ);
    angle        = _mm256_blendv_ps(angle, _mm256_sub_ps(pi, angle), x_neg);

    __m256 y_sign  = _mm256_and_ps(y, sign_mask);
    __m256 res     = _mm256_xor_ps(angle, y_sign);

    __m256 zero_pt = _mm256_and_ps(_mm256_cmp_ps(x, _mm256_setzero_ps(), _CMP_EQ_OQ),
                                   _mm256_cmp_ps(y, _mm256_setzero_ps(), _CMP_EQ_OQ));
    return _mm256_blendv_ps(res, _mm256_setzero_ps(), zero_pt);
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
    const float *f_flt = (const float *) field;

    long i = 0;
    for (; i <= npix - 8; i += 8)
    {
        long idx = 2 * i;
        __m256 f0 = _mm256_loadu_ps(&f_flt[idx]);
        __m256 f1 = _mm256_loadu_ps(&f_flt[idx + 8]);

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
        __m256 dphi = atan2_256_ps(d_im, d_re);
        _mm256_storeu_ps(&pha[i], _mm256_add_ps(vphi, dphi));
    }
    for (; i < npix; i++)
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

#else

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

#endif
