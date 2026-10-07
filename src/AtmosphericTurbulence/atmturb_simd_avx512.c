// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_simd_avx512.c
 * @brief   AVX-512 accelerated compute kernels for atmospheric turbulence simulation
 */

#include <math.h>
#if defined(__x86_64__) || defined(_M_X64)
#    include <immintrin.h>
#endif
#include "atmturb_simd.h"

#if defined(__AVX512F__)

static inline void atmturb_keys_weights(
    float f,
    float w[4])
{
    float f2 = f * f;
    float f3 = f2 * f;
    w[0] = -0.5f * f + f2 - 0.5f * f3;
    w[1] = 1.0f - 2.5f * f2 + 1.5f * f3;
    w[2] = 0.5f * f + 2.0f * f2 - 1.5f * f3;
    w[3] = -0.5f * f2 + 0.5f * f3;
}

static inline __m512 atmturb_bilinear_tap_avx512_os1(
    const float *row0,
    const float *row1,
    long         sx,
    __m512       vw00,
    __m512       vw10,
    __m512       vw01,
    __m512       vw11)
{
    __m512 v00 = _mm512_loadu_ps(&row0[sx]);
    __m512 v10 = _mm512_loadu_ps(&row0[sx + 1]);
    __m512 v01 = _mm512_loadu_ps(&row1[sx]);
    __m512 v11 = _mm512_loadu_ps(&row1[sx + 1]);
    __m512 acc = _mm512_mul_ps(v00, vw00);
    acc = _mm512_fmadd_ps(v10, vw10, acc);
    acc = _mm512_fmadd_ps(v01, vw01, acc);
    return _mm512_fmadd_ps(v11, vw11, acc);
}

static inline __m512 atmturb_bilinear_tap_avx512_os2(
    const float *row0,
    const float *row1,
    long         sx,
    __m512i      vidx,
    __m512       vw00,
    __m512       vw10,
    __m512       vw01,
    __m512       vw11)
{
    __m512 v00 = _mm512_i32gather_ps(vidx, (const void *) &row0[sx], 4);
    __m512 v10 = _mm512_i32gather_ps(vidx, (const void *) &row0[sx + 1], 4);
    __m512 v01 = _mm512_i32gather_ps(vidx, (const void *) &row1[sx], 4);
    __m512 v11 = _mm512_i32gather_ps(vidx, (const void *) &row1[sx + 1], 4);
    __m512 acc = _mm512_mul_ps(v00, vw00);
    acc = _mm512_fmadd_ps(v10, vw10, acc);
    acc = _mm512_fmadd_ps(v01, vw01, acc);
    return _mm512_fmadd_ps(v11, vw11, acc);
}

static inline __m512 atmturb_bicubic_tap_avx512_os1(
    const float *row,
    long         sx,
    __m512       vwx0,
    __m512       vwx1,
    __m512       vwx2,
    __m512       vwx3)
{
    __m512 v0 = _mm512_loadu_ps(&row[sx - 1]);
    __m512 v1 = _mm512_loadu_ps(&row[sx]);
    __m512 v2 = _mm512_loadu_ps(&row[sx + 1]);
    __m512 v3 = _mm512_loadu_ps(&row[sx + 2]);
    __m512 h = _mm512_mul_ps(v0, vwx0);
    h = _mm512_fmadd_ps(v1, vwx1, h);
    h = _mm512_fmadd_ps(v2, vwx2, h);
    return _mm512_fmadd_ps(v3, vwx3, h);
}

static inline __m512 atmturb_bicubic_tap_avx512_os2(
    const float *row,
    long         sx,
    __m512i      vidx,
    __m512       vwx0,
    __m512       vwx1,
    __m512       vwx2,
    __m512       vwx3)
{
    __m512 v0 = _mm512_i32gather_ps(vidx, (const void *) &row[sx - 1], 4);
    __m512 v1 = _mm512_i32gather_ps(vidx, (const void *) &row[sx], 4);
    __m512 v2 = _mm512_i32gather_ps(vidx, (const void *) &row[sx + 1], 4);
    __m512 v3 = _mm512_i32gather_ps(vidx, (const void *) &row[sx + 2], 4);
    __m512 h = _mm512_mul_ps(v0, vwx0);
    h = _mm512_fmadd_ps(v1, vwx1, h);
    h = _mm512_fmadd_ps(v2, vwx2, h);
    return _mm512_fmadd_ps(v3, vwx3, h);
}

/**
 * atmturb_extrude_accumulate_bilinear_avx512 - AVX-512 accelerated bilinear extrusion
 * @params: Extrusion configuration and data pointers.
 */
void atmturb_extrude_accumulate_bilinear_avx512(
    const atmturb_extrude_params_t *params)
{
    long msize = params->msize;
    long pup_size = params->pup_size;
    long os = (params->os > 1) ? params->os : 1;
    float weight = params->weight;
    const float *master = params->master;
    float *out_pha = params->out_pha;

    float fx = (float)(params->x0 - floor(params->x0));
    float fy = (float)(params->y0 - floor(params->y0));
    long base_x = (long)floor(params->x0);
    long base_y = (long)floor(params->y0);

    __m512 vw00 = _mm512_set1_ps((1.0f - fx) * (1.0f - fy) * weight);
    __m512 vw10 = _mm512_set1_ps(fx * (1.0f - fy) * weight);
    __m512 vw01 = _mm512_set1_ps((1.0f - fx) * fy * weight);
    __m512 vw11 = _mm512_set1_ps(fx * fy * weight);

    static const int s_idx[16] = {
        0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30
    };
    __m512i vidx = _mm512_loadu_si512((const __m512i *) s_idx);

    long start_x = base_x % msize;
    if (start_x < 0) start_x += msize;

    for (long jj = 0; jj < pup_size; jj++)
    {
        long iy0 = (base_y + jj * os) % msize;
        if (iy0 < 0) iy0 += msize;
        long iy1 = (iy0 + 1) % msize;

        const float *row0 = &master[iy0 * msize];
        const float *row1 = &master[iy1 * msize];
        float *out_row = &out_pha[jj * pup_size];

        if (start_x + (pup_size - 1) * os + 1 < msize)
        {
            long ii = 0;
            if (os == 1)
            {
                for (; ii <= pup_size - 16; ii += 16)
                {
                    __m512 acc = atmturb_bilinear_tap_avx512_os1(row0, row1, start_x + ii,
                                                                 vw00, vw10, vw01, vw11);
                    __m512 out = _mm512_loadu_ps(&out_row[ii]);
                    _mm512_storeu_ps(&out_row[ii], _mm512_add_ps(out, acc));
                }
            }
            else if (os == 2)
            {
                for (; ii <= pup_size - 16; ii += 16)
                {
                    __m512 acc = atmturb_bilinear_tap_avx512_os2(row0, row1, start_x + 2 * ii,
                                                                 vidx, vw00, vw10, vw01, vw11);
                    __m512 out = _mm512_loadu_ps(&out_row[ii]);
                    _mm512_storeu_ps(&out_row[ii], _mm512_add_ps(out, acc));
                }
            }
            float w00 = _mm512_cvtss_f32(vw00), w10 = _mm512_cvtss_f32(vw10);
            float w01 = _mm512_cvtss_f32(vw01), w11 = _mm512_cvtss_f32(vw11);
            for (; ii < pup_size; ii++)
            {
                long ix0 = start_x + ii * os;
                out_row[ii] += w00 * row0[ix0] + w10 * row0[ix0 + 1] +
                               w01 * row1[ix0] + w11 * row1[ix0 + 1];
            }
        }
        else
        {
            float w00 = _mm512_cvtss_f32(vw00), w10 = _mm512_cvtss_f32(vw10);
            float w01 = _mm512_cvtss_f32(vw01), w11 = _mm512_cvtss_f32(vw11);
            for (long ii = 0; ii < pup_size; ii++)
            {
                long ix0 = (start_x + ii * os) % msize;
                if (ix0 < 0) ix0 += msize;
                long ix1 = (ix0 + 1) % msize;
                out_row[ii] += w00 * row0[ix0] + w10 * row0[ix1] +
                               w01 * row1[ix0] + w11 * row1[ix1];
            }
        }
    }
}

/**
 * atmturb_extrude_accumulate_bicubic_avx512 - AVX-512 accelerated Keys bicubic extrusion
 * @params: Extrusion configuration and data pointers.
 */
void atmturb_extrude_accumulate_bicubic_avx512(
    const atmturb_extrude_params_t *params)
{
    long msize = params->msize;
    long pup_size = params->pup_size;
    long os = (params->os > 1) ? params->os : 1;
    float weight = params->weight;
    const float *master = params->master;
    float *out_pha = params->out_pha;

    float fx = (float)(params->x0 - floor(params->x0));
    float fy = (float)(params->y0 - floor(params->y0));
    long base_x = (long)floor(params->x0);
    long base_y = (long)floor(params->y0);

    float wx[4], wy[4];
    atmturb_keys_weights(fx, wx);
    atmturb_keys_weights(fy, wy);

    __m512 vwx0 = _mm512_set1_ps(wx[0]), vwx1 = _mm512_set1_ps(wx[1]);
    __m512 vwx2 = _mm512_set1_ps(wx[2]), vwx3 = _mm512_set1_ps(wx[3]);
    __m512 vwy0 = _mm512_set1_ps(wy[0] * weight), vwy1 = _mm512_set1_ps(wy[1] * weight);
    __m512 vwy2 = _mm512_set1_ps(wy[2] * weight), vwy3 = _mm512_set1_ps(wy[3] * weight);

    static const int s_idx[16] = {
        0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30
    };
    __m512i vidx = _mm512_loadu_si512((const __m512i *) s_idx);

    long start_x = base_x % msize;
    if (start_x < 0) start_x += msize;

    for (long jj = 0; jj < pup_size; jj++)
    {
        long by = base_y + jj * os;
        long iy[4];
        for (int n = 0; n < 4; n++)
        {
            long y = (by + n - 1) % msize;
            if (y < 0) y += msize;
            iy[n] = y;
        }

        const float *row0 = &master[iy[0] * msize];
        const float *row1 = &master[iy[1] * msize];
        const float *row2 = &master[iy[2] * msize];
        const float *row3 = &master[iy[3] * msize];
        float *out_row = &out_pha[jj * pup_size];

        if (start_x >= 1 && start_x + (pup_size - 1) * os + 2 < msize)
        {
            long ii = 0;
            if (os == 1)
            {
                for (; ii <= pup_size - 16; ii += 16)
                {
                    __m512 h0 = atmturb_bicubic_tap_avx512_os1(row0, start_x + ii,
                                                               vwx0, vwx1, vwx2, vwx3);
                    __m512 h1 = atmturb_bicubic_tap_avx512_os1(row1, start_x + ii,
                                                               vwx0, vwx1, vwx2, vwx3);
                    __m512 h2 = atmturb_bicubic_tap_avx512_os1(row2, start_x + ii,
                                                               vwx0, vwx1, vwx2, vwx3);
                    __m512 h3 = atmturb_bicubic_tap_avx512_os1(row3, start_x + ii,
                                                               vwx0, vwx1, vwx2, vwx3);
                    __m512 acc = _mm512_mul_ps(h0, vwy0);
                    acc = _mm512_fmadd_ps(h1, vwy1, acc);
                    acc = _mm512_fmadd_ps(h2, vwy2, acc);
                    acc = _mm512_fmadd_ps(h3, vwy3, acc);
                    __m512 out = _mm512_loadu_ps(&out_row[ii]);
                    _mm512_storeu_ps(&out_row[ii], _mm512_add_ps(out, acc));
                }
            }
            else if (os == 2)
            {
                for (; ii <= pup_size - 16; ii += 16)
                {
                    __m512 h0 = atmturb_bicubic_tap_avx512_os2(row0, start_x + 2 * ii, vidx,
                                                               vwx0, vwx1, vwx2, vwx3);
                    __m512 h1 = atmturb_bicubic_tap_avx512_os2(row1, start_x + 2 * ii, vidx,
                                                               vwx0, vwx1, vwx2, vwx3);
                    __m512 h2 = atmturb_bicubic_tap_avx512_os2(row2, start_x + 2 * ii, vidx,
                                                               vwx0, vwx1, vwx2, vwx3);
                    __m512 h3 = atmturb_bicubic_tap_avx512_os2(row3, start_x + 2 * ii, vidx,
                                                               vwx0, vwx1, vwx2, vwx3);
                    __m512 acc = _mm512_mul_ps(h0, vwy0);
                    acc = _mm512_fmadd_ps(h1, vwy1, acc);
                    acc = _mm512_fmadd_ps(h2, vwy2, acc);
                    acc = _mm512_fmadd_ps(h3, vwy3, acc);
                    __m512 out = _mm512_loadu_ps(&out_row[ii]);
                    _mm512_storeu_ps(&out_row[ii], _mm512_add_ps(out, acc));
                }
            }
            float wy0 = wy[0] * weight, wy1 = wy[1] * weight;
            float wy2 = wy[2] * weight, wy3 = wy[3] * weight;
            for (; ii < pup_size; ii++)
            {
                long bx = start_x + ii * os;
                float h0 = wx[0] * row0[bx - 1] + wx[1] * row0[bx] +
                           wx[2] * row0[bx + 1] + wx[3] * row0[bx + 2];
                float h1 = wx[0] * row1[bx - 1] + wx[1] * row1[bx] +
                           wx[2] * row1[bx + 1] + wx[3] * row1[bx + 2];
                float h2 = wx[0] * row2[bx - 1] + wx[1] * row2[bx] +
                           wx[2] * row2[bx + 1] + wx[3] * row2[bx + 2];
                float h3 = wx[0] * row3[bx - 1] + wx[1] * row3[bx] +
                           wx[2] * row3[bx + 2] + wx[3] * row3[bx + 2];
                out_row[ii] += wy0 * h0 + wy1 * h1 + wy2 * h2 + wy3 * h3;
            }
        }
        else
        {
            float wy0 = wy[0] * weight, wy1 = wy[1] * weight;
            float wy2 = wy[2] * weight, wy3 = wy[3] * weight;
            for (long ii = 0; ii < pup_size; ii++)
            {
                long bx = start_x + ii * os;
                long ix0 = (bx - 1) % msize; if (ix0 < 0) ix0 += msize;
                long ix1 = bx % msize;       if (ix1 < 0) ix1 += msize;
                long ix2 = (bx + 1) % msize; if (ix2 < 0) ix2 += msize;
                long ix3 = (bx + 2) % msize; if (ix3 < 0) ix3 += msize;

                float h0 = wx[0] * row0[ix0] + wx[1] * row0[ix1] +
                           wx[2] * row0[ix2] + wx[3] * row0[ix3];
                float h1 = wx[0] * row1[ix0] + wx[1] * row1[ix1] +
                           wx[2] * row1[ix2] + wx[3] * row1[ix3];
                float h2 = wx[0] * row2[ix0] + wx[1] * row2[ix1] +
                           wx[2] * row2[ix2] + wx[3] * row2[ix3];
                float h3 = wx[0] * row3[ix0] + wx[1] * row3[ix1] +
                           wx[2] * row3[ix2] + wx[3] * row3[ix3];

                out_row[ii] += wy0 * h0 + wy1 * h1 + wy2 * h2 + wy3 * h3;
            }
        }
    }
}

/**
 * atmturb_extrude_accumulate_avx512 - AVX-512 extrusion dispatcher (bilinear/bicubic)
 * @params: Extrusion configuration and data pointers.
 */
void atmturb_extrude_accumulate_avx512(
    const atmturb_extrude_params_t *params)
{
    if (params->interp == ATMTURB_INTERP_BICUBIC)
    {
        atmturb_extrude_accumulate_bicubic_avx512(params);
    }
    else
    {
        atmturb_extrude_accumulate_bilinear_avx512(params);
    }
}

/**
 * atmturb_scale_float_array_avx512 - Vectorized array multiplication by scalar (AVX-512)
 * @dest: Output float array.
 * @src: Input float array.
 * @scale: Scalar multiplier.
 * @n: Number of elements.
 */
void atmturb_scale_float_array_avx512(
    float       *dest,
    const float *src,
    float        scale,
    long         n)
{
    __m512 vscale = _mm512_set1_ps(scale);
    long i = 0;
    for (; i <= n - 16; i += 16)
    {
        __m512 v = _mm512_loadu_ps(&src[i]);
        _mm512_storeu_ps(&dest[i], _mm512_mul_ps(v, vscale));
    }
    for (; i < n; i++)
    {
        dest[i] = src[i] * scale;
    }
}

/**
 * atmturb_init_phase_amp_avx512 - Vectorized phase zeroing and amplitude one-filling (AVX-512)
 * @pha: Phase array.
 * @amp: Amplitude array.
 * @n: Number of elements.
 */
void atmturb_init_phase_amp_avx512(
    float *pha,
    float *amp,
    long   n)
{
    __m512 vzero = _mm512_setzero_ps();
    __m512 vone = _mm512_set1_ps(1.0f);
    long i = 0;
    for (; i <= n - 16; i += 16)
    {
        _mm512_storeu_ps(&pha[i], vzero);
        _mm512_storeu_ps(&amp[i], vone);
    }
    for (; i < n; i++)
    {
        pha[i] = 0.0f;
        amp[i] = 1.0f;
    }
}

#else

void atmturb_extrude_accumulate_avx512(
    const atmturb_extrude_params_t *params)
{
    atmturb_extrude_accumulate_scalar(params);
}

void atmturb_extrude_accumulate_bilinear_avx512(
    const atmturb_extrude_params_t *params)
{
    atmturb_extrude_accumulate_bilinear_scalar(params);
}

void atmturb_extrude_accumulate_bicubic_avx512(
    const atmturb_extrude_params_t *params)
{
    atmturb_extrude_accumulate_bicubic_scalar(params);
}

void atmturb_scale_float_array_avx512(
    float       *dest,
    const float *src,
    float        scale,
    long         n)
{
    atmturb_scale_float_array_scalar(dest, src, scale, n);
}

void atmturb_init_phase_amp_avx512(
    float *pha,
    float *amp,
    long   n)
{
    atmturb_init_phase_amp_scalar(pha, amp, n);
}

#endif
