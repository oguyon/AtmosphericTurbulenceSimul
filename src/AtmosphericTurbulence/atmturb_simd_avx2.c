// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_simd_avx2.c
 * @brief   AVX2-accelerated compute kernels for atmospheric turbulence simulation
 */

#define _GNU_SOURCE
#include <math.h>
#include <string.h>
#if defined(__x86_64__) || defined(_M_X64)
#    include <immintrin.h>
#endif
#include "atmturb_lowfreq.h"
#include "atmturb_simd.h"

#if defined(__AVX2__)

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

static inline __m256 atmturb_bilinear_tap_avx2_os1(
    const float *row0,
    const float *row1,
    long         sx,
    __m256       vw00,
    __m256       vw10,
    __m256       vw01,
    __m256       vw11)
{
    __m256 v00 = _mm256_loadu_ps(&row0[sx]);
    __m256 v10 = _mm256_loadu_ps(&row0[sx + 1]);
    __m256 v01 = _mm256_loadu_ps(&row1[sx]);
    __m256 v11 = _mm256_loadu_ps(&row1[sx + 1]);
    __m256 acc = _mm256_mul_ps(v00, vw00);
    acc = _mm256_fmadd_ps(v10, vw10, acc);
    acc = _mm256_fmadd_ps(v01, vw01, acc);
    return _mm256_fmadd_ps(v11, vw11, acc);
}

static inline __m256 atmturb_bilinear_tap_avx2_os2(
    const float *row0,
    const float *row1,
    long         sx,
    __m256i      vidx,
    __m256       vw00,
    __m256       vw10,
    __m256       vw01,
    __m256       vw11)
{
    __m256 v00 = _mm256_i32gather_ps(&row0[sx], vidx, 4);
    __m256 v10 = _mm256_i32gather_ps(&row0[sx + 1], vidx, 4);
    __m256 v01 = _mm256_i32gather_ps(&row1[sx], vidx, 4);
    __m256 v11 = _mm256_i32gather_ps(&row1[sx + 1], vidx, 4);
    __m256 acc = _mm256_mul_ps(v00, vw00);
    acc = _mm256_fmadd_ps(v10, vw10, acc);
    acc = _mm256_fmadd_ps(v01, vw01, acc);
    return _mm256_fmadd_ps(v11, vw11, acc);
}

static inline __m256 atmturb_bicubic_tap_avx2_os1(
    const float *row,
    long         sx,
    __m256       vwx0,
    __m256       vwx1,
    __m256       vwx2,
    __m256       vwx3)
{
    __m256 v0 = _mm256_loadu_ps(&row[sx - 1]);
    __m256 v1 = _mm256_loadu_ps(&row[sx]);
    __m256 v2 = _mm256_loadu_ps(&row[sx + 1]);
    __m256 v3 = _mm256_loadu_ps(&row[sx + 2]);
    __m256 h = _mm256_mul_ps(v0, vwx0);
    h = _mm256_fmadd_ps(v1, vwx1, h);
    h = _mm256_fmadd_ps(v2, vwx2, h);
    return _mm256_fmadd_ps(v3, vwx3, h);
}

static inline __m256 atmturb_bicubic_tap_avx2_os2(
    const float *row,
    long         sx,
    __m256i      vidx,
    __m256       vwx0,
    __m256       vwx1,
    __m256       vwx2,
    __m256       vwx3)
{
    __m256 v0 = _mm256_i32gather_ps(&row[sx - 1], vidx, 4);
    __m256 v1 = _mm256_i32gather_ps(&row[sx], vidx, 4);
    __m256 v2 = _mm256_i32gather_ps(&row[sx + 1], vidx, 4);
    __m256 v3 = _mm256_i32gather_ps(&row[sx + 2], vidx, 4);
    __m256 h = _mm256_mul_ps(v0, vwx0);
    h = _mm256_fmadd_ps(v1, vwx1, h);
    h = _mm256_fmadd_ps(v2, vwx2, h);
    return _mm256_fmadd_ps(v3, vwx3, h);
}

static inline void atmturb_bilinear_row_avx2_os1(
    const float *row0,
    const float *row1,
    float       *out_row,
    long         start_x,
    long         pup_size,
    __m256       vw00,
    __m256       vw10,
    __m256       vw01,
    __m256       vw11)
{
    long ii = 0;
    for (; ii <= pup_size - 16; ii += 16)
    {
        __m256 v00_0 = _mm256_loadu_ps(&row0[start_x + ii]);
        __m256 v00_1 = _mm256_loadu_ps(&row0[start_x + ii + 8]);
        __m256 v10_0 = _mm256_loadu_ps(&row0[start_x + ii + 1]);
        __m256 v10_1 = _mm256_loadu_ps(&row0[start_x + ii + 9]);

        __m256 v01_0 = _mm256_loadu_ps(&row1[start_x + ii]);
        __m256 v01_1 = _mm256_loadu_ps(&row1[start_x + ii + 8]);
        __m256 v11_0 = _mm256_loadu_ps(&row1[start_x + ii + 1]);
        __m256 v11_1 = _mm256_loadu_ps(&row1[start_x + ii + 9]);

        __m256 acc0 = _mm256_mul_ps(v00_0, vw00);
        __m256 acc1 = _mm256_mul_ps(v00_1, vw00);

        acc0 = _mm256_fmadd_ps(v10_0, vw10, acc0);
        acc1 = _mm256_fmadd_ps(v10_1, vw10, acc1);

        acc0 = _mm256_fmadd_ps(v01_0, vw01, acc0);
        acc1 = _mm256_fmadd_ps(v01_1, vw01, acc1);

        acc0 = _mm256_fmadd_ps(v11_0, vw11, acc0);
        acc1 = _mm256_fmadd_ps(v11_1, vw11, acc1);

        __m256 out0 = _mm256_loadu_ps(&out_row[ii]);
        __m256 out1 = _mm256_loadu_ps(&out_row[ii + 8]);

        _mm256_storeu_ps(&out_row[ii],     _mm256_add_ps(out0, acc0));
        _mm256_storeu_ps(&out_row[ii + 8], _mm256_add_ps(out1, acc1));
    }
    for (; ii <= pup_size - 8; ii += 8)
    {
        __m256 acc = atmturb_bilinear_tap_avx2_os1(row0, row1, start_x + ii,
                                                   vw00, vw10, vw01, vw11);
        __m256 out = _mm256_loadu_ps(&out_row[ii]);
        _mm256_storeu_ps(&out_row[ii], _mm256_add_ps(out, acc));
    }
    float w00 = _mm256_cvtss_f32(vw00), w10 = _mm256_cvtss_f32(vw10);
    float w01 = _mm256_cvtss_f32(vw01), w11 = _mm256_cvtss_f32(vw11);
    for (; ii < pup_size; ii++)
    {
        long ix0 = start_x + ii;
        out_row[ii] += w00 * row0[ix0] + w10 * row0[ix0 + 1] +
                       w01 * row1[ix0] + w11 * row1[ix0 + 1];
    }
}

/**
 * atmturb_extrude_accumulate_bilinear_avx2 - AVX2 accelerated bilinear extrusion
 * @params: Extrusion configuration and data pointers.
 */
void atmturb_extrude_accumulate_bilinear_avx2(
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

    __m256 vw00 = _mm256_set1_ps((1.0f - fx) * (1.0f - fy) * weight);
    __m256 vw10 = _mm256_set1_ps(fx * (1.0f - fy) * weight);
    __m256 vw01 = _mm256_set1_ps((1.0f - fx) * fy * weight);
    __m256 vw11 = _mm256_set1_ps(fx * fy * weight);

    static const int s_idx[8] = {0, 2, 4, 6, 8, 10, 12, 14};
    __m256i vidx = _mm256_loadu_si256((const __m256i *) s_idx);

    long start_x = base_x % msize;
    if (start_x < 0)
    {
        start_x += msize;
    }
    long iy0 = base_y % msize;
    if (iy0 < 0)
    {
        iy0 += msize;
    }

    for (long jj = 0; jj < pup_size; jj++)
    {
        long iy1 = (iy0 + 1 == msize) ? 0 : (iy0 + 1);

        const float *row0 = &master[iy0 * msize];
        const float *row1 = &master[iy1 * msize];
        float *out_row = &out_pha[jj * pup_size];

        if (start_x + (pup_size - 1) * os + 1 < msize)
        {
            if (os == 1)
            {
                atmturb_bilinear_row_avx2_os1(row0, row1, out_row, start_x, pup_size,
                                             vw00, vw10, vw01, vw11);
            }
            else if (os == 2)
            {
                long ii = 0;
                for (; ii <= pup_size - 8; ii += 8)
                {
                    __m256 acc = atmturb_bilinear_tap_avx2_os2(row0, row1, start_x + 2 * ii,
                                                               vidx, vw00, vw10, vw01, vw11);
                    __m256 out = _mm256_loadu_ps(&out_row[ii]);
                    _mm256_storeu_ps(&out_row[ii], _mm256_add_ps(out, acc));
                }
                float w00 = _mm256_cvtss_f32(vw00), w10 = _mm256_cvtss_f32(vw10);
                float w01 = _mm256_cvtss_f32(vw01), w11 = _mm256_cvtss_f32(vw11);
                for (; ii < pup_size; ii++)
                {
                    long ix0 = start_x + ii * os;
                    out_row[ii] += w00 * row0[ix0] + w10 * row0[ix0 + 1] +
                                   w01 * row1[ix0] + w11 * row1[ix0 + 1];
                }
            }
        }
        else
        {
            float w00 = _mm256_cvtss_f32(vw00), w10 = _mm256_cvtss_f32(vw10);
            float w01 = _mm256_cvtss_f32(vw01), w11 = _mm256_cvtss_f32(vw11);
            for (long ii = 0; ii < pup_size; ii++)
            {
                long ix0 = (start_x + ii * os) % msize;
                if (ix0 < 0) ix0 += msize;
                long ix1 = (ix0 + 1) % msize;
                out_row[ii] += w00 * row0[ix0] + w10 * row0[ix1] +
                               w01 * row1[ix0] + w11 * row1[ix1];
            }
        }

        iy0 += os;
        if (iy0 >= msize)
        {
            iy0 -= msize;
        }
    }
}

/**
 * atmturb_extrude_accumulate_bicubic_avx2 - AVX2 accelerated Keys bicubic extrusion
 * @params: Extrusion configuration and data pointers.
 */
void atmturb_extrude_accumulate_bicubic_avx2(
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

    __m256 vwx0 = _mm256_set1_ps(wx[0]), vwx1 = _mm256_set1_ps(wx[1]);
    __m256 vwx2 = _mm256_set1_ps(wx[2]), vwx3 = _mm256_set1_ps(wx[3]);
    __m256 vwy0 = _mm256_set1_ps(wy[0] * weight), vwy1 = _mm256_set1_ps(wy[1] * weight);
    __m256 vwy2 = _mm256_set1_ps(wy[2] * weight), vwy3 = _mm256_set1_ps(wy[3] * weight);

    static const int s_idx[8] = {0, 2, 4, 6, 8, 10, 12, 14};
    __m256i vidx = _mm256_loadu_si256((const __m256i *) s_idx);

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
                for (; ii <= pup_size - 8; ii += 8)
                {
                    __m256 h0 = atmturb_bicubic_tap_avx2_os1(row0, start_x + ii,
                                                             vwx0, vwx1, vwx2, vwx3);
                    __m256 h1 = atmturb_bicubic_tap_avx2_os1(row1, start_x + ii,
                                                             vwx0, vwx1, vwx2, vwx3);
                    __m256 h2 = atmturb_bicubic_tap_avx2_os1(row2, start_x + ii,
                                                             vwx0, vwx1, vwx2, vwx3);
                    __m256 h3 = atmturb_bicubic_tap_avx2_os1(row3, start_x + ii,
                                                             vwx0, vwx1, vwx2, vwx3);
                    __m256 acc = _mm256_mul_ps(h0, vwy0);
                    acc = _mm256_fmadd_ps(h1, vwy1, acc);
                    acc = _mm256_fmadd_ps(h2, vwy2, acc);
                    acc = _mm256_fmadd_ps(h3, vwy3, acc);
                    __m256 out = _mm256_loadu_ps(&out_row[ii]);
                    _mm256_storeu_ps(&out_row[ii], _mm256_add_ps(out, acc));
                }
            }
            else if (os == 2)
            {
                for (; ii <= pup_size - 8; ii += 8)
                {
                    __m256 h0 = atmturb_bicubic_tap_avx2_os2(row0, start_x + 2 * ii, vidx,
                                                             vwx0, vwx1, vwx2, vwx3);
                    __m256 h1 = atmturb_bicubic_tap_avx2_os2(row1, start_x + 2 * ii, vidx,
                                                             vwx0, vwx1, vwx2, vwx3);
                    __m256 h2 = atmturb_bicubic_tap_avx2_os2(row2, start_x + 2 * ii, vidx,
                                                             vwx0, vwx1, vwx2, vwx3);
                    __m256 h3 = atmturb_bicubic_tap_avx2_os2(row3, start_x + 2 * ii, vidx,
                                                             vwx0, vwx1, vwx2, vwx3);
                    __m256 acc = _mm256_mul_ps(h0, vwy0);
                    acc = _mm256_fmadd_ps(h1, vwy1, acc);
                    acc = _mm256_fmadd_ps(h2, vwy2, acc);
                    acc = _mm256_fmadd_ps(h3, vwy3, acc);
                    __m256 out = _mm256_loadu_ps(&out_row[ii]);
                    _mm256_storeu_ps(&out_row[ii], _mm256_add_ps(out, acc));
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
 * atmturb_extrude_accumulate_avx2 - AVX2 extrusion dispatcher (bilinear/bicubic)
 * @params: Extrusion configuration and data pointers.
 */
void atmturb_extrude_accumulate_avx2(
    const atmturb_extrude_params_t *params)
{
    if (params->interp == ATMTURB_INTERP_BICUBIC)
    {
        atmturb_extrude_accumulate_bicubic_avx2(params);
    }
    else
    {
        atmturb_extrude_accumulate_bilinear_avx2(params);
    }
}

/**
 * atmturb_scale_float_array_avx2 - Vectorized array multiplication by scalar (AVX2)
 * @dest: Output float array.
 * @src: Input float array.
 * @scale: Scalar multiplier.
 * @n: Number of elements.
 */
void atmturb_scale_float_array_avx2(
    float       *dest,
    const float *src,
    float        scale,
    long         n)
{
    __m256 vscale = _mm256_set1_ps(scale);
    long i = 0;
    for (; i <= n - 32; i += 32)
    {
        __m256 v0 = _mm256_loadu_ps(&src[i]);
        __m256 v1 = _mm256_loadu_ps(&src[i + 8]);
        __m256 v2 = _mm256_loadu_ps(&src[i + 16]);
        __m256 v3 = _mm256_loadu_ps(&src[i + 24]);
        _mm256_storeu_ps(&dest[i],      _mm256_mul_ps(v0, vscale));
        _mm256_storeu_ps(&dest[i + 8],  _mm256_mul_ps(v1, vscale));
        _mm256_storeu_ps(&dest[i + 16], _mm256_mul_ps(v2, vscale));
        _mm256_storeu_ps(&dest[i + 24], _mm256_mul_ps(v3, vscale));
    }
    for (; i <= n - 8; i += 8)
    {
        __m256 v = _mm256_loadu_ps(&src[i]);
        _mm256_storeu_ps(&dest[i], _mm256_mul_ps(v, vscale));
    }
    for (; i < n; i++)
    {
        dest[i] = src[i] * scale;
    }
}

/**
 * atmturb_init_phase_amp_avx2 - Vectorized phase zeroing and amplitude one-filling (AVX2)
 * @pha: Phase array.
 * @amp: Amplitude array.
 * @n: Number of elements.
 */
void atmturb_init_phase_amp_avx2(
    float *pha,
    float *amp,
    long   n)
{
    if (pha != NULL)
    {
        memset(pha, 0, sizeof(float) * (size_t) n);
    }

    if (amp == NULL)
    {
        return;
    }

    __m256 vone = _mm256_set1_ps(1.0f);
    long i = 0;
    for (; i <= n - 32; i += 32)
    {
        _mm256_storeu_ps(&amp[i],      vone);
        _mm256_storeu_ps(&amp[i + 8],  vone);
        _mm256_storeu_ps(&amp[i + 16], vone);
        _mm256_storeu_ps(&amp[i + 24], vone);
    }
    for (; i <= n - 8; i += 8)
    {
        _mm256_storeu_ps(&amp[i], vone);
    }
    for (; i < n; i++)
    {
        amp[i] = 1.0f;
    }
}

/**
 * atmturb_add_float_array_avx2 - Vectorized array accumulation (dest[i] += src[i])
 * @dest: Output/accumulator float array.
 * @src: Input float array.
 * @n: Number of elements.
 */
void atmturb_add_float_array_avx2(
    float       *dest,
    const float *src,
    long         n)
{
    long i = 0;
    for (; i <= n - 32; i += 32)
    {
        __m256 d0 = _mm256_loadu_ps(&dest[i]);
        __m256 d1 = _mm256_loadu_ps(&dest[i + 8]);
        __m256 d2 = _mm256_loadu_ps(&dest[i + 16]);
        __m256 d3 = _mm256_loadu_ps(&dest[i + 24]);

        __m256 s0 = _mm256_loadu_ps(&src[i]);
        __m256 s1 = _mm256_loadu_ps(&src[i + 8]);
        __m256 s2 = _mm256_loadu_ps(&src[i + 16]);
        __m256 s3 = _mm256_loadu_ps(&src[i + 24]);

        _mm256_storeu_ps(&dest[i],      _mm256_add_ps(d0, s0));
        _mm256_storeu_ps(&dest[i + 8],  _mm256_add_ps(d1, s1));
        _mm256_storeu_ps(&dest[i + 16], _mm256_add_ps(d2, s2));
        _mm256_storeu_ps(&dest[i + 24], _mm256_add_ps(d3, s3));
    }
    for (; i <= n - 8; i += 8)
    {
        __m256 d = _mm256_loadu_ps(&dest[i]);
        __m256 s = _mm256_loadu_ps(&src[i]);
        _mm256_storeu_ps(&dest[i], _mm256_add_ps(d, s));
    }
    for (; i < n; i++)
    {
        dest[i] += src[i];
    }
}

/**
 * atmturb_complex_mul_array_avx2 - Vectorized complex array point-wise product
 * @dest: Output complex float array (length 2 * n_complex).
 * @src1: First input complex float array.
 * @src2: Second input complex float array.
 * @n_complex: Number of complex elements.
 */
void atmturb_complex_mul_array_avx2(
    float       *dest,
    const float *src1,
    const float *src2,
    long         n_complex)
{
    long i = 0;
    for (; i <= n_complex - 8; i += 8)
    {
        long idx = 2 * i;
        __m256 va0 = _mm256_loadu_ps(&src1[idx]);
        __m256 va1 = _mm256_loadu_ps(&src1[idx + 8]);
        __m256 vb0 = _mm256_loadu_ps(&src2[idx]);
        __m256 vb1 = _mm256_loadu_ps(&src2[idx + 8]);

        __m256 a0_re = _mm256_moveldup_ps(va0);
        __m256 a0_im = _mm256_movehdup_ps(va0);
        __m256 a1_re = _mm256_moveldup_ps(va1);
        __m256 a1_im = _mm256_movehdup_ps(va1);

        __m256 b0_sw = _mm256_permute_ps(vb0, _MM_SHUFFLE(2, 3, 0, 1));
        __m256 b1_sw = _mm256_permute_ps(vb1, _MM_SHUFFLE(2, 3, 0, 1));

        __m256 r0 = _mm256_addsub_ps(_mm256_mul_ps(a0_re, vb0),
                                     _mm256_mul_ps(a0_im, b0_sw));
        __m256 r1 = _mm256_addsub_ps(_mm256_mul_ps(a1_re, vb1),
                                     _mm256_mul_ps(a1_im, b1_sw));

        _mm256_storeu_ps(&dest[idx],     r0);
        _mm256_storeu_ps(&dest[idx + 8], r1);
    }
    for (; i < n_complex; i++)
    {
        long idx = 2 * i;
        float r1 = src1[idx];
        float i1 = src1[idx + 1];
        float r2 = src2[idx];
        float i2 = src2[idx + 1];
        dest[idx]     = r1 * r2 - i1 * i2;
        dest[idx + 1] = r1 * i2 + i1 * r2;
    }
}

/**
 * atmturb_lowfreq_mode_accumulate_avx2 - AVX2 separable mode accumulation into pupil
 * @params: Low-frequency configuration bundle.
 * @kx: Mode spatial frequency along X [rad/master px].
 * @ky: Mode spatial frequency along Y [rad/master px].
 * @amp_re: Mode real amplitude.
 * @amp_im: Mode imaginary amplitude.
 */
static void atmturb_lowfreq_mode_accumulate_avx2(
    const atmturb_lowfreq_params_t *params,
    float                           kx,
    float                           ky,
    float                           amp_re,
    float                           amp_im)
{
    float h_re[1024], h_im[1024];
    long n = (params->pup_size <= 1024) ? params->pup_size : 1024;
    for (long i = 0; i < n; i++)
    {
        float th_x = kx * ((float) params->x0 + (float) (i * params->os));
        sincosf(th_x, &h_im[i], &h_re[i]);
    }

    for (long j = 0; j < params->pup_size; j++)
    {
        float th_y = ky * ((float) params->y0 + (float) (j * params->os));
        float vy_re, vy_im;
        sincosf(th_y, &vy_im, &vy_re);
        float c_re = (amp_re * vy_re - amp_im * vy_im) * params->weight;
        float c_im = (amp_re * vy_im + amp_im * vy_re) * params->weight;

        __m256 v_cre = _mm256_set1_ps(c_re);
        __m256 v_cim = _mm256_set1_ps(c_im);
        long row = j * params->pup_size;
        long i = 0;
        for (; i + 16 <= n; i += 16)
        {
            __m256 v_out0 = _mm256_loadu_ps(&params->out_pha[row + i]);
            __m256 v_out1 = _mm256_loadu_ps(&params->out_pha[row + i + 8]);
            __m256 v_hre0 = _mm256_loadu_ps(&h_re[i]);
            __m256 v_hre1 = _mm256_loadu_ps(&h_re[i + 8]);
            __m256 v_him0 = _mm256_loadu_ps(&h_im[i]);
            __m256 v_him1 = _mm256_loadu_ps(&h_im[i + 8]);
            v_out0 = _mm256_fmadd_ps(v_cre, v_hre0, v_out0);
            v_out1 = _mm256_fmadd_ps(v_cre, v_hre1, v_out1);
            v_out0 = _mm256_fnmadd_ps(v_cim, v_him0, v_out0);
            v_out1 = _mm256_fnmadd_ps(v_cim, v_him1, v_out1);
            _mm256_storeu_ps(&params->out_pha[row + i],     v_out0);
            _mm256_storeu_ps(&params->out_pha[row + i + 8], v_out1);
        }
        for (; i + 8 <= n; i += 8)
        {
            __m256 v_out = _mm256_loadu_ps(&params->out_pha[row + i]);
            __m256 v_hre = _mm256_loadu_ps(&h_re[i]);
            __m256 v_him = _mm256_loadu_ps(&h_im[i]);
            v_out = _mm256_fmadd_ps(v_cre, v_hre, v_out);
            v_out = _mm256_fnmadd_ps(v_cim, v_him, v_out);
            _mm256_storeu_ps(&params->out_pha[row + i], v_out);
        }
        for (; i < n; i++)
        {
            params->out_pha[row + i] += c_re * h_re[i] - c_im * h_im[i];
        }
    }
}

/**
 * atmturb_lowfreq_batch4_accumulate_avx2 - AVX2 separable mode accumulation for 4 modes
 * @params: Low-frequency configuration bundle.
 * @kx: Mode spatial frequencies along X (4 elements).
 * @ky: Mode spatial frequencies along Y (4 elements).
 * @amp_re: Mode real amplitudes (4 elements).
 * @amp_im: Mode imaginary amplitudes (4 elements).
 */
static void atmturb_lowfreq_batch4_accumulate_avx2(
    const atmturb_lowfreq_params_t *params,
    const float                    *kx,
    const float                    *ky,
    const float                    *amp_re,
    const float                    *amp_im)
{
    float h_re[4][1024] __attribute__((aligned(32)));
    float h_im[4][1024] __attribute__((aligned(32)));
    long n = (params->pup_size <= 1024) ? params->pup_size : 1024;

    for (int m = 0; m < 4; m++)
    {
        for (long i = 0; i < n; i++)
        {
            float th_x = kx[m] * ((float) params->x0 + (float) (i * params->os));
            sincosf(th_x, &h_im[m][i], &h_re[m][i]);
        }
    }

    for (long j = 0; j < params->pup_size; j++)
    {
        float c_re[4], c_im[4];
        for (int m = 0; m < 4; m++)
        {
            float th_y = ky[m] * ((float) params->y0 + (float) (j * params->os));
            float vy_re, vy_im;
            sincosf(th_y, &vy_im, &vy_re);
            c_re[m] = (amp_re[m] * vy_re - amp_im[m] * vy_im) * params->weight;
            c_im[m] = (amp_re[m] * vy_im + amp_im[m] * vy_re) * params->weight;
        }

        __m256 v_cre0 = _mm256_set1_ps(c_re[0]);
        __m256 v_cim0 = _mm256_set1_ps(c_im[0]);
        __m256 v_cre1 = _mm256_set1_ps(c_re[1]);
        __m256 v_cim1 = _mm256_set1_ps(c_im[1]);
        __m256 v_cre2 = _mm256_set1_ps(c_re[2]);
        __m256 v_cim2 = _mm256_set1_ps(c_im[2]);
        __m256 v_cre3 = _mm256_set1_ps(c_re[3]);
        __m256 v_cim3 = _mm256_set1_ps(c_im[3]);

        long row = j * params->pup_size;
        long i = 0;
        for (; i + 8 <= n; i += 8)
        {
            __m256 v_out = _mm256_loadu_ps(&params->out_pha[row + i]);

            v_out = _mm256_fmadd_ps(v_cre0, _mm256_load_ps(&h_re[0][i]), v_out);
            v_out = _mm256_fnmadd_ps(v_cim0, _mm256_load_ps(&h_im[0][i]), v_out);

            v_out = _mm256_fmadd_ps(v_cre1, _mm256_load_ps(&h_re[1][i]), v_out);
            v_out = _mm256_fnmadd_ps(v_cim1, _mm256_load_ps(&h_im[1][i]), v_out);

            v_out = _mm256_fmadd_ps(v_cre2, _mm256_load_ps(&h_re[2][i]), v_out);
            v_out = _mm256_fnmadd_ps(v_cim2, _mm256_load_ps(&h_im[2][i]), v_out);

            v_out = _mm256_fmadd_ps(v_cre3, _mm256_load_ps(&h_re[3][i]), v_out);
            v_out = _mm256_fnmadd_ps(v_cim3, _mm256_load_ps(&h_im[3][i]), v_out);

            _mm256_storeu_ps(&params->out_pha[row + i], v_out);
        }
        for (; i < n; i++)
        {
            float d = (c_re[0] * h_re[0][i] - c_im[0] * h_im[0][i]) +
                      (c_re[1] * h_re[1][i] - c_im[1] * h_im[1][i]) +
                      (c_re[2] * h_re[2][i] - c_im[2] * h_im[2][i]) +
                      (c_re[3] * h_re[3][i] - c_im[3] * h_im[3][i]);
            params->out_pha[row + i] += d;
        }
    }
}

/**
 * atmturb_extrude_lowfreq_avx2 - AVX2 separable low-order mode accumulation
 * @params: Low-frequency configuration and data pointers.
 */
void atmturb_extrude_lowfreq_avx2(
    const atmturb_lowfreq_params_t *params)
{
    const atmturb_lowfreq_t *lf = (const atmturb_lowfreq_t *) params->lf;
    const float *are = params->custom_are;
    const float *aim = params->custom_aim;
    if (are == NULL || aim == NULL)
    {
        int is_s1 = (params->screen_idx == 1);
        are = is_s1 ? lf->bre : lf->are;
        aim = is_s1 ? lf->bim : lf->aim;
    }

    int m = 0;
    for (; m + 4 <= lf->nmodes; m += 4)
    {
        atmturb_lowfreq_batch4_accumulate_avx2(params, &lf->kx[m], &lf->ky[m],
                                               &are[m], &aim[m]);
    }
    for (; m < lf->nmodes; m++)
    {
        atmturb_lowfreq_mode_accumulate_avx2(params, lf->kx[m], lf->ky[m], are[m], aim[m]);
    }
}

/**
 * atmturb_remove_piston_stream_avx2 - AVX2 vectorized piston removal with stream copy
 * @dst: Destination phase array.
 * @src: Source phase array.
 * @npix: Total number of pixels.
 */
void atmturb_remove_piston_stream_avx2(
    float       *restrict dst,
    const float *restrict src,
    long                  npix)
{
    if (dst == NULL || src == NULL || npix <= 0)
    {
        return;
    }

    __m256 vsum0 = _mm256_setzero_ps();
    __m256 vsum1 = _mm256_setzero_ps();
    __m256 vsum2 = _mm256_setzero_ps();
    __m256 vsum3 = _mm256_setzero_ps();

    long i = 0;
    for (; i <= npix - 32; i += 32)
    {
        vsum0 = _mm256_add_ps(vsum0, _mm256_loadu_ps(&src[i]));
        vsum1 = _mm256_add_ps(vsum1, _mm256_loadu_ps(&src[i + 8]));
        vsum2 = _mm256_add_ps(vsum2, _mm256_loadu_ps(&src[i + 16]));
        vsum3 = _mm256_add_ps(vsum3, _mm256_loadu_ps(&src[i + 24]));
    }
    __m256 vsum = _mm256_add_ps(_mm256_add_ps(vsum0, vsum1), _mm256_add_ps(vsum2, vsum3));
    for (; i <= npix - 8; i += 8)
    {
        vsum = _mm256_add_ps(vsum, _mm256_loadu_ps(&src[i]));
    }
    __m128 vlow    = _mm256_castps256_ps128(vsum);
    __m128 vhigh   = _mm256_extractf128_ps(vsum, 1);
    __m128 vsum128 = _mm_add_ps(vlow, vhigh);
    vsum128 = _mm_hadd_ps(vsum128, vsum128);
    vsum128 = _mm_hadd_ps(vsum128, vsum128);
    double sum = (double) _mm_cvtss_f32(vsum128);
    for (; i < npix; i++)
    {
        sum += (double) src[i];
    }

    float mean = (float) (sum / (double) npix);
    __m256 vmean = _mm256_set1_ps(mean);

    i = 0;
    for (; i <= npix - 32; i += 32)
    {
        __m256 s0 = _mm256_sub_ps(_mm256_loadu_ps(&src[i]), vmean);
        __m256 s1 = _mm256_sub_ps(_mm256_loadu_ps(&src[i + 8]), vmean);
        __m256 s2 = _mm256_sub_ps(_mm256_loadu_ps(&src[i + 16]), vmean);
        __m256 s3 = _mm256_sub_ps(_mm256_loadu_ps(&src[i + 24]), vmean);
        _mm256_storeu_ps(&dst[i], s0);
        _mm256_storeu_ps(&dst[i + 8], s1);
        _mm256_storeu_ps(&dst[i + 16], s2);
        _mm256_storeu_ps(&dst[i + 24], s3);
    }
    for (; i <= npix - 8; i += 8)
    {
        _mm256_storeu_ps(&dst[i], _mm256_sub_ps(_mm256_loadu_ps(&src[i]), vmean));
    }
    for (; i < npix; i++)
    {
        dst[i] = src[i] - mean;
    }
}

#else

void atmturb_extrude_accumulate_avx2(
    const atmturb_extrude_params_t *params)
{
    atmturb_extrude_accumulate_scalar(params);
}

void atmturb_extrude_accumulate_bilinear_avx2(
    const atmturb_extrude_params_t *params)
{
    atmturb_extrude_accumulate_bilinear_scalar(params);
}

void atmturb_extrude_accumulate_bicubic_avx2(
    const atmturb_extrude_params_t *params)
{
    atmturb_extrude_accumulate_bicubic_scalar(params);
}

void atmturb_scale_float_array_avx2(
    float       *dest,
    const float *src,
    float        scale,
    long         n)
{
    atmturb_scale_float_array_scalar(dest, src, scale, n);
}

void atmturb_init_phase_amp_avx2(
    float *pha,
    float *amp,
    long   n)
{
    atmturb_init_phase_amp_scalar(pha, amp, n);
}

void atmturb_extrude_lowfreq_avx2(
    const atmturb_lowfreq_params_t *params)
{
    atmturb_extrude_lowfreq_scalar(params);
}

void atmturb_add_float_array_avx2(
    float       *dest,
    const float *src,
    long         n)
{
    atmturb_add_float_array_scalar(dest, src, n);
}

void atmturb_complex_mul_array_avx2(
    float       *dest,
    const float *src1,
    const float *src2,
    long         n_complex)
{
    atmturb_complex_mul_array_scalar(dest, src1, src2, n_complex);
}

void atmturb_remove_piston_stream_avx2(
    float       *restrict dst,
    const float *restrict src,
    long                  npix)
{
    atmturb_remove_piston_stream_scalar(dst, src, npix);
}

#endif
