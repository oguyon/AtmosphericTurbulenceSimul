// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_simd_avx2.c
 * @brief   AVX2-accelerated compute kernels for atmospheric turbulence simulation
 */

#include <math.h>
#if defined(__x86_64__) || defined(_M_X64)
#    include <immintrin.h>
#endif
#include "atmturb_simd.h"

#if defined(__AVX2__)

/**
 * atmturb_extrude_row_avx2 - Vectorized bilinear interpolation along contiguous row
 * @row0: Master screen upper row.
 * @row1: Master screen lower row.
 * @start_x: Starting horizontal pixel offset.
 * @pup_size: Number of pixels in pupil row.
 * @vw00: Vector broadcast of upper-left weight.
 * @vw10: Vector broadcast of upper-right weight.
 * @vw01: Vector broadcast of lower-left weight.
 * @vw11: Vector broadcast of lower-right weight.
 * @out_row: Destination row to accumulate into.
 */
static inline void atmturb_extrude_row_avx2(const float *row0, const float *row1,
                                            long start_x, long pup_size,
                                            __m256 vw00, __m256 vw10,
                                            __m256 vw01, __m256 vw11,
                                            float *out_row)
{
    long ii = 0;
    for (; ii <= pup_size - 8; ii += 8)
    {
        __m256 v00 = _mm256_loadu_ps(&row0[start_x + ii]);
        __m256 v10 = _mm256_loadu_ps(&row0[start_x + ii + 1]);
        __m256 v01 = _mm256_loadu_ps(&row1[start_x + ii]);
        __m256 v11 = _mm256_loadu_ps(&row1[start_x + ii + 1]);

        __m256 acc = _mm256_mul_ps(v00, vw00);
        acc = _mm256_fmadd_ps(v10, vw10, acc);
        acc = _mm256_fmadd_ps(v01, vw01, acc);
        acc = _mm256_fmadd_ps(v11, vw11, acc);

        __m256 out = _mm256_loadu_ps(&out_row[ii]);
        _mm256_storeu_ps(&out_row[ii], _mm256_add_ps(out, acc));
    }
    float w00 = _mm256_cvtss_f32(vw00);
    float w10 = _mm256_cvtss_f32(vw10);
    float w01 = _mm256_cvtss_f32(vw01);
    float w11 = _mm256_cvtss_f32(vw11);
    for (; ii < pup_size; ii++)
    {
        out_row[ii] += w00 * row0[start_x + ii] + w10 * row0[start_x + ii + 1] +
                       w01 * row1[start_x + ii] + w11 * row1[start_x + ii + 1];
    }
}

/**
 * atmturb_extrude_accumulate_avx2 - AVX2 accelerated bilinear extrusion and accumulation
 * @master: Master screen array.
 * @msize: Master screen dimension.
 * @x0: Sub-pixel X coordinate offset.
 * @y0: Sub-pixel Y coordinate offset.
 * @pup_size: Pupil dimension.
 * @weight: Layer weight.
 * @out_pha: Output phase array.
 */
void atmturb_extrude_accumulate_avx2(const float *master, long msize, double x0, double y0,
                                     long pup_size, float weight, float *out_pha)
{
    float fx = (float)(x0 - floor(x0));
    float fy = (float)(y0 - floor(y0));

    long base_x = (long)floor(x0);
    long base_y = (long)floor(y0);

    float w00 = (1.0f - fx) * (1.0f - fy) * weight;
    float w10 = fx * (1.0f - fy) * weight;
    float w01 = (1.0f - fx) * fy * weight;
    float w11 = fx * fy * weight;

    __m256 vw00 = _mm256_set1_ps(w00);
    __m256 vw10 = _mm256_set1_ps(w10);
    __m256 vw01 = _mm256_set1_ps(w01);
    __m256 vw11 = _mm256_set1_ps(w11);

    long start_x = base_x % msize;
    if (start_x < 0) start_x += msize;

    for (long jj = 0; jj < pup_size; jj++)
    {
        long iy0 = (base_y + jj) % msize;
        if (iy0 < 0) iy0 += msize;
        long iy1 = (iy0 + 1) % msize;

        const float *row0 = &master[iy0 * msize];
        const float *row1 = &master[iy1 * msize];
        float *out_row = &out_pha[jj * pup_size];

        if (start_x + pup_size < msize)
        {
            atmturb_extrude_row_avx2(row0, row1, start_x, pup_size,
                                     vw00, vw10, vw01, vw11, out_row);
        }
        else
        {
            for (long ii = 0; ii < pup_size; ii++)
            {
                long ix0 = (start_x + ii) % msize;
                long ix1 = (ix0 + 1) % msize;
                out_row[ii] += w00 * row0[ix0] + w10 * row0[ix1] +
                               w01 * row1[ix0] + w11 * row1[ix1];
            }
        }
    }
}

/**
 * atmturb_scale_float_array_avx2 - Vectorized array multiplication by scalar (AVX2)
 * @dest: Output float array.
 * @src: Input float array.
 * @scale: Scalar multiplier.
 * @n: Number of elements.
 */
void atmturb_scale_float_array_avx2(float *dest, const float *src, float scale, long n)
{
    __m256 vscale = _mm256_set1_ps(scale);
    long i = 0;
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
void atmturb_init_phase_amp_avx2(float *pha, float *amp, long n)
{
    __m256 vzero = _mm256_setzero_ps();
    __m256 vone = _mm256_set1_ps(1.0f);
    long i = 0;
    for (; i <= n - 8; i += 8)
    {
        _mm256_storeu_ps(&pha[i], vzero);
        _mm256_storeu_ps(&amp[i], vone);
    }
    for (; i < n; i++)
    {
        pha[i] = 0.0f;
        amp[i] = 1.0f;
    }
}

#else

void atmturb_extrude_accumulate_avx2(const float *master, long msize, double x0, double y0,
                                     long pup_size, float weight, float *out_pha)
{
    atmturb_extrude_accumulate_scalar(master, msize, x0, y0, pup_size, weight, out_pha);
}

void atmturb_scale_float_array_avx2(float *dest, const float *src, float scale, long n)
{
    atmturb_scale_float_array_scalar(dest, src, scale, n);
}

void atmturb_init_phase_amp_avx2(float *pha, float *amp, long n)
{
    atmturb_init_phase_amp_scalar(pha, amp, n);
}

#endif
