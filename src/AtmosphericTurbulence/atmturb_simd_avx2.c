// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_simd_avx2.c
 * @brief   AVX2-accelerated 2D phase-screen extrusion kernels (bilinear and Keys bicubic)
 *
 * Implements 256-bit AVX2 + FMA sub-pixel 2D interpolation from toroidally wrapped master phase
 * screens onto pupil grids for oversampling strides `os == 1` (contiguous vector loads),
 * `os == 2` (`_mm256_i32gather_ps`), and general stride fallbacks. Array, low-frequency
 * subharmonic, and piston-removal AVX2 kernels reside in `atmturb_simd_avx2_ops.c`.
 */

#include <math.h>
#ifdef _OPENMP
#    include <omp.h>
#endif
#if defined(__x86_64__) || defined(_M_X64)
#    include <immintrin.h>
#endif
#include "atmturb_simd.h"

#if defined(__AVX2__)

/**
 * atmturb_keys_weights - Compute 4-tap Keys Catmull-Rom cubic interpolation weights (a = -0.5)
 * @f: Fractional sub-pixel offset in [0, 1).
 * @w: Output 4-element weight array for taps {-1, 0, +1, +2}.
 */
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

/**
 * atmturb_bilinear_tap_avx2_os1 - Evaluate 8 contiguous bilinear samples at unit stride (AVX2)
 * @row0: Upper master screen row pointer.
 * @row1: Lower master screen row pointer.
 * @sx: Starting column index in master screen.
 * @vw00: Broadcast top-left weight.
 * @vw10: Broadcast top-right weight.
 * @vw01: Broadcast bottom-left weight.
 * @vw11: Broadcast bottom-right weight.
 *
 * Return: 8-lane `__m256` vector of weighted bilinear phase values.
 */
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

/**
 * atmturb_bilinear_tap_avx2_os2 - Gather and evaluate 8 bilinear samples at stride 2 (AVX2)
 * @row0: Upper master screen row pointer.
 * @row1: Lower master screen row pointer.
 * @sx: Starting column index in master screen.
 * @vidx: Stride-2 gather index vector `{0, 2, 4, 6, 8, 10, 12, 14}`.
 * @vw00: Broadcast top-left weight.
 * @vw10: Broadcast top-right weight.
 * @vw01: Broadcast bottom-left weight.
 * @vw11: Broadcast bottom-right weight.
 *
 * Return: 8-lane `__m256` vector of weighted bilinear phase values.
 */
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

/**
 * atmturb_bicubic_tap_avx2_os1 - Evaluate 8 horizontal 4-tap cubic convolutions at stride 1 (AVX2)
 * @row: Master screen row pointer.
 * @sx: Starting anchor column index (`sx - 1` must be valid).
 * @vwx0: Broadcast horizontal cubic weight for tap -1.
 * @vwx1: Broadcast horizontal cubic weight for tap 0.
 * @vwx2: Broadcast horizontal cubic weight for tap +1.
 * @vwx3: Broadcast horizontal cubic weight for tap +2.
 *
 * Return: 8-lane `__m256` vector of horizontally interpolated row samples.
 */
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

/**
 * atmturb_bicubic_tap_avx2_os2 - Gather and evaluate 8 horizontal cubic taps at stride 2 (AVX2)
 * @row: Master screen row pointer.
 * @sx: Starting anchor column index (`sx - 1` must be valid).
 * @vidx: Stride-2 gather index vector `{0, 2, 4, 6, 8, 10, 12, 14}`.
 * @vwx0: Broadcast horizontal cubic weight for tap -1.
 * @vwx1: Broadcast horizontal cubic weight for tap 0.
 * @vwx2: Broadcast horizontal cubic weight for tap +1.
 * @vwx3: Broadcast horizontal cubic weight for tap +2.
 *
 * Return: 8-lane `__m256` vector of horizontally interpolated row samples.
 */
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

/**
 * atmturb_bilinear_row_avx2_os1 - Accumulate a non-wrapping stride-1 bilinear pupil row (AVX2)
 * @row0: Upper master screen row pointer.
 * @row1: Lower master screen row pointer.
 * @out_row: Output pupil row accumulator.
 * @start_x: Starting master column index (`start_x + pup_size < msize`).
 * @pup_size: Number of output pixels to accumulate.
 * @vw00: Broadcast top-left weight.
 * @vw10: Broadcast top-right weight.
 * @vw01: Broadcast bottom-left weight.
 * @vw11: Broadcast bottom-right weight.
 */
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
 * atmturb_bilinear_row_wrapped_avx2 - Accumulate a bilinear row wrapping across master seam
 * @row0: Upper master screen row pointer.
 * @row1: Lower master screen row pointer.
 * @out_row: Output pupil row accumulator.
 * @start_x: Wrapped starting column index in [0, msize).
 * @pup_size: Pupil row length in pixels.
 * @os: Oversampling stride in master pixels per pupil pixel.
 * @msize: Master screen side dimension in pixels.
 * @vw00: Broadcast top-left weight.
 * @vw10: Broadcast top-right weight.
 * @vw01: Broadcast bottom-left weight.
 * @vw11: Broadcast bottom-right weight.
 */
static inline void atmturb_bilinear_row_wrapped_avx2(
    const float *row0,
    const float *row1,
    float       *out_row,
    long         start_x,
    long         pup_size,
    long         os,
    long         msize,
    __m256       vw00,
    __m256       vw10,
    __m256       vw01,
    __m256       vw11)
{
    float w00 = _mm256_cvtss_f32(vw00);
    float w10 = _mm256_cvtss_f32(vw10);
    float w01 = _mm256_cvtss_f32(vw01);
    float w11 = _mm256_cvtss_f32(vw11);

    if (os == 1)
    {
        long seam = msize - 1 - start_x;
        if (seam >= 0 && seam < pup_size)
        {
            if (seam > 0)
            {
                atmturb_bilinear_row_avx2_os1(row0, row1, out_row, start_x, seam,
                                             vw00, vw10, vw01, vw11);
            }
            out_row[seam] += w00 * row0[msize - 1] + w10 * row0[0] +
                             w01 * row1[msize - 1] + w11 * row1[0];
            long rem = pup_size - (seam + 1);
            if (rem > 0)
            {
                atmturb_bilinear_row_avx2_os1(row0, row1, &out_row[seam + 1], 0, rem,
                                             vw00, vw10, vw01, vw11);
            }
            return;
        }
    }

    for (long ii = 0; ii < pup_size; ii++)
    {
        long ix0 = (start_x + ii * os) % msize;
        if (ix0 < 0)
        {
            ix0 += msize;
        }
        long ix1 = (ix0 + 1) % msize;
        out_row[ii] += w00 * row0[ix0] + w10 * row0[ix1] +
                       w01 * row1[ix0] + w11 * row1[ix1];
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

#ifdef _OPENMP
    #pragma omp parallel for if(!omp_in_parallel() && pup_size >= 128) schedule(static)
#endif
    for (long jj = 0; jj < pup_size; jj++)
    {
        long iy0 = (base_y + jj * os) % msize;
        if (iy0 < 0)
        {
            iy0 += msize;
        }
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
            else
            {
                long ii = 0;
                if (os == 2)
                {
                    for (; ii <= pup_size - 8; ii += 8)
                    {
                        __m256 acc = atmturb_bilinear_tap_avx2_os2(row0, row1, start_x + 2 * ii,
                                                                   vidx, vw00, vw10, vw01, vw11);
                        __m256 out = _mm256_loadu_ps(&out_row[ii]);
                        _mm256_storeu_ps(&out_row[ii], _mm256_add_ps(out, acc));
                    }
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
            atmturb_bilinear_row_wrapped_avx2(row0, row1, out_row, start_x, pup_size,
                                             os, msize, vw00, vw10, vw01, vw11);
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
    if (start_x < 0)
    {
        start_x += msize;
    }

#ifdef _OPENMP
    #pragma omp parallel for if(!omp_in_parallel() && pup_size >= 128) schedule(static)
#endif
    for (long jj = 0; jj < pup_size; jj++)
    {
        long by = base_y + jj * os;
        long iy[4];
        for (int n = 0; n < 4; n++)
        {
            long y = (by + n - 1) % msize;
            if (y < 0)
            {
                y += msize;
            }
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
                    __m256 acc = _mm256_fmadd_ps(h1, vwy1, _mm256_mul_ps(h0, vwy0));
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
                    __m256 acc = _mm256_fmadd_ps(h1, vwy1, _mm256_mul_ps(h0, vwy0));
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
                           wx[2] * row3[bx + 1] + wx[3] * row3[bx + 2];
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

#else /* !defined(__AVX2__) */

/*
 * Scalar fallback wrappers compiled when the translation unit is built without
 * AVX2 target support (e.g. non-x86_64 architectures or scalar-only builds).
 */

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

#endif /* defined(__AVX2__) */
