// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_simd_avx512.c
 * @brief   AVX-512 accelerated compute kernels for atmospheric turbulence simulation
 */

#define _GNU_SOURCE
#include <math.h>
#if defined(__x86_64__) || defined(_M_X64)
#    include <immintrin.h>
#endif
#include "atmturb_lowfreq.h"
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

    float fx = (float) (params->x0 - floor(params->x0));
    float fy = (float) (params->y0 - floor(params->y0));
    long base_x = (long) floor(params->x0);
    long base_y = (long) floor(params->y0);

    __m512 vw00 = _mm512_set1_ps((1.0f - fx) * (1.0f - fy) * weight);
    __m512 vw10 = _mm512_set1_ps(fx * (1.0f - fy) * weight);
    __m512 vw01 = _mm512_set1_ps((1.0f - fx) * fy * weight);
    __m512 vw11 = _mm512_set1_ps(fx * fy * weight);

    static const int s_idx[16] = {
        0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30
    };
    __m512i vidx = _mm512_loadu_si512((const __m512i *) s_idx);

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

    float w00 = (1.0f - fx) * (1.0f - fy) * weight;
    float w10 = fx * (1.0f - fy) * weight;
    float w01 = (1.0f - fx) * fy * weight;
    float w11 = fx * fy * weight;

    for (long jj = 0; jj < pup_size; jj++)
    {
        long iy1 = iy0 + 1;
        if (iy1 >= msize)
        {
            iy1 -= msize;
        }

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
            for (; ii < pup_size; ii++)
            {
                long ix0 = start_x + ii * os;
                out_row[ii] += w00 * row0[ix0] + w10 * row0[ix0 + 1] +
                               w01 * row1[ix0] + w11 * row1[ix0 + 1];
            }
        }
        else
        {
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
        iy0 += os;
        if (iy0 >= msize)
        {
            iy0 -= msize;
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

    float fx = (float) (params->x0 - floor(params->x0));
    float fy = (float) (params->y0 - floor(params->y0));
    long base_x = (long) floor(params->x0);
    long base_y = (long) floor(params->y0);

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
    if (start_x < 0)
    {
        start_x += msize;
    }
    long cur_y = base_y % msize;
    if (cur_y < 0)
    {
        cur_y += msize;
    }

    float wy0 = wy[0] * weight, wy1 = wy[1] * weight;
    float wy2 = wy[2] * weight, wy3 = wy[3] * weight;

    for (long jj = 0; jj < pup_size; jj++)
    {
        long iy0 = (cur_y > 0) ? (cur_y - 1) : (msize - 1);
        long iy1 = cur_y;
        long iy2 = (cur_y + 1 < msize) ? (cur_y + 1) : 0;
        long iy3 = (cur_y + 2 < msize) ? (cur_y + 2) : (cur_y + 2 - msize);

        const float *row0 = &master[iy0 * msize];
        const float *row1 = &master[iy1 * msize];
        const float *row2 = &master[iy2 * msize];
        const float *row3 = &master[iy3 * msize];
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
        cur_y += os;
        if (cur_y >= msize)
        {
            cur_y -= msize;
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
    for (; i <= n - 64; i += 64)
    {
        __m512 v0 = _mm512_loadu_ps(&src[i]);
        __m512 v1 = _mm512_loadu_ps(&src[i + 16]);
        __m512 v2 = _mm512_loadu_ps(&src[i + 32]);
        __m512 v3 = _mm512_loadu_ps(&src[i + 48]);

        _mm512_storeu_ps(&dest[i],      _mm512_mul_ps(v0, vscale));
        _mm512_storeu_ps(&dest[i + 16], _mm512_mul_ps(v1, vscale));
        _mm512_storeu_ps(&dest[i + 32], _mm512_mul_ps(v2, vscale));
        _mm512_storeu_ps(&dest[i + 48], _mm512_mul_ps(v3, vscale));
    }
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
    for (; i <= n - 32; i += 32)
    {
        _mm512_storeu_ps(&pha[i],      vzero);
        _mm512_storeu_ps(&pha[i + 16], vzero);
        _mm512_storeu_ps(&amp[i],      vone);
        _mm512_storeu_ps(&amp[i + 16], vone);
    }
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

/**
 * atmturb_lowfreq_mode_accumulate_avx512 - AVX-512 separable mode accumulation into pupil
 * @params: Low-frequency configuration bundle.
 * @kx: Mode spatial frequency along X [rad/master px].
 * @ky: Mode spatial frequency along Y [rad/master px].
 * @amp_re: Mode real amplitude.
 * @amp_im: Mode imaginary amplitude.
 */
static void atmturb_lowfreq_mode_accumulate_avx512(
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
        float th_x = kx * (float) (params->x0 + (double) (i * params->os));
        sincosf(th_x, &h_im[i], &h_re[i]);
    }

    for (long j = 0; j < params->pup_size; j++)
    {
        float th_y = ky * (float) (params->y0 + (double) (j * params->os));
        float vy_re, vy_im;
        sincosf(th_y, &vy_im, &vy_re);
        float c_re = (amp_re * vy_re - amp_im * vy_im) * params->weight;
        float c_im = (amp_re * vy_im + amp_im * vy_re) * params->weight;

        __m512 v_cre = _mm512_set1_ps(c_re);
        __m512 v_cim = _mm512_set1_ps(c_im);
        long row = j * params->pup_size;
        long i = 0;
        for (; i + 32 <= n; i += 32)
        {
            __m512 v_out0 = _mm512_loadu_ps(&params->out_pha[row + i]);
            __m512 v_out1 = _mm512_loadu_ps(&params->out_pha[row + i + 16]);
            __m512 v_hre0 = _mm512_loadu_ps(&h_re[i]);
            __m512 v_hre1 = _mm512_loadu_ps(&h_re[i + 16]);
            __m512 v_him0 = _mm512_loadu_ps(&h_im[i]);
            __m512 v_him1 = _mm512_loadu_ps(&h_im[i + 16]);

            v_out0 = _mm512_fmadd_ps(v_cre, v_hre0, v_out0);
            v_out1 = _mm512_fmadd_ps(v_cre, v_hre1, v_out1);
            v_out0 = _mm512_fnmadd_ps(v_cim, v_him0, v_out0);
            v_out1 = _mm512_fnmadd_ps(v_cim, v_him1, v_out1);

            _mm512_storeu_ps(&params->out_pha[row + i],      v_out0);
            _mm512_storeu_ps(&params->out_pha[row + i + 16], v_out1);
        }
        for (; i + 16 <= n; i += 16)
        {
            __m512 v_out = _mm512_loadu_ps(&params->out_pha[row + i]);
            __m512 v_hre = _mm512_loadu_ps(&h_re[i]);
            __m512 v_him = _mm512_loadu_ps(&h_im[i]);
            v_out = _mm512_fmadd_ps(v_cre, v_hre, v_out);
            v_out = _mm512_fnmadd_ps(v_cim, v_him, v_out);
            _mm512_storeu_ps(&params->out_pha[row + i], v_out);
        }
        for (; i < n; i++)
        {
            params->out_pha[row + i] += c_re * h_re[i] - c_im * h_im[i];
        }
    }
}

/**
 * atmturb_extrude_lowfreq_avx512 - AVX-512 separable low-order mode accumulation
 * @params: Low-frequency configuration and data pointers.
 */
void atmturb_extrude_lowfreq_avx512(
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

    for (int m = 0; m < lf->nmodes; m++)
    {
        atmturb_lowfreq_mode_accumulate_avx512(params, lf->kx[m], lf->ky[m], are[m], aim[m]);
    }
}

/**
 * atmturb_add_float_array_avx512 - Vectorized array accumulation (dest[i] += src[i])
 * @dest: Output/accumulator float array.
 * @src: Input float array.
 * @n: Number of elements.
 */
void atmturb_add_float_array_avx512(
    float       *dest,
    const float *src,
    long         n)
{
    long i = 0;
    for (; i <= n - 64; i += 64)
    {
        __m512 d0 = _mm512_loadu_ps(&dest[i]);
        __m512 d1 = _mm512_loadu_ps(&dest[i + 16]);
        __m512 d2 = _mm512_loadu_ps(&dest[i + 32]);
        __m512 d3 = _mm512_loadu_ps(&dest[i + 48]);

        __m512 s0 = _mm512_loadu_ps(&src[i]);
        __m512 s1 = _mm512_loadu_ps(&src[i + 16]);
        __m512 s2 = _mm512_loadu_ps(&src[i + 32]);
        __m512 s3 = _mm512_loadu_ps(&src[i + 48]);

        _mm512_storeu_ps(&dest[i],      _mm512_add_ps(d0, s0));
        _mm512_storeu_ps(&dest[i + 16], _mm512_add_ps(d1, s1));
        _mm512_storeu_ps(&dest[i + 32], _mm512_add_ps(d2, s2));
        _mm512_storeu_ps(&dest[i + 48], _mm512_add_ps(d3, s3));
    }
    for (; i <= n - 16; i += 16)
    {
        __m512 d = _mm512_loadu_ps(&dest[i]);
        __m512 s = _mm512_loadu_ps(&src[i]);
        _mm512_storeu_ps(&dest[i], _mm512_add_ps(d, s));
    }
    for (; i < n; i++)
    {
        dest[i] += src[i];
    }
}

/**
 * atmturb_complex_mul_array_avx512 - Vectorized complex array product (AVX-512)
 * @dest: Output complex float array (length 2 * n_complex).
 * @src1: First input complex float array.
 * @src2: Second input complex float array.
 * @n_complex: Number of complex elements.
 */
void atmturb_complex_mul_array_avx512(
    float       *dest,
    const float *src1,
    const float *src2,
    long         n_complex)
{
    long i = 0;
    for (; i <= n_complex - 16; i += 16)
    {
        long idx = 2 * i;
        __m512 va0 = _mm512_loadu_ps(&src1[idx]);
        __m512 va1 = _mm512_loadu_ps(&src1[idx + 16]);
        __m512 vb0 = _mm512_loadu_ps(&src2[idx]);
        __m512 vb1 = _mm512_loadu_ps(&src2[idx + 16]);

        __m512 a0_re = _mm512_moveldup_ps(va0);
        __m512 a0_im = _mm512_movehdup_ps(va0);
        __m512 a1_re = _mm512_moveldup_ps(va1);
        __m512 a1_im = _mm512_movehdup_ps(va1);

        __m512 b0_sw = _mm512_permute_ps(vb0, _MM_SHUFFLE(2, 3, 0, 1));
        __m512 b1_sw = _mm512_permute_ps(vb1, _MM_SHUFFLE(2, 3, 0, 1));

        __m512 p0_im = _mm512_mul_ps(a0_im, b0_sw);
        __m512 p1_im = _mm512_mul_ps(a1_im, b1_sw);

        __m512 r0 = _mm512_fmaddsub_ps(a0_re, vb0, p0_im);
        __m512 r1 = _mm512_fmaddsub_ps(a1_re, vb1, p1_im);

        _mm512_storeu_ps(&dest[idx],      r0);
        _mm512_storeu_ps(&dest[idx + 16], r1);
    }
    for (; i <= n_complex - 8; i += 8)
    {
        long idx = 2 * i;
        __m512 va = _mm512_loadu_ps(&src1[idx]);
        __m512 vb = _mm512_loadu_ps(&src2[idx]);
        __m512 a_re = _mm512_moveldup_ps(va);
        __m512 a_im = _mm512_movehdup_ps(va);
        __m512 b_sw = _mm512_permute_ps(vb, _MM_SHUFFLE(2, 3, 0, 1));
        __m512 p_im = _mm512_mul_ps(a_im, b_sw);
        __m512 res  = _mm512_fmaddsub_ps(a_re, vb, p_im);
        _mm512_storeu_ps(&dest[idx], res);
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

#else

void atmturb_extrude_accumulate_avx512(const atmturb_extrude_params_t *params)
{
    atmturb_extrude_accumulate_scalar(params);
}

void atmturb_extrude_accumulate_bilinear_avx512(const atmturb_extrude_params_t *params)
{
    atmturb_extrude_accumulate_bilinear_scalar(params);
}

void atmturb_extrude_accumulate_bicubic_avx512(const atmturb_extrude_params_t *params)
{
    atmturb_extrude_accumulate_bicubic_scalar(params);
}

void atmturb_scale_float_array_avx512(float *dest, const float *src, float scale, long n)
{
    atmturb_scale_float_array_scalar(dest, src, scale, n);
}

void atmturb_init_phase_amp_avx512(float *pha, float *amp, long n)
{
    atmturb_init_phase_amp_scalar(pha, amp, n);
}

void atmturb_extrude_lowfreq_avx512(const atmturb_lowfreq_params_t *params)
{
    atmturb_extrude_lowfreq_scalar(params);
}

void atmturb_add_float_array_avx512(float *dest, const float *src, long n)
{
    atmturb_add_float_array_scalar(dest, src, n);
}

void atmturb_complex_mul_array_avx512(
    float *dest, const float *src1, const float *src2, long n_complex)
{
    atmturb_complex_mul_array_scalar(dest, src1, src2, n_complex);
}

#endif
