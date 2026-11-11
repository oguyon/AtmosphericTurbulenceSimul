// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_simd_avx2_ops.c
 * @brief   AVX2-accelerated 1D array, low-frequency subharmonic, and piston removal kernels
 *
 * Implements 256-bit AVX2 + FMA vectorization for pupil buffer initialization, scalar scaling,
 * pointwise array accumulation, interleaved complex multiplication, separable low-frequency
 * subharmonic synthesis (batched 4-mode and single-mode passes), and two-pass piston subtraction.
 * Falls back cleanly to portable scalar implementations when compiled without `-mavx2`.
 */

#define _GNU_SOURCE
#include <math.h>
#include <stddef.h>
#include <string.h>
#ifdef _OPENMP
#    include <omp.h>
#endif
#if defined(__x86_64__) || defined(_M_X64)
#    include <immintrin.h>
#endif
#include "atmturb_lowfreq.h"
#include "atmturb_simd.h"

#if defined(__AVX2__)

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
        _mm256_storeu_ps(&dest[i],      _mm256_mul_ps(_mm256_loadu_ps(&src[i]), vscale));
        _mm256_storeu_ps(&dest[i + 8],  _mm256_mul_ps(_mm256_loadu_ps(&src[i + 8]), vscale));
        _mm256_storeu_ps(&dest[i + 16], _mm256_mul_ps(_mm256_loadu_ps(&src[i + 16]), vscale));
        _mm256_storeu_ps(&dest[i + 24], _mm256_mul_ps(_mm256_loadu_ps(&src[i + 24]), vscale));
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
 * @pha: Phase array to zero (or NULL to skip).
 * @amp: Amplitude array to fill with 1.0f (or NULL to skip).
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
        __m256 d0 = _mm256_loadu_ps(&dest[i]), s0 = _mm256_loadu_ps(&src[i]);
        __m256 d1 = _mm256_loadu_ps(&dest[i + 8]), s1 = _mm256_loadu_ps(&src[i + 8]);
        __m256 d2 = _mm256_loadu_ps(&dest[i + 16]), s2 = _mm256_loadu_ps(&src[i + 16]);
        __m256 d3 = _mm256_loadu_ps(&dest[i + 24]), s3 = _mm256_loadu_ps(&src[i + 24]);

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
 * @dest: Output interleaved complex float array (length 2 * n_complex).
 * @src1: First input interleaved complex float array.
 * @src2: Second input interleaved complex float array.
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

        __m256 r0 = _mm256_fmaddsub_ps(a0_re, vb0, _mm256_mul_ps(a0_im, b0_sw));
        __m256 r1 = _mm256_fmaddsub_ps(a1_re, vb1, _mm256_mul_ps(a1_im, b1_sw));

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

    for (long col_start = 0; col_start < params->pup_size; col_start += 1024)
    {
        long n = params->pup_size - col_start;
        if (n > 1024)
        {
            n = 1024;
        }

        for (long i = 0; i < n; i++)
        {
            long col = col_start + i;
            float th_x = kx * ((float) params->x0 + (float) (col * params->os));
            sincosf(th_x, &h_im[i], &h_re[i]);
        }

#ifdef _OPENMP
        #pragma omp parallel for if(!omp_in_parallel() && params->pup_size >= 128) schedule(static)
#endif
        for (long j = 0; j < params->pup_size; j++)
        {
            float th_y = ky * ((float) params->y0 + (float) (j * params->os));
            float vy_re, vy_im;
            sincosf(th_y, &vy_im, &vy_re);
            float c_re = (amp_re * vy_re - amp_im * vy_im) * params->weight;
            float c_im = (amp_re * vy_im + amp_im * vy_re) * params->weight;

            __m256 v_cre = _mm256_set1_ps(c_re);
            __m256 v_cim = _mm256_set1_ps(c_im);
            long row = j * params->pup_size + col_start;
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
}

/**
 * atmturb_lowfreq_batch4_accumulate_avx2 - AVX2 separable mode accumulation for 4 modes
 * @params: Low-frequency configuration bundle.
 * @kx: Mode spatial frequencies along X [rad/master px] (4 elements).
 * @ky: Mode spatial frequencies along Y [rad/master px] (4 elements).
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

    for (long col_start = 0; col_start < params->pup_size; col_start += 1024)
    {
        long n = params->pup_size - col_start;
        if (n > 1024)
        {
            n = 1024;
        }

        for (int m = 0; m < 4; m++)
        {
            for (long i = 0; i < n; i++)
            {
                long col = col_start + i;
                float th_x = kx[m] * ((float) params->x0 + (float) (col * params->os));
                sincosf(th_x, &h_im[m][i], &h_re[m][i]);
            }
        }

#ifdef _OPENMP
        #pragma omp parallel for if(!omp_in_parallel() && params->pup_size >= 128) schedule(static)
#endif
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

            long row = j * params->pup_size + col_start;
            long i = 0;
            for (; i + 16 <= n; i += 16)
            {
                __m256 v_out0 = _mm256_loadu_ps(&params->out_pha[row + i]);
                __m256 v_out1 = _mm256_loadu_ps(&params->out_pha[row + i + 8]);

                v_out0 = _mm256_fmadd_ps(v_cre0, _mm256_load_ps(&h_re[0][i]), v_out0);
                v_out1 = _mm256_fmadd_ps(v_cre0, _mm256_load_ps(&h_re[0][i + 8]), v_out1);
                v_out0 = _mm256_fnmadd_ps(v_cim0, _mm256_load_ps(&h_im[0][i]), v_out0);
                v_out1 = _mm256_fnmadd_ps(v_cim0, _mm256_load_ps(&h_im[0][i + 8]), v_out1);

                v_out0 = _mm256_fmadd_ps(v_cre1, _mm256_load_ps(&h_re[1][i]), v_out0);
                v_out1 = _mm256_fmadd_ps(v_cre1, _mm256_load_ps(&h_re[1][i + 8]), v_out1);
                v_out0 = _mm256_fnmadd_ps(v_cim1, _mm256_load_ps(&h_im[1][i]), v_out0);
                v_out1 = _mm256_fnmadd_ps(v_cim1, _mm256_load_ps(&h_im[1][i + 8]), v_out1);

                v_out0 = _mm256_fmadd_ps(v_cre2, _mm256_load_ps(&h_re[2][i]), v_out0);
                v_out1 = _mm256_fmadd_ps(v_cre2, _mm256_load_ps(&h_re[2][i + 8]), v_out1);
                v_out0 = _mm256_fnmadd_ps(v_cim2, _mm256_load_ps(&h_im[2][i]), v_out0);
                v_out1 = _mm256_fnmadd_ps(v_cim2, _mm256_load_ps(&h_im[2][i + 8]), v_out1);

                v_out0 = _mm256_fmadd_ps(v_cre3, _mm256_load_ps(&h_re[3][i]), v_out0);
                v_out1 = _mm256_fmadd_ps(v_cre3, _mm256_load_ps(&h_re[3][i + 8]), v_out1);
                v_out0 = _mm256_fnmadd_ps(v_cim3, _mm256_load_ps(&h_im[3][i]), v_out0);
                v_out1 = _mm256_fnmadd_ps(v_cim3, _mm256_load_ps(&h_im[3][i + 8]), v_out1);

                _mm256_storeu_ps(&params->out_pha[row + i], v_out0);
                _mm256_storeu_ps(&params->out_pha[row + i + 8], v_out1);
            }
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

    __m256 vsum0 = _mm256_setzero_ps(), vsum1 = _mm256_setzero_ps();
    __m256 vsum2 = _mm256_setzero_ps(), vsum3 = _mm256_setzero_ps();

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
        _mm256_storeu_ps(&dst[i],      _mm256_sub_ps(_mm256_loadu_ps(&src[i]), vmean));
        _mm256_storeu_ps(&dst[i + 8],  _mm256_sub_ps(_mm256_loadu_ps(&src[i + 8]), vmean));
        _mm256_storeu_ps(&dst[i + 16], _mm256_sub_ps(_mm256_loadu_ps(&src[i + 16]), vmean));
        _mm256_storeu_ps(&dst[i + 24], _mm256_sub_ps(_mm256_loadu_ps(&src[i + 24]), vmean));
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

#else /* !defined(__AVX2__) */

/*
 * Scalar fallback wrappers compiled when the translation unit is built without
 * AVX2 target support (e.g. non-x86_64 architectures or scalar-only builds).
 */

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

#endif /* defined(__AVX2__) */
