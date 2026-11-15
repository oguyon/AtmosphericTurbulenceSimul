// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_psf_engine.c
 * @brief   Pruned 2-stage 1D and 2D FFTW3 Fourier optics engine for Nyquist science PSF synthesis
 */

#include <stddef.h>
#include <string.h>
#include <fftw3.h>

#include "atmturb_fft_wisdom.h"
#include "atmturb_fresnel_simd.h"
#include "atmturb_psf_engine.h"
#include "atmturb_screens_batch.h"

/**
 * atmturb_ao_psf_alloc_bufs - Allocate zero-padded FFT working arrays for PSF engine
 * @eng: PSF engine descriptor with dimensions populated.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_ao_psf_alloc_bufs(
    atmturb_ao_psf_engine_t *eng)
{
    size_t nc = (size_t) eng->n_psf * (size_t) eng->col_stride;
    eng->col_in  = fftwf_alloc_complex(nc);
    eng->col_out = fftwf_alloc_complex(nc);
    if (eng->col_in == NULL || eng->col_out == NULL)
    {
        return -1;
    }
    memset(eng->col_in, 0, sizeof(fftwf_complex) * nc);

    if (eng->plan_2d == NULL && eng->col_stride != eng->n_fft)
    {
        size_t nr    = (size_t) eng->pup_size * (size_t) eng->n_fft;
        eng->row_in  = fftwf_alloc_complex(nr);
        eng->row_out = fftwf_alloc_complex(nr);
        if (eng->row_in == NULL || eng->row_out == NULL)
        {
            return -1;
        }
        memset(eng->row_in, 0, sizeof(fftwf_complex) * nr);
    }
    return 0;
}

/**
 * atmturb_ao_psf_free_bufs - Release complex working arrays of PSF engine
 * @eng: PSF engine descriptor.
 */
void atmturb_ao_psf_free_bufs(
    atmturb_ao_psf_engine_t *eng)
{
    fftwf_free(eng->row_in);
    fftwf_free(eng->row_out);
    fftwf_free(eng->col_in);
    fftwf_free(eng->col_out);
    eng->row_in  = NULL;
    eng->row_out = NULL;
    eng->col_in  = NULL;
    eng->col_out = NULL;
}

/**
 * atmturb_ao_psf_engine_free - Release FFTW plans and complex arrays
 * @eng: Initialized PSF engine.
 */
void atmturb_ao_psf_engine_free(
    atmturb_ao_psf_engine_t *eng)
{
    if (eng->plan_rows)
    {
        fftwf_destroy_plan(eng->plan_rows);
    }
    if (eng->plan_cols)
    {
        fftwf_destroy_plan(eng->plan_cols);
    }
    if (eng->plan_2d)
    {
        fftwf_destroy_plan(eng->plan_2d);
    }
    eng->plan_rows = NULL;
    eng->plan_cols = NULL;
    eng->plan_2d   = NULL;
    atmturb_ao_psf_free_bufs(eng);
}

/**
 * atmturb_ao_psf_engine_init - Allocate buffers and create cached FFTW plans
 * @eng: PSF engine structure to initialize.
 * @pup_size: Pupil dimension in pixels.
 * @n_fft: Zero-padded FFT dimension (power of 2, >= 2 * pup_size).
 * @fftw_nth: Number of FFTW worker threads per transform.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_ao_psf_engine_init(
    atmturb_ao_psf_engine_t *eng,
    long                     pup_size,
    long                     n_fft,
    int                      fftw_nth)
{
    memset(eng, 0, sizeof(*eng));
    eng->pup_size    = pup_size;
    eng->n_fft       = n_fft;
    eng->n_psf       = (n_fft <= 512) ? n_fft : 512;
    int use_2d       = (fftw_nth <= 1 && eng->n_psf == n_fft);
    eng->col_stride  = use_2d ? n_fft : (n_fft + 8);
    eng->norm_factor = 1.0 / ((double) n_fft * (double) n_fft);

    if (atmturb_ao_psf_alloc_bufs(eng) != 0)
    {
        atmturb_ao_psf_free_bufs(eng);
        return -1;
    }

    atmturb_fftwf_plan_threads(fftw_nth);
    atmturb_fftwf_wisdom_load();
    unsigned int flags = atmturb_fftwf_plan_flags();
    int          nf    = (int) n_fft;
    int          cs    = (int) eng->col_stride;

    if (use_2d)
    {
        eng->plan_2d = fftwf_plan_dft_2d(nf, nf, eng->col_in, eng->col_out, FFTW_FORWARD, flags);
        if (!eng->plan_2d)
        {
            atmturb_ao_psf_engine_free(eng);
            return -1;
        }
    }
    else
    {
        eng->plan_rows = fftwf_plan_many_dft(1, &nf, (int) pup_size, eng->row_in, NULL, 1, nf,
                                             eng->row_out, NULL, 1, nf, FFTW_FORWARD, flags);
        eng->plan_cols = fftwf_plan_many_dft(1, &nf, (int) eng->n_psf, eng->col_in, NULL, 1, cs,
                                             eng->col_out, NULL, 1, cs, FFTW_FORWARD, flags);
        if (!eng->plan_rows || !eng->plan_cols)
        {
            atmturb_ao_psf_engine_free(eng);
            return -1;
        }
    }

    if (flags == FFTW_MEASURE)
    {
        atmturb_fftwf_wisdom_save();
    }
    if (eng->row_in != NULL)
    {
        memset(eng->row_in, 0, sizeof(fftwf_complex) * (size_t) pup_size * (size_t) n_fft);
    }
    memset(eng->col_in, 0, sizeof(fftwf_complex) * (size_t) eng->n_psf * (size_t) eng->col_stride);
    return 0;
}

/**
 * atmturb_ao_transpose_pruned - L1-tiled FFT-shifted transpose of active pupil rows into columns
 * @eng: Initialized PSF engine with populated @row_out.
 */
static void atmturb_ao_transpose_pruned(
    const atmturb_ao_psf_engine_t *eng)
{
    long half_fft   = eng->n_fft / 2;
    long u_start    = (eng->n_fft - eng->n_psf) / 2;
    long n_fft      = eng->n_fft;
    long n_psf      = eng->n_psf;
    long pup_size   = eng->pup_size;
    long col_stride = eng->col_stride;

    for (long j0 = 0; j0 < pup_size; j0 += 16)
    {
        long j_end = (j0 + 16 <= pup_size) ? (j0 + 16) : pup_size;
        for (long ui0 = 0; ui0 < n_psf; ui0 += 16)
        {
            long u0  = u_start + ui0;
            long fx0 = (u0 < half_fft) ? (u0 + half_fft) : (u0 - half_fft);
            for (long j = j0; j < j_end; j++)
            {
                const fftwf_complex *restrict src = &eng->row_out[j * n_fft + fx0];
                for (long du = 0; du < 16; du++)
                {
                    fftwf_complex *dst = &eng->col_in[(ui0 + du) * col_stride + j];
                    (*dst)[0]          = src[du][0];
                    (*dst)[1]          = src[du][1];
                }
            }
        }
    }
}

/**
 * atmturb_ao_extract_2d_psf - Quadrant-swapped SIMD PSF extraction from row-major 2D DFT output
 */
static float atmturb_ao_extract_2d_psf(
    const atmturb_ao_psf_engine_t *eng,
    float               *restrict  psf_buf,
    int                            accum)
{
    long  n_fft = eng->n_fft;
    long  half  = n_fft / 2;
    float norm  = (float) eng->norm_factor;
    float peak  = 0.0f;

    for (long vj = 0; vj < n_fft; vj++)
    {
        long                          fy  = (vj < half) ? (vj + half) : (vj - half);
        const fftwf_complex *restrict row = &eng->col_out[fy * n_fft];
        float               *restrict dst = &psf_buf[vj * n_fft];
        for (int q = 0; q < 2; q++)
        {
            const fftwf_complex *restrict s = (q == 0) ? &row[half] : &row[0];
            float               *restrict d = (q == 0) ? &dst[0]    : &dst[half];
            #pragma omp simd reduction(max:peak)
            for (long u = 0; u < half; u++)
            {
                float val = (s[u][0] * s[u][0] + s[u][1] * s[u][1]) * norm;
                d[u]      = accum ? (d[u] + val) : val;
                if (val > peak)
                {
                    peak = val;
                }
            }
        }
    }
    return peak;
}

/**
 * atmturb_ao_extract_pruned_psf - L1-tiled intensity extraction from transposed column FFTs
 */
static float atmturb_ao_extract_pruned_psf(
    const atmturb_ao_psf_engine_t *eng,
    float               *restrict  psf_buf,
    int                            accum)
{
    long  half_fft   = eng->n_fft / 2;
    long  v_start    = (eng->n_fft - eng->n_psf) / 2;
    long  n_psf      = eng->n_psf;
    long  col_stride = eng->col_stride;
    float norm       = (float) eng->norm_factor;
    float peak       = 0.0f;
    float tile[16][16];

    for (long vj0 = 0; vj0 < n_psf; vj0 += 16)
    {
        long v0  = v_start + vj0;
        long fy0 = (v0 < half_fft) ? (v0 + half_fft) : (v0 - half_fft);
        for (long ui0 = 0; ui0 < n_psf; ui0 += 16)
        {
            for (int du = 0; du < 16; du++)
            {
                const fftwf_complex *restrict src = &eng->col_out[(ui0 + du) * col_stride + fy0];
                #pragma omp simd reduction(max:peak)
                for (int dv = 0; dv < 16; dv++)
                {
                    float val    = (src[dv][0] * src[dv][0] + src[dv][1] * src[dv][1]) * norm;
                    tile[du][dv] = val;
                    if (val > peak)
                    {
                        peak = val;
                    }
                }
            }
            for (int dv = 0; dv < 16; dv++)
            {
                float *restrict dst = &psf_buf[(vj0 + dv) * n_psf + ui0];
                #pragma omp simd
                for (int du = 0; du < 16; du++)
                {
                    dst[du] = accum ? (dst[du] + tile[du][dv]) : tile[du][dv];
                }
            }
        }
    }
    return peak;
}

/**
 * atmturb_ao_render_psf_amp - Form centered PSF with optional scintillation amplitude map
 */
double atmturb_ao_render_psf_amp(
    const atmturb_ao_psf_engine_t *eng,
    const float                   *phase_in,
    const float                   *amp_in,
    const float                   *pupil_mask,
    float                          phase_scale,
    float                         *psf_buf,
    int                            accum)
{
    long           pup_size = eng->pup_size;
    long           n_fft    = eng->n_fft;
    fftwf_complex *in_b     = (eng->plan_2d != NULL) ? eng->col_in : eng->row_in;

    for (long j = 0; j < pup_size; j++)
    {
        long         pup_row = j * pup_size;
        const float *p_row   = (phase_in != NULL) ? &phase_in[pup_row] : NULL;
        atmturb_fresnel_simd_phasor_masked(&in_b[j * n_fft], p_row,
                                           &pupil_mask[pup_row], phase_scale, pup_size);
        if (amp_in != NULL)
        {
            const float   *restrict a_row = &amp_in[pup_row];
            fftwf_complex *restrict dst   = &in_b[j * n_fft];
            #pragma omp simd
            for (long i = 0; i < pup_size; i++)
            {
                float a    = a_row[i];
                dst[i][0] *= a;
                dst[i][1] *= a;
            }
        }
    }

    if (eng->plan_2d != NULL)
    {
        fftwf_execute_dft(eng->plan_2d, eng->col_in, eng->col_out);
        return (double) atmturb_ao_extract_2d_psf(eng, psf_buf, accum);
    }

    fftwf_execute_dft(eng->plan_rows, eng->row_in, eng->row_out);
    atmturb_ao_transpose_pruned(eng);
    fftwf_execute_dft(eng->plan_cols, eng->col_in, eng->col_out);
    return (double) atmturb_ao_extract_pruned_psf(eng, psf_buf, accum);
}

/**
 * atmturb_ao_render_psf - Form instantaneous zero-padded centered PSF
 */
double atmturb_ao_render_psf(
    const atmturb_ao_psf_engine_t *eng,
    const float                   *phase_in,
    const float                   *pupil_mask,
    long                           pup_size,
    float                          phase_scale,
    float                         *psf_buf,
    int                            accum)
{
    (void) pup_size;
    return atmturb_ao_render_psf_amp(eng, phase_in, NULL, pupil_mask, phase_scale, psf_buf, accum);
}
