// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_simd_ops_impl.h
 * @brief   Compiler-vectorized 1D array, piston removal, and tiled subharmonic kernels
 */

#ifndef ATMTURB_SIMD_OPS_IMPL_H
#define ATMTURB_SIMD_OPS_IMPL_H

#include <math.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "atmturb_lowfreq.h"
#include "atmturb_simd.h"

#define ATMTURB_LF_TILE 256

/**
 * struct atmturb_lf_phasors_t - Precomputed 2D phasor recurrence tables for subharmonics
 * @nmodes: Active mode count (clamped to ATMTURB_LOWFREQ_NMODES).
 * @nmodes_pad4: Active mode count rounded up to a multiple of 4.
 * @c0_re: Initial mode real amplitudes at (x0, y0) scaled by layer weight.
 * @c0_im: Initial mode imaginary amplitudes at (x0, y0) scaled by layer weight.
 * @dx_re: Per-pixel X step cosine for each mode.
 * @dx_im: Per-pixel X step sine for each mode.
 * @dy_re: Per-pixel Y step cosine for each mode.
 * @dy_im: Per-pixel Y step sine for each mode.
 * @kx_os: Scaled horizontal wavenumber os * kx for tile offset.
 * @ky_os: Scaled vertical wavenumber os * ky for OpenMP row-chunk offset.
 */
typedef struct
{
    int    nmodes;
    int    nmodes_pad4;
    double c0_re[ATMTURB_LOWFREQ_NMODES];
    double c0_im[ATMTURB_LOWFREQ_NMODES];
    double dx_re[ATMTURB_LOWFREQ_NMODES];
    double dx_im[ATMTURB_LOWFREQ_NMODES];
    double dy_re[ATMTURB_LOWFREQ_NMODES];
    double dy_im[ATMTURB_LOWFREQ_NMODES];
    double kx_os[ATMTURB_LOWFREQ_NMODES];
    double ky_os[ATMTURB_LOWFREQ_NMODES];
} atmturb_lf_phasors_t;

/**
 * atmturb_scale_float_array_impl - Compiler-vectorized float array scaling
 */
static inline void atmturb_scale_float_array_impl(
    float       *dest,
    const float *src,
    float        scale,
    long         n)
{
    if (dest == NULL || src == NULL || n <= 0)
    {
        return;
    }
    if (dest == src)
    {
        float *restrict inout = dest;
        #pragma omp simd
        for (long i = 0; i < n; i++)
        {
            inout[i] *= scale;
        }
        return;
    }

    float       *restrict d = dest;
    const float *restrict s = src;
    #pragma omp simd
    for (long i = 0; i < n; i++)
    {
        d[i] = s[i] * scale;
    }
}

/**
 * atmturb_init_phase_amp_impl - Zero phase via libc memset and fill unit amplitude via SIMD
 */
static inline void atmturb_init_phase_amp_impl(
    float *pha,
    float *amp,
    long   n)
{
    if (n <= 0)
    {
        return;
    }
    if (pha != NULL)
    {
        memset(pha, 0, sizeof(float) * (size_t) n);
    }
    if (amp != NULL)
    {
        float *restrict a = amp;
        #pragma omp simd
        for (long i = 0; i < n; i++)
        {
            a[i] = 1.0f;
        }
    }
}

/**
 * atmturb_add_float_array_impl - Compiler-vectorized pointwise float accumulation
 */
static inline void atmturb_add_float_array_impl(
    float       *dest,
    const float *src,
    long         n)
{
    if (dest == NULL || src == NULL || n <= 0)
    {
        return;
    }
    if (dest == src)
    {
        float *restrict inout = dest;
        #pragma omp simd
        for (long i = 0; i < n; i++)
        {
            inout[i] += inout[i];
        }
        return;
    }

    float       *restrict d = dest;
    const float *restrict s = src;
    #pragma omp simd
    for (long i = 0; i < n; i++)
    {
        d[i] += s[i];
    }
}

/**
 * atmturb_complex_mul_inplace_impl - In-place interleaved complex float multiplication
 */
static inline void atmturb_complex_mul_inplace_impl(
    float       *restrict inout,
    const float *restrict other,
    long                  n_complex)
{
    #pragma omp simd
    for (long i = 0; i < n_complex; i++)
    {
        float r1 = inout[2 * i];
        float i1 = inout[2 * i + 1];
        float r2 = other[2 * i];
        float i2 = other[2 * i + 1];

        inout[2 * i]     = r1 * r2 - i1 * i2;
        inout[2 * i + 1] = r1 * i2 + i1 * r2;
    }
}

/**
 * atmturb_complex_mul_array_impl - Compiler-vectorized interleaved complex array product
 */
static inline void atmturb_complex_mul_array_impl(
    float       *dest,
    const float *src1,
    const float *src2,
    long         n_complex)
{
    if (dest == NULL || src1 == NULL || src2 == NULL || n_complex <= 0)
    {
        return;
    }
    if (dest == src1 && src1 != src2)
    {
        atmturb_complex_mul_inplace_impl(dest, src2, n_complex);
        return;
    }
    if (dest == src2 && src1 != src2)
    {
        atmturb_complex_mul_inplace_impl(dest, src1, n_complex);
        return;
    }
    if (dest == src1 && src1 == src2)
    {
        float *restrict inout = dest;
        #pragma omp simd
        for (long i = 0; i < n_complex; i++)
        {
            float r1 = inout[2 * i];
            float i1 = inout[2 * i + 1];
            inout[2 * i]     = r1 * r1 - i1 * i1;
            inout[2 * i + 1] = 2.0f * r1 * i1;
        }
        return;
    }

    float       *restrict d  = dest;
    const float *restrict s1 = src1;
    const float *restrict s2 = src2;

    #pragma omp simd
    for (long i = 0; i < n_complex; i++)
    {
        float r1 = s1[2 * i];
        float i1 = s1[2 * i + 1];
        float r2 = s2[2 * i];
        float i2 = s2[2 * i + 1];

        d[2 * i]     = r1 * r2 - i1 * i2;
        d[2 * i + 1] = r1 * i2 + i1 * r2;
    }
}

/**
 * atmturb_sum_float_chunked - 4-way ILP SIMD reduction with chunked FP64 accumulation
 */
static inline double atmturb_sum_float_chunked(
    const float *restrict src,
    long                  npix)
{
    double total = 0.0;
    long   base  = 0;

    for (; base + 1024 <= npix; base += 1024)
    {
        float                 s0 = 0.0f, s1 = 0.0f, s2 = 0.0f, s3 = 0.0f;
        const float *restrict p  = src + base;
        #pragma omp simd reduction(+ : s0, s1, s2, s3)
        for (long k = 0; k < 256; k++)
        {
            s0 += p[k];
            s1 += p[k + 256];
            s2 += p[k + 512];
            s3 += p[k + 768];
        }
        total += (double) s0 + (double) s1 + (double) s2 + (double) s3;
    }

    float tail = 0.0f;
    #pragma omp simd reduction(+ : tail)
    for (long i = base; i < npix; i++)
    {
        tail += src[i];
    }
    return total + (double) tail;
}

/**
 * atmturb_remove_piston_stream_impl - Vectorized FP64 piston reduction and subtraction
 */
static inline void atmturb_remove_piston_stream_impl(
    float       *dst,
    const float *src,
    long         npix)
{
    if (dst == NULL || src == NULL || npix <= 0)
    {
        return;
    }

    const double sum  = atmturb_sum_float_chunked(src, npix);
    const float  mean = (float) (sum / (double) npix);
    if (dst == src)
    {
        float *restrict inout = dst;
        #pragma omp simd
        for (long i = 0; i < npix; i++)
        {
            inout[i] -= mean;
        }
        return;
    }

    float       *restrict d = dst;
    const float *restrict s = src;
    #pragma omp simd
    for (long i = 0; i < npix; i++)
    {
        d[i] = s[i] - mean;
    }
}

/**
 * atmturb_lf_init_phasors - Precompute 2D phasor steps for active subharmonic modes
 */
static inline void atmturb_lf_init_phasors(
    const atmturb_lowfreq_params_t *params,
    const float                    *are,
    const float                    *aim,
    atmturb_lf_phasors_t           *ph)
{
    const atmturb_lowfreq_t *lf = params->lf;
    const double             os = (double) ((params->os > 0) ? params->os : 1);
    const double             w  = (double) params->weight;
    int                      nm = lf->nmodes;
    if (nm > ATMTURB_LOWFREQ_NMODES)
    {
        nm = ATMTURB_LOWFREQ_NMODES;
    }

    memset(ph, 0, sizeof(*ph));
    ph->nmodes      = nm;
    ph->nmodes_pad4 = (nm + 3) & ~3;

    for (int m = 0; m < nm; m++)
    {
        double kx = (double) lf->kx[m];
        double ky = (double) lf->ky[m];
        double p0 = params->x0 * kx + params->y0 * ky;
        double s0, c0;
        sincos(p0, &s0, &c0);

        double ar    = (double) are[m];
        double ai    = (double) aim[m];
        ph->c0_re[m] = w * (ar * c0 - ai * s0);
        ph->c0_im[m] = w * (ar * s0 + ai * c0);

        ph->kx_os[m] = os * kx;
        ph->ky_os[m] = os * ky;
        sincos(ph->kx_os[m], &ph->dx_im[m], &ph->dx_re[m]);
        sincos(ph->ky_os[m], &ph->dy_im[m], &ph->dy_re[m]);
    }
}

/**
 * atmturb_lf_build_tile_basis - Build horizontal basis table for a column tile via phasor step
 */
static inline void atmturb_lf_build_tile_basis(
    const atmturb_lf_phasors_t *ph,
    long                        col0,
    long                        tlen,
    float                       h_re[ATMTURB_LOWFREQ_NMODES][ATMTURB_LF_TILE],
    float                       h_im[ATMTURB_LOWFREQ_NMODES][ATMTURB_LF_TILE])
{
    for (int m = 0; m < ph->nmodes; m++)
    {
        double ux, vx;
        sincos((double) col0 * ph->kx_os[m], &vx, &ux);
        const double s_re = ph->dx_re[m];
        const double s_im = ph->dx_im[m];

        for (long ii = 0; ii < tlen; ii++)
        {
            h_re[m][ii] = (float) ux;
            h_im[m][ii] = (float) vx;
            double nux  = ux * s_re - vx * s_im;
            double nvx  = vx * s_re + ux * s_im;
            ux          = nux;
            vx          = nvx;
        }
    }
    for (int m = ph->nmodes; m < ph->nmodes_pad4; m++)
    {
        memset(h_re[m], 0, sizeof(float) * (size_t) tlen);
        memset(h_im[m], 0, sizeof(float) * (size_t) tlen);
    }
}

/**
 * atmturb_lf_eval_row_tile - Accumulate active modes for one row tile in L1 cache
 */
static inline void atmturb_lf_eval_row_tile(
    const float           cr[ATMTURB_LOWFREQ_NMODES],
    const float           ci[ATMTURB_LOWFREQ_NMODES],
    int                   nmodes_pad4,
    const float           h_re[ATMTURB_LOWFREQ_NMODES][ATMTURB_LF_TILE],
    const float           h_im[ATMTURB_LOWFREQ_NMODES][ATMTURB_LF_TILE],
    float       *restrict out_tile,
    long                  tlen)
{
    float acc[ATMTURB_LF_TILE];
    memcpy(acc, out_tile, sizeof(float) * (size_t) tlen);

    for (int m = 0; m < nmodes_pad4; m += 4)
    {
        const float cr0 = cr[m + 0], ci0 = ci[m + 0];
        const float cr1 = cr[m + 1], ci1 = ci[m + 1];
        const float cr2 = cr[m + 2], ci2 = ci[m + 2];
        const float cr3 = cr[m + 3], ci3 = ci[m + 3];

        const float *restrict r0 = h_re[m + 0], *restrict i0 = h_im[m + 0];
        const float *restrict r1 = h_re[m + 1], *restrict i1 = h_im[m + 1];
        const float *restrict r2 = h_re[m + 2], *restrict i2 = h_im[m + 2];
        const float *restrict r3 = h_re[m + 3], *restrict i3 = h_im[m + 3];

        #pragma omp simd
        for (long ii = 0; ii < tlen; ii++)
        {
            float v0 = cr0 * r0[ii] - ci0 * i0[ii];
            float v1 = cr1 * r1[ii] - ci1 * i1[ii];
            float v2 = cr2 * r2[ii] - ci2 * i2[ii];
            float v3 = cr3 * r3[ii] - ci3 * i3[ii];
            acc[ii] += (v0 + v1) + (v2 + v3);
        }
    }

    memcpy(out_tile, acc, sizeof(float) * (size_t) tlen);
}

/**
 * atmturb_lf_eval_row_range - Evaluate all column tiles across row range [r_beg, r_end)
 */
static inline void atmturb_lf_eval_row_range(
    const atmturb_lf_phasors_t *ph,
    float             *restrict out_pha,
    long                        n,
    long                        r_beg,
    long                        r_end)
{
    float h_re[ATMTURB_LOWFREQ_NMODES][ATMTURB_LF_TILE];
    float h_im[ATMTURB_LOWFREQ_NMODES][ATMTURB_LF_TILE];

    for (long col0 = 0; col0 < n; col0 += ATMTURB_LF_TILE)
    {
        long tlen = (col0 + ATMTURB_LF_TILE <= n) ? ATMTURB_LF_TILE : (n - col0);
        atmturb_lf_build_tile_basis(ph, col0, tlen, h_re, h_im);

        double cr_d[ATMTURB_LOWFREQ_NMODES] = {0.0};
        double ci_d[ATMTURB_LOWFREQ_NMODES] = {0.0};
        for (int m = 0; m < ph->nmodes; m++)
        {
            if (r_beg == 0)
            {
                cr_d[m] = ph->c0_re[m];
                ci_d[m] = ph->c0_im[m];
            }
            else
            {
                double sy, cy;
                sincos((double) r_beg * ph->ky_os[m], &sy, &cy);
                cr_d[m] = ph->c0_re[m] * cy - ph->c0_im[m] * sy;
                ci_d[m] = ph->c0_re[m] * sy + ph->c0_im[m] * cy;
            }
        }

        for (long row = r_beg; row < r_end; row++)
        {
            float cr[ATMTURB_LOWFREQ_NMODES] = {0.0f};
            float ci[ATMTURB_LOWFREQ_NMODES] = {0.0f};
            for (int m = 0; m < ph->nmodes; m++)
            {
                cr[m]      = (float) cr_d[m];
                ci[m]      = (float) ci_d[m];
                double ncr = cr_d[m] * ph->dy_re[m] - ci_d[m] * ph->dy_im[m];
                double nci = ci_d[m] * ph->dy_re[m] + cr_d[m] * ph->dy_im[m];
                cr_d[m]    = ncr;
                ci_d[m]    = nci;
            }
            atmturb_lf_eval_row_tile(cr, ci, ph->nmodes_pad4, h_re, h_im,
                                     &out_pha[row * n + col0], tlen);
        }
    }
}

/**
 * atmturb_extrude_lowfreq_impl - L1-tiled multi-mode fused subharmonic evaluator
 */
static inline void atmturb_extrude_lowfreq_impl(
    const atmturb_lowfreq_params_t *params)
{
    const atmturb_lowfreq_t *lf = params->lf;
    const long               n  = params->pup_size;
    if (lf == NULL || lf->nmodes <= 0 || n <= 0 || params->out_pha == NULL)
    {
        return;
    }

    const float *are = (params->custom_are != NULL && params->custom_aim != NULL)
                           ? params->custom_are
                           : ((params->screen_idx == 1) ? lf->bre : lf->are);
    const float *aim = (params->custom_are != NULL && params->custom_aim != NULL)
                           ? params->custom_aim
                           : ((params->screen_idx == 1) ? lf->bim : lf->aim);

    atmturb_lf_phasors_t ph;
    atmturb_lf_init_phasors(params, are, aim, &ph);

#ifdef _OPENMP
    if (n >= 512 && !omp_in_parallel())
    {
        int max_t = omp_get_max_threads();
        int nth   = (max_t > 8) ? 8 : max_t;
        #pragma omp parallel num_threads(nth)
        {
            int  tid   = omp_get_thread_num();
            int  nt    = omp_get_num_threads();
            long r_beg = (n * (long) tid) / (long) nt;
            long r_end = (n * (long) (tid + 1)) / (long) nt;
            atmturb_lf_eval_row_range(&ph, params->out_pha, n, r_beg, r_end);
        }
        return;
    }
#endif
    atmturb_lf_eval_row_range(&ph, params->out_pha, n, 0, n);
}

#endif // ATMTURB_SIMD_OPS_IMPL_H
