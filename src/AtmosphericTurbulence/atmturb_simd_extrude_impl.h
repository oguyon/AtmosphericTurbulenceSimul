// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_simd_extrude_impl.h
 * @brief   Compiler-vectorized 2D phase-screen extrusion kernels (bilinear & Keys bicubic)
 */

#ifndef ATMTURB_SIMD_EXTRUDE_IMPL_H
#define ATMTURB_SIMD_EXTRUDE_IMPL_H

#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "atmturb_simd.h"

/**
 * struct atmturb_bicubic_ctx_t - Row pointers and separable weights for Keys bicubic row
 * @r0: Row y - 1 pointer.
 * @r1: Row y pointer.
 * @r2: Row y + 1 pointer.
 * @r3: Row y + 2 pointer.
 * @wx: Horizontal Keys weights scaled by layer weight.
 * @wy: Vertical Keys weights.
 */
typedef struct
{
    const float *r0;
    const float *r1;
    const float *r2;
    const float *r3;
    float        wx[4];
    float        wy[4];
} atmturb_bicubic_ctx_t;

/**
 * atmturb_extrude_omp_threads - Bound OpenMP row-parallel team size on hybrid/many-core CPUs
 *
 * Return: Optimal thread count (1 to 8) for single-frame large-pupil extrusion.
 */
static inline int atmturb_extrude_omp_threads(void)
{
#ifdef _OPENMP
    int max_t = omp_get_max_threads();
    return (max_t > 8) ? 8 : max_t;
#else
    return 1;
#endif
}

/**
 * atmturb_keys_weights_impl - Compute 4-tap Keys cubic convolution weights (a = -0.5)
 * @f: Fractional coordinate in [0, 1).
 * @w: Output array of 4 weights.
 */
static inline void atmturb_keys_weights_impl(
    float  f,
    float *w)
{
    float f2 = f * f;
    float f3 = f2 * f;

    w[0] = -0.5f * f3 + f2 - 0.5f * f;
    w[1] =  1.5f * f3 - 2.5f * f2 + 1.0f;
    w[2] = -1.5f * f3 + 2.0f * f2 + 0.5f * f;
    w[3] =  0.5f * f3 - 0.5f * f2;
}

/**
 * atmturb_bilinear_span_os1 - Contiguous unit-stride bilinear row span
 */
static inline void atmturb_bilinear_span_os1(
    const float *restrict r0,
    const float *restrict r1,
    float       *restrict dst,
    long                  len,
    const float           w[4])
{
    const float w00 = w[0];
    const float w10 = w[1];
    const float w01 = w[2];
    const float w11 = w[3];

    #pragma omp simd
    for (long ii = 0; ii < len; ii++)
    {
        dst[ii] += w00 * r0[ii] + w10 * r0[ii + 1] + w01 * r1[ii] + w11 * r1[ii + 1];
    }
}

/**
 * atmturb_bilinear_span_os2 - Contiguous stride-2 bilinear row span (deinterleaved by compiler)
 */
static inline void atmturb_bilinear_span_os2(
    const float *restrict r0,
    const float *restrict r1,
    float       *restrict dst,
    long                  len,
    const float           w[4])
{
    const float w00 = w[0];
    const float w10 = w[1];
    const float w01 = w[2];
    const float w11 = w[3];

    #pragma omp simd
    for (long ii = 0; ii < len; ii++)
    {
        long sx = 2 * ii;
        dst[ii] += w00 * r0[sx] + w10 * r0[sx + 1] + w01 * r1[sx] + w11 * r1[sx + 1];
    }
}

/**
 * atmturb_bilinear_row_impl - Accumulate one pupil row with seam-splitting vectorization
 */
static inline void atmturb_bilinear_row_impl(
    const float *restrict r0,
    const float *restrict r1,
    float       *restrict dst,
    const long            geom_i[4],
    const float           w[4])
{
    const long msize    = geom_i[0];
    const long start_x  = geom_i[1];
    const long pup_size = geom_i[2];
    const long os       = geom_i[3];

    if (os == 1 && pup_size <= msize)
    {
        if (start_x + pup_size < msize)
        {
            atmturb_bilinear_span_os1(&r0[start_x], &r1[start_x], dst, pup_size, w);
            return;
        }
        long seg1 = msize - 1 - start_x;
        if (seg1 > 0)
        {
            atmturb_bilinear_span_os1(&r0[start_x], &r1[start_x], dst, seg1, w);
        }
        dst[seg1] += w[0] * r0[msize - 1] + w[1] * r0[0] + w[2] * r1[msize - 1] + w[3] * r1[0];
        if (pup_size - seg1 - 1 > 0)
        {
            atmturb_bilinear_span_os1(r0, r1, &dst[seg1 + 1], pup_size - seg1 - 1, w);
        }
        return;
    }
    if (os == 2 && 2 * pup_size <= msize)
    {
        if (start_x + 2 * (pup_size - 1) + 1 < msize)
        {
            atmturb_bilinear_span_os2(&r0[start_x], &r1[start_x], dst, pup_size, w);
            return;
        }
        long seg1 = (start_x <= msize - 2) ? ((msize - 2 - start_x) / 2 + 1) : 0;
        if (seg1 > 0)
        {
            atmturb_bilinear_span_os2(&r0[start_x], &r1[start_x], dst, seg1, w);
        }
        long ix0 = (start_x + 2 * seg1) % msize;
        long ix1 = (ix0 + 1 == msize) ? 0 : (ix0 + 1);
        dst[seg1] += w[0] * r0[ix0] + w[1] * r0[ix1] + w[2] * r1[ix0] + w[3] * r1[ix1];
        if (pup_size - seg1 - 1 > 0)
        {
            long sx2 = (start_x + 2 * (seg1 + 1)) % msize;
            atmturb_bilinear_span_os2(&r0[sx2], &r1[sx2], &dst[seg1 + 1], pup_size - seg1 - 1, w);
        }
        return;
    }

    for (long ii = 0; ii < pup_size; ii++)
    {
        long ix0 = (start_x + ii * os) % msize;
        long ix1 = (ix0 + 1 == msize) ? 0 : (ix0 + 1);
        dst[ii] += w[0] * r0[ix0] + w[1] * r0[ix1] + w[2] * r1[ix0] + w[3] * r1[ix1];
    }
}

/**
 * atmturb_extrude_bilinear_rows - Render pupil row slice [j0, j1) for bilinear extrusion
 */
static inline void atmturb_extrude_bilinear_rows(
    const atmturb_extrude_params_t *params,
    const long                      geom_i[4],
    const float                     w[4],
    long                            start_y,
    long                            j0,
    long                            j1)
{
    const float *restrict master   = params->master;
    float       *restrict out_pha  = params->out_pha;
    const long            msize    = geom_i[0];
    const long            pup_size = geom_i[2];
    const long            os       = geom_i[3];

    for (long jj = j0; jj < j1; jj++)
    {
        long         iy0 = (start_y + jj * os) % msize;
        long         iy1 = (iy0 + 1 == msize) ? 0 : (iy0 + 1);
        const float *r0  = &master[iy0 * msize];
        const float *r1  = &master[iy1 * msize];
        float       *dst = &out_pha[jj * pup_size];

        atmturb_bilinear_row_impl(r0, r1, dst, geom_i, w);
    }
}

/**
 * atmturb_extrude_bilinear_impl - Compiler-vectorized 2D bilinear phase extrusion
 */
static inline void atmturb_extrude_bilinear_impl(
    const atmturb_extrude_params_t *params)
{
    const long  msize    = params->msize;
    const long  pup_size = params->pup_size;
    const long  os       = (params->os > 0) ? params->os : 1;
    const float weight   = params->weight;

    double rx = fmod(params->x0, (double) msize);
    if (rx < 0.0)
    {
        rx += (double) msize;
    }
    double ry = fmod(params->y0, (double) msize);
    if (ry < 0.0)
    {
        ry += (double) msize;
    }

    long  start_x = (long) floor(rx);
    long  start_y = (long) floor(ry);
    float fx      = (float) (rx - (double) start_x);
    float fy      = (float) (ry - (double) start_y);
    start_x %= msize;
    start_y %= msize;

    const float w[4] = {
        weight * (1.0f - fx) * (1.0f - fy),
        weight * fx * (1.0f - fy),
        weight * (1.0f - fx) * fy,
        weight * fx * fy
    };
    const long geom_i[4] = {msize, start_x, pup_size, os};

#ifdef _OPENMP
    if (pup_size >= 512 && !omp_in_parallel())
    {
        int nth = atmturb_extrude_omp_threads();
        #pragma omp parallel num_threads(nth)
        {
            int  tid = omp_get_thread_num();
            int  nt  = omp_get_num_threads();
            long j0  = (pup_size * (long) tid) / (long) nt;
            long j1  = (pup_size * (long) (tid + 1)) / (long) nt;
            atmturb_extrude_bilinear_rows(params, geom_i, w, start_y, j0, j1);
        }
        return;
    }
#endif
    atmturb_extrude_bilinear_rows(params, geom_i, w, start_y, 0, pup_size);
}

/**
 * atmturb_bicubic_span_os1 - Contiguous unit-stride Keys bicubic row span
 */
static inline void atmturb_bicubic_span_os1(
    const atmturb_bicubic_ctx_t *bc,
    long                         off,
    float             *restrict  dst,
    long                         len)
{
    const float *restrict r0  = bc->r0 + off;
    const float *restrict r1  = bc->r1 + off;
    const float *restrict r2  = bc->r2 + off;
    const float *restrict r3  = bc->r3 + off;
    const float           wx0 = bc->wx[0], wx1 = bc->wx[1], wx2 = bc->wx[2], wx3 = bc->wx[3];
    const float           wy0 = bc->wy[0], wy1 = bc->wy[1], wy2 = bc->wy[2], wy3 = bc->wy[3];

    #pragma omp simd
    for (long ii = 0; ii < len; ii++)
    {
        float c0 = wx0 * r0[ii - 1] + wx1 * r0[ii] + wx2 * r0[ii + 1] + wx3 * r0[ii + 2];
        float c1 = wx0 * r1[ii - 1] + wx1 * r1[ii] + wx2 * r1[ii + 1] + wx3 * r1[ii + 2];
        float c2 = wx0 * r2[ii - 1] + wx1 * r2[ii] + wx2 * r2[ii + 1] + wx3 * r2[ii + 2];
        float c3 = wx0 * r3[ii - 1] + wx1 * r3[ii] + wx2 * r3[ii + 1] + wx3 * r3[ii + 2];
        dst[ii] += wy0 * c0 + wy1 * c1 + wy2 * c2 + wy3 * c3;
    }
}

/**
 * atmturb_bicubic_span_os2 - Contiguous stride-2 Keys bicubic row span
 */
static inline void atmturb_bicubic_span_os2(
    const atmturb_bicubic_ctx_t *bc,
    long                         off,
    float             *restrict  dst,
    long                         len)
{
    const float *restrict r0  = bc->r0 + off;
    const float *restrict r1  = bc->r1 + off;
    const float *restrict r2  = bc->r2 + off;
    const float *restrict r3  = bc->r3 + off;
    const float           wx0 = bc->wx[0], wx1 = bc->wx[1], wx2 = bc->wx[2], wx3 = bc->wx[3];
    const float           wy0 = bc->wy[0], wy1 = bc->wy[1], wy2 = bc->wy[2], wy3 = bc->wy[3];

    #pragma omp simd
    for (long ii = 0; ii < len; ii++)
    {
        long  sx = 2 * ii;
        float c0 = wx0 * r0[sx - 1] + wx1 * r0[sx] + wx2 * r0[sx + 1] + wx3 * r0[sx + 2];
        float c1 = wx0 * r1[sx - 1] + wx1 * r1[sx] + wx2 * r1[sx + 1] + wx3 * r1[sx + 2];
        float c2 = wx0 * r2[sx - 1] + wx1 * r2[sx] + wx2 * r2[sx + 1] + wx3 * r2[sx + 2];
        float c3 = wx0 * r3[sx - 1] + wx1 * r3[sx] + wx2 * r3[sx + 1] + wx3 * r3[sx + 2];
        dst[ii] += wy0 * c0 + wy1 * c1 + wy2 * c2 + wy3 * c3;
    }
}

/**
 * atmturb_bicubic_pixel_wrap - Evaluate one wrapped Keys bicubic pixel
 */
static inline float atmturb_bicubic_pixel_wrap(
    const atmturb_bicubic_ctx_t *bc,
    long                         sx,
    long                         msize)
{
    long xm1 = (sx == 0) ? (msize - 1) : (sx - 1);
    long xp1 = (sx + 1 == msize) ? 0 : (sx + 1);
    long xp2 = (xp1 + 1 == msize) ? 0 : (xp1 + 1);

    float c0 = bc->wx[0] * bc->r0[xm1] + bc->wx[1] * bc->r0[sx] +
               bc->wx[2] * bc->r0[xp1] + bc->wx[3] * bc->r0[xp2];
    float c1 = bc->wx[0] * bc->r1[xm1] + bc->wx[1] * bc->r1[sx] +
               bc->wx[2] * bc->r1[xp1] + bc->wx[3] * bc->r1[xp2];
    float c2 = bc->wx[0] * bc->r2[xm1] + bc->wx[1] * bc->r2[sx] +
               bc->wx[2] * bc->r2[xp1] + bc->wx[3] * bc->r2[xp2];
    float c3 = bc->wx[0] * bc->r3[xm1] + bc->wx[1] * bc->r3[sx] +
               bc->wx[2] * bc->r3[xp1] + bc->wx[3] * bc->r3[xp2];

    return bc->wy[0] * c0 + bc->wy[1] * c1 + bc->wy[2] * c2 + bc->wy[3] * c3;
}

/**
 * atmturb_bicubic_row_impl - Accumulate one Keys bicubic row with seam-peeled vectorization
 */
static inline void atmturb_bicubic_row_impl(
    const atmturb_bicubic_ctx_t *bc,
    float             *restrict  dst,
    const long                   geom_i[4])
{
    const long msize    = geom_i[0];
    const long start_x  = geom_i[1];
    const long pup_size = geom_i[2];
    const long os       = geom_i[3];

    if (start_x >= 1 && start_x + (pup_size - 1) * os + 2 < msize)
    {
        if (os == 1)
        {
            atmturb_bicubic_span_os1(bc, start_x, dst, pup_size);
            return;
        }
        if (os == 2)
        {
            atmturb_bicubic_span_os2(bc, start_x, dst, pup_size);
            return;
        }
    }

    long ii = 0;
    while (ii < pup_size)
    {
        long sx = (start_x + ii * os) % msize;
        if ((os == 1 || os == 2) && sx >= 1 && sx + 2 < msize)
        {
            long rem = pup_size - ii;
            long fit = (msize - 3 - sx) / os + 1;
            long len = (rem < fit) ? rem : fit;
            if (os == 1)
            {
                atmturb_bicubic_span_os1(bc, sx, &dst[ii], len);
            }
            else
            {
                atmturb_bicubic_span_os2(bc, sx, &dst[ii], len);
            }
            ii += len;
        }
        else
        {
            dst[ii] += atmturb_bicubic_pixel_wrap(bc, sx, msize);
            ii++;
        }
    }
}

/**
 * atmturb_extrude_bicubic_rows - Render pupil row slice [j0, j1) for Keys bicubic extrusion
 */
static inline void atmturb_extrude_bicubic_rows(
    const atmturb_extrude_params_t *params,
    const long                      geom_i[4],
    const float                     wx_w[4],
    const float                     wy_raw[4],
    long                            start_y,
    long                            j0,
    long                            j1)
{
    const float *restrict master   = params->master;
    float       *restrict out_pha  = params->out_pha;
    const long            msize    = geom_i[0];
    const long            pup_size = geom_i[2];
    const long            os       = geom_i[3];

    for (long jj = j0; jj < j1; jj++)
    {
        long iy0  = (start_y + jj * os) % msize;
        long iym1 = (iy0 == 0) ? (msize - 1) : (iy0 - 1);
        long iyp1 = (iy0 + 1 == msize) ? 0 : (iy0 + 1);
        long iyp2 = (iyp1 + 1 == msize) ? 0 : (iyp1 + 1);

        atmturb_bicubic_ctx_t bc = {
            .r0 = &master[iym1 * msize],
            .r1 = &master[iy0 * msize],
            .r2 = &master[iyp1 * msize],
            .r3 = &master[iyp2 * msize],
            .wx = {wx_w[0], wx_w[1], wx_w[2], wx_w[3]},
            .wy = {wy_raw[0], wy_raw[1], wy_raw[2], wy_raw[3]}
        };
        atmturb_bicubic_row_impl(&bc, &out_pha[jj * pup_size], geom_i);
    }
}

/**
 * atmturb_extrude_bicubic_impl - Compiler-vectorized 2D Keys bicubic phase extrusion
 */
static inline void atmturb_extrude_bicubic_impl(
    const atmturb_extrude_params_t *params)
{
    const long  msize    = params->msize;
    const long  pup_size = params->pup_size;
    const long  os       = (params->os > 0) ? params->os : 1;
    const float weight   = params->weight;

    double rx = fmod(params->x0, (double) msize);
    if (rx < 0.0)
    {
        rx += (double) msize;
    }
    double ry = fmod(params->y0, (double) msize);
    if (ry < 0.0)
    {
        ry += (double) msize;
    }

    long  start_x = (long) floor(rx);
    long  start_y = (long) floor(ry);
    float fx      = (float) (rx - (double) start_x);
    float fy      = (float) (ry - (double) start_y);
    start_x %= msize;
    start_y %= msize;

    float wx_raw[4];
    float wy_raw[4];
    atmturb_keys_weights_impl(fx, wx_raw);
    atmturb_keys_weights_impl(fy, wy_raw);
    const float wx_w[4] = {
        weight * wx_raw[0], weight * wx_raw[1], weight * wx_raw[2], weight * wx_raw[3]
    };
    const long geom_i[4] = {msize, start_x, pup_size, os};

#ifdef _OPENMP
    if (pup_size >= 512 && !omp_in_parallel())
    {
        int nth = atmturb_extrude_omp_threads();
        #pragma omp parallel num_threads(nth)
        {
            int  tid = omp_get_thread_num();
            int  nt  = omp_get_num_threads();
            long j0  = (pup_size * (long) tid) / (long) nt;
            long j1  = (pup_size * (long) (tid + 1)) / (long) nt;
            atmturb_extrude_bicubic_rows(params, geom_i, wx_w, wy_raw, start_y, j0, j1);
        }
        return;
    }
#endif
    atmturb_extrude_bicubic_rows(params, geom_i, wx_w, wy_raw, start_y, 0, pup_size);
}

#endif // ATMTURB_SIMD_EXTRUDE_IMPL_H
