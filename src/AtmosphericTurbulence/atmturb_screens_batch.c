// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_screens_batch.c
 * @brief   Multi-pair CPU/GPU batch synthesis of von Karman master turbulence screens
 */

#include <math.h>
#include <stdlib.h>
#include <fftw3.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "AtmosphericTurbulence.h"
#include "atmturb_cuda_screens.h"
#include "atmturb_fft_wisdom.h"
#include "atmturb_screens_batch.h"
#include "atmturb_simd.h"
#include "atmturb_types.h"

static long            s_single_msize = 0;
static fftwf_plan      s_single_plan  = NULL;
static fftwf_complex  *s_single_buf   = NULL;

static long            s_multi_msize  = 0;
static int             s_multi_nw     = 0;
static fftwf_plan      s_multi_plan   = NULL;
static fftwf_complex **s_multi_bufs   = NULL;

static long            s_mdbl_msize   = 0;
static int             s_mdbl_nw      = 0;
static fftw_plan       s_mdbl_plan    = NULL;
static fftw_complex  **s_mdbl_bufs    = NULL;

/**
 * atmturb_screens_free_multi_f - Release cached multi-worker FP32 FFTW buffers and plan
 */
static void atmturb_screens_free_multi_f(void)
{
    if (s_multi_plan != NULL)
    {
        fftwf_destroy_plan(s_multi_plan);
        s_multi_plan = NULL;
    }
    if (s_multi_bufs != NULL)
    {
        for (int w = 0; w < s_multi_nw; w++)
        {
            if (s_multi_bufs[w] != NULL)
            {
                fftwf_free(s_multi_bufs[w]);
            }
        }
        free(s_multi_bufs);
        s_multi_bufs = NULL;
    }
    s_multi_msize = 0;
    s_multi_nw    = 0;
}

/**
 * atmturb_screens_free_multi_d - Release cached multi-worker FP64 FFTW buffers and plan
 */
static void atmturb_screens_free_multi_d(void)
{
    if (s_mdbl_plan != NULL)
    {
        fftw_destroy_plan(s_mdbl_plan);
        s_mdbl_plan = NULL;
    }
    if (s_mdbl_bufs != NULL)
    {
        for (int w = 0; w < s_mdbl_nw; w++)
        {
            if (s_mdbl_bufs[w] != NULL)
            {
                fftw_free(s_mdbl_bufs[w]);
            }
        }
        free(s_mdbl_bufs);
        s_mdbl_bufs = NULL;
    }
    s_mdbl_msize = 0;
    s_mdbl_nw    = 0;
}

/**
 * atmturb_screens_batch_cleanup - Release all cached CPU FFTW plans and buffers at exit
 */
__attribute__((destructor)) static void atmturb_screens_batch_cleanup(void)
{
    if (s_single_plan != NULL)
    {
        fftwf_destroy_plan(s_single_plan);
        s_single_plan = NULL;
    }
    if (s_single_buf != NULL)
    {
        fftwf_free(s_single_buf);
        s_single_buf = NULL;
    }
    s_single_msize = 0;
    atmturb_screens_free_multi_f();
    atmturb_screens_free_multi_d();
}

/**
 * atmturb_deinterleave_complex_pair - Extract real and optional imaginary screens from FP32 DFT
 */
static void atmturb_deinterleave_complex_pair(
    const fftwf_complex *restrict buf,
    float               *restrict out_a,
    float               *restrict out_b,
    long                          total,
    int                           use_omp)
{
    if (out_a != NULL && out_b != NULL)
    {
#ifdef _OPENMP
        #pragma omp parallel for if (use_omp) schedule(static)
#else
        (void) use_omp;
#endif
        for (long i = 0; i < total; i++)
        {
            out_a[i] = buf[i][0];
            out_b[i] = buf[i][1];
        }
        return;
    }

    float *restrict dst = (out_a != NULL) ? out_a : out_b;
    int             ch  = (out_a != NULL) ? 0 : 1;
    if (dst != NULL)
    {
#ifdef _OPENMP
        #pragma omp parallel for if (use_omp) schedule(static)
#endif
        for (long i = 0; i < total; i++)
        {
            dst[i] = buf[i][ch];
        }
    }
}

/**
 * atmturb_deinterleave_double_pair - Extract real and optional imaginary screens from FP64 DFT
 */
static void atmturb_deinterleave_double_pair(
    const fftw_complex *restrict buf,
    float              *restrict out_a,
    float              *restrict out_b,
    long                         total,
    int                          use_omp)
{
    if (out_a != NULL && out_b != NULL)
    {
#ifdef _OPENMP
        #pragma omp parallel for if (use_omp) schedule(static)
#else
        (void) use_omp;
#endif
        for (long i = 0; i < total; i++)
        {
            out_a[i] = (float) buf[i][0];
            out_b[i] = (float) buf[i][1];
        }
        return;
    }

    float *restrict dst = (out_a != NULL) ? out_a : out_b;
    int             ch  = (out_a != NULL) ? 0 : 1;
    if (dst != NULL)
    {
#ifdef _OPENMP
        #pragma omp parallel for if (use_omp) schedule(static)
#endif
        for (long i = 0; i < total; i++)
        {
            dst[i] = (float) buf[i][ch];
        }
    }
}

/**
 * atmturb_generate_single_pair_cpu_cached - Execute single FP32 screen pair with cached FFTW plan
 */
static int atmturb_generate_single_pair_cpu_cached(
    const atmturb_screen_spec_t *spec,
    float                       *out_a,
    float                       *out_b)
{
    long msize = spec->size;
    long total = msize * msize;

#ifdef _OPENMP
    if (omp_in_parallel())
    {
        fftwf_complex *buf = (fftwf_complex *) fftwf_alloc_complex(total);
        if (buf == NULL)
        {
            return -1;
        }
        atmturb_fill_spectrum_float_rows(buf, spec, 0);
        fftwf_plan plan = NULL;
        #pragma omp critical
        {
            atmturb_fftwf_wisdom_load();
            atmturb_fftwf_plan_threads(1);
            plan = fftwf_plan_dft_2d((int) msize, (int) msize, buf, buf,
                                     FFTW_FORWARD, FFTW_ESTIMATE);
        }
        fftwf_execute(plan);
        #pragma omp critical
        {
            fftwf_destroy_plan(plan);
        }
        atmturb_deinterleave_complex_pair(buf, out_a, out_b, total, 0);
        fftwf_free(buf);
        return 0;
    }
#endif

    if (s_single_msize != msize || s_single_plan == NULL || s_single_buf == NULL)
    {
        if (s_single_plan != NULL)
        {
            fftwf_destroy_plan(s_single_plan);
            s_single_plan = NULL;
        }
        if (s_single_buf != NULL)
        {
            fftwf_free(s_single_buf);
        }
        s_single_buf = (fftwf_complex *) fftwf_alloc_complex(total);
        if (s_single_buf == NULL)
        {
            return -1;
        }
        atmturb_fftwf_wisdom_load();
        atmturb_fftwf_plan_threads(0);
        s_single_plan = fftwf_plan_dft_2d((int) msize, (int) msize, s_single_buf, s_single_buf,
                                           FFTW_FORWARD, atmturb_fftwf_plan_flags());
        atmturb_fftwf_wisdom_save();
        if (s_single_plan == NULL)
        {
            fftwf_free(s_single_buf);
            s_single_buf = NULL;
            return -1;
        }
        s_single_msize = msize;
    }

    atmturb_fill_spectrum_float_rows(s_single_buf, spec, 1);
    fftwf_execute(s_single_plan);
    atmturb_deinterleave_complex_pair(s_single_buf, out_a, out_b, total, 1);
    return 0;
}

/**
 * atmturb_screens_pick_workers - Compute worker slot count and inner FFTW thread count
 */
static void atmturb_screens_pick_workers(
    int  npairs,
    long msize,
    int *out_nw,
    int *out_fftw_nth)
{
    int max_t = 1;
#ifdef _OPENMP
    max_t = omp_get_max_threads();
#endif
    int nw = (npairs < max_t) ? npairs : max_t;
    if (msize >= 4096 && nw > 4)
    {
        nw = 4;
    }
    else if (msize >= 2048 && nw > 8)
    {
        nw = 8;
    }
    else if (nw > 16)
    {
        nw = 16;
    }
    if (nw < 1)
    {
        nw = 1;
    }
    *out_nw       = nw;
    *out_fftw_nth = (nw == 1) ? ((max_t > 8) ? 8 : max_t) : 1;
}

/**
 * atmturb_ensure_multi_cache_f - Ensure persistent FP32 worker buffers and wisdom plan
 */
static int atmturb_ensure_multi_cache_f(
    long msize,
    int  nw,
    int  fftw_nth)
{
    if (s_multi_msize == msize && s_multi_nw >= nw && s_multi_plan != NULL && s_multi_bufs != NULL)
    {
        return 0;
    }

    atmturb_screens_free_multi_f();
    long total   = msize * msize;
    s_multi_bufs = (fftwf_complex **) calloc((size_t) nw, sizeof(fftwf_complex *));
    if (s_multi_bufs == NULL)
    {
        return -1;
    }
    s_multi_nw = nw;
    for (int w = 0; w < nw; w++)
    {
        s_multi_bufs[w] = (fftwf_complex *) fftwf_alloc_complex(total);
        if (s_multi_bufs[w] == NULL)
        {
            atmturb_screens_free_multi_f();
            return -1;
        }
    }

    atmturb_fftwf_wisdom_load();
    atmturb_fftwf_plan_threads(fftw_nth);
    s_multi_plan = fftwf_plan_dft_2d((int) msize, (int) msize, s_multi_bufs[0], s_multi_bufs[0],
                                      FFTW_FORWARD, atmturb_fftwf_plan_flags());
    atmturb_fftwf_wisdom_save();
    if (s_multi_plan == NULL)
    {
        atmturb_screens_free_multi_f();
        return -1;
    }
    s_multi_msize = msize;
    return 0;
}

/**
 * atmturb_generate_multi_pairs_cpu - Parallelize FP32 screen pairs with all-thread spectrum fill
 */
static int atmturb_generate_multi_pairs_cpu(
    const atmturb_screen_batch_req_t *req,
    long                              msize)
{
    int nw       = 1;
    int fftw_nth = 1;
    atmturb_screens_pick_workers(req->npairs, msize, &nw, &fftw_nth);
    if (atmturb_ensure_multi_cache_f(msize, nw, fftw_nth) != 0)
    {
        return -1;
    }

    long        total    = msize * msize;
    const float km_scale = 5.92f / (2.0f * (float) M_PI);

    for (int p0 = 0; p0 < req->npairs; p0 += nw)
    {
        int chunk = req->npairs - p0;
        if (chunk > nw)
        {
            chunk = nw;
        }

        float    bp[16][3];
        uint64_t seeds[16];
        for (int c = 0; c < chunk; c++)
        {
            const atmturb_screen_spec_t *sp = &req->specs[p0 + c];
            float r0   = (sp->r0_pix > 0.0) ? (float) sp->r0_pix : powf(6.88f, 0.6f);
            float k0   = (sp->L0_pix > 0.0) ? ((float) msize / (float) sp->L0_pix) : 0.0f;
            float km   = (sp->l0_pix > 0.0) ? (km_scale * (float) msize / (float) sp->l0_pix)
                                            : 0.0f;
            float pref = 0.023f * powf(r0, -5.0f / 3.0f) * powf((float) msize, 5.0f / 3.0f);
            bp[c][0]   = k0 * k0;
            bp[c][1]   = sqrtf(pref);
            bp[c][2]   = (km > 0.0f) ? (0.5f / (km * km)) : 0.0f;
            seeds[c]   = atmturb_resolve_seed(sp->seed);
        }

#ifdef _OPENMP
        #pragma omp parallel for collapse(2) schedule(dynamic, 16)
#endif
        for (int c = 0; c < chunk; c++)
        {
            for (long jj = 0; jj < msize; jj++)
            {
                atmturb_fill_spectrum_row_float(&s_multi_bufs[c][jj * msize], msize, jj,
                                                seeds[c], bp[c]);
            }
        }

#ifdef _OPENMP
        #pragma omp parallel for num_threads(chunk) schedule(dynamic, 1)
#endif
        for (int c = 0; c < chunk; c++)
        {
            fftwf_execute_dft(s_multi_plan, s_multi_bufs[c], s_multi_bufs[c]);
        }

#ifdef _OPENMP
        #pragma omp parallel for collapse(2) schedule(static)
#endif
        for (int c = 0; c < chunk; c++)
        {
            for (long blk = 0; blk < total; blk += 4096)
            {
                long   len = (blk + 4096 <= total) ? 4096 : (total - blk);
                float *oa  = (req->out_a[p0 + c] != NULL) ? &req->out_a[p0 + c][blk] : NULL;
                float *ob  = (req->out_b != NULL && req->out_b[p0 + c] != NULL)
                                 ? &req->out_b[p0 + c][blk]
                                 : NULL;
                atmturb_deinterleave_complex_pair(&s_multi_bufs[c][blk], oa, ob, len, 0);
            }
        }
    }
    return 0;
}

/**
 * atmturb_ensure_multi_cache_d - Ensure persistent FP64 worker buffers and FFTW plan
 */
static int atmturb_ensure_multi_cache_d(
    long msize,
    int  nw,
    int  fftw_nth)
{
    if (s_mdbl_msize == msize && s_mdbl_nw == nw && s_mdbl_plan != NULL && s_mdbl_bufs != NULL)
    {
        return 0;
    }

    atmturb_screens_free_multi_d();
    long total  = msize * msize;
    s_mdbl_bufs = (fftw_complex **) calloc((size_t) nw, sizeof(fftw_complex *));
    if (s_mdbl_bufs == NULL)
    {
        return -1;
    }
    s_mdbl_nw = nw;
    for (int w = 0; w < nw; w++)
    {
        s_mdbl_bufs[w] = (fftw_complex *) fftw_alloc_complex(total);
        if (s_mdbl_bufs[w] == NULL)
        {
            atmturb_screens_free_multi_d();
            return -1;
        }
    }

    atmturb_fftw_plan_threads(fftw_nth);
    s_mdbl_plan = fftw_plan_dft_2d((int) msize, (int) msize, s_mdbl_bufs[0], s_mdbl_bufs[0],
                                    FFTW_FORWARD, FFTW_ESTIMATE);
    if (s_mdbl_plan == NULL)
    {
        atmturb_screens_free_multi_d();
        return -1;
    }
    s_mdbl_msize = msize;
    return 0;
}

/**
 * atmturb_generate_pairs_double_cpu - Parallelize FP64 screen pairs across OpenMP workers
 */
static int atmturb_generate_pairs_double_cpu(
    const atmturb_screen_batch_req_t *req,
    long                              msize)
{
    int nw       = 1;
    int fftw_nth = 1;
    atmturb_screens_pick_workers(req->npairs, msize, &nw, &fftw_nth);
    if (atmturb_ensure_multi_cache_d(msize, nw, fftw_nth) != 0)
    {
        return -1;
    }

    long         total    = msize * msize;
    const double km_scale = 5.92 / (2.0 * M_PI);

    for (int p0 = 0; p0 < req->npairs; p0 += nw)
    {
        int chunk = req->npairs - p0;
        if (chunk > nw)
        {
            chunk = nw;
        }

        double   bp[16][3];
        uint64_t seeds[16];
        for (int c = 0; c < chunk; c++)
        {
            const atmturb_screen_spec_t *sp = &req->specs[p0 + c];
            double r0   = (sp->r0_pix > 0.0) ? sp->r0_pix : pow(6.88, 0.6);
            double k0   = (sp->L0_pix > 0.0) ? ((double) msize / sp->L0_pix) : 0.0;
            double km   = (sp->l0_pix > 0.0) ? (km_scale * (double) msize / sp->l0_pix) : 0.0;
            double pref = 0.023 * pow(r0, -5.0 / 3.0) * pow((double) msize, 5.0 / 3.0);
            bp[c][0]    = k0 * k0;
            bp[c][1]    = sqrt(pref);
            bp[c][2]    = (km > 0.0) ? (0.5 / (km * km)) : 0.0;
            seeds[c]    = atmturb_resolve_seed(sp->seed);
        }

#ifdef _OPENMP
        #pragma omp parallel for collapse(2) schedule(dynamic, 16)
#endif
        for (int c = 0; c < chunk; c++)
        {
            for (long jj = 0; jj < msize; jj++)
            {
                atmturb_fill_spectrum_row_double(&s_mdbl_bufs[c][jj * msize], msize, jj,
                                                 seeds[c], bp[c]);
            }
        }

        if (chunk == 1)
        {
            fftw_execute_dft(s_mdbl_plan, s_mdbl_bufs[0], s_mdbl_bufs[0]);
        }
        else
        {
#ifdef _OPENMP
            #pragma omp parallel for num_threads(chunk) schedule(dynamic, 1)
#endif
            for (int c = 0; c < chunk; c++)
            {
                fftw_execute_dft(s_mdbl_plan, s_mdbl_bufs[c], s_mdbl_bufs[c]);
            }
        }

#ifdef _OPENMP
        #pragma omp parallel for collapse(2) schedule(static)
#endif
        for (int c = 0; c < chunk; c++)
        {
            for (long blk = 0; blk < total; blk += 4096)
            {
                long   len = (blk + 4096 <= total) ? 4096 : (total - blk);
                float *oa  = (req->out_a[p0 + c] != NULL) ? &req->out_a[p0 + c][blk] : NULL;
                float *ob  = (req->out_b != NULL && req->out_b[p0 + c] != NULL)
                                 ? &req->out_b[p0 + c][blk]
                                 : NULL;
                atmturb_deinterleave_double_pair(&s_mdbl_bufs[c][blk], oa, ob, len, 0);
            }
        }
    }
    return 0;
}

/**
 * atmturb_generate_screen_pairs_batch - Synthesize one or more master screen pairs via GPU or CPU
 * @req: Batch synthesis request descriptor.
 *
 * Return: 0 on success, -1 on invalid arguments or allocation failure.
 */
int atmturb_generate_screen_pairs_batch(
    const atmturb_screen_batch_req_t *req)
{
    if (req == NULL || req->specs == NULL || req->out_a == NULL || req->npairs <= 0)
    {
        return -1;
    }

    long msize = req->specs[0].size;
    int  prec  = req->specs[0].precision;
    int  unif  = (msize > 0) ? 1 : 0;
    for (int p = 1; p < req->npairs && unif; p++)
    {
        if (req->specs[p].size != msize || req->specs[p].precision != prec)
        {
            unif = 0;
        }
    }

    if (!unif)
    {
        for (int p = 0; p < req->npairs; p++)
        {
            float *ob = (req->out_b != NULL) ? req->out_b[p] : NULL;
            if (atmturb_generate_screen_pair(&req->specs[p], req->out_a[p], ob) != 0)
            {
                return -1;
            }
        }
        return 0;
    }

    if (prec != 0)
    {
        return atmturb_generate_pairs_double_cpu(req, msize);
    }

    if (atmturb_simd_is_gpu() &&
        (req->npairs >= 2 || msize >= 2048 || atmturb_simd_is_gpu_forced()))
    {
        if (atmturb_cuda_generate_screen_pairs(req) == 0)
        {
            return 0;
        }
    }

    if (req->npairs == 1)
    {
        float *ob = (req->out_b != NULL) ? req->out_b[0] : NULL;
        return atmturb_generate_single_pair_cpu_cached(&req->specs[0], req->out_a[0], ob);
    }
    return atmturb_generate_multi_pairs_cpu(req, msize);
}
