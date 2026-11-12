// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_screens_batch.c
 * @brief   Multi-pair CPU/GPU batch synthesis of von Karman master turbulence screens
 */

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

static long           s_single_msize = 0;
static fftwf_plan     s_single_plan  = NULL;
static fftwf_complex *s_single_buf   = NULL;

/**
 * atmturb_screens_batch_cleanup_single - Release cached single-pair CPU FFTW plan and buffer
 */
__attribute__((destructor)) static void atmturb_screens_batch_cleanup_single(void)
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
}

/**
 * atmturb_deinterleave_complex_pair - Extract real and optional imaginary screens from DFT output
 * @buf: Source complex spatial-domain buffer (`total` elements).
 * @out_a: Destination real screen buffer (or NULL).
 * @out_b: Destination imaginary screen buffer (or NULL).
 * @total: Number of pixels (`size * size`).
 * @use_omp: 1 to parallelize extraction via OpenMP, 0 inside an outer worker team.
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
    }
    else if (out_a != NULL)
    {
#ifdef _OPENMP
        #pragma omp parallel for if (use_omp) schedule(static)
#endif
        for (long i = 0; i < total; i++)
        {
            out_a[i] = buf[i][0];
        }
    }
    else if (out_b != NULL)
    {
#ifdef _OPENMP
        #pragma omp parallel for if (use_omp) schedule(static)
#endif
        for (long i = 0; i < total; i++)
        {
            out_b[i] = buf[i][1];
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
        atmturb_screens_batch_cleanup_single();
        s_single_buf = (fftwf_complex *) fftwf_alloc_complex(total);
        if (s_single_buf == NULL)
        {
            return -1;
        }
        atmturb_fftwf_wisdom_load();
        atmturb_fftwf_plan_threads(0);
        s_single_plan = fftwf_plan_dft_2d((int) msize, (int) msize, s_single_buf, s_single_buf,
                                          FFTW_FORWARD, FFTW_ESTIMATE);
        if (s_single_plan == NULL)
        {
            atmturb_screens_batch_cleanup_single();
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
 * atmturb_generate_multi_pairs_cpu - Parallelize multiple FP32 screen pairs across OpenMP workers
 */
static int atmturb_generate_multi_pairs_cpu(
    const atmturb_screen_batch_req_t *req,
    long                              msize)
{
    long total     = msize * msize;
    int  max_t     = 1;
#ifdef _OPENMP
    max_t = omp_get_max_threads();
#endif
    int n_workers = (req->npairs < max_t) ? req->npairs : max_t;
    if (msize >= 4096 && n_workers > 4)
    {
        n_workers = 4;
    }
    if (msize >= 2048 && n_workers > 8)
    {
        n_workers = 8;
    }
    int fftw_nth = max_t / n_workers;
    if (fftw_nth < 1)
    {
        fftw_nth = 1;
    }
    if (fftw_nth > 8)
    {
        fftw_nth = 8;
    }

    fftwf_complex **bufs = (fftwf_complex **) calloc((size_t) n_workers, sizeof(fftwf_complex *));
    if (bufs == NULL)
    {
        return -1;
    }
    for (int w = 0; w < n_workers; w++)
    {
        bufs[w] = (fftwf_complex *) fftwf_alloc_complex(total);
        if (bufs[w] == NULL)
        {
            for (int j = 0; j < w; j++)
            {
                fftwf_free(bufs[j]);
            }
            free(bufs);
            return -1;
        }
    }

    atmturb_fftwf_wisdom_load();
    atmturb_fftwf_plan_threads(fftw_nth);
    fftwf_plan plan = fftwf_plan_dft_2d((int) msize, (int) msize, bufs[0], bufs[0],
                                        FFTW_FORWARD, FFTW_ESTIMATE);
    if (plan == NULL)
    {
        for (int w = 0; w < n_workers; w++)
        {
            fftwf_free(bufs[w]);
        }
        free(bufs);
        return -1;
    }

#ifdef _OPENMP
    #pragma omp parallel for num_threads(n_workers) schedule(dynamic, 1)
#endif
    for (int p = 0; p < req->npairs; p++)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        fftwf_complex *wbuf = bufs[tid];
        atmturb_fill_spectrum_float_rows(wbuf, &req->specs[p], 0);
        fftwf_execute_dft(plan, wbuf, wbuf);
        float *ob = (req->out_b != NULL) ? req->out_b[p] : NULL;
        atmturb_deinterleave_complex_pair(wbuf, req->out_a[p], ob, total, 0);
    }

    fftwf_destroy_plan(plan);
    for (int w = 0; w < n_workers; w++)
    {
        fftwf_free(bufs[w]);
    }
    free(bufs);
    return 0;
}

int atmturb_generate_screen_pairs_batch(
    const atmturb_screen_batch_req_t *req)
{
    if (req == NULL || req->specs == NULL || req->out_a == NULL || req->npairs <= 0)
    {
        return -1;
    }

    long msize   = req->specs[0].size;
    int  uniform = (req->specs[0].precision == 0 && msize > 0) ? 1 : 0;
    for (int p = 1; p < req->npairs && uniform; p++)
    {
        if (req->specs[p].size != msize || req->specs[p].precision != 0)
        {
            uniform = 0;
        }
    }

    if (!uniform)
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
