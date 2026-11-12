// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_geom_ext.cu
 * @brief   General GPU geometric wavefront batch and stream engine (rolling, lowfreq, bicubic)
 */

#include <cuda_runtime.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

extern "C" {
#include "atmturb_cuda.h"
#include "atmturb_cuda_dma_ring.h"
#include "atmturb_cuda_geom_ext.h"
#include "atmturb_cuda_piston.h"
#include "atmturb_cuda_rytov_internal.h"
#include "atmturb_lowfreq.h"
#include "atmturb_simd.h"
#include "atmturb_superlayer.h"
#include "atmturb_wfs_cuda_common.h"
}

typedef struct
{
    cudaStream_t stream;
    char        *h_buf;
    size_t       h_cap;
    char        *d_buf;
    size_t       d_cap;
} geom_ext_cache_t;

static geom_ext_cache_t s_gbat = {0, NULL, 0, NULL, 0};

typedef struct
{
    int                           *layer_indices;
    atmturb_superlayer_t           sl;
    const float                  **masters;
    int                           *screen_offsets;
    const float                   *d_masters;
    long                           msize;
    long                           pup_size;
    long                           npix;
    int                            batch_b;
    int                            has_sec;
    int                            interp;
    int                            os;
    int                            lowfreq;
    int                            max_sub;
    int                            max_lf;
    cudaStream_t                   stream;
    atmturb_cuda_rytov_sublayer_t *h_sub[2];
    atmturb_cuda_rytov_sublayer_t *d_sub[2];
    atmturb_cuda_lowfreq_mode_t   *h_lf[2];
    atmturb_cuda_lowfreq_mode_t   *d_lf[2];
    float                         *d_pha[2];
    float                         *d_spha[2];
} geom_ext_batch_t;

struct atmturb_cuda_geom_ext_stream
{
    int                           *layer_indices;
    atmturb_superlayer_t           sl;
    const float                  **masters;
    int                           *screen_offsets;
    const float                   *d_masters;
    long                           msize;
    long                           pup_size;
    long                           npix;
    int                            interp;
    int                            os;
    int                            lowfreq;
    int                            max_sub;
    int                            max_lf;
    cudaStream_t                   streams[2];
    cudaEvent_t                    ev_done[2];
    char                          *h_buf;
    char                          *d_buf;
    atmturb_cuda_rytov_sublayer_t *h_sub[2];
    atmturb_cuda_rytov_sublayer_t *d_sub[2];
    atmturb_cuda_lowfreq_mode_t   *h_lf[2];
    atmturb_cuda_lowfreq_mode_t   *d_lf[2];
    float                         *d_pha[2];
    float                         *d_spha[2];
    float                         *d_sum[2];
};

static inline size_t align256(
    size_t n)
{
    return (n + 255) & ~(size_t) 255;
}

void atmturb_cuda_geom_ext_cleanup(void)
{
    if (s_gbat.stream != 0)
    {
        cudaStreamSynchronize(s_gbat.stream);
        cudaStreamDestroy(s_gbat.stream);
        s_gbat.stream = 0;
    }
    if (s_gbat.h_buf != NULL)
    {
        cudaFreeHost(s_gbat.h_buf);
        s_gbat.h_buf = NULL;
    }
    if (s_gbat.d_buf != NULL)
    {
        cudaFree(s_gbat.d_buf);
        s_gbat.d_buf = NULL;
    }
    s_gbat.h_cap = 0;
    s_gbat.d_cap = 0;
}

static void geom_ext_fill_amp_pair(
    float *amp,
    float *samp,
    long   npix)
{
    if (amp == NULL && samp == NULL)
    {
        return;
    }
#ifdef _OPENMP
    if (npix >= 65536 && !omp_in_parallel())
    {
        int max_t = omp_get_max_threads();
        #pragma omp parallel num_threads((max_t > 8) ? 8 : max_t)
        {
            int  tid = omp_get_thread_num();
            int  nt  = omp_get_num_threads();
            long i0  = (npix * (long) tid) / (long) nt;
            long i1  = (npix * (long) (tid + 1)) / (long) nt;
            if (i1 > i0)
            {
                if (amp != NULL)
                {
                    atmturb_init_phase_amp(NULL, amp + i0, i1 - i0);
                }
                if (samp != NULL)
                {
                    atmturb_init_phase_amp(NULL, samp + i0, i1 - i0);
                }
            }
        }
        return;
    }
#endif
    if (amp != NULL)
    {
        atmturb_init_phase_amp(NULL, amp, npix);
    }
    if (samp != NULL)
    {
        atmturb_init_phase_amp(NULL, samp, npix);
    }
}

static int geom_ext_batch_ensure_cache(
    geom_ext_batch_t *c)
{
    if (s_gbat.stream == 0 &&
        cudaStreamCreateWithFlags(&s_gbat.stream, cudaStreamNonBlocking) != cudaSuccess)
    {
        return -1;
    }
    c->stream = s_gbat.stream;

    size_t fbytes  = align256(sizeof(float) * (size_t) c->npix * (size_t) c->batch_b);
    size_t sub_byt = align256(sizeof(atmturb_cuda_rytov_sublayer_t) *
                              (size_t) (c->batch_b * c->max_sub));
    size_t lf_byt  = align256(sizeof(atmturb_cuda_lowfreq_mode_t) *
                              (size_t) (c->batch_b * c->max_lf));
    size_t h_need  = 2 * (sub_byt + lf_byt);
    size_t d_need  = 2 * (sub_byt + lf_byt + fbytes + (c->has_sec ? fbytes : 0));

    if (atmturb_cuda_dma_ring_prepare(sizeof(float) * (size_t) c->npix * (size_t) c->batch_b) != 0)
    {
        return -1;
    }
    if (s_gbat.h_cap < h_need)
    {
        cudaFreeHost(s_gbat.h_buf);
        s_gbat.h_buf = NULL;
        if (cudaMallocHost((void **) &s_gbat.h_buf, h_need) != cudaSuccess)
        {
            s_gbat.h_cap = 0;
            return -1;
        }
        s_gbat.h_cap = h_need;
    }
    if (s_gbat.d_cap < d_need)
    {
        cudaFree(s_gbat.d_buf);
        s_gbat.d_buf = NULL;
        if (cudaMalloc((void **) &s_gbat.d_buf, d_need) != cudaSuccess)
        {
            s_gbat.d_cap = 0;
            return -1;
        }
        s_gbat.d_cap = d_need;
    }

    char *hp = s_gbat.h_buf;
    char *dp = s_gbat.d_buf;
    for (int s = 0; s < 2; s++)
    {
        c->h_sub[s]  = (atmturb_cuda_rytov_sublayer_t *) hp;
        hp          += sub_byt;
        c->h_lf[s]   = (lf_byt > 0) ? (atmturb_cuda_lowfreq_mode_t *) hp : NULL;
        hp          += lf_byt;
        c->d_sub[s]  = (atmturb_cuda_rytov_sublayer_t *) dp;
        dp          += sub_byt;
        c->d_lf[s]   = (lf_byt > 0) ? (atmturb_cuda_lowfreq_mode_t *) dp : NULL;
        dp          += lf_byt;
        c->d_pha[s]  = (float *) dp;
        dp          += fbytes;
        c->d_spha[s] = c->has_sec ? (float *) dp : NULL;
        if (c->has_sec)
        {
            dp += fbytes;
        }
    }
    return 0;
}

static void geom_ext_batch_free(
    geom_ext_batch_t *c)
{
    free(c->layer_indices);
    free(c->masters);
    free(c->screen_offsets);
}

static int geom_ext_batch_init(
    geom_ext_batch_t        *c,
    const atmturb_rolling_t *r,
    const atmturb_geom_t    *geom,
    long                     master_size,
    long                     pup_size,
    long                     nbframes,
    int                      has_sec)
{
    memset(c, 0, sizeof(*c));
    c->msize    = master_size;
    c->pup_size = pup_size;
    c->npix     = pup_size * pup_size;
    c->has_sec  = has_sec;
    c->interp   = geom->interp;
    c->os       = (geom->oversample > 1) ? (int) geom->oversample : 1;
    c->lowfreq  = r->lowfreq;
    c->max_sub  = (geom->nlayers > 0) ? (2 * geom->nlayers) : 1;
    c->max_lf   = (r->lowfreq && geom->nlayers > 0) ? ATMTURB_LOWFREQ_NMODES : 0;

    int max_b  = (pup_size <= 256) ? 32 : 16;
    c->batch_b = (nbframes < (long) max_b) ? (int) nbframes : max_b;

    c->layer_indices = (int *) malloc(sizeof(int) * (size_t) geom->nlayers);
    if (c->layer_indices == NULL)
    {
        return -1;
    }
    for (int k = 0; k < geom->nlayers; k++)
    {
        c->layer_indices[k] = k;
    }
    c->sl.nlayers       = geom->nlayers;
    c->sl.layer_indices = c->layer_indices;

    int total_screens = 0;
    if (atmturb_wfs_cuda_init_screen_pool(r, geom->nlayers, &c->masters,
                                          &c->screen_offsets, &total_screens) != 0)
    {
        return -1;
    }
    c->d_masters = atmturb_cuda_sync_device_masters(master_size, total_screens, c->masters);
    if (c->d_masters == NULL)
    {
        return -1;
    }
    return geom_ext_batch_ensure_cache(c);
}

static int geom_ext_batch_run_chunk(
    geom_ext_batch_t        *c,
    const atmturb_rolling_t *r,
    const atmturb_geom_t    *geom,
    double                   time_step_s,
    int                      slot,
    long                     t0,
    int                      count,
    float                   *pha,
    float                   *spha)
{
    if (atmturb_cuda_dma_ring_drain(slot) != 0)
    {
        return -1;
    }

#ifdef _OPENMP
    #pragma omp parallel for if (count >= 4) schedule(static)
#endif
    for (int b = 0; b < count; b++)
    {
        atmturb_wfs_cuda_build_frame_sublayers(
            r, geom, 1, &c->sl, c->screen_offsets, t0 + (long) b, time_step_s, 0.0,
            c->max_sub, c->h_sub[slot] + (size_t) (b * c->max_sub));
        if (c->max_lf > 0)
        {
            atmturb_wfs_cuda_build_frame_lf_modes(
                r, geom, 1, &c->sl, t0 + (long) b, time_step_s, 0.0,
                c->max_lf, c->h_lf[slot] + (size_t) (b * c->max_lf));
        }
    }

    size_t sub_byt = sizeof(atmturb_cuda_rytov_sublayer_t) * (size_t) (count * c->max_sub);
    cudaMemcpyAsync(c->d_sub[slot], c->h_sub[slot], sub_byt, cudaMemcpyHostToDevice, c->stream);
    if (c->max_lf > 0)
    {
        size_t lf_byt = sizeof(atmturb_cuda_lowfreq_mode_t) * (size_t) (count * c->max_lf);
        cudaMemcpyAsync(c->d_lf[slot], c->h_lf[slot], lf_byt, cudaMemcpyHostToDevice, c->stream);
    }

    atmturb_cuda_extrude_batch_args_t ext = {
        .interp_mode  = c->interp,
        .d_masters    = c->d_masters,
        .msize        = (int) c->msize,
        .pad_size     = (int) c->pup_size,
        .pup_size     = (int) c->pup_size,
        .guard        = 0,
        .n_sub        = 2 * geom->nlayers,
        .sub_stride   = c->max_sub,
        .has_sec      = c->has_sec,
        .os           = c->os,
        .batch_count  = count,
        .d_sublayers  = c->d_sub[slot],
        .n_lf_modes   = c->max_lf,
        .lf_stride    = c->max_lf,
        .d_lf_modes   = c->d_lf[slot],
        .d_super_pha  = c->d_pha[slot],
        .d_super_spha = c->has_sec ? c->d_spha[slot] : NULL,
        .d_frame_pha  = NULL,
        .d_frame_spha = NULL
    };
    atmturb_cuda_extrude_sl_batch_launch(&ext, c->stream);

    size_t             dst_off  = (size_t) t0 * (size_t) c->npix;
    size_t             cbytes   = sizeof(float) * (size_t) count * (size_t) c->npix;
    const float *const d_src[4] = {
        c->d_pha[slot], NULL, c->has_sec ? c->d_spha[slot] : NULL, NULL
    };
    float *const h_dst[4] = {
        (pha != NULL) ? (pha + dst_off) : NULL, NULL,
        (c->has_sec && spha != NULL) ? (spha + dst_off) : NULL, NULL
    };
    return atmturb_cuda_dma_ring_enqueue(slot, c->stream, d_src, h_dst, cbytes);
}

int atmturb_cuda_geom_ext_render_batch(
    const atmturb_rolling_t *r,
    const atmturb_geom_t    *geom,
    long                     master_size,
    long                     pup_size,
    long                     nbframes,
    double                   time_step_s,
    float                   *pha,
    float                   *amp,
    float                   *spha,
    float                   *samp)
{
    if (r == NULL || geom == NULL || geom->nlayers <= 0 || nbframes <= 0 ||
        !atmturb_cuda_device_available())
    {
        return -1;
    }

    geom_ext_batch_t ctx;
    if (geom_ext_batch_init(&ctx, r, geom, master_size, pup_size, nbframes, spha != NULL) != 0)
    {
        geom_ext_batch_free(&ctx);
        return -1;
    }

    long tot_pix  = nbframes * ctx.npix;
    int  ret      = 0;
    int  chunk    = 0;
    int  amp_done = 0;
    for (long t0 = 0; t0 < nbframes; t0 += (long) ctx.batch_b, chunk++)
    {
        if (chunk == 2 && !amp_done)
        {
            geom_ext_fill_amp_pair(amp, samp, tot_pix);
            amp_done = 1;
        }
        long rem   = nbframes - t0;
        int  count = (rem < (long) ctx.batch_b) ? (int) rem : ctx.batch_b;
        if (geom_ext_batch_run_chunk(&ctx, r, geom, time_step_s, chunk & 1, t0, count,
                                     pha, spha) != 0)
        {
            ret = -1;
            break;
        }
    }

    if (!amp_done)
    {
        geom_ext_fill_amp_pair(amp, samp, tot_pix);
    }
    if (atmturb_cuda_dma_ring_flush() != 0)
    {
        ret = -1;
    }
    geom_ext_batch_free(&ctx);
    return ret;
}

atmturb_cuda_geom_ext_stream_t *atmturb_cuda_geom_ext_stream_init(
    const atmturb_rolling_t *r,
    const atmturb_geom_t    *geom,
    long                     master_size,
    long                     pup_size)
{
    if (r == NULL || geom == NULL || geom->nlayers <= 0 || !atmturb_cuda_device_available())
    {
        return NULL;
    }

    atmturb_cuda_geom_ext_stream_t *c =
        (atmturb_cuda_geom_ext_stream_t *) calloc(1, sizeof(*c));
    if (c == NULL)
    {
        return NULL;
    }

    c->msize         = master_size;
    c->pup_size      = pup_size;
    c->npix          = pup_size * pup_size;
    c->interp        = geom->interp;
    c->os            = (geom->oversample > 1) ? (int) geom->oversample : 1;
    c->lowfreq       = r->lowfreq;
    c->max_sub       = 2 * geom->nlayers;
    c->max_lf        = r->lowfreq ? ATMTURB_LOWFREQ_NMODES : 0;
    c->layer_indices = (int *) malloc(sizeof(int) * (size_t) geom->nlayers);
    if (c->layer_indices == NULL)
    {
        free(c);
        return NULL;
    }
    for (int k = 0; k < geom->nlayers; k++)
    {
        c->layer_indices[k] = k;
    }
    c->sl.nlayers       = geom->nlayers;
    c->sl.layer_indices = c->layer_indices;

    int total_screens = 0;
    if (atmturb_wfs_cuda_init_screen_pool(r, geom->nlayers, &c->masters,
                                          &c->screen_offsets, &total_screens) != 0)
    {
        atmturb_cuda_geom_ext_stream_free(c);
        return NULL;
    }
    c->d_masters = atmturb_cuda_sync_device_masters(master_size, total_screens, c->masters);
    if (c->d_masters == NULL)
    {
        atmturb_cuda_geom_ext_stream_free(c);
        return NULL;
    }

    size_t nbytes  = align256(sizeof(float) * (size_t) c->npix);
    size_t sub_byt = align256(sizeof(atmturb_cuda_rytov_sublayer_t) * (size_t) c->max_sub);
    size_t lf_byt  = align256(sizeof(atmturb_cuda_lowfreq_mode_t) * (size_t) c->max_lf);
    size_t h_tot   = 2 * (sub_byt + lf_byt);
    size_t d_tot   = 2 * (sub_byt + lf_byt + 2 * nbytes + 256);

    if (cudaMallocHost((void **) &c->h_buf, h_tot) != cudaSuccess ||
        cudaMalloc((void **) &c->d_buf, d_tot) != cudaSuccess)
    {
        atmturb_cuda_geom_ext_stream_free(c);
        return NULL;
    }

    char *hp = c->h_buf;
    char *dp = c->d_buf;
    for (int s = 0; s < 2; s++)
    {
        if (cudaStreamCreateWithFlags(&c->streams[s], cudaStreamNonBlocking) != cudaSuccess ||
            cudaEventCreateWithFlags(&c->ev_done[s], cudaEventDisableTiming) != cudaSuccess)
        {
            atmturb_cuda_geom_ext_stream_free(c);
            return NULL;
        }
        c->h_sub[s]  = (atmturb_cuda_rytov_sublayer_t *) hp;
        hp          += sub_byt;
        c->h_lf[s]   = (lf_byt > 0) ? (atmturb_cuda_lowfreq_mode_t *) hp : NULL;
        hp          += lf_byt;
        c->d_sub[s]  = (atmturb_cuda_rytov_sublayer_t *) dp;
        dp          += sub_byt;
        c->d_lf[s]   = (lf_byt > 0) ? (atmturb_cuda_lowfreq_mode_t *) dp : NULL;
        dp          += lf_byt;
        c->d_pha[s]  = (float *) dp;
        dp          += nbytes;
        c->d_spha[s] = (float *) dp;
        dp          += nbytes;
        c->d_sum[s]  = (float *) dp;
        dp          += 256;
    }
    return c;
}

int atmturb_cuda_geom_ext_stream_launch(
    atmturb_cuda_geom_ext_stream_t *c,
    const atmturb_rolling_t        *r,
    const atmturb_geom_t           *geom,
    long                            t,
    double                          time_step_s,
    float                          *pha,
    float                          *spha,
    int                             slot)
{
    if (c == NULL || pha == NULL)
    {
        return -1;
    }
    int          s       = slot & 1;
    cudaStream_t st      = c->streams[s];
    int          has_sec = (spha != NULL) ? 1 : 0;
    size_t       sub_byt = sizeof(atmturb_cuda_rytov_sublayer_t) * (size_t) c->max_sub;

    atmturb_wfs_cuda_build_frame_sublayers(
        r, geom, 1, &c->sl, c->screen_offsets, t, time_step_s, 0.0, c->max_sub, c->h_sub[s]);
    cudaMemcpyAsync(c->d_sub[s], c->h_sub[s], sub_byt, cudaMemcpyHostToDevice, st);

    if (c->max_lf > 0)
    {
        size_t lf_byt = sizeof(atmturb_cuda_lowfreq_mode_t) * (size_t) c->max_lf;
        atmturb_wfs_cuda_build_frame_lf_modes(
            r, geom, 1, &c->sl, t, time_step_s, 0.0, c->max_lf, c->h_lf[s]);
        cudaMemcpyAsync(c->d_lf[s], c->h_lf[s], lf_byt, cudaMemcpyHostToDevice, st);
    }

    atmturb_cuda_extrude_sl_launch(
        c->interp, c->d_masters, (int) c->msize, (int) c->pup_size, (int) c->pup_size,
        0, c->max_sub, has_sec, c->os, c->d_sub[s], c->max_lf, c->d_lf[s],
        c->d_pha[s], has_sec ? c->d_spha[s] : NULL, NULL, NULL, st);

    atmturb_cuda_remove_piston_launch(c->d_pha[s], c->d_sum[s], (int) c->npix, st);
    if (has_sec)
    {
        atmturb_cuda_remove_piston_launch(c->d_spha[s], c->d_sum[s], (int) c->npix, st);
    }

    size_t nbytes = sizeof(float) * (size_t) c->npix;
    cudaMemcpyAsync(pha, c->d_pha[s], nbytes, cudaMemcpyDeviceToHost, st);
    if (has_sec)
    {
        cudaMemcpyAsync(spha, c->d_spha[s], nbytes, cudaMemcpyDeviceToHost, st);
    }
    cudaEventRecord(c->ev_done[s], st);
    return 0;
}

int atmturb_cuda_geom_ext_stream_wait(
    atmturb_cuda_geom_ext_stream_t *c,
    int                             slot)
{
    if (c == NULL)
    {
        return -1;
    }
    return (cudaEventSynchronize(c->ev_done[slot & 1]) == cudaSuccess) ? 0 : -1;
}

void atmturb_cuda_geom_ext_stream_free(
    atmturb_cuda_geom_ext_stream_t *c)
{
    if (c == NULL)
    {
        return;
    }
    for (int s = 0; s < 2; s++)
    {
        if (c->streams[s] != 0)
        {
            cudaStreamSynchronize(c->streams[s]);
            cudaStreamDestroy(c->streams[s]);
        }
        if (c->ev_done[s] != 0)
        {
            cudaEventDestroy(c->ev_done[s]);
        }
    }
    cudaFreeHost(c->h_buf);
    cudaFree(c->d_buf);
    free(c->layer_indices);
    free(c->masters);
    free(c->screen_offsets);
    free(c);
}
