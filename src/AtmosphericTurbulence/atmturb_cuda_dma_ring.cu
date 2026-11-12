// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_dma_ring.cu
 * @brief   Persistent double-buffered pinned host DMA ring for 3D wavefront cubes
 */

#include <cuda_runtime.h>
#include <stddef.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "atmturb_cuda_dma_ring.h"

static float       *s_h_stage[2]    = {NULL, NULL};
static size_t       s_stage_cap     = 0;
static cudaStream_t s_dma_stream    = 0;
static cudaEvent_t  s_ev_comp[2]    = {0, 0};
static cudaEvent_t  s_ev_dma[2]     = {0, 0};
static int          s_pending[2]    = {0, 0};
static float       *s_pending_dst[2][4];
static size_t       s_pending_bytes[2] = {0, 0};

/**
 * atmturb_cuda_dma_ring_cleanup - Free persistent pinned staging buffers, events, and stream
 */
void atmturb_cuda_dma_ring_cleanup(void)
{
    if (s_dma_stream != 0)
    {
        cudaStreamSynchronize(s_dma_stream);
        cudaStreamDestroy(s_dma_stream);
        s_dma_stream = 0;
    }
    for (int s = 0; s < 2; s++)
    {
        if (s_ev_comp[s] != 0)
        {
            cudaEventDestroy(s_ev_comp[s]);
            s_ev_comp[s] = 0;
        }
        if (s_ev_dma[s] != 0)
        {
            cudaEventDestroy(s_ev_dma[s]);
            s_ev_dma[s] = 0;
        }
        if (s_h_stage[s] != NULL)
        {
            cudaFreeHost(s_h_stage[s]);
            s_h_stage[s] = NULL;
        }
        s_pending[s]       = 0;
        s_pending_bytes[s] = 0;
    }
    s_stage_cap = 0;
}

/**
 * atmturb_cuda_dma_ring_prepare - Ensure persistent 2-slot pinned D2H ring capacity
 * @max_ch_bytes: Maximum byte size of a single channel chunk per slot.
 *
 * Return: 0 on success, -1 on CUDA allocation failure.
 */
int atmturb_cuda_dma_ring_prepare(
    size_t max_ch_bytes)
{
    if (max_ch_bytes == 0)
    {
        return -1;
    }
    size_t need = 4 * max_ch_bytes;

    if (s_dma_stream == 0)
    {
        if (cudaStreamCreateWithFlags(&s_dma_stream, cudaStreamNonBlocking) != cudaSuccess)
        {
            return -1;
        }
        for (int s = 0; s < 2; s++)
        {
            if (cudaEventCreateWithFlags(&s_ev_comp[s], cudaEventDisableTiming) != cudaSuccess ||
                cudaEventCreateWithFlags(&s_ev_dma[s], cudaEventDisableTiming) != cudaSuccess)
            {
                atmturb_cuda_dma_ring_cleanup();
                return -1;
            }
        }
    }

    if (s_h_stage[0] != NULL && s_h_stage[1] != NULL && s_stage_cap >= need)
    {
        s_pending[0] = 0;
        s_pending[1] = 0;
        return 0;
    }

    for (int s = 0; s < 2; s++)
    {
        if (s_h_stage[s] != NULL)
        {
            cudaFreeHost(s_h_stage[s]);
            s_h_stage[s] = NULL;
        }
        if (cudaMallocHost((void **) &s_h_stage[s], need) != cudaSuccess)
        {
            atmturb_cuda_dma_ring_cleanup();
            return -1;
        }
        s_pending[s] = 0;
    }
    s_stage_cap = need;
    return 0;
}

/**
 * atmturb_cuda_dma_ring_drain - Synchronize and unpack one ring slot to pageable host cubes
 * @slot: Ring slot index (0 or 1).
 *
 * Return: 0 on success, -1 on CUDA synchronization error.
 */
int atmturb_cuda_dma_ring_drain(
    int slot)
{
    int s = slot & 1;
    if (!s_pending[s])
    {
        return 0;
    }
    if (cudaEventSynchronize(s_ev_dma[s]) != cudaSuccess)
    {
        s_pending[s] = 0;
        return -1;
    }

    size_t ch_bytes = s_pending_bytes[s];

#ifdef _OPENMP
    if (ch_bytes >= 262144 && !omp_in_parallel())
    {
        size_t ch_floats = ch_bytes / sizeof(float);
        int    max_t     = omp_get_max_threads();
        #pragma omp parallel num_threads((max_t > 8) ? 8 : max_t)
        {
            int    tid = omp_get_thread_num();
            int    nt  = omp_get_num_threads();
            size_t b0  = ((ch_floats * (size_t) tid) / (size_t) nt) * sizeof(float);
            size_t b1  = ((ch_floats * (size_t) (tid + 1)) / (size_t) nt) * sizeof(float);
            if (b1 > b0)
            {
                for (int c = 0; c < 4; c++)
                {
                    if (s_pending_dst[s][c] != NULL)
                    {
                        const char *src = (const char *) s_h_stage[s] + (size_t) c * ch_bytes + b0;
                        char       *dst = (char *) s_pending_dst[s][c] + b0;
                        memcpy(dst, src, b1 - b0);
                    }
                }
            }
        }
        s_pending[s] = 0;
        return 0;
    }
#endif
    for (int c = 0; c < 4; c++)
    {
        if (s_pending_dst[s][c] != NULL)
        {
            const char *src = (const char *) s_h_stage[s] + (size_t) c * ch_bytes;
            memcpy(s_pending_dst[s][c], src, ch_bytes);
        }
    }
    s_pending[s] = 0;
    return 0;
}

/**
 * atmturb_cuda_dma_ring_enqueue - Hand off device chunk from compute stream to async D2H DMA
 * @slot: Ring slot index (0 or 1); caller must drain @slot before reuse.
 * @compute_stream: Stream producing @d_src buffers.
 * @d_src: 4-element array of device channel pointers (pha, amp, spha, samp).
 * @h_dst: 4-element array of pageable host destination pointers (NULL entries skipped).
 * @ch_bytes: Byte count per active channel in this chunk.
 *
 * Return: 0 on success, -1 on CUDA error.
 */
int atmturb_cuda_dma_ring_enqueue(
    int                slot,
    cudaStream_t       compute_stream,
    const float *const d_src[4],
    float *const       h_dst[4],
    size_t             ch_bytes)
{
    int s = slot & 1;
    if (s_dma_stream == 0 || s_h_stage[s] == NULL || 4 * ch_bytes > s_stage_cap)
    {
        return -1;
    }

    if (cudaEventRecord(s_ev_comp[s], compute_stream) != cudaSuccess ||
        cudaStreamWaitEvent(s_dma_stream, s_ev_comp[s], 0) != cudaSuccess)
    {
        return -1;
    }

    for (int c = 0; c < 4; c++)
    {
        s_pending_dst[s][c] = h_dst[c];
        if (h_dst[c] != NULL && d_src[c] != NULL)
        {
            void *stage_c = (char *) s_h_stage[s] + (size_t) c * ch_bytes;
            if (cudaMemcpyAsync(stage_c, d_src[c], ch_bytes,
                                cudaMemcpyDeviceToHost, s_dma_stream) != cudaSuccess)
            {
                return -1;
            }
        }
    }

    if (cudaEventRecord(s_ev_dma[s], s_dma_stream) != cudaSuccess)
    {
        return -1;
    }
    s_pending_bytes[s] = ch_bytes;
    s_pending[s]       = 1;
    return 0;
}

/**
 * atmturb_cuda_dma_ring_flush - Drain both ring slots at batch completion
 *
 * Return: 0 on success, -1 on error.
 */
int atmturb_cuda_dma_ring_flush(void)
{
    int r0 = atmturb_cuda_dma_ring_drain(0);
    int r1 = atmturb_cuda_dma_ring_drain(1);
    return (r0 == 0 && r1 == 0) ? 0 : -1;
}
