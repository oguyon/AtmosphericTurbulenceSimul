// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_dma_ring.h
 * @brief   Persistent double-buffered pinned host DMA ring for 3D wavefront cubes
 */

#ifndef ATMTURB_CUDA_DMA_RING_H
#define ATMTURB_CUDA_DMA_RING_H

#include <cuda_runtime.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * atmturb_cuda_dma_ring_prepare - Ensure persistent 2-slot pinned D2H ring capacity
 * @max_ch_bytes: Maximum byte size of a single channel chunk per slot.
 *
 * Return: 0 on success, -1 on CUDA allocation failure.
 */
int atmturb_cuda_dma_ring_prepare(
    size_t max_ch_bytes);

/**
 * atmturb_cuda_dma_ring_drain - Synchronize and unpack one ring slot to pageable host cubes
 * @slot: Ring slot index (0 or 1).
 *
 * Return: 0 on success, -1 on CUDA synchronization error.
 */
int atmturb_cuda_dma_ring_drain(
    int slot);

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
    size_t             ch_bytes);

/**
 * atmturb_cuda_dma_ring_flush - Drain both ring slots at batch completion
 *
 * Return: 0 on success, -1 on error.
 */
int atmturb_cuda_dma_ring_flush(void);

/**
 * atmturb_cuda_dma_ring_cleanup - Free persistent pinned staging buffers, events, and stream
 */
void atmturb_cuda_dma_ring_cleanup(void);

#ifdef __cplusplus
}
#endif

#endif /* ATMTURB_CUDA_DMA_RING_H */
