// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_screens_batch.h
 * @brief   Multi-pair CPU/GPU batch synthesis of von Karman master turbulence screens
 */

#ifndef ATMTURB_SCREENS_BATCH_H
#define ATMTURB_SCREENS_BATCH_H

#include <fftw3.h>

#include "AtmosphericTurbulence.h"
#include "atmturb_cuda_screens.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * atmturb_fill_spectrum_float_rows - Populate 2D single-precision von Karman Fourier spectrum
 * @buf: Destination FFTW single-precision complex spectrum buffer.
 * @spec: Screen generation specification.
 * @use_omp: 1 to parallelize rows via OpenMP, 0 when called inside an outer worker team.
 */
void atmturb_fill_spectrum_float_rows(
    fftwf_complex               *buf,
    const atmturb_screen_spec_t *spec,
    int                          use_omp);

/**
 * atmturb_fftwf_plan_threads - Configure FFTW single-precision planner thread count
 * @req_nthreads: Requested thread count (>0 explicit, <=0 uses omp_get_max_threads()).
 */
void atmturb_fftwf_plan_threads(
    int req_nthreads);

/**
 * atmturb_generate_screen_pairs_batch - Synthesize one or more master screen pairs via GPU or CPU
 * @req: Batch synthesis request descriptor.
 *
 * Dispatches multi-layer / large-grid FP32 requests to CUDA `cuFFT` (populating the persistent
 * GPU master screen and 2D Morton texture cache directly in VRAM when requested), with automatic
 * fallback to OpenMP-parallelized CPU synthesis sharing a single read-only FFTW3 plan and
 * per-worker scratch buffers.
 *
 * Return: 0 on success, -1 on invalid arguments or allocation failure.
 */
int atmturb_generate_screen_pairs_batch(
    const atmturb_screen_batch_req_t *req);

#ifdef __cplusplus
}
#endif

#endif /* ATMTURB_SCREENS_BATCH_H */
