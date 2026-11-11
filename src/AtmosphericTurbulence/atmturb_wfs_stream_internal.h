// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_wfs_stream_internal.h
 * @brief   Internal runtime state containers for 2D shared-memory wavefront streaming
 */

#ifndef ATMTURB_WFS_STREAM_INTERNAL_H
#define ATMTURB_WFS_STREAM_INTERNAL_H

#include "atmturb_cuda.h"
#include "atmturb_cuda_fresnel.h"
#include "atmturb_cuda_rytov.h"
#include "atmturb_fresnel.h"
#include "atmturb_rytov.h"
#include "atmturb_types.h"
#include "atmturb_wfs_render_cuda.h"

/**
 * struct atmturb_stream_state_t - Active runtime state for 2D SHM wavefront streaming
 * @id_pha:       Primary phase 2D shared-memory image ID (-1 if inactive).
 * @id_amp:       Primary amplitude 2D shared-memory image ID (-1 if inactive).
 * @id_spha:      Secondary phase 2D shared-memory image ID (-1 if inactive).
 * @id_samp:      Secondary amplitude 2D shared-memory image ID (-1 if inactive).
 * @pha:          Host scratch buffer for primary phase frame (pup_size^2 floats).
 * @amp:          Host scratch buffer for primary amplitude frame (or NULL).
 * @spha:         Host scratch buffer for secondary phase frame (or NULL).
 * @samp:         Host scratch buffer for secondary amplitude frame (or NULL).
 * @fplan:        Active CPU Fresnel propagation plan (or NULL).
 * @fctx:         Active CPU Fresnel thread scratch context (or NULL).
 * @rplan:        Active CPU Rytov propagation plan (or NULL).
 * @rctx:         Active CPU Rytov thread scratch context (or NULL).
 * @cuda_fstream: Persistent GPU Fresnel streaming context (or NULL).
 * @cuda_rstream: Persistent GPU Rytov streaming context (or NULL).
 * @cuda_gstream: Persistent GPU geometric streaming context (or NULL).
 * @is_cuda:      1 when a GPU streaming backend owns frame generation, 0 for CPU.
 * @dma_pha:      Double-buffered pinned host DMA slots for primary phase.
 * @dma_amp:      Double-buffered pinned host DMA slots for primary amplitude.
 * @dma_spha:     Double-buffered pinned host DMA slots for secondary phase.
 * @dma_samp:     Double-buffered pinned host DMA slots for secondary amplitude.
 */
typedef struct
{
    imageID                        id_pha;
    imageID                        id_amp;
    imageID                        id_spha;
    imageID                        id_samp;
    float                         *pha;
    float                         *amp;
    float                         *spha;
    float                         *samp;
    const atmturb_fresnel_plan_t  *fplan;
    atmturb_fresnel_ctx_t         *fctx;
    const atmturb_rytov_plan_t    *rplan;
    atmturb_rytov_ctx_t           *rctx;
    atmturb_cuda_fresnel_stream_t *cuda_fstream;
    atmturb_cuda_rytov_stream_t   *cuda_rstream;
    atmturb_cuda_geom_stream_t    *cuda_gstream;
    int                            is_cuda;
    float                         *dma_pha[2];
    float                         *dma_amp[2];
    float                         *dma_spha[2];
    float                         *dma_samp[2];
} atmturb_stream_state_t;

/**
 * struct atmturb_stream_engines_t - Pre-allocated CPU diffractive propagation engines
 * @fplan: Fresnel split-step propagation plan.
 * @fctx:  Fresnel split-step execution scratchpad and FFTW plans.
 * @rplan: Fourier Rytov propagation plan.
 * @rctx:  Fourier Rytov execution scratchpad and FFTW plans.
 */
typedef struct
{
    atmturb_fresnel_plan_t fplan;
    atmturb_fresnel_ctx_t  fctx;
    atmturb_rytov_plan_t   rplan;
    atmturb_rytov_ctx_t    rctx;
} atmturb_stream_engines_t;

#endif // ATMTURB_WFS_STREAM_INTERNAL_H
