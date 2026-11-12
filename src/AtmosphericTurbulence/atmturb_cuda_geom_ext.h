// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_geom_ext.h
 * @brief   General GPU geometric wavefront batch and stream engine (rolling, lowfreq, bicubic)
 */

#ifndef ATMTURB_CUDA_GEOM_EXT_H
#define ATMTURB_CUDA_GEOM_EXT_H

#include "atmturb_geometry.h"
#include "atmturb_rolling.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef HAVE_CUDA

typedef struct atmturb_cuda_geom_ext_stream atmturb_cuda_geom_ext_stream_t;

/**
 * atmturb_cuda_geom_ext_render_batch - Render 3D geometric cube on GPU with rolling/lowfreq/bicubic
 */
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
    float                   *samp);

/**
 * atmturb_cuda_geom_ext_stream_init - Allocate general GPU geometric 2D streaming context
 */
atmturb_cuda_geom_ext_stream_t *atmturb_cuda_geom_ext_stream_init(
    const atmturb_rolling_t *r,
    const atmturb_geom_t    *geom,
    long                     master_size,
    long                     pup_size);

/**
 * atmturb_cuda_geom_ext_stream_launch - Asynchronously dispatch one general geometric frame
 */
int atmturb_cuda_geom_ext_stream_launch(
    atmturb_cuda_geom_ext_stream_t *ctx,
    const atmturb_rolling_t        *r,
    const atmturb_geom_t           *geom,
    long                            t,
    double                          time_step_s,
    float                          *pha,
    float                          *spha,
    int                             slot);

/**
 * atmturb_cuda_geom_ext_stream_wait - Synchronize on an in-flight general geometric frame slot
 */
int atmturb_cuda_geom_ext_stream_wait(
    atmturb_cuda_geom_ext_stream_t *ctx,
    int                             slot);

/**
 * atmturb_cuda_geom_ext_stream_free - Release general GPU geometric streaming context
 */
void atmturb_cuda_geom_ext_stream_free(
    atmturb_cuda_geom_ext_stream_t *ctx);

/**
 * atmturb_cuda_geom_ext_cleanup - Release cached batch workspace buffers and stream
 */
void atmturb_cuda_geom_ext_cleanup(void);

#endif /* HAVE_CUDA */

#ifdef __cplusplus
}
#endif

#endif /* ATMTURB_CUDA_GEOM_EXT_H */
