// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_wfs_render_cuda.h
 * @brief   CUDA GPU dispatch wrappers for wavefront series rendering
 */

#ifndef ATMTURB_WFS_RENDER_CUDA_H
#define ATMTURB_WFS_RENDER_CUDA_H

#include "atmturb_geometry.h"
#include "atmturb_profile.h"
#include "atmturb_rolling.h"
#include "atmturb_wfs_render.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef HAVE_CUDA

typedef struct atmturb_cuda_rytov_stream atmturb_cuda_rytov_stream_t;

/**
 * atmturb_wfs_render_cuda - Render geometric simulation frames on GPU
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @master_size: Master screen dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 * @nbframes: Number of simulation frames.
 * @imgs: Container of output 3D image handles.
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_wfs_render_cuda(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    long                        master_size,
    long                        pup_size,
    long                        nbframes,
    const atmturb_wfs_images_t *imgs);

/**
 * atmturb_wfs_render_rytov_cuda - Render Rytov diffractive simulation frames on GPU
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @prof: Active turbulence profile.
 * @params: Observation parameters container.
 * @master_size: Master screen dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 * @nbframes: Number of simulation frames.
 * @imgs: Container of output 3D image handles.
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_wfs_render_rytov_cuda(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    const atmturb_profile_t    *prof,
    const atmturb_obs_params_t *params,
    long                        master_size,
    long                        pup_size,
    long                        nbframes,
    const atmturb_wfs_images_t *imgs);

/**
 * atmturb_cuda_rytov_stream_init - Initialize persistent GPU Rytov streaming context
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @prof: Active turbulence profile.
 * @params: Observation parameters container.
 * @master_size: Master screen dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 *
 * Return: Pointer to allocated stream context, or NULL on failure.
 */
atmturb_cuda_rytov_stream_t *atmturb_cuda_rytov_stream_init(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    const atmturb_profile_t    *prof,
    const atmturb_obs_params_t *params,
    long                        master_size,
    long                        pup_size);

/**
 * atmturb_cuda_rytov_stream_render_step - Render one stream frame on GPU
 * @ctx: Persistent stream context.
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @t: Frame index.
 * @time_step_s: Time step in seconds.
 * @pha: Destination primary phase frame buffer.
 * @amp: Destination primary amplitude frame buffer.
 * @spha: Destination secondary phase frame buffer (or NULL).
 * @samp: Destination secondary amplitude frame buffer (or NULL).
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_cuda_rytov_stream_render_step(
    atmturb_cuda_rytov_stream_t *ctx,
    const atmturb_rolling_t     *r,
    const atmturb_geom_t        *geom,
    long                         t,
    double                       time_step_s,
    float                       *pha,
    float                       *amp,
    float                       *spha,
    float                       *samp);

/**
 * atmturb_cuda_rytov_stream_free - Free persistent GPU Rytov streaming context
 * @ctx: Stream context to release.
 */
void atmturb_cuda_rytov_stream_free(
    atmturb_cuda_rytov_stream_t *ctx);

#else /* !HAVE_CUDA */

typedef void atmturb_cuda_rytov_stream_t;

static inline atmturb_cuda_rytov_stream_t *atmturb_cuda_rytov_stream_init(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    const atmturb_profile_t    *prof,
    const atmturb_obs_params_t *params,
    long                        master_size,
    long                        pup_size)
{
    (void) r;
    (void) geom;
    (void) prof;
    (void) params;
    (void) master_size;
    (void) pup_size;
    return NULL;
}

static inline int atmturb_cuda_rytov_stream_render_step(
    atmturb_cuda_rytov_stream_t *ctx,
    const atmturb_rolling_t     *r,
    const atmturb_geom_t        *geom,
    long                         t,
    double                       time_step_s,
    float                       *pha,
    float                       *amp,
    float                       *spha,
    float                       *samp)
{
    (void) ctx;
    (void) r;
    (void) geom;
    (void) t;
    (void) time_step_s;
    (void) pha;
    (void) amp;
    (void) spha;
    (void) samp;
    return -1;
}

static inline void atmturb_cuda_rytov_stream_free(
    atmturb_cuda_rytov_stream_t *ctx)
{
    (void) ctx;
}

#endif /* HAVE_CUDA */

#ifdef __cplusplus
}
#endif

#endif /* ATMTURB_WFS_RENDER_CUDA_H */
