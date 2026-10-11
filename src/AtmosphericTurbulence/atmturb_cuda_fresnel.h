// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_fresnel.h
 * @brief   CUDA and cuFFT accelerated multi-layer Fresnel diffractive propagation
 */

#ifndef ATMTURB_CUDA_FRESNEL_H
#define ATMTURB_CUDA_FRESNEL_H

#include "atmturb_fresnel.h"
#include "atmturb_geometry.h"
#include "atmturb_profile.h"
#include "atmturb_rolling.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef HAVE_CUDA

typedef void atmturb_cuda_fresnel_stream_t;

static inline atmturb_cuda_fresnel_stream_t *atmturb_cuda_fresnel_stream_init(
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

static inline int atmturb_cuda_fresnel_stream_render_step(
    atmturb_cuda_fresnel_stream_t *ctx,
    const atmturb_rolling_t       *r,
    const atmturb_geom_t          *geom,
    long                           t,
    double                         time_step_s,
    float                         *pha,
    float                         *amp,
    float                         *spha,
    float                         *samp)
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

static inline int atmturb_cuda_fresnel_stream_render_launch(
    atmturb_cuda_fresnel_stream_t *ctx,
    const atmturb_rolling_t       *r,
    const atmturb_geom_t          *geom,
    long                           t,
    double                         time_step_s,
    float                         *pha,
    float                         *amp,
    float                         *spha,
    float                         *samp)
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

static inline int atmturb_cuda_fresnel_stream_render_wait(
    atmturb_cuda_fresnel_stream_t *ctx)
{
    (void) ctx;
    return -1;
}

static inline void atmturb_cuda_fresnel_stream_free(
    atmturb_cuda_fresnel_stream_t *ctx)
{
    (void) ctx;
}

#else

typedef struct atmturb_cuda_fresnel_stream atmturb_cuda_fresnel_stream_t;

/**
 * atmturb_cuda_fresnel_stream_init - Initialize persistent GPU Fresnel streaming context
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @prof: Active turbulence profile.
 * @params: Observation parameters container.
 * @master_size: Master screen linear dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 *
 * Return: Allocated streaming context, or NULL on failure.
 */
atmturb_cuda_fresnel_stream_t *atmturb_cuda_fresnel_stream_init(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    const atmturb_profile_t    *prof,
    const atmturb_obs_params_t *params,
    long                        master_size,
    long                        pup_size);

/**
 * atmturb_cuda_fresnel_stream_render_step - Render one Fresnel stream frame on GPU
 * @ctx: Persistent stream context.
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @t: Frame index.
 * @time_step_s: Time step in seconds.
 * @pha: Destination primary phase frame buffer.
 * @amp: Destination primary amplitude frame buffer (or NULL).
 * @spha: Destination secondary phase frame buffer (or NULL).
 * @samp: Destination secondary amplitude frame buffer (or NULL).
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_cuda_fresnel_stream_render_step(
    atmturb_cuda_fresnel_stream_t *ctx,
    const atmturb_rolling_t       *r,
    const atmturb_geom_t          *geom,
    long                           t,
    double                         time_step_s,
    float                         *pha,
    float                         *amp,
    float                         *spha,
    float                         *samp);

/**
 * atmturb_cuda_fresnel_stream_render_launch - Enqueue GPU ops for Fresnel frame asynchronously
 * @ctx: Persistent stream context.
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @t: Frame index.
 * @time_step_s: Time step in seconds.
 * @pha: Destination primary phase frame buffer.
 * @amp: Destination primary amplitude frame buffer (or NULL).
 * @spha: Destination secondary phase frame buffer (or NULL).
 * @samp: Destination secondary amplitude frame buffer (or NULL).
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_cuda_fresnel_stream_render_launch(
    atmturb_cuda_fresnel_stream_t *ctx,
    const atmturb_rolling_t       *r,
    const atmturb_geom_t          *geom,
    long                           t,
    double                         time_step_s,
    float                         *pha,
    float                         *amp,
    float                         *spha,
    float                         *samp);

/**
 * atmturb_cuda_fresnel_stream_render_wait - Wait for enqueued Fresnel frame to finish on GPU
 * @ctx: Persistent stream context.
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_cuda_fresnel_stream_render_wait(
    atmturb_cuda_fresnel_stream_t *ctx);

/**
 * atmturb_cuda_fresnel_stream_free - Free persistent GPU Fresnel streaming context
 * @ctx: Stream context to release.
 */
void atmturb_cuda_fresnel_stream_free(
    atmturb_cuda_fresnel_stream_t *ctx);

#endif /* HAVE_CUDA */

#ifdef __cplusplus
}
#endif

#endif /* ATMTURB_CUDA_FRESNEL_H */
