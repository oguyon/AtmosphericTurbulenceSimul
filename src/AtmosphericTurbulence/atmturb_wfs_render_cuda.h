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

#endif /* HAVE_CUDA */

#ifdef __cplusplus
}
#endif

#endif /* ATMTURB_WFS_RENDER_CUDA_H */
