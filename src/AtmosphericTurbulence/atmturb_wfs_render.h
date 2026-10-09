// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_wfs_render.h
 * @brief   Wavefront series rendering drivers (geometric, diffractive, CUDA)
 */

#ifndef ATMTURB_WFS_RENDER_H
#define ATMTURB_WFS_RENDER_H

#include "CLIcore.h"
#include "atmturb_geometry.h"
#include "atmturb_profile.h"
#include "atmturb_rolling.h"
#include "ImageStreamIO/ImageStreamIO.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * struct atmturb_wfs_images_t - Output 3D image handles container
 * @ID_pha: Primary phase 3D image handle.
 * @ID_amp: Primary amplitude 3D image handle.
 * @ID_spha: Secondary phase 3D image handle.
 * @ID_samp: Secondary amplitude 3D image handle.
 */
typedef struct
{
    imageID ID_pha;
    imageID ID_amp;
    imageID ID_spha;
    imageID ID_samp;
} atmturb_wfs_images_t;

/**
 * struct atmturb_wfs_render_target_t - Phase slice destination with guard margin
 * @pup_size: Linear dimension of pupil.
 * @guard_pix: Guard band margin in pixels (0 for unpadded).
 * @pha: Destination primary phase array.
 * @spha: Destination secondary phase array.
 * @weight_scale: Fractional extrusion weight multiplier (<=0 for 1.0 default).
 */
typedef struct
{
    long   pup_size;
    long   guard_pix;
    float *pha;
    float *spha;
    float  weight_scale;
} atmturb_wfs_render_target_t;

/**
 * atmturb_wfs_render_layer_target - Render one turbulence layer into target buffer
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @k: Layer index.
 * @t: Frame index.
 * @time_step_s: Time step in seconds.
 * @master_size: Master screen dimension.
 * @target: Render target specifications and buffers.
 */
void atmturb_wfs_render_layer_target(
    const atmturb_rolling_t           *r,
    const atmturb_geom_t              *geom,
    int                                k,
    long                               t,
    double                             time_step_s,
    long                               master_size,
    const atmturb_wfs_render_target_t *target);

/**
 * atmturb_wfs_render_layer - Render one turbulence layer into phase slices for frame t
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @k: Layer index.
 * @t: Frame index.
 * @time_step_s: Time step in seconds.
 * @master_size: Master screen dimension.
 * @pup_size: Pupil dimension.
 * @pha_slice: Primary phase frame slice.
 * @spha_slice: Secondary phase frame slice.
 */
void atmturb_wfs_render_layer(
    const atmturb_rolling_t *r,
    const atmturb_geom_t    *geom,
    int                      k,
    long                     t,
    double                   time_step_s,
    long                     master_size,
    long                     pup_size,
    float                   *pha_slice,
    float                   *spha_slice);

/**
 * atmturb_wfs_render_frames - Dispatch rendering to diffractive or geometric engine
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @prof: Active turbulence profile.
 * @params: Observation parameters container.
 * @master_size: Master screen dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 * @nbframes: Number of simulation frames.
 * @imgs: Container of output 3D image handles.
 */
void atmturb_wfs_render_frames(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    const atmturb_profile_t    *prof,
    const atmturb_obs_params_t *params,
    long                        master_size,
    long                        pup_size,
    long                        nbframes,
    const atmturb_wfs_images_t *imgs);

#ifdef __cplusplus
}
#endif

#endif // ATMTURB_WFS_RENDER_H
