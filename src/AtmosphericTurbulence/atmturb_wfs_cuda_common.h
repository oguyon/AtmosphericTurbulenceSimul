// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_wfs_cuda_common.h
 * @brief   Shared helper declarations for CUDA wavefront rendering and streaming
 */

#ifndef ATMTURB_WFS_CUDA_COMMON_H
#define ATMTURB_WFS_CUDA_COMMON_H

#include "atmturb_cuda_rytov.h"
#include "atmturb_cuda_rytov_internal.h"
#include "atmturb_geometry.h"
#include "atmturb_rolling.h"
#include "atmturb_superlayer.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * atmturb_wfs_cuda_init_screen_pool - Assemble pointers to master phase screens
 * @r: Rolling simulation context.
 * @nlayers: Number of simulation layers.
 * @out_masters: Pointer to store allocated screen pointer array.
 * @out_offsets: Pointer to store allocated layer screen offsets.
 * @out_total: Pointer to store total screen count.
 *
 * Flattens master phase screen arrays across all turbulent layers into a contiguous pointer list
 * for upload and caching in GPU device memory.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_wfs_cuda_init_screen_pool(
    const atmturb_rolling_t  *r,
    int                       nlayers,
    const float            ***out_masters,
    int                     **out_offsets,
    int                      *out_total);

/**
 * atmturb_wfs_cuda_build_frame_sublayers - Populate sublayer parameters for a single frame
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @nsuper: Number of superlayers.
 * @supers: Array of superlayer definitions.
 * @screen_offsets: Screen index offsets per layer.
 * @t: Frame index.
 * @time_step_s: Time step in seconds.
 * @offset_os: Offset in oversampled master pixels for guard band alignment.
 * @max_sublayers: Maximum allocated sublayers per superlayer.
 * @sublayers: Output array of sublayer interpolation parameters.
 *
 * Computes wind-displaced coordinates, chromatic secondary shifts, and rolling screen weights.
 */
void atmturb_wfs_cuda_build_frame_sublayers(
    const atmturb_rolling_t       *r,
    const atmturb_geom_t          *geom,
    int                            nsuper,
    const atmturb_superlayer_t    *supers,
    const int                     *screen_offsets,
    long                           t,
    double                         time_step_s,
    double                         offset_os,
    int                            max_sublayers,
    atmturb_cuda_rytov_sublayer_t *sublayers);

/**
 * atmturb_wfs_cuda_build_frame_lf_modes - Populate subharmonic low-frequency modes for a frame
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @nsuper: Number of superlayers.
 * @supers: Array of superlayer definitions.
 * @t: Frame index.
 * @time_step_s: Time step in seconds.
 * @offset_os: Offset in oversampled master pixels for guard band alignment.
 * @max_lf_modes: Maximum allocated modes per superlayer.
 * @lf_modes: Output array of subharmonic Fourier mode parameters.
 *
 * Precomputes spatial wavevectors, primary/secondary initial phases, and temporal amplitudes.
 */
void atmturb_wfs_cuda_build_frame_lf_modes(
    const atmturb_rolling_t     *r,
    const atmturb_geom_t        *geom,
    int                          nsuper,
    const atmturb_superlayer_t  *supers,
    long                         t,
    double                       time_step_s,
    double                       offset_os,
    int                          max_lf_modes,
    atmturb_cuda_lowfreq_mode_t *lf_modes);

#ifdef __cplusplus
}
#endif

#endif /* ATMTURB_WFS_CUDA_COMMON_H */
