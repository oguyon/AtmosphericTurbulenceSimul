// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_wfs_stream.h
 * @brief   Header for frame-by-frame 2D shared memory wavefront streaming
 */

#ifndef ATMTURB_WFS_STREAM_H
#define ATMTURB_WFS_STREAM_H

#include "AtmosphericTurbulence.h"
#include "atmturb_geometry.h"
#include "atmturb_profile.h"
#include "atmturb_rolling.h"
#include "atmturb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * make_AtmosphericTurbulence_wavefront_stream - Stream 2D wavefront frames to SHM
 * @slambdaum: Secondary observing wavelength in um.
 * @WFprecision: Precision mode flag (0=single, 1=double).
 * @stream_mode: Streaming mode (1=continuous stream, 2=finite stream, 3=unpaced continuous).
 *
 * Return: 0 on success, -1 on failure.
 */
int make_AtmosphericTurbulence_wavefront_stream(
    float slambdaum,
    long  WFprecision,
    int   stream_mode);

/**
 * atmturb_wfs_validate_config - Validate configuration before wavefront simulation
 *
 * Return: 0 if valid, -1 on configuration error.
 */
int atmturb_wfs_validate_config(void);

/**
 * atmturb_wfs_init_obs_params - Populate observation parameters from global configuration
 * @slambdaum: Secondary observing wavelength in um.
 * @params: Observation parameters container to populate.
 */
void atmturb_wfs_init_obs_params(
    float                 slambdaum,
    atmturb_obs_params_t *params);

/**
 * atmturb_wfs_setup_sim - Initialize atmospheric simulation structures
 * @slambdaum: Secondary observing wavelength in um.
 * @WFprecision: Precision mode flag.
 * @nbframes: Number of frames (or estimated pool frames).
 * @prof: Destination profile struct.
 * @params: Destination observation params struct.
 * @geom: Destination geometry struct.
 * @rsim: Destination rolling screens struct.
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_wfs_setup_sim(
    float                 slambdaum,
    long                  WFprecision,
    long                  nbframes,
    atmturb_profile_t    *prof,
    atmturb_obs_params_t *params,
    atmturb_geom_t       *geom,
    atmturb_rolling_t    *rsim);

/**
 * atmturb_wfs_teardown_sim - Free atmospheric simulation structures
 * @prof: Active profile struct.
 * @geom: Active geometry struct.
 * @rsim: Active rolling screens struct.
 */
void atmturb_wfs_teardown_sim(
    atmturb_profile_t *prof,
    atmturb_geom_t    *geom,
    atmturb_rolling_t *rsim);

#ifdef __cplusplus
}
#endif

#endif // ATMTURB_WFS_STREAM_H
