// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_geometry.h
 * @brief   Observing geometry, differential refraction, and layer parameters
 */

#ifndef ATMTURB_GEOMETRY_H
#define ATMTURB_GEOMETRY_H

#include <stdint.h>
#include "atmturb_profile.h"

/**
 * struct atmturb_layer_geom_t - Per-layer geometric parameters for extrusion
 * @weight: Screen amplitude multiplier at reference wavelength.
 * @weight_s: Screen amplitude multiplier at secondary wavelength.
 * @x0: Initial primary x coordinate in master pixels.
 * @y0: Initial primary y coordinate in master pixels.
 * @xs0: Initial secondary x coordinate in master pixels (includes chromatic shift).
 * @ys0: Initial secondary y coordinate in master pixels.
 * @vx_pix: Pupil-plane velocity in x [master pixels / frame].
 * @vy_pix: Pupil-plane velocity in y [master pixels / frame].
 * @dist_m: Line-of-sight slant distance from telescope pupil [meters].
 * @t_dec_s: Rolling screen epoch duration [seconds].
 */
typedef struct
{
    double weight;
    double weight_s;
    double x0;
    double y0;
    double xs0;
    double ys0;
    double vx_pix;
    double vy_pix;
    double dist_m;
    double t_dec_s;
} atmturb_layer_geom_t;

/**
 * struct atmturb_geom_t - Full observing geometry context
 * @r0_ref_m: Kolmogorov Fried parameter at lambda_ref at zenith [m].
 * @r0_ref_pix: Reference r0 in master grid pixel units.
 * @dx_master_m: Physical size of master pixel [m/pixel].
 * @site_alt_m: Telescope altitude above sea level [m].
 * @cos_z: Cosine of zenith angle.
 * @sin_z: Sine of zenith angle.
 * @ez_x: Unit vector toward zenith in pupil frame (x component).
 * @ez_y: Unit vector toward zenith in pupil frame (y component).
 * @nlayers: Number of layers.
 * @layers: Array of per-layer geometric properties.
 */
typedef struct
{
    double               r0_ref_m;
    double               r0_ref_pix;
    double               dx_master_m;
    double               site_alt_m;
    double               cos_z;
    double               sin_z;
    double               ez_x;
    double               ez_y;
    int                  nlayers;
    atmturb_layer_geom_t *layers;
} atmturb_geom_t;

/**
 * struct atmturb_obs_params_t - Input parameters for geometry calculation
 * @lambda_ref_m: Reference wavelength [m] (e.g. 0.5e-6).
 * @lambda_s_m: Secondary observing wavelength [m] (e.g. 1.65e-6).
 * @seeing_arcsec: Seeing at zenith at reference wavelength [arcseconds].
 * @zenith_rad: Zenith angle [radians].
 * @parallactic_rad: Parallactic angle [radians].
 * @site_alt_m: Telescope altitude ASL [m] (-1 for auto-detect).
 * @pupil_scale_m: Pupil sampling scale [m/pixel].
 * @oversample: Master grid oversampling factor (1 or 2).
 * @master_size: Linear dimension of master screen [pixels].
 * @time_step_s: Time step between frames [seconds].
 * @source_x_rad: Off-axis source x angular position [radians].
 * @source_y_rad: Off-axis source y angular position [radians].
 * @seed: Master RNG seed.
 */
typedef struct
{
    double   lambda_ref_m;
    double   lambda_s_m;
    double   seeing_arcsec;
    double   zenith_rad;
    double   parallactic_rad;
    double   site_alt_m;
    double   pupil_scale_m;
    int      oversample;
    long     master_size;
    double   time_step_s;
    double   source_x_rad;
    double   source_y_rad;
    uint64_t seed;
} atmturb_obs_params_t;

/**
 * atmturb_geometry_compute - Derive all layer geometric parameters and scales
 * @prof: Active turbulence profile.
 * @params: Observing conditions and simulation parameters.
 * @geom: Output geometry structure to populate.
 *
 * Return: 0 on success, -1 on invalid configuration or allocation failure.
 */
int atmturb_geometry_compute(
    const atmturb_profile_t    *prof,
    const atmturb_obs_params_t *params,
    atmturb_geom_t             *geom);

/**
 * atmturb_geometry_free - Release geometry layer array
 * @geom: Pointer to geometry structure.
 */
void atmturb_geometry_free(
    atmturb_geom_t *geom);

/**
 * atmturb_refraction_ray_shift - Compute cumulative atmospheric refraction displacement
 * @h_layer: Layer altitude above sea level [m].
 * @h_site: Telescope site altitude above sea level [m].
 * @zenith_angle: Apparent zenith angle [rad].
 * @lambda_m: Optical wavelength [m].
 *
 * Return: Lateral deflection in meters relative to an unrefracted straight ray.
 */
double atmturb_refraction_ray_shift(
    double h_layer,
    double h_site,
    double zenith_angle,
    double lambda_m);

#endif // ATMTURB_GEOMETRY_H
