// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_profile.h
 * @brief   Atmospheric turbulence layer profile data types and loader interface
 */

#ifndef ATMTURB_PROFILE_H
#define ATMTURB_PROFILE_H

#include <stdint.h>

/**
 * struct atmturb_layer_t - Parameters of one atmospheric turbulence layer
 * @alt_m: Layer altitude above sea level in meters.
 * @cn2_frac: Normalized Cn2 fraction (sum over all layers = 1.0).
 * @speed_mps: Wind velocity magnitude in meters per second.
 * @dir_rad: Wind velocity direction in radians.
 * @L0_m: Turbulence outer scale in meters (<= 0 for infinite).
 * @l0_m: Turbulence inner scale in meters (<= 0 for disabled).
 * @sigma_wind_mps: Wind speed turbulence standard deviation [m/s].
 * @L_wind_m: Wind turbulence outer scale [m].
 */
typedef struct
{
    double alt_m;
    double cn2_frac;
    double speed_mps;
    double dir_rad;
    double L0_m;
    double l0_m;
    double sigma_wind_mps;
    double L_wind_m;
} atmturb_layer_t;

/**
 * struct atmturb_profile_t - Multi-layer atmospheric profile container
 * @nlayers: Total count of active turbulence layers.
 * @total_cn2_raw: Raw integral sum of Cn2 values before normalization.
 * @layers: Dynamic array of turbulence layer descriptions.
 */
typedef struct
{
    int             nlayers;
    double          total_cn2_raw;
    atmturb_layer_t *layers;
} atmturb_profile_t;

/**
 * atmturb_profile_load - Load and normalize turbulence profile from disk
 * @fname: Path to text profile file (or NULL to load built-in default).
 * @prof: Output profile container to initialize.
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_profile_load(
    const char        *fname,
    atmturb_profile_t *prof);

/**
 * atmturb_profile_free - Deallocate turbulence layer array
 * @prof: Pointer to profile container.
 */
void atmturb_profile_free(
    atmturb_profile_t *prof);

/**
 * atmturb_profile_init_default - Populate standard 7-layer atmospheric profile
 * @prof: Output profile container to initialize.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_profile_init_default(
    atmturb_profile_t *prof);

#endif // ATMTURB_PROFILE_H
