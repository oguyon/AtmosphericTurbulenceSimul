// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_superlayer_internal.h
 * @brief   Internal declarations for super-layer binning and node interpolation
 */

#ifndef ATMTURB_SUPERLAYER_INTERNAL_H
#define ATMTURB_SUPERLAYER_INTERNAL_H

#include "atmturb_geometry.h"
#include "atmturb_profile.h"
#include "atmturb_superlayer.h"

/**
 * atmturb_superlayer_sort_layers - Sort layer indices descending by LOS distance
 * @order:   Output array of @nlayers layer indices sorted top-to-bottom.
 * @geom:    Computed observing geometry.
 * @nlayers: Total number of active turbulence layers.
 */
void atmturb_superlayer_sort_layers(
    int                  *order,
    const atmturb_geom_t *geom,
    int                   nlayers);

/**
 * atmturb_superlayer_calc_chromatic - Compute weighted chromatic shifts and spreads
 * @sl:    Super-layer record to populate with chromatic statistics.
 * @prof:  Active turbulence profile.
 * @geom:  Computed observing geometry.
 * @sum_w: Sum of effective Cn2 weights in @sl.
 */
void atmturb_superlayer_calc_chromatic(
    atmturb_superlayer_t    *sl,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    double                   sum_w);

/**
 * atmturb_superlayer_count_nodes - Count linear z-interpolation nodes across clusters
 * @order:   Layer index array sorted descending by distance.
 * @geom:    Computed observing geometry.
 * @nlayers: Total number of active turbulence layers.
 * @z_bin_m: Distance interval threshold [m].
 *
 * Return: Total number of interpolated super-layer nodes.
 */
int atmturb_superlayer_count_nodes(
    const int            *order,
    const atmturb_geom_t *geom,
    int                   nlayers,
    double                z_bin_m);

/**
 * atmturb_superlayer_populate_nodes - Build linear z-interpolated super-layer nodes
 * @supers:  Pre-allocated super-layer node array to populate.
 * @order:   Layer index array sorted descending by distance.
 * @prof:    Active turbulence profile.
 * @geom:    Computed observing geometry.
 * @nlayers: Total number of active turbulence layers.
 * @z_bin_m: Distance interval threshold [m].
 *
 * Return: Number of populated nodes on success, or -1 on allocation failure.
 */
int atmturb_superlayer_populate_nodes(
    atmturb_superlayer_t    *supers,
    const int               *order,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    int                      nlayers,
    double                   z_bin_m);

#endif // ATMTURB_SUPERLAYER_INTERNAL_H
