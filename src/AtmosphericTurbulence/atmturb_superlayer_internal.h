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

void atmturb_superlayer_sort_layers(
    int                  *order,
    const atmturb_geom_t *geom,
    int                   nlayers);

void atmturb_superlayer_calc_chromatic(
    atmturb_superlayer_t    *sl,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    double                   sum_w);

int atmturb_superlayer_count_nodes(
    const int            *order,
    const atmturb_geom_t *geom,
    int                   nlayers,
    double                z_bin_m);

int atmturb_superlayer_populate_nodes(
    atmturb_superlayer_t    *supers,
    const int               *order,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    int                      nlayers,
    double                   z_bin_m);

#endif // ATMTURB_SUPERLAYER_INTERNAL_H
