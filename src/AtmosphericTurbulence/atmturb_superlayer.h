// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_superlayer.h
 * @brief   Super-layer grouping and chromatic binning for diffractive propagation
 */

#ifndef ATMTURB_SUPERLAYER_H
#define ATMTURB_SUPERLAYER_H

#include "atmturb_geometry.h"
#include "atmturb_profile.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * struct atmturb_superlayer_t - Grouped atmospheric layers for Fresnel diffraction
 * @nlayers: Count of profile layers binned into this super-layer.
 * @layer_indices: Array of indices in geom->layers / rsim->layers.
 * @layer_weights: Fractional layer weights (NULL for centroid binning, 1.0 default).
 * @dist_m: Line-of-sight centroid distance from telescope pupil [meters].
 * @dist_s_m: Secondary wavelength centroid distance from telescope pupil [meters].
 * @step_dist_m: Propagation distance to the next super-layer or ground [meters].
 * @chrom_dx_px: Weighted chromatic offset along pupil X [pixels].
 * @chrom_dy_px: Weighted chromatic offset along pupil Y [pixels].
 * @chrom_spread_px: Maximum chromatic offset deviation within bin [pixels].
 * @weight_ratio: Weighted ratio of secondary to primary screen weights.
 * @weight_ratio_spread: Maximum relative deviation of weight ratio within bin.
 */
typedef struct
{
    int     nlayers;
    int    *layer_indices;
    float  *layer_weights;
    double  dist_m;
    double  dist_s_m;
    double  step_dist_m;
    double  chrom_dx_px;
    double  chrom_dy_px;
    double  chrom_spread_px;
    double  weight_ratio;
    double  weight_ratio_spread;
} atmturb_superlayer_t;

/**
 * atmturb_superlayer_build - Bin profile layers into super-layers with chromatic statistics
 * @supers_out: Output pointer to allocated super-layer array.
 * @nsuper_out: Output pointer to number of super-layers.
 * @prof: Active turbulence profile.
 * @geom: Computed observing geometry.
 * @z_bin_m: Altitude binning distance threshold [meters].
 * @pixscale_m: Grid physical pixel size [meters/pixel].
 *
 * Return: 0 on success, -1 on allocation failure or invalid input.
 */
int atmturb_superlayer_build(
    atmturb_superlayer_t   **supers_out,
    int                     *nsuper_out,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    double                   z_bin_m,
    double                   pixscale_m);

/**
 * atmturb_superlayer_free - Release memory allocated for super-layer array
 * @supers: Array of super-layers to release.
 * @nsuper: Number of super-layers in array.
 */
void atmturb_superlayer_free(
    atmturb_superlayer_t *supers,
    int                   nsuper);

#ifdef __cplusplus
}
#endif

#endif // ATMTURB_SUPERLAYER_H
