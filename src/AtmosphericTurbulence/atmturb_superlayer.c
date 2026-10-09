// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_superlayer.c
 * @brief   Super-layer grouping and chromatic binning implementation
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "atmturb_geometry.h"
#include "atmturb_profile.h"
#include "atmturb_superlayer.h"

/**
 * atmturb_superlayer_sort_layers - Sort layer indices descending by line-of-sight distance
 * @order: Output array of sorted layer indices.
 * @geom: Computed observing geometry.
 * @nlayers: Total count of active turbulence layers.
 */
static void atmturb_superlayer_sort_layers(
    int                  *order,
    const atmturb_geom_t *geom,
    int                   nlayers)
{
    for (int i = 0; i < nlayers; i++)
    {
        order[i] = i;
    }
    for (int i = 0; i < nlayers - 1; i++)
    {
        for (int j = i + 1; j < nlayers; j++)
        {
            if (geom->layers[order[j]].dist_m > geom->layers[order[i]].dist_m)
            {
                int tmp  = order[i];
                order[i] = order[j];
                order[j] = tmp;
            }
        }
    }
}

/**
 * atmturb_superlayer_count_bins - Determine number of super-layers based on distance threshold
 * @order: Sorted layer index array.
 * @geom: Computed observing geometry.
 * @nlayers: Total count of active turbulence layers.
 * @z_bin_m: Altitude binning distance threshold [meters].
 *
 * Return: Number of super-layers.
 */
static int atmturb_superlayer_count_bins(
    const int            *order,
    const atmturb_geom_t *geom,
    int                   nlayers,
    double                z_bin_m)
{
    if (nlayers <= 0)
    {
        return 0;
    }
    if (z_bin_m <= 0.0)
    {
        return nlayers;
    }

    int    nsuper      = 1;
    double anchor_dist = geom->layers[order[0]].dist_m;

    for (int i = 1; i < nlayers; i++)
    {
        double d = geom->layers[order[i]].dist_m;
        if ((anchor_dist - d) > z_bin_m)
        {
            nsuper++;
            anchor_dist = d;
        }
    }

    return nsuper;
}

/**
 * atmturb_superlayer_calc_chromatic - Derive chromatic offsets and spreads for a bin
 * @sl: Pointer to superlayer structure being populated.
 * @prof: Active turbulence profile.
 * @geom: Computed observing geometry.
 * @sum_w: Sum of cn2 weights for this superlayer.
 */
static void atmturb_superlayer_calc_chromatic(
    atmturb_superlayer_t    *sl,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    double                   sum_w)
{
    double sum_dx = 0.0;
    double sum_dy = 0.0;
    double sum_wr = 0.0;
    double os     = (geom->oversample > 0) ? (double) geom->oversample : 1.0;

    for (int j = 0; j < sl->nlayers; j++)
    {
        int k = sl->layer_indices[j];
        double w = prof->layers[k].cn2_frac;
        double dx_pup = geom->layers[k].d_chrom_x / os;
        double dy_pup = geom->layers[k].d_chrom_y / os;
        double wr = (geom->layers[k].weight > 0.0)
                    ? (geom->layers[k].weight_s / geom->layers[k].weight)
                    : 1.0;

        sum_dx += w * dx_pup;
        sum_dy += w * dy_pup;
        sum_wr += w * wr;
    }

    sl->chrom_dx_px  = (sum_w > 0.0) ? (sum_dx / sum_w) : 0.0;
    sl->chrom_dy_px  = (sum_w > 0.0) ? (sum_dy / sum_w) : 0.0;
    sl->weight_ratio = (sum_w > 0.0) ? (sum_wr / sum_w) : 1.0;

    double max_dev_px = 0.0;
    double max_dev_wr = 0.0;

    for (int j = 0; j < sl->nlayers; j++)
    {
        int k = sl->layer_indices[j];
        double dx_pup = geom->layers[k].d_chrom_x / os;
        double dy_pup = geom->layers[k].d_chrom_y / os;
        double wr = (geom->layers[k].weight > 0.0)
                    ? (geom->layers[k].weight_s / geom->layers[k].weight)
                    : 1.0;

        double ddx = dx_pup - sl->chrom_dx_px;
        double ddy = dy_pup - sl->chrom_dy_px;
        double dev_px = sqrt(ddx * ddx + ddy * ddy);
        if (dev_px > max_dev_px)
        {
            max_dev_px = dev_px;
        }

        double dev_wr = fabs(wr - sl->weight_ratio);
        if (sl->weight_ratio > 0.0)
        {
            dev_wr /= sl->weight_ratio;
        }
        if (dev_wr > max_dev_wr)
        {
            max_dev_wr = dev_wr;
        }
    }

    sl->chrom_spread_px     = max_dev_px;
    sl->weight_ratio_spread = max_dev_wr;
}

/**
 * atmturb_superlayer_populate_bins - Group sorted layers into bins and compute centroids
 * @supers: Array of super-layers to populate.
 * @order: Sorted layer index array.
 * @prof: Active turbulence profile.
 * @geom: Computed observing geometry.
 * @nlayers: Total count of active turbulence layers.
 * @z_bin_m: Altitude binning distance threshold [meters].
 *
 * Return: Number of successfully built super-layers, or -1 on allocation failure.
 */
static int atmturb_superlayer_populate_bins(
    atmturb_superlayer_t    *supers,
    const int               *order,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    int                      nlayers,
    double                   z_bin_m)
{
    int    s           = 0;
    int    idx_start   = 0;
    double anchor_dist = geom->layers[order[0]].dist_m;

    for (int i = 0; i < nlayers; i++)
    {
        int is_last = (i == nlayers - 1);
        int split   = 0;

        if (!is_last && z_bin_m > 0.0)
        {
            double next_dist = geom->layers[order[i + 1]].dist_m;
            if ((anchor_dist - next_dist) > z_bin_m)
            {
                split = 1;
            }
        }

        if (is_last || split)
        {
            int count = i - idx_start + 1;
            supers[s].nlayers = count;
            supers[s].layer_indices = (int *) malloc(sizeof(int) * (size_t) count);
            if (supers[s].layer_indices == NULL)
            {
                return -1;
            }

            double sum_w  = 0.0;
            double sum_wd = 0.0;
            for (int j = 0; j < count; j++)
            {
                int k = order[idx_start + j];
                supers[s].layer_indices[j] = k;
                double w = prof->layers[k].cn2_frac;
                double d = geom->layers[k].dist_m;
                sum_w  += w;
                sum_wd += w * d;
            }

            supers[s].dist_m = (sum_w > 0.0) ? (sum_wd / sum_w)
                                             : geom->layers[order[idx_start]].dist_m;

            atmturb_superlayer_calc_chromatic(&supers[s], prof, geom, sum_w);

            s++;
            idx_start = i + 1;
            if (!is_last)
            {
                anchor_dist = geom->layers[order[idx_start]].dist_m;
            }
        }
    }

    for (int m = 0; m < s - 1; m++)
    {
        double step = supers[m].dist_m - supers[m + 1].dist_m;
        supers[m].step_dist_m = (step > 0.0) ? step : 0.0;
    }
    supers[s - 1].step_dist_m = (supers[s - 1].dist_m > 0.0) ? supers[s - 1].dist_m : 0.0;

    return s;
}

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
    double                   pixscale_m)
{
    (void) pixscale_m;
    if (supers_out == NULL || nsuper_out == NULL || prof == NULL ||
        geom == NULL || prof->nlayers <= 0)
    {
        return -1;
    }

    int nlayers = prof->nlayers;
    int *order  = (int *) malloc(sizeof(int) * (size_t) nlayers);
    if (order == NULL)
    {
        return -1;
    }

    atmturb_superlayer_sort_layers(order, geom, nlayers);
    int nsuper = atmturb_superlayer_count_bins(order, geom, nlayers, z_bin_m);

    atmturb_superlayer_t *supers = (atmturb_superlayer_t *) calloc((size_t) nsuper,
                                                                   sizeof(atmturb_superlayer_t));
    if (supers == NULL)
    {
        free(order);
        return -1;
    }

    int s_built = atmturb_superlayer_populate_bins(supers, order, prof, geom, nlayers, z_bin_m);
    free(order);
    if (s_built < 0)
    {
        atmturb_superlayer_free(supers, nsuper);
        return -1;
    }

    *supers_out = supers;
    *nsuper_out = s_built;
    return 0;
}

/**
 * atmturb_superlayer_free - Release memory allocated for super-layer array
 * @supers: Array of super-layers to release.
 * @nsuper: Number of super-layers in array.
 */
void atmturb_superlayer_free(
    atmturb_superlayer_t *supers,
    int                   nsuper)
{
    if (supers == NULL)
    {
        return;
    }

    for (int m = 0; m < nsuper; m++)
    {
        free(supers[m].layer_indices);
        supers[m].layer_indices = NULL;
    }
    free(supers);
}
