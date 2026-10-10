// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_superlayer_nodes.c
 * @brief   Distance node generation and linear interpolation for diffractive propagation
 */

#include <math.h>
#include <stdlib.h>

#include "atmturb_geometry.h"
#include "atmturb_profile.h"
#include "atmturb_superlayer.h"
#include "atmturb_superlayer_internal.h"

/**
 * atmturb_superlayer_count_nodes - Count total distance nodes across all layer clusters
 * @order: Sorted layer index array.
 * @geom: Computed observing geometry.
 * @nlayers: Total count of active turbulence layers.
 * @z_bin_m: Altitude binning distance threshold [meters].
 *
 * Return: Total number of distance nodes.
 */
int atmturb_superlayer_count_nodes(
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

    int total_nodes = 0;
    int idx_start   = 0;

    for (int i = 0; i < nlayers; i++)
    {
        int is_last = (i == nlayers - 1);
        int split   = 0;

        if (!is_last)
        {
            double d_curr = geom->layers[order[i]].dist_m;
            double d_next = geom->layers[order[i + 1]].dist_m;
            if ((d_curr - d_next) > z_bin_m)
            {
                split = 1;
            }
        }

        if (is_last || split)
        {
            double z_top = geom->layers[order[idx_start]].dist_m;
            double z_bot = geom->layers[order[i]].dist_m;
            double L     = z_top - z_bot;
            if (L < 1e-3)
            {
                total_nodes += 1;
            }
            else
            {
                int n_int = (int) ceil(L / z_bin_m);
                if (n_int < 1)
                {
                    n_int = 1;
                }
                total_nodes += (n_int + 1);
            }
            idx_start = i + 1;
        }
    }

    return total_nodes;
}

/**
 * atmturb_superlayer_calc_weight - Compute layer interpolation weight for a candidate node
 * @d: Distance of physical layer [meters].
 * @z_top: Top distance of cluster [meters].
 * @dz: Interval step size [meters].
 * @n_int: Number of intervals in cluster.
 * @node_idx: Candidate node index within cluster.
 *
 * Return: Fractional weight in [0, 1].
 */
static double atmturb_superlayer_calc_weight(
    double d,
    double z_top,
    double dz,
    int    n_int,
    int    node_idx)
{
    int m = (int) floor((z_top - d) / dz);
    if (m < 0)
    {
        m = 0;
    }
    if (m >= n_int)
    {
        m = n_int - 1;
    }

    double z_lower = z_top - (double) (m + 1) * dz;
    double alpha   = (d - z_lower) / dz;
    if (alpha < 0.0)
    {
        alpha = 0.0;
    }
    if (alpha > 1.0)
    {
        alpha = 1.0;
    }

    if (m == node_idx)
    {
        return alpha;
    }
    if (m + 1 == node_idx)
    {
        return 1.0 - alpha;
    }
    return 0.0;
}

/**
 * atmturb_superlayer_populate_single_node - Populate a single node for an isolated cluster
 * @sl: Superlayer node structure to populate.
 * @order: Sorted layer index array.
 * @idx_start: Starting layer index in order.
 * @count: Number of layers in cluster.
 * @prof: Active turbulence profile.
 * @geom: Computed observing geometry.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
static int atmturb_superlayer_populate_single_node(
    atmturb_superlayer_t    *sl,
    const int               *order,
    int                      idx_start,
    int                      count,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom)
{
    sl->nlayers       = count;
    sl->dist_m        = geom->layers[order[idx_start]].dist_m;
    sl->dist_s_m      = (geom->layers[order[idx_start]].path_s_m > 0.0)
                            ? geom->layers[order[idx_start]].path_s_m
                            : sl->dist_m;
    sl->layer_indices = (int *) malloc(sizeof(int) * (size_t) count);
    sl->layer_weights = (float *) malloc(sizeof(float) * (size_t) count);
    if (sl->layer_indices == NULL || sl->layer_weights == NULL)
    {
        return -1;
    }

    double sum_w = 0.0;
    for (int j = 0; j < count; j++)
    {
        int k = order[idx_start + j];
        sl->layer_indices[j] = k;
        sl->layer_weights[j] = 1.0f;
        sum_w += prof->layers[k].cn2_frac;
    }

    atmturb_superlayer_calc_chromatic(sl, prof, geom, sum_w);
    return 0;
}

/**
 * atmturb_superlayer_init_node - Populate one node in an interpolated cluster
 * @sl: Superlayer node structure to populate.
 * @node_dist: Physical distance assigned to this node [meters].
 * @order: Sorted layer index array.
 * @idx_start: Starting layer index in order.
 * @count: Number of layers in cluster.
 * @prof: Active turbulence profile.
 * @geom: Computed observing geometry.
 * @z_top: Top distance of cluster [meters].
 * @dz: Step size [meters].
 * @n_int: Number of intervals.
 * @node_idx: Index of node within cluster.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
static int atmturb_superlayer_init_node(
    atmturb_superlayer_t    *sl,
    double                   node_dist,
    const int               *order,
    int                      idx_start,
    int                      count,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    double                   z_top,
    double                   dz,
    int                      n_int,
    int                      node_idx)
{
    sl->dist_m   = node_dist;
    sl->dist_s_m = node_dist;

    int n_contrib = 0;
    for (int l = 0; l < count; l++)
    {
        double d = geom->layers[order[idx_start + l]].dist_m;
        double w = atmturb_superlayer_calc_weight(d, z_top, dz, n_int, node_idx);
        if (w > 1e-6)
        {
            n_contrib++;
        }
    }

    if (n_contrib == 0)
    {
        sl->nlayers       = 0;
        sl->layer_indices = NULL;
        sl->layer_weights = NULL;
        return 0;
    }

    sl->nlayers       = n_contrib;
    sl->layer_indices = (int *) malloc(sizeof(int) * (size_t) n_contrib);
    sl->layer_weights = (float *) malloc(sizeof(float) * (size_t) n_contrib);
    if (sl->layer_indices == NULL || sl->layer_weights == NULL)
    {
        return -1;
    }

    int idx = 0;
    double sum_w = 0.0;
    for (int l = 0; l < count; l++)
    {
        int k = order[idx_start + l];
        double d = geom->layers[k].dist_m;
        double w = atmturb_superlayer_calc_weight(d, z_top, dz, n_int, node_idx);
        if (w > 1e-6)
        {
            if (w > 1.0)
            {
                w = 1.0;
            }
            sl->layer_indices[idx] = k;
            sl->layer_weights[idx] = (float) w;
            sum_w += w * prof->layers[k].cn2_frac;
            idx++;
        }
    }

    atmturb_superlayer_calc_chromatic(sl, prof, geom, sum_w);
    return 0;
}

/**
 * atmturb_superlayer_build_cluster_nodes - Populate interpolated nodes for one cluster
 * @supers: Destination superlayer array.
 * @s_curr: Current index in supers array.
 * @order: Sorted layer index array.
 * @idx_start: Starting layer index in order.
 * @count: Number of layers in cluster.
 * @prof: Active turbulence profile.
 * @geom: Computed observing geometry.
 * @z_bin_m: Altitude binning distance threshold [meters].
 *
 * Return: Number of nodes created, or -1 on allocation failure.
 */
static int atmturb_superlayer_build_cluster_nodes(
    atmturb_superlayer_t    *supers,
    int                      s_curr,
    const int               *order,
    int                      idx_start,
    int                      count,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    double                   z_bin_m)
{
    double z_top = geom->layers[order[idx_start]].dist_m;
    double z_bot = geom->layers[order[idx_start + count - 1]].dist_m;
    double L     = z_top - z_bot;

    if (L < 1e-3)
    {
        if (atmturb_superlayer_populate_single_node(&supers[s_curr], order, idx_start,
                                                    count, prof, geom) != 0)
        {
            return -1;
        }
        return 1;
    }

    int n_int = (int) ceil(L / z_bin_m);
    if (n_int < 1)
    {
        n_int = 1;
    }
    double dz      = L / (double) n_int;
    int    n_nodes = n_int + 1;

    for (int j = 0; j < n_nodes; j++)
    {
        double z_node = (j == n_nodes - 1) ? z_bot : (z_top - (double) j * dz);
        if (atmturb_superlayer_init_node(&supers[s_curr + j], z_node, order,
                                         idx_start, count, prof, geom,
                                         z_top, dz, n_int, j) != 0)
        {
            return -1;
        }
    }
    return n_nodes;
}

/**
 * atmturb_superlayer_populate_nodes - Build all distance nodes across clusters
 * @supers: Array of superlayer node structures to populate.
 * @order: Sorted layer index array.
 * @prof: Active turbulence profile.
 * @geom: Computed observing geometry.
 * @nlayers: Total count of active turbulence layers.
 * @z_bin_m: Altitude binning distance threshold [meters].
 *
 * Return: Total number of successfully built nodes, or -1 on allocation failure.
 */
int atmturb_superlayer_populate_nodes(
    atmturb_superlayer_t    *supers,
    const int               *order,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    int                      nlayers,
    double                   z_bin_m)
{
    int s         = 0;
    int idx_start = 0;

    for (int i = 0; i < nlayers; i++)
    {
        int is_last = (i == nlayers - 1);
        int split   = 0;

        if (!is_last && z_bin_m > 0.0)
        {
            double d_curr = geom->layers[order[i]].dist_m;
            double d_next = geom->layers[order[i + 1]].dist_m;
            if ((d_curr - d_next) > z_bin_m)
            {
                split = 1;
            }
        }

        if (is_last || split)
        {
            int count = i - idx_start + 1;
            int n_built = atmturb_superlayer_build_cluster_nodes(supers, s, order, idx_start,
                                                                count, prof, geom, z_bin_m);
            if (n_built < 0)
            {
                return -1;
            }
            s += n_built;
            idx_start = i + 1;
        }
    }

    for (int m = 0; m < s - 1; m++)
    {
        double step = supers[m].dist_m - supers[m + 1].dist_m;
        supers[m].step_dist_m = (step > 0.0) ? step : 0.0;
    }
    if (s > 0)
    {
        supers[s - 1].step_dist_m = (supers[s - 1].dist_m > 0.0) ? supers[s - 1].dist_m : 0.0;
    }

    return s;
}
