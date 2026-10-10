// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_geometry.c
 * @brief   Observing geometry, differential atmospheric refraction, and layer parameters
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "AtmosphereModel/AtmosphereModel.h"
#include "AtmosphericTurbulence.h"
#include "atmturb_geometry.h"
#include "atmturb_types.h"

/**
 * atmturb_geometry_free - Release geometry layer array
 * @geom: Pointer to geometry structure.
 */
void atmturb_geometry_free(
    atmturb_geom_t *geom)
{
    if (geom == NULL)
    {
        return;
    }
    if (geom->layers != NULL)
    {
        for (int k = 0; k < geom->nlayers; k++)
        {
            free(geom->layers[k].traj_x);
            free(geom->layers[k].traj_y);
            geom->layers[k].traj_x = NULL;
            geom->layers[k].traj_y = NULL;
        }
    }
    free(geom->layers);
    geom->layers = NULL;
    geom->nlayers = 0;
}

/**
 * atmturb_refraction_ray_shift_path - Compute refraction shift and curved ray path length
 * @h_layer: Layer altitude above sea level [m].
 * @h_site: Telescope site altitude above sea level [m].
 * @zenith_angle: Apparent zenith angle [rad].
 * @lambda_m: Optical wavelength [m].
 * @path_len_out: Output pointer for curved ray path length in meters (or NULL).
 *
 * Return: Lateral deflection in meters relative to an unrefracted straight ray.
 */
double atmturb_refraction_ray_shift_path(
    double  h_layer,
    double  h_site,
    double  zenith_angle,
    double  lambda_m,
    double *path_len_out)
{
    double cos_z0 = cos(zenith_angle);
    double straight_dist = (zenith_angle <= 1e-6)
                               ? (h_layer - h_site)
                               : ((cos_z0 > 1e-9) ? (h_layer - h_site) / cos_z0 : 0.0);
    if (straight_dist < 0.0)
    {
        straight_dist = 0.0;
    }

    if (zenith_angle <= 1e-6 || h_layer <= h_site)
    {
        if (path_len_out != NULL)
        {
            *path_len_out = straight_dist;
        }
        return 0.0;
    }

    double sin_z0 = sin(zenith_angle);
    double tan_z0 = tan(zenith_angle);
    double h_curr = h_site;
    double shift = 0.0;
    double path_len = 0.0;
    const double step = 10.0;

    while (h_curr < h_layer)
    {
        double next_h = h_curr + step;
        if (next_h > h_layer)
        {
            next_h = h_layer;
        }
        double dh = next_h - h_curr;
        double h_mid = 0.5 * (h_curr + next_h);
        double N = (double) AtmosphereModel_stdAtmModel_N((float) h_mid, (float) lambda_m, 0);
        double sin_z = sin_z0 / (1.0 + N);
        if (sin_z > 0.9999)
        {
            sin_z = 0.9999;
        }
        double cos_z = sqrt(1.0 - sin_z * sin_z);
        double tan_z = (cos_z > 1e-9) ? (sin_z / cos_z) : tan_z0;
        shift += (tan_z - tan_z0) * dh;
        path_len += (cos_z > 1e-9) ? (dh / cos_z) : (dh / cos_z0);
        h_curr = next_h;
    }

    if (path_len_out != NULL)
    {
        *path_len_out = path_len;
    }
    return shift;
}

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
    double lambda_m)
{
    return atmturb_refraction_ray_shift_path(h_layer, h_site, zenith_angle, lambda_m, NULL);
}

/**
 * atmturb_compute_layer_velocity - Project wind onto pupil plane with foreshortening
 * @speed_mps: Horizontal wind speed in m/s.
 * @dir_rad: Wind heading direction in radians.
 * @ez_x: Unit vector toward zenith (x component).
 * @ez_y: Unit vector toward zenith (y component).
 * @cos_z: Cosine of zenith angle.
 * @dt_s: Time step in seconds.
 * @dx_m: Master grid pixel scale in meters.
 * @vx_pix: Output pointer for x velocity in master pixels per frame.
 * @vy_pix: Output pointer for y velocity in master pixels per frame.
 */
static void atmturb_compute_layer_velocity(
    double  speed_mps,
    double  dir_rad,
    double  ez_x,
    double  ez_y,
    double  cos_z,
    double  dt_s,
    double  dx_m,
    double *vx_pix,
    double *vy_pix)
{
    double v_raw_x = speed_mps * cos(dir_rad);
    double v_raw_y = speed_mps * sin(dir_rad);

    double ep_x = ez_y;
    double ep_y = -ez_x;

    double v_parallel = v_raw_x * ez_x + v_raw_y * ez_y;
    double v_perp     = v_raw_x * ep_x + v_raw_y * ep_y;

    double v_proj_x = (v_parallel * cos_z) * ez_x + v_perp * ep_x;
    double v_proj_y = (v_parallel * cos_z) * ez_y + v_perp * ep_y;

    *vx_pix = v_proj_x * dt_s / dx_m;
    *vy_pix = v_proj_y * dt_s / dx_m;
}

/**
 * atmturb_resolve_site_alt - Determine site altitude ASL from config or profile
 * @configured_alt: Configured site altitude in meters (-1 for auto).
 * @prof: Loaded turbulence profile.
 *
 * Return: Site altitude in meters ASL.
 */
static double atmturb_resolve_site_alt(
    double                   configured_alt,
    const atmturb_profile_t *prof)
{
    if (configured_alt >= 0.0)
    {
        return configured_alt;
    }
    double lowest = prof->layers[0].alt_m;
    for (int k = 1; k < prof->nlayers; k++)
    {
        if (prof->layers[k].alt_m < lowest)
        {
            lowest = prof->layers[k].alt_m;
        }
    }
    return lowest;
}

/**
 * atmturb_wrap_coord - Wrap coordinate modulo screen size into [0, size)
 * @val: Real coordinate in pixel units.
 * @size: Master grid size in pixels.
 *
 * Return: Wrapped positive coordinate in [0, size).
 */
static inline double atmturb_wrap_coord(
    double val,
    long   size)
{
    double m = fmod(val, (double) size);
    if (m < 0.0)
    {
        m += (double) size;
    }
    return m;
}

/**
 * atmturb_init_layer_turbulent_wind - Initialize wind trajectory table if layer has turbulent wind
 * @layer: Layer physical specification.
 * @params: Observing simulation parameters.
 * @geom: Observing geometry context.
 * @k: Layer index.
 * @base_seed: Master PRNG seed.
 * @lg: Layer geometry to populate with trajectory.
 */
static void atmturb_init_layer_turbulent_wind(
    const atmturb_layer_t      *layer,
    const atmturb_obs_params_t *params,
    const atmturb_geom_t       *geom,
    int                         k,
    uint64_t                    base_seed,
    atmturb_layer_geom_t       *lg)
{
    lg->traj_x   = NULL;
    lg->traj_y   = NULL;
    lg->nbframes = 0;

    if (params->nbframes <= 0 || !(layer->sigma_wind_mps > 0.0) || !(layer->L_wind_m > 0.0))
    {
        return;
    }

    lg->traj_x = (double *) malloc(sizeof(double) * (size_t) params->nbframes);
    lg->traj_y = (double *) malloc(sizeof(double) * (size_t) params->nbframes);
    if (lg->traj_x == NULL || lg->traj_y == NULL)
    {
        free(lg->traj_x);
        free(lg->traj_y);
        lg->traj_x = NULL;
        lg->traj_y = NULL;
        return;
    }
    lg->nbframes = params->nbframes;

    atmturb_wind_traj_params_t tp = {
        .nbframes       = params->nbframes,
        .dt_s           = params->time_step_s,
        .vx_pix         = lg->vx_pix,
        .vy_pix         = lg->vy_pix,
        .dx_master_m    = geom->dx_master_m,
        .sigma_wind_mps = layer->sigma_wind_mps,
        .L_wind_m       = layer->L_wind_m,
        .seed           = atmturb_rng_stream_seed(base_seed, (uint64_t) (5000 + k))
    };
    atmturb_wind_synthesize_trajectory(&tp, lg->traj_x, lg->traj_y);
}

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
    atmturb_geom_t             *geom)
{
    if (!prof || prof->nlayers <= 0 || !params || !geom)
    {
        return -1;
    }

    memset(geom, 0, sizeof(*geom));
    geom->site_alt_m = atmturb_resolve_site_alt(params->site_alt_m, prof);
    int os = (params->oversample > 1) ? params->oversample : 1;
    geom->oversample = os;
    geom->interp = (params->interp == 0) ? 0 : 1;
    geom->lowfreq = (params->lowfreq != 0) ? 1 : 0;
    geom->rolling = (params->rolling != 0) ? 1 : 0;
    geom->boil_time_s = (params->boil_time_s > 0.0) ? params->boil_time_s : 0.0;
    geom->dx_master_m = params->pupil_scale_m / (double) os;

    geom->cos_z = cos(params->zenith_rad);
    if (geom->cos_z < 0.05)
    {
        geom->cos_z = 0.05;
    }
    geom->sin_z = sin(params->zenith_rad);

    geom->ez_x = sin(params->parallactic_rad);
    geom->ez_y = cos(params->parallactic_rad);

    double seeing_rad = params->seeing_arcsec * (M_PI / (180.0 * 3600.0));
    geom->r0_ref_m = 0.98 * params->lambda_ref_m / seeing_rad;
    geom->r0_ref_pix = geom->r0_ref_m / geom->dx_master_m;

    geom->layers = (atmturb_layer_geom_t *) calloc((size_t) prof->nlayers,
                                                   sizeof(atmturb_layer_geom_t));
    if (geom->layers == NULL)
    {
        return -1;
    }
    geom->nlayers = prof->nlayers;

    uint64_t base_seed = atmturb_resolve_seed(params->seed);

    for (int k = 0; k < prof->nlayers; k++)
    {
        const atmturb_layer_t *layer = &prof->layers[k];
        atmturb_layer_geom_t *lg = &geom->layers[k];

        double h_layer = (layer->alt_m >= geom->site_alt_m) ? layer->alt_m : geom->site_alt_m;
        double path_ref = 0.0;
        double path_s   = 0.0;
        double r_ref = atmturb_refraction_ray_shift_path(h_layer, geom->site_alt_m,
                                                         params->zenith_rad,
                                                         params->lambda_ref_m,
                                                         &path_ref);
        double r_s   = atmturb_refraction_ray_shift_path(h_layer, geom->site_alt_m,
                                                         params->zenith_rad,
                                                         params->lambda_s_m,
                                                         &path_s);
        lg->path_m   = path_ref;
        lg->path_s_m = path_s;

        double straight_dist = (h_layer - geom->site_alt_m) / geom->cos_z;
        lg->dist_m = (CONF_FRESNEL_REFRACT_PATH == 1 && path_ref > 0.0)
                         ? path_ref
                         : straight_dist;
        lg->weight = sqrt(layer->cn2_frac / geom->cos_z);

        double n_ref = (double) AtmosphereModel_stdAtmModel_N((float) h_layer,
                                                              (float) params->lambda_ref_m, 0);
        double n_s   = (double) AtmosphereModel_stdAtmModel_N((float) h_layer,
                                                              (float) params->lambda_s_m, 0);
        double scoeff = (n_ref != 0.0)
                            ? (params->lambda_ref_m / params->lambda_s_m) * (n_s / n_ref)
                            : (params->lambda_ref_m / params->lambda_s_m);
        lg->weight_s = lg->weight * scoeff;
        double delta_r = r_s - r_ref;
        double d_chrom_x = (delta_r / geom->dx_master_m) * geom->ez_x;
        double d_chrom_y = (delta_r / geom->dx_master_m) * geom->ez_y;
        lg->d_chrom_x = d_chrom_x;
        lg->d_chrom_y = d_chrom_y;

        double d_src_x = (params->source_x_rad * lg->dist_m) / geom->dx_master_m;
        double d_src_y = (params->source_y_rad * lg->dist_m) / geom->dx_master_m;

        uint64_t rng = atmturb_rng_stream_seed(base_seed, (uint64_t) (1000 + k));
        double u1 = (atmturb_rng_splitmix64(&rng) >> 11) * (1.0 / 9007199254740992.0);
        double u2 = (atmturb_rng_splitmix64(&rng) >> 11) * (1.0 / 9007199254740992.0);
        double x0_raw = floor(u1 * (double) params->master_size);
        double y0_raw = floor(u2 * (double) params->master_size);

        lg->x0  = atmturb_wrap_coord(x0_raw + d_src_x, params->master_size);
        lg->y0  = atmturb_wrap_coord(y0_raw + d_src_y, params->master_size);
        lg->xs0 = atmturb_wrap_coord(x0_raw + d_src_x + d_chrom_x, params->master_size);
        lg->ys0 = atmturb_wrap_coord(y0_raw + d_src_y + d_chrom_y, params->master_size);

        atmturb_compute_layer_velocity(layer->speed_mps, layer->dir_rad,
                                       geom->ez_x, geom->ez_y, geom->cos_z,
                                       params->time_step_s, geom->dx_master_m,
                                       &lg->vx_pix, &lg->vy_pix);

        double l_screen_m = (double) params->master_size * geom->dx_master_m;
        double v_eff = sqrt(lg->vx_pix * lg->vx_pix + lg->vy_pix * lg->vy_pix)
                       * (geom->dx_master_m / params->time_step_s);
        double t_dec = (v_eff > 0.1) ? (0.5 * l_screen_m / v_eff) : 10.0;
        if (params->boil_time_s > 0.0 && params->boil_time_s < t_dec)
        {
            t_dec = params->boil_time_s;
        }
        lg->t_dec_s = t_dec;

        atmturb_init_layer_turbulent_wind(layer, params, geom, k, base_seed, lg);
    }

    printf("[milkatmturb] Geometry: r0_ref = %.4f m (%.2f pix), site_alt = %.1f m, cos(z) = %.4f\n",
           geom->r0_ref_m, geom->r0_ref_pix, geom->site_alt_m, geom->cos_z);
    fflush(stdout);

    return 0;
}

/**
 * atmturb_geom_get_layer_dx - Evaluate cumulative x translation for layer at frame t
 * @lg: Layer geometry structure.
 * @t: Simulation frame index.
 *
 * Return: Cumulative x offset in master pixels.
 */
double atmturb_geom_get_layer_dx(
    const atmturb_layer_geom_t *lg,
    long                        t)
{
    if (lg == NULL)
    {
        return 0.0;
    }
    if (lg->traj_x == NULL || lg->nbframes <= 1)
    {
        return (double) t * lg->vx_pix;
    }
    long span  = lg->nbframes - 1;
    long cycle = t / span;
    long rem   = t % span;
    return (double) cycle * lg->traj_x[span] + lg->traj_x[rem];
}

/**
 * atmturb_geom_get_layer_dy - Evaluate cumulative y translation for layer at frame t
 * @lg: Layer geometry structure.
 * @t: Simulation frame index.
 *
 * Return: Cumulative y offset in master pixels.
 */
double atmturb_geom_get_layer_dy(
    const atmturb_layer_geom_t *lg,
    long                        t)
{
    if (lg == NULL)
    {
        return 0.0;
    }
    if (lg->traj_y == NULL || lg->nbframes <= 1)
    {
        return (double) t * lg->vy_pix;
    }
    long span  = lg->nbframes - 1;
    long cycle = t / span;
    long rem   = t % span;
    return (double) cycle * lg->traj_y[span] + lg->traj_y[rem];
}
