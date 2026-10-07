// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_wind.c
 * @brief   1D von Karman turbulent wind velocity profile generation
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>

#include "AtmosphericTurbulence.h"
#include "atmturb_types.h"

/**
 * struct atmturb_wind_spec_t - Parameters shared by the three wind components
 * @vksize: Number of samples in the 1D series.
 * @pixscale: Spatial sampling step in meters.
 * @sigmawind: Target RMS velocity standard deviation in m/s.
 * @Lwind: Wind velocity outer scale in meters.
 * @seed: Resolved (non-zero) base RNG seed.
 */
typedef struct
{
    long     vksize;
    float    pixscale;
    float    sigmawind;
    float    Lwind;
    uint64_t seed;
} atmturb_wind_spec_t;

/**
 * atmturb_synthesize_wind_component - Synthesize 1D von Karman wind velocity component
 * @spec: Shared wind synthesis parameters.
 * @component: 0 = longitudinal (u), 1 = transverse (v), 2 = vertical (w).
 * @out: Output pointer to destination float buffer (length spec->vksize).
 */
static void atmturb_synthesize_wind_component(
    const atmturb_wind_spec_t *spec,
    int                        component,
    float                     *out)
{
    long vksize = spec->vksize;
    fftwf_complex *buf = (fftwf_complex *)fftwf_alloc_complex(vksize);
    if (!buf)
    {
        return;
    }

    double length_total = (double)vksize * spec->pixscale;
    double Lwind = spec->Lwind;
    // independent RNG stream per component
    uint64_t rng = atmturb_rng_stream_seed(spec->seed, (uint64_t)component);

    for (long ii = 0; ii < vksize; ii++)
    {
        long fx = (ii < vksize / 2) ? ii : (ii - vksize);
        double r = fabs((double)fx);

        if (fx == 0)
        {
            buf[ii][0] = 0.0f;
            buf[ii][1] = 0.0f;
            continue;
        }

        double amp = 0.0;
        if (component == 0)
        {
            double k = 1.339 * 2.0 * M_PI * r / length_total * Lwind;
            amp = sqrt(1.0 / pow(1.0 + k * k, 5.0 / 6.0));
        }
        else
        {
            double k = 1.339 * 2.0 * M_PI * r / length_total * Lwind;
            double num = 1.0 + (8.0 / 3.0) * (k * k);
            double den = pow(1.0 + k * k, 11.0 / 6.0);
            amp = sqrt(num / den);
        }

        double g0, g1;
        atmturb_rng_gaussian_pair(&rng, &g0, &g1);
        buf[ii][0] = (float)(amp * g0);
        buf[ii][1] = (float)(amp * g1);
    }

    fftwf_plan plan = fftwf_plan_dft_1d((int)vksize, buf, buf, FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    double sum_sq = 0.0;
    for (long ii = 0; ii < vksize; ii++)
    {
        sum_sq += (double)buf[ii][0] * (double)buf[ii][0];
    }
    double rms = sqrt(sum_sq / (double)vksize);
    double scale = (rms > 0.0) ? ((double)spec->sigmawind / rms) : 1.0;

    for (long ii = 0; ii < vksize; ii++)
    {
        out[ii] = (float)((double)buf[ii][0] * scale);
    }

    fftwf_free(buf);
}

/**
 * make_AtmosphericTurbulence_vonKarmanWind - Generate 3D velocity cube [u, v, w]
 * @vKsize: Sample length of 1D series.
 * @pixscale: Physical sampling step in meters.
 * @sigmawind: Velocity standard deviation in m/s.
 * @Lwind: Wind velocity turbulence outer scale in meters.
 * @seed: RNG seed (0 = time-based). The three components use independent streams.
 * @IDout_name: Output 3D image name (vKsize x 1 x 3).
 *
 * Return: Output image ID on success.
 */
long make_AtmosphericTurbulence_vonKarmanWind(
    long        vKsize,
    float       pixscale,
    float       sigmawind,
    float       Lwind,
    long        seed,
    const char *IDout_name)
{
    delete_image_ID(IDout_name);
    imageID IDc = create_3Dimage_ID(IDout_name, vKsize, 1, 3);

    atmturb_wind_spec_t spec = {
        .vksize = vKsize,
        .pixscale = pixscale,
        .sigmawind = sigmawind,
        .Lwind = Lwind,
        .seed = atmturb_resolve_seed((uint64_t)seed)
    };

    printf("vK wind outer scale = %f m\n", Lwind);
    printf("pixscale            = %f m\n", pixscale);
    printf("Image size          = %f m\n", vKsize * pixscale);
    printf("seed                = %llu\n", (unsigned long long)spec.seed);

    // Longitudinal (u), tangential (v) and vertical (w) components
    for (int comp = 0; comp < 3; comp++)
    {
        atmturb_synthesize_wind_component(&spec, comp, &dcimg[IDc].array.F[comp * vKsize]);
    }

    return IDc;
}

/**
 * atmturb_wind_synthesize_trajectory - Synthesize cumulative 2D trajectory with turbulent wind
 * @params: Trajectory generation input parameters.
 * @traj_x: Output buffer of size nbframes for cumulative x offset [pixels].
 * @traj_y: Output buffer of size nbframes for cumulative y offset [pixels].
 *
 * Return: 0 on success, -1 on invalid arguments or memory allocation failure.
 */
int atmturb_wind_synthesize_trajectory(
    const atmturb_wind_traj_params_t *params,
    double                           *traj_x,
    double                           *traj_y)
{
    if (params == NULL || traj_x == NULL || traj_y == NULL)
    {
        return -1;
    }
    if (params->nbframes <= 0 || !(params->dt_s > 0.0) || !(params->dx_master_m > 0.0))
    {
        return -1;
    }

    traj_x[0] = 0.0;
    traj_y[0] = 0.0;
    if (params->nbframes == 1)
    {
        return 0;
    }

    if (!(params->sigma_wind_mps > 0.0) || !(params->L_wind_m > 0.0))
    {
        for (long t = 1; t < params->nbframes; t++)
        {
            traj_x[t] = (double) t * params->vx_pix;
            traj_y[t] = (double) t * params->vy_pix;
        }
        return 0;
    }

    float *u_fluc = (float *) malloc(sizeof(float) * (size_t) params->nbframes);
    float *v_fluc = (float *) malloc(sizeof(float) * (size_t) params->nbframes);
    if (u_fluc == NULL || v_fluc == NULL)
    {
        free(u_fluc);
        free(v_fluc);
        return -1;
    }

    double v_mean_pix = sqrt(params->vx_pix * params->vx_pix + params->vy_pix * params->vy_pix);
    double v_mean_mps = v_mean_pix * params->dx_master_m / params->dt_s;
    double ev_x = (v_mean_pix > 1e-6) ? (params->vx_pix / v_mean_pix) : 1.0;
    double ev_y = (v_mean_pix > 1e-6) ? (params->vy_pix / v_mean_pix) : 0.0;
    double ep_x = -ev_y;
    double ep_y =  ev_x;
    if (v_mean_mps < 0.1)
    {
        v_mean_mps = 0.1;
    }

    double dx_step = v_mean_mps * params->dt_s;
    atmturb_wind_spec_t spec = {
        .vksize    = params->nbframes,
        .pixscale  = (float) dx_step,
        .sigmawind = (float) params->sigma_wind_mps,
        .Lwind     = (float) params->L_wind_m,
        .seed      = atmturb_resolve_seed(params->seed)
    };
    atmturb_synthesize_wind_component(&spec, 0, u_fluc);
    atmturb_synthesize_wind_component(&spec, 1, v_fluc);

    double scale_pix = params->dt_s / params->dx_master_m;
    for (long t = 1; t < params->nbframes; t++)
    {
        double u = v_mean_mps + (double) u_fluc[t - 1];
        double v = (double) v_fluc[t - 1];
        double step_x = (u * ev_x + v * ep_x) * scale_pix;
        double step_y = (u * ev_y + v * ep_y) * scale_pix;
        traj_x[t] = traj_x[t - 1] + step_x;
        traj_y[t] = traj_y[t - 1] + step_y;
    }

    free(u_fluc);
    free(v_fluc);
    return 0;
}

