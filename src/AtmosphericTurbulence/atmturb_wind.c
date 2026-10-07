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
            double k = 2.678 * 2.0 * M_PI * r / length_total * Lwind;
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
