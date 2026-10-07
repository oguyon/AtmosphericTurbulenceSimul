// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_hvturb.c
 * @brief   Hufnagel-Valley Cn2 atmospheric turbulence profile generator
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "AtmosphericTurbulence.h"
#include "atmturb_types.h"

/**
 * atmturb_hv_cn2 - Hufnagel-Valley 5/7 refractive index structure constant
 * @h: Altitude above sea level [m].
 * @wspeed: High-altitude (200 mbar) wind speed [m/s].
 * @A: Ground boundary layer amplitude [m^-2/3].
 * @sitealt: Observatory altitude [m] (origin of the ground layer term).
 *
 * Return: Cn2(h) [m^-2/3].
 */
static double atmturb_hv_cn2(
    double h,
    double wspeed,
    double A,
    double sitealt)
{
    return 5.94e-53 * pow(wspeed / 27.0, 2.0) * pow(h, 10.0) * exp(-h / 1000.0) +
           2.7e-16 * exp(-h / 1500.0) + A * exp(-(h - sitealt) / 100.0);
}

/**
 * atmturb_solve_hv_ground_amplitude - Analytically solve HV ground layer amplitude to match r0
 * @wspeed: High-altitude wind speed [m/s].
 * @r0: Target Fried parameter in meters.
 * @sitealt: Observatory altitude in meters.
 * @hmax: Upper integration altitude in meters.
 * @cn2sum_out: Output integrated Cn2 value [m^1/3].
 * @r0val_out: Output computed r0 value [m].
 *
 * Return: Solved boundary layer amplitude A, or -1.0 if r0 is unreachable.
 */
static double atmturb_solve_hv_ground_amplitude(
    double  wspeed,
    double  r0,
    double  sitealt,
    double  hmax,
    double *cn2sum_out,
    double *r0val_out)
{
    const double lambda = 0.55e-6;
    const double hstep  = 1.0;
    double k            = 2.0 * M_PI / lambda;
    double target_cn2   = pow(r0, -5.0 / 3.0) / (0.423 * k * k);

    double jupper = 0.0;
    for (double h = sitealt; h < hmax; h += hstep)
    {
        double term1 = 5.94e-53 * pow(wspeed / 27.0, 2.0) * pow(h, 10.0) * exp(-h / 1000.0);
        double term2 = 2.7e-16 * exp(-h / 1500.0);
        jupper += (term1 + term2) * hstep;
    }

    double iground = 100.0 * (1.0 - exp(-(hmax - sitealt) / 100.0));
    double A       = (target_cn2 - jupper) / iground;

    if (A < 0.0)
    {
        double r0_max = pow(0.423 * k * k * jupper, -3.0 / 5.0);
        printf("ERROR: requested r0 = %.4f m exceeds maximum reachable r0 = %.4f m "
               "from upper atmosphere alone\n", r0, r0_max);
        return -1.0;
    }

    *cn2sum_out = target_cn2;
    *r0val_out  = r0;
    return A;
}

/**
 * atmturb_hv_compute_moments - Calculate isoplanatic angle theta0 and coherence time tau0
 * @wspeed: High-altitude wind speed [m/s].
 * @A: Ground boundary layer amplitude.
 * @sitealt: Observatory altitude [m].
 * @hmax: Upper integration altitude [m].
 * @wind_model: Wind model selector (0 = legacy, 1 = Bufton).
 * @target_cn2: Total integrated Cn2.
 * @r0: Target Fried parameter in meters.
 * @theta0_out: Output isoplanatic angle [rad].
 * @tau0_out: Output coherence time [s].
 */
static void atmturb_hv_compute_moments(
    double  wspeed,
    double  A,
    double  sitealt,
    double  hmax,
    int     wind_model,
    double  target_cn2,
    double  r0,
    double *theta0_out,
    double *tau0_out)
{
    const double hstep = 1.0;
    double sum_h       = 0.0;
    double sum_v       = 0.0;

    for (double h = sitealt; h < hmax; h += hstep)
    {
        double cn2   = atmturb_hv_cn2(h, wspeed, A, sitealt) * hstep;
        double h_rel = h - sitealt;
        sum_h += cn2 * pow(h_rel, 5.0 / 3.0);

        double v;
        if (wind_model == 1)
        {
            v = 5.0 + wspeed * exp(-pow((h - 9400.0) / 4800.0, 2.0));
        }
        else
        {
            v = wspeed * (0.3 + 0.8 * sqrt((h - sitealt) / (hmax - sitealt)));
        }
        sum_v += cn2 * pow(v, 5.0 / 3.0);
    }

    double h_eff = pow(sum_h / target_cn2, 3.0 / 5.0);
    double v_eff = pow(sum_v / target_cn2, 3.0 / 5.0);

    *theta0_out = (h_eff > 0.0) ? (0.314 * r0 / h_eff) : 0.0;
    *tau0_out   = (v_eff > 0.0) ? (0.314 * r0 / v_eff) : 0.0;
}

/**
 * atmturb_hv_bin_layers - Integrate the HV profile into discrete layers
 * @wspeed: High-altitude wind speed [m/s].
 * @A: Ground boundary layer amplitude.
 * @sitealt: Observatory altitude [m].
 * @hmax: Upper integration altitude [m].
 * @NBlayer: Number of layers (>= 1).
 * @layer_h: Output layer altitudes [m] (length NBlayer).
 * @layer_cn2: Output integrated Cn2 per layer [m^1/3] (length NBlayer, zeroed).
 *
 * Layers are spaced quadratically in height above the site. With a single layer, its altitude
 * is the Cn2-weighted mean altitude.
 */
static void atmturb_hv_bin_layers(
    double  wspeed,
    double  A,
    double  sitealt,
    double  hmax,
    long    NBlayer,
    double *layer_h,
    double *layer_cn2)
{
    const double hstep   = 1.0;
    double       hmoment = 0.0;

    for (long k = 0; k < NBlayer; k++)
    {
        double frac = (NBlayer > 1) ? (1.0 * (double) k / (double) (NBlayer - 1)) : 0.0;
        layer_h[k]  = sitealt + frac * frac * (hmax - sitealt);
    }

    for (double h = sitealt; h < hmax; h += hstep)
    {
        double cn2 = atmturb_hv_cn2(h, wspeed, A, sitealt) * hstep;
        long k = (long) (sqrt((h - sitealt) / (hmax - sitealt)) * (double) (NBlayer - 1) + 0.5);
        if (k >= 0 && k < NBlayer)
        {
            layer_cn2[k] += cn2;
        }
        hmoment += cn2 * h;
    }

    if (NBlayer == 1 && layer_cn2[0] > 0.0)
    {
        layer_h[0] = hmoment / layer_cn2[0];
    }
}

/**
 * atmturb_hv_write_profile - Write discrete HV layers to a turbulence profile file
 * @outfile: Destination profile file path.
 * @wspeed: High-altitude wind speed [m/s].
 * @NBlayer: Number of layers.
 * @layer_h: Layer altitudes [m].
 * @layer_cn2: Relative Cn2 per layer (normalized).
 * @wind_model: Wind model selector (0 = legacy, 1 = Bufton).
 * @seed: Resolved PRNG seed.
 *
 * Return: 0 on success, -1 if the file cannot be opened.
 */
static int atmturb_hv_write_profile(
    const char   *outfile,
    double        wspeed,
    long          NBlayer,
    const double *layer_h,
    const double *layer_cn2,
    int           wind_model,
    uint64_t      seed)
{
    FILE *fp = fopen(outfile, "w");
    if (fp == NULL)
    {
        printf("ERROR: cannot write profile \"%s\"\n", outfile);
        return -1;
    }

    fprintf(fp, "# altitude(m)   relativeCN2     speed(m/s)   direction(rad) "
                "outerscale[m] innerscale[m] sigmaWsp[m/s] Lwind[m]\n\n");

    for (long k = 0; k < NBlayer; k++)
    {
        double l0 = 0.008 + 0.072 * pow(layer_h[k] / 20000.0, 1.6);
        double L0 = (layer_h[k] < 14000.0)
                        ? pow(10.0, 2.0 - 0.9 * (layer_h[k] / 14000.0))
                        : pow(10.0, 1.1 + 0.3 * (layer_h[k] - 14000.0) / 6000.0);

        double wsp;
        double dir;
        if (wind_model == 1)
        {
            wsp = 5.0 + wspeed * exp(-pow((layer_h[k] - 9400.0) / 4800.0, 2.0));
            uint64_t rng = atmturb_rng_stream_seed(seed, (uint64_t) (2000 + k));
            dir = (atmturb_rng_splitmix64(&rng) >> 11) * (2.0 * M_PI / 9007199254740992.0);
        }
        else
        {
            wsp = wspeed * (0.3 + 0.8 * sqrt((double) k / (1.0 + (double) NBlayer)));
            dir = 2.0 * M_PI * (double) k / (1.0 + (double) NBlayer);
        }
        double sigma_wsp = 0.1 * wsp;

        fprintf(fp, "%12f  %12f  %12f  %12f  %12f  %12f  %12f  %12f\n",
                layer_h[k], layer_cn2[k], wsp, dir, L0, l0, sigma_wsp, 500.0);
    }
    fclose(fp);
    return 0;
}

/**
 * AtmosphericTurbulence_makeHV_CN2prof_opt - Generate HV profile with model options
 * @wspeed: Upper atmospheric wind speed [m/s].
 * @r0: Target Fried parameter in meters.
 * @sitealt: Observatory elevation in meters.
 * @NBlayer: Number of vertical discrete layers (>= 1).
 * @outfile: Destination profile file path.
 * @wind_model: Wind velocity model (0 = legacy, 1 = Bufton).
 * @seed: Master RNG seed value.
 *
 * Return: 0 on success, -1 on failure.
 */
int AtmosphericTurbulence_makeHV_CN2prof_opt(
    double      wspeed,
    double      r0,
    double      sitealt,
    long        NBlayer,
    const char *outfile,
    int         wind_model,
    uint64_t    seed)
{
    const double hmax   = 30000.0;
    const double lambda = 0.55e-6;

    if (NBlayer < 1 || !(r0 > 0.0) || !(sitealt < hmax))
    {
        printf("ERROR: need NBlayer >= 1, r0 > 0 and sitealt < %.0f m\n", hmax);
        return -1;
    }

    double cn2sum = 0.0;
    double r0val  = 0.0;
    double A = atmturb_solve_hv_ground_amplitude(wspeed, r0, sitealt, hmax, &cn2sum, &r0val);
    if (A < 0.0)
    {
        return -1;
    }

    double theta0 = 0.0;
    double tau0   = 0.0;
    atmturb_hv_compute_moments(wspeed, A, sitealt, hmax, wind_model, cn2sum, r0val,
                               &theta0, &tau0);

    FILE *fp = fopen("conf_turb.txt", "w");
    if (fp != NULL)
    {
        fprintf(fp, "%f\n", 0.98 * (lambda / r0val) / M_PI * 180.0 * 3600.0);
        fclose(fp);
    }

    printf("[milkatmturb] HV profile: r0 = %.4f m (seeing = %.3f\"), "
           "theta0 = %.2f\" (%.1f urad), tau0 = %.2f ms\n",
           r0val, 0.98 * (lambda / r0val) / M_PI * 180.0 * 3600.0,
           theta0 * 180.0 * 3600.0 / M_PI, theta0 * 1e6, tau0 * 1000.0);

    int     ret       = -1;
    double *layer_h   = (double *) malloc(sizeof(double) * (size_t) NBlayer);
    double *layer_cn2 = (double *) calloc((size_t) NBlayer, sizeof(double));
    if (layer_h == NULL || layer_cn2 == NULL)
    {
        goto cleanup;
    }

    atmturb_hv_bin_layers(wspeed, A, sitealt, hmax, NBlayer, layer_h, layer_cn2);
    for (long k = 0; k < NBlayer; k++)
    {
        layer_cn2[k] /= cn2sum;
    }
    ret = atmturb_hv_write_profile(outfile, wspeed, NBlayer, layer_h, layer_cn2,
                                   wind_model, seed);

cleanup:
    free(layer_h);
    free(layer_cn2);
    return ret;
}

/**
 * AtmosphericTurbulence_makeHV_CN2prof - Generate standard Hufnagel-Valley Cn2 profile file
 * @wspeed: Upper atmospheric wind speed [m/s].
 * @r0: Target Fried parameter in meters.
 * @sitealt: Observatory elevation in meters.
 * @NBlayer: Number of vertical discrete layers (>= 1).
 * @outfile: Destination profile file path.
 *
 * Return: 0 on success, -1 on failure.
 */
int AtmosphericTurbulence_makeHV_CN2prof(
    double      wspeed,
    double      r0,
    double      sitealt,
    long        NBlayer,
    const char *outfile)
{
    return AtmosphericTurbulence_makeHV_CN2prof_opt(wspeed, r0, sitealt, NBlayer,
                                                    outfile, 1, 1);
}
