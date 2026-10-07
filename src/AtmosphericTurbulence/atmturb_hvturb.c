// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_hvturb.c
 * @brief   Hufnagel-Valley Cn2 atmospheric turbulence profile generator
 */

#include "AtmosphericTurbulence.h"
#include "atmturb_types.h"

/**
 * atmturb_solve_hv_ground_amplitude - Fit HV ground boundary layer coefficient to match r0
 * @wspeed: High-altitude wind speed [m/s].
 * @r0: Target Fried parameter in meters.
 * @sitealt: Observatory altitude in meters.
 * @hmax: Upper integration altitude in meters.
 * @cn2sum_out: Output integrated Cn2 value.
 * @r0val_out: Output computed r0 value.
 *
 * Return: Solved boundary layer amplitude A.
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
    const double A0 = 1.7e-14;
    const double hstep = 1.0;

    double Acoeff = 1.0;
    double Astep = 1.0;
    double A = A0;
    double r0val = 0.0;
    double cn2sum = 0.0;

    for (long iter = 0; iter < 30; iter++)
    {
        A = A0 * Acoeff;
        cn2sum = 0.0;
        for (double h = sitealt; h < hmax; h += hstep)
        {
            double cn2 = 5.94e-53 * pow(wspeed / 27.0, 2.0) * pow(h, 10.0) * exp(-h / 1000.0) +
                         2.7e-16 * exp(-h / 1500.0) + A * exp(-(h - sitealt) / 100.0);
            cn2sum += cn2 / hstep;
        }
        r0val = 1.0 / pow(0.423 * pow(2.0 * M_PI / lambda, 2.0) * cn2sum, 3.0 / 5.0);

        if (r0val > r0)
        {
            Acoeff *= 1.0 + Astep;
        }
        else
        {
            Acoeff /= 1.0 + Astep;
        }
        Astep *= 0.8;
    }

    *cn2sum_out = cn2sum;
    *r0val_out = r0val;
    return A;
}

/**
 * AtmosphericTurbulence_makeHV_CN2prof - Generate Hufnagel-Valley Cn2 profile file
 * @wspeed: Upper atmospheric wind speed [m/s].
 * @r0: Target Fried parameter in meters.
 * @sitealt: Observatory elevation in meters.
 * @NBlayer: Number of vertical discrete layers.
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
    const double hmax = 30000.0;
    const double lambda = 0.55e-6;

    double cn2sum = 0.0;
    double r0val = 0.0;
    double A = atmturb_solve_hv_ground_amplitude(wspeed, r0, sitealt, hmax, &cn2sum, &r0val);

    FILE *fp = fopen("conf_turb.txt", "w");
    if (fp != NULL)
    {
        fprintf(fp, "%f\n", (lambda / r0val) / M_PI * 180.0 * 3600.0);
        fclose(fp);
    }

    double *layer_h = malloc(sizeof(double) * NBlayer);
    double *layer_cn2 = calloc(NBlayer, sizeof(double));

    for (long k = 0; k < NBlayer; k++)
    {
        layer_h[k] = sitealt + pow(1.0 * k / (NBlayer - 1), 2.0) * (hmax - sitealt);
    }

    for (double h = sitealt; h < hmax; h += 1.0)
    {
        double cn2 = 5.94e-53 * pow(wspeed / 27.0, 2.0) * pow(h, 10.0) * exp(-h / 1000.0) +
                     2.7e-16 * exp(-h / 1500.0) + A * exp(-(h - sitealt) / 100.0);
        long k = (long)(sqrt((h - sitealt) / (hmax - sitealt)) * (1.0 * NBlayer - 1.0) + 0.5);
        if (k >= 0 && k < NBlayer)
        {
            layer_cn2[k] += cn2;
        }
    }

    fp = fopen(outfile, "w");
    if (fp == NULL)
    {
        free(layer_h);
        free(layer_cn2);
        return -1;
    }

    fprintf(fp, "# altitude(m)   relativeCN2     speed(m/s)   direction(rad) "
                "outerscale[m] innerscale[m] sigmaWsp[m/s] Lwind[m]\n\n");

    for (long k = 0; k < NBlayer; k++)
    {
        layer_cn2[k] /= cn2sum;
        double l0 = 0.008 + 0.072 * pow(layer_h[k] / 20000.0, 1.6);
        double L0 = (layer_h[k] < 14000.0)
                        ? pow(10.0, 2.0 - 0.9 * (layer_h[k] / 14000.0))
                        : pow(10.0, 1.1 + 0.3 * (layer_h[k] - 14000.0) / 6000.0);

        double wsp = wspeed * (0.3 + 0.8 * sqrt(1.0 * k / (1.0 + NBlayer)));
        double sigma_wsp = 0.1 * wsp;

        fprintf(fp, "%12f  %12f  %12f  %12f  %12f  %12f  %12f  %12f\n",
                layer_h[k], layer_cn2[k], wsp, 2.0 * M_PI * k / (1.0 + NBlayer),
                L0, l0, sigma_wsp, 500.0);
    }
    fclose(fp);

    free(layer_h);
    free(layer_cn2);
    return 0;
}
