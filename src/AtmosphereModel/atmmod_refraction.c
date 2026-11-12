// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmmod_refraction.c
 * @brief   Atmospheric refraction ray tracing and transmission calculation
 */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>

#include "AtmosphereModel.h"
#include "atmmod_types.h"

/**
 * atmmod_trace_upward_ray - Trace single upward ray through atmospheric layers
 * @lambda:     Wavelength in meters.
 * @Zangle0:    Ground zenith angle trial in radians.
 * @write_path: 1 to write trajectory to refractpath.txt, 0 otherwise.
 * @flux_out:   Output transmission flux.
 *
 * Return: Deflection offset angle (alpha - Zangle0) in radians.
 */
static double atmmod_trace_upward_ray(
    double  lambda,
    double  Zangle0,
    int     write_path,
    double *flux_out)
{
    const double lstep   = 20.0;
    const double Re      = 6371000.0;
    const double site_re = (double) SiteAlt + Re;
    const float  lam_f   = (float) lambda;

    double x0    = 0.0;
    double y0    = 0.0;
    double alpha = Zangle0;
    double h0    = SiteAlt;
    double n0    = 1.0 + AtmosphereModel_stdAtmModel_N((float) h0, lam_f, 0);

    FILE *fp = write_path ? fopen("refractpath.txt", "w") : NULL;
    const double tan_z0      = (fp != NULL) ? tan(Zangle0) : 0.0;
    const double rad2arcsec  = 180.0 * 3600.0 / M_PI;

    double h1      = 0.0;
    double abs_sum = 0.0;
    while (h1 < 99000.0)
    {
        double sa, ca;
        sincos(alpha, &sa, &ca);
        double x1 = x0 + lstep * sa;
        double y1 = y0 + lstep * ca;
        double yr = y1 + site_re;
        h1        = sqrt(x1 * x1 + yr * yr) - Re;
        double n1 = 1.0 + AtmosphereModel_stdAtmModel_N((float) h1, lam_f, 0);

        double alphae = atan2(x1, yr);
        abs_sum += v_ABSCOEFF;

        double alpha0 = alpha - alphae;
        alpha         = asin(n0 * sin(alpha0) / n1) + alphae;

        n0 = n1;
        x0 = x1;
        y0 = y1;
        h0 = h1;

        if (fp != NULL)
        {
            double da = alpha - Zangle0;
            fprintf(fp, "%10.3f %10.3f %10.3f  %12.9f  %20g    %20g   %20g\n",
                    h0, x0, y0, alpha, da, da * rad2arcsec, y0 * tan_z0 - x0);
        }
    }

    if (fp != NULL)
    {
        fclose(fp);
    }

    *flux_out = exp(-lstep * abs_sum);
    return alpha - Zangle0;
}

/**
 * AtmosphereModel_RefractionPath - Compute ray trajectory through atmospheric layers
 * @lambda:    Optical wavelength in meters.
 * @Zangle:    Zenith angle at ground in radians.
 * @WritePath: Flag (1 to write refractpath.txt, 0 otherwise).
 *
 * Return: Atmospheric refraction deflection in arcseconds.
 */
double AtmosphereModel_RefractionPath(
    double lambda,
    double Zangle,
    int    WritePath)
{
    double       Zangle0    = Zangle;
    double       errV       = 100000.0;
    double       flux       = 1.0;
    long         iter       = 0;
    const double rad2arcsec = 180.0 * 3600.0 / M_PI;

    while (iter < 20 && errV > 0.00001)
    {
        int    should_write = (WritePath && (iter == 0 || errV <= 0.00002));
        double offsetangle  = atmmod_trace_upward_ray(lambda, Zangle0, should_write, &flux);
        double resid        = (Zangle0 + offsetangle) - Zangle;
        Zangle0 -= resid;
        errV = fabs(resid * rad2arcsec);
        iter++;
    }

    double defl_arcsec = (Zangle - Zangle0) * rad2arcsec;
    printf("%10.6f um   %12.10f     Atmospheric Refraction = %10.6f arcsec   ",
           lambda * 1e6, 1.0 + AtmosphereModel_stdAtmModel_N(SiteAlt, (float) lambda, 0),
           defl_arcsec);
    printf("TRANSMISSION = %lf\n", flux);

    v_TRANSM = flux;
    return defl_arcsec;
}
