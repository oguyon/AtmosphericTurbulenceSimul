// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmmod_refraction.c
 * @brief   Atmospheric refraction ray tracing and transmission calculation
 */

#include "AtmosphereModel.h"
#include "atmmod_types.h"

/**
 * atmmod_trace_upward_ray - Trace single upward ray through atmospheric layers
 * @lambda: Wavelength in meters.
 * @Zangle0: Ground zenith angle trial in radians.
 * @write_path: 1 to write trajectory to refractpath.txt, 0 otherwise.
 * @flux_out: Output transmission flux.
 *
 * Return: Deflection offset angle (alpha - Zangle0) in radians.
 */
static double atmmod_trace_upward_ray(double lambda, double Zangle0, int write_path,
                                      double *flux_out)
{
    const double lstep = 20.0;         // step [m]
    const double Re = 6371000.0;       // Earth radius [m]

    double x0 = 0.0, y0 = 0.0;
    double alpha = Zangle0;
    double h0 = SiteAlt;
    double n0 = 1.0 + AtmosphereModel_stdAtmModel_N((float)h0, (float)lambda, 0);

    FILE *fp = NULL;
    if (write_path)
    {
        fp = fopen("refractpath.txt", "w");
    }

    double h1 = 0.0;
    double flux = 1.0;
    while (h1 < 99000.0)
    {
        double x1 = x0 + lstep * sin(alpha);
        double y1 = y0 + lstep * cos(alpha);
        h1 = sqrt(x1 * x1 + (y1 + SiteAlt + Re) * (y1 + SiteAlt + Re)) - Re;
        double n1 = 1.0 + AtmosphereModel_stdAtmModel_N((float)h1, (float)lambda, 0);

        double alphae = atan2(x1, y1 + SiteAlt + Re);
        flux *= exp(-lstep * v_ABSCOEFF);

        double alpha0 = alpha - alphae;
        double alpha1 = asin(n0 * sin(alpha0) / n1);
        alpha = alpha1 + alphae;

        n0 = n1;
        x0 = x1;
        y0 = y1;
        h0 = h1;

        if (fp != NULL)
        {
            fprintf(fp, "%10.3f %10.3f %10.3f  %12.9f  %20g    %20g   %20g\n",
                    h0, x0, y0, alpha, alpha - Zangle0,
                    (alpha - Zangle0) / M_PI * 180.0 * 3600.0, y0 * tan(Zangle0) - x0);
        }
    }

    if (fp != NULL)
    {
        fclose(fp);
    }

    *flux_out = flux;
    return alpha - Zangle0;
}

/**
 * AtmosphereModel_RefractionPath - Compute ray trajectory through atmospheric layers
 * @lambda: Optical wavelength in meters.
 * @Zangle: Zenith angle at ground in radians.
 * @WritePath: Flag (1 to write refractpath.txt, 0 otherwise).
 *
 * Return: Atmospheric refraction deflection in arcseconds.
 */
double AtmosphereModel_RefractionPath(double lambda, double Zangle, int WritePath)
{
    double Zangle0 = Zangle;
    double errV = 100000.0;
    double flux = 1.0;
    long iter = 0;

    while (iter < 20 && errV > 0.00001)
    {
        int should_write = (WritePath && (iter == 0 || errV <= 0.00002));
        double offsetangle = atmmod_trace_upward_ray(lambda, Zangle0, should_write, &flux);
        Zangle0 -= offsetangle;
        errV = fabs(offsetangle / M_PI * 180.0 * 3600.0);
        iter++;
    }

    double defl_arcsec = (Zangle - Zangle0) / M_PI * 180.0 * 3600.0;
    printf("%10.6f um   %12.10f     Atmospheric Refraction = %10.6f arcsec   ",
           lambda * 1e6, 1.0 + AtmosphereModel_stdAtmModel_N(SiteAlt, (float)lambda, 0),
           defl_arcsec);
    printf("TRANSMISSION = %lf\n", flux);

    v_TRANSM = flux;
    return defl_arcsec;
}
