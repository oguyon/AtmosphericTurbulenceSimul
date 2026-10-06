// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_wind.c
 * @brief   1D von Karman turbulent wind velocity profile generation
 */

#include "AtmosphericTurbulence.h"
#include "atmturb_types.h"

/**
 * atmturb_compute_wind_psd_amplitude - Compute 1D von Karman spectrum amplitude
 * @vksize: Array size in samples.
 * @pixscale: Physical sampling step in meters.
 * @Lwind: Wind velocity outer scale in meters.
 * @is_transverse: 0 for longitudinal component, 1 for transverse / vertical.
 */
static void atmturb_compute_wind_psd_amplitude(long vksize, float pixscale,
                                               float Lwind, int is_transverse)
{
    imageID ID = create_2Dimage_ID("tmpamp0", vksize, 1);
    double length_total = vksize * pixscale;

    for (long ii = 0; ii < vksize; ii++)
    {
        double dx = 1.0 * ii - vksize / 2;
        double r = fabs(dx);

        if (!is_transverse)
        {
            double k = 1.339 * 2.0 * M_PI * r / length_total * Lwind;
            dcimg[ID].array.F[ii] = (float)sqrt(1.0 / pow(1.0 + k * k, 5.0 / 6.0));
        }
        else
        {
            double k = 2.678 * 2.0 * M_PI * r / length_total * Lwind;
            double num = 1.0 + (8.0 / 3.0) * (k * k);
            double den = pow(1.0 + k * k, 11.0 / 6.0);
            dcimg[ID].array.F[ii] = (float)sqrt(num / den);
        }
    }
}

/**
 * atmturb_synthesize_wind_series - Transform PSD into time series and scale to sigma
 * @vksize: Number of samples.
 * @sigma: Target RMS velocity standard deviation.
 * @idc: Destination 3D image ID.
 * @channel_offset: Output offset in 3D array (e.g. 0, vksize, 2*vksize).
 */
static void atmturb_synthesize_wind_series(long vksize, float sigma, imageID idc,
                                          long channel_offset)
{
    make_rnd("tmppha0", vksize, 1, "");
    arith_image_cstmult("tmppha0", 2.0 * PI, "tmppha");
    delete_image_ID("tmppha0");

    make_rnd("tmpg", vksize, 1, "-gauss");
    arith_image_mult("tmpg", "tmpamp0", "tmpamp");
    delete_image_ID("tmpamp0");
    delete_image_ID("tmpg");

    arith_set_pixel("tmpamp", 0.0, vksize / 2, 0);
    mk_complex_from_amph("tmpamp", "tmppha", "tmpc", 0);
    delete_image_ID("tmpamp");
    delete_image_ID("tmppha");

    permut("tmpc");
    do2dfft("tmpc", "tmpcf");
    delete_image_ID("tmpc");

    mk_reim_from_complex("tmpcf", "tmpo1", "tmpo2", 0);
    delete_image_ID("tmpcf");
    delete_image_ID("tmpo2");

    imageID ID = image_ID("tmpo1");
    double rms = 0.0;
    for (long ii = 0; ii < vksize; ii++)
    {
        rms += dcimg[ID].array.F[ii] * dcimg[ID].array.F[ii];
    }
    rms = sqrt(rms / vksize);

    for (long ii = 0; ii < vksize; ii++)
    {
        dcimg[idc].array.F[channel_offset + ii] =
            (float)(dcimg[ID].array.F[ii] / rms * sigma);
    }
    delete_image_ID("tmpo1");
}

/**
 * make_AtmosphericTurbulence_vonKarmanWind - Generate 3D velocity cube [u, v, w]
 * @vKsize: Sample length of 1D series.
 * @pixscale: Physical sampling step in meters.
 * @sigmawind: Velocity standard deviation in m/s.
 * @Lwind: Wind velocity turbulence outer scale in meters.
 * @size: Unused legacy size parameter.
 * @IDout_name: Output 3D image name (vKsize x 1 x 3).
 *
 * Return: Output image ID on success.
 */
long make_AtmosphericTurbulence_vonKarmanWind(long vKsize, float pixscale,
                                             float sigmawind, float Lwind,
                                             long size, char *IDout_name)
{
    (void)size;
    imageID IDc = create_3Dimage_ID(IDout_name, vKsize, 1, 3);

    printf("vK wind outer scale = %f m\n", Lwind);
    printf("pixscale            = %f m\n", pixscale);
    printf("Image size          = %f m\n", vKsize * pixscale);

    // Longitudinal component (u)
    atmturb_compute_wind_psd_amplitude(vKsize, pixscale, Lwind, 0);
    atmturb_synthesize_wind_series(vKsize, sigmawind, IDc, 0);

    // Tangential component (v)
    atmturb_compute_wind_psd_amplitude(vKsize, pixscale, Lwind, 1);
    atmturb_synthesize_wind_series(vKsize, sigmawind, IDc, vKsize);

    // Vertical component (w)
    atmturb_compute_wind_psd_amplitude(vKsize, pixscale, Lwind, 1);
    atmturb_synthesize_wind_series(vKsize, sigmawind, IDc, 2 * vKsize);

    return IDc;
}
