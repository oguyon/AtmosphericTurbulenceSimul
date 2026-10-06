// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_screens.c
 * @brief   Master atmospheric turbulence phase screen generation
 */

#include "AtmosphericTurbulence.h"
#include "atmturb_types.h"

/**
 * atmturb_generate_spatial_distance_map - Build 2D spatial frequency distance map
 * @outname: Output image name.
 * @size: Grid linear dimension in pixels.
 * @outer_f0: Spatial frequency corresponding to outer scale.
 * @rlim: Inner frequency cutoff radius.
 * @rlim_mode: 1 to zero out frequencies below rlim.
 * @precision: 0 for single precision, 1 for double precision.
 */
static void atmturb_generate_spatial_distance_map(const char *outname, long size,
                                                  double outer_f0, double rlim,
                                                  int rlim_mode, long precision)
{
    imageID ID = (precision == 0) ? create_2Dimage_ID(outname, size, size)
                                  : create_2Dimage_ID_double(outname, size, size);

    double rlim2 = rlim * rlim;
    double f0_sq = outer_f0 * outer_f0;

    #pragma omp parallel for
    for (long jj = 0; jj < size; jj++)
    {
        double dy = 1.0 * jj - size / 2;
        double dy2 = dy * dy;
        long row_idx = jj * size;

        for (long ii = 0; ii < size; ii++)
        {
            double dx = 1.0 * ii - size / 2;
            double r2 = dx * dx + dy2;
            double val = (rlim_mode == 1 && r2 < rlim2) ? 0.0 : sqrt(r2 + f0_sq);

            if (precision == 0)
            {
                dcimg[ID].array.F[row_idx + ii] = (float)val;
            }
            else
            {
                dcimg[ID].array.D[row_idx + ii] = val;
            }
        }
    }
}

/**
 * atmturb_apply_innerscale_cutoff - Apply Gaussian inner scale attenuation in Fourier domain
 * @imname: Target image name.
 * @size: Grid linear dimension in pixels.
 * @inner_f0: Inner scale cutoff frequency.
 * @precision: 0 for single precision, 1 for double precision.
 */
static void atmturb_apply_innerscale_cutoff(const char *imname, long size,
                                           double inner_f0, long precision)
{
    imageID ID = image_ID(imname);
    double inv_two_inner2 = 0.5 / (inner_f0 * inner_f0);

    #pragma omp parallel for
    for (long jj = 0; jj < size; jj++)
    {
        double dy = 1.0 * jj - size / 2;
        double dy2 = dy * dy;
        long row_idx = jj * size;

        for (long ii = 0; ii < size; ii++)
        {
            double dx = 1.0 * ii - size / 2;
            double factor = exp(-(dx * dx + dy2) * inv_two_inner2);

            if (precision == 0)
            {
                dcimg[ID].array.F[row_idx + ii] *= (float)factor;
            }
            else
            {
                dcimg[ID].array.D[row_idx + ii] *= factor;
            }
        }
    }
}

/**
 * atmturb_measure_structure_constant - Compute Kolmogorov structure function scaling constant
 * @imname: Input phase screen image name.
 * @size: Grid linear dimension in pixels.
 * @power_exponent: Power law exponent (e.g. 5/3 for Kolmogorov).
 *
 * Return: Structure function scaling constant C.
 */
static double atmturb_measure_structure_constant(const char *imname, long size,
                                                 double power_exponent)
{
    fft_structure_function((char *)imname, "strf");
    imageID ID = image_ID("strf");

    double value = 0.0;
    long cnt = 0;
    long Dlim = 3;

    if (dcimg[ID].md[0].atype == FLOAT)
    {
        for (long ii = 1; ii < Dlim; ii++)
        {
            for (long jj = 1; jj < Dlim; jj++)
            {
                value += log10(dcimg[ID].array.F[jj * size + ii]) -
                         power_exponent * log10(sqrt(ii * ii + jj * jj));
                cnt++;
            }
        }
    }
    else
    {
        for (long ii = 1; ii < Dlim; ii++)
        {
            for (long jj = 1; jj < Dlim; jj++)
            {
                value += log10(dcimg[ID].array.D[jj * size + ii]) -
                         power_exponent * log10(sqrt(ii * ii + jj * jj));
                cnt++;
            }
        }
    }
    delete_image_ID("strf");
    return pow(10.0, value / cnt);
}

/**
 * make_master_turbulence_screen - Generate von Karman master turbulence screens
 * @ID_name1: Output name for screen 1.
 * @ID_name2: Output name for screen 2.
 * @size: Grid dimension in pixels.
 * @outerscale: Outer scale in pixels.
 * @innerscale: Inner scale in pixels.
 * @WFprecision: Precision flag (0=single, 1=double).
 *
 * Return: 0 on success.
 */
int make_master_turbulence_screen(char *ID_name1, char *ID_name2, long size,
                                 float outerscale, float innerscale, long WFprecision)
{
    printf("Make turbulence screen, precision = %ld\n", WFprecision);
    fflush(stdout);

    double rlim = 0.0;
    int RLIMMODE = 0;
    imageID IDv = variable_ID("RLIM");
    if (IDv != -1)
    {
        RLIMMODE = 1;
        rlim = dcvar[IDv].value.f;
        printf("R limit = %f pix\n", rlim);
    }

    double outer_f0 = 1.0 * size / outerscale;
    double inner_f0 = (5.92 / (2.0 * M_PI)) * size / innerscale;

    if (WFprecision == 0)
    {
        make_rnd("tmppha", size, size, "");
        make_rnd("tmpg", size, size, "-gauss");
    }
    else
    {
        make_rnd_double("tmppha", size, size, "");
        make_rnd_double("tmpg", size, size, "-gauss");
    }

    arith_image_cstmult("tmppha", 2.0 * PI, "tmppha1");
    delete_image_ID("tmppha");

    atmturb_generate_spatial_distance_map("tmpd", size, outer_f0, rlim, RLIMMODE, WFprecision);
    atmturb_apply_innerscale_cutoff("tmpg", size, inner_f0, WFprecision);

    arith_image_cstpow("tmpd", 11.0 / 6.0, "tmpd1");
    delete_image_ID("tmpd");
    arith_image_div("tmpg", "tmpd1", "tmpamp");
    delete_image_ID("tmpg");
    delete_image_ID("tmpd1");

    arith_set_pixel("tmpamp", 0.0, size / 2, size / 2);
    mk_complex_from_amph("tmpamp", "tmppha1", "tmpc", 0);
    delete_image_ID("tmpamp");
    delete_image_ID("tmppha1");

    permut("tmpc");
    do2dfft("tmpc", "tmpcf");
    delete_image_ID("tmpc");
    mk_reim_from_complex("tmpcf", "tmpo1", "tmpo2", 0);
    delete_image_ID("tmpcf");

    double C1 = atmturb_measure_structure_constant("tmpo1", size, 5.0 / 3.0);
    double C2 = atmturb_measure_structure_constant("tmpo2", size, 5.0 / 3.0);
    printf("C1, C2 =   %f %f\n", C1, C2);

    arith_image_cstmult("tmpo1", 1.0 / sqrt(C1), ID_name1);
    arith_image_cstmult("tmpo2", 1.0 / sqrt(C2), ID_name2);
    delete_image_ID("tmpo1");
    delete_image_ID("tmpo2");

    return 0;
}

/**
 * make_master_turbulence_screen_pow - Generate phase screen with arbitrary power-law PSD
 * @ID_name1: Output name for screen 1.
 * @ID_name2: Output name for screen 2.
 * @size: Grid dimension in pixels.
 * @power: Spatial frequency power-law exponent.
 *
 * Return: 0 on success.
 */
int make_master_turbulence_screen_pow(char *ID_name1, char *ID_name2, long size, float power)
{
    make_rnd("tmppha", size, size, "");
    arith_image_cstmult("tmppha", 2.0 * PI, "tmppha1");
    delete_image_ID("tmppha");

    make_dist("tmpd", size, size, size / 2, size / 2);
    make_rnd("tmpg", size, size, "-gauss");

    arith_image_cstpow("tmpd", power, "tmpd1");
    delete_image_ID("tmpd");
    arith_image_div("tmpg", "tmpd1", "tmpamp");
    delete_image_ID("tmpg");
    delete_image_ID("tmpd1");

    arith_set_pixel("tmpamp", 0.0, size / 2, size / 2);
    mk_complex_from_amph("tmpamp", "tmppha1", "tmpc", 0);
    delete_image_ID("tmpamp");
    delete_image_ID("tmppha1");

    permut("tmpc");
    do2dfft("tmpc", "tmpcf");
    delete_image_ID("tmpc");
    mk_reim_from_complex("tmpcf", "tmpo1", "tmpo2", 0);
    delete_image_ID("tmpcf");

    double C1 = atmturb_measure_structure_constant("tmpo1", size, power);
    double C2 = atmturb_measure_structure_constant("tmpo2", size, power);

    arith_image_cstmult("tmpo1", 1.0 / sqrt(C1), ID_name1);
    arith_image_cstmult("tmpo2", 1.0 / sqrt(C2), ID_name2);
    delete_image_ID("tmpo1");
    delete_image_ID("tmpo2");

    return 0;
}
