// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cube_ops.c
 * @brief   Wavefront cube spatial contraction and binning operations
 */

#include "AtmosphericTurbulence.h"
#include "atmturb_types.h"

/**
 * atmturb_bin_complex_wavefront - Spatial 2D boxcar binning of complex wavefront cube
 * @in_amp: Input amplitude buffer.
 * @in_pha: Input phase buffer.
 * @out_amp: Output binned amplitude buffer.
 * @out_pha: Output binned phase buffer.
 * @nx: Input X dimension.
 * @ny: Input Y dimension.
 * @nz: Number of Z slices / frames.
 * @pfactor: Decimation factor along X and Y.
 */
static void atmturb_bin_complex_wavefront(
    const float *in_amp,
    const float *in_pha,
    float       *out_amp,
    float       *out_pha,
    long         nx,
    long         ny,
    long         nz,
    int          pfactor)
{
    long out_nx = nx / pfactor;
    long out_ny = ny / pfactor;
    long LARGE = 10000;

    for (long kk = 0; kk < nz; kk++)
    {
        for (long ii = 0; ii < out_nx; ii++)
        {
            for (long jj = 0; jj < out_ny; jj++)
            {
                float re = 0.0f, im = 0.0f;
                float pharef = 0.0f, ampref = 0.0f;

                for (long i = 0; i < pfactor; i++)
                {
                    for (long j = 0; j < pfactor; j++)
                    {
                        long in_idx = kk * nx * ny + (pfactor * jj + j) * nx + (pfactor * ii + i);
                        float a = in_amp[in_idx];
                        float p = in_pha[in_idx];
                        pharef += a * p;
                        ampref += a;
                        re += a * cosf(p);
                        im += a * sinf(p);
                    }
                }

                float a_out = sqrtf(re * re + im * im);
                float p_out = atan2f(im, re);
                if (ampref > 0.0001f)
                {
                    pharef /= ampref;
                }
                float P = 2.0f * (float)PI *
                          ((long)(0.5f + 1.0f * LARGE
                                  + (pharef - p_out) / (2.0f * (float)PI)) - LARGE);
                if (ampref < 0.01f)
                {
                    P = 0.0f;
                }

                long out_idx = kk * out_nx * out_ny + jj * out_nx + ii;
                out_pha[out_idx] = p_out + P;
                out_amp[out_idx] = a_out / (pfactor * pfactor);
            }
        }
    }
}

/**
 * atmturb_bin_phaseonly_wavefront - Spatial 2D boxcar binning of phase-only wavefront cube
 * @in_pha: Input phase buffer.
 * @out_pha: Output binned phase buffer.
 * @nx: Input X dimension.
 * @ny: Input Y dimension.
 * @nz: Number of Z slices / frames.
 * @pfactor: Decimation factor along X and Y.
 */
static void atmturb_bin_phaseonly_wavefront(
    const float *in_pha,
    float       *out_pha,
    long         nx,
    long         ny,
    long         nz,
    int          pfactor)
{
    long out_nx = nx / pfactor;
    long out_ny = ny / pfactor;
    long LARGE = 10000;

    for (long kk = 0; kk < nz; kk++)
    {
        for (long ii = 0; ii < out_nx; ii++)
        {
            for (long jj = 0; jj < out_ny; jj++)
            {
                float re = 0.0f, im = 0.0f;
                float pharef = 0.0f;

                for (long i = 0; i < pfactor; i++)
                {
                    for (long j = 0; j < pfactor; j++)
                    {
                        long in_idx = kk * nx * ny + (pfactor * jj + j) * nx + (pfactor * ii + i);
                        float p = in_pha[in_idx];
                        pharef += p;
                        re += cosf(p);
                        im += sinf(p);
                    }
                }

                float p_out = atan2f(im, re);
                pharef /= (pfactor * pfactor);
                float P = 2.0f * (float)PI *
                          ((long)(0.5f + 1.0f * LARGE
                                  + (pharef - p_out) / (2.0f * (float)PI)) - LARGE);

                long out_idx = kk * out_nx * out_ny + jj * out_nx + ii;
                out_pha[out_idx] = p_out + P;
            }
        }
    }
}

/**
 * contract_wavefront_cube - Decimate complex wavefront cube by 2^factor
 * @ina_file: Input amplitude FITS file.
 * @inp_file: Input phase FITS file.
 * @outa_file: Output amplitude FITS file.
 * @outp_file: Output phase FITS file.
 * @factor: Power of 2 decimation exponent.
 *
 * Return: 0 on success.
 */
int contract_wavefront_cube(
    const char *ina_file,
    const char *inp_file,
    const char *outa_file,
    const char *outp_file,
    int         factor)
{
    int pfactor = 1 << factor;

    load_fits(inp_file, "tmpwfp", 1);
    imageID IDpha = image_ID("tmpwfp");
    load_fits(ina_file, "tmpwfa", 1);
    imageID IDamp = image_ID("tmpwfa");

    long nx = dcimg[IDpha].md[0].size[0];
    long ny = dcimg[IDpha].md[0].size[1];
    long nz = dcimg[IDpha].md[0].size[2];

    imageID IDoutpha = create_3Dimage_ID("tmpwfop", nx / pfactor, ny / pfactor, nz);
    imageID IDoutamp = create_3Dimage_ID("tmpwfoa", nx / pfactor, ny / pfactor, nz);

    atmturb_bin_complex_wavefront(dcimg[IDamp].array.F, dcimg[IDpha].array.F,
                                 dcimg[IDoutamp].array.F, dcimg[IDoutpha].array.F,
                                 nx, ny, nz, pfactor);

    save_fl_fits("tmpwfop", outp_file);
    save_fl_fits("tmpwfoa", outa_file);

    delete_image_ID("tmpwfa");
    delete_image_ID("tmpwfp");
    delete_image_ID("tmpwfoa");
    delete_image_ID("tmpwfop");

    return 0;
}

/**
 * contract_wavefront_cube_phaseonly - Decimate phase-only wavefront cube by 2^factor
 * @inp_file: Input phase FITS file.
 * @outp_file: Output phase FITS file.
 * @factor: Power of 2 decimation exponent.
 *
 * Return: 0 on success.
 */
int contract_wavefront_cube_phaseonly(
    const char *inp_file,
    const char *outp_file,
    int         factor)
{
    int pfactor = 1 << factor;

    load_fits(inp_file, "tmpwfp", 1);
    imageID IDpha = image_ID("tmpwfp");

    long nx = dcimg[IDpha].md[0].size[0];
    long ny = dcimg[IDpha].md[0].size[1];
    long nz = dcimg[IDpha].md[0].size[2];

    imageID IDoutpha = create_3Dimage_ID("tmpwfop", nx / pfactor, ny / pfactor, nz);

    atmturb_bin_phaseonly_wavefront(dcimg[IDpha].array.F, dcimg[IDoutpha].array.F,
                                    nx, ny, nz, pfactor);

    save_fl_fits("tmpwfop", outp_file);
    delete_image_ID("tmpwfp");
    delete_image_ID("tmpwfop");

    return 0;
}

/**
 * contract_wavefront_series - Decimate series of complex wavefront cubes by factor of 2
 * @in_prefix: Input filename prefix.
 * @out_prefix: Output filename prefix.
 * @NB_files: Number of files in sequence.
 *
 * Return: 0 on success.
 */
int contract_wavefront_series(
    const char *in_prefix,
    const char *out_prefix,
    long        NB_files)
{
    char fname_p[200], fname_a[200];
    const double SLAMBDA = 1.65e-6;

    for (long index = 0; index < NB_files; index++)
    {
        printf("INDEX = %ld/%ld\n", index, NB_files);
        snprintf(fname_p, sizeof(fname_p), "%s%08ld.%09ld.pha.fits",
                 in_prefix, index, (long)(1.0e12 * SLAMBDA + 0.5));
        snprintf(fname_a, sizeof(fname_a), "%s%08ld.%09ld.amp.fits",
                 in_prefix, index, (long)(1.0e12 * SLAMBDA + 0.5));

        load_fits(fname_p, "tmpwfp", 1);
        imageID IDpha = image_ID("tmpwfp");
        load_fits(fname_a, "tmpwfa", 1);
        imageID IDamp = image_ID("tmpwfa");

        long nx = dcimg[IDpha].md[0].size[0];
        long ny = dcimg[IDpha].md[0].size[1];
        long nz = dcimg[IDpha].md[0].size[2];

        imageID IDoutpha = create_3Dimage_ID("tmpwfop", nx / 2, ny / 2, nz);
        imageID IDoutamp = create_3Dimage_ID("tmpwfoa", nx / 2, ny / 2, nz);

        atmturb_bin_complex_wavefront(dcimg[IDamp].array.F, dcimg[IDpha].array.F,
                                     dcimg[IDoutamp].array.F, dcimg[IDoutpha].array.F,
                                     nx, ny, nz, 2);

        snprintf(fname_p, sizeof(fname_p), "%s%08ld.%09ld.pha.fits",
                 out_prefix, index, (long)(1.0e12 * SLAMBDA + 0.1));
        snprintf(fname_a, sizeof(fname_a), "%s%08ld.%09ld.amp.fits",
                 out_prefix, index, (long)(1.0e12 * SLAMBDA + 0.5));

        save_fl_fits("tmpwfop", fname_p);
        save_fl_fits("tmpwfoa", fname_a);

        delete_image_ID("tmpwfa");
        delete_image_ID("tmpwfp");
        delete_image_ID("tmpwfoa");
        delete_image_ID("tmpwfop");
    }

    return 0;
}
