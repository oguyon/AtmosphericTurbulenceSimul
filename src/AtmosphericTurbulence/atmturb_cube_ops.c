// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cube_ops.c
 * @brief   Wavefront cube spatial contraction and binning operations
 */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "AtmosphericTurbulence.h"
#include "atmturb_fresnel_simd.h"
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
    long        out_nx    = nx / pfactor;
    long        out_ny    = ny / pfactor;
    long        LARGE     = 10000;
    const float inv_pf2   = 1.0f / (float) (pfactor * pfactor);
    const float inv_twopi = 1.0f / (2.0f * (float) PI);
    const float twopi     = 2.0f * (float) PI;

    #pragma omp parallel
    {
        fftwf_complex *row_z = fftwf_alloc_complex((size_t) nx);
        float         *acc   = (float *) malloc(sizeof(float) * (size_t) (4 * out_nx));
        if (row_z != NULL && acc != NULL)
        {
            float *restrict acc_re  = acc;
            float *restrict acc_im  = acc + out_nx;
            float *restrict acc_pr  = acc + 2 * out_nx;
            float *restrict acc_ar  = acc + 3 * out_nx;

            #pragma omp for schedule(static)
            for (long kk = 0; kk < nz; kk++)
            {
                long in_base  = kk * nx * ny;
                long out_base = kk * out_nx * out_ny;
                for (long jj = 0; jj < out_ny; jj++)
                {
                    memset(acc, 0, sizeof(float) * (size_t) (4 * out_nx));
                    for (int j = 0; j < pfactor; j++)
                    {
                        long         roff = in_base + (pfactor * jj + j) * nx;
                        const float *arow = &in_amp[roff];
                        const float *prow = &in_pha[roff];
                        atmturb_fresnel_simd_phasor_masked(row_z, prow, arow, 1.0f, nx);

                        for (long ii = 0; ii < out_nx; ii++)
                        {
                            long  x0 = (long) pfactor * ii;
                            float re = 0.0f;
                            float im = 0.0f;
                            float pr = 0.0f;
                            float ar = 0.0f;
                            for (int i = 0; i < pfactor; i++)
                            {
                                float a = arow[x0 + i];
                                pr += a * prow[x0 + i];
                                ar += a;
                                re += row_z[x0 + i][0];
                                im += row_z[x0 + i][1];
                            }
                            acc_re[ii] += re;
                            acc_im[ii] += im;
                            acc_pr[ii] += pr;
                            acc_ar[ii] += ar;
                        }
                    }

                    long ooff = out_base + jj * out_nx;
                    for (long ii = 0; ii < out_nx; ii++)
                    {
                        float re     = acc_re[ii];
                        float im     = acc_im[ii];
                        float ampref = acc_ar[ii];
                        float pharef = (ampref > 0.0001f) ? (acc_pr[ii] / ampref) : acc_pr[ii];
                        float a_out  = sqrtf(re * re + im * im);
                        float p_out  = atan2f(im, re);
                        float P      = 0.0f;
                        if (ampref >= 0.01f)
                        {
                            long kwrap = (long) (0.5f + (float) LARGE +
                                                 (pharef - p_out) * inv_twopi) - LARGE;
                            P          = twopi * (float) kwrap;
                        }
                        out_pha[ooff + ii] = p_out + P;
                        out_amp[ooff + ii] = a_out * inv_pf2;
                    }
                }
            }
        }
        free(acc);
        fftwf_free(row_z);
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
    long        out_nx    = nx / pfactor;
    long        out_ny    = ny / pfactor;
    long        LARGE     = 10000;
    const float inv_pf2   = 1.0f / (float) (pfactor * pfactor);
    const float inv_twopi = 1.0f / (2.0f * (float) PI);
    const float twopi     = 2.0f * (float) PI;

    #pragma omp parallel
    {
        fftwf_complex *row_z = fftwf_alloc_complex((size_t) nx);
        float         *acc   = (float *) malloc(sizeof(float) * (size_t) (3 * out_nx));
        if (row_z != NULL && acc != NULL)
        {
            float *restrict acc_re = acc;
            float *restrict acc_im = acc + out_nx;
            float *restrict acc_pr = acc + 2 * out_nx;

            #pragma omp for schedule(static)
            for (long kk = 0; kk < nz; kk++)
            {
                long in_base  = kk * nx * ny;
                long out_base = kk * out_nx * out_ny;
                for (long jj = 0; jj < out_ny; jj++)
                {
                    memset(acc, 0, sizeof(float) * (size_t) (3 * out_nx));
                    for (int j = 0; j < pfactor; j++)
                    {
                        const float *prow = &in_pha[in_base + (pfactor * jj + j) * nx];
                        atmturb_fresnel_simd_init_modulate_field(row_z, prow, nx);

                        for (long ii = 0; ii < out_nx; ii++)
                        {
                            long  x0 = (long) pfactor * ii;
                            float re = 0.0f;
                            float im = 0.0f;
                            float pr = 0.0f;
                            for (int i = 0; i < pfactor; i++)
                            {
                                pr += prow[x0 + i];
                                re += row_z[x0 + i][0];
                                im += row_z[x0 + i][1];
                            }
                            acc_re[ii] += re;
                            acc_im[ii] += im;
                            acc_pr[ii] += pr;
                        }
                    }

                    long ooff = out_base + jj * out_nx;
                    for (long ii = 0; ii < out_nx; ii++)
                    {
                        float pharef = acc_pr[ii] * inv_pf2;
                        float p_out  = atan2f(acc_im[ii], acc_re[ii]);
                        long  kwrap  = (long) (0.5f + (float) LARGE +
                                               (pharef - p_out) * inv_twopi) - LARGE;
                        out_pha[ooff + ii] = p_out + twopi * (float) kwrap;
                    }
                }
            }
        }
        free(acc);
        fftwf_free(row_z);
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
    if (factor < 0 || factor > 6)
    {
        return -1;
    }
    int pfactor = 1 << factor;

    delete_image_ID("tmpwfp");
    delete_image_ID("tmpwfa");
    load_fits(inp_file, "tmpwfp", 1);
    imageID IDpha = image_ID("tmpwfp");
    load_fits(ina_file, "tmpwfa", 1);
    imageID IDamp = image_ID("tmpwfa");

    if (IDpha < 0 || IDamp < 0)
    {
        delete_image_ID("tmpwfa");
        delete_image_ID("tmpwfp");
        return -1;
    }

    long nx = dcimg[IDpha].md[0].size[0];
    long ny = dcimg[IDpha].md[0].size[1];
    long nz = dcimg[IDpha].md[0].size[2];

    if (nx < pfactor || ny < pfactor || nz < 1)
    {
        delete_image_ID("tmpwfa");
        delete_image_ID("tmpwfp");
        return -1;
    }

    delete_image_ID("tmpwfop");
    delete_image_ID("tmpwfoa");
    imageID IDoutpha = create_3Dimage_ID("tmpwfop", nx / pfactor, ny / pfactor, nz);
    imageID IDoutamp = create_3Dimage_ID("tmpwfoa", nx / pfactor, ny / pfactor, nz);
    if (IDoutpha < 0 || IDoutamp < 0)
    {
        delete_image_ID("tmpwfa");
        delete_image_ID("tmpwfp");
        delete_image_ID("tmpwfoa");
        delete_image_ID("tmpwfop");
        return -1;
    }

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
    if (factor < 0 || factor > 6)
    {
        return -1;
    }
    int pfactor = 1 << factor;

    delete_image_ID("tmpwfp");
    load_fits(inp_file, "tmpwfp", 1);
    imageID IDpha = image_ID("tmpwfp");

    if (IDpha < 0)
    {
        return -1;
    }

    long nx = dcimg[IDpha].md[0].size[0];
    long ny = dcimg[IDpha].md[0].size[1];
    long nz = dcimg[IDpha].md[0].size[2];

    if (nx < pfactor || ny < pfactor || nz < 1)
    {
        delete_image_ID("tmpwfp");
        return -1;
    }

    delete_image_ID("tmpwfop");
    imageID IDoutpha = create_3Dimage_ID("tmpwfop", nx / pfactor, ny / pfactor, nz);
    if (IDoutpha < 0)
    {
        delete_image_ID("tmpwfp");
        return -1;
    }

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

        delete_image_ID("tmpwfp");
        delete_image_ID("tmpwfa");
        load_fits(fname_p, "tmpwfp", 1);
        imageID IDpha = image_ID("tmpwfp");
        load_fits(fname_a, "tmpwfa", 1);
        imageID IDamp = image_ID("tmpwfa");

        if (IDpha < 0 || IDamp < 0)
        {
            delete_image_ID("tmpwfa");
            delete_image_ID("tmpwfp");
            return -1;
        }

        long nx = dcimg[IDpha].md[0].size[0];
        long ny = dcimg[IDpha].md[0].size[1];
        long nz = dcimg[IDpha].md[0].size[2];

        if (nx < 2 || ny < 2 || nz < 1)
        {
            delete_image_ID("tmpwfa");
            delete_image_ID("tmpwfp");
            return -1;
        }

        delete_image_ID("tmpwfop");
        delete_image_ID("tmpwfoa");
        imageID IDoutpha = create_3Dimage_ID("tmpwfop", nx / 2, ny / 2, nz);
        imageID IDoutamp = create_3Dimage_ID("tmpwfoa", nx / 2, ny / 2, nz);
        if (IDoutpha < 0 || IDoutamp < 0)
        {
            delete_image_ID("tmpwfa");
            delete_image_ID("tmpwfp");
            delete_image_ID("tmpwfoa");
            delete_image_ID("tmpwfop");
            return -1;
        }

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
