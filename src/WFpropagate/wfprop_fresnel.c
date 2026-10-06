// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    wfprop_fresnel.c
 * @brief   2D Fresnel diffractive wavefront propagation algorithms
 */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "fft/fft.h"
#include "WFpropagate.h"
#include "atmturb_compat.h"

#ifndef PI
#    define PI 3.14159265358979323846264338328
#endif

#define SBUFFERSIZE 2000

/**
 * Fresnel_propagate_wavefront - Fresnel propagate complex optical field
 * @in: Name of input complex image
 * @out: Name of output complex image
 * @PUPIL_SCALE: Physical pixel scale [m/pixel]
 * @z: Propagation distance [m]
 * @lambda: Optical wavelength [m]
 *
 * Propagates a 2D optical field across distance z using Fourier Fresnel kernel.
 *
 * Return: 0 on success.
 */
int Fresnel_propagate_wavefront(char *in, char *out, double PUPIL_SCALE, double z, double lambda)
{
    do2dfft(in, "tmp");
    permut("tmp");
    long ID    = image_ID("tmp");
    int  atype = dcimg[ID].md[0].atype;

    long   naxes[2];
    naxes[0]     = dcimg[ID].md[0].size[0];
    naxes[1]     = dcimg[ID].md[0].size[1];
    double coeff = PI * z * lambda / (PUPIL_SCALE * naxes[0]) / (PUPIL_SCALE * naxes[0]);
    double co1   = 1.0 * naxes[0] * naxes[1];
    long   n0h   = naxes[0] / 2;

    if (atype == COMPLEX_FLOAT)
    {
        #pragma omp parallel for
        for (long jj = 0; jj < naxes[1]; jj++)
        {
            long jj1 = naxes[0] * jj;
            long jj2 = (jj - naxes[1] / 2) * (jj - naxes[1] / 2);
            for (long ii = 0; ii < naxes[0]; ii++)
            {
                long   ii1    = jj1 + ii;
                long   ii2    = ii - n0h;
                double sqdist = (double) (ii2 * ii2 + jj2);
                float  angle  = (float) (-coeff * sqdist);
                float  s, c;
                sincosf(angle, &s, &c);
                float re = (float) (dcimg[ID].array.CF[ii1].re / co1);
                float im = (float) (dcimg[ID].array.CF[ii1].im / co1);
                dcimg[ID].array.CF[ii1].re = re * c - im * s;
                dcimg[ID].array.CF[ii1].im = re * s + im * c;
            }
        }
    }
    else
    {
        #pragma omp parallel for
        for (long jj = 0; jj < naxes[1]; jj++)
        {
            long jj1 = naxes[0] * jj;
            long jj2 = (jj - naxes[1] / 2) * (jj - naxes[1] / 2);
            for (long ii = 0; ii < naxes[0]; ii++)
            {
                long   ii1    = jj1 + ii;
                long   ii2    = ii - n0h;
                double sqdist = (double) (ii2 * ii2 + jj2);
                double angle  = -coeff * sqdist;
                double s, c;
                sincos(angle, &s, &c);
                double re = dcimg[ID].array.CD[ii1].re / co1;
                double im = dcimg[ID].array.CD[ii1].im / co1;
                dcimg[ID].array.CD[ii1].re = re * c - im * s;
                dcimg[ID].array.CD[ii1].im = re * s + im * c;
            }
        }
    }

    permut("tmp");
    do2dffti("tmp", out);
    delete_image_ID("tmp");

    return 0;
}

/**
 * Init_Fresnel_propagate_wavefront - Initialize anti-aliased Fresnel transfer kernel
 * @Cim: Output complex kernel image name
 * @size: Grid dimension in pixels
 * @PUPIL_SCALE: Physical pixel scale [m/pixel]
 * @z: Propagation distance [m]
 * @lambda: Optical wavelength [m]
 * @FPMASKRAD: Focal plane mask radius cutoff
 * @Precision: 0 for single precision float, 1 for double precision
 *
 * Pre-computes the quadratic phase transfer function for Fresnel propagation.
 *
 * Return: 0 on success.
 */
int Init_Fresnel_propagate_wavefront(char *Cim, long size, double PUPIL_SCALE, double z,
                                    double lambda, double FPMASKRAD, int Precision)
{
    long ID;
    if (Precision == 0)
    {
        create_2DCimage_ID(Cim, size, size);
    }
    else
    {
        create_2DCimage_ID_double(Cim, size, size);
    }

    ID = image_ID(Cim);
    double coeff = PI * z * lambda / (PUPIL_SCALE * size) / (PUPIL_SCALE * size);
    double co1   = 1.0 * size * size;
    long   n0h   = size / 2;

    if (Precision == 0)
    {
        for (long jj = 0; jj < size; jj++)
        {
            long jj2 = (jj - n0h) * (jj - n0h);
            for (long ii = 0; ii < size; ii++)
            {
                long   ii2    = ii - n0h;
                double sqdist = (double) (ii2 * ii2 + jj2);
                double Pha    = -coeff * sqdist;
                double dist   = sqrt(sqdist);
                double Amp;

                if (dist < FPMASKRAD - 1.0)
                {
                    Amp = 1.0;
                }
                else if (dist < FPMASKRAD)
                {
                    Amp = 0.5 * (1.0 + cos(PI * (dist - FPMASKRAD + 1.0)));
                }
                else
                {
                    Amp = 0.0;
                }

                dcimg[ID].array.CF[jj * size + ii].re = (float) (Amp * cos(Pha) / co1);
                dcimg[ID].array.CF[jj * size + ii].im = (float) (Amp * sin(Pha) / co1);
            }
        }
    }
    else
    {
        for (long jj = 0; jj < size; jj++)
        {
            long jj2 = (jj - n0h) * (jj - n0h);
            for (long ii = 0; ii < size; ii++)
            {
                long   ii2    = ii - n0h;
                double sqdist = (double) (ii2 * ii2 + jj2);
                double Pha    = -coeff * sqdist;
                double dist   = sqrt(sqdist);
                double Amp;

                if (dist < FPMASKRAD - 1.0)
                {
                    Amp = 1.0;
                }
                else if (dist < FPMASKRAD)
                {
                    Amp = 0.5 * (1.0 + cos(PI * (dist - FPMASKRAD + 1.0)));
                }
                else
                {
                    Amp = 0.0;
                }

                dcimg[ID].array.CD[jj * size + ii].re = Amp * cos(Pha) / co1;
                dcimg[ID].array.CD[jj * size + ii].im = Amp * sin(Pha) / co1;
            }
        }
    }

    permut(Cim);
    return 0;
}

/**
 * Fresnel_propagate_wavefront1 - Apply precomputed transfer function kernel
 * @in: Name of input complex optical field
 * @out: Name of output complex optical field
 * @Cin: Name of precomputed complex transfer kernel
 *
 * Propagates field by multiplying Fourier transform with precomputed transfer function.
 *
 * Return: 0 on success.
 */
int Fresnel_propagate_wavefront1(char *in, char *out, char *Cin)
{
    char fname[SBUFFERSIZE];
    long ID     = image_ID(in);
    long sizein = dcimg[ID].md[0].size[0];
    snprintf(fname, sizeof(fname), "tmpfpw%ld", sizein);
    int atype = dcimg[ID].md[0].atype;

    do2dfft(in, fname);

    ID = image_ID(fname);
    long naxes[2];
    naxes[0]   = dcimg[ID].md[0].size[0];
    naxes[1]   = dcimg[ID].md[0].size[1];
    long nbelem = naxes[0] * naxes[1];
    long IDref  = image_ID(Cin);

    if (atype == COMPLEX_FLOAT)
    {
        for (long ii = 0; ii < nbelem; ii++)
        {
            double re    = dcimg[ID].array.CF[ii].re;
            double im    = dcimg[ID].array.CF[ii].im;
            double reref = dcimg[IDref].array.CF[ii].re;
            double imref = dcimg[IDref].array.CF[ii].im;

            dcimg[ID].array.CF[ii].re = (float) (re * reref - im * imref);
            dcimg[ID].array.CF[ii].im = (float) (re * imref + im * reref);
        }
    }
    else
    {
        for (long ii = 0; ii < nbelem; ii++)
        {
            double re    = dcimg[ID].array.CD[ii].re;
            double im    = dcimg[ID].array.CD[ii].im;
            double reref = dcimg[IDref].array.CD[ii].re;
            double imref = dcimg[IDref].array.CD[ii].im;

            dcimg[ID].array.CD[ii].re = re * reref - im * imref;
            dcimg[ID].array.CD[ii].im = re * imref + im * reref;
        }
    }

    do2dffti(fname, out);
    delete_image_ID(fname);

    return 0;
}
