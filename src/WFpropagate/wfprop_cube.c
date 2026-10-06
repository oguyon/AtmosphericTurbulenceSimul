// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    wfprop_cube.c
 * @brief   Multi-distance 3D Fresnel propagation cube generator
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "WFpropagate.h"
#include "atmturb_compat.h"

/**
 * Fresnel_propagate_cube - Propagate optical field through a series of distances
 * @IDcin_name: Input complex field name
 * @IDout_name_amp: Output amplitude 3D cube image name
 * @IDout_name_pha: Output phase 3D cube image name
 * @PUPIL_SCALE: Physical pixel scale [m/pixel]
 * @zstart: Starting propagation distance [m]
 * @zend: Ending propagation distance [m]
 * @NBzpts: Number of distance sampling points
 * @lambda: Optical wavelength [m]
 *
 * Computes 3D amplitude and phase cubes sampling propagation over a range of z.
 *
 * Return: 0 on success.
 */
long Fresnel_propagate_cube(char *IDcin_name, char *IDout_name_amp, char *IDout_name_pha,
                           double PUPIL_SCALE, double zstart, double zend, long NBzpts,
                           double lambda)
{
    long IDcin = image_ID(IDcin_name);
    long xsize = dcimg[IDcin].md[0].size[0];
    long ysize = dcimg[IDcin].md[0].size[1];
    int  atype = dcimg[IDcin].md[0].atype;

    long IDouta, IDoutp;
    if (atype == COMPLEX_FLOAT)
    {
        IDouta = create_3Dimage_ID(IDout_name_amp, xsize, ysize, NBzpts);
        IDoutp = create_3Dimage_ID(IDout_name_pha, xsize, ysize, NBzpts);
    }
    else
    {
        IDouta = create_3Dimage_ID_double(IDout_name_amp, xsize, ysize, NBzpts);
        IDoutp = create_3Dimage_ID_double(IDout_name_pha, xsize, ysize, NBzpts);
    }

    for (long kk = 0; kk < NBzpts; kk++)
    {
        double zprop = zstart + (zend - zstart) * (double) kk / (double) NBzpts;
        printf("[%ld] propagating by %f m\n", kk, zprop);
        Fresnel_propagate_wavefront(IDcin_name, "_propim", PUPIL_SCALE, zprop, lambda);
        long IDtmp = image_ID("_propim");

        if (atype == COMPLEX_FLOAT)
        {
            for (long ii = 0; ii < xsize; ii++)
            {
                for (long jj = 0; jj < ysize; jj++)
                {
                    double re  = dcimg[IDtmp].array.CF[jj * xsize + ii].re;
                    double im  = dcimg[IDtmp].array.CF[jj * xsize + ii].im;
                    float  amp = (float) sqrt(re * re + im * im);
                    float  pha = (float) atan2(im, re);
                    dcimg[IDouta].array.F[kk * xsize * ysize + jj * xsize + ii] = amp;
                    dcimg[IDoutp].array.F[kk * xsize * ysize + jj * xsize + ii] = pha;
                }
            }
        }
        else
        {
            for (long ii = 0; ii < xsize; ii++)
            {
                for (long jj = 0; jj < ysize; jj++)
                {
                    double re  = dcimg[IDtmp].array.CD[jj * xsize + ii].re;
                    double im  = dcimg[IDtmp].array.CD[jj * xsize + ii].im;
                    double amp = sqrt(re * re + im * im);
                    double pha = atan2(im, re);
                    dcimg[IDouta].array.D[kk * xsize * ysize + jj * xsize + ii] = amp;
                    dcimg[IDoutp].array.D[kk * xsize * ysize + jj * xsize + ii] = pha;
                }
            }
        }

        delete_image_ID("_propim");
    }

    return 0;
}
