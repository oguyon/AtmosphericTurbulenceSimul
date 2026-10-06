// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    wfprop_lyot.c
 * @brief   Lyot coronagraph propagation and contrast testing routines
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CLIcore.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "fft/fft.h"
#include "image_gen/image_gen.h"
#include "WFpropagate.h"
#include "atmturb_compat.h"

/**
 * WFpropagate_TestLyot - Test multi-mask Lyot coronagraph propagation
 * @NBmask: Number of sequential masks
 * @maskpos: Array of axial mask positions [m]
 *
 * Propagates pupil through multiple masks and computes focal plane contrast.
 *
 * Return: Average contrast inside the focal plane dark hole.
 */
double WFpropagate_TestLyot(long NBmask, double *maskpos)
{
    const double lambda   = 0.55e-6;
    const double pixscale = 5.0693766e-5;
    const double rin      = 3.0;  // in l/D
    const double rout     = 10.0; // in l/D
    double       z        = 0.0;
    double       value    = 0.0;
    double       valuecnt = 0.0;
    char         fname[200];

    printf("Testing Lyot Masks\n");
    copy_image_ID("imc", "imc0", 0);

    for (long k = 0; k < NBmask; k++)
    {
        snprintf(fname, sizeof(fname), "mask%ld", k);
        long IDm = image_ID(fname);
        Fresnel_propagate_wavefront("imc0", "imc1", pixscale, maskpos[k] - z, lambda);
        z = maskpos[k];
        delete_image_ID("imc0");
        mk_amph_from_complex("imc1", "ima1", "imp1", 0);
        delete_image_ID("imc1");

        long ID   = image_ID("ima1");
        long size = dcimg[ID].md[0].size[0];
        snprintf(fname, sizeof(fname), "!ima_%ld_0.fits", k);
        save_fl_fits("ima1", fname);
        snprintf(fname, sizeof(fname), "!imp_%ld_0.fits", k);
        save_fl_fits("imp1", fname);

        for (long ii = 0; ii < size * size; ii++)
        {
            dcimg[ID].array.F[ii] *= dcimg[IDm].array.F[ii];
        }

        snprintf(fname, sizeof(fname), "!ima_%ld_1.fits", k);
        save_fl_fits("ima1", fname);
        snprintf(fname, sizeof(fname), "!imp_%ld_1.fits", k);
        save_fl_fits("imp1", fname);
        mk_complex_from_amph("ima1", "imp1", "imc0", 0);
        delete_image_ID("ima1");
        delete_image_ID("imp1");
    }

    mk_amph_from_complex("imc0", "pup0a", "pup0p", 0);
    delete_image_ID("pup0p");
    save_fl_fits("pup0a", "!pup0a.fits");
    delete_image_ID("pup0a");

    permut("imc0");
    do2dfft("imc0", "imc1");
    delete_image_ID("imc0");
    permut("imc1");
    mk_amph_from_complex("imc1", "foca", "focp", 0);
    delete_image_ID("imc1");
    execute_arith("foci=foca*foca/98130");

    long IDa  = image_ID("foca");
    long size = dcimg[IDa].md[0].size[0];
    for (long ii = 0; ii < size; ii++)
    {
        for (long jj = 0; jj < size; jj++)
        {
            double x = ((double) ii - 0.5 * size) / 5.12;
            double y = ((double) jj - 0.5 * size) / 5.12;
            double r = sqrt(x * x + y * y);
            if ((r > 5.0 * rout) || (r < rin))
            {
                dcimg[IDa].array.F[jj * size + ii] = 0.0f;
            }
        }
    }

    mk_complex_from_amph("foca", "focp", "focc", 0);
    delete_image_ID("focp");
    delete_image_ID("foca");
    permut("focc");
    do2dfft("focc", "pupc1");
    delete_image_ID("focc");
    permut("pupc1");
    mk_amph_from_complex("pupc1", "pupa1", "pupp1", 0);
    delete_image_ID("pupc1");
    save_fl_fits("pupa1", "!pupa_res.fits");
    delete_image_ID("pupa1");
    delete_image_ID("pupp1");

    long IDf = image_ID("foci");
    size     = dcimg[IDf].md[0].size[0];
    for (long ii = 0; ii < size; ii++)
    {
        for (long jj = 0; jj < size; jj++)
        {
            double x = ((double) ii - 0.5 * size) / 5.12;
            double y = ((double) jj - 0.5 * size) / 5.12;
            double r = sqrt(x * x + y * y);
            if ((r > rin) && (r < rout))
            {
                value += dcimg[IDf].array.F[jj * size + ii];
                valuecnt += 1.0;
            }
        }
    }

    return (valuecnt > 0.0) ? (value / valuecnt) : 0.0;
}

/**
 * WFpropagate_run - Standalone test harness for Lyot propagation
 *
 * Sets up 4-mask coronagraph geometry and computes focal plane contrast.
 *
 * Return: 0 on success.
 */
long WFpropagate_run(void)
{
    long    NBmask  = 4;
    double *maskpos = (double *) malloc(sizeof(double) * NBmask);
    if (!maskpos)
    {
        return -1;
    }

    load_fits("pa1a_post2.fits", "pa1a", 1);
    execute_arith("refpup=pa1a*pa1a");
    double tot0 = arith_image_total("refpup");

    make_disk("mask0", 2048, 2048, 1024, 1024, 190.0);
    save_fl_fits("mask0", "!mask0.fits");
    maskpos[0] = -1.35;

    load_fits("mask_i100_o5.fits", "mask1", 1);
    maskpos[1] = 0.0;

    load_fits("mask_i40_o5.fits", "mask2", 1);
    maskpos[2] = 0.08;

    load_fits("mask_o5_r48.fits", "mask3", 1);
    maskpos[3] = 0.14;

    execute_arith("refpup1=refpup*mask0*mask1*mask2*mask3");
    double tot1 = arith_image_total("refpup1");

    FILE *fp = fopen("result.txt", "w");
    if (fp)
    {
        fclose(fp);
    }

    double x     = 0.0;
    double value = WFpropagate_TestLyot(NBmask, maskpos);
    save_fl_fits("foci", "!foci.fits");
    delete_image_ID("foci");

    printf("AVERAGE CONTRAST = %g\n", value);
    printf("MASK THROUGHPUT = %g\n", (tot0 > 0.0) ? (tot1 / tot0) : 0.0);

    fp = fopen("result.txt", "a");
    if (fp)
    {
        fprintf(fp, "%f %g %g\n", x, value, (tot0 > 0.0) ? (tot1 / tot0) : 0.0);
        fclose(fp);
    }

    free(maskpos);
    return 0;
}
