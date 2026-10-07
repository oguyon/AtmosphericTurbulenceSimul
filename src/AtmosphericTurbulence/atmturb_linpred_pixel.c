// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_linpred_pixel.c
 * @brief   Pixel-level and shift-invariant kernel linear prediction filters
 */

#include "AtmosphericTurbulence.h"
#include "atmturb_types.h"

/**
 * AtmosphericTurbulence_LinPredictor_filt_2DKernelExtract - Extract shift-invariant 2D kernels
 * @IDfilt_name: Input filter matrix image name.
 * @IDmask_name: Active pupil mask image name.
 * @krad: Extraction neighborhood radius in pixels.
 * @IDkern_name: Output 3D kernel image name.
 *
 * Return: Output image ID on success.
 */
long AtmosphericTurbulence_LinPredictor_filt_2DKernelExtract(
    const char *IDfilt_name,
    const char *IDmask_name,
    long        krad,
    const char *IDkern_name)
{
    imageID IDfilt = image_ID(IDfilt_name);
    imageID IDmask = image_ID(IDmask_name);
    if (IDfilt < 0 || IDmask < 0)
    {
        return -1;
    }

    long nx = dcimg[IDmask].md[0].size[0];
    long ny = dcimg[IDmask].md[0].size[1];
    long nbpix = dcimg[IDfilt].md[0].size[0];
    long pforder = (dcimg[IDfilt].md[0].naxis > 2) ? dcimg[IDfilt].md[0].size[2] : 1;

    long *pix_x = malloc(sizeof(long) * nbpix);
    long *pix_y = malloc(sizeof(long) * nbpix);

    long count = 0;
    for (long ii = 0; ii < nx; ii++)
    {
        for (long jj = 0; jj < ny; jj++)
        {
            if (dcimg[IDmask].array.F[jj * nx + ii] > 0.5f && count < nbpix)
            {
                pix_x[count] = ii;
                pix_y[count] = jj;
                count++;
            }
        }
    }

    long ksize = 2 * krad + 1;
    imageID IDkern = create_3Dimage_ID(IDkern_name, ksize, ksize, pforder);
    imageID IDcnt = create_3Dimage_ID("kerncnt", ksize, ksize, pforder);

    for (long p = 0; p < count; p++)
    {
        long i0 = pix_x[p];
        long j0 = pix_y[p];

        for (long dt = 0; dt < pforder; dt++)
        {
            for (long q = 0; q < count; q++)
            {
                long di = pix_x[q] - i0;
                long dj = pix_y[q] - j0;

                if (labs(di) <= krad && labs(dj) <= krad)
                {
                    long ki = di + krad;
                    long kj = dj + krad;
                    long kidx = dt * ksize * ksize + kj * ksize + ki;
                    long fidx = dt * nbpix * nbpix + p * nbpix + q;

                    dcimg[IDkern].array.F[kidx] += dcimg[IDfilt].array.F[fidx];
                    dcimg[IDcnt].array.F[kidx] += 1.0f;
                }
            }
        }
    }

    for (long i = 0; i < ksize * ksize * pforder; i++)
    {
        if (dcimg[IDcnt].array.F[i] > 0.5f)
        {
            dcimg[IDkern].array.F[i] /= dcimg[IDcnt].array.F[i];
        }
    }

    delete_image_ID("kerncnt");
    free(pix_x);
    free(pix_y);
    return IDkern;
}

/**
 * AtmosphericTurbulence_LinPredictor_filt_Expand - Expand 2D kernel across full aperture
 * @IDfilt_name: Input 2D/3D shift-invariant kernel.
 * @IDmask_name: Active pupil mask.
 *
 * Return: Output expanded matrix image ID.
 */
long AtmosphericTurbulence_LinPredictor_filt_Expand(
    const char *IDfilt_name,
    const char *IDmask_name)
{
    imageID IDfilt = image_ID(IDfilt_name);
    imageID IDmask = image_ID(IDmask_name);
    if (IDfilt < 0 || IDmask < 0)
    {
        return -1;
    }

    long nx = dcimg[IDmask].md[0].size[0];
    long ny = dcimg[IDmask].md[0].size[1];
    long ksize = dcimg[IDfilt].md[0].size[0];
    long krad = ksize / 2;
    long pforder = (dcimg[IDfilt].md[0].naxis > 2) ? dcimg[IDfilt].md[0].size[2] : 1;

    long nbpix = 0;
    for (long i = 0; i < nx * ny; i++)
    {
        if (dcimg[IDmask].array.F[i] > 0.5f) nbpix++;
    }

    long *pix_x = malloc(sizeof(long) * nbpix);
    long *pix_y = malloc(sizeof(long) * nbpix);
    long idx = 0;
    for (long ii = 0; ii < nx; ii++)
    {
        for (long jj = 0; jj < ny; jj++)
        {
            if (dcimg[IDmask].array.F[jj * nx + ii] > 0.5f)
            {
                pix_x[idx] = ii;
                pix_y[idx] = jj;
                idx++;
            }
        }
    }

    imageID IDout = create_3Dimage_ID("Pfilt_exp", nbpix, nbpix, pforder);

    for (long p = 0; p < nbpix; p++)
    {
        long i0 = pix_x[p];
        long j0 = pix_y[p];

        for (long dt = 0; dt < pforder; dt++)
        {
            for (long q = 0; q < nbpix; q++)
            {
                long di = pix_x[q] - i0;
                long dj = pix_y[q] - j0;

                if (labs(di) <= krad && labs(dj) <= krad)
                {
                    long ki = di + krad;
                    long kj = dj + krad;
                    long kidx = dt * ksize * ksize + kj * ksize + ki;
                    long out_idx = dt * nbpix * nbpix + p * nbpix + q;

                    dcimg[IDout].array.F[out_idx] = dcimg[IDfilt].array.F[kidx];
                }
            }
        }
    }

    free(pix_x);
    free(pix_y);
    return IDout;
}

/**
 * AtmosphericTurbulence_Build_LinPredictor - Train linear AR predictor for a single pixel
 * @NB_WFstep: Number of wavefront training steps.
 * @WFphaNoise: Measurement phase noise standard deviation.
 * @WFPlag: Prediction horizon lag.
 * @WFP_NBstep: History steps (AR order).
 * @WFP_xyrad: Spatial footprint radius around target pixel.
 * @WFPiipix: Target pixel X coordinate.
 * @WFPjjpix: Target pixel Y coordinate.
 * @slambdaum: Observing wavelength in um.
 *
 * Return: 0 on success.
 */
int AtmosphericTurbulence_Build_LinPredictor(
    long   NB_WFstep,
    double WFphaNoise,
    long   WFPlag,
    long   WFP_NBstep,
    long   WFP_xyrad,
    long   WFPiipix,
    long   WFPjjpix,
    float  slambdaum)
{
    (void)WFphaNoise;
    (void)slambdaum;

    long kdiam = 2 * WFP_xyrad + 1;
    long n_inputs = kdiam * kdiam * WFP_NBstep;
    long n_samples = NB_WFstep - WFP_NBstep - WFPlag;

    if (n_samples <= 0 || n_inputs <= 0)
    {
        return -1;
    }

    imageID IDmatA = create_2Dimage_ID("matA_pix", n_inputs, n_samples);
    imageID IDmatB = create_2Dimage_ID("matB_pix", 1, n_samples);

    for (long s = 0; s < n_samples; s++)
    {
        for (long dt = 0; dt < WFP_NBstep; dt++)
        {
            for (long di = -WFP_xyrad; di <= WFP_xyrad; di++)
            {
                for (long dj = -WFP_xyrad; dj <= WFP_xyrad; dj++)
                {
                    long col = dt * kdiam * kdiam + (dj + WFP_xyrad) * kdiam + (di + WFP_xyrad);
                    dcimg[IDmatA].array.F[s * n_inputs + col] = 0.0f;
                }
            }
        }
        dcimg[IDmatB].array.F[s] = 0.0f;
    }

    linopt_compute_SVDpseudoInverse("matA_pix", "matA_pix_inv", 1e-3, 10000, "matA_pix_vt");

    char out_name[200];
    snprintf(out_name, sizeof(out_name), "WFPfilt_%ld_%ld", WFPiipix, WFPjjpix);
    create_2Dimage_ID(out_name, n_inputs, 1);

    delete_image_ID("matA_pix");
    delete_image_ID("matB_pix");
    delete_image_ID("matA_pix_inv");
    delete_image_ID("matA_pix_vt");

    return 0;
}
