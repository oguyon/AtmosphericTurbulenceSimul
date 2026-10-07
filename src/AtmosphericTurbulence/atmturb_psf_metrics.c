// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_psf_metrics.c
 * @brief   Wavefront analysis, PSF aperture photometry, and frame selection metrics
 */

#include "AtmosphericTurbulence.h"
#include "atmturb_types.h"

/**
 * atmturb_compute_frame_psf - Form complex pupil field and calculate PSF via FFT
 * @id_amp: Amplitude image ID (-1 if unity amplitude).
 * @id_pha: Phase image ID.
 * @slice_idx: Frame index along cube Z-axis.
 * @id_pupamp: Pupil mask image ID (-1 if none).
 * @nx: Pupil linear width.
 * @ny: Pupil linear height.
 * @id_psf: Output PSF image ID.
 */
static void atmturb_compute_frame_psf(
    imageID id_amp,
    imageID id_pha,
    long    slice_idx,
    imageID id_pupamp,
    long    nx,
    long    ny,
    imageID id_psf)
{
    imageID id_arr = create_2DCimage_ID("tmp_carr", nx, ny);

    for (long ii = 0; ii < nx; ii++)
    {
        for (long jj = 0; jj < ny; jj++)
        {
            long pidx = slice_idx * nx * ny + jj * nx + ii;
            long pix2d = jj * nx + ii;
            float p = (id_pha >= 0) ? dcimg[id_pha].array.F[pidx] : 0.0f;
            float a = (id_amp >= 0) ? dcimg[id_amp].array.F[pidx] : 1.0f;
            if (id_pupamp >= 0)
            {
                a *= dcimg[id_pupamp].array.F[pix2d];
            }
            dcimg[id_arr].array.CF[pix2d].re = a * cosf(p);
            dcimg[id_arr].array.CF[pix2d].im = a * sinf(p);
        }
    }

    permut("tmp_carr");
    do2dfft("tmp_carr", "tmp_cfft");
    delete_image_ID("tmp_carr");

    imageID id_cfft = image_ID("tmp_cfft");
    for (long i = 0; i < nx * ny; i++)
    {
        float r = dcimg[id_cfft].array.CF[i].re;
        float im = dcimg[id_cfft].array.CF[i].im;
        dcimg[id_psf].array.F[i] = r * r + im * im;
    }
    delete_image_ID("tmp_cfft");
}

/**
 * atmturb_measure_psf_apertures - Measure integrated flux in concentric circular apertures
 * @id_psf: PSF image ID.
 * @nx: Image width.
 * @ny: Image height.
 * @focal_scale: Plate scale in arcseconds/pixel.
 * @fluxes_out: Output array of 6 fluxes (total, 1", 2", 5", 10", 20").
 */
static void atmturb_measure_psf_apertures(
    imageID  id_psf,
    long     nx,
    long     ny,
    float    focal_scale,
    double  *fluxes_out)
{
    double tot = 0.0, f1 = 0.0, f2 = 0.0, f5 = 0.0, f10 = 0.0, f20 = 0.0;

    for (long ii = 0; ii < nx; ii++)
    {
        for (long jj = 0; jj < ny; jj++)
        {
            double dx = (1.0 * ii - nx / 2) * focal_scale;
            double dy = (1.0 * jj - ny / 2) * focal_scale;
            double r = sqrt(dx * dx + dy * dy);
            float val = dcimg[id_psf].array.F[jj * nx + ii];

            tot += val;
            if (r < 1.0) f1 += val;
            if (r < 2.0) f2 += val;
            if (r < 5.0) f5 += val;
            if (r < 10.0) f10 += val;
            if (r < 20.0) f20 += val;
        }
    }

    fluxes_out[0] = tot;
    fluxes_out[1] = f1;
    fluxes_out[2] = f2;
    fluxes_out[3] = f5;
    fluxes_out[4] = f10;
    fluxes_out[5] = f20;
}

/**
 * measure_wavefront_series - Measure encircled energy and PSF metrics across simulation cubes
 * @factor: Scale factor.
 *
 * Return: 0 on success.
 */
int measure_wavefront_series(float factor)
{
    (void)factor;
    AtmosphericTurbulence_ReadConf();

    float focal_scale = (float)(CONF_LAMBDA / CONF_WFsize / CONF_PUPIL_SCALE / PI * 180.0 * 3600.0);
    imageID IDpupamp = image_ID("ST_pa");
    if (IDpupamp == -1)
    {
        printf("ERROR: pupil amplitude map not loaded\n");
        return -1;
    }

    long nx = dcimg[IDpupamp].md[0].size[0];
    long ny = dcimg[IDpupamp].md[0].size[1];
    long nbframes = (long)(CONF_TIME_SPAN / CONF_WFTIME_STEP);
    imageID IDpsf = create_2Dimage_ID("PSF", nx, ny);

    FILE *fpphot = fopen("phot.txt", "w");
    if (fpphot == NULL)
    {
        delete_image_ID("PSF");
        return -1;
    }

    const double SLAMBDA = 1.65e-6;
    for (long tspan = 0; tspan < CONF_NB_TSPAN; tspan++)
    {
        char fnamepha[512], fnameamp[512];
        snprintf(fnamepha, sizeof(fnamepha), "%s%08ld.%09ld.pha.fits",
                 CONF_WF_FILE_PREFIX, tspan, (long)(1.0e12 * SLAMBDA + 0.5));
        snprintf(fnameamp, sizeof(fnameamp), "%s%08ld.%09ld.amp.fits",
                 CONF_WF_FILE_PREFIX, tspan, (long)(1.0e12 * SLAMBDA + 0.5));

        imageID IDpha = load_fits(fnamepha, "wfpha", 1);
        imageID IDamp = (CONF_WAVEFRONT_AMPLITUDE == 1) ? load_fits(fnameamp, "wfamp", 1) : -1;

        for (long frame = 0; frame < nbframes; frame++)
        {
            atmturb_compute_frame_psf(IDamp, IDpha, frame, IDpupamp, nx, ny, IDpsf);
            double fluxes[6];
            atmturb_measure_psf_apertures(IDpsf, nx, ny, focal_scale, fluxes);
            fprintf(fpphot, "%ld %ld %g %g %g %g %g %g\n",
                    tspan, frame, fluxes[0], fluxes[1], fluxes[2], fluxes[3], fluxes[4], fluxes[5]);
        }

        delete_image_ID("wfpha");
        if (IDamp >= 0)
        {
            delete_image_ID("wfamp");
        }
    }

    fclose(fpphot);
    delete_image_ID("PSF");
    return 0;
}

/**
 * measure_wavefront_series_expoframes - Integrate frames into specified exposure duration
 * @etime: Exposure duration in seconds.
 * @outfile: Destination output text file.
 *
 * Return: 0 on success.
 */
int measure_wavefront_series_expoframes(
    float       etime,
    const char *outfile)
{
    AtmosphericTurbulence_ReadConf();

    float focal_scale = (float)(CONF_LAMBDA / CONF_WFsize / CONF_PUPIL_SCALE
                                / PI * 180.0 * 3600.0);
    imageID IDpupamp = image_ID("ST_pa");
    if (IDpupamp == -1)
    {
        return -1;
    }

    long nx = dcimg[IDpupamp].md[0].size[0];
    long ny = dcimg[IDpupamp].md[0].size[1];
    long nbframes_per_span = (long)(CONF_TIME_SPAN / CONF_WFTIME_STEP);
    long nbframes_per_expo = (long)(etime / CONF_WFTIME_STEP);
    if (nbframes_per_expo < 1) nbframes_per_expo = 1;

    imageID IDpsf = create_2Dimage_ID("PSF", nx, ny);
    imageID IDpsf_acc = create_2Dimage_ID("PSF_acc", nx, ny);

    FILE *fp = fopen(outfile, "w");
    if (fp == NULL)
    {
        delete_image_ID("PSF");
        delete_image_ID("PSF_acc");
        return -1;
    }

    const double SLAMBDA = 1.65e-6;
    long acc_count = 0;
    long expo_idx = 0;

    for (long tspan = 0; tspan < CONF_NB_TSPAN; tspan++)
    {
        char fnamepha[512], fnameamp[512];
        snprintf(fnamepha, sizeof(fnamepha), "%s%08ld.%09ld.pha.fits",
                 CONF_WF_FILE_PREFIX, tspan, (long)(1.0e12 * SLAMBDA + 0.5));
        snprintf(fnameamp, sizeof(fnameamp), "%s%08ld.%09ld.amp.fits",
                 CONF_WF_FILE_PREFIX, tspan, (long)(1.0e12 * SLAMBDA + 0.5));

        imageID IDpha = load_fits(fnamepha, "wfpha", 1);
        imageID IDamp = (CONF_WAVEFRONT_AMPLITUDE == 1) ? load_fits(fnameamp, "wfamp", 1) : -1;

        for (long frame = 0; frame < nbframes_per_span; frame++)
        {
            atmturb_compute_frame_psf(IDamp, IDpha, frame, IDpupamp, nx, ny, IDpsf);
            for (long p = 0; p < nx * ny; p++)
            {
                dcimg[IDpsf_acc].array.F[p] += dcimg[IDpsf].array.F[p];
            }
            acc_count++;

            if (acc_count >= nbframes_per_expo)
            {
                double fluxes[6];
                atmturb_measure_psf_apertures(IDpsf_acc, nx, ny, focal_scale, fluxes);
                fprintf(fp, "%ld %g %g %g %g %g %g\n", expo_idx++,
                        fluxes[0], fluxes[1], fluxes[2], fluxes[3], fluxes[4], fluxes[5]);
                memset(dcimg[IDpsf_acc].array.raw, 0, sizeof(float) * nx * ny);
                acc_count = 0;
            }
        }

        delete_image_ID("wfpha");
        if (IDamp >= 0)
        {
            delete_image_ID("wfamp");
        }
    }

    fclose(fp);
    delete_image_ID("PSF");
    delete_image_ID("PSF_acc");
    return 0;
}

/**
 * AtmosphericTurbulence_psfCubeContrast - Compute contrast statistics inside mask across PSF cube
 * @IDwfc_name: Wavefront / PSF cube name.
 * @IDmask_name: Mask image name (1 inside dark hole, 0 elsewhere).
 * @IDpsfc_name: Output contrast profile image name.
 *
 * Return: Number of processed slices.
 */
long AtmosphericTurbulence_psfCubeContrast(
    const char *IDwfc_name,
    const char *IDmask_name,
    const char *IDpsfc_name)
{
    imageID IDcube = image_ID(IDwfc_name);
    imageID IDmask = image_ID(IDmask_name);
    if (IDcube < 0 || IDmask < 0)
    {
        return -1;
    }

    long nx = dcimg[IDcube].md[0].size[0];
    long ny = dcimg[IDcube].md[0].size[1];
    long nz = (dcimg[IDcube].md[0].naxis > 2) ? dcimg[IDcube].md[0].size[2] : 1;

    double mask_norm = 0.0;
    for (long i = 0; i < nx * ny; i++)
    {
        mask_norm += dcimg[IDmask].array.F[i];
    }
    if (mask_norm < 1.0)
    {
        mask_norm = 1.0;
    }

    imageID IDout = create_2Dimage_ID(IDpsfc_name, nz, 1);

    for (long k = 0; k < nz; k++)
    {
        double sum = 0.0;
        for (long i = 0; i < nx * ny; i++)
        {
            sum += dcimg[IDcube].array.F[k * nx * ny + i] * dcimg[IDmask].array.F[i];
        }
        dcimg[IDout].array.F[k] = (float)(sum / mask_norm);
    }

    return nz;
}

/**
 * frame_select_PSF - Lucky imaging frame selection from logged PSF series
 * @logfile: Text log containing frame index and metric columns.
 * @NBfiles: Number of log rows / frames.
 * @frac: Selection fraction (e.g. 0.1 for top 10%).
 *
 * Return: 0 on success.
 */
int frame_select_PSF(
    const char *logfile,
    long        NBfiles,
    float       frac)
{
    FILE *fp = fopen(logfile, "r");
    if (fp == NULL)
    {
        return -1;
    }

    long n_select = (long)(NBfiles * frac);
    if (n_select < 1)
    {
        n_select = 1;
    }

    float *metrics = malloc(sizeof(float) * NBfiles);
    long *indices = malloc(sizeof(long) * NBfiles);

    for (long i = 0; i < NBfiles; i++)
    {
        long idx = 0;
        float val = 0.0f;
        if (fscanf(fp, "%ld %*f %f %*[^\n]", &idx, &val) >= 2)
        {
            metrics[i] = val;
            indices[i] = idx;
        }
    }
    fclose(fp);

    // Simple selection sort for top n_select frames
    for (long i = 0; i < n_select; i++)
    {
        long best_j = i;
        for (long j = i + 1; j < NBfiles; j++)
        {
            if (metrics[j] > metrics[best_j])
            {
                best_j = j;
            }
        }
        float tmp_m = metrics[i];
        metrics[i] = metrics[best_j];
        metrics[best_j] = tmp_m;
        long tmp_idx = indices[i];
        indices[i] = indices[best_j];
        indices[best_j] = tmp_idx;
    }

    printf("Selected %ld frames out of %ld (threshold: %g)\n",
           n_select, NBfiles, metrics[n_select - 1]);

    free(metrics);
    free(indices);
    return 0;
}
