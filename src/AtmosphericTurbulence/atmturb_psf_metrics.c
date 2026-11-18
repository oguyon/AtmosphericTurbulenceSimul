// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_psf_metrics.c
 * @brief   Wavefront analysis, PSF aperture photometry, and frame selection metrics
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
#include "atmturb_fft_wisdom.h"
#include "atmturb_fresnel_simd.h"
#include "atmturb_screens_batch.h"
#include "atmturb_types.h"

typedef struct
{
    long           nx;
    long           ny;
    fftwf_complex *in_buf;
    fftwf_complex *out_buf;
    float         *amp_row;
    fftwf_plan     plan;
} atmturb_psf_metrics_fft_t;

/**
 * atmturb_psf_metrics_fft_init - Allocate persistent buffers and cached 2D FFTW plan
 */
static int atmturb_psf_metrics_fft_init(
    atmturb_psf_metrics_fft_t *ctx,
    long                       nx,
    long                       ny)
{
    memset(ctx, 0, sizeof(*ctx));
    ctx->nx      = nx;
    ctx->ny      = ny;
    ctx->in_buf  = fftwf_alloc_complex((size_t) nx * (size_t) ny);
    ctx->out_buf = fftwf_alloc_complex((size_t) nx * (size_t) ny);
    ctx->amp_row = (float *) malloc(sizeof(float) * (size_t) nx);
    if (ctx->in_buf == NULL || ctx->out_buf == NULL || ctx->amp_row == NULL)
    {
        fftwf_free(ctx->in_buf);
        fftwf_free(ctx->out_buf);
        free(ctx->amp_row);
        memset(ctx, 0, sizeof(*ctx));
        return -1;
    }

    atmturb_fftwf_plan_threads(0);
    atmturb_fftwf_wisdom_load();
    unsigned int flags = atmturb_fftwf_plan_flags();
    ctx->plan          = fftwf_plan_dft_2d((int) ny, (int) nx, ctx->in_buf, ctx->out_buf,
                                           FFTW_FORWARD, flags);
    if (ctx->plan == NULL)
    {
        fftwf_free(ctx->in_buf);
        fftwf_free(ctx->out_buf);
        free(ctx->amp_row);
        memset(ctx, 0, sizeof(*ctx));
        return -1;
    }
    if (flags == FFTW_MEASURE)
    {
        atmturb_fftwf_wisdom_save();
    }
    return 0;
}

/**
 * atmturb_psf_metrics_fft_free - Release persistent 2D FFTW plan and buffers
 */
static void atmturb_psf_metrics_fft_free(
    atmturb_psf_metrics_fft_t *ctx)
{
    if (ctx->plan != NULL)
    {
        fftwf_destroy_plan(ctx->plan);
    }
    fftwf_free(ctx->in_buf);
    fftwf_free(ctx->out_buf);
    free(ctx->amp_row);
    memset(ctx, 0, sizeof(*ctx));
}

/**
 * atmturb_psf_metrics_accum_span - Write or accumulate squared modulus over a contiguous span
 */
static inline void atmturb_psf_metrics_accum_span(
    const fftwf_complex *restrict src,
    float               *restrict dst,
    long                          len,
    int                           accum)
{
    if (accum)
    {
        #pragma omp simd
        for (long i = 0; i < len; i++)
        {
            float r  = src[i][0];
            float im = src[i][1];
            dst[i] += r * r + im * im;
        }
    }
    else
    {
        #pragma omp simd
        for (long i = 0; i < len; i++)
        {
            float r  = src[i][0];
            float im = src[i][1];
            dst[i]   = r * r + im * im;
        }
    }
}

/**
 * atmturb_compute_frame_psf - Form complex pupil field and calculate centered PSF via cached FFT
 */
static void atmturb_compute_frame_psf(
    const atmturb_psf_metrics_fft_t *ctx,
    imageID                          id_amp,
    imageID                          id_pha,
    long                             slice_idx,
    imageID                          id_pupamp,
    float                           *psf_dst,
    int                              accum)
{
    long         nx        = ctx->nx;
    long         ny        = ctx->ny;
    long         slice_off = slice_idx * nx * ny;
    const float *pha_base  = (id_pha >= 0) ? &dcimg[id_pha].array.F[slice_off] : NULL;
    const float *amp_base  = (id_amp >= 0) ? &dcimg[id_amp].array.F[slice_off] : NULL;
    const float *pup_base  = (id_pupamp >= 0) ? dcimg[id_pupamp].array.F : NULL;

    if (amp_base == NULL && pup_base == NULL)
    {
        #pragma omp simd
        for (long ii = 0; ii < nx; ii++)
        {
            ctx->amp_row[ii] = 1.0f;
        }
    }

    for (long jj = 0; jj < ny; jj++)
    {
        long         roff = jj * nx;
        const float *prow = (pha_base != NULL) ? &pha_base[roff] : NULL;
        const float *mask = ctx->amp_row;

        if (amp_base != NULL && pup_base != NULL)
        {
            const float *restrict arow = &amp_base[roff];
            const float *restrict mrow = &pup_base[roff];
            float       *restrict dst  = ctx->amp_row;
            #pragma omp simd
            for (long ii = 0; ii < nx; ii++)
            {
                dst[ii] = arow[ii] * mrow[ii];
            }
        }
        else if (amp_base != NULL)
        {
            mask = &amp_base[roff];
        }
        else if (pup_base != NULL)
        {
            mask = &pup_base[roff];
        }

        atmturb_fresnel_simd_phasor_masked(&ctx->in_buf[roff], prow, mask, 1.0f, nx);
    }

    fftwf_execute_dft(ctx->plan, ctx->in_buf, ctx->out_buf);

    long xhalf = nx / 2;
    long yhalf = ny / 2;
    for (long jj = 0; jj < ny; jj++)
    {
        long sy = (jj < yhalf) ? (jj + yhalf) : (jj - yhalf);
        if (sy >= ny)
        {
            sy = jj;
        }
        const fftwf_complex *srow = &ctx->out_buf[sy * nx];
        float               *drow = &psf_dst[jj * nx];
        atmturb_psf_metrics_accum_span(&srow[xhalf], drow, xhalf, accum);
        atmturb_psf_metrics_accum_span(srow, &drow[xhalf], xhalf, accum);
        for (long ii = 2 * xhalf; ii < nx; ii++)
        {
            float r    = srow[ii][0];
            float im   = srow[ii][1];
            float val  = r * r + im * im;
            drow[ii]   = accum ? (drow[ii] + val) : val;
        }
    }
}

/**
 * atmturb_measure_psf_apertures - Measure integrated flux in concentric circular apertures
 */
static void atmturb_measure_psf_apertures(
    const float *restrict psf,
    long                  nx,
    long                  ny,
    float                 focal_scale,
    double               *fluxes_out)
{
    double tot = 0.0;
    double f1  = 0.0;
    double f2  = 0.0;
    double f5  = 0.0;
    double f10 = 0.0;
    double f20 = 0.0;
    long   np  = nx * ny;

    #pragma omp simd reduction(+:tot)
    for (long i = 0; i < np; i++)
    {
        tot += (double) psf[i];
    }

    double cx    = 0.5 * (double) nx;
    double cy    = 0.5 * (double) ny;
    double fs    = (double) focal_scale;
    double rmax  = (fs > 0.0) ? (20.0 / fs) : (double) ((nx > ny) ? nx : ny);
    long   j_min = (long) floor(cy - rmax);
    long   j_max = (long) ceil(cy + rmax);
    long   i_min = (long) floor(cx - rmax);
    long   i_max = (long) ceil(cx + rmax);
    if (j_min < 0)
    {
        j_min = 0;
    }
    if (j_max >= ny)
    {
        j_max = ny - 1;
    }
    if (i_min < 0)
    {
        i_min = 0;
    }
    if (i_max >= nx)
    {
        i_max = nx - 1;
    }

    for (long jj = j_min; jj <= j_max; jj++)
    {
        double               dy   = ((double) jj - cy) * fs;
        double               dy2  = dy * dy;
        const float *restrict row = &psf[jj * nx];
        for (long ii = i_min; ii <= i_max; ii++)
        {
            double dx = ((double) ii - cx) * fs;
            double r2 = dx * dx + dy2;
            if (r2 < 400.0)
            {
                double val = (double) row[ii];
                f20 += val;
                f10 += (r2 < 100.0) ? val : 0.0;
                f5  += (r2 < 25.0)  ? val : 0.0;
                f2  += (r2 < 4.0)   ? val : 0.0;
                f1  += (r2 < 1.0)   ? val : 0.0;
            }
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
int measure_wavefront_series(
    float factor)
{
    (void) factor;
    AtmosphericTurbulence_ReadConf();

    const double SLAMBDA  = 1.65e-6;
    imageID      IDpupamp = image_ID("ST_pa");
    if (IDpupamp == -1)
    {
        printf("ERROR: pupil amplitude map not loaded\n");
        return -1;
    }

    long  nx          = dcimg[IDpupamp].md[0].size[0];
    long  ny          = dcimg[IDpupamp].md[0].size[1];
    float focal_scale = (float) (SLAMBDA / ((double) nx * (double) CONF_PUPIL_SCALE)
                                 / M_PI * 180.0 * 3600.0);
    long  nbframes    = (long) (CONF_TIME_SPAN / CONF_WFTIME_STEP);

    atmturb_psf_metrics_fft_t fft_ctx;
    if (atmturb_psf_metrics_fft_init(&fft_ctx, nx, ny) != 0)
    {
        return -1;
    }

    imageID IDpsf  = create_2Dimage_ID("PSF", nx, ny);
    FILE   *fpphot = fopen("phot.txt", "w");
    if (fpphot == NULL || IDpsf < 0)
    {
        if (fpphot != NULL)
        {
            fclose(fpphot);
        }
        delete_image_ID("PSF");
        atmturb_psf_metrics_fft_free(&fft_ctx);
        return -1;
    }

    for (long tspan = 0; tspan < CONF_NB_TSPAN; tspan++)
    {
        char fnamepha[512];
        char fnameamp[512];
        long lam_pm = (long) (1.0e12 * SLAMBDA + 0.5);
        snprintf(fnamepha, sizeof(fnamepha), "%s%08ld.%09ld.pha.fits",
                 CONF_WF_FILE_PREFIX, tspan, lam_pm);
        snprintf(fnameamp, sizeof(fnameamp), "%s%08ld.%09ld.amp.fits",
                 CONF_WF_FILE_PREFIX, tspan, lam_pm);

        imageID IDpha = load_fits(fnamepha, "wfpha", 1);
        imageID IDamp = (CONF_WAVEFRONT_AMPLITUDE == 1) ? load_fits(fnameamp, "wfamp", 1) : -1;

        for (long frame = 0; frame < nbframes; frame++)
        {
            atmturb_compute_frame_psf(&fft_ctx, IDamp, IDpha, frame, IDpupamp,
                                      dcimg[IDpsf].array.F, 0);
            double fluxes[6];
            atmturb_measure_psf_apertures(dcimg[IDpsf].array.F, nx, ny, focal_scale, fluxes);
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
    atmturb_psf_metrics_fft_free(&fft_ctx);
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

    const double SLAMBDA  = 1.65e-6;
    imageID      IDpupamp = image_ID("ST_pa");
    if (IDpupamp == -1)
    {
        return -1;
    }

    long  nx                = dcimg[IDpupamp].md[0].size[0];
    long  ny                = dcimg[IDpupamp].md[0].size[1];
    float focal_scale       = (float) (SLAMBDA / ((double) nx * (double) CONF_PUPIL_SCALE)
                                       / M_PI * 180.0 * 3600.0);
    long  nbframes_per_span = (long) (CONF_TIME_SPAN / CONF_WFTIME_STEP);
    long  nbframes_per_expo = (long) (etime / CONF_WFTIME_STEP);
    if (nbframes_per_expo < 1)
    {
        nbframes_per_expo = 1;
    }

    atmturb_psf_metrics_fft_t fft_ctx;
    if (atmturb_psf_metrics_fft_init(&fft_ctx, nx, ny) != 0)
    {
        return -1;
    }

    imageID IDpsf_acc = create_2Dimage_ID("PSF_acc", nx, ny);
    FILE   *fp        = fopen(outfile, "w");
    if (fp == NULL || IDpsf_acc < 0)
    {
        if (fp != NULL)
        {
            fclose(fp);
        }
        delete_image_ID("PSF_acc");
        atmturb_psf_metrics_fft_free(&fft_ctx);
        return -1;
    }

    long acc_count = 0;
    long expo_idx  = 0;

    for (long tspan = 0; tspan < CONF_NB_TSPAN; tspan++)
    {
        char fnamepha[512];
        char fnameamp[512];
        long lam_pm = (long) (1.0e12 * SLAMBDA + 0.5);
        snprintf(fnamepha, sizeof(fnamepha), "%s%08ld.%09ld.pha.fits",
                 CONF_WF_FILE_PREFIX, tspan, lam_pm);
        snprintf(fnameamp, sizeof(fnameamp), "%s%08ld.%09ld.amp.fits",
                 CONF_WF_FILE_PREFIX, tspan, lam_pm);

        imageID IDpha = load_fits(fnamepha, "wfpha", 1);
        imageID IDamp = (CONF_WAVEFRONT_AMPLITUDE == 1) ? load_fits(fnameamp, "wfamp", 1) : -1;

        for (long frame = 0; frame < nbframes_per_span; frame++)
        {
            atmturb_compute_frame_psf(&fft_ctx, IDamp, IDpha, frame, IDpupamp,
                                      dcimg[IDpsf_acc].array.F, 1);
            acc_count++;

            if (acc_count >= nbframes_per_expo)
            {
                double fluxes[6];
                atmturb_measure_psf_apertures(dcimg[IDpsf_acc].array.F, nx, ny,
                                              focal_scale, fluxes);
                fprintf(fp, "%ld %g %g %g %g %g %g\n", expo_idx++,
                        fluxes[0], fluxes[1], fluxes[2], fluxes[3], fluxes[4], fluxes[5]);
                memset(dcimg[IDpsf_acc].array.F, 0, sizeof(float) * (size_t) nx * (size_t) ny);
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
    delete_image_ID("PSF_acc");
    atmturb_psf_metrics_fft_free(&fft_ctx);
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

    long                  nx  = dcimg[IDcube].md[0].size[0];
    long                  ny  = dcimg[IDcube].md[0].size[1];
    long                  nz  =
        (dcimg[IDcube].md[0].naxis > 2) ? dcimg[IDcube].md[0].size[2] : 1;
    long                  nxy = nx * ny;
    const float *restrict cube_data = dcimg[IDcube].array.F;
    const float *restrict mask_data = dcimg[IDmask].array.F;
    double                mask_norm = 0.0;

    #pragma omp simd reduction(+:mask_norm)
    for (long i = 0; i < nxy; i++)
    {
        mask_norm += (double) mask_data[i];
    }
    if (mask_norm < 1.0)
    {
        mask_norm = 1.0;
    }

    imageID         IDout     = create_2Dimage_ID(IDpsfc_name, nz, 1);
    float *restrict out_data  = dcimg[IDout].array.F;
    double          inv_norm  = 1.0 / mask_norm;

    #pragma omp parallel for schedule(static) if (nz >= 4)
    for (long k = 0; k < nz; k++)
    {
        const float *restrict slice = &cube_data[k * nxy];
        double                sum   = 0.0;
        #pragma omp simd reduction(+:sum)
        for (long i = 0; i < nxy; i++)
        {
            sum += (double) slice[i] * (double) mask_data[i];
        }
        out_data[k] = (float) (sum * inv_norm);
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

    long n_select = (long) ((float) NBfiles * frac);
    if (n_select < 1)
    {
        n_select = 1;
    }

    float *metrics = (float *) calloc((size_t) NBfiles, sizeof(float));
    long  *indices = (long *) calloc((size_t) NBfiles, sizeof(long));
    if (metrics == NULL || indices == NULL)
    {
        free(metrics);
        free(indices);
        fclose(fp);
        return -1;
    }

    for (long i = 0; i < NBfiles; i++)
    {
        long  idx = 0;
        float val = 0.0f;
        if (fscanf(fp, "%ld %*f %f %*[^\n]", &idx, &val) >= 2)
        {
            metrics[i] = val;
            indices[i] = idx;
        }
    }
    fclose(fp);

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
        float tmp_m     = metrics[i];
        metrics[i]      = metrics[best_j];
        metrics[best_j] = tmp_m;
        long tmp_idx    = indices[i];
        indices[i]      = indices[best_j];
        indices[best_j] = tmp_idx;
    }

    printf("Selected %ld frames out of %ld (threshold: %g)\n",
           n_select, NBfiles, metrics[n_select - 1]);

    free(metrics);
    free(indices);
    return 0;
}
