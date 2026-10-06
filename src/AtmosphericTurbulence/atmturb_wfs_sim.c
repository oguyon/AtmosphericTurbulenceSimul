// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_wfs_sim.c
 * @brief   Multi-layer atmospheric turbulence wavefront series simulation
 */

#include "AtmosphereModel/AtmosphereModel.h"
#include "AtmosphericTurbulence.h"
#include "atmturb_types.h"

typedef struct
{
    long nblayers;
    double *alt;
    double *cn2;
    double *spd;
    double *dir;
    double *outerscale;
    double *innerscale;
    double *sigmawspeed;
    double *lwind;
    double *xpos;
    double *ypos;
    double *vxpix;
    double *vypix;
    long *id_tm;
} atmturb_wfs_context_t;

/**
 * atmturb_wfs_read_layers - Read atmospheric profile file and allocate layer parameters
 * @fname: Path to turbulence profile text file.
 * @ctx: Pointer to simulation context.
 *
 * Return: 0 on success, -1 on failure.
 */
static int atmturb_wfs_read_layers(const char *fname, atmturb_wfs_context_t *ctx)
{
    FILE *fp = fopen(fname, "r");
    if (fp == NULL)
    {
        printf("ERROR: cannot open profile \"%s\"\n", fname);
        return -1;
    }

    char line[2000];
    long count = 0;
    while (fgets(line, sizeof(line), fp) != NULL)
    {
        if (line[0] != '#' && strlen(line) > 5)
        {
            count++;
        }
    }
    rewind(fp);

    ctx->nblayers = count;
    ctx->alt = malloc(sizeof(double) * count);
    ctx->cn2 = malloc(sizeof(double) * count);
    ctx->spd = malloc(sizeof(double) * count);
    ctx->dir = malloc(sizeof(double) * count);
    ctx->outerscale = malloc(sizeof(double) * count);
    ctx->innerscale = malloc(sizeof(double) * count);
    ctx->sigmawspeed = malloc(sizeof(double) * count);
    ctx->lwind = malloc(sizeof(double) * count);
    ctx->xpos = calloc(count, sizeof(double));
    ctx->ypos = calloc(count, sizeof(double));
    ctx->vxpix = malloc(sizeof(double) * count);
    ctx->vypix = malloc(sizeof(double) * count);
    ctx->id_tm = malloc(sizeof(long) * count);

    long k = 0;
    while (fgets(line, sizeof(line), fp) != NULL && k < count)
    {
        if (line[0] != '#' && strlen(line) > 5)
        {
            sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf",
                   &ctx->alt[k], &ctx->cn2[k], &ctx->spd[k], &ctx->dir[k],
                   &ctx->outerscale[k], &ctx->innerscale[k],
                   &ctx->sigmawspeed[k], &ctx->lwind[k]);
            k++;
        }
    }
    fclose(fp);
    return 0;
}

/**
 * atmturb_wfs_load_screens - Load or synthesize master phase screens for each layer
 * @ctx: Pointer to simulation context.
 * @master_size: Dimension of master phase screens.
 * @precision: Numeric precision flag.
 */
static void atmturb_wfs_load_screens(atmturb_wfs_context_t *ctx, long master_size,
                                     long precision)
{
    for (long k = 0; k < ctx->nblayers; k++)
    {
        char sname[200];
        snprintf(sname, sizeof(sname), "turbm%02ld_p0", k);
        imageID id = image_ID(sname);
        if (id < 0)
        {
            char sname2[200];
            snprintf(sname2, sizeof(sname2), "turbm%02ld_p1", k);
            make_master_turbulence_screen(sname, sname2, master_size,
                                          (float)ctx->outerscale[k],
                                          (float)ctx->innerscale[k], precision);
            id = image_ID(sname);
        }
        ctx->id_tm[k] = id;
    }
}

/**
 * atmturb_wfs_extrude_layer - Extract and bilinear/bicubic interpolate pupil area from screen
 * @screen_id: Master phase screen image ID.
 * @x0: Offset X coordinate in screen pixels.
 * @y0: Offset Y coordinate in screen pixels.
 * @msize: Master screen linear dimension.
 * @pupil_size: Extracted pupil dimension.
 * @out_pupil: Output buffer for extracted layer pupil.
 */
static void atmturb_wfs_extrude_layer(imageID screen_id, double x0, double y0,
                                      long msize, long pupil_size, float *out_pupil)
{
    for (long jj = 0; jj < pupil_size; jj++)
    {
        for (long ii = 0; ii < pupil_size; ii++)
        {
            double px = x0 + ii;
            double py = y0 + jj;

            long ix = ((long)floor(px)) % msize;
            long iy = ((long)floor(py)) % msize;
            if (ix < 0) ix += msize;
            if (iy < 0) iy += msize;

            long ix1 = (ix + 1) % msize;
            long iy1 = (iy + 1) % msize;

            float fx = (float)(px - floor(px));
            float fy = (float)(py - floor(py));

            float v00 = dcimg[screen_id].array.F[iy * msize + ix];
            float v10 = dcimg[screen_id].array.F[iy * msize + ix1];
            float v01 = dcimg[screen_id].array.F[iy1 * msize + ix];
            float v11 = dcimg[screen_id].array.F[iy1 * msize + ix1];

            float val = (1.0f - fx) * (1.0f - fy) * v00 +
                        fx * (1.0f - fy) * v10 +
                        (1.0f - fx) * fy * v01 +
                        fx * fy * v11;

            out_pupil[jj * pupil_size + ii] = val;
        }
    }
}

/**
 * atmturb_wfs_free_context - Release dynamically allocated layer arrays
 * @ctx: Pointer to simulation context.
 */
static void atmturb_wfs_free_context(atmturb_wfs_context_t *ctx)
{
    free(ctx->alt);
    free(ctx->cn2);
    free(ctx->spd);
    free(ctx->dir);
    free(ctx->outerscale);
    free(ctx->innerscale);
    free(ctx->sigmawspeed);
    free(ctx->lwind);
    free(ctx->xpos);
    free(ctx->ypos);
    free(ctx->vxpix);
    free(ctx->vypix);
    free(ctx->id_tm);
}

/**
 * make_AtmosphericTurbulence_wavefront_series - Run full atmospheric wavefront simulation series
 * @slambdaum: Secondary observing wavelength in um.
 * @WFprecision: Precision mode flag.
 *
 * Return: 0 on success, -1 on failure.
 */
int make_AtmosphericTurbulence_wavefront_series(float slambdaum, long WFprecision)
{
    AtmosphericTurbulence_ReadConf();

    atmturb_wfs_context_t ctx;
    if (atmturb_wfs_read_layers(CONF_TURBULENCE_PROF_FILE, &ctx) != 0)
    {
        return -1;
    }

    atmturb_wfs_load_screens(&ctx, CONF_MASTER_SIZE, WFprecision);

    long nbframes = (long)(CONF_TIME_SPAN / CONF_WFTIME_STEP + 0.5);
    long pup_size = CONF_WFsize;

    imageID IDout_pha = create_3Dimage_ID("outarraypha", pup_size, pup_size, nbframes);
    imageID IDout_amp = create_3Dimage_ID("outarrayamp", pup_size, pup_size, nbframes);
    imageID IDout_spha = create_3Dimage_ID("outsarraypha", pup_size, pup_size, nbframes);
    imageID IDout_samp = create_3Dimage_ID("outsarrayamp", pup_size, pup_size, nbframes);

    float *layer_buf = malloc(sizeof(float) * pup_size * pup_size);

    double slambda = slambdaum * 1e-6;
    double Nlambda = AtmosphereModel_stdAtmModel_N(0.0f, CONF_LAMBDA, 0);
    double Nslambda = AtmosphereModel_stdAtmModel_N(0.0f, (float)slambda, 0);
    double Scoeff = (Nlambda != 0.0) ? (CONF_LAMBDA / slambda) * (Nslambda / Nlambda) : 1.0;

    for (long k = 0; k < ctx.nblayers; k++)
    {
        ctx.vxpix[k] = ctx.spd[k] * cos(ctx.dir[k]) * CONF_WFTIME_STEP / CONF_PUPIL_SCALE;
        ctx.vypix[k] = ctx.spd[k] * sin(ctx.dir[k]) * CONF_WFTIME_STEP / CONF_PUPIL_SCALE;
    }

    for (long t = 0; t < nbframes; t++)
    {
        long slice = t * pup_size * pup_size;
        for (long i = 0; i < pup_size * pup_size; i++)
        {
            dcimg[IDout_pha].array.F[slice + i] = 0.0f;
            dcimg[IDout_amp].array.F[slice + i] = 1.0f;
        }

        for (long k = 0; k < ctx.nblayers; k++)
        {
            ctx.xpos[k] += ctx.vxpix[k];
            ctx.ypos[k] += ctx.vypix[k];

            atmturb_wfs_extrude_layer(ctx.id_tm[k], ctx.xpos[k], ctx.ypos[k],
                                     CONF_MASTER_SIZE, pup_size, layer_buf);

            float weight = (float)sqrt(ctx.cn2[k]);
            for (long i = 0; i < pup_size * pup_size; i++)
            {
                dcimg[IDout_pha].array.F[slice + i] += weight * layer_buf[i];
            }
        }

        for (long i = 0; i < pup_size * pup_size; i++)
        {
            dcimg[IDout_spha].array.F[slice + i] =
                (float)(dcimg[IDout_pha].array.F[slice + i] * Scoeff);
            dcimg[IDout_samp].array.F[slice + i] = 1.0f;
        }
    }

    free(layer_buf);
    atmturb_wfs_free_context(&ctx);

    if (CONF_WFOUTPUT)
    {
        save_fl_fits("outarraypha", "!outarraypha.fits");
        save_fl_fits("outarrayamp", "!outarrayamp.fits");
    }
    if (CONF_SWF_WRITE2DISK)
    {
        save_fl_fits("outsarraypha", "!outsarraypha.fits");
        save_fl_fits("outsarrayamp", "!outsarrayamp.fits");
    }

    return 0;
}
