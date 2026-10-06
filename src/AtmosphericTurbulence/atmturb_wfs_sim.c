// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_wfs_sim.c
 * @brief   Multi-layer atmospheric turbulence wavefront series simulation
 */

#include <omp.h>
#include "AtmosphereModel/AtmosphereModel.h"
#include "AtmosphericTurbulence.h"
#include "atmturb_types.h"
#include "atmturb_simd.h"
#if defined(HAVE_CUDA)
#    include "atmturb_cuda.h"
#endif

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
 * atmturb_wfs_load_default_layers - Load built-in 7-layer atmospheric turbulence profile
 * @ctx: Pointer to simulation context.
 *
 * Return: 0 on success.
 */
static int atmturb_wfs_load_default_layers(atmturb_wfs_context_t *ctx)
{
    static const double def_alt[7] = {4215.0, 4230.0, 4349.0, 5007.0, 12000.0, 16200.0, 23701.0};
    static const double def_cn2[7] = {5.32,   1.47,   1.08,   2.11,   1.83,    1.48,    0.697};
    static const double def_spd[7] = {6.5,    6.55,   6.6,    6.7,   22.0,     9.5,     5.6};
    static const double def_dir[7] = {1.47,   1.57,   1.67,   1.77,   3.10,    3.20,    3.30};

    long count = 7;
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

    for (long k = 0; k < count; k++)
    {
        ctx->alt[k] = def_alt[k];
        ctx->cn2[k] = def_cn2[k];
        ctx->spd[k] = def_spd[k];
        ctx->dir[k] = def_dir[k];
        ctx->outerscale[k] = 50.0;
        ctx->innerscale[k] = 0.01;
        ctx->sigmawspeed[k] = 0.0;
        ctx->lwind[k] = 500.0;
    }
    return 0;
}

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
        if (strcmp(fname, "turbul.prof") == 0)
        {
            printf("[milkatmturb] Notice: Profile \"%s\" not found, using built-in 7-layer profile.\n",
                   fname);
            return atmturb_wfs_load_default_layers(ctx);
        }
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
            ctx->outerscale[k] = 50.0;
            ctx->innerscale[k] = 0.01;
            ctx->sigmawspeed[k] = 0.0;
            ctx->lwind[k] = 500.0;
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
    float pscale = (CONF_PUPIL_SCALE > 1e-6f) ? CONF_PUPIL_SCALE : 0.18f;
    for (long k = 0; k < ctx->nblayers; k++)
    {
        char sname[200];
        snprintf(sname, sizeof(sname), "turbm%02ld_p0", k);
        imageID id = image_ID(sname);
        if (id < 0)
        {
            char sname2[200];
            snprintf(sname2, sizeof(sname2), "turbm%02ld_p1", k);
            float osc = (float)(ctx->outerscale[k] / pscale);
            float isc = (float)(ctx->innerscale[k] / pscale);
            if (osc < 1.0f)
            {
                osc = 100.0f;
            }
            if (isc < 0.1f)
            {
                isc = 1.0f;
            }
            make_master_turbulence_screen(sname, sname2, master_size,
                                          osc, isc, precision);
            id = image_ID(sname);
        }
        ctx->id_tm[k] = id;
    }
}

/**
 * atmturb_wfs_render_frames - Multi-threaded SIMD rendering of simulation time steps
 * @ctx: Pointer to simulation context.
 * @pup_size: Linear dimension of the pupil grid.
 * @nbframes: Number of frames to synthesize.
 * @Scoeff: Chromatic dispersion scaling coefficient for secondary wavelength.
 * @IDout_pha: Primary phase 3D image ID.
 * @IDout_amp: Primary amplitude 3D image ID.
 * @IDout_spha: Secondary phase 3D image ID.
 * @IDout_samp: Secondary amplitude 3D image ID.
 */
static void atmturb_wfs_render_frames(const atmturb_wfs_context_t *ctx, long pup_size,
                                      long nbframes, double Scoeff, imageID IDout_pha,
                                      imageID IDout_amp, imageID IDout_spha,
                                      imageID IDout_samp)
{
    long frame_pixels = pup_size * pup_size;

    #pragma omp parallel for schedule(dynamic)
    for (long t = 0; t < nbframes; t++)
    {
        long slice = t * frame_pixels;
        float *pha_slice  = &dcimg[IDout_pha].array.F[slice];
        float *amp_slice  = &dcimg[IDout_amp].array.F[slice];
        float *spha_slice = &dcimg[IDout_spha].array.F[slice];
        float *samp_slice = &dcimg[IDout_samp].array.F[slice];

        atmturb_init_phase_amp(pha_slice, amp_slice, frame_pixels);

        for (long k = 0; k < ctx->nblayers; k++)
        {
            double cur_x = (double)(t + 1) * ctx->vxpix[k];
            double cur_y = (double)(t + 1) * ctx->vypix[k];
            float weight = (float)sqrt(ctx->cn2[k]);

            atmturb_extrude_accumulate(dcimg[ctx->id_tm[k]].array.F, CONF_MASTER_SIZE,
                                       cur_x, cur_y, pup_size, weight, pha_slice);
        }

        atmturb_scale_float_array(spha_slice, pha_slice, (float)Scoeff, frame_pixels);
        atmturb_scale_float_array(samp_slice, amp_slice, 1.0f, frame_pixels);
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
 * atmturb_wfs_dispatch_render - Dispatch rendering to CUDA GPU or multi-threaded CPU SIMD
 * @ctx: Pointer to simulation context.
 * @pup_size: Linear dimension of the pupil grid.
 * @nbframes: Number of frames to synthesize.
 * @Scoeff: Chromatic dispersion scaling coefficient for secondary wavelength.
 * @IDout_pha: Primary phase 3D image ID.
 * @IDout_amp: Primary amplitude 3D image ID.
 * @IDout_spha: Secondary phase 3D image ID.
 * @IDout_samp: Secondary amplitude 3D image ID.
 */
static void atmturb_wfs_dispatch_render(const atmturb_wfs_context_t *ctx, long pup_size,
                                        long nbframes, double Scoeff, imageID IDout_pha,
                                        imageID IDout_amp, imageID IDout_spha,
                                        imageID IDout_samp)
{
#if defined(HAVE_CUDA)
    if (atmturb_cuda_device_available())
    {
        const float **h_masters = (const float **)malloc(sizeof(const float *) * ctx->nblayers);
        for (long k = 0; k < ctx->nblayers; k++)
        {
            h_masters[k] = dcimg[ctx->id_tm[k]].array.F;
        }

        atmturb_cuda_sim_params_t params = {
            .nblayers = ctx->nblayers,
            .msize = CONF_MASTER_SIZE,
            .pup_size = pup_size,
            .nbframes = nbframes,
            .Scoeff = Scoeff,
            .h_masters = h_masters,
            .vxpix = ctx->vxpix,
            .vypix = ctx->vypix,
            .cn2 = ctx->cn2
        };

        atmturb_cuda_sim_outputs_t outputs = {
            .pha = dcimg[IDout_pha].array.F,
            .amp = dcimg[IDout_amp].array.F,
            .spha = dcimg[IDout_spha].array.F,
            .samp = dcimg[IDout_samp].array.F
        };

        printf("Synthesizing %ld wavefront frames [CUDA GPU]\n", nbframes);
        int res = atmturb_wfs_render_frames_cuda(&params, &outputs);
        free(h_masters);
        if (res == 0)
        {
            return;
        }
    }
#endif

    printf("Synthesizing %ld wavefront frames [%s SIMD]\n", nbframes,
           atmturb_simd_active_isa());

    atmturb_wfs_render_frames(ctx, pup_size, nbframes, Scoeff,
                              IDout_pha, IDout_amp, IDout_spha, IDout_samp);
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
    if (CONFFILE[0] != '\0')
    {
        if (AtmosphericTurbulence_ReadConf() != 0)
        {
            return -1;
        }
    }

    atmturb_wfs_context_t ctx;
    if (atmturb_wfs_read_layers(CONF_TURBULENCE_PROF_FILE, &ctx) != 0)
    {
        return -1;
    }

    atmturb_wfs_load_screens(&ctx, CONF_MASTER_SIZE, WFprecision);

    long nbframes = (long)(CONF_TIME_SPAN / CONF_WFTIME_STEP + 0.5);
    long pup_size = CONF_WFsize;

    const char *pha_name = (CONF_WF_PHASE_NAME[0] != '\0') ? CONF_WF_PHASE_NAME : "outarraypha";
    const char *amp_name = (CONF_WF_AMPL_NAME[0] != '\0') ? CONF_WF_AMPL_NAME : "outarrayamp";

    imageID IDout_pha = create_3Dimage_ID(pha_name, pup_size, pup_size, nbframes);
    imageID IDout_amp = create_3Dimage_ID(amp_name, pup_size, pup_size, nbframes);
    imageID IDout_spha = create_3Dimage_ID("outsarraypha", pup_size, pup_size, nbframes);
    imageID IDout_samp = create_3Dimage_ID("outsarrayamp", pup_size, pup_size, nbframes);

    double slambda = slambdaum * 1e-6;
    double Nlambda = AtmosphereModel_stdAtmModel_N(0.0f, CONF_LAMBDA, 0);
    double Nslambda = AtmosphereModel_stdAtmModel_N(0.0f, (float)slambda, 0);
    double Scoeff = (Nlambda != 0.0) ? (CONF_LAMBDA / slambda) * (Nslambda / Nlambda) : 1.0;

    for (long k = 0; k < ctx.nblayers; k++)
    {
        ctx.vxpix[k] = ctx.spd[k] * cos(ctx.dir[k]) * CONF_WFTIME_STEP / CONF_PUPIL_SCALE;
        ctx.vypix[k] = ctx.spd[k] * sin(ctx.dir[k]) * CONF_WFTIME_STEP / CONF_PUPIL_SCALE;
    }

    atmturb_wfs_dispatch_render(&ctx, pup_size, nbframes, Scoeff,
                                IDout_pha, IDout_amp, IDout_spha, IDout_samp);

    atmturb_wfs_free_context(&ctx);

    if (CONF_WFOUTPUT)
    {
        char fname_pha[200], fname_amp[200];
        snprintf(fname_pha, sizeof(fname_pha), "%s.fits", pha_name);
        snprintf(fname_amp, sizeof(fname_amp), "%s.fits", amp_name);
        save_fl_fits(pha_name, fname_pha);
        save_fl_fits(amp_name, fname_amp);
    }
    if (CONF_SWF_WRITE2DISK)
    {
        save_fl_fits("outsarraypha", "outsarraypha.fits");
        save_fl_fits("outsarrayamp", "outsarrayamp.fits");
    }

    return 0;
}
