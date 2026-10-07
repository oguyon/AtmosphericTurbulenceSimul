// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_wfs_sim.c
 * @brief   Multi-layer atmospheric turbulence wavefront series simulation
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "AtmosphereModel/AtmosphereModel.h"
#include "AtmosphericTurbulence.h"
#include "atmturb_geometry.h"
#include "atmturb_profile.h"
#include "atmturb_simd.h"
#include "atmturb_types.h"

/**
 * struct atmturb_wfs_images_t - Output 3D image handles container
 * @ID_pha: Primary phase 3D image handle.
 * @ID_amp: Primary amplitude 3D image handle.
 * @ID_spha: Secondary phase 3D image handle.
 * @ID_samp: Secondary amplitude 3D image handle.
 */
typedef struct
{
    imageID ID_pha;
    imageID ID_amp;
    imageID ID_spha;
    imageID ID_samp;
} atmturb_wfs_images_t;

/**
 * atmturb_wfs_ensure_float_screen - Validate a master screen and convert it to FP32 if needed
 * @name: Master screen image name.
 * @msize: Expected linear dimension in pixels.
 *
 * Return: Image ID of the FP32 screen, or -1 on size/type mismatch or allocation failure.
 */
static imageID atmturb_wfs_ensure_float_screen(
    const char *name,
    long        msize)
{
    imageID id = image_ID(name);
    if (id < 0)
    {
        return -1;
    }
    if ((long) dcimg[id].md[0].size[0] != msize || (long) dcimg[id].md[0].size[1] != msize)
    {
        printf("ERROR: master screen \"%s\" is %ld x %ld, expected %ld x %ld\n", name,
               (long) dcimg[id].md[0].size[0], (long) dcimg[id].md[0].size[1], msize, msize);
        return -1;
    }
    if (dcimg[id].md[0].datatype == _DATATYPE_FLOAT)
    {
        return id;
    }
    if (dcimg[id].md[0].datatype != _DATATYPE_DOUBLE)
    {
        printf("ERROR: master screen \"%s\" must be FLOAT or DOUBLE\n", name);
        return -1;
    }

    long ntot = msize * msize;
    float *tmp = (float *) malloc(sizeof(float) * (size_t) ntot);
    if (tmp == NULL)
    {
        return -1;
    }
    for (long ii = 0; ii < ntot; ii++)
    {
        tmp[ii] = (float) dcimg[id].array.D[ii];
    }
    delete_image_ID(name);
    id = create_2Dimage_ID(name, msize, msize);
    memcpy(dcimg[id].array.F, tmp, sizeof(float) * (size_t) ntot);
    free(tmp);
    return id;
}

/**
 * atmturb_wfs_load_screens - Load or synthesize FP32 master phase screens for each layer
 * @prof: Active turbulence profile.
 * @geom: Computed observing geometry.
 * @master_size: Master screen dimension in pixels.
 * @precision: FFT precision flag (0=single, 1=double).
 * @seed: Master PRNG seed.
 * @id_tm: Output array of image IDs per layer.
 *
 * Return: 0 on success, -1 on failure.
 */
static int atmturb_wfs_load_screens(
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    long                     master_size,
    long                     precision,
    uint64_t                 seed,
    long                    *id_tm)
{
    for (int k = 0; k < prof->nlayers; k++)
    {
        char sname[200];
        snprintf(sname, sizeof(sname), "turbm%02d_p0", k);
        if (!CONF_SKIP_EXISTING || image_ID(sname) < 0)
        {
            char sname2[200];
            snprintf(sname2, sizeof(sname2), "turbm%02d_p1", k);

            atmturb_screen_spec_t spec;
            memset(&spec, 0, sizeof(spec));
            spec.size = master_size;
            spec.r0_pix = geom->r0_ref_pix;
            spec.L0_pix = (prof->layers[k].L0_m > 0.0)
                              ? (prof->layers[k].L0_m / geom->dx_master_m) : 0.0;
            spec.l0_pix = (prof->layers[k].l0_m > 0.0)
                              ? (prof->layers[k].l0_m / geom->dx_master_m) : 0.0;
            spec.seed = atmturb_rng_stream_seed(seed, (uint64_t) k);
            spec.precision = (int) precision;

            delete_image_ID(sname);
            delete_image_ID(sname2);
            imageID id0 = create_2Dimage_ID(sname, master_size, master_size);
            imageID id1 = create_2Dimage_ID(sname2, master_size, master_size);
            if (id0 < 0 || id1 < 0)
            {
                return -1;
            }

            int ret = atmturb_generate_screen_pair(&spec, dcimg[id0].array.F, dcimg[id1].array.F);
            if (ret != 0)
            {
                return -1;
            }
        }
        id_tm[k] = atmturb_wfs_ensure_float_screen(sname, master_size);
        if (id_tm[k] < 0)
        {
            return -1;
        }
    }
    return 0;
}

/**
 * atmturb_wfs_render_frames - Multi-threaded SIMD rendering of simulation time steps
 * @geom: Computed observing geometry.
 * @id_tm: Array of master screen image IDs.
 * @master_size: Master screen dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 * @nbframes: Number of simulation frames.
 * @imgs: Container of output 3D image handles.
 */
static void atmturb_wfs_render_frames(
    const atmturb_geom_t       *geom,
    const long                 *id_tm,
    long                        master_size,
    long                        pup_size,
    long                        nbframes,
    const atmturb_wfs_images_t *imgs)
{
    long frame_pixels = pup_size * pup_size;

    #pragma omp parallel for schedule(dynamic)
    for (long t = 0; t < nbframes; t++)
    {
        long slice = t * frame_pixels;
        float *pha_slice  = &dcimg[imgs->ID_pha].array.F[slice];
        float *amp_slice  = &dcimg[imgs->ID_amp].array.F[slice];
        float *spha_slice = &dcimg[imgs->ID_spha].array.F[slice];
        float *samp_slice = &dcimg[imgs->ID_samp].array.F[slice];

        atmturb_init_phase_amp(pha_slice, amp_slice, frame_pixels);
        atmturb_init_phase_amp(spha_slice, samp_slice, frame_pixels);

        for (int k = 0; k < geom->nlayers; k++)
        {
            const atmturb_layer_geom_t *lg = &geom->layers[k];
            const float *scr = dcimg[id_tm[k]].array.F;

            double cur_x = fmod(lg->x0 + (double) t * lg->vx_pix, (double) master_size);
            if (cur_x < 0.0)
            {
                cur_x += (double) master_size;
            }
            double cur_y = fmod(lg->y0 + (double) t * lg->vy_pix, (double) master_size);
            if (cur_y < 0.0)
            {
                cur_y += (double) master_size;
            }

            atmturb_extrude_accumulate(scr, master_size, cur_x, cur_y,
                                       pup_size, (float) lg->weight, pha_slice);

            double cur_sx = fmod(lg->xs0 + (double) t * lg->vx_pix, (double) master_size);
            if (cur_sx < 0.0)
            {
                cur_sx += (double) master_size;
            }
            double cur_sy = fmod(lg->ys0 + (double) t * lg->vy_pix, (double) master_size);
            if (cur_sy < 0.0)
            {
                cur_sy += (double) master_size;
            }

            atmturb_extrude_accumulate(scr, master_size, cur_sx, cur_sy,
                                       pup_size, (float) lg->weight_s, spha_slice);
        }
    }
}

/**
 * atmturb_wfs_save_outputs - Write simulated phase and amplitude cubes to disk
 * @pha_name: Primary phase image stream name.
 * @amp_name: Primary amplitude image stream name.
 */
static void atmturb_wfs_save_outputs(
    const char *pha_name,
    const char *amp_name)
{
    if (CONF_WFOUTPUT == 1)
    {
        char fname_pha[200], fname_amp[200];
        snprintf(fname_pha, sizeof(fname_pha), "%s.fits", pha_name);
        save_fl_fits(pha_name, fname_pha);

        snprintf(fname_amp, sizeof(fname_amp), "%s.fits", amp_name);
        save_fl_fits(amp_name, fname_amp);

        save_fl_fits("outsarraypha", "outsarraypha.fits");
        save_fl_fits("outsarrayamp", "outsarrayamp.fits");
    }
}

/**
 * atmturb_wfs_validate_config - Check simulation geometry parameters before allocation
 *
 * Return: 0 if the configuration is usable, -1 otherwise (diagnostic printed).
 */
static int atmturb_wfs_validate_config(void)
{
    if (!(CONF_PUPIL_SCALE > 0.0f))
    {
        printf("ERROR: PUPIL_SCALE must be > 0 (got %g m/pix)\n", (double) CONF_PUPIL_SCALE);
        return -1;
    }
    if (!(CONF_WFTIME_STEP > 0.0f) || !(CONF_TIME_SPAN > 0.0f))
    {
        printf("ERROR: WFTIME_STEP and TIME_SPAN must be > 0 (got %g s, %g s)\n",
               (double) CONF_WFTIME_STEP, (double) CONF_TIME_SPAN);
        return -1;
    }
    if (CONF_WFsize < 1 || CONF_MASTER_SIZE < CONF_WFsize)
    {
        printf("ERROR: need 1 <= WFsize <= MASTER_SIZE (got WFsize=%ld, MASTER_SIZE=%ld)\n",
               CONF_WFsize, CONF_MASTER_SIZE);
        return -1;
    }
    if (!(CONF_LAMBDA > 0.0f))
    {
        printf("ERROR: TURBULENCE_REF_WAVEL must be > 0\n");
        return -1;
    }
    return 0;
}

/**
 * make_AtmosphericTurbulence_wavefront_series - Run full atmospheric wavefront simulation series
 * @slambdaum: Secondary observing wavelength in um.
 * @WFprecision: Precision mode flag (0=single, 1=double).
 *
 * Return: 0 on success, -1 on failure.
 */
int make_AtmosphericTurbulence_wavefront_series(
    float slambdaum,
    long  WFprecision)
{
    if (CONFFILE[0] != '\0' && AtmosphericTurbulence_ReadConf() != 0)
    {
        return -1;
    }
    if (atmturb_wfs_validate_config() != 0 || !(slambdaum > 0.0f))
    {
        return -1;
    }

    if (slambdaum > 20.0f)
    {
        slambdaum *= 1e-3f; // Convert nm to um
    }

    atmturb_profile_t prof;
    memset(&prof, 0, sizeof(prof));
    if (atmturb_profile_load(CONF_TURBULENCE_PROF_FILE, &prof) != 0)
    {
        return -1;
    }

    atmturb_obs_params_t params;
    memset(&params, 0, sizeof(params));
    params.lambda_ref_m = (double) CONF_LAMBDA;
    params.lambda_s_m = (double) slambdaum * 1e-6;
    params.seeing_arcsec = (double) CONF_SEEING;
    params.zenith_rad = (double) CONF_ZANGLE;
    params.parallactic_rad = (double) CONF_PARALLACTIC_ANGLE;
    params.site_alt_m = (double) CONF_SITE_ALT;
    params.pupil_scale_m = (double) CONF_PUPIL_SCALE;
    params.oversample = 1;
    params.master_size = CONF_MASTER_SIZE;
    params.time_step_s = (double) CONF_WFTIME_STEP;
    params.source_x_rad = (double) CONF_SOURCE_Xpos;
    params.source_y_rad = (double) CONF_SOURCE_Ypos;
    params.seed = CONF_SEED;

    atmturb_geom_t geom;
    memset(&geom, 0, sizeof(geom));
    if (atmturb_geometry_compute(&prof, &params, &geom) != 0)
    {
        atmturb_profile_free(&prof);
        return -1;
    }

    long *id_tm = (long *) calloc((size_t) prof.nlayers, sizeof(long));
    if (id_tm == NULL || atmturb_wfs_load_screens(&prof, &geom, CONF_MASTER_SIZE,
                                                  WFprecision, CONF_SEED, id_tm) != 0)
    {
        free(id_tm);
        atmturb_geometry_free(&geom);
        atmturb_profile_free(&prof);
        return -1;
    }

    long nbframes = (long) (CONF_TIME_SPAN / CONF_WFTIME_STEP + 0.5);
    nbframes = (nbframes < 1) ? 1 : nbframes;
    long pup_size = CONF_WFsize;

    const char *pha_name = (CONF_WF_PHASE_NAME[0] != '\0') ? CONF_WF_PHASE_NAME : "outarraypha";
    const char *amp_name = (CONF_WF_AMPL_NAME[0] != '\0') ? CONF_WF_AMPL_NAME : "outarrayamp";

    atmturb_wfs_images_t imgs;
    imgs.ID_pha  = create_3Dimage_ID(pha_name, pup_size, pup_size, nbframes);
    imgs.ID_amp  = create_3Dimage_ID(amp_name, pup_size, pup_size, nbframes);
    imgs.ID_spha = create_3Dimage_ID("outsarraypha", pup_size, pup_size, nbframes);
    imgs.ID_samp = create_3Dimage_ID("outsarrayamp", pup_size, pup_size, nbframes);

    printf("Synthesizing %ld wavefront frames [%s]\n", nbframes,
           atmturb_simd_active_isa());
    fflush(stdout);

    atmturb_wfs_render_frames(&geom, id_tm, CONF_MASTER_SIZE, pup_size, nbframes, &imgs);
    atmturb_wfs_save_outputs(pha_name, amp_name);

    free(id_tm);
    atmturb_geometry_free(&geom);
    atmturb_profile_free(&prof);
    return 0;
}
