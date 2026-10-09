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
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "AtmosphereModel/AtmosphereModel.h"
#include "AtmosphericTurbulence.h"
#include "atmturb_fresnel.h"
#include "atmturb_geometry.h"
#include "atmturb_lowfreq.h"
#include "atmturb_profile.h"
#include "atmturb_rolling.h"
#include "atmturb_simd.h"
#include "atmturb_types.h"
#include "atmturb_wfs_render.h"
#include "atmturb_wfs_stream.h"

/**
 * atmturb_wfs_notify_streams - Signal shared memory semaphores and increment write counters
 * @imgs: Container of output 3D image handles.
 */
static void atmturb_wfs_notify_streams(
    const atmturb_wfs_images_t *imgs)
{
    const imageID ids[4] = {imgs->ID_pha, imgs->ID_amp, imgs->ID_spha, imgs->ID_samp};
    for (int i = 0; i < 4; i++)
    {
        if (ids[i] >= 0)
        {
            dcimg[ids[i]].md[0].cnt0++;
            COREMOD_MEMORY_image_set_sempost_byID(ids[i], -1);
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
    int save_ref = (CONF_WFOUTPUT & 1);
    int save_sci = (CONF_WFOUTPUT & 2);

    if (!save_ref && !save_sci)
    {
        return;
    }

    char fname[200];

    printf("[milkatmturb] Writing output files to disk:\n");
    if (save_ref)
    {
        snprintf(fname, sizeof(fname), "%s.fits", pha_name);
        printf("[milkatmturb]   - file:   \"%s\" (FITS)\n", fname);
        save_fl_fits(pha_name, fname);

        snprintf(fname, sizeof(fname), "%s.fits", amp_name);
        printf("[milkatmturb]   - file:   \"%s\" (FITS)\n", fname);
        save_fl_fits(amp_name, fname);
    }

    if (save_sci)
    {
        printf("[milkatmturb]   - file:   \"outsarraypha.fits\" (FITS)\n");
        save_fl_fits("outsarraypha", "outsarraypha.fits");

        printf("[milkatmturb]   - file:   \"outsarrayamp.fits\" (FITS)\n");
        save_fl_fits("outsarrayamp", "outsarrayamp.fits");
    }
    fflush(stdout);
}

/**
 * atmturb_wfs_validate_config - Check simulation geometry parameters before allocation
 *
 * Return: 0 if the configuration is usable, -1 otherwise (diagnostic printed).
 */
int atmturb_wfs_validate_config(void)
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
    long os = (CONF_OVERSAMPLE > 1) ? (long) CONF_OVERSAMPLE : 1L;
    if (CONF_WFsize < 1 || CONF_MASTER_SIZE < CONF_WFsize * os)
    {
        printf("ERROR: need 1 <= WFsize * os <= MASTER_SIZE "
               "(got WFsize=%ld, os=%ld, MASTER_SIZE=%ld)\n",
               CONF_WFsize, os, CONF_MASTER_SIZE);
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
 * atmturb_wfs_init_obs_params - Populate observation parameters from global configuration
 * @slambdaum: Secondary observing wavelength in um.
 * @params: Observation parameters container to populate.
 */
void atmturb_wfs_init_obs_params(
    float                 slambdaum,
    atmturb_obs_params_t *params)
{
    memset(params, 0, sizeof(*params));
    params->lambda_ref_m    = (double) CONF_LAMBDA;
    params->lambda_s_m      = (double) slambdaum * 1e-6;
    params->seeing_arcsec   = (double) CONF_SEEING;
    params->zenith_rad      = (double) CONF_ZANGLE;
    params->parallactic_rad = (double) CONF_PARALLACTIC_ANGLE;
    params->site_alt_m      = (double) CONF_SITE_ALT;
    params->pupil_scale_m   = (double) CONF_PUPIL_SCALE;
    params->oversample      = (CONF_OVERSAMPLE > 1) ? CONF_OVERSAMPLE : 1;
    params->interp          = CONF_INTERP;
    params->lowfreq         = CONF_LOWFREQ;
    params->rolling         = CONF_ROLLING;
    params->boil_time_s     = (double) CONF_BOIL_TIME;
    params->master_size     = CONF_MASTER_SIZE;
    long nbframes           = (long) (CONF_TIME_SPAN / CONF_WFTIME_STEP + 0.5);
    params->nbframes        = (nbframes < 1) ? 1 : nbframes;
    params->time_step_s     = (double) CONF_WFTIME_STEP;
    params->source_x_rad    = (double) CONF_SOURCE_Xpos;
    params->source_y_rad    = (double) CONF_SOURCE_Ypos;
    params->seed            = CONF_SEED;
}

/**
 * atmturb_wfs_setup_sim - Initialize atmospheric simulation structures
 * @slambdaum: Secondary observing wavelength in um.
 * @WFprecision: Precision mode flag.
 * @nbframes: Number of frames (or estimated pool frames).
 * @prof: Destination profile struct.
 * @params: Destination observation params struct.
 * @geom: Destination geometry struct.
 * @rsim: Destination rolling screens struct.
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_wfs_setup_sim(
    float                 slambdaum,
    long                  WFprecision,
    long                  nbframes,
    atmturb_profile_t    *prof,
    atmturb_obs_params_t *params,
    atmturb_geom_t       *geom,
    atmturb_rolling_t    *rsim)
{
    if ((CONFFILE[0] != '\0' && AtmosphericTurbulence_ReadConf() != 0) ||
        atmturb_wfs_validate_config() != 0 || !(slambdaum > 0.0f))
    {
        return -1;
    }
    if (slambdaum > 20.0f)
    {
        slambdaum *= 1e-3f;
    }

    memset(prof, 0, sizeof(*prof));
    if (atmturb_profile_load(CONF_TURBULENCE_PROF_FILE, prof) != 0)
    {
        return -1;
    }

    atmturb_wfs_init_obs_params(slambdaum, params);

    memset(geom, 0, sizeof(*geom));
    if (atmturb_geometry_compute(prof, params, geom) != 0)
    {
        atmturb_profile_free(prof);
        return -1;
    }

    long pool_frames = (nbframes > 0) ? nbframes : params->nbframes;
    if (atmturb_rolling_init(rsim, prof, geom, CONF_MASTER_SIZE,
                             pool_frames, params->time_step_s, WFprecision, CONF_SEED) != 0)
    {
        atmturb_geometry_free(geom);
        atmturb_profile_free(prof);
        return -1;
    }
    return 0;
}

/**
 * atmturb_wfs_teardown_sim - Free atmospheric simulation structures
 * @prof: Active profile struct.
 * @geom: Active geometry struct.
 * @rsim: Active rolling screens struct.
 */
void atmturb_wfs_teardown_sim(
    atmturb_profile_t *prof,
    atmturb_geom_t    *geom,
    atmturb_rolling_t *rsim)
{
    atmturb_rolling_free(rsim);
    atmturb_geometry_free(geom);
    atmturb_profile_free(prof);
}

/**
 * atmturb_wfs_print_output_targets - Display planned simulation output targets and types
 * @pha_name: Primary phase image stream name.
 * @amp_name: Primary amplitude image stream name.
 * @pup_size: Linear dimension of square pupil grid in pixels.
 * @nbframes: Number of simulated 3D cube frames.
 * @precision: Floating point precision flag (0=float32, 1=float64).
 */
static void atmturb_wfs_print_output_targets(
    const char *pha_name,
    const char *amp_name,
    long        pup_size,
    long        nbframes,
    long        precision)
{
    const char *prec_str = (precision == 1) ? "float64" : "float32";

    printf("[milkatmturb] Outputs:\n");
    printf("[milkatmturb]   - stream: \"%s\" (3D SHM, %ldx%ldx%ld %s)\n",
           pha_name, pup_size, pup_size, nbframes, prec_str);
    printf("[milkatmturb]   - stream: \"%s\" (3D SHM, %ldx%ldx%ld %s)\n",
           amp_name, pup_size, pup_size, nbframes, prec_str);
    printf("[milkatmturb]   - stream: \"outsarraypha\" (3D SHM, %ldx%ldx%ld %s)\n",
           pup_size, pup_size, nbframes, prec_str);
    printf("[milkatmturb]   - stream: \"outsarrayamp\" (3D SHM, %ldx%ldx%ld %s)\n",
           pup_size, pup_size, nbframes, prec_str);

    int save_ref = (CONF_WFOUTPUT & 1);
    int save_sci = (CONF_WFOUTPUT & 2);

    if (save_ref)
    {
        printf("[milkatmturb]   - file:   \"%s.fits\" (FITS 3D cube)\n", pha_name);
        printf("[milkatmturb]   - file:   \"%s.fits\" (FITS 3D cube)\n", amp_name);
    }
    if (save_sci)
    {
        printf("[milkatmturb]   - file:   \"outsarraypha.fits\" (FITS 3D cube)\n");
        printf("[milkatmturb]   - file:   \"outsarrayamp.fits\" (FITS 3D cube)\n");
    }
    if (!save_ref && !save_sci)
    {
        printf("[milkatmturb]   - file:   none (save_fits = 0)\n");
    }
    fflush(stdout);
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
    if (CONF_STREAM_MODE > 0)
    {
        return make_AtmosphericTurbulence_wavefront_stream(slambdaum, WFprecision,
                                                           CONF_STREAM_MODE);
    }

    atmturb_profile_t    prof;
    atmturb_obs_params_t params;
    atmturb_geom_t       geom;
    atmturb_rolling_t    rsim;

    if (atmturb_wfs_setup_sim(slambdaum, WFprecision, 0, &prof, &params, &geom, &rsim) != 0)
    {
        return -1;
    }

    long nbframes = params.nbframes;
    long pup_size = CONF_WFsize;

    const char *pha_name = (CONF_WF_PHASE_NAME[0] != '\0') ? CONF_WF_PHASE_NAME : "outarraypha";
    const char *amp_name = (CONF_WF_AMPL_NAME[0] != '\0') ? CONF_WF_AMPL_NAME : "outarrayamp";

    atmturb_wfs_images_t imgs;
    imgs.ID_pha  = create_3Dimage_ID(pha_name, pup_size, pup_size, nbframes);
    imgs.ID_amp  = create_3Dimage_ID(amp_name, pup_size, pup_size, nbframes);
    imgs.ID_spha = create_3Dimage_ID("outsarraypha", pup_size, pup_size, nbframes);
    imgs.ID_samp = create_3Dimage_ID("outsarrayamp", pup_size, pup_size, nbframes);

    atmturb_wfs_print_output_targets(pha_name, amp_name, pup_size, nbframes, WFprecision);

    printf("Synthesizing %ld wavefront frames [%s]\n", nbframes,
           atmturb_simd_active_isa());
    fflush(stdout);

    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    atmturb_wfs_render_frames(&rsim, &geom, &prof, &params, CONF_MASTER_SIZE, pup_size,
                              nbframes, &imgs);
    atmturb_wfs_notify_streams(&imgs);
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    double render_time = (double) (ts1.tv_sec - ts0.tv_sec) +
                         (double) (ts1.tv_nsec - ts0.tv_nsec) * 1e-9;
    double fps = (render_time > 0.0) ? ((double) nbframes / render_time) : 0.0;
    double mps = fps * (double) (pup_size * pup_size) * 1e-6;
    printf("[milkatmturb] Rendered %ld frames (%ldx%ld) in %.3f s (%.1f fps, %.2f MP/s)\n",
           nbframes, pup_size, pup_size, render_time, fps, mps);
    fflush(stdout);
    atmturb_wfs_save_outputs(pha_name, amp_name);
    atmturb_wfs_teardown_sim(&prof, &geom, &rsim);
    return 0;
}
