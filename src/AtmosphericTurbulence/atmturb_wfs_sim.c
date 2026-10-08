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
#include "atmturb_wfs_stream.h"

#ifdef HAVE_CUDA
#include "atmturb_cuda.h"
#endif

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
 * atmturb_wfs_extrude_channel - Extrude a single turbulence layer phase into pupil
 * @scr: Master screen pointer.
 * @master_size: Master screen linear dimension.
 * @x: Unwrapped X position in master pixels.
 * @y: Unwrapped Y position in master pixels.
 * @pup_size: Linear dimension of pupil.
 * @geom: Observing geometry context.
 * @weight: Scaled layer weight.
 * @out_pha: Destination phase array to accumulate into.
 */
static inline void atmturb_wfs_extrude_channel(
    const float          *scr,
    long                  master_size,
    double                x,
    double                y,
    long                  pup_size,
    const atmturb_geom_t *geom,
    float                 weight,
    float                *out_pha)
{
    double cur_x = fmod(x, (double) master_size);
    if (cur_x < 0.0)
    {
        cur_x += (double) master_size;
    }
    double cur_y = fmod(y, (double) master_size);
    if (cur_y < 0.0)
    {
        cur_y += (double) master_size;
    }

    atmturb_extrude_params_t ep;
    ep.master   = scr;
    ep.msize    = master_size;
    ep.x0       = cur_x;
    ep.y0       = cur_y;
    ep.pup_size = pup_size;
    ep.os       = (long) geom->oversample;
    ep.interp   = geom->interp;
    ep.weight   = weight;
    ep.out_pha  = out_pha;
    atmturb_extrude_accumulate(&ep);
}

/**
 * atmturb_wfs_render_layer - Render one turbulence layer into phase slices for frame t
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @k: Layer index.
 * @t: Frame index.
 * @time_step_s: Time step in seconds.
 * @master_size: Master screen dimension.
 * @pup_size: Pupil dimension.
 * @pha_slice: Primary phase frame slice.
 * @spha_slice: Secondary phase frame slice.
 */
void atmturb_wfs_render_layer(
    const atmturb_rolling_t *r,
    const atmturb_geom_t    *geom,
    int                      k,
    long                     t,
    double                   time_step_s,
    long                     master_size,
    long                     pup_size,
    float                   *pha_slice,
    float                   *spha_slice)
{
    const atmturb_layer_geom_t *lg = &geom->layers[k];
    atmturb_rolling_eval_t rev;
    atmturb_rolling_get_frame(r, k, t, time_step_s, &rev);

    double dx = (lg->traj_x != NULL) ? lg->traj_x[t] : ((double) t * lg->vx_pix);
    double dy = (lg->traj_y != NULL) ? lg->traj_y[t] : ((double) t * lg->vy_pix);

    double x = lg->x0 + dx;
    double y = lg->y0 + dy;
    atmturb_wfs_extrude_channel(rev.scrA, master_size, x, y, pup_size, geom,
                                (float) (lg->weight * (double) rev.wA), pha_slice);
    if (rev.wB > 0.0f && rev.scrB != NULL)
    {
        atmturb_wfs_extrude_channel(rev.scrB, master_size, x, y, pup_size, geom,
                                    (float) (lg->weight * (double) rev.wB), pha_slice);
    }

    double xs = lg->xs0 + dx;
    double ys = lg->ys0 + dy;
    atmturb_wfs_extrude_channel(rev.scrA, master_size, xs, ys, pup_size, geom,
                                (float) (lg->weight_s * (double) rev.wA), spha_slice);
    if (rev.wB > 0.0f && rev.scrB != NULL)
    {
        atmturb_wfs_extrude_channel(rev.scrB, master_size, xs, ys, pup_size, geom,
                                    (float) (lg->weight_s * (double) rev.wB), spha_slice);
    }

    if (r->lowfreq)
    {
        atmturb_lowfreq_accumulate_custom(&r->layers[k].lf_base, rev.are_eff, rev.aim_eff,
                                          x, y, pup_size, (long) geom->oversample,
                                          (float) lg->weight, pha_slice);

        atmturb_lowfreq_accumulate_custom(&r->layers[k].lf_base, rev.are_eff, rev.aim_eff,
                                          xs, ys, pup_size, (long) geom->oversample,
                                          (float) lg->weight_s, spha_slice);
    }
}

/**
 * atmturb_wfs_render_geometric - Multi-threaded geometric rendering of simulation time steps
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @master_size: Master screen dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 * @nbframes: Number of simulation frames.
 * @time_step_s: Time step between frames in seconds.
 * @imgs: Container of output 3D image handles.
 */
static void atmturb_wfs_render_geometric(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    long                        master_size,
    long                        pup_size,
    long                        nbframes,
    double                      time_step_s,
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
            atmturb_wfs_render_layer(r, geom, k, t, time_step_s, master_size,
                                     pup_size, pha_slice, spha_slice);
        }
    }
}

/**
 * atmturb_wfs_render_diffractive - Multi-threaded diffractive rendering of simulation time steps
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @prof: Active turbulence profile.
 * @params: Observation parameters container.
 * @master_size: Master screen dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 * @nbframes: Number of simulation frames.
 * @imgs: Container of output 3D image handles.
 *
 * Return: 0 on success, -1 on plan or context allocation failure.
 */
static int atmturb_wfs_render_diffractive(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    const atmturb_profile_t    *prof,
    const atmturb_obs_params_t *params,
    long                        master_size,
    long                        pup_size,
    long                        nbframes,
    const atmturb_wfs_images_t *imgs)
{
    atmturb_fresnel_plan_t plan;
    double z_bin = (double) CONF_FRESNEL_PROPAGATION_BIN;
    if (atmturb_fresnel_plan_init(&plan, prof, geom, pup_size, params->pupil_scale_m,
                                  params->lambda_ref_m, params->lambda_s_m, z_bin) != 0)
    {
        return -1;
    }

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    atmturb_fresnel_ctx_t *ctxs = (atmturb_fresnel_ctx_t *) calloc((size_t) nthreads,
                                                                   sizeof(atmturb_fresnel_ctx_t));
    if (ctxs == NULL)
    {
        atmturb_fresnel_plan_free(&plan);
        return -1;
    }

    for (int tid = 0; tid < nthreads; tid++)
    {
        if (atmturb_fresnel_ctx_init(&ctxs[tid], pup_size) != 0)
        {
            for (int j = 0; j < tid; j++)
            {
                atmturb_fresnel_ctx_free(&ctxs[j]);
            }
            free(ctxs);
            atmturb_fresnel_plan_free(&plan);
            return -1;
        }
    }

    long frame_pixels = pup_size * pup_size;

    #pragma omp parallel for schedule(dynamic)
    for (long t = 0; t < nbframes; t++)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        long slice = t * frame_pixels;
        float *pha_slice  = &dcimg[imgs->ID_pha].array.F[slice];
        float *amp_slice  = &dcimg[imgs->ID_amp].array.F[slice];
        float *spha_slice = &dcimg[imgs->ID_spha].array.F[slice];
        float *samp_slice = &dcimg[imgs->ID_samp].array.F[slice];

        atmturb_fresnel_render_step(&ctxs[tid], &plan, r, geom, t, params->time_step_s,
                                    master_size, pup_size, pha_slice, amp_slice,
                                    spha_slice, samp_slice);
    }

    for (int tid = 0; tid < nthreads; tid++)
    {
        atmturb_fresnel_ctx_free(&ctxs[tid]);
    }
    free(ctxs);
    atmturb_fresnel_plan_free(&plan);

    return 0;
}

#ifdef HAVE_CUDA
/**
 * atmturb_wfs_render_cuda - Render simulation frames using CUDA GPU acceleration
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @master_size: Master screen dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 * @nbframes: Number of simulation frames.
 * @imgs: Container of output 3D image handles.
 *
 * Return: 0 on success, -1 on failure.
 */
static int atmturb_wfs_render_cuda(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    long                        master_size,
    long                        pup_size,
    long                        nbframes,
    const atmturb_wfs_images_t *imgs)
{
    int nl = geom->nlayers;
    size_t m_sz = (size_t) nl * sizeof(float *);
    size_t d_sz = (size_t) nl * sizeof(double);
    void *buf = malloc(m_sz + 5 * d_sz);
    if (!buf)
    {
        return -1;
    }
    const float **masters = (const float **) buf;
    double *vx  = (double *) ((char *) buf + m_sz);
    double *vy  = vx + nl;
    double *cn2 = vy + nl;
    double *x0  = cn2 + nl;
    double *y0  = x0 + nl;

    for (int k = 0; k < nl; k++)
    {
        masters[k] = r->layers[k].screens[0].data;
        vx[k]      = geom->layers[k].vx_pix;
        vy[k]      = geom->layers[k].vy_pix;
        cn2[k]     = geom->layers[k].weight * geom->layers[k].weight;
        x0[k]      = geom->layers[k].x0;
        y0[k]      = geom->layers[k].y0;
    }

    double w0 = geom->layers[0].weight;
    double scoeff = (w0 > 0.0) ? (geom->layers[0].weight_s / w0) : 1.0;
    atmturb_cuda_sim_params_t cparams = {
        .nblayers  = nl,
        .msize     = master_size,
        .pup_size  = pup_size,
        .nbframes  = nbframes,
        .Scoeff    = scoeff,
        .h_masters = (const float *const *) masters,
        .vxpix     = vx,
        .vypix     = vy,
        .cn2       = cn2,
        .x0        = x0,
        .y0        = y0
    };
    long total_pixels = nbframes * pup_size * pup_size;
    float *amp_ptr  = &dcimg[imgs->ID_amp].array.F[0];
    float *samp_ptr = &dcimg[imgs->ID_samp].array.F[0];

    #pragma omp parallel for
    for (long i = 0; i < total_pixels; i++)
    {
        amp_ptr[i]  = 1.0f;
        samp_ptr[i] = 1.0f;
    }

    atmturb_cuda_sim_outputs_t outputs = {
        .pha  = &dcimg[imgs->ID_pha].array.F[0],
        .amp  = NULL,
        .spha = &dcimg[imgs->ID_spha].array.F[0],
        .samp = NULL
    };

    int ret = atmturb_wfs_render_frames_cuda(&cparams, &outputs);
    free(buf);
    return ret;
}
#endif

/**
 * atmturb_wfs_render_frames - Dispatch rendering to diffractive or geometric engine
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @prof: Active turbulence profile.
 * @params: Observation parameters container.
 * @master_size: Master screen dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 * @nbframes: Number of simulation frames.
 * @imgs: Container of output 3D image handles.
 */
static void atmturb_wfs_render_frames(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    const atmturb_profile_t    *prof,
    const atmturb_obs_params_t *params,
    long                        master_size,
    long                        pup_size,
    long                        nbframes,
    const atmturb_wfs_images_t *imgs)
{
    if (CONF_FRESNEL_PROPAGATION == 1)
    {
        if (atmturb_wfs_render_diffractive(r, geom, prof, params, master_size,
                                           pup_size, nbframes, imgs) == 0)
        {
            return;
        }
    }

#ifdef HAVE_CUDA
    if (atmturb_simd_is_gpu())
    {
        if (atmturb_wfs_render_cuda(r, geom, master_size, pup_size, nbframes, imgs) == 0)
        {
            return;
        }
    }
#endif

    atmturb_wfs_render_geometric(r, geom, master_size, pup_size, nbframes,
                                 params->time_step_s, imgs);
}

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
    if (CONF_WFOUTPUT == 1)
    {
        char fname[200];
        snprintf(fname, sizeof(fname), "%s.fits", pha_name);
        save_fl_fits(pha_name, fname);
        snprintf(fname, sizeof(fname), "%s.fits", amp_name);
        save_fl_fits(amp_name, fname);
        save_fl_fits("outsarraypha", "outsarraypha.fits");
        save_fl_fits("outsarrayamp", "outsarrayamp.fits");
    }
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
