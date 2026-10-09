// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_wfs_stream.c
 * @brief   Frame-by-frame 2D shared memory wavefront streaming engine
 */

#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "AtmosphereModel/AtmosphereModel.h"
#include "AtmosphericTurbulence.h"
#include "atmturb_compat.h"
#include "atmturb_fresnel.h"
#include "atmturb_geometry.h"
#include "atmturb_profile.h"
#include "atmturb_rolling.h"
#include "atmturb_rytov.h"
#include "atmturb_simd.h"
#include "atmturb_types.h"
#include "atmturb_wfs_render.h"
#include "atmturb_wfs_stream.h"
#include "ImageStreamIO/ImageStreamIO.h"

static volatile sig_atomic_t g_stream_stop = 0;

/**
 * atmturb_wfs_stream_sig_handler - Signal handler to request graceful stream loop stop
 * @sig: Signal number.
 */
static void atmturb_wfs_stream_sig_handler(int sig)
{
    (void) sig;
    g_stream_stop = 1;
}

/**
 * atmturb_wfs_stream_ensure_image - Validate or create a 2D float shared memory stream
 * @name: Shared memory image stream name.
 * @pup_size: Linear dimension of square image.
 *
 * Return: Valid Image ID, or -1 on creation failure.
 */
static imageID atmturb_wfs_stream_ensure_image(
    const char *name,
    long        pup_size)
{
    if (name == NULL || name[0] == '\0')
    {
        return -1;
    }
    imageID id = image_ID(name);
    if (id >= 0 && dcimg[id].used == 1)
    {
        if (dcimg[id].md[0].naxis == 2 &&
            (long) dcimg[id].md[0].size[0] == pup_size &&
            (long) dcimg[id].md[0].size[1] == pup_size &&
            dcimg[id].md[0].datatype == _DATATYPE_FLOAT)
        {
            return id;
        }
        delete_image_ID(name);
        id = -1;
    }

    uint32_t sz[2] = {(uint32_t) pup_size, (uint32_t) pup_size};
    create_image_ID(name, 2, sz, _DATATYPE_FLOAT, 1, 10, 0, &id);
    return id;
}

/**
 * atmturb_wfs_stream_post - Update stream metadata and post semaphores
 * @id: Image ID to update.
 * @ts: Frame timestamp.
 */
static inline void atmturb_wfs_stream_post(
    imageID                id,
    const struct timespec *ts)
{
    if (id < 0 || id >= dcnimg || dcimg[id].used != 1)
    {
        return;
    }
    dcimg[id].md[0].cnt1 = 0;
    if (dcimg[id].md[0].shared == 1)
    {
        ImageStreamIO_UpdateIm_atime(&dcimg[id], (struct timespec *) ts);
    }
    else
    {
        dcimg[id].md[0].writetime = *ts;
        dcimg[id].md[0].atime     = *ts;
        dcimg[id].md[0].cnt0++;
        dcimg[id].md[0].write     = 0;
        ImageStreamIO_sempost(&dcimg[id], -1);
    }
}

/**
 * struct atmturb_stream_state_t - Container for stream execution state
 * @id_pha: Phase image ID.
 * @id_amp: Amplitude image ID.
 * @id_spha: Secondary phase image ID.
 * @id_samp: Secondary amplitude image ID.
 * @pha: Phase pointer.
 * @amp: Amplitude pointer.
 * @spha: Secondary phase pointer.
 * @samp: Secondary amplitude pointer.
 * @fplan: Fresnel plan pointer.
 * @fctx: Fresnel context pointer.
 * @rplan: Rytov plan pointer.
 * @rctx: Rytov context pointer.
 */
typedef struct
{
    imageID                       id_pha;
    imageID                       id_amp;
    imageID                       id_spha;
    imageID                       id_samp;
    float                        *pha;
    float                        *amp;
    float                        *spha;
    float                        *samp;
    const atmturb_fresnel_plan_t *fplan;
    atmturb_fresnel_ctx_t        *fctx;
    const atmturb_rytov_plan_t   *rplan;
    atmturb_rytov_ctx_t          *rctx;
} atmturb_stream_state_t;

/**
 * atmturb_wfs_stream_render_step - Compute a single simulation frame
 * @rsim: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @st: Stream state container.
 * @t: Frame index.
 * @time_step_s: Time step in seconds.
 * @master_size: Master screen linear dimension.
 * @pup_size: Linear pupil dimension.
 */
static void atmturb_wfs_stream_render_step(
    const atmturb_rolling_t *rsim,
    const atmturb_geom_t    *geom,
    atmturb_stream_state_t  *st,
    long                     t,
    double                   time_step_s,
    long                     master_size,
    long                     pup_size)
{
    if (CONF_FRESNEL_PROPAGATION == 1 && st->fplan != NULL && st->fctx != NULL)
    {
        atmturb_fresnel_render_step(st->fctx, st->fplan, rsim, geom, t, time_step_s,
                                    master_size, pup_size, st->pha, st->amp,
                                    st->spha, st->samp);
        return;
    }

    if (CONF_FRESNEL_PROPAGATION == 2 && st->rplan != NULL && st->rctx != NULL)
    {
        atmturb_rytov_render_step(st->rctx, st->rplan, rsim, geom, t, time_step_s,
                                  master_size, pup_size, st->pha, st->amp,
                                  st->spha, st->samp);
        return;
    }

    long frame_pixels = pup_size * pup_size;
    atmturb_init_phase_amp(st->pha, st->amp, frame_pixels);
    atmturb_init_phase_amp(st->spha, st->samp, frame_pixels);

    for (int k = 0; k < geom->nlayers; k++)
    {
        atmturb_wfs_render_layer(rsim, geom, k, t, time_step_s, master_size,
                                 pup_size, st->pha, st->spha);
    }
}

/**
 * atmturb_wfs_stream_pace - Delay until the target frame time
 * @start_ts: Simulation start monotonic timestamp.
 * @next_frame_idx: Next frame index.
 * @time_step_s: Time step between frames in seconds.
 */
static inline void atmturb_wfs_stream_pace(
    struct timespec start_ts,
    long            next_frame_idx,
    double          time_step_s)
{
    double target_sec = (double) next_frame_idx * time_step_s;
    int64_t target_nsec = (int64_t) start_ts.tv_nsec + (int64_t) (target_sec * 1e9);
    struct timespec target_ts;
    target_ts.tv_sec  = start_ts.tv_sec + (time_t) (target_nsec / 1000000000LL);
    target_ts.tv_nsec = (long) (target_nsec % 1000000000LL);
    clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &target_ts, NULL);
}

/**
 * atmturb_wfs_stream_remove_piston - Zero pupil-averaged phase across frame
 * @pha: Phase array to zero mean.
 * @npix: Total pixel count.
 */
static void atmturb_wfs_stream_remove_piston(
    float *pha,
    long   npix)
{
    if (pha == NULL || npix <= 0)
    {
        return;
    }
    double sum = 0.0;
    for (long i = 0; i < npix; i++)
    {
        sum += (double) pha[i];
    }
    float mean = (float) (sum / (double) npix);
    for (long i = 0; i < npix; i++)
    {
        pha[i] -= mean;
    }
}

/**
 * atmturb_wfs_stream_publish - Commit scratch buffers to shared memory and notify readers
 * @st: Stream state container.
 * @frame_pixels: Total pixels in frame.
 * @ts: Frame acquisition timestamp.
 */
static void atmturb_wfs_stream_publish(
    const atmturb_stream_state_t *st,
    long                          frame_pixels,
    const struct timespec        *ts)
{
    size_t nbytes = (size_t) frame_pixels * sizeof(float);

    if (st->id_pha >= 0)
    {
        atmturb_wfs_stream_remove_piston(st->pha, frame_pixels);
        dcimg[st->id_pha].md[0].write = 1;
        memcpy(dcimg[st->id_pha].array.F, st->pha, nbytes);
        atmturb_wfs_stream_post(st->id_pha, ts);
    }
    if (st->id_amp >= 0)
    {
        dcimg[st->id_amp].md[0].write = 1;
        memcpy(dcimg[st->id_amp].array.F, st->amp, nbytes);
        atmturb_wfs_stream_post(st->id_amp, ts);
    }
    if (st->id_spha >= 0)
    {
        atmturb_wfs_stream_remove_piston(st->spha, frame_pixels);
        dcimg[st->id_spha].md[0].write = 1;
        memcpy(dcimg[st->id_spha].array.F, st->spha, nbytes);
        atmturb_wfs_stream_post(st->id_spha, ts);
    }
    if (st->id_samp >= 0)
    {
        dcimg[st->id_samp].md[0].write = 1;
        memcpy(dcimg[st->id_samp].array.F, st->samp, nbytes);
        atmturb_wfs_stream_post(st->id_samp, ts);
    }
}

/**
 * atmturb_wfs_stream_loop - Inner frame-by-frame streaming loop
 * @rsim: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @params: Observation parameters container.
 * @st: Stream state container.
 * @pup_size: Output pupil linear dimension.
 * @max_frames: Maximum frames to render (0 = infinite).
 * @pace: Flag to pace with simulated time_step.
 *
 * Return: Total frames streamed.
 */
static long atmturb_wfs_stream_loop(
    const atmturb_rolling_t    *rsim,
    const atmturb_geom_t       *geom,
    const atmturb_obs_params_t *params,
    atmturb_stream_state_t     *st,
    long                        pup_size,
    long                        max_frames,
    int                         pace)
{
    long frame_pixels = pup_size * pup_size;
    long t = 0;
    struct timespec start_ts;
    clock_gettime(CLOCK_MONOTONIC, &start_ts);

    for (t = 0; ; t++)
    {
        if (g_stream_stop != 0 || (max_frames > 0 && t >= max_frames))
        {
            break;
        }

        atmturb_wfs_stream_render_step(rsim, geom, st, t, params->time_step_s,
                                       CONF_MASTER_SIZE, pup_size);

        struct timespec ts_now;
        clock_gettime(CLOCK_REALTIME, &ts_now);

        atmturb_wfs_stream_publish(st, frame_pixels, &ts_now);

        if (pace && params->time_step_s > 0.0)
        {
            atmturb_wfs_stream_pace(start_ts, t + 1, params->time_step_s);
        }
    }
    return t;
}

/**
 * atmturb_wfs_stream_init_buffers - Allocate and map 2D streaming buffers
 * @st: Stream state container.
 * @pup_size: Linear dimension of square pupil.
 * @pha_name: Phase stream name.
 * @amp_name: Amplitude stream name.
 *
 * Return: Scratch buffer pointer, or NULL on creation failure.
 */
static float *atmturb_wfs_stream_init_buffers(
    atmturb_stream_state_t *st,
    long                    pup_size,
    const char             *pha_name,
    const char             *amp_name)
{
    st->id_pha  = atmturb_wfs_stream_ensure_image(pha_name, pup_size);
    st->id_amp  = (CONF_WAVEFRONT_AMPLITUDE == 1)
                      ? atmturb_wfs_stream_ensure_image(amp_name, pup_size) : -1;
    st->id_spha = (CONF_MAKE_SWAVEFRONT == 1)
                      ? atmturb_wfs_stream_ensure_image("outsarraypha", pup_size) : -1;
    st->id_samp = (CONF_MAKE_SWAVEFRONT == 1 && CONF_WAVEFRONT_AMPLITUDE == 1)
                      ? atmturb_wfs_stream_ensure_image("outsarrayamp", pup_size) : -1;

    if (st->id_pha < 0)
    {
        return NULL;
    }

    long frame_pixels = pup_size * pup_size;
    float *scratch = (float *) calloc((size_t) (frame_pixels * 4), sizeof(float));
    if (scratch == NULL)
    {
        return NULL;
    }
    st->pha  = scratch;
    st->amp  = scratch + frame_pixels;
    st->spha = scratch + 2 * frame_pixels;
    st->samp = scratch + 3 * frame_pixels;
    return scratch;
}

/**
 * atmturb_wfs_stream_print_output_targets - Display streaming output targets and types
 * @pha_name: Phase stream name.
 * @amp_name: Amplitude stream name.
 * @pup_size: Linear pupil dimension.
 * @st: Stream state container.
 */
static void atmturb_wfs_stream_print_output_targets(
    const char                   *pha_name,
    const char                   *amp_name,
    long                          pup_size,
    const atmturb_stream_state_t *st)
{
    printf("[milkatmturb] Outputs:\n");
    printf("[milkatmturb]   - stream: \"%s\" (2D SHM, %ldx%ld float32)\n",
           pha_name, pup_size, pup_size);
    if (st->id_amp >= 0)
    {
        printf("[milkatmturb]   - stream: \"%s\" (2D SHM, %ldx%ld float32)\n",
               amp_name, pup_size, pup_size);
    }
    if (st->id_spha >= 0)
    {
        printf("[milkatmturb]   - stream: \"outsarraypha\" (2D SHM, %ldx%ld float32)\n",
               pup_size, pup_size);
    }
    if (st->id_samp >= 0)
    {
        printf("[milkatmturb]   - stream: \"outsarrayamp\" (2D SHM, %ldx%ld float32)\n",
               pup_size, pup_size);
    }
    printf("[milkatmturb]   - file:   none (streaming mode)\n");
    fflush(stdout);
}

/**
 * make_AtmosphericTurbulence_wavefront_stream - Stream 2D wavefront frames to SHM
 * @slambdaum: Secondary observing wavelength in um.
 * @WFprecision: Precision mode flag (0=single, 1=double).
 * @stream_mode: Streaming mode (1=continuous stream, 2=finite stream, 3=unpaced continuous).
 *
 * Return: 0 on success, -1 on failure.
 */
int make_AtmosphericTurbulence_wavefront_stream(
    float slambdaum,
    long  WFprecision,
    int   stream_mode)
{
    atmturb_profile_t    prof;
    atmturb_obs_params_t params;
    atmturb_geom_t       geom;
    atmturb_rolling_t    rsim;

    long max_frames = (stream_mode == 2) ? (long) (CONF_TIME_SPAN / CONF_WFTIME_STEP + 0.5) : 0;
    long pool_frames = (max_frames > 0) ? max_frames : 200;

    if (atmturb_wfs_setup_sim(slambdaum, WFprecision, pool_frames, &prof, &params,
                              &geom, &rsim) != 0)
    {
        return -1;
    }

    long pup_size = CONF_WFsize;
    const char *pha_name = (CONF_WF_PHASE_NAME[0] != '\0') ? CONF_WF_PHASE_NAME : "outarraypha";
    const char *amp_name = (CONF_WF_AMPL_NAME[0] != '\0') ? CONF_WF_AMPL_NAME : "outarrayamp";

    atmturb_stream_state_t st;
    memset(&st, 0, sizeof(st));
    float *scratch = atmturb_wfs_stream_init_buffers(&st, pup_size, pha_name, amp_name);
    if (scratch == NULL)
    {
        atmturb_wfs_teardown_sim(&prof, &geom, &rsim);
        return -1;
    }

    atmturb_wfs_stream_print_output_targets(pha_name, amp_name, pup_size, &st);

    atmturb_fresnel_plan_t fplan;
    atmturb_fresnel_ctx_t  fctx;
    atmturb_rytov_plan_t   rplan;
    atmturb_rytov_ctx_t    rctx;
    if (CONF_FRESNEL_PROPAGATION == 1)
    {
        double z_bin = (double) CONF_FRESNEL_PROPAGATION_BIN;
        atmturb_fresnel_plan_init(&fplan, &prof, &geom, pup_size, params.pupil_scale_m,
                                  params.lambda_ref_m, params.lambda_s_m, z_bin);
        atmturb_fresnel_ctx_init(&fctx, pup_size);
        st.fplan = &fplan;
        st.fctx  = &fctx;
    }
    else if (CONF_FRESNEL_PROPAGATION == 2)
    {
        double z_bin = (double) CONF_FRESNEL_PROPAGATION_BIN;
        atmturb_rytov_plan_init(&rplan, &prof, &geom, pup_size, params.pupil_scale_m,
                                params.lambda_ref_m, params.lambda_s_m, z_bin,
                                CONF_FRESNEL_RYTOV_SEC_EXACT);
        atmturb_rytov_ctx_init(&rctx, rplan.pad_size);
        st.rplan = &rplan;
        st.rctx  = &rctx;
    }

    g_stream_stop = 0;
    struct sigaction sa, old_sa_int, old_sa_term, old_sa_quit;
    memset(&sa, 0, sizeof(sa));
    sa.sa_handler = atmturb_wfs_stream_sig_handler;
    sigaction(SIGINT, &sa, &old_sa_int);
    sigaction(SIGTERM, &sa, &old_sa_term);
    sigaction(SIGQUIT, &sa, &old_sa_quit);

    printf("[milkatmturb] Streaming 2D wavefront frames (%ldx%ld pix) to SHM '%s'%s\n",
           pup_size, pup_size, pha_name, (st.id_amp >= 0) ? " and amplitude" : "");
    fflush(stdout);

    int pace = (stream_mode == 1 || stream_mode == 2);
    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    long n_streamed = atmturb_wfs_stream_loop(&rsim, &geom, &params, &st, pup_size,
                                             max_frames, pace);
    clock_gettime(CLOCK_MONOTONIC, &ts1);

    sigaction(SIGINT, &old_sa_int, NULL);
    sigaction(SIGTERM, &old_sa_term, NULL);
    sigaction(SIGQUIT, &old_sa_quit, NULL);

    double dt = (double) (ts1.tv_sec - ts0.tv_sec) + (double) (ts1.tv_nsec - ts0.tv_nsec) * 1e-9;
    double fps = (dt > 0.0) ? ((double) n_streamed / dt) : 0.0;
    printf("[milkatmturb] Streamed %ld frames in %.3f s (%.1f fps)\n", n_streamed, dt, fps);
    fflush(stdout);

    if (CONF_FRESNEL_PROPAGATION == 1)
    {
        atmturb_fresnel_ctx_free(&fctx);
        atmturb_fresnel_plan_free(&fplan);
    }
    else if (CONF_FRESNEL_PROPAGATION == 2)
    {
        atmturb_rytov_ctx_free(&rctx);
        atmturb_rytov_plan_free(&rplan);
    }
    free(scratch);
    atmturb_wfs_teardown_sim(&prof, &geom, &rsim);
    return 0;
}
