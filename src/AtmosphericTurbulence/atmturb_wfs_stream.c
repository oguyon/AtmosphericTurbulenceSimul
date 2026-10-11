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
#include "atmturb_cuda.h"
#include "atmturb_cuda_fresnel.h"
#include "atmturb_cuda_rytov.h"
#include "atmturb_wfs_render_cuda.h"
#include "atmturb_wfs_stream.h"
#include "ImageStreamIO/ImageStreamIO.h"
#include "fps.h"
#include "processtools.h"

static volatile sig_atomic_t g_stream_stop = 0;
static PROCESSINFO          *g_stream_processinfo = NULL;

/**
 * atmturb_wfs_stream_set_processinfo - Register active FPS processinfo for rate control
 * @procinfo: Active processinfo handle, or NULL if processinfo is disabled.
 */
void atmturb_wfs_stream_set_processinfo(
    PROCESSINFO *procinfo)
{
    g_stream_processinfo = procinfo;
}

/**
 * atmturb_wfs_stream_get_processinfo - Retrieve active FPS processinfo handle
 *
 * Return: Active processinfo handle or NULL.
 */
PROCESSINFO *atmturb_wfs_stream_get_processinfo(void)
{
    return g_stream_processinfo;
}

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
    uint32_t req_cbsize = (CONF_STREAM_CBSIZE > 0) ? (uint32_t) CONF_STREAM_CBSIZE : 0;
    imageID  id         = image_ID(name);
    if (id >= 0 && dcimg[id].used == 1)
    {
        if (dcimg[id].md[0].naxis == 2 &&
            (long) dcimg[id].md[0].size[0] == pup_size &&
            (long) dcimg[id].md[0].size[1] == pup_size &&
            dcimg[id].md[0].datatype == _DATATYPE_FLOAT &&
            dcimg[id].md[0].CBsize == req_cbsize)
        {
            return id;
        }
        delete_image_ID(name);
        id = -1;
    }

    uint32_t sz[2] = {(uint32_t) pup_size, (uint32_t) pup_size};
    create_image_ID(name, 2, sz, _DATATYPE_FLOAT, 1, 10, req_cbsize, &id);
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
 * @cuda_rstream: GPU Rytov streaming context pointer.
 * @cuda_gstream: GPU geometric streaming context pointer.
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
    atmturb_cuda_fresnel_stream_t *cuda_fstream;
    atmturb_cuda_rytov_stream_t  *cuda_rstream;
    atmturb_cuda_geom_stream_t   *cuda_gstream;
    int                           is_cuda;
    float                        *dma_pha[2];
    float                        *dma_amp[2];
    float                        *dma_spha[2];
    float                        *dma_samp[2];
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
    float *pha_dst = (st->is_cuda && st->id_pha >= 0) ? dcimg[st->id_pha].array.F : st->pha;
    float *amp_dst = (st->id_amp >= 0)
                         ? ((st->is_cuda) ? dcimg[st->id_amp].array.F : st->amp)
                         : NULL;
    float *spha_dst = (st->id_spha >= 0)
                          ? ((st->is_cuda) ? dcimg[st->id_spha].array.F : st->spha)
                          : NULL;
    float *samp_dst = (st->id_samp >= 0)
                          ? ((st->is_cuda) ? dcimg[st->id_samp].array.F : st->samp)
                          : NULL;

    if (CONF_FRESNEL_PROPAGATION == 1 && st->cuda_fstream != NULL)
    {
        if (atmturb_cuda_fresnel_stream_render_step(st->cuda_fstream, rsim, geom, t,
                                                    time_step_s, pha_dst, amp_dst,
                                                    spha_dst, samp_dst) == 0)
        {
            return;
        }
    }

    if (CONF_FRESNEL_PROPAGATION == 1 && st->fplan != NULL && st->fctx != NULL)
    {
        atmturb_fresnel_render_step(st->fctx, st->fplan, rsim, geom, t, time_step_s,
                                    master_size, pup_size, st->pha, amp_dst,
                                    spha_dst, samp_dst);
        return;
    }

    if (CONF_FRESNEL_PROPAGATION == 2 && st->cuda_rstream != NULL)
    {
        if (atmturb_cuda_rytov_stream_render_step(st->cuda_rstream, rsim, geom, t,
                                                 time_step_s, pha_dst, amp_dst,
                                                 spha_dst, samp_dst) == 0)
        {
            return;
        }
    }

    if (CONF_FRESNEL_PROPAGATION == 2 && st->rplan != NULL && st->rctx != NULL)
    {
        atmturb_rytov_render_step(st->rctx, st->rplan, rsim, geom, t, time_step_s,
                                  master_size, pup_size, st->pha, amp_dst,
                                  spha_dst, samp_dst);
        return;
    }

    if (CONF_FRESNEL_PROPAGATION == 0 && st->cuda_gstream != NULL)
    {
        if (atmturb_wfs_cuda_geom_stream_render_step(
                st->cuda_gstream, t, pha_dst, amp_dst, spha_dst, samp_dst) == 0)
        {
            return;
        }
    }

    long frame_pixels = pup_size * pup_size;
    atmturb_init_phase_amp(st->pha, amp_dst, frame_pixels);
    if (spha_dst != NULL || samp_dst != NULL)
    {
        atmturb_init_phase_amp(spha_dst, samp_dst, frame_pixels);
    }

    for (int k = 0; k < geom->nlayers; k++)
    {
        atmturb_wfs_render_layer(rsim, geom, k, t, time_step_s, master_size,
                                 pup_size, st->pha, spha_dst);
    }
}

/**
 * atmturb_wfs_stream_notify_channel - Post semaphores on an updated stream channel
 * @id: Image ID to notify.
 * @ts: Frame timestamp.
 */
static inline void atmturb_wfs_stream_notify_channel(
    imageID                id,
    const struct timespec *ts)
{
    if (id >= 0)
    {
        dcimg[id].md[0].write = 1;
        atmturb_wfs_stream_post(id, ts);
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
    if (st->is_cuda)
    {
        atmturb_wfs_stream_notify_channel(st->id_pha, ts);
        atmturb_wfs_stream_notify_channel(st->id_amp, ts);
        atmturb_wfs_stream_notify_channel(st->id_spha, ts);
        atmturb_wfs_stream_notify_channel(st->id_samp, ts);
        return;
    }

    size_t nbytes = (size_t) frame_pixels * sizeof(float);

    if (st->id_pha >= 0)
    {
        dcimg[st->id_pha].md[0].write = 1;
        atmturb_remove_piston_stream(dcimg[st->id_pha].array.F, st->pha, frame_pixels);
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
        dcimg[st->id_spha].md[0].write = 1;
        atmturb_remove_piston_stream(dcimg[st->id_spha].array.F, st->spha, frame_pixels);
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
 * atmturb_wfs_stream_render_wait - Wait for in-flight GPU frame completion
 * @st: Stream state container.
 * @buf_idx: Double-buffer slot index (0 or 1).
 *
 * Return: 0 on success, -1 on failure.
 */
static inline int atmturb_wfs_stream_render_wait(
    atmturb_stream_state_t *st,
    int                     buf_idx)
{
    if (CONF_FRESNEL_PROPAGATION == 1 && st->cuda_fstream != NULL)
    {
        return atmturb_cuda_fresnel_stream_render_wait(st->cuda_fstream, buf_idx);
    }
    if (CONF_FRESNEL_PROPAGATION == 2 && st->cuda_rstream != NULL)
    {
        return atmturb_cuda_rytov_stream_render_wait(st->cuda_rstream, buf_idx);
    }
    return 0;
}

/**
 * atmturb_wfs_stream_pipeline_publish - Copy double-buffer frame to SHM and post semaphores
 * @st: Stream state container.
 * @buf_idx: Double-buffer index (0 or 1).
 * @nbytes: Frame size in bytes.
 * @ts: Timestamp for frame.
 */
static inline void atmturb_wfs_stream_write_channel(
    imageID                id,
    const void            *src,
    size_t                 nbytes,
    const struct timespec *ts)
{
    if (id >= 0 && src != NULL)
    {
        memcpy(dcimg[id].array.F, src, nbytes);
        dcimg[id].md[0].write = 1;
        atmturb_wfs_stream_post(id, ts);
    }
}

static void atmturb_wfs_stream_pipeline_publish(
    const atmturb_stream_state_t *st,
    int                           buf_idx,
    size_t                        nbytes,
    const struct timespec        *ts)
{
    atmturb_wfs_stream_write_channel(st->id_pha, st->dma_pha[buf_idx], nbytes, ts);
    atmturb_wfs_stream_write_channel(st->id_amp, st->dma_amp[buf_idx], nbytes, ts);
    atmturb_wfs_stream_write_channel(st->id_spha, st->dma_spha[buf_idx], nbytes, ts);
    atmturb_wfs_stream_write_channel(st->id_samp, st->dma_samp[buf_idx], nbytes, ts);
}

/**
 * atmturb_wfs_stream_pipeline_launch_frame - Dispatch next frame to GPU asynchronously
 * @rsim: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @params: Observation parameters container.
 * @st: Stream state container.
 * @t: Frame index.
 * @buf_idx: Double-buffer slot index.
 *
 * Return: 0 on success, -1 on failure.
 */
static inline int atmturb_wfs_stream_pipeline_launch_frame(
    const atmturb_rolling_t    *rsim,
    const atmturb_geom_t       *geom,
    const atmturb_obs_params_t *params,
    atmturb_stream_state_t     *st,
    long                        t,
    int                         buf_idx)
{
    float *pha_dst  = st->dma_pha[buf_idx];
    float *amp_dst  = (st->id_amp >= 0) ? st->dma_amp[buf_idx] : NULL;
    float *spha_dst = (st->id_spha >= 0) ? st->dma_spha[buf_idx] : NULL;
    float *samp_dst = (st->id_samp >= 0) ? st->dma_samp[buf_idx] : NULL;

    if (CONF_FRESNEL_PROPAGATION == 1 && st->cuda_fstream != NULL)
    {
        return atmturb_cuda_fresnel_stream_render_launch(
            st->cuda_fstream, rsim, geom, t, params->time_step_s,
            pha_dst, amp_dst, spha_dst, samp_dst, buf_idx);
    }
    if (CONF_FRESNEL_PROPAGATION == 2 && st->cuda_rstream != NULL)
    {
        return atmturb_cuda_rytov_stream_render_launch(
            st->cuda_rstream, rsim, geom, t, params->time_step_s,
            pha_dst, amp_dst, spha_dst, samp_dst, buf_idx);
    }
    return -1;
}

/**
 * atmturb_wfs_stream_loop_pipelined - Asynchronous double-buffered GPU streaming loop
 * @rsim: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @params: Observation parameters container.
 * @st: Stream state container.
 * @pup_size: Linear dimension of square pupil.
 * @max_frames: Maximum frames to render (0 = infinite).
 *
 * Return: Total frames streamed.
 */
static long atmturb_wfs_stream_loop_pipelined(
    const atmturb_rolling_t    *rsim,
    const atmturb_geom_t       *geom,
    const atmturb_obs_params_t *params,
    atmturb_stream_state_t     *st,
    long                        pup_size,
    long                        max_frames)
{
    size_t nbytes = (size_t) (pup_size * pup_size) * sizeof(float);
    PROCESSINFO *pinfo = g_stream_processinfo;
    long t = 0;
    int in_flight = 1;

    atmturb_wfs_stream_pipeline_launch_frame(rsim, geom, params, st, 0, 0);
    if (g_stream_stop == 0 && (max_frames <= 0 || 1 < max_frames))
    {
        atmturb_wfs_stream_pipeline_launch_frame(rsim, geom, params, st, 1, 1);
        in_flight++;
    }

    for (t = 0; ; t++)
    {
        int cur_buf = (int) (t & 1);

        atmturb_wfs_stream_render_wait(st, cur_buf);
        in_flight--;

        struct timespec ts_now;
        clock_gettime(CLOCK_REALTIME, &ts_now);
        atmturb_wfs_stream_pipeline_publish(st, cur_buf, nbytes, &ts_now);

        long next_t = t + 2;
        if (g_stream_stop == 0 && (max_frames <= 0 || next_t < max_frames))
        {
            atmturb_wfs_stream_pipeline_launch_frame(
                rsim, geom, params, st, next_t, cur_buf);
            in_flight++;
        }

        if (pinfo != NULL && t > 0)
        {
            processinfo_exec_end(pinfo);
            if (processinfo_loopstep(pinfo) == 0 || g_stream_stop != 0)
            {
                t++;
                break;
            }
            if (dcfpsptr != NULL && dcfpsptr->md != NULL &&
                (dcfpsptr->md->signal & FUNCTION_PARAMETER_STRUCT_SIGNAL_UPDATE))
            {
                fps_to_processinfo(dcfpsptr, pinfo);
                dcfpsptr->md->signal &= ~FUNCTION_PARAMETER_STRUCT_SIGNAL_UPDATE;
            }
            processinfo_waitoninputstream(pinfo);
            processinfo_exec_start(pinfo);
        }

        if (g_stream_stop != 0 || (max_frames > 0 && (t + 1) >= max_frames))
        {
            t++;
            break;
        }
    }

    for (int k = 0; k < in_flight; k++)
    {
        atmturb_wfs_stream_render_wait(st, (int) ((t + k) & 1));
    }

    if (pinfo != NULL)
    {
        pinfo->CTRLval = 3;
    }
    return t;
}

/**
 * atmturb_wfs_stream_loop - Inner frame-by-frame streaming loop governed by processinfo
 * @rsim: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @params: Observation parameters container.
 * @st: Stream state container.
 * @pup_size: Output pupil linear dimension.
 * @max_frames: Maximum frames to render (0 = infinite).
 *
 * Return: Total frames streamed.
 */
static long atmturb_wfs_stream_loop(
    const atmturb_rolling_t    *rsim,
    const atmturb_geom_t       *geom,
    const atmturb_obs_params_t *params,
    atmturb_stream_state_t     *st,
    long                        pup_size,
    long                        max_frames)
{
    if (st->dma_pha[0] != NULL && (st->cuda_fstream != NULL || st->cuda_rstream != NULL))
    {
        return atmturb_wfs_stream_loop_pipelined(
            rsim, geom, params, st, pup_size, max_frames);
    }

    long frame_pixels = pup_size * pup_size;
    long t = 0;
    PROCESSINFO *pinfo = g_stream_processinfo;

    for (t = 0; ; t++)
    {
        if (g_stream_stop != 0 || (max_frames > 0 && t >= max_frames))
        {
            break;
        }

        if (pinfo != NULL && t > 0)
        {
            processinfo_exec_end(pinfo);
            if (processinfo_loopstep(pinfo) == 0)
            {
                break;
            }
            if (dcfpsptr != NULL && dcfpsptr->md != NULL &&
                (dcfpsptr->md->signal & FUNCTION_PARAMETER_STRUCT_SIGNAL_UPDATE))
            {
                fps_to_processinfo(dcfpsptr, pinfo);
                dcfpsptr->md->signal &= ~FUNCTION_PARAMETER_STRUCT_SIGNAL_UPDATE;
            }
            processinfo_waitoninputstream(pinfo);
            if (pinfo->triggerstatus == PROCESSINFO_TRIGGERSTATUS_TIMEDOUT &&
                pinfo->triggermode == PROCESSINFO_TRIGGERMODE_SEMAPHORE)
            {
                continue;
            }
            processinfo_exec_start(pinfo);
        }

        atmturb_wfs_stream_render_step(rsim, geom, st, t, params->time_step_s,
                                       CONF_MASTER_SIZE, pup_size);

        struct timespec ts_now;
        clock_gettime(CLOCK_REALTIME, &ts_now);

        atmturb_wfs_stream_publish(st, frame_pixels, &ts_now);
    }

    if (pinfo != NULL)
    {
        pinfo->CTRLval = 3;
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
 * struct atmturb_stream_engines_t - Propagation plans and execution contexts
 * @fplan: Fresnel propagation plan.
 * @fctx: Fresnel execution context.
 * @rplan: Rytov propagation plan.
 * @rctx: Rytov execution context.
 */
typedef struct
{
    atmturb_fresnel_plan_t fplan;
    atmturb_fresnel_ctx_t  fctx;
    atmturb_rytov_plan_t   rplan;
    atmturb_rytov_ctx_t    rctx;
} atmturb_stream_engines_t;

/**
 * atmturb_wfs_stream_init_engine - Initialize GPU or CPU propagation engine
 * @st: Stream state container.
 * @rsim: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @prof: Active atmospheric profile.
 * @params: Observation parameters container.
 * @pup_size: Linear dimension of pupil.
 * @eng: CPU propagation engine container.
 */
static void atmturb_wfs_stream_init_engine(
    atmturb_stream_state_t          *st,
    const atmturb_rolling_t         *rsim,
    const atmturb_geom_t            *geom,
    const atmturb_profile_t         *prof,
    const atmturb_obs_params_t      *params,
    long                             pup_size,
    atmturb_stream_engines_t        *eng)
{
    if (CONF_FRESNEL_PROPAGATION == 0)
    {
        if (atmturb_simd_is_gpu())
        {
            st->cuda_gstream = atmturb_wfs_cuda_geom_stream_init(
                rsim, geom, CONF_MASTER_SIZE, pup_size);
        }
    }
    else if (CONF_FRESNEL_PROPAGATION == 1)
    {
        if (atmturb_simd_is_gpu())
        {
            st->cuda_fstream = atmturb_cuda_fresnel_stream_init(
                rsim, geom, prof, params, CONF_MASTER_SIZE, pup_size);
        }
        if (st->cuda_fstream == NULL)
        {
            double z_bin = (double) CONF_FRESNEL_PROPAGATION_BIN;
            atmturb_fresnel_plan_init(&eng->fplan, prof, geom, pup_size, params->pupil_scale_m,
                                      params->lambda_ref_m, params->lambda_s_m, z_bin);
            atmturb_fresnel_ctx_init(&eng->fctx, pup_size);
            st->fplan = &eng->fplan;
            st->fctx  = &eng->fctx;
        }
    }
    else if (CONF_FRESNEL_PROPAGATION == 2)
    {
        if (atmturb_simd_is_gpu())
        {
            st->cuda_rstream = atmturb_cuda_rytov_stream_init(
                rsim, geom, prof, params, CONF_MASTER_SIZE, pup_size);
        }
        if (st->cuda_rstream == NULL)
        {
            double z_bin = (double) CONF_FRESNEL_PROPAGATION_BIN;
            atmturb_rytov_plan_init(&eng->rplan, prof, geom, pup_size, params->pupil_scale_m,
                                    params->lambda_ref_m, params->lambda_s_m, z_bin,
                                    CONF_FRESNEL_RYTOV_SEC_EXACT);
            atmturb_rytov_ctx_init(&eng->rctx, eng->rplan.pad_size);
            st->rplan = &eng->rplan;
            st->rctx  = &eng->rctx;
        }
    }
}

/**
 * atmturb_wfs_stream_free_engine - Release allocated GPU and CPU propagation engines
 * @st: Stream state container.
 * @eng: CPU propagation engine container.
 */
static void atmturb_wfs_stream_free_engine(
    atmturb_stream_state_t   *st,
    atmturb_stream_engines_t *eng)
{
    if (CONF_FRESNEL_PROPAGATION == 0)
    {
        if (st->cuda_gstream != NULL)
        {
            atmturb_wfs_cuda_geom_stream_free(st->cuda_gstream);
            st->cuda_gstream = NULL;
        }
    }
    else if (CONF_FRESNEL_PROPAGATION == 1)
    {
        if (st->cuda_fstream != NULL)
        {
            atmturb_cuda_fresnel_stream_free(st->cuda_fstream);
            st->cuda_fstream = NULL;
        }
        if (st->fctx != NULL)
        {
            atmturb_fresnel_ctx_free(&eng->fctx);
            atmturb_fresnel_plan_free(&eng->fplan);
        }
    }
    else if (CONF_FRESNEL_PROPAGATION == 2)
    {
        if (st->cuda_rstream != NULL)
        {
            atmturb_cuda_rytov_stream_free(st->cuda_rstream);
            st->cuda_rstream = NULL;
        }
        if (st->rctx != NULL)
        {
            atmturb_rytov_ctx_free(&eng->rctx);
            atmturb_rytov_plan_free(&eng->rplan);
        }
    }

    atmturb_cuda_cleanup();
    atmturb_cuda_rytov_cleanup();
}

/**
 * atmturb_wfs_dma_alloc_pair - Allocate a pair of pinned host buffers
 * @p: Array of 2 pointers to allocate.
 * @nbytes: Size in bytes of each buffer.
 */
static void atmturb_wfs_dma_alloc_pair(
    float **p,
    size_t  nbytes)
{
    atmturb_cuda_host_alloc((void **) &p[0], nbytes);
    atmturb_cuda_host_alloc((void **) &p[1], nbytes);
}

/**
 * atmturb_wfs_dma_free_pair - Free a pair of pinned host buffers
 * @p: Array of 2 pointers to release.
 */
static void atmturb_wfs_dma_free_pair(
    float **p)
{
    for (int b = 0; b < 2; b++)
    {
        if (p[b] != NULL)
        {
            atmturb_cuda_host_free(p[b]);
            p[b] = NULL;
        }
    }
}

/**
 * atmturb_wfs_stream_setup_dma - Allocate pinned host buffers or register SHM for GPU DMA
 * @st: Stream state container.
 * @pup_size: Linear dimension of square pupil.
 */
static inline void atmturb_wfs_host_reg_channel(
    imageID id,
    size_t  nbytes)
{
    if (id >= 0)
    {
        atmturb_cuda_host_register(dcimg[id].array.F, nbytes);
    }
}

static inline void atmturb_wfs_host_unreg_channel(imageID id)
{
    if (id >= 0)
    {
        atmturb_cuda_host_unregister(dcimg[id].array.F);
    }
}

static void atmturb_wfs_stream_setup_dma(
    atmturb_stream_state_t *st,
    long                    pup_size)
{
    if (!st->is_cuda)
    {
        return;
    }

    size_t nbytes = (size_t) (pup_size * pup_size) * sizeof(float);
    if (st->cuda_gstream != NULL)
    {
        atmturb_wfs_host_reg_channel(st->id_pha, nbytes);
        atmturb_wfs_host_reg_channel(st->id_amp, nbytes);
        atmturb_wfs_host_reg_channel(st->id_spha, nbytes);
        atmturb_wfs_host_reg_channel(st->id_samp, nbytes);
        return;
    }

    atmturb_wfs_dma_alloc_pair(st->dma_pha, nbytes);
    if (st->id_amp >= 0)
    {
        atmturb_wfs_dma_alloc_pair(st->dma_amp, nbytes);
    }
    if (st->id_spha >= 0)
    {
        atmturb_wfs_dma_alloc_pair(st->dma_spha, nbytes);
    }
    if (st->id_samp >= 0)
    {
        atmturb_wfs_dma_alloc_pair(st->dma_samp, nbytes);
    }
}

/**
 * atmturb_wfs_stream_cleanup_dma - Release pinned host buffers and unregister SHM
 * @st: Stream state container.
 */
static void atmturb_wfs_stream_cleanup_dma(
    atmturb_stream_state_t *st)
{
    if (!st->is_cuda)
    {
        return;
    }

    if (st->cuda_gstream != NULL)
    {
        atmturb_wfs_host_unreg_channel(st->id_pha);
        atmturb_wfs_host_unreg_channel(st->id_amp);
        atmturb_wfs_host_unreg_channel(st->id_spha);
        atmturb_wfs_host_unreg_channel(st->id_samp);
        return;
    }

    atmturb_wfs_dma_free_pair(st->dma_pha);
    atmturb_wfs_dma_free_pair(st->dma_amp);
    atmturb_wfs_dma_free_pair(st->dma_spha);
    atmturb_wfs_dma_free_pair(st->dma_samp);
}

/**
 * make_AtmosphericTurbulence_wavefront_stream - Stream 2D wavefront frames to SHM
 * @slambdaum: Secondary observing wavelength in um.
 * @WFprecision: Precision mode flag (0=single, 1=double).
 * @stream_mode: Streaming mode (1=continuous stream, 2=finite stream).
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

    atmturb_stream_engines_t eng;
    memset(&eng, 0, sizeof(eng));
    atmturb_wfs_stream_init_engine(&st, &rsim, &geom, &prof, &params, pup_size, &eng);

    st.is_cuda = (st.cuda_rstream != NULL || st.cuda_gstream != NULL || st.cuda_fstream != NULL);
    atmturb_wfs_stream_setup_dma(&st, pup_size);

    g_stream_stop = 0;
    struct sigaction sa, old_sa_int, old_sa_term, old_sa_quit;
    memset(&sa, 0, sizeof(sa));
    sa.sa_handler = atmturb_wfs_stream_sig_handler;
    sigaction(SIGINT, &sa, &old_sa_int);
    sigaction(SIGTERM, &sa, &old_sa_term);
    sigaction(SIGQUIT, &sa, &old_sa_quit);

    const char *isa_str = "AVX2";
    if (st.cuda_rstream != NULL || st.cuda_gstream != NULL || st.cuda_fstream != NULL)
    {
        isa_str = "CUDA GPU";
    }
    else
    {
        const char *active = atmturb_simd_active_isa();
        isa_str = (strcmp(active, "CUDA GPU") == 0) ? "AVX2" : active;
    }
    printf("[milkatmturb] Streaming 2D wavefront frames (%ldx%ld pix) to SHM '%s'%s [%s]\n",
           pup_size, pup_size, pha_name, (st.id_amp >= 0) ? " and amplitude" : "",
           isa_str);
    fflush(stdout);

    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    long n_streamed = atmturb_wfs_stream_loop(&rsim, &geom, &params, &st, pup_size,
                                             max_frames);
    clock_gettime(CLOCK_MONOTONIC, &ts1);

    sigaction(SIGINT, &old_sa_int, NULL);
    sigaction(SIGTERM, &old_sa_term, NULL);
    sigaction(SIGQUIT, &old_sa_quit, NULL);

    double dt = (double) (ts1.tv_sec - ts0.tv_sec) + (double) (ts1.tv_nsec - ts0.tv_nsec) * 1e-9;
    double fps = (dt > 0.0) ? ((double) n_streamed / dt) : 0.0;
    printf("[milkatmturb] Streamed %ld frames in %.3f s (%.1f fps)\n", n_streamed, dt, fps);
    fflush(stdout);

    atmturb_wfs_stream_cleanup_dma(&st);

    atmturb_wfs_stream_free_engine(&st, &eng);
    free(scratch);
    atmturb_wfs_teardown_sim(&prof, &geom, &rsim);
    return 0;
}
