// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_psf_sim.c
 * @brief   End-to-end adaptive optics closed-loop simulation and PSF synthesis
 */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fitsio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "AtmosphericTurbulence.h"
#include "atmturb_cuda_psf.h"
#include "atmturb_psf_engine.h"
#include "atmturb_simd.h"
#include "atmturb_types.h"

/**
 * struct atmturb_ao_controller_t - Closed-loop AO state and latency ring buffer
 * @delay: Latency delay in frames (0 to 3).
 * @depth: Ring buffer depth (delay + 1).
 * @head: Current write pointer in ring buffer.
 * @ntot: Total pupil pixels (pup_size * pup_size).
 * @dm_opd: Current deformable mirror shape in meters.
 * @pid_int: Accumulated integral error for PID mode.
 * @prev_err: Previous error buffer for derivative calculation.
 * @delay_ring: Ring buffer slots holding delayed error frames.
 */
typedef struct
{
    int    delay;
    int    depth;
    int    head;
    long   ntot;
    float *dm_opd;
    float *pid_int;
    float *prev_err;
    float *delay_ring[4];
} atmturb_ao_controller_t;

/**
 * struct atmturb_ao_wf_t - Wavefront data container
 * @data: Float pointer to 3D cube [pup_size * pup_size * nbframes].
 * @pup_size: Linear dimension of pupil.
 * @nbframes: Number of frames in cube.
 * @is_allocated: Flag if data was allocated via malloc.
 * @is_shm: Flag if data was connected via ImageStreamIO.
 * @shm_img: ImageStreamIO image structure.
 */
typedef struct
{
    float *data;
    long   pup_size;
    long   nbframes;
    int    is_allocated;
    int    is_shm;
    IMAGE  shm_img;
} atmturb_ao_wf_t;

/**
 * atmturb_ao_init_params - Populate default parameters for AO simulation
 * @params: Structure to initialize.
 */
void atmturb_ao_init_params(
    atmturb_ao_params_t *params)
{
    if (params == NULL)
    {
        return;
    }
    memset(params, 0, sizeof(*params));
    params->in_wfname     = "outarraypha";
    params->in_ampname    = NULL;
    params->out_psfname   = "PSFcumul";
    params->out_fitsname  = NULL;
    params->pup_size      = 0;
    params->nbframes      = 0;
    params->lambda_ref_m  = 0.50e-6;
    params->lambda_sci_m  = 1.65e-6;
    params->tel_diam_m    = 8.0;
    params->pupil_scale_m = 0.0;
    params->loop_mode     = 1;
    params->gain          = 0.5;
    params->leak          = 0.001;
    params->loop_delay    = 1;
    params->Kp            = 0.5;
    params->Ki            = 0.0;
    params->Kd            = 0.0;
    params->time_step_s   = 0.001;
    params->save_psfcube  = 0;
}

/**
 * atmturb_ao_init_pupil_mask - Construct centered circular telescope pupil aperture
 * @pup_size: Linear dimension of pupil in pixels.
 * @tel_diam_m: Primary aperture diameter in meters.
 * @pupil_scale_m: Pixel scale in meters/pixel (<= 0 for pupil fitting grid).
 * @mask: Output buffer of size pup_size * pup_size.
 *
 * Return: Total number of active pupil pixels.
 */
static double atmturb_ao_init_pupil_mask(
    long    pup_size,
    double  tel_diam_m,
    double  pupil_scale_m,
    float  *mask)
{
    double r_pix = (pupil_scale_m > 0.0 && tel_diam_m > 0.0)
                       ? (0.5 * tel_diam_m / pupil_scale_m)
                       : (0.5 * (double) pup_size);
    if (r_pix > 0.5 * (double) pup_size)
    {
        r_pix = 0.5 * (double) pup_size;
    }
    double r2   = r_pix * r_pix;
    double cx   = 0.5 * ((double) pup_size - 1.0);
    double cy   = 0.5 * ((double) pup_size - 1.0);
    double npup = 0.0;

    for (long j = 0; j < pup_size; j++)
    {
        double dy = (double) j - cy;
        for (long i = 0; i < pup_size; i++)
        {
            double dx  = (double) i - cx;
            long   idx = j * pup_size + i;
            if (dx * dx + dy * dy <= r2)
            {
                mask[idx] = 1.0f;
                npup += 1.0;
            }
            else
            {
                mask[idx] = 0.0f;
            }
        }
    }

    return (npup > 0.0) ? npup : 1.0;
}

/**
 * atmturb_ao_controller_init - Allocate DM state and latency ring buffer
 * @ctrl: Controller state structure to initialize.
 * @ntot: Number of pupil pixels.
 * @delay: Latency delay in frames (0 to 3).
 *
 * Return: 0 on success, -1 on allocation failure.
 */
static int atmturb_ao_controller_init(
    atmturb_ao_controller_t *ctrl,
    long                     ntot,
    int                      delay)
{
    memset(ctrl, 0, sizeof(*ctrl));
    ctrl->delay = (delay < 0) ? 0 : ((delay > 3) ? 3 : delay);
    ctrl->depth = ctrl->delay + 1;
    ctrl->head  = 0;
    ctrl->ntot  = ntot;

    ctrl->dm_opd   = (float *) calloc((size_t) ntot, sizeof(float));
    ctrl->pid_int  = (float *) calloc((size_t) ntot, sizeof(float));
    ctrl->prev_err = (float *) calloc((size_t) ntot, sizeof(float));
    if (ctrl->dm_opd == NULL || ctrl->pid_int == NULL || ctrl->prev_err == NULL)
    {
        return -1;
    }

    for (int d = 0; d < ctrl->depth; d++)
    {
        ctrl->delay_ring[d] = (float *) calloc((size_t) ntot, sizeof(float));
        if (ctrl->delay_ring[d] == NULL)
        {
            return -1;
        }
    }
    return 0;
}

/**
 * atmturb_ao_controller_free - Deallocate controller memory
 * @ctrl: Controller state pointer.
 */
static void atmturb_ao_controller_free(
    atmturb_ao_controller_t *ctrl)
{
    free(ctrl->dm_opd);
    free(ctrl->pid_int);
    free(ctrl->prev_err);
    for (int d = 0; d < 4; d++)
    {
        free(ctrl->delay_ring[d]);
    }
    memset(ctrl, 0, sizeof(*ctrl));
}

/**
 * atmturb_ao_controller_step - Update DM and compute science residual for current frame
 * @ctrl: Active controller state.
 * @params: Simulation parameters.
 * @wf_slice: Input WFS phase slice in reference radians [ntot], or NULL for zero input.
 * @swf_slice: Optional science phase slice in science radians [ntot], or NULL to scale WFS OPD.
 * @pupil_mask: Circular pupil mask [ntot].
 * @sci_res_out: Output science residual (science radians if @swf_slice != NULL, else OPD [m]).
 */
static void atmturb_ao_controller_step(
    atmturb_ao_controller_t   *ctrl,
    const atmturb_ao_params_t *params,
    const float               *wf_slice,
    const float               *swf_slice,
    const float               *pupil_mask,
    float                     *sci_res_out)
{
    long               ntot      = ctrl->ntot;
    float              opd_scale = (float) (params->lambda_ref_m / (2.0 * M_PI));
    float              k_sci     = (float) (2.0 * M_PI / params->lambda_sci_m);
    float *restrict    slot      = ctrl->delay_ring[ctrl->head];
    const float *const dm_r      = ctrl->dm_opd;

    #pragma omp simd
    for (long i = 0; i < ntot; i++)
    {
        float wf_v     = (wf_slice != NULL) ? wf_slice[i] : 0.0f;
        float res      = wf_v * opd_scale - dm_r[i];
        slot[i]        = res * pupil_mask[i];
        sci_res_out[i] = (swf_slice != NULL) ? (swf_slice[i] - dm_r[i] * k_sci) : res;
    }

    int          read_idx    = (ctrl->head - ctrl->delay + ctrl->depth) % ctrl->depth;
    const float *delayed_err = ctrl->delay_ring[read_idx];

    if (params->loop_mode == 1)
    {
        float           keep = 1.0f - (float) params->leak;
        float           gain = (float) params->gain;
        float *restrict dm   = ctrl->dm_opd;
        #pragma omp simd
        for (long i = 0; i < ntot; i++)
        {
            dm[i] = keep * dm[i] + gain * delayed_err[i];
        }
    }
    else if (params->loop_mode == 2)
    {
        float           kp    = (float) params->Kp;
        float           ki    = (float) params->Ki;
        float           dt    = (params->time_step_s > 0.0) ? (float) params->time_step_s : 0.001f;
        float           kd_dt = (float) params->Kd / dt;
        float           keep  = 1.0f - (float) params->leak;
        float *restrict dm    = ctrl->dm_opd;
        float *restrict pint  = ctrl->pid_int;
        float *restrict perr  = ctrl->prev_err;
        #pragma omp simd
        for (long i = 0; i < ntot; i++)
        {
            float e   = delayed_err[i];
            float der = (e - perr[i]) * kd_dt;
            float pi  = keep * pint[i] + e * dt;
            pint[i]   = pi;
            dm[i]     = kp * e + ki * pi + der;
            perr[i]   = e;
        }
    }

    ctrl->head = (ctrl->head + 1) % ctrl->depth;
}

/**
 * atmturb_ao_wf_load_fits - Load 2D/3D single-precision wavefront cube from FITS file
 */
static int atmturb_ao_wf_load_fits(
    const char      *target,
    atmturb_ao_wf_t *wf)
{
    fitsfile *fptr   = NULL;
    int       status = 0;
    fits_open_file(&fptr, target, READONLY, &status);
    if (status != 0)
    {
        return -1;
    }

    int  naxis    = 0;
    long naxes[3] = {1, 1, 1};
    fits_get_img_dim(fptr, &naxis, &status);
    fits_get_img_size(fptr, (naxis > 3) ? 3 : naxis, naxes, &status);
    long   pup_size = naxes[0];
    long   nz       = (naxis >= 3) ? naxes[2] : 1;
    long   ntot     = pup_size * pup_size * nz;
    float *arr      = (float *) malloc(sizeof(float) * (size_t) ntot);
    if (arr == NULL)
    {
        fits_close_file(fptr, &status);
        return -1;
    }

    long fpixel[3] = {1, 1, 1};
    int  anynul    = 0;
    fits_read_pix(fptr, TFLOAT, fpixel, ntot, NULL, arr, &anynul, &status);
    fits_close_file(fptr, &status);
    if (status != 0)
    {
        free(arr);
        return -1;
    }

    wf->data         = arr;
    wf->pup_size     = pup_size;
    wf->nbframes     = nz;
    wf->is_allocated = 1;
    return 0;
}

/**
 * atmturb_ao_wf_load - Resolve wavefront input from memory stream or FITS file
 * @wfname: Input stream name or FITS file path.
 * @default_pup_size: Fallback pupil dimension if flat.
 * @wf: Wavefront structure to populate.
 *
 * Return: 0 on success, -1 on failure.
 */
static int atmturb_ao_wf_load(
    const char      *wfname,
    long             default_pup_size,
    atmturb_ao_wf_t *wf)
{
    memset(wf, 0, sizeof(*wf));
    if (wfname == NULL || wfname[0] == '\0' || strcmp(wfname, "none") == 0)
    {
        wf->pup_size = (default_pup_size > 0) ? default_pup_size : 128;
        wf->nbframes = 1;
        return 0;
    }

    imageID id = image_ID(wfname, dcimg, dcnimg);
    if (id >= 0 && dcimg[id].used == 1)
    {
        wf->data     = dcimg[id].array.F;
        wf->pup_size = dcimg[id].md[0].size[0];
        wf->nbframes = (dcimg[id].md[0].naxis >= 3) ? dcimg[id].md[0].size[2] : 1;
        return 0;
    }

    if (ImageStreamIO_openIm(&wf->shm_img, wfname) == IMAGESTREAMIO_SUCCESS)
    {
        wf->data     = wf->shm_img.array.F;
        wf->pup_size = wf->shm_img.md->size[0];
        wf->nbframes = (wf->shm_img.md->naxis >= 3) ? wf->shm_img.md->size[2] : 1;
        wf->is_shm   = 1;
        return 0;
    }

    char        tmpname[512];
    const char *target = NULL;
    if (access(wfname, R_OK) == 0)
    {
        target = wfname;
    }
    else
    {
        snprintf(tmpname, sizeof(tmpname), "%s.fits", wfname);
        if (access(tmpname, R_OK) == 0)
        {
            target = tmpname;
        }
    }

    if (target != NULL && atmturb_ao_wf_load_fits(target, wf) == 0)
    {
        return 0;
    }

    printf("ERROR: cannot find or load wavefront stream or file '%s'\n", wfname);
    return -1;
}

/**
 * atmturb_ao_wf_free - Release heap or shared-memory handles for loaded wavefront cube
 * @wf: Wavefront source handle initialized by atmturb_ao_wf_load().
 */
static void atmturb_ao_wf_free(
    atmturb_ao_wf_t *wf)
{
    if (wf->is_allocated && wf->data != NULL)
    {
        free(wf->data);
    }
    if (wf->is_shm)
    {
        ImageStreamIO_closeIm(&wf->shm_img);
    }
    memset(wf, 0, sizeof(*wf));
}

/**
 * atmturb_ao_save_fits - Write a 2D single-precision PSF array to a FITS file
 */
static int atmturb_ao_save_fits(
    const char  *filename,
    const float *data,
    long         nx,
    long         ny)
{
    fitsfile *fptr   = NULL;
    int       status = 0;
    remove(filename);
    fits_create_file(&fptr, filename, &status);
    if (status != 0)
    {
        return -1;
    }
    long naxes[2]  = {nx, ny};
    long fpixel[2] = {1, 1};
    fits_create_img(fptr, FLOAT_IMG, 2, naxes, &status);
    fits_write_pix(fptr, TFLOAT, fpixel, nx * ny, (void *) data, &status);
    fits_close_file(fptr, &status);
    return (status == 0) ? 0 : -1;
}

/**
 * atmturb_ao_opt_wf_load - Load and validate optional amplitude or chromatic phase cube/map
 */
static int atmturb_ao_opt_wf_load(
    const char      *name,
    long            *io_pup_size,
    long            *io_req_frames,
    atmturb_ao_wf_t *out,
    long            *out_stride)
{
    memset(out, 0, sizeof(*out));
    *out_stride = 0;
    if (name == NULL || name[0] == '\0' || strcmp(name, "none") == 0)
    {
        return 0;
    }
    if (atmturb_ao_wf_load(name, *io_pup_size, out) != 0 || out->data == NULL)
    {
        atmturb_ao_wf_free(out);
        return -1;
    }
    if (*io_pup_size <= 0)
    {
        *io_pup_size = out->pup_size;
    }
    if (*io_req_frames <= 0)
    {
        *io_req_frames = out->nbframes;
    }
    if (out->pup_size != *io_pup_size || (out->nbframes != 1 && out->nbframes < *io_req_frames))
    {
        printf("ERROR: channel '%s' size (%ldx%ldx%ld) incompatible with (%ldx%ldx%ld)\n",
               name, out->pup_size, out->pup_size, out->nbframes,
               *io_pup_size, *io_pup_size, *io_req_frames);
        atmturb_ao_wf_free(out);
        return -1;
    }
    *out_stride = (out->nbframes == 1) ? 0 : (*io_pup_size * *io_pup_size);
    return 0;
}

/**
 * atmturb_ao_run_open_loop_mt - Parallel open-loop PSF accumulation across frames
 */
static void atmturb_ao_run_open_loop_mt(
    const atmturb_ao_psf_engine_t *base_eng,
    const atmturb_cuda_ao_req_t   *req,
    float                          phase_scale,
    double                         dl_peak)
{
    long         n_pix      = base_eng->n_psf * base_eng->n_psf;
    long         ntot       = req->pup_size * req->pup_size;
    const float *pha_src    = (req->swf_data != NULL) ? req->swf_data : req->wf_data;
    long         pha_stride = (req->swf_data != NULL) ? req->swf_stride : ntot;
    float        p_scale    = (req->swf_data != NULL) ? 1.0f : phase_scale;

#ifdef _OPENMP
    int nth = omp_get_max_threads();
    if (nth > 16)
    {
        nth = 16;
    }
    if (nth > 1 && req->nbframes >= 8 && req->psf_out != NULL)
    {
        #pragma omp parallel num_threads(nth)
        {
            atmturb_ao_psf_engine_t teng = *base_eng;
            teng.row_in = teng.row_out = teng.col_in = teng.col_out = NULL;
            float *lpsf = (float *) calloc((size_t) n_pix, sizeof(float));
            if (atmturb_ao_psf_alloc_bufs(&teng) == 0 && lpsf != NULL)
            {
                #pragma omp for schedule(dynamic, 2)
                for (long t = 0; t < req->nbframes; t++)
                {
                    const float *sl  = (pha_src != NULL) ? &pha_src[t * pha_stride] : NULL;
                    const float *asl = (req->amp_data != NULL)
                                           ? &req->amp_data[t * req->amp_stride]
                                           : NULL;
                    double       sp  = atmturb_ao_render_psf_amp(&teng, sl, asl, req->pupil_mask,
                                                                 p_scale, lpsf, 1);
                    double       sr  = (dl_peak > 0.0) ? (sp / dl_peak) : 0.0;
                    if (t == 0)
                    {
                        req->s_endpoints[0] = sr;
                    }
                    if (t == req->nbframes - 1)
                    {
                        req->s_endpoints[1] = sr;
                    }
                }
                #pragma omp critical
                {
                    atmturb_add_float_array(req->psf_out, lpsf, n_pix);
                }
            }
            free(lpsf);
            atmturb_ao_psf_free_bufs(&teng);
        }
        return;
    }
#endif

    for (long t = 0; t < req->nbframes; t++)
    {
        const float *sl  = (pha_src != NULL) ? &pha_src[t * pha_stride] : NULL;
        const float *asl = (req->amp_data != NULL) ? &req->amp_data[t * req->amp_stride] : NULL;
        double       sp  = atmturb_ao_render_psf_amp(base_eng, sl, asl, req->pupil_mask,
                                                     p_scale, req->psf_out, 1);
        double       sr  = (dl_peak > 0.0) ? (sp / dl_peak) : 0.0;
        if (t == 0)
        {
            req->s_endpoints[0] = sr;
        }
        if (t == req->nbframes - 1)
        {
            req->s_endpoints[1] = sr;
        }
    }
}

/**
 * atmturb_ao_run_closed_loop_mt - Chunked controller rollout + parallel PSF accumulation
 */
static int atmturb_ao_run_closed_loop_mt(
    const atmturb_ao_psf_engine_t *eng,
    const atmturb_cuda_ao_req_t   *req,
    atmturb_ao_controller_t       *ctrl,
    float                          k_sci,
    double                         dl_peak)
{
    const atmturb_ao_params_t *p         = req->params;
    long                       ntot      = req->pup_size * req->pup_size;
    size_t                     frame_b   = (size_t) ntot * sizeof(float);
    long                       chunk_cap = (long) ((16ULL * 1024ULL * 1024ULL) / frame_b);
    if (chunk_cap < 16)
    {
        chunk_cap = 16;
    }
    if (chunk_cap > req->nbframes)
    {
        chunk_cap = req->nbframes;
    }

    float *opd_chunk = (float *) malloc((size_t) chunk_cap * frame_b);
    if (opd_chunk == NULL)
    {
        return -1;
    }

    float  chunk_scale = (req->swf_data != NULL) ? 1.0f : k_sci;
    double s_first     = 0.0;
    double s_last      = 0.0;

    for (long t0 = 0; t0 < req->nbframes; t0 += chunk_cap)
    {
        long n_cur = req->nbframes - t0;
        if (n_cur > chunk_cap)
        {
            n_cur = chunk_cap;
        }
        for (long c = 0; c < n_cur; c++)
        {
            long         t  = t0 + c;
            const float *sl = (req->wf_data != NULL) ? &req->wf_data[t * ntot] : NULL;
            const float *sw = (req->swf_data != NULL) ? &req->swf_data[t * req->swf_stride] : NULL;
            atmturb_ao_controller_step(ctrl, p, sl, sw, req->pupil_mask, &opd_chunk[c * ntot]);
        }

        double                chunk_ep[2] = {0.0, 0.0};
        atmturb_cuda_ao_req_t sub         = *req;
        sub.wf_data                       = opd_chunk;
        sub.swf_data                      = NULL;
        sub.swf_stride                    = 0;
        sub.amp_data                      = (req->amp_data != NULL)
                                                ? (req->amp_data + t0 * req->amp_stride)
                                                : NULL;
        sub.nbframes                      = n_cur;
        sub.s_endpoints                   = chunk_ep;
        atmturb_ao_run_open_loop_mt(eng, &sub, chunk_scale, dl_peak);

        if (t0 == 0)
        {
            s_first = chunk_ep[0];
        }
        if (t0 + n_cur == req->nbframes)
        {
            s_last = chunk_ep[1];
        }
    }

    req->s_endpoints[0] = s_first;
    req->s_endpoints[1] = s_last;
    free(opd_chunk);
    return 0;
}

/**
 * atmturb_ao_run_cpu - Execute open-loop or closed-loop AO PSF synthesis on CPU
 */
static int atmturb_ao_run_cpu(
    const atmturb_cuda_ao_req_t *req)
{
    const atmturb_ao_params_t *p        = req->params;
    int                        use_mt   = (req->do_accum && req->nbframes >= 8);
    int                        fftw_nth = 1;
#ifdef _OPENMP
    if (!use_mt && !omp_in_parallel())
    {
        fftw_nth = omp_get_max_threads();
        if (fftw_nth > 8)
        {
            fftw_nth = 8;
        }
    }
#endif

    atmturb_ao_psf_engine_t eng;
    if (atmturb_ao_psf_engine_init(&eng, req->pup_size, req->n_fft, fftw_nth) != 0)
    {
        return -1;
    }

    long                    ntot    = req->pup_size * req->pup_size;
    float                  *tmp_psf = (float *) malloc(sizeof(float) *
                                                       (size_t) (eng.n_psf * eng.n_psf));
    float                  *opd_res = NULL;
    atmturb_ao_controller_t ctrl;
    memset(&ctrl, 0, sizeof(ctrl));
    if (tmp_psf == NULL ||
        (p->loop_mode != 0 &&
         ((!use_mt && (opd_res = (float *) malloc(sizeof(float) * (size_t) ntot)) == NULL) ||
          atmturb_ao_controller_init(&ctrl, ntot, p->loop_delay) != 0)))
    {
        free(tmp_psf);
        free(opd_res);
        atmturb_ao_controller_free(&ctrl);
        atmturb_ao_psf_engine_free(&eng);
        return -1;
    }

    float  k_sci   = (float) (2.0 * M_PI / p->lambda_sci_m);
    double dl_peak = atmturb_ao_render_psf(&eng, NULL, req->pupil_mask, req->pup_size,
                                           k_sci, tmp_psf, 0);
    if (p->loop_mode == 0 && use_mt)
    {
        float open_scale = (float) (p->lambda_ref_m / p->lambda_sci_m);
        atmturb_ao_run_open_loop_mt(&eng, req, open_scale, dl_peak);
    }
    else if (p->loop_mode != 0 && use_mt &&
             atmturb_ao_run_closed_loop_mt(&eng, req, &ctrl, k_sci, dl_peak) == 0)
    {
        /* Completed via chunked multi-threaded rollout */
    }
    else
    {
        float open_scale = (float) (p->lambda_ref_m / p->lambda_sci_m);
        float cl_scale   = (req->swf_data != NULL) ? 1.0f : k_sci;
        for (long t = 0; t < req->nbframes; t++)
        {
            const float *slice = (req->wf_data != NULL) ? &req->wf_data[t * ntot] : NULL;
            const float *sw    = (req->swf_data != NULL)
                                     ? &req->swf_data[t * req->swf_stride]
                                     : NULL;
            const float *asl   = (req->amp_data != NULL)
                                     ? &req->amp_data[t * req->amp_stride]
                                     : NULL;
            double       sp    = 0.0;
            if (p->loop_mode == 0)
            {
                const float *psl = (req->swf_data != NULL) ? sw : slice;
                float        psc = (req->swf_data != NULL) ? 1.0f : open_scale;
                sp               = atmturb_ao_render_psf_amp(&eng, psl, asl, req->pupil_mask,
                                                             psc, req->psf_out, req->do_accum);
            }
            else
            {
                atmturb_ao_controller_step(&ctrl, p, slice, sw, req->pupil_mask, opd_res);
                sp = atmturb_ao_render_psf_amp(&eng, opd_res, asl, req->pupil_mask,
                                               cl_scale, req->psf_out, req->do_accum);
            }
            double sr = (dl_peak > 0.0) ? (sp / dl_peak) : 0.0;
            if (t == 0)
            {
                req->s_endpoints[0] = sr;
            }
            if (t == req->nbframes - 1)
            {
                req->s_endpoints[1] = sr;
            }
        }
    }

    *req->dl_peak_out = dl_peak;
    free(tmp_psf);
    free(opd_res);
    atmturb_ao_controller_free(&ctrl);
    atmturb_ao_psf_engine_free(&eng);
    return 0;
}

/**
 * atmturb_ao_sim_run - Execute closed-loop AO simulation and generate science PSF
 * @params: Simulation parameters.
 * @results: Output results container (optional, can be NULL).
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_ao_sim_run(
    const atmturb_ao_params_t *params,
    atmturb_ao_results_t      *results)
{
    atmturb_ao_params_t p;
    if (params != NULL)
    {
        p = *params;
    }
    else
    {
        atmturb_ao_init_params(&p);
    }

    atmturb_ao_wf_t wf, amp_wf, swf, samp_wf;
    memset(&amp_wf, 0, sizeof(amp_wf));
    memset(&swf, 0, sizeof(swf));
    memset(&samp_wf, 0, sizeof(samp_wf));
    if (atmturb_ao_wf_load(p.in_wfname, p.pup_size, &wf) != 0)
    {
        return -1;
    }

    long pup_size = (wf.data != NULL) ? wf.pup_size : p.pup_size;
    long nbframes = (wf.data != NULL)
                        ? ((p.nbframes > 0 && p.nbframes < wf.nbframes) ? p.nbframes : wf.nbframes)
                        : p.nbframes;
    long amp_str  = 0;
    long swf_str  = 0;
    long samp_str = 0;
    if (atmturb_ao_opt_wf_load(p.in_ampname, &pup_size, &nbframes, &amp_wf, &amp_str) != 0 ||
        atmturb_ao_opt_wf_load(p.in_swfname, &pup_size, &nbframes, &swf, &swf_str) != 0 ||
        atmturb_ao_opt_wf_load(p.in_sampname, &pup_size, &nbframes, &samp_wf, &samp_str) != 0)
    {
        atmturb_ao_wf_free(&samp_wf);
        atmturb_ao_wf_free(&swf);
        atmturb_ao_wf_free(&amp_wf);
        atmturb_ao_wf_free(&wf);
        return -1;
    }
    if (pup_size <= 0)
    {
        pup_size = wf.pup_size;
    }
    if (nbframes <= 0)
    {
        nbframes = wf.nbframes;
    }
    long ntot = pup_size * pup_size;

    long n_fft = 256;
    while (n_fft < 2 * pup_size)
    {
        n_fft *= 2;
    }
    long n_psf = (n_fft <= 512) ? n_fft : 512;

    if (p.pupil_scale_m <= 0.0 && CONF_PUPIL_SCALE > 0.0f)
    {
        p.pupil_scale_m = (double) CONF_PUPIL_SCALE;
    }

    float *pupil_mask = (float *) malloc(sizeof(float) * (size_t) ntot);
    float *fallback   = NULL;
    int    ret        = -1;
    if (pupil_mask == NULL)
    {
        goto cleanup;
    }
    atmturb_ao_init_pupil_mask(pup_size, p.tel_diam_m, p.pupil_scale_m, pupil_mask);

    {
        delete_image_ID(p.out_psfname);
        imageID id_psf   = create_2Dimage_ID(p.out_psfname, n_psf, n_psf);
        int     do_accum = (id_psf >= 0) ? 1 : 0;
        float  *psf_out  = (id_psf >= 0) ? dcimg[id_psf].array.F : NULL;
        if (psf_out == NULL)
        {
            fallback = (float *) calloc((size_t) (n_psf * n_psf), sizeof(float));
            psf_out  = fallback;
            if (psf_out == NULL)
            {
                goto cleanup;
            }
        }

        double                dl_peak        = 0.0;
        double                s_endpoints[2] = {0.0, 0.0};
        atmturb_cuda_ao_req_t creq           = {
            .params      = &p,
            .wf_data     = wf.data,
            .amp_data    = (samp_wf.data != NULL) ? samp_wf.data : amp_wf.data,
            .swf_data    = swf.data,
            .pupil_mask  = pupil_mask,
            .pup_size    = pup_size,
            .amp_stride  = (samp_wf.data != NULL) ? samp_str : amp_str,
            .swf_stride  = swf_str,
            .n_fft       = n_fft,
            .n_psf       = n_psf,
            .nbframes    = nbframes,
            .do_accum    = do_accum,
            .psf_out     = psf_out,
            .dl_peak_out = &dl_peak,
            .s_endpoints = s_endpoints
        };
        int ran_gpu = (atmturb_simd_is_gpu() && (nbframes >= 4 || atmturb_simd_is_gpu_forced()) &&
                       p.save_psfcube == 0 && atmturb_cuda_ao_sim_run(&creq) == 0);
        if (!ran_gpu && atmturb_ao_run_cpu(&creq) != 0)
        {
            goto cleanup;
        }

        float               max_cumul_f = 0.0f;
        const float *const  psf_r       = psf_out;
        #pragma omp simd reduction(max:max_cumul_f)
        for (long i = 0; i < n_psf * n_psf; i++)
        {
            if (psf_r[i] > max_cumul_f)
            {
                max_cumul_f = psf_r[i];
            }
        }

        double d_denom = do_accum ? ((double) nbframes * dl_peak) : dl_peak;
        double s_cumul = (d_denom > 0.0) ? ((double) max_cumul_f / d_denom) : 0.0;
        if (p.out_fitsname != NULL && p.out_fitsname[0] != '\0')
        {
            atmturb_ao_save_fits(p.out_fitsname, psf_out, n_psf, n_psf);
        }
        if (results != NULL)
        {
            results->strehl_cumul = s_cumul;
            results->strehl_first = s_endpoints[0];
            results->strehl_last  = s_endpoints[1];
            results->dl_peak      = dl_peak;
        }
    }
    ret = 0;

cleanup:
    free(fallback);
    free(pupil_mask);
    atmturb_ao_wf_free(&samp_wf);
    atmturb_ao_wf_free(&swf);
    atmturb_ao_wf_free(&amp_wf);
    atmturb_ao_wf_free(&wf);
    return ret;
}

/**
 * AtmosphericTurbulence_makePSF - Legacy wrapper for AO closed-loop PSF generation
 * @Kp: Proportional feedback gain.
 * @Ki: Integral feedback gain.
 * @Kd: Derivative feedback gain.
 * @Kdgain: Dynamic gain multiplier (unused).
 *
 * Return: Cumulative Strehl ratio on success, -1.0 on failure.
 */
double AtmosphericTurbulence_makePSF(
    double Kp,
    double Ki,
    double Kd,
    double Kdgain)
{
    (void) Kdgain;
    AtmosphericTurbulence_ReadConf();

    atmturb_ao_params_t p;
    atmturb_ao_init_params(&p);

    p.in_wfname   = "outarraypha";
    p.out_psfname = "PSFcumul";
    p.out_fitsname = "!PSFcumul.fits";
    p.lambda_ref_m = (double) CONF_LAMBDA;
    p.lambda_sci_m = 1.65e-6;
    p.pup_size     = CONF_WFsize;
    p.gain         = (Kp > 0.0) ? Kp : 0.5;
    p.loop_mode    = (Ki > 0.0 || Kd > 0.0) ? 2 : 1;
    p.Kp           = Kp;
    p.Ki           = Ki;
    p.Kd           = Kd;

    atmturb_ao_results_t res;
    if (atmturb_ao_sim_run(&p, &res) != 0)
    {
        return -1.0;
    }

    return res.strehl_cumul;
}
