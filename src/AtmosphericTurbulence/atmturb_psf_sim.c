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
#include <fftw3.h>

#include "AtmosphericTurbulence.h"
#include "atmturb_types.h"

/**
 * struct atmturb_ao_psf_engine_t - Fourier optics engine for Nyquist PSF formation
 * @n_fft: Zero-padded 2D grid dimension (power of 2, >= 2 * pup_size).
 * @n_psf: Output cropped PSF dimension (min(n_fft, 512)).
 * @off_x: Pupil embedding X offset in zero-padded grid.
 * @off_y: Pupil embedding Y offset in zero-padded grid.
 * @norm_factor: Parseval intensity normalization factor (1 / n_fft^2).
 * @in_grid: Pre-allocated complex pupil input array.
 * @out_grid: Pre-allocated complex focal output array.
 * @plan: Cached FFTW 2D DFT forward execution plan.
 */
typedef struct
{
    long           n_fft;
    long           n_psf;
    long           off_x;
    long           off_y;
    double         norm_factor;
    fftwf_complex *in_grid;
    fftwf_complex *out_grid;
    fftwf_plan     plan;
} atmturb_ao_psf_engine_t;

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
    params->lambda_ref_m  = 0.55e-6;
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
    double r_pix = (pupil_scale_m > 0.0)
                       ? (0.5 * tel_diam_m / pupil_scale_m)
                       : (0.5 * (double) pup_size);
    double r2 = r_pix * r_pix;
    double cx = 0.5 * ((double) pup_size - 1.0);
    double cy = 0.5 * ((double) pup_size - 1.0);
    double npup = 0.0;

    for (long j = 0; j < pup_size; j++)
    {
        double dy = (double) j - cy;
        for (long i = 0; i < pup_size; i++)
        {
            double dx = (double) i - cx;
            long idx = j * pup_size + i;
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
 * atmturb_ao_psf_engine_init - Allocate zero-padded grid and create cached FFTW plan
 * @eng: PSF engine structure to initialize.
 * @pup_size: Pupil dimension in pixels.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
static int atmturb_ao_psf_engine_init(
    atmturb_ao_psf_engine_t *eng,
    long                     pup_size)
{
    long n_fft = 256;
    while (n_fft < 2 * pup_size)
    {
        n_fft *= 2;
    }

    eng->n_fft = n_fft;
    eng->n_psf = (n_fft <= 512) ? n_fft : 512;
    eng->off_x = (n_fft - pup_size) / 2;
    eng->off_y = (n_fft - pup_size) / 2;
    eng->norm_factor = 1.0 / ((double) n_fft * (double) n_fft);

    eng->in_grid = (fftwf_complex *) fftwf_alloc_complex((size_t) (n_fft * n_fft));
    eng->out_grid = (fftwf_complex *) fftwf_alloc_complex((size_t) (n_fft * n_fft));
    if (eng->in_grid == NULL || eng->out_grid == NULL)
    {
        fftwf_free(eng->in_grid);
        fftwf_free(eng->out_grid);
        return -1;
    }

    eng->plan = fftwf_plan_dft_2d((int) n_fft, (int) n_fft,
                                  eng->in_grid, eng->out_grid,
                                  FFTW_FORWARD, FFTW_ESTIMATE);
    if (!eng->plan)
    {
        fftwf_free(eng->in_grid);
        fftwf_free(eng->out_grid);
        return -1;
    }

    return 0;
}

/**
 * atmturb_ao_psf_engine_free - Release FFTW plan and complex arrays
 * @eng: Initialized PSF engine.
 */
static void atmturb_ao_psf_engine_free(
    atmturb_ao_psf_engine_t *eng)
{
    if (eng->plan)
    {
        fftwf_destroy_plan(eng->plan);
        eng->plan = NULL;
    }
    fftwf_free(eng->in_grid);
    fftwf_free(eng->out_grid);
    eng->in_grid = NULL;
    eng->out_grid = NULL;
}

/**
 * atmturb_ao_render_psf - Form instantaneous zero-padded centered PSF
 * @eng: Initialized PSF engine.
 * @opd_res: Residual optical path difference in meters [pup_size * pup_size].
 * @pupil_mask: Circular pupil mask [pup_size * pup_size].
 * @pup_size: Pupil grid linear dimension.
 * @lambda_sci: Science observing wavelength in meters.
 * @psf_out: Output 2D PSF buffer [n_psf * n_psf].
 *
 * Return: Peak pixel intensity of the rendered PSF slice.
 */
static double atmturb_ao_render_psf(
    const atmturb_ao_psf_engine_t *eng,
    const float                   *opd_res,
    const float                   *pupil_mask,
    long                           pup_size,
    double                         lambda_sci,
    float                         *psf_out)
{
    long n_fft = eng->n_fft;
    long n_psf = eng->n_psf;
    memset(eng->in_grid, 0, sizeof(fftwf_complex) * (size_t) (n_fft * n_fft));

    float k_sci = (float) (2.0 * M_PI / lambda_sci);
    for (long j = 0; j < pup_size; j++)
    {
        long in_row = (eng->off_y + j) * n_fft + eng->off_x;
        long pup_row = j * pup_size;
        for (long i = 0; i < pup_size; i++)
        {
            float m = pupil_mask[pup_row + i];
            if (m > 0.0f)
            {
                float p = (opd_res != NULL) ? (opd_res[pup_row + i] * k_sci) : 0.0f;
                eng->in_grid[in_row + i][0] = m * cosf(p);
                eng->in_grid[in_row + i][1] = m * sinf(p);
            }
        }
    }

    fftwf_execute(eng->plan);

    long u_start = (n_fft - n_psf) / 2;
    long v_start = (n_fft - n_psf) / 2;
    double peak = 0.0;
    float norm = (float) eng->norm_factor;

    for (long vj = 0; vj < n_psf; vj++)
    {
        long v = v_start + vj;
        long fy = (v < n_fft / 2) ? (v + n_fft / 2) : (v - n_fft / 2);
        long out_row = fy * n_fft;
        long psf_row = vj * n_psf;

        for (long ui = 0; ui < n_psf; ui++)
        {
            long u = u_start + ui;
            long fx = (u < n_fft / 2) ? (u + n_fft / 2) : (u - n_fft / 2);
            long idx = out_row + fx;
            float re = eng->out_grid[idx][0];
            float im = eng->out_grid[idx][1];
            float val = (re * re + im * im) * norm;
            psf_out[psf_row + ui] = val;
            if ((double) val > peak)
            {
                peak = (double) val;
            }
        }
    }

    return peak;
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
    ctrl->head = 0;
    ctrl->ntot = ntot;

    ctrl->dm_opd = (float *) calloc((size_t) ntot, sizeof(float));
    ctrl->pid_int = (float *) calloc((size_t) ntot, sizeof(float));
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
 * atmturb_ao_controller_step - Update DM and compute residual OPD for current frame
 * @ctrl: Active controller state.
 * @params: Simulation parameters.
 * @opd_atm: Atmospheric optical path difference in meters [ntot].
 * @pupil_mask: Circular pupil mask [ntot].
 * @opd_res_out: Output residual optical path difference in meters [ntot].
 */
static void atmturb_ao_controller_step(
    atmturb_ao_controller_t   *ctrl,
    const atmturb_ao_params_t *params,
    const float               *opd_atm,
    const float               *pupil_mask,
    float                     *opd_res_out)
{
    long ntot = ctrl->ntot;
    float *slot = ctrl->delay_ring[ctrl->head];

    for (long i = 0; i < ntot; i++)
    {
        float res = (opd_atm != NULL) ? (opd_atm[i] - ctrl->dm_opd[i]) : (-ctrl->dm_opd[i]);
        slot[i] = res * pupil_mask[i];
        opd_res_out[i] = res;
    }

    int read_idx = (ctrl->head - ctrl->delay + ctrl->depth) % ctrl->depth;
    const float *delayed_err = ctrl->delay_ring[read_idx];

    if (params->loop_mode == 1)
    {
        float leak = (float) params->leak;
        float gain = (float) params->gain;
        for (long i = 0; i < ntot; i++)
        {
            if (pupil_mask[i] > 0.0f)
            {
                ctrl->dm_opd[i] = (1.0f - leak) * ctrl->dm_opd[i] + gain * delayed_err[i];
            }
            else
            {
                ctrl->dm_opd[i] = (1.0f - leak) * ctrl->dm_opd[i];
            }
        }
    }
    else if (params->loop_mode == 2)
    {
        float kp = (float) params->Kp;
        float ki = (float) params->Ki;
        float kd = (float) params->Kd;
        float dt = 0.001f;
        for (long i = 0; i < ntot; i++)
        {
            if (pupil_mask[i] > 0.0f)
            {
                float e = delayed_err[i];
                float der = (e - ctrl->prev_err[i]) / dt;
                ctrl->pid_int[i] += e * dt;
                float delta = kp * e + ki * ctrl->pid_int[i] + kd * der;
                ctrl->dm_opd[i] += delta;
                ctrl->prev_err[i] = e;
            }
        }
    }

    ctrl->head = (ctrl->head + 1) % ctrl->depth;
}

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
        wf->data = dcimg[id].array.F;
        wf->pup_size = dcimg[id].md[0].size[0];
        wf->nbframes = (dcimg[id].md[0].naxis >= 3) ? dcimg[id].md[0].size[2] : 1;
        return 0;
    }

    if (ImageStreamIO_openIm(&wf->shm_img, wfname) == IMAGESTREAMIO_SUCCESS)
    {
        wf->data = wf->shm_img.array.F;
        wf->pup_size = wf->shm_img.md->size[0];
        wf->nbframes = (wf->shm_img.md->naxis >= 3) ? wf->shm_img.md->size[2] : 1;
        wf->is_shm = 1;
        return 0;
    }

    char tmpname[512];
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

    if (target != NULL)
    {
        fitsfile *fptr = NULL;
        int status = 0;
        fits_open_file(&fptr, target, READONLY, &status);
        if (status == 0)
        {
            int naxis = 0;
            long naxes[3] = {1, 1, 1};
            fits_get_img_dim(fptr, &naxis, &status);
            fits_get_img_size(fptr, (naxis > 3) ? 3 : naxis, naxes, &status);
            long pup_size = naxes[0];
            long nz = (naxis >= 3) ? naxes[2] : 1;
            long ntot = pup_size * pup_size * nz;
            float *arr = (float *) malloc(sizeof(float) * (size_t) ntot);
            if (arr != NULL)
            {
                long fpixel[3] = {1, 1, 1};
                int anynul = 0;
                fits_read_pix(fptr, TFLOAT, fpixel, ntot, NULL, arr, &anynul, &status);
                fits_close_file(fptr, &status);
                if (status == 0)
                {
                    wf->data = arr;
                    wf->pup_size = pup_size;
                    wf->nbframes = nz;
                    wf->is_allocated = 1;
                    return 0;
                }
                free(arr);
            }
            else
            {
                fits_close_file(fptr, &status);
            }
        }
    }

    printf("ERROR: cannot find or load wavefront stream or file '%s'\n", wfname);
    return -1;
}

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

static int atmturb_ao_save_fits(
    const char  *filename,
    const float *data,
    long         nx,
    long         ny)
{
    fitsfile *fptr = NULL;
    int status = 0;
    remove(filename);
    fits_create_file(&fptr, filename, &status);
    if (status != 0)
    {
        return -1;
    }
    long naxes[2] = {nx, ny};
    fits_create_img(fptr, FLOAT_IMG, 2, naxes, &status);
    long fpixel[2] = {1, 1};
    fits_write_pix(fptr, TFLOAT, fpixel, nx * ny, (void *) data, &status);
    fits_close_file(fptr, &status);
    return (status == 0) ? 0 : -1;
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

    atmturb_ao_wf_t wf;
    if (atmturb_ao_wf_load(p.in_wfname, p.pup_size, &wf) != 0)
    {
        return -1;
    }

    long pup_size = wf.pup_size;
    long nz = wf.nbframes;
    long nbframes = (p.nbframes > 0 && p.nbframes < nz) ? p.nbframes : nz;
    long ntot = pup_size * pup_size;

    float *pupil_mask = (float *) malloc(sizeof(float) * (size_t) ntot);
    if (pupil_mask == NULL)
    {
        atmturb_ao_wf_free(&wf);
        return -1;
    }
    atmturb_ao_init_pupil_mask(pup_size, p.tel_diam_m, p.pupil_scale_m, pupil_mask);

    atmturb_ao_psf_engine_t eng;
    if (atmturb_ao_psf_engine_init(&eng, pup_size) != 0)
    {
        free(pupil_mask);
        atmturb_ao_wf_free(&wf);
        return -1;
    }

    long n_psf = eng.n_psf;
    float *tmp_psf = (float *) malloc(sizeof(float) * (size_t) (n_psf * n_psf));
    float *opd_frame = (float *) malloc(sizeof(float) * (size_t) ntot);
    float *opd_res = (float *) malloc(sizeof(float) * (size_t) ntot);
    atmturb_ao_controller_t ctrl;
    if (tmp_psf == NULL || opd_frame == NULL || opd_res == NULL ||
        atmturb_ao_controller_init(&ctrl, ntot, p.loop_delay) != 0)
    {
        free(pupil_mask);
        free(tmp_psf);
        free(opd_frame);
        free(opd_res);
        atmturb_ao_psf_engine_free(&eng);
        atmturb_ao_wf_free(&wf);
        return -1;
    }

    double dl_peak = atmturb_ao_render_psf(&eng, NULL, pupil_mask, pup_size,
                                          p.lambda_sci_m, tmp_psf);

    delete_image_ID(p.out_psfname);
    imageID id_psf = create_2Dimage_ID(p.out_psfname, n_psf, n_psf);

    float opd_scale = (float) (p.lambda_ref_m / (2.0 * M_PI));
    double s_first = 0.0, s_last = 0.0;

    for (long t = 0; t < nbframes; t++)
    {
        if (wf.data != NULL)
        {
            const float *slice = &wf.data[t * ntot];
            for (long i = 0; i < ntot; i++)
            {
                opd_frame[i] = slice[i] * opd_scale;
            }
        }
        else
        {
            memset(opd_frame, 0, sizeof(float) * (size_t) ntot);
        }

        atmturb_ao_controller_step(&ctrl, &p, opd_frame, pupil_mask, opd_res);
        double slice_peak = atmturb_ao_render_psf(&eng, opd_res, pupil_mask, pup_size,
                                                 p.lambda_sci_m, tmp_psf);
        if (id_psf >= 0)
        {
            for (long i = 0; i < n_psf * n_psf; i++)
            {
                dcimg[id_psf].array.F[i] += tmp_psf[i];
            }
        }

        double s_slice = (dl_peak > 0.0) ? (slice_peak / dl_peak) : 0.0;
        if (t == 0) s_first = s_slice;
        if (t == nbframes - 1) s_last = s_slice;
    }

    double max_cumul = 0.0;
    if (id_psf >= 0)
    {
        for (long i = 0; i < n_psf * n_psf; i++)
        {
            if ((double) dcimg[id_psf].array.F[i] > max_cumul)
            {
                max_cumul = (double) dcimg[id_psf].array.F[i];
            }
        }
    }

    double s_cumul = (dl_peak > 0.0 && nbframes > 0)
                         ? (max_cumul / ((double) nbframes * dl_peak))
                         : 0.0;

    if (p.out_fitsname != NULL && p.out_fitsname[0] != '\0')
    {
        const float *psf_buf = (id_psf >= 0) ? dcimg[id_psf].array.F : tmp_psf;
        atmturb_ao_save_fits(p.out_fitsname, psf_buf, n_psf, n_psf);
    }

    if (results != NULL)
    {
        results->strehl_cumul = s_cumul;
        results->strehl_first = s_first;
        results->strehl_last  = s_last;
        results->dl_peak      = dl_peak;
    }

    free(pupil_mask);
    free(tmp_psf);
    free(opd_frame);
    free(opd_res);
    atmturb_ao_controller_free(&ctrl);
    atmturb_ao_psf_engine_free(&eng);
    atmturb_ao_wf_free(&wf);

    return 0;
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

    p.in_wfname = "outarraypha";
    p.out_psfname = "PSFcumul";
    p.out_fitsname = "!PSFcumul.fits";
    p.lambda_ref_m = (double) CONF_LAMBDA;
    p.lambda_sci_m = 1.65e-6;
    p.pup_size = CONF_WFsize;
    p.gain = (Kp > 0.0) ? Kp : 0.5;
    p.loop_mode = (Ki > 0.0 || Kd > 0.0) ? 2 : 1;
    p.Kp = Kp;
    p.Ki = Ki;
    p.Kd = Kd;

    atmturb_ao_results_t res;
    if (atmturb_ao_sim_run(&p, &res) != 0)
    {
        return -1.0;
    }

    return res.strehl_cumul;
}
