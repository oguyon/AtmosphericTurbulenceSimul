// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_psf.h
 * @brief   CUDA/cuFFT open-loop and closed-loop AO science PSF simulation engine
 */

#ifndef ATMTURB_CUDA_PSF_H
#define ATMTURB_CUDA_PSF_H

#include "AtmosphericTurbulence.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * struct atmturb_cuda_ao_req - GPU AO closed-loop / open-loop PSF simulation request
 * @params: Active AO controller and optical wavelength parameters.
 * @wf_data: Host 3D input phase cube (`nbframes * pup_size * pup_size` floats, or NULL for flat).
 * @pupil_mask: Host 2D circular pupil mask (`pup_size * pup_size` floats).
 * @pup_size: Pupil linear dimension in pixels.
 * @n_fft: Zero-padded 2D FFT linear dimension (power of 2, >= 2 * pup_size).
 * @n_psf: Cropped output focal PSF dimension (`min(n_fft, 512)`).
 * @nbframes: Number of frames to simulate.
 * @do_accum: 1 to accumulate frames into `psf_out`, 0 to overwrite each frame.
 * @psf_out: Destination host 2D PSF array (`n_psf * n_psf` floats).
 * @dl_peak_out: Output diffraction-limited flat-wavefront peak intensity.
 * @s_endpoints: Output 2-element array `{strehl_first, strehl_last}`.
 */
typedef struct
{
    const atmturb_ao_params_t *params;
    const float               *wf_data;
    const float               *pupil_mask;
    long                       pup_size;
    long                       n_fft;
    long                       n_psf;
    long                       nbframes;
    int                        do_accum;
    float                     *psf_out;
    double                    *dl_peak_out;
    double                    *s_endpoints;
} atmturb_cuda_ao_req_t;

#ifdef HAVE_CUDA

/**
 * atmturb_cuda_ao_sim_run - Execute open-loop or closed-loop AO PSF synthesis on CUDA GPU
 * @req: Simulation request descriptor.
 *
 * Return: 0 on success, -1 on unsupported configuration or CUDA/cuFFT failure.
 */
int atmturb_cuda_ao_sim_run(
    const atmturb_cuda_ao_req_t *req);

/**
 * atmturb_cuda_psf_cleanup - Release cached cuFFT plans and stream for GPU AO PSF engine
 */
void atmturb_cuda_psf_cleanup(void);

#else /* !HAVE_CUDA */

static inline int atmturb_cuda_ao_sim_run(
    const atmturb_cuda_ao_req_t *req)
{
    (void) req;
    return -1;
}

static inline void atmturb_cuda_psf_cleanup(void)
{
}

#endif /* HAVE_CUDA */

#ifdef __cplusplus
}
#endif

#endif /* ATMTURB_CUDA_PSF_H */
