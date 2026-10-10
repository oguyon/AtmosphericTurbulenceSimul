// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda.h
 * @brief   CUDA-accelerated multi-layer wavefront extrusion
 */

#ifndef ATMTURB_CUDA_H
#define ATMTURB_CUDA_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * struct atmturb_cuda_sim_params - Simulation configuration for CUDA renderer
 * @nblayers: Number of turbulent atmospheric layers.
 * @msize: Master screen linear dimension.
 * @pup_size: Linear dimension of the synthesized pupil.
 * @nbframes: Number of frames to synthesize.
 * @Scoeff: Chromatic dispersion scaling factor for secondary wavelength.
 * @h_masters: Array of host pointers to master screens.
 * @vxpix: Array of X velocities in pupil pixels/frame.
 * @vypix: Array of Y velocities in pupil pixels/frame.
 * @cn2: Array of layer Cn2 weights.
 * @x0: Optional array of initial X offsets in master pixels (NULL for vxpix).
 * @y0: Optional array of initial Y offsets in master pixels (NULL for vypix).
 */
typedef struct
{
    long                nblayers;
    long                msize;
    long                pup_size;
    long                nbframes;
    double              Scoeff;
    const float *const *h_masters;
    const double       *vxpix;
    const double       *vypix;
    const double       *cn2;
    const double       *x0;
    const double       *y0;
} atmturb_cuda_sim_params_t;

/**
 * struct atmturb_cuda_sim_outputs - Host output buffers for synthesized wavefronts
 * @pha: Primary phase 3D buffer.
 * @amp: Primary amplitude 3D buffer.
 * @spha: Secondary phase 3D buffer.
 * @samp: Secondary amplitude 3D buffer.
 */
typedef struct
{
    float *pha;
    float *amp;
    float *spha;
    float *samp;
} atmturb_cuda_sim_outputs_t;

#ifndef HAVE_CUDA

static inline int atmturb_cuda_device_available(void)
{
    return 0;
}

static inline int atmturb_wfs_render_frames_cuda(
    const atmturb_cuda_sim_params_t *params,
    atmturb_cuda_sim_outputs_t      *outputs)
{
    (void) params;
    (void) outputs;
    return -1;
}

static inline void atmturb_cuda_cleanup(void)
{
}

#else

/**
 * atmturb_cuda_device_available - Check if CUDA GPU is present and ready
 *
 * Return: 1 if CUDA device available, 0 otherwise.
 */
int atmturb_cuda_device_available(void);

/**
 * atmturb_wfs_render_frames_cuda - CUDA GPU accelerated wavefront time-series rendering
 * @params: Simulation parameters and input screen pointers.
 * @outputs: Destination host arrays for synthesized frames.
 *
 * Return: 0 on success, -1 on CUDA runtime error.
 */
int atmturb_wfs_render_frames_cuda(
    const atmturb_cuda_sim_params_t *params,
    atmturb_cuda_sim_outputs_t      *outputs);

/**
 * atmturb_cuda_sync_device_masters - Upload and cache master screens in GPU device memory
 * @msize: Master screen dimension in pixels.
 * @nblayers: Number of simulation layers.
 * @h_masters: Array of host screen pointers.
 *
 * Return: Pointer to device master screens array, or NULL on failure.
 */
float *atmturb_cuda_sync_device_masters(
    long                msize,
    long                nblayers,
    const float *const *h_masters);

/**
 * atmturb_cuda_cleanup - Release persistent GPU buffers and context
 */
void atmturb_cuda_cleanup(void);

#endif /* HAVE_CUDA */

#ifdef __cplusplus
}
#endif

#endif /* ATMTURB_CUDA_H */
