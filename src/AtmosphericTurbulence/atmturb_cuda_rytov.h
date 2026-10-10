// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_rytov.h
 * @brief   CUDA and cuFFT accelerated multi-layer Rytov propagation
 */

#ifndef ATMTURB_CUDA_RYTOV_H
#define ATMTURB_CUDA_RYTOV_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * struct atmturb_cuda_rytov_sublayer_t - Extrusion parameters for one sublayer
 * @k: Master screen layer index.
 * @w_pri: Primary wavelength weight factor.
 * @w_sec: Secondary wavelength weight factor.
 * @x: Primary X coordinate on master screen.
 * @y: Primary Y coordinate on master screen.
 * @xs: Secondary X coordinate on master screen.
 * @ys: Secondary Y coordinate on master screen.
 */
typedef struct
{
    int   k;
    float w_pri;
    float w_sec;
    float x;
    float y;
    float xs;
    float ys;
} atmturb_cuda_rytov_sublayer_t;

/**
 * struct atmturb_cuda_rytov_params_t - Parameters for CUDA Rytov execution
 * @msize: Master screen linear dimension in pixels.
 * @pup_size: Synthesized pupil linear dimension in pixels.
 * @pad_size: Padded compute grid linear dimension (pup_size + 2 * guard_pix).
 * @guard_pix: Guard band padding margin in pixels.
 * @nbframes: Number of simulation frames.
 * @nblayers: Total number of physical turbulence layers.
 * @nsuper: Number of superlayers.
 * @has_sec: 1 if secondary wavelength enabled, 0 otherwise.
 * @sec_shared: 1 if Option B shared forward FFT enabled, 0 for Option C.
 * @os: Master screen oversampling factor.
 * @interp: Sub-pixel interpolation mode (0: bilinear, 1: Keys bicubic).
 * @h_masters: Array of host master screen pointers.
 * @h_laplace_inv: Moisan Poisson inverse Laplace multipliers.
 * @h_exp_x: Moisan X boundary phase factor (complex float[2]).
 * @h_exp_y: Moisan Y boundary phase factor (complex float[2]).
 * @h_filt_a_pri: Primary phase filters per superlayer.
 * @h_filt_b_pri: Primary amplitude filters per superlayer.
 * @h_filt_a_sec: Secondary phase filters per superlayer.
 * @h_filt_b_sec: Secondary amplitude filters per superlayer.
 * @h_chrom_ramp: Option B chromatic phase ramp per superlayer (complex float[2]).
 * @frame_sublayers: Precomputed sublayer extrusion parameters [nbframes * nsuper * max_sublayers].
 * @super_nlayers: Number of physical layers in each superlayer [nsuper].
 * @super_dist_m: Line-of-sight distance for each superlayer in meters [nsuper].
 * @max_sublayers: Stride of sublayers array per (frame, superlayer).
 * @pha: Output primary phase 3D array (nbframes * pup_size * pup_size floats).
 * @amp: Output primary amplitude 3D array (nbframes * pup_size * pup_size floats).
 * @spha: Output secondary phase 3D array (nbframes * pup_size * pup_size floats).
 * @samp: Output secondary amplitude 3D array (nbframes * pup_size * pup_size floats).
 */
typedef struct
{
    long                                msize;
    long                                pup_size;
    long                                pad_size;
    long                                guard_pix;
    long                                nbframes;
    long                                nblayers;
    int                                 nsuper;
    int                                 has_sec;
    int                                 sec_shared;
    int                                 os;
    int                                 interp;
    int                                 use_moisan;

    const float *const                 *h_masters;

    const float                        *h_laplace_inv;
    const void                         *h_exp_x;
    const void                         *h_exp_y;
    const float *const                 *h_filt_a_pri;
    const float *const                 *h_filt_b_pri;
    const float *const                 *h_filt_a_sec;
    const float *const                 *h_filt_b_sec;
    const void *const                  *h_chrom_ramp;

    const atmturb_cuda_rytov_sublayer_t *frame_sublayers;
    const int                          *super_nlayers;
    const double                       *super_dist_m;
    int                                 max_sublayers;

    float                              *pha;
    float                              *amp;
    float                              *spha;
    float                              *samp;
} atmturb_cuda_rytov_params_t;

#ifndef HAVE_CUDA

static inline int atmturb_cuda_rytov_render(
    const atmturb_cuda_rytov_params_t *params)
{
    (void) params;
    return -1;
}

static inline void atmturb_cuda_rytov_cleanup(void)
{
}

#else

/**
 * atmturb_cuda_rytov_render - GPU-accelerated Rytov diffractive propagation
 * @params: Rytov simulation and filter parameters.
 *
 * Return: 0 on success, -1 on CUDA or cuFFT error.
 */
int atmturb_cuda_rytov_render(
    const atmturb_cuda_rytov_params_t *params);

/**
 * atmturb_cuda_rytov_cleanup - Free cached GPU buffers and cuFFT plans
 */
void atmturb_cuda_rytov_cleanup(void);

#endif /* HAVE_CUDA */

#ifdef __cplusplus
}
#endif

#endif /* ATMTURB_CUDA_RYTOV_H */
