// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_rytov_internal.h
 * @brief   Internal declarations for CUDA Rytov propagation compute kernels
 */

#ifndef ATMTURB_CUDA_RYTOV_INTERNAL_H
#define ATMTURB_CUDA_RYTOV_INTERNAL_H

#include <cuda_runtime.h>
#include <cufft.h>

#include "atmturb_cuda_piston.h"
#include "atmturb_cuda_rytov.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * atmturb_cuda_extrude_sl_launch - Launch superlayer phase screen extrusion and subharmonics
 * @interp_mode: Interpolation scheme (0=bilinear, 1=Keys bicubic).
 * @d_masters: Array of device pointers to master phase screens.
 * @msize: Master phase screen linear dimension in pixels.
 * @pad_size: Linear dimension of padded simulation grid in pixels.
 * @pup_size: Linear dimension of synthesized pupil in pixels.
 * @guard: Guard band margin in pixels surrounding pupil.
 * @n_sub: Number of active sublayers in current superlayer.
 * @has_sec: 1 if secondary wavelength is enabled, 0 otherwise.
 * @os: Master screen oversampling factor (1 or 2).
 * @d_sublayers: Array of device sublayer coordinates and interpolation weights.
 * @n_lf_modes: Number of active low-frequency subharmonic modes.
 * @d_lf_modes: Array of device subharmonic mode parameters.
 * @d_super_pha: Destination superlayer primary phase array.
 * @d_super_spha: Destination superlayer secondary phase array.
 * @d_frame_pha: Cumulative pupil primary phase array.
 * @d_frame_spha: Cumulative pupil secondary phase array.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_extrude_sl_launch(
    int                                  interp_mode,
    const float                         *d_masters,
    int                                  msize,
    int                                  pad_size,
    int                                  pup_size,
    int                                  guard,
    int                                  n_sub,
    int                                  has_sec,
    int                                  os,
    const atmturb_cuda_rytov_sublayer_t *d_sublayers,
    int                                  n_lf_modes,
    const atmturb_cuda_lowfreq_mode_t   *d_lf_modes,
    float                               *d_super_pha,
    float                               *d_super_spha,
    float                               *d_frame_pha,
    float                               *d_frame_spha,
    cudaStream_t                         stream);

/**
 * atmturb_cuda_lowfreq_accumulate_launch - Launch analytic low-frequency mode accumulation
 * @pad_size: Linear dimension of padded simulation grid in pixels.
 * @pup_size: Linear dimension of synthesized pupil in pixels.
 * @guard: Guard margin in pixels.
 * @n_modes: Number of active low-frequency modes to evaluate.
 * @has_sec: 1 if secondary wavelength enabled, 0 otherwise.
 * @d_modes: Array of device low-frequency mode parameters.
 * @d_super_pha: Destination superlayer primary phase array.
 * @d_super_spha: Destination superlayer secondary phase array.
 * @d_frame_pha: Cumulative pupil primary phase array.
 * @d_frame_spha: Cumulative pupil secondary phase array.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_lowfreq_accumulate_launch(
    int                                pad_size,
    int                                pup_size,
    int                                guard,
    int                                n_modes,
    int                                has_sec,
    const atmturb_cuda_lowfreq_mode_t *d_modes,
    float                             *d_super_pha,
    float                             *d_super_spha,
    float                             *d_frame_pha,
    float                             *d_frame_spha,
    cudaStream_t                       stream);

/**
 * atmturb_cuda_extract_bounds_launch - Extract boundary pixel gradients for Moisan decomposition
 * @d_super_pha: Input superlayer phase screen on padded grid.
 * @d_bounds: Output boundary gradient array.
 * @pad_size: Linear dimension of padded grid in pixels.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_extract_bounds_launch(
    const float *d_super_pha,
    float       *d_bounds,
    int          pad_size,
    cudaStream_t stream);

/**
 * atmturb_cuda_moisan_adjust_launch - Solve Poisson boundary correction in frequency space
 * @d_spec: In-place 2D Fourier spectrum of phase screen.
 * @d_hat_bounds: 1D Fourier transforms of boundary row/col differences.
 * @d_exp_x: Precomputed complex exponential roots along X.
 * @d_exp_y: Precomputed complex exponential roots along Y.
 * @d_laplace_inv: Inverse discrete Laplacian grid values.
 * @pad_size: Linear dimension of padded grid.
 * @n_half: Half-complex dimension along X (pad_size/2 + 1).
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_moisan_adjust_launch(
    cufftComplex       *d_spec,
    const cufftComplex *d_hat_bounds,
    const cufftComplex *d_exp_x,
    const cufftComplex *d_exp_y,
    const float        *d_laplace_inv,
    int                 pad_size,
    int                 n_half,
    cudaStream_t        stream);

/**
 * struct atmturb_cuda_moisan_mac - Parameter container for fused Moisan adjustment and filter MAC
 * @d_acc_dphi: Cumulative phase perturbation spectrum buffer.
 * @d_acc_chi: Cumulative log-amplitude perturbation spectrum buffer.
 * @d_spec: Input uncorrected 2D DFT spectrum of superlayer phase screen.
 * @d_hat_bounds: 1D Fourier transforms of boundary row/column differences.
 * @d_exp_x: Precomputed complex roots along X.
 * @d_exp_y: Precomputed complex roots along Y.
 * @d_laplace_inv: Precomputed inverse discrete Laplacian values.
 * @d_filt_a: Rytov phase filter transfer function kernel.
 * @d_filt_b: Rytov log-amplitude filter transfer function kernel.
 * @pad_size: Padded grid linear dimension.
 * @n_half: Half-complex linear dimension along X.
 */
typedef struct
{
    cufftComplex       *d_acc_dphi;
    cufftComplex       *d_acc_chi;
    const cufftComplex *d_spec;
    const cufftComplex *d_hat_bounds;
    const cufftComplex *d_exp_x;
    const cufftComplex *d_exp_y;
    const float        *d_laplace_inv;
    const float        *d_filt_a;
    const float        *d_filt_b;
    int                 pad_size;
    int                 n_half;
} atmturb_cuda_moisan_mac_t;

/**
 * struct atmturb_cuda_moisan_dual_mac - Container for dual-channel fused Moisan adjustment and MAC
 * @d_acc_dphi_pri: Cumulative primary phase perturbation spectrum.
 * @d_acc_chi_pri: Cumulative primary log-amplitude perturbation spectrum.
 * @d_acc_dphi_sec: Cumulative secondary phase perturbation spectrum.
 * @d_acc_chi_sec: Cumulative secondary log-amplitude perturbation spectrum.
 * @d_spec: Input uncorrected 2D DFT spectrum.
 * @d_ramp: Chromatic differential spatial shift phase ramp.
 * @d_hat_bounds: 1D Fourier transforms of boundary row/column differences.
 * @d_exp_x: Precomputed complex roots along X.
 * @d_exp_y: Precomputed complex roots along Y.
 * @d_laplace_inv: Precomputed inverse discrete Laplacian values.
 * @d_filt_a_pri: Primary wavelength Rytov phase filter transfer function.
 * @d_filt_b_pri: Primary wavelength Rytov log-amplitude filter transfer function.
 * @d_filt_a_sec: Secondary wavelength Rytov phase filter transfer function.
 * @d_filt_b_sec: Secondary wavelength Rytov log-amplitude filter transfer function.
 * @pad_size: Padded grid linear dimension.
 * @n_half: Half-complex linear dimension along X.
 */
typedef struct
{
    cufftComplex       *d_acc_dphi_pri;
    cufftComplex       *d_acc_chi_pri;
    cufftComplex       *d_acc_dphi_sec;
    cufftComplex       *d_acc_chi_sec;
    const cufftComplex *d_spec;
    const cufftComplex *d_ramp;
    const cufftComplex *d_hat_bounds;
    const cufftComplex *d_exp_x;
    const cufftComplex *d_exp_y;
    const float        *d_laplace_inv;
    const float        *d_filt_a_pri;
    const float        *d_filt_b_pri;
    const float        *d_filt_a_sec;
    const float        *d_filt_b_sec;
    int                 pad_size;
    int                 n_half;
} atmturb_cuda_moisan_dual_mac_t;

/**
 * atmturb_cuda_moisan_filter_mac_launch - Launch single-channel fused Moisan adjustment and MAC
 * @args: Parameter container with device buffer pointers and grid dimensions.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_moisan_filter_mac_launch(
    const atmturb_cuda_moisan_mac_t *args,
    cudaStream_t                     stream);

/**
 * atmturb_cuda_moisan_dual_mac_launch - Launch dual-channel fused Moisan adjustment and MAC
 * @args: Parameter container with device buffer pointers and grid dimensions.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_moisan_dual_mac_launch(
    const atmturb_cuda_moisan_dual_mac_t *args,
    cudaStream_t                          stream);

/**
 * atmturb_cuda_filter_mac_launch - Launch Rytov spectral multiply-accumulate pass
 * @d_acc_dphi: Cumulative phase perturbation spectrum buffer.
 * @d_acc_chi: Cumulative log-amplitude perturbation spectrum buffer.
 * @d_spec: Corrected periodic 2D DFT spectrum of superlayer.
 * @d_filt_a: Rytov phase filter kernel array.
 * @d_filt_b: Rytov log-amplitude filter kernel array.
 * @ntot: Total half-complex spectrum elements (pad_size * n_half).
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_filter_mac_launch(
    cufftComplex       *d_acc_dphi,
    cufftComplex       *d_acc_chi,
    const cufftComplex *d_spec,
    const float        *d_filt_a,
    const float        *d_filt_b,
    int                 ntot,
    cudaStream_t        stream);

/**
 * atmturb_cuda_filter_mac_dual_launch - Launch dual-wavelength Rytov spectral MAC pass
 * @d_acc_dphi_pri: Cumulative primary phase perturbation spectrum.
 * @d_acc_chi_pri: Cumulative primary log-amplitude perturbation spectrum.
 * @d_acc_dphi_sec: Cumulative secondary phase perturbation spectrum.
 * @d_acc_chi_sec: Cumulative secondary log-amplitude perturbation spectrum.
 * @d_spec: Corrected periodic 2D DFT spectrum of superlayer.
 * @d_ramp: Chromatic phase shift ramp.
 * @d_filt_a_pri: Primary Rytov phase filter array.
 * @d_filt_b_pri: Primary Rytov log-amplitude filter array.
 * @d_filt_a_sec: Secondary Rytov phase filter array.
 * @d_filt_b_sec: Secondary Rytov log-amplitude filter array.
 * @ntot: Total half-complex spectrum elements.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_filter_mac_dual_launch(
    cufftComplex       *d_acc_dphi_pri,
    cufftComplex       *d_acc_chi_pri,
    cufftComplex       *d_acc_dphi_sec,
    cufftComplex       *d_acc_chi_sec,
    const cufftComplex *d_spec,
    const cufftComplex *d_ramp,
    const float        *d_filt_a_pri,
    const float        *d_filt_b_pri,
    const float        *d_filt_a_sec,
    const float        *d_filt_b_sec,
    int                 ntot,
    cudaStream_t        stream);

/**
 * atmturb_cuda_assemble_launch - Assemble pupil phase and log-amplitude into complex wavefront
 * @pup_size: Linear dimension of synthesized pupil in pixels.
 * @guard: Guard margin in pixels.
 * @pad_size: Linear dimension of padded simulation grid.
 * @d_frame_pha: Cumulative geometric pupil phase (in-place updated).
 * @d_frame_amp: Output pupil amplitude array (exp(chi)).
 * @d_dphi_out: Synthesized phase perturbation from inverse 2D DFT.
 * @d_chi_out: Synthesized log-amplitude perturbation from inverse 2D DFT.
 * @d_sum_I: Accumulator buffer for mean intensity over pupil aperture.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_assemble_launch(
    int          pup_size,
    int          guard,
    int          pad_size,
    float       *d_frame_pha,
    float       *d_frame_amp,
    const float *d_dphi_out,
    const float *d_chi_out,
    float       *d_sum_I,
    cudaStream_t stream);

/**
 * atmturb_cuda_normalize_launch - Normalize pupil amplitude array to enforce energy conservation
 * @d_frame_amp: In-place pupil amplitude array.
 * @npix: Total pupil aperture pixels (pup_size * pup_size).
 * @d_sum_I: Pointer to total intensity sum across pupil.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_normalize_launch(
    float       *d_frame_amp,
    int          npix,
    const float *d_sum_I,
    cudaStream_t stream);

#ifdef __cplusplus
}
#endif

#endif // ATMTURB_CUDA_RYTOV_INTERNAL_H
