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

#include "atmturb_cuda_rytov.h"

#ifdef __cplusplus
extern "C" {
#endif

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
    float                               *d_super_pha,
    float                               *d_super_spha,
    float                               *d_frame_pha,
    float                               *d_frame_spha,
    cudaStream_t                         stream);

void atmturb_cuda_extract_bounds_launch(
    const float *d_super_pha,
    float       *d_bounds,
    int          pad_size,
    cudaStream_t stream);

void atmturb_cuda_moisan_adjust_launch(
    cufftComplex       *d_spec,
    const cufftComplex *d_hat_bounds,
    const cufftComplex *d_exp_x,
    const cufftComplex *d_exp_y,
    const float        *d_laplace_inv,
    int                 pad_size,
    int                 n_half,
    cudaStream_t        stream);

void atmturb_cuda_filter_mac_launch(
    cufftComplex       *d_acc_dphi,
    cufftComplex       *d_acc_chi,
    const cufftComplex *d_spec,
    const float        *d_filt_a,
    const float        *d_filt_b,
    int                 ntot,
    cudaStream_t        stream);

void atmturb_cuda_filter_mac_rotated_launch(
    cufftComplex       *d_acc_dphi_sec,
    cufftComplex       *d_acc_chi_sec,
    const cufftComplex *d_spec,
    const cufftComplex *d_ramp,
    const float        *d_filt_a_sec,
    const float        *d_filt_b_sec,
    int                 ntot,
    cudaStream_t        stream);

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

void atmturb_cuda_normalize_launch(
    float       *d_frame_amp,
    int          npix,
    const float *d_sum_I,
    cudaStream_t stream);

#ifdef __cplusplus
}
#endif

#endif // ATMTURB_CUDA_RYTOV_INTERNAL_H
