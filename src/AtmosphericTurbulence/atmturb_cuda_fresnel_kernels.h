// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_fresnel_kernels.h
 * @brief   GPU device kernel prototypes and launch wrappers for Fresnel propagation
 */

#ifndef ATMTURB_CUDA_FRESNEL_KERNELS_H
#define ATMTURB_CUDA_FRESNEL_KERNELS_H

#include <cuda_runtime.h>
#include <cufft.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * struct atmturb_cuda_fresnel_extract_args_t - Parameters for Fresnel pupil extraction kernel
 * @d_field: Input propagated complex electric field (`count * pad_size * pad_size`).
 * @d_pha: Input geometric pupil phase and output unwrapped phase (`count * pup_size * pup_size`).
 * @d_amp: Destination pupil amplitude array (or NULL if disabled).
 * @d_sum: Optional device scalar accumulator for piston sum reduction (or NULL).
 * @pad_size: Linear dimension of padded propagation grid in pixels.
 * @pup_size: Linear dimension of output pupil grid in pixels.
 * @guard: Guard band margin around pupil in pixels (`0` when unpadded).
 * @count: Number of batched frames (`1` for streaming).
 */
typedef struct
{
    const cufftComplex *d_field;
    float              *d_pha;
    float              *d_amp;
    float              *d_sum;
    int                 pad_size;
    int                 pup_size;
    int                 guard;
    int                 count;
} atmturb_cuda_fresnel_extract_args_t;

/**
 * atmturb_cuda_fresnel_init_modulate_launch - Initialize complex field from initial phase screen
 * @d_field: Destination complex electric field array (`count * pad_size * pad_size`).
 * @d_super_pha: Source phase screen on padded grid (`count * pad_size * pad_size`).
 * @d_sponge_1d: Optional 1D rim absorber window (`pad_size` floats, or NULL).
 * @pad_size: Linear dimension of padded propagation grid.
 * @count: Number of batched frames.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_fresnel_init_modulate_launch(
    cufftComplex *d_field,
    const float  *d_super_pha,
    const float  *d_sponge_1d,
    int           pad_size,
    int           count,
    cudaStream_t  stream);

/**
 * atmturb_cuda_fresnel_modulate_launch - Modulate propagating complex field by superlayer phase
 * @d_field: In-place complex electric field array (`count * pad_size * pad_size`).
 * @d_super_pha: Superlayer phase screen to apply (`count * pad_size * pad_size`).
 * @d_sponge_1d: Optional 1D rim absorber window (`pad_size` floats, or NULL).
 * @pad_size: Linear dimension of padded propagation grid.
 * @count: Number of batched frames.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_fresnel_modulate_launch(
    cufftComplex *d_field,
    const float  *d_super_pha,
    const float  *d_sponge_1d,
    int           pad_size,
    int           count,
    cudaStream_t  stream);

/**
 * atmturb_cuda_fresnel_tf1d_mul_launch - Multiply Fourier field by separable 1D Fresnel chirp
 * @d_field: Complex spatial frequency field array (`count * pad_size * pad_size`, in place).
 * @d_tf1d: 1D complex Fresnel chirp factor (`pad_size` elements).
 * @pad_size: Linear dimension of padded propagation grid.
 * @count: Number of batched frames.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_fresnel_tf1d_mul_launch(
    cufftComplex       *d_field,
    const cufftComplex *d_tf1d,
    int                 pad_size,
    int                 count,
    cudaStream_t        stream);

/**
 * atmturb_cuda_fresnel_extract_launch - Extract pupil amplitude, unwrapped phase, and piston sum
 * @args: Extraction descriptor with device pointers and grid dimensions.
 * @stream: CUDA execution stream.
 */
void atmturb_cuda_fresnel_extract_launch(
    const atmturb_cuda_fresnel_extract_args_t *args,
    cudaStream_t                               stream);

#ifdef __cplusplus
}
#endif

#endif /* ATMTURB_CUDA_FRESNEL_KERNELS_H */
