// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    wfprop_fresnel_cuda.h
 * @brief   CUDA-accelerated Fresnel diffractive wavefront propagation
 */

#ifndef WFPROP_FRESNEL_CUDA_H
#define WFPROP_FRESNEL_CUDA_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * wfprop_fresnel_device_available - Check if CUDA GPU is present and ready
 *
 * Return: 1 if CUDA device available, 0 otherwise.
 */
int wfprop_fresnel_device_available(void);

/**
 * wfprop_fresnel_propagate_cuda - CUDA GPU accelerated Fresnel propagation
 * @h_in: Pointer to input complex host array.
 * @h_out: Pointer to output complex host array.
 * @nx: Horizontal grid dimension.
 * @ny: Vertical grid dimension.
 * @pupil_scale: Physical pixel sampling scale [m/pixel].
 * @z: Propagation distance [m].
 * @lambda: Optical wavelength [m].
 * @is_double: 1 for double precision complex, 0 for single precision.
 *
 * Return: 0 on success, -1 on CUDA runtime error.
 */
int wfprop_fresnel_propagate_cuda(const void *h_in, void *h_out, long nx, long ny,
                                  double pupil_scale, double z, double lambda,
                                  int is_double);

#ifdef __cplusplus
}
#endif

#endif /* WFPROP_FRESNEL_CUDA_H */
