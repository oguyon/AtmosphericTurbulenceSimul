// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_fresnel_simd_scalar.c
 * @brief   Scalar fallback implementations for Fresnel propagation kernels
 */

#define _GNU_SOURCE
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "atmturb_fresnel_simd.h"

/**
 * atmturb_fresnel_modulate_field_scalar - Scalar fallback for complex field modulation
 * @field: In-out complex field array.
 * @phase: Real phase screen array in radians.
 * @npix: Total number of pixels.
 */
void atmturb_fresnel_modulate_field_scalar(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    npix)
{
#ifdef _OPENMP
    #pragma omp parallel for if(!omp_in_parallel() && npix >= 65536) schedule(static)
#endif
    for (long i = 0; i < npix; i++)
    {
        float p = phase[i];
        float s, c;
        sincosf(p, &s, &c);
        float re    = field[i][0];
        float im    = field[i][1];
        field[i][0] = re * c - im * s;
        field[i][1] = re * s + im * c;
    }
}

/**
 * atmturb_fresnel_init_modulate_field_scalar - Scalar fallback for unit plane wave init
 * @field: Output complex field array.
 * @phase: Real phase screen array in radians.
 * @npix: Total number of pixels.
 */
void atmturb_fresnel_init_modulate_field_scalar(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    npix)
{
#ifdef _OPENMP
    #pragma omp parallel for if(!omp_in_parallel() && npix >= 65536) schedule(static)
#endif
    for (long i = 0; i < npix; i++)
    {
        float p = phase[i];
        float s, c;
        sincosf(p, &s, &c);
        field[i][0] = c;
        field[i][1] = s;
    }
}

/**
 * atmturb_fresnel_extract_unwrapped_scalar - Scalar fallback for unwrapped phase extraction
 * @field: Complex optical field at pupil.
 * @pha: Phase array holding accumulated geometric phase (updated in-place).
 * @amp: Destination amplitude array (can be NULL).
 * @npix: Total number of pixels.
 */
void atmturb_fresnel_extract_unwrapped_scalar(
    const fftwf_complex *restrict field,
    float               *restrict pha,
    float               *restrict amp,
    long                          npix)
{
#ifdef _OPENMP
    #pragma omp parallel for if(!omp_in_parallel() && npix >= 65536) schedule(static)
#endif
    for (long i = 0; i < npix; i++)
    {
        float re = field[i][0];
        float im = field[i][1];
        if (amp != NULL)
        {
            amp[i] = sqrtf(re * re + im * im);
        }

        float phi_geo = pha[i];
        float s_geo, c_geo;
        sincosf(phi_geo, &s_geo, &c_geo);

        float d_re = re * c_geo + im * s_geo;
        float d_im = im * c_geo - re * s_geo;
        float dphi = atan2f(d_im, d_re);
        pha[i]     = phi_geo + dphi;
    }
}
