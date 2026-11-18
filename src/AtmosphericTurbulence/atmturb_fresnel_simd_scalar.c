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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * atmturb_fresnel_phasor_masked_scalar - Scalar fallback for masked pupil complex field
 * @field: Output complex field array.
 * @opd: Input phase or OPD array (NULL for flat wavefront).
 * @mask: Real aperture mask array.
 * @scale: Scale factor converting @opd to radians.
 * @npix: Number of contiguous pixels.
 */
void atmturb_fresnel_phasor_masked_scalar(
    fftwf_complex *restrict field,
    const float   *restrict opd,
    const float   *restrict mask,
    float                   scale,
    long                    npix)
{
    if (opd == NULL)
    {
        for (long i = 0; i < npix; i++)
        {
            field[i][0] = mask[i];
            field[i][1] = 0.0f;
        }
        return;
    }

    for (long i = 0; i < npix; i++)
    {
        float m = mask[i];
        if (m > 0.0f)
        {
            float s, c;
            sincosf(opd[i] * scale, &s, &c);
            field[i][0] = m * c;
            field[i][1] = m * s;
        }
        else
        {
            field[i][0] = 0.0f;
            field[i][1] = 0.0f;
        }
    }
}

/**
 * atmturb_fresnel_spectrum_strip_scalar - Scalar fallback for 8-mode von Karman spectrum strip
 * @out: Output complex spectrum strip (8 elements).
 * @u1: Uniform (0, 1) Box-Muller radial samples (8 floats).
 * @u2: Uniform [0, 1) Box-Muller angular samples (8 floats).
 * @fx0: Starting X spatial frequency index.
 * @spec_p: 4-float parameter tuple {fy2, k0_sq, sqrt_pref, inv_two_km2}.
 */
void atmturb_fresnel_spectrum_strip_scalar(
    fftwf_complex *restrict out,
    const float   *restrict u1,
    const float   *restrict u2,
    float                   fx0,
    const float   *restrict spec_p)
{
    const float fy2         = spec_p[0];
    const float k0_sq       = spec_p[1];
    const float sqrt_pref   = spec_p[2];
    const float inv_two_km2 = spec_p[3];
    const float two_pi      = 2.0f * (float) M_PI;

    for (int lane = 0; lane < 8; lane++)
    {
        float fx    = fx0 + (float) lane;
        float r2    = fx * fx + fy2;
        float dist2 = r2 + k0_sq;
        float r_bm  = sqrtf(-2.0f * logf(u1[lane]));
        float s, c;
        sincosf(two_pi * u2[lane], &s, &c);

        float amp = sqrt_pref * powf(dist2, -11.0f / 12.0f);
        if (inv_two_km2 > 0.0f)
        {
            amp *= expf(-r2 * inv_two_km2);
        }
        out[lane][0] = amp * r_bm * c;
        out[lane][1] = amp * r_bm * s;
    }
}
