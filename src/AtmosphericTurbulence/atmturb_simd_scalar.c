// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_simd_scalar.c
 * @brief   Scalar reference kernels for atmospheric turbulence simulation
 */

#include <math.h>
#include "atmturb_simd.h"

/**
 * atmturb_extrude_accumulate_scalar - Scalar bilinear extrusion and accumulation
 * @master: Master screen array.
 * @msize: Master screen dimension.
 * @x0: Sub-pixel X coordinate offset.
 * @y0: Sub-pixel Y coordinate offset.
 * @pup_size: Pupil dimension.
 * @weight: Layer weight.
 * @out_pha: Output phase array.
 */
void atmturb_extrude_accumulate_scalar(const float *master, long msize, double x0, double y0,
                                       long pup_size, float weight, float *out_pha)
{
    float fx = (float)(x0 - floor(x0));
    float fy = (float)(y0 - floor(y0));

    long base_x = (long)floor(x0);
    long base_y = (long)floor(y0);

    float w00 = (1.0f - fx) * (1.0f - fy) * weight;
    float w10 = fx * (1.0f - fy) * weight;
    float w01 = (1.0f - fx) * fy * weight;
    float w11 = fx * fy * weight;

    for (long jj = 0; jj < pup_size; jj++)
    {
        long iy0 = (base_y + jj) % msize;
        if (iy0 < 0) iy0 += msize;
        long iy1 = (iy0 + 1) % msize;

        const float *row0 = &master[iy0 * msize];
        const float *row1 = &master[iy1 * msize];
        float *out_row = &out_pha[jj * pup_size];

        for (long ii = 0; ii < pup_size; ii++)
        {
            long ix0 = (base_x + ii) % msize;
            if (ix0 < 0) ix0 += msize;
            long ix1 = (ix0 + 1) % msize;

            float val = w00 * row0[ix0] + w10 * row0[ix1] +
                        w01 * row1[ix0] + w11 * row1[ix1];
            out_row[ii] += val;
        }
    }
}

/**
 * atmturb_scale_float_array_scalar - Scalar float array scaling
 * @dest: Destination array.
 * @src: Source array.
 * @scale: Scaling multiplier.
 * @n: Array length.
 */
void atmturb_scale_float_array_scalar(float *dest, const float *src, float scale, long n)
{
    for (long i = 0; i < n; i++)
    {
        dest[i] = src[i] * scale;
    }
}

/**
 * atmturb_init_phase_amp_scalar - Scalar phase and amplitude initialization
 * @pha: Phase array (set to 0.0).
 * @amp: Amplitude array (set to 1.0).
 * @n: Array length.
 */
void atmturb_init_phase_amp_scalar(float *pha, float *amp, long n)
{
    for (long i = 0; i < n; i++)
    {
        pha[i] = 0.0f;
        amp[i] = 1.0f;
    }
}
