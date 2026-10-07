// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_simd_scalar.c
 * @brief   Scalar reference kernels for atmospheric turbulence simulation
 */

#include <math.h>
#include "atmturb_lowfreq.h"
#include "atmturb_simd.h"

/**
 * atmturb_keys_weights - Keys bicubic convolution tap weights for a = -0.5
 * @f: Fractional sub-pixel position in [0, 1).
 * @w: Output 4-element tap weight array for offsets [-1, 0, 1, 2].
 */
static inline void atmturb_keys_weights(
    float f,
    float w[4])
{
    float f2 = f * f;
    float f3 = f2 * f;
    w[0] = -0.5f * f + f2 - 0.5f * f3;
    w[1] = 1.0f - 2.5f * f2 + 1.5f * f3;
    w[2] = 0.5f * f + 2.0f * f2 - 1.5f * f3;
    w[3] = -0.5f * f2 + 0.5f * f3;
}

/**
 * atmturb_extrude_accumulate_bilinear_scalar - Scalar bilinear extrusion
 * @params: Extrusion configuration and data pointers.
 */
void atmturb_extrude_accumulate_bilinear_scalar(
    const atmturb_extrude_params_t *params)
{
    long msize = params->msize;
    long pup_size = params->pup_size;
    long os = (params->os > 1) ? params->os : 1;
    float weight = params->weight;
    const float *master = params->master;
    float *out_pha = params->out_pha;

    float fx = (float)(params->x0 - floor(params->x0));
    float fy = (float)(params->y0 - floor(params->y0));

    long base_x = (long)floor(params->x0);
    long base_y = (long)floor(params->y0);

    float w00 = (1.0f - fx) * (1.0f - fy) * weight;
    float w10 = fx * (1.0f - fy) * weight;
    float w01 = (1.0f - fx) * fy * weight;
    float w11 = fx * fy * weight;

    for (long jj = 0; jj < pup_size; jj++)
    {
        long iy0 = (base_y + jj * os) % msize;
        if (iy0 < 0)
        {
            iy0 += msize;
        }
        long iy1 = (iy0 + 1) % msize;

        const float *row0 = &master[iy0 * msize];
        const float *row1 = &master[iy1 * msize];
        float *out_row = &out_pha[jj * pup_size];

        for (long ii = 0; ii < pup_size; ii++)
        {
            long ix0 = (base_x + ii * os) % msize;
            if (ix0 < 0)
            {
                ix0 += msize;
            }
            long ix1 = (ix0 + 1) % msize;

            float val = w00 * row0[ix0] + w10 * row0[ix1] +
                        w01 * row1[ix0] + w11 * row1[ix1];
            out_row[ii] += val;
        }
    }
}

/**
 * atmturb_extrude_accumulate_bicubic_scalar - Scalar Keys bicubic extrusion
 * @params: Extrusion configuration and data pointers.
 */
void atmturb_extrude_accumulate_bicubic_scalar(
    const atmturb_extrude_params_t *params)
{
    long msize = params->msize;
    long pup_size = params->pup_size;
    long os = (params->os > 1) ? params->os : 1;
    float weight = params->weight;
    const float *master = params->master;
    float *out_pha = params->out_pha;

    float fx = (float)(params->x0 - floor(params->x0));
    float fy = (float)(params->y0 - floor(params->y0));

    long base_x = (long)floor(params->x0);
    long base_y = (long)floor(params->y0);

    float wx[4];
    float wy[4];
    atmturb_keys_weights(fx, wx);
    atmturb_keys_weights(fy, wy);

    float wy_wt[4];
    for (int n = 0; n < 4; n++)
    {
        wy_wt[n] = wy[n] * weight;
    }

    for (long jj = 0; jj < pup_size; jj++)
    {
        long by = base_y + jj * os;
        long iy[4];
        for (int n = 0; n < 4; n++)
        {
            long y = (by + n - 1) % msize;
            if (y < 0)
            {
                y += msize;
            }
            iy[n] = y;
        }

        const float *row0 = &master[iy[0] * msize];
        const float *row1 = &master[iy[1] * msize];
        const float *row2 = &master[iy[2] * msize];
        const float *row3 = &master[iy[3] * msize];
        float *out_row = &out_pha[jj * pup_size];

        for (long ii = 0; ii < pup_size; ii++)
        {
            long bx = base_x + ii * os;
            long ix[4];
            for (int m = 0; m < 4; m++)
            {
                long x = (bx + m - 1) % msize;
                if (x < 0)
                {
                    x += msize;
                }
                ix[m] = x;
            }

            float h0 = wx[0] * row0[ix[0]] + wx[1] * row0[ix[1]] +
                       wx[2] * row0[ix[2]] + wx[3] * row0[ix[3]];
            float h1 = wx[0] * row1[ix[0]] + wx[1] * row1[ix[1]] +
                       wx[2] * row1[ix[2]] + wx[3] * row1[ix[3]];
            float h2 = wx[0] * row2[ix[0]] + wx[1] * row2[ix[1]] +
                       wx[2] * row2[ix[2]] + wx[3] * row2[ix[3]];
            float h3 = wx[0] * row3[ix[0]] + wx[1] * row3[ix[1]] +
                       wx[2] * row3[ix[2]] + wx[3] * row3[ix[3]];

            out_row[ii] += wy_wt[0] * h0 + wy_wt[1] * h1 +
                           wy_wt[2] * h2 + wy_wt[3] * h3;
        }
    }
}

/**
 * atmturb_extrude_accumulate_scalar - Scalar extrusion dispatcher (bilinear/bicubic)
 * @params: Extrusion configuration and data pointers.
 */
void atmturb_extrude_accumulate_scalar(
    const atmturb_extrude_params_t *params)
{
    if (params->interp == ATMTURB_INTERP_BICUBIC)
    {
        atmturb_extrude_accumulate_bicubic_scalar(params);
    }
    else
    {
        atmturb_extrude_accumulate_bilinear_scalar(params);
    }
}

/**
 * atmturb_scale_float_array_scalar - Scalar float array scaling
 * @dest: Destination array.
 * @src: Source array.
 * @scale: Scaling multiplier.
 * @n: Array length.
 */
void atmturb_scale_float_array_scalar(
    float       *dest,
    const float *src,
    float        scale,
    long         n)
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
void atmturb_init_phase_amp_scalar(
    float *pha,
    float *amp,
    long   n)
{
    for (long i = 0; i < n; i++)
    {
        pha[i] = 0.0f;
        amp[i] = 1.0f;
    }
}

/**
 * atmturb_lowfreq_mode_accumulate_scalar - Accumulate one harmonic mode into pupil
 * @params: Low-frequency configuration bundle.
 * @kx: Mode spatial frequency along X [rad/master px].
 * @ky: Mode spatial frequency along Y [rad/master px].
 * @amp_re: Mode real amplitude.
 * @amp_im: Mode imaginary amplitude.
 */
static void atmturb_lowfreq_mode_accumulate_scalar(
    const atmturb_lowfreq_params_t *params,
    float                           kx,
    float                           ky,
    float                           amp_re,
    float                           amp_im)
{
    float h_re[1024], h_im[1024];
    long n = (params->pup_size <= 1024) ? params->pup_size : 1024;
    for (long i = 0; i < n; i++)
    {
        double th_x = (double) kx * (params->x0 + (double) (i * params->os));
        h_re[i] = (float) cos(th_x);
        h_im[i] = (float) sin(th_x);
    }

    for (long j = 0; j < params->pup_size; j++)
    {
        double th_y = (double) ky * (params->y0 + (double) (j * params->os));
        float vy_re = (float) cos(th_y);
        float vy_im = (float) sin(th_y);
        float c_re = (amp_re * vy_re - amp_im * vy_im) * params->weight;
        float c_im = (amp_re * vy_im + amp_im * vy_re) * params->weight;

        long row = j * params->pup_size;
        for (long i = 0; i < n; i++)
        {
            params->out_pha[row + i] += c_re * h_re[i] - c_im * h_im[i];
        }
    }
}

/**
 * atmturb_extrude_lowfreq_scalar - Scalar separable low-order mode accumulation
 * @params: Low-frequency configuration and data pointers.
 */
void atmturb_extrude_lowfreq_scalar(
    const atmturb_lowfreq_params_t *params)
{
    const atmturb_lowfreq_t *lf = (const atmturb_lowfreq_t *) params->lf;
    const float *are = params->custom_are;
    const float *aim = params->custom_aim;
    if (are == NULL || aim == NULL)
    {
        int is_s1 = (params->screen_idx == 1);
        are = is_s1 ? lf->bre : lf->are;
        aim = is_s1 ? lf->bim : lf->aim;
    }

    for (int m = 0; m < lf->nmodes; m++)
    {
        atmturb_lowfreq_mode_accumulate_scalar(params, lf->kx[m], lf->ky[m], are[m], aim[m]);
    }
}

