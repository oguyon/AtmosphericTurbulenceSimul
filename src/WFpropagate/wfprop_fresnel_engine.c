// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    wfprop_fresnel_engine.c
 * @brief   Thread-safe 2D Fresnel diffractive propagation engine implementation
 */

#define _GNU_SOURCE
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "wfprop_fresnel_engine.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

/**
 * wfprop_fresnel_engine_init - Initialize thread-local Fresnel propagation engine
 * @eng: Pointer to engine structure to initialize.
 * @n: Grid dimension in pixels.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int wfprop_fresnel_engine_init(
    wfprop_fresnel_engine_t *eng,
    long                     n)
{
    if (eng == NULL || n <= 0)
    {
        return -1;
    }

    eng->n        = n;
    long ntot     = n * n;
    eng->buf      = (fftwf_complex *) fftwf_alloc_complex((size_t) ntot);
    if (eng->buf == NULL)
    {
        return -1;
    }

    eng->fwd_plan = fftwf_plan_dft_2d((int) n, (int) n, eng->buf, eng->buf,
                                      FFTW_FORWARD, FFTW_ESTIMATE);
    eng->inv_plan = fftwf_plan_dft_2d((int) n, (int) n, eng->buf, eng->buf,
                                      FFTW_BACKWARD, FFTW_ESTIMATE);

    if (eng->fwd_plan == NULL || eng->inv_plan == NULL)
    {
        wfprop_fresnel_engine_free(eng);
        return -1;
    }

    return 0;
}

/**
 * wfprop_fresnel_engine_free - Release resources held by Fresnel propagation engine
 * @eng: Pointer to engine structure to tear down.
 */
void wfprop_fresnel_engine_free(
    wfprop_fresnel_engine_t *eng)
{
    if (eng == NULL)
    {
        return;
    }

    if (eng->fwd_plan != NULL)
    {
        fftwf_destroy_plan(eng->fwd_plan);
        eng->fwd_plan = NULL;
    }

    if (eng->inv_plan != NULL)
    {
        fftwf_destroy_plan(eng->inv_plan);
        eng->inv_plan = NULL;
    }

    if (eng->buf != NULL)
    {
        fftwf_free(eng->buf);
        eng->buf = NULL;
    }

    eng->n = 0;
}

/**
 * wfprop_fresnel_tf_build - Precompute 2D angular spectrum transfer function
 * @tf: Output complex transfer function array (n * n elements).
 * @n: Grid dimension in pixels.
 * @pixscale_m: Grid physical sampling scale [m/pixel].
 * @z_m: Propagation distance along optical path [m].
 * @lambda_m: Optical wavelength [m].
 * @cutoff_rad: Frequency mask cutoff radius in pixels (<= 0 for Nyquist disc).
 *
 * Return: 0 on success, -1 on invalid argument.
 */
int wfprop_fresnel_tf_build(
    fftwf_complex *tf,
    long           n,
    double         pixscale_m,
    double         z_m,
    double         lambda_m,
    double         cutoff_rad)
{
    if (tf == NULL || n <= 0 || pixscale_m <= 0.0 || lambda_m <= 0.0)
    {
        return -1;
    }

    double l_grid   = (double) n * pixscale_m;
    double coeff    = M_PI * z_m * lambda_m / (l_grid * l_grid);
    double inv_norm = 1.0 / (double) (n * n);

    #pragma omp parallel for schedule(static)
    for (long j = 0; j < n; j++)
    {
        long fy = (j < n / 2) ? j : (j - n);
        double fy2 = (double) (fy * fy);
        long row = j * n;

        for (long i = 0; i < n; i++)
        {
            long fx = (i < n / 2) ? i : (i - n);
            double sqdist = (double) (fx * fx) + fy2;
            double pha = -coeff * sqdist;
            double amp = 1.0;

            if (cutoff_rad > 0.0)
            {
                double dist = sqrt(sqdist);
                if (dist >= cutoff_rad)
                {
                    amp = 0.0;
                }
                else if (dist >= cutoff_rad - 1.0)
                {
                    amp = 0.5 * (1.0 + cos(M_PI * (dist - cutoff_rad + 1.0)));
                }
            }

            tf[row + i][0] = (float) (amp * cos(pha) * inv_norm);
            tf[row + i][1] = (float) (amp * sin(pha) * inv_norm);
        }
    }

    return 0;
}

/**
 * wfprop_fresnel_engine_apply - Apply optical diffraction step using precomputed transfer function
 * @eng: Initialized Fresnel engine.
 * @field: In-out 2D complex optical field (n * n elements).
 * @tf: Precomputed 2D complex transfer function (n * n elements).
 */
void wfprop_fresnel_engine_apply(
    wfprop_fresnel_engine_t *eng,
    fftwf_complex           *field,
    const fftwf_complex     *tf)
{
    long ntot = eng->n * eng->n;

    memcpy(eng->buf, field, sizeof(fftwf_complex) * (size_t) ntot);

    fftwf_execute(eng->fwd_plan);

    for (long k = 0; k < ntot; k++)
    {
        float r_b = eng->buf[k][0];
        float i_b = eng->buf[k][1];
        float r_t = tf[k][0];
        float i_t = tf[k][1];

        eng->buf[k][0] = r_b * r_t - i_b * i_t;
        eng->buf[k][1] = r_b * i_t + i_b * r_t;
    }

    fftwf_execute(eng->inv_plan);

    memcpy(field, eng->buf, sizeof(fftwf_complex) * (size_t) ntot);
}
