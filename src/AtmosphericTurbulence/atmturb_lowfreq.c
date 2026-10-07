// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_lowfreq.c
 * @brief   Analytic Johansson-Gavel subharmonic low-order mode synthesis
 */

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "atmturb_lowfreq.h"
#include "atmturb_simd.h"
#include "atmturb_types.h"

/**
 * atmturb_lowfreq_init_kvec - Calculate spatial frequency wavevectors for one mode
 * @lf: Target low-order modes container.
 * @idx: Mode index (0..23).
 * @p: Subharmonic level (1..3).
 * @nx: Integer harmonic offset along X (-1..1).
 * @ny: Integer harmonic offset along Y (-1..1).
 * @msize: Master screen linear dimension.
 */
static void atmturb_lowfreq_init_kvec(
    atmturb_lowfreq_t *lf,
    int                idx,
    int                p,
    int                nx,
    int                ny,
    long               msize)
{
    double pow3 = pow(3.0, (double) p);
    double k_scale = 2.0 * M_PI / (pow3 * (double) msize);
    lf->kx[idx] = (float) ((double) nx * k_scale);
    lf->ky[idx] = (float) ((double) ny * k_scale);
}

/**
 * atmturb_lowfreq_draw_amplitudes - Draw random complex Gaussian amplitudes for 24 modes
 */
void atmturb_lowfreq_draw_amplitudes(
    float    *are,
    float    *aim,
    long      msize,
    double    r0_pix,
    double    L0_pix,
    double    l0_pix,
    uint64_t  seed)
{
    double r0 = (r0_pix > 0.0) ? r0_pix : pow(6.88, 0.6);
    double k0 = (L0_pix > 0.0) ? ((double) msize / L0_pix) : 0.0;
    double km = (l0_pix > 0.0) ? ((5.92 / (2.0 * M_PI)) * (double) msize / l0_pix) : 0.0;
    double pref = 0.023 * pow(r0, -5.0 / 3.0) * pow((double) msize, 5.0 / 3.0);
    double sqrt_pref = sqrt(pref);
    double inv_two_km2 = (km > 0.0) ? (0.5 / (km * km)) : 0.0;
    double k0_sq = k0 * k0;
    uint64_t rng = atmturb_resolve_seed(seed);

    int idx = 0;
    for (int p = 1; p <= 3; p++)
    {
        double pow3 = pow(3.0, (double) p);
        for (int ny = -1; ny <= 1; ny++)
        {
            for (int nx = -1; nx <= 1; nx++)
            {
                if (nx == 0 && ny == 0)
                {
                    continue;
                }
                double r2 = (double) (nx * nx + ny * ny) / (pow3 * pow3);
                double dist2 = r2 + k0_sq;
                double amp = sqrt_pref * (1.0 / pow3) * pow(dist2, -11.0 / 12.0);
                if (inv_two_km2 > 0.0)
                {
                    amp *= exp(-r2 * inv_two_km2);
                }
                double g0, g1;
                atmturb_rng_gaussian_pair(&rng, &g0, &g1);
                are[idx] = (float) (amp * g0);
                aim[idx] = (float) (amp * g1);
                idx++;
            }
        }
    }
}

/**
 * atmturb_lowfreq_init - Initialize 24 subharmonic Fourier modes for a screen pair
 */
void atmturb_lowfreq_init(
    atmturb_lowfreq_t *lf,
    long               msize,
    double             r0_pix,
    double             L0_pix,
    double             l0_pix,
    uint64_t           seed)
{
    memset(lf, 0, sizeof(*lf));
    int idx = 0;
    for (int p = 1; p <= 3; p++)
    {
        for (int ny = -1; ny <= 1; ny++)
        {
            for (int nx = -1; nx <= 1; nx++)
            {
                if (nx == 0 && ny == 0)
                {
                    continue;
                }
                atmturb_lowfreq_init_kvec(lf, idx++, p, nx, ny, msize);
            }
        }
    }
    lf->nmodes = idx;

    atmturb_lowfreq_draw_amplitudes(lf->are, lf->aim, msize, r0_pix, L0_pix, l0_pix, seed);
    atmturb_lowfreq_draw_amplitudes(lf->bre, lf->bim, msize, r0_pix, L0_pix, l0_pix,
                                    atmturb_rng_stream_seed(seed, 1));
}

/**
 * atmturb_lowfreq_accumulate_custom - Evaluate and accumulate modes with custom amplitudes
 */
void atmturb_lowfreq_accumulate_custom(
    const atmturb_lowfreq_t *lf,
    const float             *are,
    const float             *aim,
    double                   x0,
    double                   y0,
    long                     pup_size,
    long                     os,
    float                    weight,
    float                   *out_pha)
{
    if (lf == NULL || lf->nmodes <= 0 || weight == 0.0f)
    {
        return;
    }
    atmturb_lowfreq_params_t ep = {
        .lf         = lf,
        .screen_idx = 0,
        .custom_are = are,
        .custom_aim = aim,
        .x0         = x0,
        .y0         = y0,
        .pup_size   = pup_size,
        .os         = os,
        .weight     = weight,
        .out_pha    = out_pha
    };
    atmturb_extrude_lowfreq(&ep);
}

/**
 * atmturb_lowfreq_accumulate - Analytically evaluate and accumulate low-order modes
 */
void atmturb_lowfreq_accumulate(
    const atmturb_lowfreq_t *lf,
    int                      screen_idx,
    double                   x0,
    double                   y0,
    long                     pup_size,
    long                     os,
    float                    weight,
    float                   *out_pha)
{
    if (lf == NULL)
    {
        return;
    }
    const float *are = (screen_idx == 1) ? lf->bre : lf->are;
    const float *aim = (screen_idx == 1) ? lf->bim : lf->aim;
    atmturb_lowfreq_accumulate_custom(lf, are, aim, x0, y0, pup_size, os, weight, out_pha);
}
