// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_boiling_kernels.c
 * @brief   SIMD/OpenMP Fourier AR(1)/AR(2) and subharmonic compute kernels for boiling
 */

#define _GNU_SOURCE
#include <math.h>
#include <stdint.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "atmturb_boiling.h"
#include "atmturb_fresnel_simd.h"
#include "atmturb_types.h"

/**
 * atmturb_boil_omr2 - Evaluate 1 - exp(-2 * eta) without FP32 cancellation at small eta
 */
static inline float atmturb_boil_omr2(
    float eta,
    float r2_pole)
{
    float v = (eta < 0.02f) ? (2.0f * eta * (1.0f - eta * (1.0f - 0.6666667f * eta)))
                            : (1.0f - r2_pole);
    return (v > 0.0f) ? v : 0.0f;
}

/**
 * atmturb_boil_draw_strip8 - Counter-unrolled 8-lane SplitMix64 uniform pair draw
 */
static inline void atmturb_boil_draw_strip8(
    uint64_t *rng_io,
    float     u1[8],
    float     u2[8])
{
    const uint64_t gamma    = 0x9e3779b97f4a7c15ULL;
    const double   inv_2p53 = 1.0 / 9007199254740992.0;
    uint64_t       s        = *rng_io;

    for (int lane = 0; lane < 8; lane++)
    {
        uint64_t z1 = s + (uint64_t) (2 * lane + 1) * gamma;
        uint64_t z2 = s + (uint64_t) (2 * lane + 2) * gamma;
        z1 = (z1 ^ (z1 >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z1 = (z1 ^ (z1 >> 27)) * 0x94d049bb133111ebULL;
        z1 = z1 ^ (z1 >> 31);
        z2 = (z2 ^ (z2 >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z2 = (z2 ^ (z2 >> 27)) * 0x94d049bb133111ebULL;
        z2 = z2 ^ (z2 >> 31);
        float f1 = (float) (((double) (z1 >> 11) + 1.0) * inv_2p53);
        u1[lane] = (f1 < 0.99999994f) ? f1 : 0.99999994f;
        u2[lane] = (float) ((double) (z2 >> 11) * inv_2p53);
    }
    *rng_io = s + 16ULL * gamma;
}

/**
 * atmturb_boil_mode_amp - Compute von Karman mode standard deviation at squared radius r2
 */
static inline float atmturb_boil_mode_amp(
    float r2,
    float k0_sq,
    float sqrt_pref,
    float inv_two_km2)
{
    float amp = sqrt_pref * powf(r2 + k0_sq, -11.0f / 12.0f);
    if (inv_two_km2 > 0.0f)
    {
        amp *= expf(-r2 * inv_two_km2);
    }
    return amp;
}

/**
 * atmturb_boil_apply_init - Stationary t=0 AR(1) or cascaded AR(2) mode initialization
 */
static inline void atmturb_boil_apply_init(
    fftwf_complex       *restrict cur,
    fftwf_complex       *restrict aux,
    fftwf_complex       *restrict wrk,
    const fftwf_complex           z0,
    const fftwf_complex           z1,
    float                         eta,
    int                           ar_order)
{
    if (ar_order == 2 && aux != NULL)
    {
        float rp       = expf(-eta);
        float rp2      = rp * rp;
        float omr2     = atmturb_boil_omr2(eta, rp2);
        float inv_1pr2 = 1.0f / (1.0f + rp2);
        float c0       = 2.0f * rp * inv_1pr2;
        float c1       = omr2 * inv_1pr2;
        float cre      = c0 * z0[0] + c1 * z1[0];
        float cim      = c0 * z0[1] + c1 * z1[1];
        (*cur)[0] = cre;
        (*cur)[1] = cim;
        (*aux)[0] = c1 * (rp * z0[0] + z1[0]);
        (*aux)[1] = c1 * (rp * z0[1] + z1[1]);
        (*wrk)[0] = cre;
        (*wrk)[1] = cim;
    }
    else
    {
        (*cur)[0] = z0[0];
        (*cur)[1] = z0[1];
        (*wrk)[0] = z0[0];
        (*wrk)[1] = z0[1];
    }
}

/**
 * atmturb_boil_apply_step - Advance one Fourier mode by one AR(1) or cascaded AR(2) step
 */
static inline void atmturb_boil_apply_step(
    fftwf_complex       *restrict cur,
    fftwf_complex       *restrict aux,
    fftwf_complex       *restrict wrk,
    const fftwf_complex           wn,
    float                         eta,
    int                           ar_order)
{
    float rp   = expf(-eta);
    float rp2  = rp * rp;
    float omr2 = atmturb_boil_omr2(eta, rp2);

    if (ar_order == 2 && aux != NULL)
    {
        float b   = omr2 * sqrtf(omr2 / (1.0f + rp2));
        float ure = rp * (*aux)[0] + b * wn[0];
        float uim = rp * (*aux)[1] + b * wn[1];
        float cre = rp * (*cur)[0] + ure;
        float cim = rp * (*cur)[1] + uim;
        (*aux)[0] = ure;
        (*aux)[1] = uim;
        (*cur)[0] = cre;
        (*cur)[1] = cim;
        (*wrk)[0] = cre;
        (*wrk)[1] = cim;
    }
    else
    {
        float b   = sqrtf(omr2);
        float cre = rp * (*cur)[0] + b * wn[0];
        float cim = rp * (*cur)[1] + b * wn[1];
        (*cur)[0] = cre;
        (*cur)[1] = cim;
        (*wrk)[0] = cre;
        (*wrk)[1] = cim;
    }
}

/**
 * atmturb_boil_init_row - Initialize interior Hermitian columns 1..half_end-1 of row jj
 */
static void atmturb_boil_init_row(
    atmturb_boil_layer_t *bl,
    fftwf_complex        *c2r_work,
    long                  msize,
    long                  jj,
    int                   ar_order,
    const uint64_t        seeds[2])
{
    const long     hcols    = msize / 2 + 1;
    const long     half_end = (msize % 2 == 0) ? (msize / 2) : hcols;
    const long     fy       = (jj <= msize / 2) ? jj : (jj - msize);
    const float    fy2      = (float) (fy * fy);
    const float    sp[4]    = {fy2, bl->k0_sq, bl->sqrt_pref_pair, bl->inv_two_km2};
    uint64_t       rng0     = atmturb_rng_stream_seed(seeds[0], (uint64_t) jj);
    uint64_t       rng1     = atmturb_rng_stream_seed(seeds[1], (uint64_t) jj);
    fftwf_complex *rcur     = &bl->spec_cur[jj * hcols];
    fftwf_complex *raux     = (bl->spec_aux != NULL) ? &bl->spec_aux[jj * hcols] : NULL;
    fftwf_complex *rwrk     = &c2r_work[jj * hcols];

    for (long ii = 1; ii < half_end; ii++)
    {
        if (ii + 8 <= half_end && (fy2 >= 256.0f || ii >= 16))
        {
            float         u1a[8], u2a[8], u1b[8], u2b[8];
            fftwf_complex z0[8], z1[8];
            atmturb_boil_draw_strip8(&rng0, u1a, u2a);
            atmturb_fresnel_simd_spectrum_strip(z0, u1a, u2a, (float) ii, sp);
            if (ar_order == 2)
            {
                atmturb_boil_draw_strip8(&rng1, u1b, u2b);
                atmturb_fresnel_simd_spectrum_strip(z1, u1b, u2b, (float) ii, sp);
            }
            for (int l = 0; l < 8; l++)
            {
                float fx  = (float) (ii + l);
                float eta = bl->alpha_dt * sqrtf(fx * fx + fy2);
                atmturb_boil_apply_init(&rcur[ii + l], raux ? &raux[ii + l] : NULL,
                                        &rwrk[ii + l], z0[l], z1[l], eta, ar_order);
            }
            ii += 7;
            continue;
        }

        float  r2  = (float) (ii * ii) + fy2;
        float  amp = atmturb_boil_mode_amp(r2, sp[1], sp[2], sp[3]);
        double g0, g1, h0 = 0.0, h1 = 0.0;
        atmturb_rng_gaussian_pair(&rng0, &g0, &g1);
        if (ar_order == 2)
        {
            atmturb_rng_gaussian_pair(&rng1, &h0, &h1);
        }
        fftwf_complex z0 = {amp * (float) g0, amp * (float) g1};
        fftwf_complex z1 = {amp * (float) h0, amp * (float) h1};
        atmturb_boil_apply_init(&rcur[ii], raux ? &raux[ii] : NULL, &rwrk[ii],
                                z0, z1, bl->alpha_dt * sqrtf(r2), ar_order);
    }
}

/**
 * atmturb_boil_step_row - Advance interior Hermitian columns 1..half_end-1 of row jj
 */
static void atmturb_boil_step_row(
    atmturb_boil_layer_t *bl,
    fftwf_complex        *c2r_work,
    long                  msize,
    long                  jj,
    int                   ar_order,
    uint64_t              step_seed)
{
    const long     hcols    = msize / 2 + 1;
    const long     half_end = (msize % 2 == 0) ? (msize / 2) : hcols;
    const long     fy       = (jj <= msize / 2) ? jj : (jj - msize);
    const float    fy2      = (float) (fy * fy);
    const float    sp[4]    = {fy2, bl->k0_sq, bl->sqrt_pref_pair, bl->inv_two_km2};
    uint64_t       rng      = atmturb_rng_stream_seed(step_seed, (uint64_t) jj);
    fftwf_complex *rcur     = &bl->spec_cur[jj * hcols];
    fftwf_complex *raux     = (bl->spec_aux != NULL) ? &bl->spec_aux[jj * hcols] : NULL;
    fftwf_complex *rwrk     = &c2r_work[jj * hcols];

    for (long ii = 1; ii < half_end; ii++)
    {
        if (ii + 8 <= half_end && (fy2 >= 256.0f || ii >= 16))
        {
            float         u1[8], u2[8];
            fftwf_complex wn[8];
            atmturb_boil_draw_strip8(&rng, u1, u2);
            atmturb_fresnel_simd_spectrum_strip(wn, u1, u2, (float) ii, sp);
            for (int l = 0; l < 8; l++)
            {
                float fx  = (float) (ii + l);
                float eta = bl->alpha_dt * sqrtf(fx * fx + fy2);
                atmturb_boil_apply_step(&rcur[ii + l], raux ? &raux[ii + l] : NULL,
                                        &rwrk[ii + l], wn[l], eta, ar_order);
            }
            ii += 7;
            continue;
        }

        float  r2  = (float) (ii * ii) + fy2;
        float  amp = atmturb_boil_mode_amp(r2, sp[1], sp[2], sp[3]);
        double g0, g1;
        atmturb_rng_gaussian_pair(&rng, &g0, &g1);
        fftwf_complex wn = {amp * (float) g0, amp * (float) g1};
        atmturb_boil_apply_step(&rcur[ii], raux ? &raux[ii] : NULL, &rwrk[ii],
                                wn, bl->alpha_dt * sqrtf(r2), ar_order);
    }
}

/**
 * atmturb_boil_axis_col - Update self-conjugate column ii in {0, msize/2} and mirror rows
 */
static void atmturb_boil_axis_col(
    atmturb_boil_layer_t *bl,
    fftwf_complex        *c2r_work,
    long                  msize,
    long                  ii,
    int                   ar_order,
    int                   is_init,
    const uint64_t        seeds[2])
{
    const long hcols   = msize / 2 + 1;
    const long half_lim = msize / 2;
    const int  col_tag = (ii == 0) ? 0 : 1;
    uint64_t   rng0    = atmturb_rng_stream_seed(seeds[0], (uint64_t) (msize + col_tag));
    uint64_t   rng1    = atmturb_rng_stream_seed(seeds[1], (uint64_t) (msize + col_tag));

    for (long jj = 0; jj <= half_lim; jj++)
    {
        long idx = jj * hcols + ii;
        if (ii == 0 && jj == 0)
        {
            bl->spec_cur[idx][0] = 0.0f;
            bl->spec_cur[idx][1] = 0.0f;
            c2r_work[idx][0]     = 0.0f;
            c2r_work[idx][1]     = 0.0f;
            continue;
        }

        int    is_nyq = (jj == 0) || ((msize % 2 == 0) && (jj == half_lim));
        float  r2     = (float) (ii * ii + jj * jj);
        float  spref  = is_nyq ? bl->sqrt_pref_real : bl->sqrt_pref_pair;
        float  amp    = atmturb_boil_mode_amp(r2, bl->k0_sq, spref, bl->inv_two_km2);
        double g0, g1, h0 = 0.0, h1 = 0.0;
        atmturb_rng_gaussian_pair(&rng0, &g0, &g1);
        fftwf_complex z0 = {amp * (float) g0, is_nyq ? 0.0f : (amp * (float) g1)};
        fftwf_complex *aux_p = bl->spec_aux ? &bl->spec_aux[idx] : NULL;
        float eta = bl->alpha_dt * sqrtf(r2);

        if (is_init)
        {
            if (ar_order == 2)
            {
                atmturb_rng_gaussian_pair(&rng1, &h0, &h1);
            }
            fftwf_complex z1 = {amp * (float) h0, is_nyq ? 0.0f : (amp * (float) h1)};
            atmturb_boil_apply_init(&bl->spec_cur[idx], aux_p, &c2r_work[idx],
                                    z0, z1, eta, ar_order);
        }
        else
        {
            atmturb_boil_apply_step(&bl->spec_cur[idx], aux_p, &c2r_work[idx],
                                    z0, eta, ar_order);
        }

        if (!is_nyq && jj > 0)
        {
            long mir = (msize - jj) * hcols + ii;
            bl->spec_cur[mir][0] = c2r_work[idx][0];
            bl->spec_cur[mir][1] = -c2r_work[idx][1];
            c2r_work[mir][0]     = c2r_work[idx][0];
            c2r_work[mir][1]     = -c2r_work[idx][1];
        }
    }
}

/**
 * atmturb_boil_init_spectrum_cpu - Stationary t=0 half-spectrum draw across all modes
 */
void atmturb_boil_init_spectrum_cpu(
    atmturb_boil_layer_t *bl,
    fftwf_complex        *c2r_work,
    long                  msize,
    int                   ar_order,
    const uint64_t        seeds[2])
{
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (long jj = 0; jj < msize; jj++)
    {
        atmturb_boil_init_row(bl, c2r_work, msize, jj, ar_order, seeds);
    }

    atmturb_boil_axis_col(bl, c2r_work, msize, 0, ar_order, 1, seeds);
    if (msize % 2 == 0)
    {
        atmturb_boil_axis_col(bl, c2r_work, msize, msize / 2, ar_order, 1, seeds);
    }
}

/**
 * atmturb_boil_step_spectrum_cpu - Advance half-spectrum by one AR(1)/AR(2) keyframe step
 */
void atmturb_boil_step_spectrum_cpu(
    atmturb_boil_layer_t *bl,
    fftwf_complex        *c2r_work,
    long                  msize,
    int                   ar_order,
    uint64_t              step_seed)
{
    uint64_t seeds[2] = {step_seed, 0ULL};

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (long jj = 0; jj < msize; jj++)
    {
        atmturb_boil_step_row(bl, c2r_work, msize, jj, ar_order, step_seed);
    }

    atmturb_boil_axis_col(bl, c2r_work, msize, 0, ar_order, 0, seeds);
    if (msize % 2 == 0)
    {
        atmturb_boil_axis_col(bl, c2r_work, msize, msize / 2, ar_order, 0, seeds);
    }
}

/**
 * atmturb_boil_compute_rho1 - Compute PSD-weighted 1-step keyframe correlation rho_bar
 */
double atmturb_boil_compute_rho1(
    const atmturb_boil_layer_t *bl,
    long                        msize,
    int                         ar_order)
{
    const long hcols   = msize / 2 + 1;
    double     num_sum = 0.0;
    double     den_sum = 0.0;

#ifdef _OPENMP
    #pragma omp parallel for reduction(+:num_sum, den_sum) schedule(static)
#endif
    for (long jj = 0; jj < msize; jj++)
    {
        long   fy      = (jj <= msize / 2) ? jj : (jj - msize);
        double fy2     = (double) (fy * fy);
        double row_num = 0.0;
        double row_den = 0.0;

        for (long ii = 0; ii < hcols; ii++)
        {
            if (ii == 0 && jj == 0)
            {
                continue;
            }
            int is_axis = (ii == 0) || ((msize % 2 == 0) && (ii == msize / 2));
            if (is_axis && jj > msize / 2)
            {
                continue;
            }
            int    is_nyq = is_axis && ((jj == 0) || ((msize % 2 == 0) && (jj == msize / 2)));
            double mult   = is_nyq ? 1.0 : 2.0;
            double r2     = (double) (ii * ii) + fy2;
            double w_psd  = mult * pow(r2 + (double) bl->k0_sq, -11.0 / 6.0);
            if (bl->inv_two_km2 > 0.0f)
            {
                w_psd *= exp(-2.0 * r2 * (double) bl->inv_two_km2);
            }
            double rp   = exp(-(double) bl->alpha_dt * sqrt(r2));
            double rho1 = (ar_order == 2) ? ((2.0 * rp) / (1.0 + rp * rp)) : rp;
            row_num += w_psd * rho1;
            row_den += w_psd;
        }
        num_sum += row_num;
        den_sum += row_den;
    }

    return (den_sum > 0.0) ? (num_sum / den_sum) : 1.0;
}

/**
 * atmturb_boil_init_lowfreq_layer - Stationary t=0 AR(1)/AR(2) setup for 24 subharmonics
 */
void atmturb_boil_init_lowfreq_layer(
    atmturb_boil_layer_t         *bl,
    const atmturb_lowfreq_t      *lf_base,
    const atmturb_boil_lf_spec_t *spec)
{
    double   r0    = (spec->r0_pix > 0.0) ? spec->r0_pix : pow(6.88, 0.6);
    double   k0    = (spec->L0_pix > 0.0) ? ((double) spec->msize / spec->L0_pix) : 0.0;
    double   km    = (spec->l0_pix > 0.0)
                         ? ((5.92 / (2.0 * M_PI)) * (double) spec->msize / spec->l0_pix) : 0.0;
    double   spref = sqrt(0.023 * pow(r0, -5.0 / 3.0) * pow((double) spec->msize, 5.0 / 3.0));
    double   ikm2  = (km > 0.0) ? (0.5 / (km * km)) : 0.0;
    uint64_t rng0  = atmturb_resolve_seed(spec->seed0);
    uint64_t rng1  = atmturb_resolve_seed(spec->seed1);

    int m = 0;
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
                double r2  = (double) (nx * nx + ny * ny) / (pow3 * pow3);
                double amp = spref * (1.0 / pow3) * pow(r2 + k0 * k0, -11.0 / 12.0);
                if (ikm2 > 0.0)
                {
                    amp *= exp(-r2 * ikm2);
                }
                double kx   = (double) lf_base->kx[m];
                double ky   = (double) lf_base->ky[m];
                float  eta  = (float) (spec->gamma_scale * sqrt(kx * kx + ky * ky));
                float  rp   = expf(-eta);
                float  rp2  = rp * rp;
                float  omr2 = atmturb_boil_omr2(eta, rp2);

                bl->lf_amp[m]  = (float) amp;
                bl->lf_r[m]    = rp;
                bl->lf_omr2[m] = omr2;
                bl->lf_b[m]    = (spec->ar_order == 2)
                                     ? (omr2 * sqrtf(omr2 / (1.0f + rp2))) : sqrtf(omr2);

                double g0, g1, h0 = 0.0, h1 = 0.0;
                atmturb_rng_gaussian_pair(&rng0, &g0, &g1);
                if (spec->ar_order == 2)
                {
                    atmturb_rng_gaussian_pair(&rng1, &h0, &h1);
                }
                fftwf_complex z0 = {(float) (amp * g0), (float) (amp * g1)};
                fftwf_complex z1 = {(float) (amp * h0), (float) (amp * h1)};
                fftwf_complex cur, aux = {0.0f, 0.0f}, wrk;
                atmturb_boil_apply_init(&cur, &aux, &wrk, z0, z1, eta, spec->ar_order);
                bl->lf_cur_re[m]  = cur[0];
                bl->lf_cur_im[m]  = cur[1];
                bl->lf_aux_re[m]  = aux[0];
                bl->lf_aux_im[m]  = aux[1];
                bl->lf_ring_re[m] = cur[0];
                bl->lf_ring_im[m] = cur[1];
                m++;
            }
        }
    }
    bl->max_lf_ready = 0;
}

/**
 * atmturb_boil_step_lowfreq_layer - Advance 24 subharmonic modes by one frame step
 */
void atmturb_boil_step_lowfreq_layer(
    atmturb_boil_layer_t *bl,
    int                   ar_order,
    long                  frame_t,
    uint64_t              step_seed)
{
    uint64_t rng  = atmturb_resolve_seed(step_seed);
    long     slot = (frame_t % bl->lf_cap) * ATMTURB_LOWFREQ_NMODES;

    for (int m = 0; m < ATMTURB_LOWFREQ_NMODES; m++)
    {
        double g0, g1;
        atmturb_rng_gaussian_pair(&rng, &g0, &g1);
        float rp  = bl->lf_r[m];
        float b   = bl->lf_b[m];
        float wre = bl->lf_amp[m] * (float) g0;
        float wim = bl->lf_amp[m] * (float) g1;

        if (ar_order == 2)
        {
            float ure = rp * bl->lf_aux_re[m] + b * wre;
            float uim = rp * bl->lf_aux_im[m] + b * wim;
            bl->lf_aux_re[m] = ure;
            bl->lf_aux_im[m] = uim;
            bl->lf_cur_re[m] = rp * bl->lf_cur_re[m] + ure;
            bl->lf_cur_im[m] = rp * bl->lf_cur_im[m] + uim;
        }
        else
        {
            bl->lf_cur_re[m] = rp * bl->lf_cur_re[m] + b * wre;
            bl->lf_cur_im[m] = rp * bl->lf_cur_im[m] + b * wim;
        }

        bl->lf_ring_re[slot + m] = bl->lf_cur_re[m];
        bl->lf_ring_im[slot + m] = bl->lf_cur_im[m];
    }
    bl->max_lf_ready = frame_t;
}
