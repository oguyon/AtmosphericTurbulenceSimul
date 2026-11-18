// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_validate_simd.c
 * @brief   Scalar vs dispatched SIMD parity checks across extrusion, Fresnel, and Rytov kernels
 */

#define _GNU_SOURCE
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "atmturb_fresnel_simd.h"
#include "atmturb_lowfreq.h"
#include "atmturb_rytov_simd.h"
#include "atmturb_simd.h"
#include "atmturb_validate_math.h"

/**
 * val_simd_compare_buffers - Record maximum absolute element-wise difference between float buffers
 */
static void val_simd_compare_buffers(
    const float *out_s,
    const float *out_v,
    long         count,
    double      *max_diff)
{
    for (long i = 0; i < count; i++)
    {
        double diff = fabs((double) out_s[i] - (double) out_v[i]);
        if (diff > *max_diff)
        {
            *max_diff = diff;
        }
    }
}

/**
 * val_simd_compare_cplx - Record maximum absolute element-wise difference between complex buffers
 */
static void val_simd_compare_cplx(
    const fftwf_complex *out_s,
    const fftwf_complex *out_v,
    long                 count,
    double              *max_diff)
{
    val_simd_compare_buffers((const float *) out_s, (const float *) out_v, 2 * count, max_diff);
}

/**
 * val_check_extrude_cases - Test extrusion SIMD parity across schemes, oversampling, and offsets
 */
static void val_check_extrude_cases(
    const float *master,
    long         msize,
    long         pup_size,
    float       *out_s,
    float       *out_v,
    double      *max_diff)
{
    long ptot = pup_size * pup_size;
    const double test_x[] = {
        15.3,
        (double) (msize - pup_size) - 0.5,
        (double) (msize - pup_size),
        (double) (msize - pup_size / 2) + 0.5,
        -15.3
    };
    const double test_y[] = {
        22.7,
        (double) (msize - pup_size) - 0.5,
        (double) (msize - pup_size),
        (double) (msize - pup_size / 2) + 0.5,
        -22.7
    };
    int ncases = (int) (sizeof(test_x) / sizeof(test_x[0]));

    for (int interp = 0; interp <= 1; interp++)
    {
        for (int os = 1; os <= 3; os++)
        {
            if (pup_size * os >= msize)
            {
                continue;
            }
            for (int c = 0; c < ncases; c++)
            {
                atmturb_extrude_params_t ep_s = {
                    .master   = master,
                    .msize    = msize,
                    .x0       = test_x[c],
                    .y0       = test_y[c],
                    .pup_size = pup_size,
                    .os       = os,
                    .interp   = interp,
                    .weight   = 1.25f,
                    .out_pha  = out_s
                };
                atmturb_extrude_params_t ep_v = ep_s;
                ep_v.out_pha = out_v;

                memset(out_s, 0, sizeof(float) * (size_t) ptot);
                memset(out_v, 0, sizeof(float) * (size_t) ptot);

                atmturb_extrude_accumulate_scalar(&ep_s);
                atmturb_extrude_accumulate(&ep_v);
                val_simd_compare_buffers(out_s, out_v, ptot, max_diff);
            }
        }
    }
}

/**
 * val_check_lowfreq_cases - Test low-frequency subharmonic SIMD parity across strides and wraps
 */
static void val_check_lowfreq_cases(
    long    msize,
    long    pup_size,
    float  *out_s,
    float  *out_v,
    double *max_diff)
{
    long ptot = pup_size * pup_size;
    atmturb_lowfreq_t lf;
    atmturb_lowfreq_init(&lf, msize, 25.0, 100.0, 0.1, 42);

    const double test_x[] = {12.34, (double) (msize - pup_size) + 0.5, -12.34};
    const double test_y[] = {56.78, (double) (msize - pup_size) + 0.5, -56.78};
    int ncases = (int) (sizeof(test_x) / sizeof(test_x[0]));

    for (int sidx = 0; sidx <= 1; sidx++)
    {
        for (int os = 1; os <= 3; os++)
        {
            for (int c = 0; c < ncases; c++)
            {
                atmturb_lowfreq_params_t lp_s = {
                    .lf         = &lf,
                    .screen_idx = sidx,
                    .x0         = test_x[c],
                    .y0         = test_y[c],
                    .pup_size   = pup_size,
                    .os         = os,
                    .weight     = 1.15f,
                    .out_pha    = out_s
                };
                atmturb_lowfreq_params_t lp_v = lp_s;
                lp_v.out_pha = out_v;

                memset(out_s, 0, sizeof(float) * (size_t) ptot);
                memset(out_v, 0, sizeof(float) * (size_t) ptot);

                atmturb_extrude_lowfreq_scalar(&lp_s);
                atmturb_extrude_lowfreq(&lp_v);
                val_simd_compare_buffers(out_s, out_v, ptot, max_diff);
            }
        }
    }
}

/**
 * val_check_fresnel_simd_cases - Test Fresnel modulation, extraction, phasor, and spectrum SIMD
 */
static int val_check_fresnel_simd_cases(
    long    npix,
    double *max_diff)
{
    fftwf_complex *fld_s = fftwf_alloc_complex((size_t) npix);
    fftwf_complex *fld_v = fftwf_alloc_complex((size_t) npix);
    float *pha_in = (float *) calloc((size_t) npix, sizeof(float));
    float *mask   = (float *) calloc((size_t) npix, sizeof(float));
    float *pha_s  = (float *) calloc((size_t) npix, sizeof(float));
    float *pha_v  = (float *) calloc((size_t) npix, sizeof(float));
    float *amp_s  = (float *) calloc((size_t) npix, sizeof(float));
    float *amp_v  = (float *) calloc((size_t) npix, sizeof(float));

    int ok = (fld_s && fld_v && pha_in && mask && pha_s && pha_v && amp_s && amp_v);
    if (!ok)
    {
        goto cleanup;
    }

    for (long i = 0; i < npix; i++)
    {
        pha_in[i] = 3.5f * sinf((float) i * 0.13f) - 1.7f * cosf((float) i * 0.07f);
        mask[i]   = 0.5f + 0.5f * cosf((float) i * 0.05f);
        pha_s[i]  = pha_in[i];
        pha_v[i]  = pha_in[i];
    }

    atmturb_fresnel_init_modulate_field_scalar(fld_s, pha_in, npix);
    atmturb_fresnel_simd_init_modulate_field(fld_v, pha_in, npix);
    val_simd_compare_cplx(fld_s, fld_v, npix, max_diff);

    atmturb_fresnel_modulate_field_scalar(fld_s, pha_in, npix);
    atmturb_fresnel_simd_modulate_field(fld_v, pha_in, npix);
    val_simd_compare_cplx(fld_s, fld_v, npix, max_diff);

    atmturb_fresnel_extract_unwrapped_scalar(fld_s, pha_s, amp_s, npix);
    atmturb_fresnel_simd_extract_unwrapped(fld_v, pha_v, amp_v, npix);
    val_simd_compare_buffers(pha_s, pha_v, npix, max_diff);
    val_simd_compare_buffers(amp_s, amp_v, npix, max_diff);

    atmturb_fresnel_phasor_masked_scalar(fld_s, pha_in, mask, 1.35f, npix);
    atmturb_fresnel_simd_phasor_masked(fld_v, pha_in, mask, 1.35f, npix);
    val_simd_compare_cplx(fld_s, fld_v, npix, max_diff);

    float u1[8], u2[8];
    const float spec_p[4] = {4.0f, 0.01f, 100.0f, 0.001f};
    for (int k = 0; k < 8; k++)
    {
        u1[k] = 0.08f + 0.11f * (float) k;
        u2[k] = 0.05f + 0.12f * (float) k;
    }
    atmturb_fresnel_spectrum_strip_scalar(fld_s, u1, u2, 16.0f, spec_p);
    atmturb_fresnel_simd_spectrum_strip(fld_v, u1, u2, 16.0f, spec_p);
    val_simd_compare_cplx(fld_s, fld_v, 8, max_diff);

cleanup:
    fftwf_free(fld_s);
    fftwf_free(fld_v);
    free(pha_in);
    free(mask);
    free(pha_s);
    free(pha_v);
    free(amp_s);
    free(amp_v);
    return ok;
}

/**
 * val_check_rytov_simd_cases - Test Rytov filter accumulation, Moisan row, and assembly SIMD
 */
static int val_check_rytov_simd_cases(
    long    n,
    double *max_diff)
{
    fftwf_complex *dp_s  = fftwf_alloc_complex((size_t) n);
    fftwf_complex *dp_v  = fftwf_alloc_complex((size_t) n);
    fftwf_complex *ch_s  = fftwf_alloc_complex((size_t) n);
    fftwf_complex *ch_v  = fftwf_alloc_complex((size_t) n);
    fftwf_complex *spec  = fftwf_alloc_complex((size_t) n);
    fftwf_complex *ramp  = fftwf_alloc_complex((size_t) n);
    float *fa = (float *) calloc((size_t) n, sizeof(float));
    float *fb = (float *) calloc((size_t) n, sizeof(float));

    int ok = (dp_s && dp_v && ch_s && ch_v && spec && ramp && fa && fb);
    if (!ok)
    {
        goto cleanup;
    }

    memset(dp_s, 0, sizeof(fftwf_complex) * (size_t) n);
    memset(dp_v, 0, sizeof(fftwf_complex) * (size_t) n);
    memset(ch_s, 0, sizeof(fftwf_complex) * (size_t) n);
    memset(ch_v, 0, sizeof(fftwf_complex) * (size_t) n);

    for (long i = 0; i < n; i++)
    {
        spec[i][0] = 0.7f * cosf((float) i * 0.11f);
        spec[i][1] = 0.6f * sinf((float) i * 0.09f);
        ramp[i][0] = cosf((float) i * 0.05f);
        ramp[i][1] = sinf((float) i * 0.05f);
        fa[i]      = -0.25f + 0.4f * sinf((float) i * 0.03f);
        fb[i]      = 0.35f * cosf((float) i * 0.04f);
    }

    atmturb_rytov_accumulate_filters_scalar(dp_s, ch_s, spec, fa, fb, n);
    atmturb_rytov_simd_accumulate_filters(dp_v, ch_v, spec, fa, fb, n);
    atmturb_rytov_accumulate_filters_rotated_scalar(dp_s, ch_s, spec, ramp, fa, fb, n);
    atmturb_rytov_simd_accumulate_filters_rotated(dp_v, ch_v, spec, ramp, fa, fb, n);
    val_simd_compare_cplx(dp_s, dp_v, n, max_diff);
    val_simd_compare_cplx(ch_s, ch_v, n, max_diff);

    memcpy(dp_s, spec, sizeof(fftwf_complex) * (size_t) n);
    memcpy(dp_v, spec, sizeof(fftwf_complex) * (size_t) n);
    atmturb_rytov_decompose_row_scalar(dp_s, ch_s, ramp, fb, 0.8f, -0.6f, 0.3f, 0.4f, n);
    atmturb_rytov_simd_decompose_row(dp_v, ch_v, ramp, fb, 0.8f, -0.6f, 0.3f, 0.4f, n);
    val_simd_compare_cplx(dp_s, dp_v, n, max_diff);

    float *ps = (float *) dp_s, *pv = (float *) dp_v;
    float *as = (float *) ch_s, *av = (float *) ch_v;
    for (long i = 0; i < n; i++)
    {
        ps[i] = 0.5f * sinf((float) i * 0.07f);
        pv[i] = ps[i];
    }
    double ss = atmturb_rytov_assemble_row_scalar(ps, as, fa, fb, n);
    double sv = atmturb_rytov_simd_assemble_row(pv, av, fa, fb, n);
    val_simd_compare_buffers(ps, pv, n, max_diff);
    val_simd_compare_buffers(as, av, n, max_diff);
    double sdiff = fabs(ss - sv) / (double) n;
    if (sdiff > *max_diff)
    {
        *max_diff = sdiff;
    }

cleanup:
    fftwf_free(dp_s);
    fftwf_free(dp_v);
    fftwf_free(ch_s);
    fftwf_free(ch_v);
    fftwf_free(spec);
    fftwf_free(ramp);
    free(fa);
    free(fb);
    return ok;
}

/**
 * val_check_simd_parity - Compare scalar and vectorized SIMD kernels across modules
 */
int val_check_simd_parity(
    double  tol,
    char   *msg_buf,
    size_t  msg_size)
{
    long msize    = 4096;
    long max_pup  = 1100;
    long mtot     = msize * msize;
    long max_ptot = max_pup * max_pup;

    float *master = (float *) malloc(sizeof(float) * (size_t) mtot);
    float *out_s  = (float *) calloc((size_t) max_ptot, sizeof(float));
    float *out_v  = (float *) calloc((size_t) max_ptot, sizeof(float));

    if (!master || !out_s || !out_v)
    {
        free(master);
        free(out_s);
        free(out_v);
        snprintf(msg_buf, msg_size, "allocation failure in SIMD parity test");
        return 0;
    }

    for (long i = 0; i < mtot; i++)
    {
        master[i] = sinf((float) i * 0.1f) + cosf((float) i * 0.03f);
    }

    double max_diff = 0.0;
    long pup_sizes[] = {64, 1100};
    int npup = (int) (sizeof(pup_sizes) / sizeof(pup_sizes[0]));

    for (int p = 0; p < npup; p++)
    {
        val_check_extrude_cases(master, msize, pup_sizes[p], out_s, out_v, &max_diff);
        val_check_lowfreq_cases(msize, pup_sizes[p], out_s, out_v, &max_diff);
    }

    int ok_prop = val_check_fresnel_simd_cases(263, &max_diff)
                  && val_check_rytov_simd_cases(263, &max_diff);

    free(master);
    free(out_s);
    free(out_v);

    if (!ok_prop)
    {
        snprintf(msg_buf, msg_size, "allocation failure in propagation SIMD parity test");
        return 0;
    }

    snprintf(msg_buf, msg_size, "SIMD parity max diff = %.3e (tol %.3e)", max_diff, tol);
    return (max_diff <= tol);
}
