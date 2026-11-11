// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_rytov_kernels.c
 * @brief   Compute kernels and filter synthesis for Rytov propagation
 */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "atmturb_rytov.h"
#include "atmturb_rytov_internal.h"
#include "atmturb_rytov_simd.h"
#include "atmturb_simd.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

/**
 * atmturb_rytov_assemble_pixel - Evaluate one clamped boundary pixel for Rytov output
 */
static inline double atmturb_rytov_assemble_pixel(
    float       *pha,
    float       *amp,
    const float *dphi,
    const float *chi,
    long         dst_idx,
    long         src_idx)
{
    if (pha != NULL)
    {
        pha[dst_idx] += dphi[src_idx];
    }
    if (amp == NULL)
    {
        return 0.0;
    }
    float c = chi[src_idx];
    c = (c < -3.0f) ? -3.0f : ((c > 3.0f) ? 3.0f : c);
    float a      = expf(c);
    amp[dst_idx] = a;
    return (double) (a * a);
}

/**
 * atmturb_rytov_assemble_one_row - Assemble one pupil row with seam-peeled SIMD interior
 */
static inline double atmturb_rytov_assemble_one_row(
    long         y,
    long         pup_size,
    long         guard_pix,
    long         pad_size,
    float       *pha,
    float       *amp,
    const float *dphi,
    const float *chi)
{
    long max_clamp = (pup_size > 2) ? (pup_size - 2) : 0;
    long sy        = (guard_pix == 0 && (y == 0 || y == pup_size - 1))
                         ? ((y == 0) ? ((pup_size > 1) ? 1 : 0) : max_clamp)
                         : y;
    long src_row   = (sy + guard_pix) * pad_size + guard_pix;
    long dst_row   = y * pup_size;

    if (guard_pix > 0)
    {
        return atmturb_rytov_simd_assemble_row(
            (pha != NULL) ? &pha[dst_row] : NULL,
            (amp != NULL) ? &amp[dst_row] : NULL,
            &dphi[src_row], &chi[src_row], pup_size);
    }

    if (pup_size <= 2)
    {
        double s = 0.0;
        for (long x = 0; x < pup_size; x++)
        {
            long sx = (x == 0) ? ((pup_size > 1) ? 1 : 0) : max_clamp;
            s += atmturb_rytov_assemble_pixel(pha, amp, dphi, chi, dst_row + x, src_row + sx);
        }
        return s;
    }

    double sum_i = atmturb_rytov_assemble_pixel(pha, amp, dphi, chi, dst_row, src_row + 1);
    sum_i += atmturb_rytov_simd_assemble_row(
        (pha != NULL) ? &pha[dst_row + 1] : NULL,
        (amp != NULL) ? &amp[dst_row + 1] : NULL,
        &dphi[src_row + 1], &chi[src_row + 1], pup_size - 2);
    sum_i += atmturb_rytov_assemble_pixel(
        pha, amp, dphi, chi, dst_row + pup_size - 1, src_row + pup_size - 2);
    return sum_i;
}

/**
 * atmturb_rytov_assemble_output - Add diffractive phase and normalize amplitude
 * @pup_size: Linear dimension of pupil.
 * @guard_pix: Guard band margin in pixels.
 * @pad_size: Linear dimension of padded compute grid.
 * @pha: Accumulated geometric phase (updated in-place with diffractive correction).
 * @amp: Destination amplitude array.
 * @dphi: Reconstructed diffractive phase correction.
 * @chi: Reconstructed log-amplitude array.
 */
void atmturb_rytov_assemble_output(
    long         pup_size,
    long         guard_pix,
    long         pad_size,
    float       *pha,
    float       *amp,
    const float *dphi,
    const float *chi)
{
    long   npix  = pup_size * pup_size;
    double sum_i = 0.0;

#ifdef _OPENMP
    if (pup_size >= 512 && !omp_in_parallel())
    {
        int max_t = omp_get_max_threads();
        int nth   = (max_t > 8) ? 8 : max_t;
        #pragma omp parallel for num_threads(nth) reduction(+ : sum_i) schedule(static)
        for (long y = 0; y < pup_size; y++)
        {
            sum_i += atmturb_rytov_assemble_one_row(y, pup_size, guard_pix, pad_size,
                                                    pha, amp, dphi, chi);
        }
    }
    else
#endif
    {
        for (long y = 0; y < pup_size; y++)
        {
            sum_i += atmturb_rytov_assemble_one_row(y, pup_size, guard_pix, pad_size,
                                                    pha, amp, dphi, chi);
        }
    }

    if (amp != NULL)
    {
        float norm = (sum_i > 0.0) ? (float) (1.0 / sqrt(sum_i / (double) npix)) : 1.0f;
        atmturb_scale_float_array(amp, amp, norm, npix);
    }
}

/**
 * atmturb_rytov_decompose_rows - Apply Moisan harmonic subtraction over row range [j0, j1)
 */
static inline void atmturb_rytov_decompose_rows(
    atmturb_rytov_ctx_t        *ctx,
    const atmturb_rytov_plan_t *plan,
    long                        n,
    long                        n_half,
    long                        j0,
    long                        j1)
{
    for (long j = j0; j < j1; j++)
    {
        float b_re;
        float b_im;
        if (j <= n / 2)
        {
            b_re = ctx->hat_b[j][0];
            b_im = ctx->hat_b[j][1];
        }
        else
        {
            long j_sym = n - j;
            b_re       =  ctx->hat_b[j_sym][0];
            b_im       = -ctx->hat_b[j_sym][1];
        }

        float ey_re   = plan->exp_y[j][0];
        float ey_im   = plan->exp_y[j][1];
        long  row_idx = j * n_half;

        atmturb_rytov_simd_decompose_row(&ctx->spec[row_idx], ctx->hat_a,
                                         plan->exp_x, &plan->laplace_inv[row_idx],
                                         ey_re, ey_im, b_re, b_im, n_half);
    }
}

/**
 * atmturb_rytov_decompose_periodic - Subtract harmonic smooth component from 2D spectrum
 * @ctx: Thread execution context.
 * @plan: Propagation plan.
 * @phi: Input 2D real phase screen (32/64-byte FFTW-aligned).
 */
void atmturb_rytov_decompose_periodic(
    atmturb_rytov_ctx_t        *ctx,
    const atmturb_rytov_plan_t *plan,
    float                      *phi)
{
    long                  n      = plan->grid_size;
    long                  n_half = n / 2 + 1;
    float       *restrict ba     = ctx->bound_a;
    float       *restrict bb     = ctx->bound_b;
    const float *restrict r_top  = phi;
    const float *restrict r_bot  = &phi[(n - 1) * n];

    #pragma omp simd
    for (long x = 0; x < n; x++)
    {
        ba[x] = r_top[x] - r_bot[x];
    }
    for (long y = 0; y < n; y++)
    {
        bb[y] = phi[y * n] - phi[y * n + (n - 1)];
    }

    fftwf_execute(ctx->plan_1d_a);
    fftwf_execute(ctx->plan_1d_b);
    fftwf_execute_dft_r2c(ctx->plan_r2c, phi, ctx->spec);

#ifdef _OPENMP
    if (n >= 512 && !omp_in_parallel())
    {
        int max_t = omp_get_max_threads();
        int nth   = (max_t > 8) ? 8 : max_t;
        #pragma omp parallel num_threads(nth)
        {
            int  tid = omp_get_thread_num();
            int  nt  = omp_get_num_threads();
            long j0  = (n * (long) tid) / (long) nt;
            long j1  = (n * (long) (tid + 1)) / (long) nt;
            atmturb_rytov_decompose_rows(ctx, plan, n, n_half, j0, j1);
        }
        return;
    }
#endif
    atmturb_rytov_decompose_rows(ctx, plan, n, n_half, 0, n);
}

/**
 * atmturb_rytov_accumulate_filters - Multiply spectrum by filters and accumulate
 * @acc_dphi: Phase correction accumulator.
 * @acc_chi: Log-amplitude accumulator.
 * @spec: Periodic spectrum input.
 * @filt_a: Phase filter array.
 * @filt_b: Amplitude filter array.
 * @ntot: Elements in half-spectrum.
 */
void atmturb_rytov_accumulate_filters(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          ntot)
{
    atmturb_rytov_simd_accumulate_filters(acc_dphi, acc_chi, spec, filt_a, filt_b, ntot);
}

/**
 * atmturb_rytov_accumulate_filters_rotated - Apply chromatic ramp and accumulate filters
 * @acc_dphi: Phase correction accumulator.
 * @acc_chi: Log-amplitude accumulator.
 * @spec: Periodic spectrum input.
 * @ramp: Complex chromatic phase factor array.
 * @filt_a: Phase filter array.
 * @filt_b: Amplitude filter array.
 * @ntot: Elements in half-spectrum.
 */
void atmturb_rytov_accumulate_filters_rotated(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const fftwf_complex *restrict ramp,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          ntot)
{
    atmturb_rytov_simd_accumulate_filters_rotated(acc_dphi, acc_chi, spec, ramp,
                                                  filt_a, filt_b, ntot);
}

/**
 * atmturb_rytov_init_tables - Build Poisson solver lookup tables
 * @plan: Propagation plan to populate.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_rytov_init_tables(
    atmturb_rytov_plan_t *plan)
{
    long n      = plan->grid_size;
    long n_half = n / 2 + 1;
    size_t npix = (size_t) (n * n_half);

    plan->laplace_inv = (float *) malloc(sizeof(float) * npix);
    plan->exp_y       = (fftwf_complex *) fftwf_alloc_complex((size_t) n);
    plan->exp_x       = (fftwf_complex *) fftwf_alloc_complex((size_t) n_half);

    if (plan->laplace_inv == NULL || plan->exp_y == NULL || plan->exp_x == NULL)
    {
        return -1;
    }

    for (long j = 0; j < n; j++)
    {
        double theta = 2.0 * M_PI * (double) j / (double) n;
        plan->exp_y[j][0] = (float) (1.0 - cos(theta));
        plan->exp_y[j][1] = (float) (-sin(theta));
    }

    for (long i = 0; i < n_half; i++)
    {
        double theta = 2.0 * M_PI * (double) i / (double) n;
        plan->exp_x[i][0] = (float) (1.0 - cos(theta));
        plan->exp_x[i][1] = (float) (-sin(theta));
    }

    for (long j = 0; j < n; j++)
    {
        double theta_y = 2.0 * M_PI * (double) j / (double) n;
        long row_idx = j * n_half;
        for (long i = 0; i < n_half; i++)
        {
            double theta_x = 2.0 * M_PI * (double) i / (double) n;
            double denom = 2.0 * cos(theta_x) + 2.0 * cos(theta_y) - 4.0;
            plan->laplace_inv[row_idx + i] = (j == 0 && i == 0)
                                            ? 0.0f
                                            : (float) (1.0 / denom);
        }
    }
    return 0;
}

/**
 * atmturb_rytov_build_layer_filters - Compute transfer function arrays for one super-layer
 * @plan: Propagation plan container.
 * @m: Super-layer index.
 * @npix: Total pixels in half-spectrum.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_rytov_build_layer_filters(
    atmturb_rytov_plan_t *plan,
    int                   m,
    size_t                npix)
{
    plan->filter_a_pri[m] = (float *) malloc(sizeof(float) * npix);
    plan->filter_b_pri[m] = (float *) malloc(sizeof(float) * npix);
    if (plan->filter_a_pri[m] == NULL || plan->filter_b_pri[m] == NULL)
    {
        return -1;
    }

    if (plan->lambda_s_m > 0.0)
    {
        plan->filter_a_sec[m] = (float *) malloc(sizeof(float) * npix);
        plan->filter_b_sec[m] = (float *) malloc(sizeof(float) * npix);
        if (plan->filter_a_sec[m] == NULL || plan->filter_b_sec[m] == NULL)
        {
            return -1;
        }

        if (plan->sec_shared && plan->chrom_ramp != NULL)
        {
            plan->chrom_ramp[m] = (fftwf_complex *) fftwf_alloc_complex(npix);
            if (plan->chrom_ramp[m] == NULL)
            {
                return -1;
            }
        }
    }

    long n = plan->grid_size;
    long n_half = n / 2 + 1;
    double l_grid = (double) n * plan->pixscale_m;
    double inv_n2 = 1.0 / (double) (n * n);
    double dist_s = (plan->supers[m].dist_s_m > 0.0) ? plan->supers[m].dist_s_m
                                                     : plan->supers[m].dist_m;
    double coeff_pri = M_PI * plan->supers[m].dist_m * plan->lambda_ref_m / (l_grid * l_grid);
    double coeff_sec = M_PI * dist_s * plan->lambda_s_m / (l_grid * l_grid);

    double dx = plan->supers[m].chrom_dx_px;
    double dy = plan->supers[m].chrom_dy_px;
    double wr = plan->supers[m].weight_ratio;

    for (long j = 0; j < n; j++)
    {
        long fy = (j < n / 2) ? j : (j - n);
        double fy2 = (double) (fy * fy);
        long row = j * n_half;

        for (long i = 0; i < n_half; i++)
        {
            long fx = i;
            double sqdist = (double) (fx * fx) + fy2;
            double alpha_pri = coeff_pri * sqdist;
            double a_pri = (cos(alpha_pri) - 1.0) * inv_n2;
            double b_pri = sin(alpha_pri) * inv_n2;

            plan->filter_a_pri[m][row + i] = (float) a_pri;
            plan->filter_b_pri[m][row + i] = (float) b_pri;

            if (plan->lambda_s_m > 0.0)
            {
                double alpha_sec = coeff_sec * sqdist;
                double a_s = (cos(alpha_sec) - 1.0) * inv_n2;
                double b_s = sin(alpha_sec) * inv_n2;

                plan->filter_a_sec[m][row + i] = (float) a_s;
                plan->filter_b_sec[m][row + i] = (float) b_s;

                if (plan->sec_shared && plan->chrom_ramp != NULL && plan->chrom_ramp[m] != NULL)
                {
                    double beta = 2.0 * M_PI * ((double) fx * dx + (double) fy * dy) / (double) n;
                    plan->chrom_ramp[m][row + i][0] = (float) (wr * cos(beta));
                    plan->chrom_ramp[m][row + i][1] = (float) (wr * sin(beta));
                }
            }
        }
    }
    return 0;
}

/**
 * atmturb_rytov_get_wisdom_path - Determine persistent FFTW wisdom file path
 * @path: Output string buffer.
 * @maxlen: Size of output buffer in bytes.
 *
 * Return: 0 on success.
 */
static int atmturb_rytov_get_wisdom_path(
    char  *path,
    size_t maxlen)
{
    const char *env_path = getenv("ATMTURB_FFTW_WISDOM");
    if (env_path != NULL && env_path[0] != '\0')
    {
        snprintf(path, maxlen, "%s", env_path);
        return 0;
    }

    const char *home = getenv("HOME");
    if (home != NULL && home[0] != '\0')
    {
        char milk_dir[512];
        snprintf(milk_dir, sizeof(milk_dir), "%s/.milk", home);
        mkdir(milk_dir, 0755);
        snprintf(path, maxlen, "%s/.milk/atmturb_fftw_wisdom.txt", home);
        return 0;
    }

    snprintf(path, maxlen, ".atmturb_fftw_wisdom.txt");
    return 0;
}

/**
 * atmturb_rytov_wisdom_load - Import pre-computed FFTW wisdom from disk
 */
void atmturb_rytov_wisdom_load(void)
{
    static int s_loaded = 0;
    if (s_loaded != 0)
    {
        return;
    }
    s_loaded = 1;

    char path[512];
    if (atmturb_rytov_get_wisdom_path(path, sizeof(path)) == 0)
    {
        fftwf_import_wisdom_from_filename(path);
    }
}

/**
 * atmturb_rytov_wisdom_save - Export updated FFTW wisdom to disk
 */
void atmturb_rytov_wisdom_save(void)
{
    char path[512];
    if (atmturb_rytov_get_wisdom_path(path, sizeof(path)) == 0)
    {
        fftwf_export_wisdom_to_filename(path);
    }
}
