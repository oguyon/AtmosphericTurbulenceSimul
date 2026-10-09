// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_rytov_kernels.c
 * @brief   Compute kernels and filter synthesis for Rytov propagation
 */

#define _GNU_SOURCE
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "atmturb_rytov.h"
#include "atmturb_rytov_internal.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

/**
 * atmturb_rytov_assemble_output - Add diffractive phase and normalize amplitude
 * @pup_size: Linear dimension of pupil.
 * @pha: Accumulated geometric phase (updated in-place with diffractive correction).
 * @amp: Destination amplitude array.
 * @dphi: Reconstructed diffractive phase correction.
 * @chi: Reconstructed log-amplitude array.
 */
void atmturb_rytov_assemble_output(
    long         pup_size,
    float       *pha,
    float       *amp,
    const float *dphi,
    const float *chi)
{
    long npix = pup_size * pup_size;
    double sum_i = 0.0;

    for (long i = 0; i < npix; i++)
    {
        pha[i] += dphi[i];
        float a = expf(chi[i]);
        amp[i]  = a;
        sum_i  += (double) (a * a);
    }

    float norm = (sum_i > 0.0) ? (float) (1.0 / sqrt(sum_i / (double) npix)) : 1.0f;
    for (long i = 0; i < npix; i++)
    {
        amp[i] *= norm;
    }
}

/**
 * atmturb_rytov_decompose_periodic - Subtract harmonic smooth component from 2D spectrum
 * @ctx: Thread execution context.
 * @plan: Propagation plan.
 * @phi: Input 2D real phase screen.
 */
void atmturb_rytov_decompose_periodic(
    atmturb_rytov_ctx_t        *ctx,
    const atmturb_rytov_plan_t *plan,
    const float                *phi)
{
    long n    = plan->grid_size;
    long n_half = n / 2 + 1;

    for (long x = 0; x < n; x++)
    {
        ctx->bound_a[x] = phi[0 * n + x] - phi[(n - 1) * n + x];
    }
    for (long y = 0; y < n; y++)
    {
        ctx->bound_b[y] = phi[y * n + 0] - phi[y * n + (n - 1)];
    }

    memcpy(ctx->real_in, phi, sizeof(float) * (size_t) (n * n));
    fftwf_execute(ctx->plan_1d_a);
    fftwf_execute(ctx->plan_1d_b);
    fftwf_execute(ctx->plan_r2c);

    for (long j = 0; j < n; j++)
    {
        float b_re, b_im;
        if (j <= n / 2)
        {
            b_re = ctx->hat_b[j][0];
            b_im = ctx->hat_b[j][1];
        }
        else
        {
            long j_sym = n - j;
            b_re =  ctx->hat_b[j_sym][0];
            b_im = -ctx->hat_b[j_sym][1];
        }

        float ey_re = plan->exp_y[j][0];
        float ey_im = plan->exp_y[j][1];
        long row_idx = j * n_half;

        for (long i = 0; i < n_half; i++)
        {
            float a_re = ctx->hat_a[i][0];
            float a_im = ctx->hat_a[i][1];
            float ex_re = plan->exp_x[i][0];
            float ex_im = plan->exp_x[i][1];

            float v_re = (a_re * ey_re - a_im * ey_im) + (b_re * ex_re - b_im * ex_im);
            float v_im = (a_re * ey_im + a_im * ey_re) + (b_re * ex_im + b_im * ex_re);

            float linv = plan->laplace_inv[row_idx + i];
            ctx->spec[row_idx + i][0] += v_re * linv;
            ctx->spec[row_idx + i][1] += v_im * linv;
        }
    }
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
    fftwf_complex       *acc_dphi,
    fftwf_complex       *acc_chi,
    const fftwf_complex *spec,
    const float         *filt_a,
    const float         *filt_b,
    long                 ntot)
{
    for (long k = 0; k < ntot; k++)
    {
        float p_re = spec[k][0];
        float p_im = spec[k][1];
        float fa   = filt_a[k];
        float fb   = filt_b[k];

        acc_dphi[k][0] += p_re * fa;
        acc_dphi[k][1] += p_im * fa;
        acc_chi[k][0]  += p_re * fb;
        acc_chi[k][1]  += p_im * fb;
    }
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
    fftwf_complex       *acc_dphi,
    fftwf_complex       *acc_chi,
    const fftwf_complex *spec,
    const fftwf_complex *ramp,
    const float         *filt_a,
    const float         *filt_b,
    long                 ntot)
{
    for (long k = 0; k < ntot; k++)
    {
        float p_re = spec[k][0];
        float p_im = spec[k][1];
        float cb   = ramp[k][0];
        float sb   = ramp[k][1];

        float ps_re = p_re * cb - p_im * sb;
        float ps_im = p_re * sb + p_im * cb;

        float fa = filt_a[k];
        float fb = filt_b[k];

        acc_dphi[k][0] += ps_re * fa;
        acc_dphi[k][1] += ps_im * fa;
        acc_chi[k][0]  += ps_re * fb;
        acc_chi[k][1]  += ps_im * fb;
    }
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
    double coeff_pri = M_PI * plan->supers[m].dist_m * plan->lambda_ref_m / (l_grid * l_grid);
    double coeff_sec = M_PI * plan->supers[m].dist_m * plan->lambda_s_m / (l_grid * l_grid);

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
