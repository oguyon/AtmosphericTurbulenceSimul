// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_rytov.c
 * @brief   Fourier-space Rytov diffractive propagation engine implementation
 */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "atmturb_geometry.h"
#include "atmturb_profile.h"
#include "atmturb_rolling.h"
#include "atmturb_rytov.h"
#include "atmturb_rytov_internal.h"
#include "atmturb_simd.h"
#include "atmturb_superlayer.h"
#include "atmturb_types.h"
#include "atmturb_wfs_render.h"

/**
 * atmturb_rytov_select_sec_mode - Determine whether Option B is valid for secondary wavelength
 * @plan: Rytov propagation plan.
 * @lambda_s_m: Secondary wavelength in meters.
 * @force_exact_sec: User flag forcing Option C.
 * @pup_size: Linear wavefront dimension in pixels.
 *
 * Return: 1 if Option B is valid, 0 if Option C fallback is required.
 */
static int atmturb_rytov_select_sec_mode(
    const atmturb_rytov_plan_t *plan,
    double                      lambda_s_m,
    int                         force_exact_sec,
    long                        pup_size)
{
    if (lambda_s_m <= 0.0 || force_exact_sec != 0)
    {
        return 0;
    }
    double max_offset_lim = (double) pup_size / 16.0;
    for (int m = 0; m < plan->nsuper; m++)
    {
        if (plan->supers[m].chrom_spread_px > 0.25 ||
            plan->supers[m].weight_ratio_spread > 0.02 ||
            fabs(plan->supers[m].chrom_dx_px) > max_offset_lim ||
            fabs(plan->supers[m].chrom_dy_px) > max_offset_lim)
        {
            return 0;
        }
    }
    return 1;
}

/**
 * atmturb_rytov_check_scintillation - Warn if turbulence is in strong scintillation regime
 * @plan: Initialized Rytov propagation plan.
 * @geom: Computed observing geometry.
 */
static void atmturb_rytov_check_scintillation(
    const atmturb_rytov_plan_t *plan,
    const atmturb_geom_t       *geom)
{
    if (geom->r0_ref_m <= 0.0 || plan->lambda_ref_m <= 0.0 || geom->cos_z <= 0.0)
    {
        return;
    }

    double lam = plan->lambda_ref_m;
    double int_cn2 = 0.060 * (lam * lam) * pow(geom->r0_ref_m, -5.0 / 3.0);
    double geom_factor = 19.12 * pow(lam, -7.0 / 6.0) * pow(geom->cos_z, -11.0 / 6.0);
    double sigma2 = 0.0;

    for (int k = 0; k < geom->nlayers; k++)
    {
        double dh = geom->layers[k].dist_m * geom->cos_z;
        if (dh > 0.0)
        {
            sigma2 += geom_factor * (geom->layers[k].weight * int_cn2) * pow(dh, 5.0 / 6.0);
        }
    }

    if (sigma2 > 0.30)
    {
        printf("[milkatmturb] Rytov warning: Strong scintillation regime "
               "(sigma_R^2 = %.2f > 0.30)\n", sigma2);
    }
}

/**
 * atmturb_rytov_print_diagnostics - Log Rytov physical validity diagnostics and warnings
 * @plan: Initialized Rytov propagation plan.
 * @geom: Computed observing geometry.
 */
static void atmturb_rytov_print_diagnostics(
    const atmturb_rytov_plan_t *plan,
    const atmturb_geom_t       *geom)
{
    atmturb_rytov_check_scintillation(plan, geom);

    if (plan->lambda_s_m > 0.0)
    {
        printf("[milkatmturb] Rytov: Secondary wavelength using Option %s\n",
               plan->sec_shared ? "B (shared spectrum + chromatic ramp)"
                                : "C (exact extrusion and FFT)");
    }

    double max_z = 0.0;
    for (int m = 0; m < plan->nsuper; m++)
    {
        if (plan->supers[m].dist_m > max_z)
        {
            max_z = plan->supers[m].dist_m;
        }
    }

    double dx = plan->pixscale_m;
    if (max_z > 0.0 && dx > 0.0)
    {
        double fresnel_scale = sqrt(plan->lambda_ref_m * max_z);
        if (fresnel_scale < 2.0 * dx)
        {
            printf("[milkatmturb] Rytov warning: Fresnel scale (%.3f m) < 2 dx (%.3f m)\n",
                   fresnel_scale, 2.0 * dx);
        }
        double chirp_limit = (double) plan->grid_size * dx / 4.0;
        if (plan->lambda_ref_m * max_z / dx > chirp_limit)
        {
            printf("[milkatmturb] Rytov warning: Potential chirp aliasing\n");
        }
    }

    if (geom->cos_z < cos(70.0 * M_PI / 180.0))
    {
        printf("[milkatmturb] Rytov warning: High zenith angle (>70 deg)\n");
    }
}

/**
 * atmturb_rytov_plan_init - Initialize Rytov diffractive filters and lookup tables
 * @plan: Propagation plan container to populate.
 * @prof: Active turbulence profile.
 * @geom: Computed observing geometry.
 * @pup_size: Linear dimension of wavefront in pixels.
 * @pixscale_m: Physical pixel scale [m/pixel].
 * @lambda_ref_m: Primary reference wavelength [m].
 * @lambda_s_m: Secondary observing wavelength [m].
 * @z_bin_m: Altitude binning distance threshold [m].
 * @force_exact_sec: Force Option C (separate secondary extrusion and FFT).
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_rytov_plan_init(
    atmturb_rytov_plan_t    *plan,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    long                     pup_size,
    double                   pixscale_m,
    double                   lambda_ref_m,
    double                   lambda_s_m,
    double                   z_bin_m,
    int                      force_exact_sec)
{
    if (plan == NULL || prof == NULL || geom == NULL || pup_size <= 0)
    {
        return -1;
    }

    memset(plan, 0, sizeof(*plan));
    plan->guard_pix    = (CONF_FRESNEL_GUARD_PIX > 0) ? (long) CONF_FRESNEL_GUARD_PIX : 0;
    plan->pad_size     = pup_size + 2 * plan->guard_pix;
    plan->grid_size    = plan->pad_size;
    plan->pixscale_m   = pixscale_m;
    plan->lambda_ref_m = lambda_ref_m;
    plan->lambda_s_m   = lambda_s_m;

    if (atmturb_superlayer_build(&plan->supers, &plan->nsuper, prof, geom,
                                 z_bin_m, pixscale_m) != 0)
    {
        return -1;
    }

    plan->sec_shared = atmturb_rytov_select_sec_mode(plan, lambda_s_m, force_exact_sec, pup_size);
    atmturb_rytov_print_diagnostics(plan, geom);

    if (atmturb_rytov_init_tables(plan) != 0)
    {
        atmturb_rytov_plan_free(plan);
        return -1;
    }

    int ns = plan->nsuper;
    plan->filter_a_pri = (float **) calloc((size_t) ns, sizeof(float *));
    plan->filter_b_pri = (float **) calloc((size_t) ns, sizeof(float *));
    plan->filter_a_sec = (float **) calloc((size_t) ns, sizeof(float *));
    plan->filter_b_sec = (float **) calloc((size_t) ns, sizeof(float *));

    if (plan->sec_shared)
    {
        plan->chrom_ramp = (fftwf_complex **) calloc((size_t) ns, sizeof(fftwf_complex *));
    }

    if (!plan->filter_a_pri || !plan->filter_b_pri ||
        !plan->filter_a_sec || !plan->filter_b_sec ||
        (plan->sec_shared && !plan->chrom_ramp))
    {
        atmturb_rytov_plan_free(plan);
        return -1;
    }

    size_t npix = (size_t) (plan->grid_size * (plan->grid_size / 2 + 1));
    for (int m = 0; m < ns; m++)
    {
        if (atmturb_rytov_build_layer_filters(plan, m, npix) != 0)
        {
            atmturb_rytov_plan_free(plan);
            return -1;
        }
    }

    return 0;
}

/**
 * atmturb_rytov_plan_free - Release resources held by Rytov propagation plan
 * @plan: Propagation plan container to tear down.
 */
void atmturb_rytov_plan_free(
    atmturb_rytov_plan_t *plan)
{
    if (plan == NULL)
    {
        return;
    }

    if (plan->supers != NULL)
    {
        atmturb_superlayer_free(plan->supers, plan->nsuper);
        plan->supers = NULL;
    }

    if (plan->filter_a_pri != NULL)
    {
        for (int m = 0; m < plan->nsuper; m++)
        {
            free(plan->filter_a_pri[m]);
        }
        free(plan->filter_a_pri);
        plan->filter_a_pri = NULL;
    }

    if (plan->filter_b_pri != NULL)
    {
        for (int m = 0; m < plan->nsuper; m++)
        {
            free(plan->filter_b_pri[m]);
        }
        free(plan->filter_b_pri);
        plan->filter_b_pri = NULL;
    }

    if (plan->filter_a_sec != NULL)
    {
        for (int m = 0; m < plan->nsuper; m++)
        {
            free(plan->filter_a_sec[m]);
        }
        free(plan->filter_a_sec);
        plan->filter_a_sec = NULL;
    }

    if (plan->filter_b_sec != NULL)
    {
        for (int m = 0; m < plan->nsuper; m++)
        {
            free(plan->filter_b_sec[m]);
        }
        free(plan->filter_b_sec);
        plan->filter_b_sec = NULL;
    }

    if (plan->chrom_ramp != NULL)
    {
        for (int m = 0; m < plan->nsuper; m++)
        {
            if (plan->chrom_ramp[m] != NULL)
            {
                fftwf_free(plan->chrom_ramp[m]);
            }
        }
        free(plan->chrom_ramp);
        plan->chrom_ramp = NULL;
    }

    free(plan->laplace_inv);
    plan->laplace_inv = NULL;
    if (plan->exp_y != NULL)
    {
        fftwf_free(plan->exp_y);
        plan->exp_y = NULL;
    }
    if (plan->exp_x != NULL)
    {
        fftwf_free(plan->exp_x);
        plan->exp_x = NULL;
    }
    plan->nsuper = 0;
}

/**
 * atmturb_rytov_ctx_init - Initialize thread-local Rytov render context and FFTW plans
 * @ctx: Thread context structure to initialize.
 * @pup_size: Linear dimension of wavefront in pixels.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_rytov_ctx_init(
    atmturb_rytov_ctx_t *ctx,
    long                 pup_size)
{
    if (ctx == NULL || pup_size <= 0)
    {
        return -1;
    }

    memset(ctx, 0, sizeof(*ctx));
    ctx->grid_size = pup_size;
    long n_half    = pup_size / 2 + 1;
    size_t npix_2d = (size_t) (pup_size * pup_size);
    size_t npix_spec = (size_t) (pup_size * n_half);

    ctx->real_in    = (float *) malloc(sizeof(float) * npix_2d);
    ctx->spec       = (fftwf_complex *) fftwf_alloc_complex(npix_spec);
    ctx->bound_a    = (float *) malloc(sizeof(float) * (size_t) pup_size);
    ctx->bound_b    = (float *) malloc(sizeof(float) * (size_t) pup_size);
    ctx->hat_a      = (fftwf_complex *) fftwf_alloc_complex((size_t) n_half);
    ctx->hat_b      = (fftwf_complex *) fftwf_alloc_complex((size_t) n_half);

    ctx->acc_dphi_pri = (fftwf_complex *) fftwf_alloc_complex(npix_spec);
    ctx->acc_chi_pri  = (fftwf_complex *) fftwf_alloc_complex(npix_spec);
    ctx->acc_dphi_sec = (fftwf_complex *) fftwf_alloc_complex(npix_spec);
    ctx->acc_chi_sec  = (fftwf_complex *) fftwf_alloc_complex(npix_spec);

    ctx->super_pha  = (float *) malloc(sizeof(float) * npix_2d);
    ctx->super_spha = (float *) malloc(sizeof(float) * npix_2d);
    ctx->dphi_out   = (float *) malloc(sizeof(float) * npix_2d);
    ctx->chi_out    = (float *) malloc(sizeof(float) * npix_2d);

    if (!ctx->real_in || !ctx->spec || !ctx->bound_a || !ctx->bound_b ||
        !ctx->hat_a || !ctx->hat_b || !ctx->acc_dphi_pri || !ctx->acc_chi_pri ||
        !ctx->acc_dphi_sec || !ctx->acc_chi_sec || !ctx->super_pha ||
        !ctx->super_spha || !ctx->dphi_out || !ctx->chi_out)
    {
        atmturb_rytov_ctx_free(ctx);
        return -1;
    }

    ctx->plan_r2c = fftwf_plan_dft_r2c_2d((int) pup_size, (int) pup_size,
                                          ctx->real_in, ctx->spec, FFTW_ESTIMATE);
    ctx->plan_c2r = fftwf_plan_dft_c2r_2d((int) pup_size, (int) pup_size,
                                          ctx->spec, ctx->real_in, FFTW_ESTIMATE);
    ctx->plan_1d_a = fftwf_plan_dft_r2c_1d((int) pup_size, ctx->bound_a,
                                           ctx->hat_a, FFTW_ESTIMATE);
    ctx->plan_1d_b = fftwf_plan_dft_r2c_1d((int) pup_size, ctx->bound_b,
                                           ctx->hat_b, FFTW_ESTIMATE);

    if (!ctx->plan_r2c || !ctx->plan_c2r || !ctx->plan_1d_a || !ctx->plan_1d_b)
    {
        atmturb_rytov_ctx_free(ctx);
        return -1;
    }
    return 0;
}

/**
 * atmturb_rytov_ctx_free - Release resources held by thread context
 * @ctx: Thread context structure to tear down.
 */
void atmturb_rytov_ctx_free(
    atmturb_rytov_ctx_t *ctx)
{
    if (ctx == NULL)
    {
        return;
    }

    if (ctx->plan_r2c)  fftwf_destroy_plan(ctx->plan_r2c);
    if (ctx->plan_c2r)  fftwf_destroy_plan(ctx->plan_c2r);
    if (ctx->plan_1d_a) fftwf_destroy_plan(ctx->plan_1d_a);
    if (ctx->plan_1d_b) fftwf_destroy_plan(ctx->plan_1d_b);

    free(ctx->real_in);
    if (ctx->spec)          fftwf_free(ctx->spec);
    free(ctx->bound_a);
    free(ctx->bound_b);
    if (ctx->hat_a)         fftwf_free(ctx->hat_a);
    if (ctx->hat_b)         fftwf_free(ctx->hat_b);
    if (ctx->acc_dphi_pri)  fftwf_free(ctx->acc_dphi_pri);
    if (ctx->acc_chi_pri)   fftwf_free(ctx->acc_chi_pri);
    if (ctx->acc_dphi_sec)  fftwf_free(ctx->acc_dphi_sec);
    if (ctx->acc_chi_sec)   fftwf_free(ctx->acc_chi_sec);

    free(ctx->super_pha);
    free(ctx->super_spha);
    free(ctx->dphi_out);
    free(ctx->chi_out);
}

/**
 * atmturb_rytov_accumulate_geom - Accumulate geometric phase from superlayer scratchpad
 * @pha: Destination primary phase slice.
 * @spha: Destination secondary phase slice.
 * @super_pha: Source padded superlayer primary phase.
 * @super_spha: Source padded superlayer secondary phase.
 * @pup_size: Linear dimension of pupil.
 * @guard: Guard band margin in pixels.
 * @pad_size: Linear dimension of padded compute grid.
 */
static void atmturb_rytov_accumulate_geom(
    float       *pha,
    float       *spha,
    const float *super_pha,
    const float *super_spha,
    long         pup_size,
    long         guard,
    long         pad_size)
{
    long pup_pixels = pup_size * pup_size;
    if (guard == 0)
    {
        atmturb_add_float_array(pha, super_pha, pup_pixels);
        atmturb_add_float_array(spha, super_spha, pup_pixels);
        return;
    }

    for (long y = 0; y < pup_size; y++)
    {
        long src_row = (y + guard) * pad_size + guard;
        long dst_row = y * pup_size;
        atmturb_add_float_array(&pha[dst_row], &super_pha[src_row], pup_size);
        atmturb_add_float_array(&spha[dst_row], &super_spha[src_row], pup_size);
    }
}

/**
 * atmturb_rytov_render_step - Render one frame using Fourier-space Rytov accumulation
 * @ctx: Thread-local execution context.
 * @plan: Precomputed Rytov propagation plan.
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @t: Simulation frame index.
 * @time_step_s: Time step between frames in seconds.
 * @master_size: Master screen dimension in pixels.
 * @pup_size: Linear dimension of wavefront in pixels.
 * @pha_slice: Destination primary phase slice.
 * @amp_slice: Destination primary amplitude slice.
 * @spha_slice: Destination secondary phase slice.
 * @samp_slice: Destination secondary amplitude slice.
 */
void atmturb_rytov_render_step(
    atmturb_rytov_ctx_t        *ctx,
    const atmturb_rytov_plan_t *plan,
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    long                        t,
    double                      time_step_s,
    long                        master_size,
    long                        pup_size,
    float                      *pha_slice,
    float                      *amp_slice,
    float                      *spha_slice,
    float                      *samp_slice)
{
    long pup_pixels = pup_size * pup_size;
    long guard = plan->guard_pix;
    long pad_size = plan->pad_size;
    long pad_pixels = pad_size * pad_size;
    long n_half = pad_size / 2 + 1;
    long n_spec = pad_size * n_half;

    memset(pha_slice, 0, sizeof(float) * (size_t) pup_pixels);
    memset(spha_slice, 0, sizeof(float) * (size_t) pup_pixels);
    memset(ctx->acc_dphi_pri, 0, sizeof(fftwf_complex) * (size_t) n_spec);
    memset(ctx->acc_chi_pri,  0, sizeof(fftwf_complex) * (size_t) n_spec);
    memset(ctx->acc_dphi_sec, 0, sizeof(fftwf_complex) * (size_t) n_spec);
    memset(ctx->acc_chi_sec,  0, sizeof(fftwf_complex) * (size_t) n_spec);

    atmturb_wfs_render_target_t target;
    target.pup_size     = pup_size;
    target.guard_pix    = guard;
    target.pha          = ctx->super_pha;
    target.spha         = ctx->super_spha;
    target.weight_scale = 1.0f;

    for (int m = 0; m < plan->nsuper; m++)
    {
        memset(ctx->super_pha, 0, sizeof(float) * (size_t) pad_pixels);
        memset(ctx->super_spha, 0, sizeof(float) * (size_t) pad_pixels);

        const atmturb_superlayer_t *sl = &plan->supers[m];
        for (int j = 0; j < sl->nlayers; j++)
        {
            int k = sl->layer_indices[j];
            target.weight_scale = (sl->layer_weights != NULL) ? sl->layer_weights[j] : 1.0f;
            atmturb_wfs_render_layer_target(r, geom, k, t, time_step_s, master_size, &target);
        }

        atmturb_rytov_accumulate_geom(pha_slice, spha_slice, ctx->super_pha,
                                      ctx->super_spha, pup_size, guard, pad_size);

        if (sl->dist_m > 0.0)
        {
            atmturb_rytov_decompose_periodic(ctx, plan, ctx->super_pha);
            atmturb_rytov_accumulate_filters(ctx->acc_dphi_pri, ctx->acc_chi_pri,
                                             ctx->spec, plan->filter_a_pri[m],
                                             plan->filter_b_pri[m], n_spec);

            if (plan->lambda_s_m > 0.0)
            {
                if (plan->sec_shared && plan->chrom_ramp != NULL && plan->chrom_ramp[m] != NULL)
                {
                    atmturb_rytov_accumulate_filters_rotated(ctx->acc_dphi_sec, ctx->acc_chi_sec,
                                                             ctx->spec, plan->chrom_ramp[m],
                                                             plan->filter_a_sec[m],
                                                             plan->filter_b_sec[m], n_spec);
                }
                else
                {
                    atmturb_rytov_decompose_periodic(ctx, plan, ctx->super_spha);
                    atmturb_rytov_accumulate_filters(ctx->acc_dphi_sec, ctx->acc_chi_sec,
                                                     ctx->spec, plan->filter_a_sec[m],
                                                     plan->filter_b_sec[m], n_spec);
                }
            }
        }
    }

    fftwf_execute_dft_c2r(ctx->plan_c2r, ctx->acc_dphi_pri, ctx->dphi_out);
    fftwf_execute_dft_c2r(ctx->plan_c2r, ctx->acc_chi_pri,  ctx->chi_out);
    atmturb_rytov_assemble_output(pup_size, guard, pad_size,
                                  pha_slice, amp_slice, ctx->dphi_out, ctx->chi_out);

    if (plan->lambda_s_m > 0.0)
    {
        fftwf_execute_dft_c2r(ctx->plan_c2r, ctx->acc_dphi_sec, ctx->dphi_out);
        fftwf_execute_dft_c2r(ctx->plan_c2r, ctx->acc_chi_sec,  ctx->chi_out);
        atmturb_rytov_assemble_output(pup_size, guard, pad_size,
                                      spha_slice, samp_slice, ctx->dphi_out, ctx->chi_out);
    }
}
