// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_fresnel.c
 * @brief   Super-layer binning and diffractive propagation engine implementation
 */

#define _GNU_SOURCE
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "atmturb_fresnel.h"
#include "atmturb_simd.h"
#include "atmturb_superlayer.h"

/**
 * atmturb_fresnel_plan_init - Initialize super-layer binning and transfer functions
 * @plan: Propagation plan container to populate.
 * @prof: Active turbulence profile.
 * @geom: Computed observing geometry.
 * @pup_size: Linear dimension of wavefront in pixels.
 * @pixscale_m: Physical pixel scale [m/pixel].
 * @lambda_ref_m: Primary reference wavelength [m].
 * @lambda_s_m: Secondary observing wavelength [m].
 * @z_bin_m: Altitude binning distance threshold [m].
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_fresnel_plan_init(
    atmturb_fresnel_plan_t  *plan,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    long                     pup_size,
    double                   pixscale_m,
    double                   lambda_ref_m,
    double                   lambda_s_m,
    double                   z_bin_m)
{
    if (plan == NULL || prof == NULL || geom == NULL || pup_size <= 0)
    {
        return -1;
    }

    memset(plan, 0, sizeof(*plan));
    plan->grid_size    = pup_size;
    plan->pixscale_m   = pixscale_m;
    plan->lambda_ref_m = lambda_ref_m;
    plan->lambda_s_m   = lambda_s_m;

    if (atmturb_superlayer_build(&plan->supers, &plan->nsuper, prof, geom,
                                 z_bin_m, pixscale_m) != 0)
    {
        return -1;
    }
    int nsuper = plan->nsuper;

    size_t npix = (size_t) (pup_size * pup_size);
    plan->tf_pri = (fftwf_complex **) calloc((size_t) nsuper, sizeof(fftwf_complex *));
    plan->tf_sec = (fftwf_complex **) calloc((size_t) nsuper, sizeof(fftwf_complex *));
    if (plan->tf_pri == NULL || plan->tf_sec == NULL)
    {
        atmturb_fresnel_plan_free(plan);
        return -1;
    }

    for (int m = 0; m < nsuper; m++)
    {
        plan->tf_pri[m] = (fftwf_complex *) fftwf_alloc_complex(npix);
        if (plan->tf_pri[m] == NULL)
        {
            atmturb_fresnel_plan_free(plan);
            return -1;
        }
        wfprop_fresnel_tf_build(plan->tf_pri[m], pup_size, pixscale_m,
                                plan->supers[m].step_dist_m, lambda_ref_m, 0.0);

        if (lambda_s_m > 0.0)
        {
            plan->tf_sec[m] = (fftwf_complex *) fftwf_alloc_complex(npix);
            if (plan->tf_sec[m] == NULL)
            {
                atmturb_fresnel_plan_free(plan);
                return -1;
            }
            wfprop_fresnel_tf_build(plan->tf_sec[m], pup_size, pixscale_m,
                                    plan->supers[m].step_dist_m, lambda_s_m, 0.0);
        }
    }

    return 0;
}

/**
 * atmturb_fresnel_plan_free - Release resources held by Fresnel propagation plan
 * @plan: Propagation plan container to tear down.
 */
void atmturb_fresnel_plan_free(
    atmturb_fresnel_plan_t *plan)
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

    if (plan->tf_pri != NULL)
    {
        for (int m = 0; m < plan->nsuper; m++)
        {
            if (plan->tf_pri[m] != NULL)
            {
                fftwf_free(plan->tf_pri[m]);
            }
        }
        free(plan->tf_pri);
        plan->tf_pri = NULL;
    }

    if (plan->tf_sec != NULL)
    {
        for (int m = 0; m < plan->nsuper; m++)
        {
            if (plan->tf_sec[m] != NULL)
            {
                fftwf_free(plan->tf_sec[m]);
            }
        }
        free(plan->tf_sec);
        plan->tf_sec = NULL;
    }

    plan->nsuper = 0;
}

/**
 * atmturb_fresnel_ctx_init - Initialize thread-local Fresnel render context
 * @ctx: Thread context structure to initialize.
 * @pup_size: Linear dimension of wavefront in pixels.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_fresnel_ctx_init(
    atmturb_fresnel_ctx_t *ctx,
    long                   pup_size)
{
    if (ctx == NULL || pup_size <= 0)
    {
        return -1;
    }

    memset(ctx, 0, sizeof(*ctx));
    if (wfprop_fresnel_engine_init(&ctx->eng, pup_size) != 0)
    {
        return -1;
    }

    size_t npix    = (size_t) (pup_size * pup_size);
    ctx->field_pri = (fftwf_complex *) fftwf_alloc_complex(npix);
    ctx->field_sec = (fftwf_complex *) fftwf_alloc_complex(npix);
    ctx->super_pha = (float *) malloc(sizeof(float) * npix);
    ctx->super_spha = (float *) malloc(sizeof(float) * npix);

    if (ctx->field_pri == NULL || ctx->field_sec == NULL ||
        ctx->super_pha == NULL || ctx->super_spha == NULL)
    {
        atmturb_fresnel_ctx_free(ctx);
        return -1;
    }

    return 0;
}

/**
 * atmturb_fresnel_ctx_free - Release resources held by thread context
 * @ctx: Thread context structure to tear down.
 */
void atmturb_fresnel_ctx_free(
    atmturb_fresnel_ctx_t *ctx)
{
    if (ctx == NULL)
    {
        return;
    }

    wfprop_fresnel_engine_free(&ctx->eng);
    if (ctx->field_pri != NULL)
    {
        fftwf_free(ctx->field_pri);
        ctx->field_pri = NULL;
    }
    if (ctx->field_sec != NULL)
    {
        fftwf_free(ctx->field_sec);
        ctx->field_sec = NULL;
    }
    free(ctx->super_pha);
    free(ctx->super_spha);
    ctx->super_pha  = NULL;
    ctx->super_spha = NULL;
}

/**
 * atmturb_fresnel_extract_unwrapped - Extract amplitude and unwrapped phase using reference
 * @field: Diffracted complex optical field at pupil.
 * @pha: Phase slice holding accumulated geometric phase (updated in-place).
 * @amp: Destination amplitude slice.
 * @npix: Total number of pixels.
 */
static void atmturb_fresnel_extract_unwrapped(
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
        amp[i]   = sqrtf(re * re + im * im);

        float phi_geo = pha[i];
        float s_geo, c_geo;
        sincosf(phi_geo, &s_geo, &c_geo);

        float d_re = re * c_geo + im * s_geo;
        float d_im = im * c_geo - re * s_geo;
        float dphi = atan2f(d_im, d_re);
        pha[i]     = phi_geo + dphi;
    }
}

/**
 * atmturb_fresnel_modulate_field - Modulate complex optical field with phase screen
 * @field: In-out complex field array.
 * @phase: Phase screen array in radians.
 * @npix: Total number of pixels.
 */
static void atmturb_fresnel_modulate_field(
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
 * atmturb_fresnel_render_step - Render one frame using multi-layer diffractive propagation
 * @ctx: Thread-local execution context.
 * @plan: Precomputed diffractive propagation plan.
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @t: Simulation frame index.
 * @time_step_s: Time step between frames in seconds.
 * @master_size: Master screen dimension in pixels.
 * @pup_size: Linear dimension of wavefront in pixels.
 * @pha_slice: Destination primary phase slice (accumulates unwrapped diffractive phase).
 * @amp_slice: Destination primary amplitude slice.
 * @spha_slice: Destination secondary phase slice (accumulates unwrapped diffractive phase).
 * @samp_slice: Destination secondary amplitude slice.
 */
void atmturb_fresnel_render_step(
    atmturb_fresnel_ctx_t        *ctx,
    const atmturb_fresnel_plan_t *plan,
    const atmturb_rolling_t      *r,
    const atmturb_geom_t         *geom,
    long                          t,
    double                        time_step_s,
    long                          master_size,
    long                          pup_size,
    float                        *pha_slice,
    float                        *amp_slice,
    float                        *spha_slice,
    float                        *samp_slice)
{
    long npix = pup_size * pup_size;

    for (long i = 0; i < npix; i++)
    {
        ctx->field_pri[i][0] = 1.0f;
        ctx->field_pri[i][1] = 0.0f;
        ctx->field_sec[i][0] = 1.0f;
        ctx->field_sec[i][1] = 0.0f;
        pha_slice[i]         = 0.0f;
        spha_slice[i]        = 0.0f;
    }

    for (int m = 0; m < plan->nsuper; m++)
    {
        memset(ctx->super_pha, 0, sizeof(float) * (size_t) npix);
        memset(ctx->super_spha, 0, sizeof(float) * (size_t) npix);

        const atmturb_superlayer_t *sl = &plan->supers[m];
        for (int j = 0; j < sl->nlayers; j++)
        {
            int k = sl->layer_indices[j];
            atmturb_wfs_render_layer(r, geom, k, t, time_step_s, master_size,
                                     pup_size, ctx->super_pha, ctx->super_spha);
        }

        atmturb_add_float_array(pha_slice, ctx->super_pha, npix);
        atmturb_add_float_array(spha_slice, ctx->super_spha, npix);

        atmturb_fresnel_modulate_field(ctx->field_pri, ctx->super_pha, npix);
        wfprop_fresnel_engine_apply(&ctx->eng, ctx->field_pri, plan->tf_pri[m]);

        if (plan->tf_sec != NULL && plan->tf_sec[m] != NULL)
        {
            atmturb_fresnel_modulate_field(ctx->field_sec, ctx->super_spha, npix);
            wfprop_fresnel_engine_apply(&ctx->eng, ctx->field_sec, plan->tf_sec[m]);
        }
    }

    atmturb_fresnel_extract_unwrapped(ctx->field_pri, pha_slice, amp_slice, npix);
    if (plan->tf_sec != NULL && plan->tf_sec[0] != NULL)
    {
        atmturb_fresnel_extract_unwrapped(ctx->field_sec, spha_slice, samp_slice, npix);
    }
}
