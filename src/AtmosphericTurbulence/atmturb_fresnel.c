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
#include "atmturb_fresnel_simd.h"
#include "atmturb_simd.h"
#include "atmturb_superlayer.h"
#include "atmturb_types.h"

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
 * atmturb_fresnel_init_fftw_threads - Initialize OpenMP thread pool for FFTW plans
 * @pup_size: Linear dimension of wavefront in pixels.
 */
static void atmturb_fresnel_init_fftw_threads(
    long pup_size)
{
#ifdef _OPENMP
    static int s_threads_init = 0;
    if (!s_threads_init)
    {
        #pragma omp critical
        {
            if (!s_threads_init)
            {
                if (fftwf_init_threads())
                {
                    s_threads_init = 1;
                }
            }
        }
    }
    if (s_threads_init)
    {
        int nth = 1;
#ifdef _OPENMP
        if (!omp_in_parallel())
        {
            int max_th = omp_get_max_threads();
            if (pup_size > 128 && pup_size <= 256)
            {
                nth = (max_th < 4) ? max_th : 4;
            }
            else if (pup_size > 256 && pup_size <= 512)
            {
                nth = (max_th < 8) ? max_th : 8;
            }
            else if (pup_size > 512)
            {
                nth = (max_th < 16) ? max_th : 16;
            }
        }
#endif
        fftwf_plan_with_nthreads(nth);
    }
#else
    (void) pup_size;
#endif
}

/**
 * atmturb_fresnel_ctx_init - Allocate thread-local scratchpad and init FFTW engine
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
    ctx->eng.n = pup_size;

    size_t npix    = (size_t) (pup_size * pup_size);
    ctx->field_pri = (fftwf_complex *) fftwf_alloc_complex(npix);
    ctx->super_pha = (float *) malloc(sizeof(float) * npix);

    if (ctx->field_pri == NULL || ctx->super_pha == NULL)
    {
        atmturb_fresnel_ctx_free(ctx);
        return -1;
    }

    int has_sec = (CONF_MAKE_SWAVEFRONT == 1);
    if (has_sec)
    {
        ctx->eng_sec.n  = pup_size;
        ctx->field_sec  = (fftwf_complex *) fftwf_alloc_complex(npix);
        ctx->super_spha = (float *) malloc(sizeof(float) * npix);
        if (ctx->field_sec == NULL || ctx->super_spha == NULL)
        {
            atmturb_fresnel_ctx_free(ctx);
            return -1;
        }
    }

    atmturb_fresnel_init_fftw_threads(pup_size);
    ctx->eng.fwd_plan = fftwf_plan_dft_2d((int) pup_size, (int) pup_size,
                                          ctx->field_pri, ctx->field_pri,
                                          FFTW_FORWARD, FFTW_MEASURE);
    ctx->eng.inv_plan = fftwf_plan_dft_2d((int) pup_size, (int) pup_size,
                                          ctx->field_pri, ctx->field_pri,
                                          FFTW_BACKWARD, FFTW_MEASURE);

    if (has_sec)
    {
        ctx->eng_sec.fwd_plan = fftwf_plan_dft_2d((int) pup_size, (int) pup_size,
                                                  ctx->field_sec, ctx->field_sec,
                                                  FFTW_FORWARD, FFTW_MEASURE);
        ctx->eng_sec.inv_plan = fftwf_plan_dft_2d((int) pup_size, (int) pup_size,
                                                  ctx->field_sec, ctx->field_sec,
                                                  FFTW_BACKWARD, FFTW_MEASURE);
    }

    if (ctx->eng.fwd_plan == NULL || ctx->eng.inv_plan == NULL ||
        (has_sec && (ctx->eng_sec.fwd_plan == NULL || ctx->eng_sec.inv_plan == NULL)))
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
    wfprop_fresnel_engine_free(&ctx->eng_sec);
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
 * atmturb_fresnel_propagate_field - In-place diffractive propagation using FFTW and SIMD TF
 * @eng: Initialized FFTW Fresnel propagation engine.
 * @field: Optical complex field buffer.
 * @tf: Precomputed 2D complex transfer function.
 * @npix: Total number of pixels.
 */
static inline void atmturb_fresnel_propagate_field(
    const wfprop_fresnel_engine_t *eng,
    fftwf_complex                 *field,
    const fftwf_complex           *tf,
    long                           npix)
{
    fftwf_execute_dft(eng->fwd_plan, field, field);
    atmturb_complex_mul_array((float *) field, (const float *) field, (const float *) tf, npix);
    fftwf_execute_dft(eng->inv_plan, field, field);
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
    int has_sec = (plan->lambda_s_m > 0.0 && plan->tf_sec != NULL &&
                   plan->tf_sec[0] != NULL && (spha_slice != NULL || samp_slice != NULL));

    memset(pha_slice, 0, sizeof(float) * (size_t) npix);
    if (has_sec && spha_slice != NULL)
    {
        memset(spha_slice, 0, sizeof(float) * (size_t) npix);
    }

    float *sec_target = has_sec ? ctx->super_spha : NULL;

    for (int m = 0; m < plan->nsuper; m++)
    {
        memset(ctx->super_pha, 0, sizeof(float) * (size_t) npix);
        if (has_sec)
        {
            memset(ctx->super_spha, 0, sizeof(float) * (size_t) npix);
        }

        const atmturb_superlayer_t *sl = &plan->supers[m];
        for (int j = 0; j < sl->nlayers; j++)
        {
            int k = sl->layer_indices[j];
            atmturb_wfs_render_layer(r, geom, k, t, time_step_s, master_size,
                                     pup_size, ctx->super_pha, sec_target);
        }

        atmturb_add_float_array(pha_slice, ctx->super_pha, npix);
        if (m == 0)
        {
            atmturb_fresnel_simd_init_modulate_field(ctx->field_pri, ctx->super_pha, npix);
        }
        else
        {
            atmturb_fresnel_simd_modulate_field(ctx->field_pri, ctx->super_pha, npix);
        }

        int do_sec_prop = (has_sec && plan->tf_sec[m] != NULL);
        if (do_sec_prop)
        {
            if (spha_slice != NULL)
            {
                atmturb_add_float_array(spha_slice, ctx->super_spha, npix);
            }
            if (m == 0)
            {
                atmturb_fresnel_simd_init_modulate_field(ctx->field_sec, ctx->super_spha, npix);
            }
            else
            {
                atmturb_fresnel_simd_modulate_field(ctx->field_sec, ctx->super_spha, npix);
            }
        }

        int par_sec = (do_sec_prop && pup_size < 512);
#ifdef _OPENMP
        par_sec = par_sec && !omp_in_parallel();
#endif
        if (par_sec)
        {
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    atmturb_fresnel_propagate_field(&ctx->eng, ctx->field_pri,
                                                    plan->tf_pri[m], npix);
                }
                #pragma omp section
                {
                    atmturb_fresnel_propagate_field(&ctx->eng_sec, ctx->field_sec,
                                                    plan->tf_sec[m], npix);
                }
            }
        }
        else
        {
            atmturb_fresnel_propagate_field(&ctx->eng, ctx->field_pri,
                                            plan->tf_pri[m], npix);
            if (do_sec_prop)
            {
                atmturb_fresnel_propagate_field(&ctx->eng_sec, ctx->field_sec,
                                                plan->tf_sec[m], npix);
            }
        }
    }

    int par_extract = (has_sec && pup_size < 512);
#ifdef _OPENMP
    par_extract = par_extract && !omp_in_parallel();
#endif
    if (par_extract)
    {
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                atmturb_fresnel_simd_extract_unwrapped(ctx->field_pri, pha_slice, amp_slice, npix);
            }
            #pragma omp section
            {
                atmturb_fresnel_simd_extract_unwrapped(ctx->field_sec, spha_slice,
                                                       samp_slice, npix);
            }
        }
    }
    else
    {
        atmturb_fresnel_simd_extract_unwrapped(ctx->field_pri, pha_slice, amp_slice, npix);
        if (has_sec)
        {
            atmturb_fresnel_simd_extract_unwrapped(ctx->field_sec, spha_slice,
                                                   samp_slice, npix);
        }
    }
}
