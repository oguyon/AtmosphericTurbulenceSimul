// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_fresnel.c
 * @brief   Super-layer binning and separable 1D-chirp Fresnel propagation implementation
 */

#define _GNU_SOURCE
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "atmturb_fft_wisdom.h"
#include "atmturb_fresnel.h"
#include "atmturb_fresnel_simd.h"
#include "atmturb_simd.h"
#include "atmturb_superlayer.h"
#include "atmturb_types.h"
#include "atmturb_wfs_render.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * atmturb_fresnel_build_tf1d - Build 1D separable paraxial Fresnel chirp factor
 * @tf_1d: Destination complex vector of length @n.
 * @n: Grid linear dimension in pixels.
 * @pixscale_m: Physical pixel scale [m/pixel].
 * @z_m: Propagation step distance [m].
 * @lambda_m: Optical wavelength [m].
 *
 * Factors the 2D paraxial angular-spectrum kernel into H_2D[j, i] = H_1D[j] * H_1D[i]
 * with H_1D[k] = (1 / n) * exp(-i * pi * lambda * z * f_k^2 / L^2).
 */
static void atmturb_fresnel_build_tf1d(
    fftwf_complex *tf_1d,
    long           n,
    double         pixscale_m,
    double         z_m,
    double         lambda_m)
{
    double l_grid = (double) n * pixscale_m;
    double coeff  = (l_grid > 0.0) ? (M_PI * z_m * lambda_m / (l_grid * l_grid)) : 0.0;
    double inv_n  = 1.0 / (double) n;
    long   half   = (n + 1) / 2;

    for (long k = 0; k < n; k++)
    {
        long   fk    = (k < half) ? k : (k - n);
        double theta = -coeff * (double) (fk * fk);
        tf_1d[k][0]  = (float) (inv_n * cos(theta));
        tf_1d[k][1]  = (float) (inv_n * sin(theta));
    }
}

/**
 * atmturb_fresnel_build_sponge1d - Build 1D quartic super-Gaussian rim absorber
 * @n: Padded grid linear dimension in pixels.
 * @w_sp: Absorbing rim width in pixels (> 0).
 *
 * Return: Newly allocated float array of length @n, or NULL on allocation failure.
 */
static float *atmturb_fresnel_build_sponge1d(
    long n,
    long w_sp)
{
    float *sp = (float *) malloc(sizeof(float) * (size_t) n);
    if (sp == NULL)
    {
        return NULL;
    }

    for (long i = 0; i < n; i++)
    {
        long d_lo = i;
        long d_hi = (n - 1) - i;
        long d    = (d_lo < d_hi) ? d_lo : d_hi;
        if (d < w_sp)
        {
            double u  = 1.0 - (double) d / (double) w_sp;
            double u2 = u * u;
            sp[i]     = (float) exp(-3.0 * u2 * u2);
        }
        else
        {
            sp[i] = 1.0f;
        }
    }
    return sp;
}

/**
 * atmturb_fresnel_plan_init - Initialize super-layer binning and 1D chirp transfer functions
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
    long guard         = (CONF_FRESNEL_GUARD_PIX > 0) ? (long) CONF_FRESNEL_GUARD_PIX : 0;
    long pad_size      = pup_size + 2 * guard;
    plan->grid_size    = pup_size;
    plan->guard_pix    = guard;
    plan->pad_size     = pad_size;
    plan->sponge_width = (guard >= 12) ? (guard / 3) : 0;
    plan->pixscale_m   = pixscale_m;
    plan->lambda_ref_m = lambda_ref_m;
    plan->lambda_s_m   = lambda_s_m;

    if (atmturb_superlayer_build(&plan->supers, &plan->nsuper, prof, geom,
                                 z_bin_m, pixscale_m) != 0)
    {
        return -1;
    }

    if (plan->sponge_width > 0)
    {
        plan->sponge_1d = atmturb_fresnel_build_sponge1d(pad_size, plan->sponge_width);
        if (plan->sponge_1d == NULL)
        {
            atmturb_fresnel_plan_free(plan);
            return -1;
        }
    }

    int nsuper   = plan->nsuper;
    plan->tf_pri = (fftwf_complex **) calloc((size_t) nsuper, sizeof(fftwf_complex *));
    plan->tf_sec = (fftwf_complex **) calloc((size_t) nsuper, sizeof(fftwf_complex *));
    if (plan->tf_pri == NULL || plan->tf_sec == NULL)
    {
        atmturb_fresnel_plan_free(plan);
        return -1;
    }

    for (int m = 0; m < nsuper; m++)
    {
        plan->tf_pri[m] = fftwf_alloc_complex((size_t) pad_size);
        if (plan->tf_pri[m] == NULL)
        {
            atmturb_fresnel_plan_free(plan);
            return -1;
        }
        atmturb_fresnel_build_tf1d(plan->tf_pri[m], pad_size, pixscale_m,
                                   plan->supers[m].step_dist_m, lambda_ref_m);

        if (lambda_s_m > 0.0)
        {
            plan->tf_sec[m] = fftwf_alloc_complex((size_t) pad_size);
            if (plan->tf_sec[m] == NULL)
            {
                atmturb_fresnel_plan_free(plan);
                return -1;
            }
            atmturb_fresnel_build_tf1d(plan->tf_sec[m], pad_size, pixscale_m,
                                       plan->supers[m].step_dist_s_m, lambda_s_m);
        }
    } // for (int m = 0; m < nsuper; m++)

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
            fftwf_free(plan->tf_pri[m]);
        }
        free(plan->tf_pri);
        plan->tf_pri = NULL;
    }

    if (plan->tf_sec != NULL)
    {
        for (int m = 0; m < plan->nsuper; m++)
        {
            fftwf_free(plan->tf_sec[m]);
        }
        free(plan->tf_sec);
        plan->tf_sec = NULL;
    }

    free(plan->sponge_1d);
    plan->sponge_1d = NULL;
    plan->nsuper    = 0;
}

/**
 * atmturb_fresnel_init_fftw_threads - Configure FFTW thread count for Fresnel propagation
 * @grid_size: Padded grid linear dimension in pixels.
 * @req_nthreads: Explicit thread count (>0), or 0 for automatic selection.
 */
static void atmturb_fresnel_init_fftw_threads(
    long grid_size,
    int  req_nthreads)
{
#ifdef _OPENMP
    static int s_threads_init = 0;
    if (!s_threads_init)
    {
        #pragma omp critical(atmturb_fresnel_fftw_init)
        {
            if (!s_threads_init && fftwf_init_threads())
            {
                s_threads_init = 1;
            }
        }
    }
    if (s_threads_init)
    {
        int nth = (req_nthreads > 0) ? req_nthreads : 1;
        if (req_nthreads <= 0 && !omp_in_parallel())
        {
            int max_th = omp_get_max_threads();
            if (grid_size > 128 && grid_size <= 256)
            {
                nth = (max_th < 4) ? max_th : 4;
            }
            else if (grid_size > 256 && grid_size <= 512)
            {
                nth = (max_th < 8) ? max_th : 8;
            }
            else if (grid_size > 512)
            {
                nth = (max_th < 16) ? max_th : 16;
            }
        }
        fftwf_plan_with_nthreads(nth);
    }
#else
    (void) grid_size;
    (void) req_nthreads;
#endif
}

/**
 * atmturb_fresnel_ctx_init_mt - Initialize Fresnel context with explicit threading/shared plans
 * @ctx: Thread context structure to initialize.
 * @pup_size: Linear dimension of propagation grid (`plan->pad_size`).
 * @fftw_nthreads: Explicit FFTW thread count (>0), or 0 for automatic sizing.
 * @shared_src: Optional donor context whose FFTW plans are reused read-only (or NULL).
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_fresnel_ctx_init_mt(
    atmturb_fresnel_ctx_t       *ctx,
    long                         pup_size,
    int                          fftw_nthreads,
    const atmturb_fresnel_ctx_t *shared_src)
{
    if (ctx == NULL || pup_size <= 0)
    {
        return -1;
    }

    memset(ctx, 0, sizeof(*ctx));
    ctx->eng.n     = pup_size;
    ctx->eng_sec.n = pup_size;

    size_t npix     = (size_t) pup_size * (size_t) pup_size;
    ctx->field_pri  = fftwf_alloc_complex(npix);
    ctx->field_sec  = fftwf_alloc_complex(npix);
    ctx->super_pha  = fftwf_alloc_real(npix);
    ctx->super_spha = fftwf_alloc_real(npix);

    if (ctx->field_pri == NULL || ctx->field_sec == NULL ||
        ctx->super_pha == NULL || ctx->super_spha == NULL)
    {
        atmturb_fresnel_ctx_free(ctx);
        return -1;
    }

    if (shared_src != NULL)
    {
        ctx->eng        = shared_src->eng;
        ctx->eng_sec    = shared_src->eng_sec;
        ctx->owns_plans = 0;
        return 0;
    }

    atmturb_fftwf_wisdom_load();
    atmturb_fresnel_init_fftw_threads(pup_size, fftw_nthreads);
    unsigned int flags = atmturb_fftwf_plan_flags();
    ctx->eng.fwd_plan  = fftwf_plan_dft_2d((int) pup_size, (int) pup_size,
                                           ctx->field_pri, ctx->field_pri,
                                           FFTW_FORWARD, flags);
    ctx->eng.inv_plan  = fftwf_plan_dft_2d((int) pup_size, (int) pup_size,
                                           ctx->field_pri, ctx->field_pri,
                                           FFTW_BACKWARD, flags);
    ctx->owns_plans    = 1;

    if (ctx->eng.fwd_plan == NULL || ctx->eng.inv_plan == NULL)
    {
        atmturb_fresnel_ctx_free(ctx);
        return -1;
    }

    ctx->eng_sec.fwd_plan = ctx->eng.fwd_plan;
    ctx->eng_sec.inv_plan = ctx->eng.inv_plan;

    if (flags == FFTW_MEASURE)
    {
        atmturb_fftwf_wisdom_save();
    }

    return 0;
}

/**
 * atmturb_fresnel_ctx_init - Allocate thread-local scratchpad and init FFTW engine
 * @ctx: Thread context structure to initialize.
 * @pup_size: Linear dimension of propagation grid (`plan->pad_size`).
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_fresnel_ctx_init(
    atmturb_fresnel_ctx_t *ctx,
    long                   pup_size)
{
    return atmturb_fresnel_ctx_init_mt(ctx, pup_size, 0, NULL);
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

    if (ctx->owns_plans)
    {
        wfprop_fresnel_engine_free(&ctx->eng);
    }
    memset(&ctx->eng, 0, sizeof(ctx->eng));
    memset(&ctx->eng_sec, 0, sizeof(ctx->eng_sec));
    ctx->owns_plans = 0;
    fftwf_free(ctx->field_pri);
    fftwf_free(ctx->field_sec);
    fftwf_free(ctx->super_pha);
    fftwf_free(ctx->super_spha);
    ctx->field_pri  = NULL;
    ctx->field_sec  = NULL;
    ctx->super_pha  = NULL;
    ctx->super_spha = NULL;
}

/**
 * atmturb_fresnel_apply_sponge - Attenuate outer guard rim of complex field via 1D window
 */
static void atmturb_fresnel_apply_sponge(
    fftwf_complex       *restrict field,
    const float         *restrict sp_1d,
    long                          pad_size,
    long                          w_sp)
{
    for (long y = 0; y < pad_size; y++)
    {
        fftwf_complex *restrict row = &field[y * pad_size];
        if (y < w_sp || y >= pad_size - w_sp)
        {
            float wy = sp_1d[y];
            #pragma omp simd
            for (long x = 0; x < pad_size; x++)
            {
                float w   = wy * sp_1d[x];
                row[x][0] *= w;
                row[x][1] *= w;
            }
        }
        else
        {
            for (long x = 0; x < w_sp; x++)
            {
                float wl          = sp_1d[x];
                long  xr          = pad_size - 1 - x;
                float wr          = sp_1d[xr];
                row[x][0]        *= wl;
                row[x][1]        *= wl;
                row[xr][0]       *= wr;
                row[xr][1]       *= wr;
            }
        }
    } // for (long y = 0; y < pad_size; y++)
}

/**
 * atmturb_fresnel_mul_tf1d_rows - Multiply 2D spectrum rows [y0, y1) by separable 1D chirp
 */
static inline void atmturb_fresnel_mul_tf1d_rows(
    fftwf_complex               *restrict field,
    const fftwf_complex         *restrict tf_1d,
    long                                  n,
    long                                  y0,
    long                                  y1)
{
    for (long y = y0; y < y1; y++)
    {
        float                   hy_re = tf_1d[y][0];
        float                   hy_im = tf_1d[y][1];
        fftwf_complex *restrict row   = &field[y * n];
        #pragma omp simd
        for (long x = 0; x < n; x++)
        {
            float hx_re = tf_1d[x][0];
            float hx_im = tf_1d[x][1];
            float t_re  = hy_re * hx_re - hy_im * hx_im;
            float t_im  = hy_re * hx_im + hy_im * hx_re;
            float f_re  = row[x][0];
            float f_im  = row[x][1];
            row[x][0]   = f_re * t_re - f_im * t_im;
            row[x][1]   = f_re * t_im + f_im * t_re;
        }
    }
}

/**
 * atmturb_fresnel_propagate_field - In-place diffractive propagation using FFTW and 1D chirp
 */
static inline void atmturb_fresnel_propagate_field(
    const wfprop_fresnel_engine_t *eng,
    fftwf_complex                 *field,
    const fftwf_complex           *tf_1d,
    long                           pad_size)
{
    fftwf_execute_dft(eng->fwd_plan, field, field);
#ifdef _OPENMP
    if (pad_size * pad_size >= 262144 && !omp_in_parallel())
    {
        int max_t = omp_get_max_threads();
        int nth   = (max_t > 8) ? 8 : max_t;
        #pragma omp parallel num_threads(nth)
        {
            int  tid = omp_get_thread_num();
            int  nt  = omp_get_num_threads();
            long y0  = (pad_size * (long) tid) / (long) nt;
            long y1  = (pad_size * (long) (tid + 1)) / (long) nt;
            atmturb_fresnel_mul_tf1d_rows(field, tf_1d, pad_size, y0, y1);
        }
    }
    else
#endif
    {
        atmturb_fresnel_mul_tf1d_rows(field, tf_1d, pad_size, 0, pad_size);
    }
    fftwf_execute_dft(eng->inv_plan, field, field);
}

/**
 * atmturb_fresnel_step_propagate - Execute single or dual-channel Fresnel FFT step
 */
static void atmturb_fresnel_step_propagate(
    atmturb_fresnel_ctx_t        *ctx,
    const atmturb_fresnel_plan_t *plan,
    int                           m,
    int                           do_sec_prop,
    long                          pad_size)
{
    int par_sec = (do_sec_prop && pad_size < 512);
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
                                                plan->tf_pri[m], pad_size);
            }
            #pragma omp section
            {
                atmturb_fresnel_propagate_field(&ctx->eng_sec, ctx->field_sec,
                                                plan->tf_sec[m], pad_size);
            }
        }
        return;
    }

    atmturb_fresnel_propagate_field(&ctx->eng, ctx->field_pri, plan->tf_pri[m], pad_size);
    if (do_sec_prop)
    {
        atmturb_fresnel_propagate_field(&ctx->eng_sec, ctx->field_sec,
                                        plan->tf_sec[m], pad_size);
    }
}

/**
 * atmturb_fresnel_extract_channel - Extract unwrapped phase/amp from padded or unpadded field
 */
static void atmturb_fresnel_extract_channel(
    const fftwf_complex *field,
    float               *pha_slice,
    float               *amp_slice,
    long                 pad_size,
    long                 pup_size,
    long                 guard)
{
    if (guard <= 0)
    {
        atmturb_fresnel_simd_extract_unwrapped(field, pha_slice, amp_slice, pup_size * pup_size);
        return;
    }

    for (long y = 0; y < pup_size; y++)
    {
        const fftwf_complex *f_row = &field[(y + guard) * pad_size + guard];
        float               *p_row = &pha_slice[y * pup_size];
        float               *a_row = (amp_slice != NULL) ? &amp_slice[y * pup_size] : NULL;
        atmturb_fresnel_simd_extract_unwrapped(f_row, p_row, a_row, pup_size);
    }
}

/**
 * atmturb_fresnel_accum_crop - Accumulate cropped interior of padded screen into pupil slice
 */
static void atmturb_fresnel_accum_crop(
    float       *dst_pup,
    const float *src_pad,
    long         pad_size,
    long         pup_size,
    long         guard)
{
    if (dst_pup == NULL)
    {
        return;
    }
    for (long y = 0; y < pup_size; y++)
    {
        atmturb_add_float_array(&dst_pup[y * pup_size],
                                &src_pad[(y + guard) * pad_size + guard], pup_size);
    }
}

/**
 * atmturb_fresnel_modulate_step - Modulate primary/secondary fields for superlayer m
 */
static void atmturb_fresnel_modulate_step(
    atmturb_fresnel_ctx_t        *ctx,
    const atmturb_fresnel_plan_t *plan,
    int                           m,
    int                           do_sec_prop,
    long                          pad_size,
    long                          pup_size,
    long                          guard,
    float                        *pha_slice,
    float                        *spha_slice)
{
    long   pad_npix = pad_size * pad_size;
    float *cur_pha  = (m == 0 && guard == 0) ? pha_slice : ctx->super_pha;
    float *cur_sec  = (m == 0 && guard == 0 && spha_slice != NULL) ? spha_slice : ctx->super_spha;

    if (guard > 0)
    {
        atmturb_fresnel_accum_crop(pha_slice, cur_pha, pad_size, pup_size, guard);
    }
    else if (m > 0)
    {
        atmturb_add_float_array(pha_slice, cur_pha, pad_npix);
    }

    if (m == 0)
    {
        atmturb_fresnel_simd_init_modulate_field(ctx->field_pri, cur_pha, pad_npix);
    }
    else
    {
        atmturb_fresnel_simd_modulate_field(ctx->field_pri, cur_pha, pad_npix);
    }
    if (plan->sponge_1d != NULL && guard > 0)
    {
        atmturb_fresnel_apply_sponge(ctx->field_pri, plan->sponge_1d,
                                     pad_size, plan->sponge_width);
    }

    if (!do_sec_prop)
    {
        return;
    }

    if (guard > 0)
    {
        atmturb_fresnel_accum_crop(spha_slice, cur_sec, pad_size, pup_size, guard);
    }
    else if (m > 0 && spha_slice != NULL)
    {
        atmturb_add_float_array(spha_slice, cur_sec, pad_npix);
    }

    if (m == 0)
    {
        atmturb_fresnel_simd_init_modulate_field(ctx->field_sec, cur_sec, pad_npix);
    }
    else
    {
        atmturb_fresnel_simd_modulate_field(ctx->field_sec, cur_sec, pad_npix);
    }
    if (plan->sponge_1d != NULL && guard > 0)
    {
        atmturb_fresnel_apply_sponge(ctx->field_sec, plan->sponge_1d,
                                     pad_size, plan->sponge_width);
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
    (void) master_size;
    long guard    = (ctx->eng.n == plan->pad_size) ? plan->guard_pix : 0;
    long pad_size = pup_size + 2 * guard;
    long pup_npix = pup_size * pup_size;
    long pad_npix = pad_size * pad_size;
    int  has_sec  = (plan->lambda_s_m > 0.0 && plan->tf_sec != NULL &&
                     plan->tf_sec[0] != NULL && ctx->field_sec != NULL &&
                     ctx->super_spha != NULL && (spha_slice != NULL || samp_slice != NULL));

    memset(pha_slice, 0, sizeof(float) * (size_t) pup_npix);
    if (has_sec && spha_slice != NULL)
    {
        memset(spha_slice, 0, sizeof(float) * (size_t) pup_npix);
    }

    int has_trail_zero = (plan->nsuper > 0 &&
                          fabs(plan->supers[plan->nsuper - 1].step_dist_m) <= 0.0);
    int n_prop         = has_trail_zero ? (plan->nsuper - 1) : plan->nsuper;

    if (n_prop <= 0)
    {
        if (has_trail_zero)
        {
            atmturb_wfs_render_target_t zt = {
                .pup_size = pup_size, .guard_pix = 0,
                .pha = pha_slice, .spha = (has_sec ? spha_slice : NULL), .weight_scale = 1.0f
            };
            atmturb_wfs_render_superlayer(r, geom, &plan->supers[0], t, time_step_s, &zt);
        }
        atmturb_init_phase_amp(NULL, amp_slice, pup_npix);
        if (has_sec)
        {
            atmturb_init_phase_amp(NULL, samp_slice, pup_npix);
        }
        return;
    }

    for (int m = 0; m < n_prop; m++)
    {
        float *cur_pha = (m == 0 && guard == 0) ? pha_slice : ctx->super_pha;
        float *cur_sec = !has_sec ? NULL
                                  : ((m == 0 && guard == 0 && spha_slice != NULL)
                                         ? spha_slice
                                         : ctx->super_spha);
        if (cur_pha == ctx->super_pha)
        {
            memset(ctx->super_pha, 0, sizeof(float) * (size_t) pad_npix);
        }
        if (cur_sec == ctx->super_spha)
        {
            memset(ctx->super_spha, 0, sizeof(float) * (size_t) pad_npix);
        }

        atmturb_wfs_render_target_t target = {
            .pup_size = pup_size, .guard_pix = guard,
            .pha = cur_pha, .spha = cur_sec, .weight_scale = 1.0f
        };
        atmturb_wfs_render_superlayer(r, geom, &plan->supers[m], t, time_step_s, &target);

        int do_sec_prop = (has_sec && plan->tf_sec[m] != NULL);
        atmturb_fresnel_modulate_step(ctx, plan, m, do_sec_prop, pad_size, pup_size,
                                      guard, pha_slice, spha_slice);
        if (fabs(plan->supers[m].step_dist_m) > 0.0)
        {
            atmturb_fresnel_step_propagate(ctx, plan, m, do_sec_prop, pad_size);
        }
    } // for (int m = 0; m < n_prop; m++)

    float *sec_pha_ext = (has_sec && spha_slice == NULL) ? ctx->super_spha : spha_slice;
    if (has_sec && spha_slice == NULL)
    {
        memset(ctx->super_spha, 0, sizeof(float) * (size_t) pup_npix);
    }
    atmturb_fresnel_extract_channel(ctx->field_pri, pha_slice, amp_slice,
                                    pad_size, pup_size, guard);
    if (has_sec)
    {
        atmturb_fresnel_extract_channel(ctx->field_sec, sec_pha_ext, samp_slice,
                                        pad_size, pup_size, guard);
    }

    if (has_trail_zero)
    {
        atmturb_wfs_render_target_t zt = {
            .pup_size     = pup_size,
            .guard_pix    = 0,
            .pha          = pha_slice,
            .spha         = (has_sec && spha_slice != NULL) ? spha_slice : NULL,
            .weight_scale = 1.0f
        };
        atmturb_wfs_render_superlayer(r, geom, &plan->supers[plan->nsuper - 1],
                                      t, time_step_s, &zt);
    }
}
