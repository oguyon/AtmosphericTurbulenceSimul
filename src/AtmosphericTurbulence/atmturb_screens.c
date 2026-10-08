// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_screens.c
 * @brief   Master atmospheric turbulence phase screen generation
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "AtmosphericTurbulence.h"
#include "atmturb_types.h"

/**
 * atmturb_fill_spectrum_float - Fill 2D Fourier spectrum with shaped Gaussian noise (float)
 * @buf: Pointer to FFTW single-precision complex frequency buffer.
 * @spec: Pointer to screen generation specifications.
 */
static void atmturb_fill_spectrum_float(
    fftwf_complex               *buf,
    const atmturb_screen_spec_t *spec)
{
    long size = spec->size;
    float r0_pix = (spec->r0_pix > 0.0) ? (float) spec->r0_pix : powf(6.88f, 0.6f);
    float k0 = (spec->L0_pix > 0.0) ? ((float) size / (float) spec->L0_pix) : 0.0f;
    float km = (spec->l0_pix > 0.0)
                    ? ((5.92f / (2.0f * (float) M_PI)) * (float) size / (float) spec->l0_pix)
                    : 0.0f;
    float pref = 0.023f * powf(r0_pix, -5.0f / 3.0f) * powf((float) size, 5.0f / 3.0f);
    float sqrt_pref = sqrtf(pref);
    float inv_two_km2 = (km > 0.0f) ? (0.5f / (km * km)) : 0.0f;
    float k0_sq = k0 * k0;
    uint64_t base_seed = atmturb_resolve_seed(spec->seed);

    #pragma omp parallel for schedule(static)
    for (long jj = 0; jj < size; jj++)
    {
        uint64_t rng = atmturb_rng_stream_seed(base_seed, (uint64_t) jj);
        long fy = (jj < size / 2) ? jj : (jj - size);
        float fy2 = (float) (fy * fy);
        long row = jj * size;

        for (long ii = 0; ii < size; ii++)
        {
            long fx = (ii < size / 2) ? ii : (ii - size);
            float r2 = (float) (fx * fx) + fy2;
            double g0, g1;
            atmturb_rng_gaussian_pair(&rng, &g0, &g1);

            if (fx == 0 && fy == 0)
            {
                buf[row + ii][0] = 0.0f;
                buf[row + ii][1] = 0.0f;
            }
            else
            {
                float dist2 = r2 + k0_sq;
                float amp = sqrt_pref * powf(dist2, -11.0f / 12.0f);
                if (inv_two_km2 > 0.0f)
                {
                    amp *= expf(-r2 * inv_two_km2);
                }
                buf[row + ii][0] = amp * (float) g0;
                buf[row + ii][1] = amp * (float) g1;
            }
        }
    }
}

/**
 * atmturb_fill_spectrum_double - Fill 2D Fourier spectrum with shaped Gaussian noise (double)
 * @buf: Pointer to FFTW double-precision complex frequency buffer.
 * @spec: Pointer to screen generation specifications.
 */
static void atmturb_fill_spectrum_double(
    fftw_complex                *buf,
    const atmturb_screen_spec_t *spec)
{
    long size = spec->size;
    double r0_pix = (spec->r0_pix > 0.0) ? spec->r0_pix : pow(6.88, 0.6);
    double k0 = (spec->L0_pix > 0.0) ? ((double) size / spec->L0_pix) : 0.0;
    double km = (spec->l0_pix > 0.0)
                    ? ((5.92 / (2.0 * M_PI)) * (double) size / spec->l0_pix)
                    : 0.0;
    double pref = 0.023 * pow(r0_pix, -5.0 / 3.0) * pow((double) size, 5.0 / 3.0);
    double sqrt_pref = sqrt(pref);
    double inv_two_km2 = (km > 0.0) ? (0.5 / (km * km)) : 0.0;
    double k0_sq = k0 * k0;
    uint64_t base_seed = atmturb_resolve_seed(spec->seed);

    #pragma omp parallel for schedule(static)
    for (long jj = 0; jj < size; jj++)
    {
        uint64_t rng = atmturb_rng_stream_seed(base_seed, (uint64_t) jj);
        long fy = (jj < size / 2) ? jj : (jj - size);
        double fy2 = (double) (fy * fy);
        long row = jj * size;

        for (long ii = 0; ii < size; ii++)
        {
            long fx = (ii < size / 2) ? ii : (ii - size);
            double r2 = (double) (fx * fx) + fy2;
            double g0, g1;
            atmturb_rng_gaussian_pair(&rng, &g0, &g1);

            if (fx == 0 && fy == 0)
            {
                buf[row + ii][0] = 0.0;
                buf[row + ii][1] = 0.0;
            }
            else
            {
                double dist2 = r2 + k0_sq;
                double amp = sqrt_pref * pow(dist2, -11.0 / 12.0);
                if (inv_two_km2 > 0.0)
                {
                    amp *= exp(-r2 * inv_two_km2);
                }
                buf[row + ii][0] = amp * g0;
                buf[row + ii][1] = amp * g1;
            }
        }
    }
}

/**
 * atmturb_fftwf_ensure_threads - Initialize FFTW single-precision multi-threading
 */
static void atmturb_fftwf_ensure_threads(void)
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
        int nth = omp_get_max_threads();
        if (nth > 1)
        {
            fftwf_plan_with_nthreads(nth);
        }
    }
#endif
}

/**
 * atmturb_fftw_ensure_threads - Initialize FFTW double-precision multi-threading
 */
static void atmturb_fftw_ensure_threads(void)
{
#ifdef _OPENMP
    static int s_threads_init = 0;
    if (!s_threads_init)
    {
        #pragma omp critical
        {
            if (!s_threads_init)
            {
                if (fftw_init_threads())
                {
                    s_threads_init = 1;
                }
            }
        }
    }
    if (s_threads_init)
    {
        int nth = omp_get_max_threads();
        if (nth > 1)
        {
            fftw_plan_with_nthreads(nth);
        }
    }
#endif
}

/**
 * atmturb_generate_screen_pair_float - Generate screen pair using single-precision FFT
 * @spec: Pointer to screen generation specifications.
 * @screen_a: Output buffer for screen 0 (size * size floats).
 * @screen_b: Output buffer for screen 1 (size * size floats).
 *
 * Return: 0 on success, -1 on allocation failure.
 */
static int atmturb_generate_screen_pair_float(
    const atmturb_screen_spec_t *spec,
    float                       *screen_a,
    float                       *screen_b)
{
    long size = spec->size;
    long total = size * size;
    fftwf_complex *buf = (fftwf_complex *) fftwf_alloc_complex(total);
    if (!buf)
    {
        return -1;
    }

    atmturb_fill_spectrum_float(buf, spec);

    atmturb_fftwf_ensure_threads();
    fftwf_plan plan = fftwf_plan_dft_2d((int) size, (int) size, buf, buf,
                                        FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    #pragma omp parallel for schedule(static)
    for (long i = 0; i < total; i++)
    {
        if (screen_a)
        {
            screen_a[i] = buf[i][0];
        }
        if (screen_b)
        {
            screen_b[i] = buf[i][1];
        }
    }

    fftwf_free(buf);
    return 0;
}

/**
 * atmturb_generate_screen_pair_double - Generate screen pair using double-precision FFT
 * @spec: Pointer to screen generation specifications.
 * @screen_a: Output buffer for screen 0 (size * size floats).
 * @screen_b: Output buffer for screen 1 (size * size floats).
 *
 * Return: 0 on success, -1 on allocation failure.
 */
static int atmturb_generate_screen_pair_double(
    const atmturb_screen_spec_t *spec,
    float                       *screen_a,
    float                       *screen_b)
{
    long size = spec->size;
    long total = size * size;
    fftw_complex *buf = (fftw_complex *) fftw_alloc_complex(total);
    if (!buf)
    {
        return -1;
    }

    atmturb_fill_spectrum_double(buf, spec);

    atmturb_fftw_ensure_threads();
    fftw_plan plan = fftw_plan_dft_2d((int) size, (int) size, buf, buf,
                                      FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    #pragma omp parallel for schedule(static)
    for (long i = 0; i < total; i++)
    {
        if (screen_a)
        {
            screen_a[i] = (float) buf[i][0];
        }
        if (screen_b)
        {
            screen_b[i] = (float) buf[i][1];
        }
    }

    fftw_free(buf);
    return 0;
}

/**
 * atmturb_generate_screen_pair - Generate two independent phase screens from spectrum
 * @spec: Generation parameters and scales.
 * @screen_a: Output buffer of size * size floats (or NULL to skip).
 * @screen_b: Output buffer of size * size floats (or NULL to skip).
 *
 * Return: 0 on success, non-zero error code otherwise.
 */
int atmturb_generate_screen_pair(
    const atmturb_screen_spec_t *spec,
    float                       *screen_a,
    float                       *screen_b)
{
    if (!spec || spec->size <= 0)
    {
        return -1;
    }

    if (spec->precision == 0)
    {
        return atmturb_generate_screen_pair_float(spec, screen_a, screen_b);
    }
    return atmturb_generate_screen_pair_double(spec, screen_a, screen_b);
}

/**
 * atmturb_measure_r0_pix - Estimate effective r0 in pixels from structure function at lag 1
 * @data: Pointer to 2D float screen array.
 * @size: Grid dimension in pixels.
 *
 * Return: Estimated r0 in pixel units.
 */
double atmturb_measure_r0_pix(
    const float *data,
    long         size)
{
    if (!data || size < 2)
    {
        return 0.0;
    }

    double sum_sq = 0.0;
    long count = 0;

    #pragma omp parallel for reduction(+:sum_sq,count) schedule(static)
    for (long y = 0; y < size; y++)
    {
        long y_next = (y + 1) % size;
        long row = y * size;
        long row_next = y_next * size;

        for (long x = 0; x < size; x++)
        {
            long x_next = (x + 1) % size;
            double dx = (double) data[row + x_next] - (double) data[row + x];
            double dy = (double) data[row_next + x] - (double) data[row + x];
            sum_sq += dx * dx + dy * dy;
            count += 2;
        }
    }

    double d1 = sum_sq / (double) count;
    if (d1 <= 1e-30)
    {
        return 0.0;
    }
    return pow(6.88 / d1, 0.6);
}

/**
 * make_master_turbulence_screen_seeded - Generate von Karman master screens with seed
 * @ID_name1: Output name for screen 1.
 * @ID_name2: Output name for screen 2.
 * @size: Grid dimension in pixels.
 * @outerscale: Outer scale in pixels.
 * @innerscale: Inner scale in pixels.
 * @WFprecision: Precision flag (0=single, 1=double).
 * @seed: PRNG seed value (0 = time-based).
 *
 * Return: 0 on success.
 */
int make_master_turbulence_screen_seeded(
    const char *ID_name1,
    const char *ID_name2,
    long        size,
    float       outerscale,
    float       innerscale,
    long        WFprecision,
    uint64_t    seed)
{
    if (size <= 0)
    {
        return -1;
    }

    delete_image_ID(ID_name1);
    delete_image_ID(ID_name2);
    imageID ID1 = create_2Dimage_ID(ID_name1, size, size);
    imageID ID2 = create_2Dimage_ID(ID_name2, size, size);
    if (ID1 < 0 || ID2 < 0)
    {
        return -1;
    }

    atmturb_screen_spec_t spec;
    memset(&spec, 0, sizeof(spec));
    spec.size = size;
    spec.r0_pix = pow(6.88, 0.6);
    spec.L0_pix = (double) outerscale;
    spec.l0_pix = (double) innerscale;
    spec.seed = seed;
    spec.precision = (int) WFprecision;

    int ret = atmturb_generate_screen_pair(&spec, dcimg[ID1].array.F, dcimg[ID2].array.F);
    if (ret != 0)
    {
        delete_image_ID(ID_name1);
        delete_image_ID(ID_name2);
        return ret;
    }

    double r0_1 = atmturb_measure_r0_pix(dcimg[ID1].array.F, size);
    double r0_2 = atmturb_measure_r0_pix(dcimg[ID2].array.F, size);
    printf("Master screens generated: r0_eff = %.3f px, %.3f px (target %.3f px)\n",
           r0_1, r0_2, spec.r0_pix);
    fflush(stdout);

    return 0;
}

/**
 * make_master_turbulence_screen - Generate von Karman master turbulence screens
 * @ID_name1: Output name for screen 1.
 * @ID_name2: Output name for screen 2.
 * @size: Grid dimension in pixels.
 * @outerscale: Outer scale in pixels.
 * @innerscale: Inner scale in pixels.
 * @WFprecision: Precision flag (0=single, 1=double).
 *
 * Return: 0 on success.
 */
int make_master_turbulence_screen(
    const char *ID_name1,
    const char *ID_name2,
    long        size,
    float       outerscale,
    float       innerscale,
    long        WFprecision)
{
    return make_master_turbulence_screen_seeded(ID_name1, ID_name2, size, outerscale,
                                                innerscale, WFprecision, 1);
}

/**
 * make_master_turbulence_screen_pow - Generate phase screen with arbitrary power-law PSD
 * @ID_name1: Output name for screen 1.
 * @ID_name2: Output name for screen 2.
 * @size: Grid dimension in pixels.
 * @power: Spatial frequency power-law exponent.
 *
 * Return: 0 on success.
 */
int make_master_turbulence_screen_pow(
    const char *ID_name1,
    const char *ID_name2,
    long        size,
    float       power)
{
    if (size <= 0)
    {
        return -1;
    }

    delete_image_ID(ID_name1);
    delete_image_ID(ID_name2);
    imageID ID1 = create_2Dimage_ID(ID_name1, size, size);
    imageID ID2 = create_2Dimage_ID(ID_name2, size, size);
    if (ID1 < 0 || ID2 < 0)
    {
        return -1;
    }

    long total = size * size;
    fftwf_complex *buf = (fftwf_complex *) fftwf_alloc_complex(total);
    if (!buf)
    {
        delete_image_ID(ID_name1);
        delete_image_ID(ID_name2);
        return -1;
    }

    uint64_t base_seed = atmturb_resolve_seed(0);
    double half_p = 0.5 * (double) power;

    #pragma omp parallel for schedule(static)
    for (long jj = 0; jj < size; jj++)
    {
        uint64_t rng = atmturb_rng_stream_seed(base_seed, (uint64_t) jj);
        long fy = (jj < size / 2) ? jj : (jj - size);
        double fy2 = (double) (fy * fy);
        long row = jj * size;

        for (long ii = 0; ii < size; ii++)
        {
            long fx = (ii < size / 2) ? ii : (ii - size);
            double r2 = (double) (fx * fx) + fy2;
            double g0, g1;
            atmturb_rng_gaussian_pair(&rng, &g0, &g1);

            if (fx == 0 && fy == 0)
            {
                buf[row + ii][0] = 0.0f;
                buf[row + ii][1] = 0.0f;
            }
            else
            {
                double amp = pow(r2, -half_p);
                buf[row + ii][0] = (float) (amp * g0);
                buf[row + ii][1] = (float) (amp * g1);
            }
        }
    }

    fftwf_plan plan = fftwf_plan_dft_2d((int) size, (int) size, buf, buf,
                                        FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    #pragma omp parallel for schedule(static)
    for (long i = 0; i < total; i++)
    {
        dcimg[ID1].array.F[i] = buf[i][0];
        dcimg[ID2].array.F[i] = buf[i][1];
    }
    fftwf_free(buf);

    return 0;
}
