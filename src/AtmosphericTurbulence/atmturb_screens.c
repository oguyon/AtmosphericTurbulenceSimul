// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_screens.c
 * @brief   Master atmospheric turbulence phase screen generation
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <fftw3.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "AtmosphericTurbulence.h"
#include "atmturb_types.h"

/**
 * struct atmturb_screen_params_t - Parameter container for phase screen generation
 * @size: Grid linear dimension in pixels.
 * @outer_f0: Spatial frequency cutoff corresponding to outer scale.
 * @inner_f0: Spatial frequency cutoff corresponding to inner scale.
 * @rlim: Inner frequency cutoff radius.
 * @rlim_mode: 1 to zero out frequencies below rlim.
 * @psd_exponent: Amplitude spectrum roll-off exponent.
 * @power_exponent: Structure function scaling exponent.
 */
typedef struct {
    long size;
    double outer_f0;
    double inner_f0;
    double rlim;
    int rlim_mode;
    double psd_exponent;
    double power_exponent;
} atmturb_screen_params_t;

/**
 * atmturb_next_rng_seed - Generate unique base seed for parallel PRNG streams
 *
 * Return: Unique 64-bit seed value.
 */
static uint64_t atmturb_next_rng_seed(void)
{
    static uint64_t g_seed = 0x853c49e6748fea9bULL;
    uint64_t s;
#ifdef _OPENMP
    #pragma omp atomic capture
    {
        s = g_seed;
        g_seed += 0x9e3779b97f4a7c15ULL;
    }
#else
    s = g_seed;
    g_seed += 0x9e3779b97f4a7c15ULL;
#endif
    return s;
}

/**
 * atmturb_fill_spectrum_grid_float - Fill 2D Fourier spectrum with shaped Gaussian noise (float)
 * @buf: Pointer to FFTW complex frequency buffer.
 * @p: Pointer to screen generation parameters.
 */
static void atmturb_fill_spectrum_grid_float(fftwf_complex *buf,
                                            const atmturb_screen_params_t *p)
{
    long size = p->size;
    double f0_sq = p->outer_f0 * p->outer_f0;
    double inv_two_inner2 = (p->inner_f0 > 0.0) ? (0.5 / (p->inner_f0 * p->inner_f0)) : 0.0;
    double rlim2 = p->rlim * p->rlim;
    uint64_t base_seed = atmturb_next_rng_seed();

    #pragma omp parallel
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        uint64_t rng = base_seed + (uint64_t)tid * 0x517cc1b727220a95ULL;
        atmturb_rng_splitmix64(&rng);

        #pragma omp for
        for (long jj = 0; jj < size; jj++)
        {
            long fy = (jj < size / 2) ? jj : (jj - size);
            double fy2 = (double)(fy * fy);
            long row = jj * size;

            for (long ii = 0; ii < size; ii++)
            {
                long fx = (ii < size / 2) ? ii : (ii - size);
                double r2 = (double)(fx * fx) + fy2;
                if (r2 == 0.0 || (p->rlim_mode && r2 < rlim2))
                {
                    buf[row + ii][0] = 0.0f;
                    buf[row + ii][1] = 0.0f;
                }
                else
                {
                    double dist = sqrt(r2 + f0_sq);
                    double factor = (inv_two_inner2 > 0.0) ? exp(-r2 * inv_two_inner2) : 1.0;
                    double amp = factor / pow(dist, p->psd_exponent);
                    double g0, g1;
                    atmturb_rng_gaussian_pair(&rng, &g0, &g1);
                    buf[row + ii][0] = (float)(amp * g0);
                    buf[row + ii][1] = (float)(amp * g1);
                }
            }
        }
    }
}

/**
 * atmturb_fill_spectrum_grid_double - Fill 2D Fourier spectrum with shaped Gaussian noise (double)
 * @buf: Pointer to FFTW double complex frequency buffer.
 * @p: Pointer to screen generation parameters.
 */
static void atmturb_fill_spectrum_grid_double(fftw_complex *buf,
                                             const atmturb_screen_params_t *p)
{
    long size = p->size;
    double f0_sq = p->outer_f0 * p->outer_f0;
    double inv_two_inner2 = (p->inner_f0 > 0.0) ? (0.5 / (p->inner_f0 * p->inner_f0)) : 0.0;
    double rlim2 = p->rlim * p->rlim;
    uint64_t base_seed = atmturb_next_rng_seed();

    #pragma omp parallel
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        uint64_t rng = base_seed + (uint64_t)tid * 0x517cc1b727220a95ULL;
        atmturb_rng_splitmix64(&rng);

        #pragma omp for
        for (long jj = 0; jj < size; jj++)
        {
            long fy = (jj < size / 2) ? jj : (jj - size);
            double fy2 = (double)(fy * fy);
            long row = jj * size;

            for (long ii = 0; ii < size; ii++)
            {
                long fx = (ii < size / 2) ? ii : (ii - size);
                double r2 = (double)(fx * fx) + fy2;
                if (r2 == 0.0 || (p->rlim_mode && r2 < rlim2))
                {
                    buf[row + ii][0] = 0.0;
                    buf[row + ii][1] = 0.0;
                }
                else
                {
                    double dist = sqrt(r2 + f0_sq);
                    double factor = (inv_two_inner2 > 0.0) ? exp(-r2 * inv_two_inner2) : 1.0;
                    double amp = factor / pow(dist, p->psd_exponent);
                    double g0, g1;
                    atmturb_rng_gaussian_pair(&rng, &g0, &g1);
                    buf[row + ii][0] = amp * g0;
                    buf[row + ii][1] = amp * g1;
                }
            }
        }
    }
}

/**
 * atmturb_calc_structure_diff_float - Mean squared difference for displacement (dx, dy)
 * @data: Pointer to linear 2D phase screen buffer.
 * @size: Grid dimension in pixels.
 * @dx: Horizontal displacement in pixels.
 * @dy: Vertical displacement in pixels.
 *
 * Return: Mean squared difference D.
 */
static double atmturb_calc_structure_diff_float(const float *data, long size,
                                               long dx, long dy)
{
    double sum_sq = 0.0;
    #pragma omp parallel for reduction(+:sum_sq)
    for (long y = 0; y < size; y++)
    {
        long y2 = (y + dy) % size;
        long row1 = y * size;
        long row2 = y2 * size;
        for (long x = 0; x < size; x++)
        {
            long x2 = (x + dx) % size;
            double diff = (double)data[row2 + x2] - (double)data[row1 + x];
            sum_sq += diff * diff;
        }
    }
    return sum_sq / (double)(size * size);
}

/**
 * atmturb_calc_structure_diff_double - Mean squared difference for displacement (dx, dy)
 * @data: Pointer to linear 2D double phase screen buffer.
 * @size: Grid dimension in pixels.
 * @dx: Horizontal displacement in pixels.
 * @dy: Vertical displacement in pixels.
 *
 * Return: Mean squared difference D.
 */
static double atmturb_calc_structure_diff_double(const double *data, long size,
                                                long dx, long dy)
{
    double sum_sq = 0.0;
    #pragma omp parallel for reduction(+:sum_sq)
    for (long y = 0; y < size; y++)
    {
        long y2 = (y + dy) % size;
        long row1 = y * size;
        long row2 = y2 * size;
        for (long x = 0; x < size; x++)
        {
            long x2 = (x + dx) % size;
            double diff = data[row2 + x2] - data[row1 + x];
            sum_sq += diff * diff;
        }
    }
    return sum_sq / (double)(size * size);
}

/**
 * atmturb_measure_structure_constant_float - Compute structure constant for float array
 * @data: Pointer to linear 2D phase screen buffer.
 * @size: Grid dimension in pixels.
 * @power_exponent: Power law scaling exponent.
 *
 * Return: Measured structure function constant C.
 */
static double atmturb_measure_structure_constant_float(const float *data, long size,
                                                      double power_exponent)
{
    double value = 0.0;
    long cnt = 0;
    long Dlim = 3;

    for (long ii = 1; ii < Dlim; ii++)
    {
        for (long jj = 1; jj < Dlim; jj++)
        {
            double D = atmturb_calc_structure_diff_float(data, size, ii, jj);
            if (D < 1e-30)
            {
                D = 1e-30;
            }
            value += log10(D) - power_exponent * log10(sqrt((double)(ii * ii + jj * jj)));
            cnt++;
        }
    }
    return (cnt > 0) ? pow(10.0, value / cnt) : 1.0;
}

/**
 * atmturb_measure_structure_constant_double - Compute structure constant for double array
 * @data: Pointer to linear 2D double phase screen buffer.
 * @size: Grid dimension in pixels.
 * @power_exponent: Power law scaling exponent.
 *
 * Return: Measured structure function constant C.
 */
static double atmturb_measure_structure_constant_double(const double *data, long size,
                                                       double power_exponent)
{
    double value = 0.0;
    long cnt = 0;
    long Dlim = 3;

    for (long ii = 1; ii < Dlim; ii++)
    {
        for (long jj = 1; jj < Dlim; jj++)
        {
            double D = atmturb_calc_structure_diff_double(data, size, ii, jj);
            if (D < 1e-30)
            {
                D = 1e-30;
            }
            value += log10(D) - power_exponent * log10(sqrt((double)(ii * ii + jj * jj)));
            cnt++;
        }
    }
    return (cnt > 0) ? pow(10.0, value / cnt) : 1.0;
}

/**
 * atmturb_synthesize_screens_float - Generate and normalize a pair of single-precision screens
 * @ID_name1: Name of first screen image.
 * @ID_name2: Name of second screen image.
 * @p: Pointer to screen generation parameters.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
static int atmturb_synthesize_screens_float(
    const char                    *ID_name1,
    const char                    *ID_name2,
    const atmturb_screen_params_t *p)
{
    long size = p->size;
    long ntot = size * size;
    fftwf_complex *buf = (fftwf_complex *)fftwf_alloc_complex(ntot);
    if (!buf)
    {
        return -1;
    }

    atmturb_fill_spectrum_grid_float(buf, p);

    fftwf_plan plan = fftwf_plan_dft_2d((int)size, (int)size, buf, buf,
                                        FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    delete_image_ID(ID_name1);
    delete_image_ID(ID_name2);
    imageID ID1 = create_2Dimage_ID(ID_name1, size, size);
    imageID ID2 = create_2Dimage_ID(ID_name2, size, size);

    #pragma omp parallel for
    for (long i = 0; i < ntot; i++)
    {
        dcimg[ID1].array.F[i] = buf[i][0];
        dcimg[ID2].array.F[i] = buf[i][1];
    }
    fftwf_free(buf);

    double C1 = atmturb_measure_structure_constant_float(dcimg[ID1].array.F,
                                                         size, p->power_exponent);
    double C2 = atmturb_measure_structure_constant_float(dcimg[ID2].array.F,
                                                         size, p->power_exponent);
    printf("C1, C2 =   %f %f\n", C1, C2);
    fflush(stdout);

    double scale1 = (C1 > 0.0 && !isnan(C1)) ? (1.0 / sqrt(C1)) : 1.0;
    double scale2 = (C2 > 0.0 && !isnan(C2)) ? (1.0 / sqrt(C2)) : 1.0;

    #pragma omp parallel for
    for (long i = 0; i < ntot; i++)
    {
        dcimg[ID1].array.F[i] *= (float)scale1;
        dcimg[ID2].array.F[i] *= (float)scale2;
    }

    return 0;
}

/**
 * atmturb_synthesize_screens_double - Generate and normalize a pair of double-precision screens
 * @ID_name1: Name of first screen image.
 * @ID_name2: Name of second screen image.
 * @p: Pointer to screen generation parameters.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
static int atmturb_synthesize_screens_double(
    const char                    *ID_name1,
    const char                    *ID_name2,
    const atmturb_screen_params_t *p)
{
    long size = p->size;
    long ntot = size * size;
    fftw_complex *buf = (fftw_complex *)fftw_alloc_complex(ntot);
    if (!buf)
    {
        return -1;
    }

    atmturb_fill_spectrum_grid_double(buf, p);

    fftw_plan plan = fftw_plan_dft_2d((int)size, (int)size, buf, buf,
                                      FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    delete_image_ID(ID_name1);
    delete_image_ID(ID_name2);
    imageID ID1 = create_2Dimage_ID_double(ID_name1, size, size);
    imageID ID2 = create_2Dimage_ID_double(ID_name2, size, size);

    #pragma omp parallel for
    for (long i = 0; i < ntot; i++)
    {
        dcimg[ID1].array.D[i] = buf[i][0];
        dcimg[ID2].array.D[i] = buf[i][1];
    }
    fftw_free(buf);

    double C1 = atmturb_measure_structure_constant_double(dcimg[ID1].array.D,
                                                          size, p->power_exponent);
    double C2 = atmturb_measure_structure_constant_double(dcimg[ID2].array.D,
                                                          size, p->power_exponent);
    printf("C1, C2 =   %f %f\n", C1, C2);
    fflush(stdout);

    double scale1 = (C1 > 0.0 && !isnan(C1)) ? (1.0 / sqrt(C1)) : 1.0;
    double scale2 = (C2 > 0.0 && !isnan(C2)) ? (1.0 / sqrt(C2)) : 1.0;

    #pragma omp parallel for
    for (long i = 0; i < ntot; i++)
    {
        dcimg[ID1].array.D[i] *= scale1;
        dcimg[ID2].array.D[i] *= scale2;
    }

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
    printf("Make turbulence screen, precision = %ld\n", WFprecision);
    fflush(stdout);

    atmturb_screen_params_t p;
    memset(&p, 0, sizeof(p));
    p.size = size;
    p.psd_exponent = 11.0 / 6.0;
    p.power_exponent = 5.0 / 3.0;

    imageID IDv = variable_ID("RLIM");
    if (IDv != -1)
    {
        p.rlim_mode = 1;
        p.rlim = dcvar[IDv].value.f;
        printf("R limit = %f pix\n", p.rlim);
    }

    p.outer_f0 = (outerscale > 0.0f) ? (1.0 * size / outerscale) : 0.0;
    p.inner_f0 = (innerscale > 0.0f) ? ((5.92 / (2.0 * M_PI)) * size / innerscale) : 0.0;

    if (WFprecision == 0)
    {
        return atmturb_synthesize_screens_float(ID_name1, ID_name2, &p);
    }
    return atmturb_synthesize_screens_double(ID_name1, ID_name2, &p);
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
    atmturb_screen_params_t p;
    memset(&p, 0, sizeof(p));
    p.size = size;
    p.psd_exponent = power;
    p.power_exponent = power;
    p.outer_f0 = 0.0;
    p.inner_f0 = 0.0;

    return atmturb_synthesize_screens_float(ID_name1, ID_name2, &p);
}
