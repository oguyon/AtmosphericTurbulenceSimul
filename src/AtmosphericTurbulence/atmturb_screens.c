// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_screens.c
 * @brief   Master atmospheric turbulence phase screen generation
 */

#define _GNU_SOURCE
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
#include "atmturb_fft_wisdom.h"
#include "atmturb_fresnel_simd.h"
#include "atmturb_screens_batch.h"
#include "atmturb_types.h"

/**
 * atmturb_rng_draw_strip8 - Counter-unrolled 8-lane SplitMix64 uniform pair draw
 */
static inline void atmturb_rng_draw_strip8(
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
 * atmturb_fill_spectrum_row_float - Populate one row of 2D single-precision von Karman spectrum
 * @row_buf: Destination row pointer (`size` complex elements).
 * @size: Grid linear dimension in pixels.
 * @jj: Row index in `[0, size)`.
 * @base_seed: Resolved PRNG base seed for the screen pair.
 * @base_p: Precomputed spectral constants `[k0_sq, sqrt_pref, inv_two_km2]`.
 */
void atmturb_fill_spectrum_row_float(
    fftwf_complex *row_buf,
    long           size,
    long           jj,
    uint64_t       base_seed,
    const float   *base_p)
{
    fftwf_complex *restrict rbuf = row_buf;
    const long              half = size / 2;
    uint64_t                rng  = atmturb_rng_stream_seed(base_seed, (uint64_t) jj);
    long                    fy   = (jj < half) ? jj : (jj - size);
    float                   fy2  = (float) (fy * fy);
    float                   spec_p[4] = {fy2, base_p[0], base_p[1], base_p[2]};

    for (long ii = 0; ii < size; ii++)
    {
        int in_half = (ii + 8 <= half) || (ii >= half && ii + 8 <= size);
        int far_dc  = (fy2 >= 256.0f) || (ii >= 16 && ii + 8 <= size - 16);
        if (in_half && far_dc)
        {
            float u1[8];
            float u2[8];
            atmturb_rng_draw_strip8(&rng, u1, u2);
            float fx0 = (float) ((ii < half) ? ii : (ii - size));
            atmturb_fresnel_simd_spectrum_strip(&rbuf[ii], u1, u2, fx0, spec_p);
            ii += 7;
            continue;
        }

        long   fx = (ii < half) ? ii : (ii - size);
        float  r2 = (float) (fx * fx) + fy2;
        double g0, g1;
        atmturb_rng_gaussian_pair(&rng, &g0, &g1);
        if (fx == 0 && fy == 0)
        {
            rbuf[ii][0] = 0.0f;
            rbuf[ii][1] = 0.0f;
        }
        else
        {
            float amp = spec_p[2] * powf(r2 + spec_p[1], -11.0f / 12.0f);
            if (spec_p[3] > 0.0f)
            {
                amp *= expf(-r2 * spec_p[3]);
            }
            rbuf[ii][0] = amp * (float) g0;
            rbuf[ii][1] = amp * (float) g1;
        }
    }
}

/**
 * atmturb_fill_spectrum_float_rows - Fill 2D Fourier spectrum with shaped Gaussian noise (float)
 * @buf: Pointer to FFTW single-precision complex frequency buffer.
 * @spec: Pointer to screen generation specifications.
 * @use_omp: 1 to parallelize rows via OpenMP, 0 inside an outer worker team.
 */
void atmturb_fill_spectrum_float_rows(
    fftwf_complex               *buf,
    const atmturb_screen_spec_t *spec,
    int                          use_omp)
{
    long     size      = spec->size;
    float    r0_pix    = (spec->r0_pix > 0.0) ? (float) spec->r0_pix : powf(6.88f, 0.6f);
    float    k0        = (spec->L0_pix > 0.0) ? ((float) size / (float) spec->L0_pix) : 0.0f;
    float    km_scale  = 5.92f / (2.0f * (float) M_PI);
    float    km        = (spec->l0_pix > 0.0)
                             ? (km_scale * (float) size / (float) spec->l0_pix)
                             : 0.0f;
    float    pref      = 0.023f * powf(r0_pix, -5.0f / 3.0f) * powf((float) size, 5.0f / 3.0f);
    float    base_p[3] = {k0 * k0, sqrtf(pref), (km > 0.0f) ? (0.5f / (km * km)) : 0.0f};
    uint64_t base_seed = atmturb_resolve_seed(spec->seed);

#ifdef _OPENMP
    #pragma omp parallel for if (use_omp) schedule(static)
#else
    (void) use_omp;
#endif
    for (long jj = 0; jj < size; jj++)
    {
        atmturb_fill_spectrum_row_float(&buf[jj * size], size, jj, base_seed, base_p);
    }
}

/**
 * atmturb_fill_spectrum_row_double - Populate one row of 2D double-precision von Karman spectrum
 * @row_buf: Destination row pointer (`size` double complex elements).
 * @size: Grid linear dimension in pixels.
 * @jj: Row index in `[0, size)`.
 * @base_seed: Resolved PRNG base seed for the screen pair.
 * @base_p: Precomputed spectral constants `[k0_sq, sqrt_pref, inv_two_km2]`.
 */
void atmturb_fill_spectrum_row_double(
    fftw_complex *row_buf,
    long          size,
    long          jj,
    uint64_t      base_seed,
    const double *base_p)
{
    fftw_complex *restrict rbuf = row_buf;
    const long             half = size / 2;
    uint64_t               rng  = atmturb_rng_stream_seed(base_seed, (uint64_t) jj);
    long                   fy   = (jj < half) ? jj : (jj - size);
    double                 fy2  = (double) (fy * fy);
    float                  spec_pf[4] = {
        (float) fy2, (float) base_p[0], (float) base_p[1], (float) base_p[2]
    };

    for (long ii = 0; ii < size; ii++)
    {
        int in_half = (ii + 8 <= half) || (ii >= half && ii + 8 <= size);
        int far_dc  = (fy2 >= 1024.0) || (ii >= 32 && ii + 8 <= size - 32);
        if (in_half && far_dc)
        {
            float         u1[8];
            float         u2[8];
            fftwf_complex c8[8];
            atmturb_rng_draw_strip8(&rng, u1, u2);
            float fx0 = (float) ((ii < half) ? ii : (ii - size));
            atmturb_fresnel_simd_spectrum_strip(c8, u1, u2, fx0, spec_pf);
            for (int lane = 0; lane < 8; lane++)
            {
                rbuf[ii + lane][0] = (double) c8[lane][0];
                rbuf[ii + lane][1] = (double) c8[lane][1];
            }
            ii += 7;
            continue;
        }

        long   fx = (ii < half) ? ii : (ii - size);
        double r2 = (double) (fx * fx) + fy2;
        double g0, g1;
        atmturb_rng_gaussian_pair(&rng, &g0, &g1);
        if (fx == 0 && fy == 0)
        {
            rbuf[ii][0] = 0.0;
            rbuf[ii][1] = 0.0;
        }
        else
        {
            double l_exp = -11.0 / 12.0 * log(r2 + base_p[0]) - r2 * base_p[2];
            double amp   = base_p[1] * exp(l_exp);
            rbuf[ii][0]  = amp * g0;
            rbuf[ii][1]  = amp * g1;
        }
    }
}

/**
 * atmturb_fftwf_plan_threads - Configure FFTW single-precision planner thread count
 * @req_nthreads: Requested thread count (>0 explicit, <=0 uses omp_get_max_threads()).
 */
void atmturb_fftwf_plan_threads(
    int req_nthreads)
{
#ifdef _OPENMP
    static int s_threads_init = 0;
    if (!s_threads_init)
    {
        #pragma omp critical
        {
            if (!s_threads_init && fftwf_init_threads())
            {
                s_threads_init = 1;
            }
        }
    }
    if (s_threads_init)
    {
        int nth = (req_nthreads > 0) ? req_nthreads : omp_get_max_threads();
        fftwf_plan_with_nthreads((nth < 1) ? 1 : nth);
    }
#else
    (void) req_nthreads;
#endif
}

/**
 * atmturb_fftw_plan_threads - Configure FFTW double-precision planner thread count
 * @req_nthreads: Requested thread count (>0 explicit, <=0 uses omp_get_max_threads()).
 */
void atmturb_fftw_plan_threads(
    int req_nthreads)
{
#ifdef _OPENMP
    static int s_threads_init = 0;
    if (!s_threads_init)
    {
        #pragma omp critical
        {
            if (!s_threads_init && fftw_init_threads())
            {
                s_threads_init = 1;
            }
        }
    }
    if (s_threads_init)
    {
        int nth = (req_nthreads > 0) ? req_nthreads : omp_get_max_threads();
        fftw_plan_with_nthreads((nth < 1) ? 1 : nth);
    }
#else
    (void) req_nthreads;
#endif
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

    float *oa[1] = {screen_a};
    float *ob[1] = {screen_b};
    atmturb_screen_batch_req_t req;
    memset(&req, 0, sizeof(req));
    req.specs  = spec;
    req.out_a  = oa;
    req.out_b  = ob;
    req.npairs = 1;
    return atmturb_generate_screen_pairs_batch(&req);
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
    long   count  = 2 * size * size;

    #pragma omp parallel for reduction(+:sum_sq) schedule(static)
    for (long y = 0; y < size; y++)
    {
        long         y_next = (y + 1 < size) ? (y + 1) : 0;
        const float *r_cur  = &data[y * size];
        const float *r_nxt  = &data[y_next * size];
        double       row_sq = 0.0;

        #pragma omp simd reduction(+:row_sq)
        for (long x = 0; x < size - 1; x++)
        {
            double cur = (double) r_cur[x];
            double dx  = (double) r_cur[x + 1] - cur;
            double dy  = (double) r_nxt[x] - cur;
            row_sq += dx * dx + dy * dy;
        }
        double last_cur = (double) r_cur[size - 1];
        double last_dx  = (double) r_cur[0] - last_cur;
        double last_dy  = (double) r_nxt[size - 1] - last_cur;
        sum_sq += row_sq + last_dx * last_dx + last_dy * last_dy;
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

    atmturb_fftwf_wisdom_load();
    atmturb_fftwf_plan_threads(0);
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
