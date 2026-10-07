// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_validate_math.c
 * @brief   Core math, I/O, and statistical primitives for milkatmturb validation
 */

#define _GNU_SOURCE
#include <fitsio.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "atmturb_simd.h"
#include "atmturb_validate_math.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

/**
 * @brief Load a 2D image or 3D cube from a FITS file.
 */
int val_cube_load(
    const char *filename,
    val_cube_t *cube)
{
    fitsfile *fptr = NULL;
    int status = 0;
    int naxis = 0;
    long naxes[3] = {1, 1, 1};

    fits_open_file(&fptr, filename, READONLY, &status);
    if (status != 0)
    {
        fprintf(stderr, "Error: cannot open FITS file %s (status %d)\n", filename, status);
        return -1;
    }

    fits_get_img_dim(fptr, &naxis, &status);
    fits_get_img_size(fptr, naxis, naxes, &status);
    if (status != 0 || naxis < 2)
    {
        fits_close_file(fptr, &status);
        return -1;
    }

    cube->naxis = naxis;
    cube->nx = naxes[0];
    cube->ny = naxes[1];
    cube->nframes = (naxis >= 3) ? naxes[2] : 1;
    cube->nelements = cube->nx * cube->ny * cube->nframes;
    cube->data = (double *) malloc(cube->nelements * sizeof(double));
    if (cube->data == NULL)
    {
        fits_close_file(fptr, &status);
        return -1;
    }

    long fpixel[3] = {1, 1, 1};
    fits_read_pix(fptr, TDOUBLE, fpixel, cube->nelements, NULL, cube->data, NULL, &status);
    fits_close_file(fptr, &status);
    if (status != 0)
    {
        free(cube->data);
        cube->data = NULL;
        return -1;
    }
    return 0;
}

/**
 * @brief Free resources associated with a val_cube_t.
 */
void val_cube_free(
    val_cube_t *cube)
{
    if (cube != NULL && cube->data != NULL)
    {
        free(cube->data);
        cube->data = NULL;
        cube->nelements = 0;
    }
}

/**
 * @brief Subtract mean per frame to remove global piston.
 */
void val_cube_subtract_piston(
    val_cube_t *cube)
{
    long frame_size = cube->nx * cube->ny;

    for (long f = 0; f < cube->nframes; f++)
    {
        double *frame = &cube->data[f * frame_size];
        double sum = 0.0;
        for (long i = 0; i < frame_size; i++)
        {
            sum += frame[i];
        }
        double mean = sum / (double) frame_size;
        for (long i = 0; i < frame_size; i++)
        {
            frame[i] -= mean;
        }
    }
}

/**
 * @brief Check if all elements are finite and calculate standard deviation.
 */
int val_cube_check_finite(
    const val_cube_t *cube,
    double           *out_std)
{
    int all_finite = 1;
    double sum = 0.0;
    double sum2 = 0.0;

    for (long i = 0; i < cube->nelements; i++)
    {
        double v = cube->data[i];
        if (!isfinite(v))
        {
            all_finite = 0;
            break;
        }
        sum += v;
        sum2 += v * v;
    }

    if (!all_finite || cube->nelements <= 1)
    {
        *out_std = 0.0;
        return 0;
    }

    double n = (double) cube->nelements;
    double var = (sum2 - (sum * sum) / n) / (n - 1.0);
    *out_std = (var > 0.0) ? sqrt(var) : 0.0;
    return 1;
}

/**
 * @brief Compute standard deviation across all pixels.
 */
double val_cube_std(
    const val_cube_t *cube)
{
    double std_val = 0.0;
    val_cube_check_finite(cube, &std_val);
    return std_val;
}

/**
 * @brief Accurate numerical quadrature of modified Bessel K_{5/6}(x).
 */
static double val_kv_5_6(
    double x)
{
    if (x <= 1e-12)
    {
        return 0.0;
    }
    double arg = (40.0 / x > 1.0) ? (40.0 / x) : 1.0;
    double t_max = acosh(arg);
    int nsteps = 400;
    double dt = t_max / (double) nsteps;
    double sum = 0.5 * (exp(-x) + exp(-x * cosh(t_max)) * cosh(5.0 / 6.0 * t_max));

    for (int i = 1; i < nsteps; i++)
    {
        double t = i * dt;
        sum += exp(-x * cosh(t)) * cosh(5.0 / 6.0 * t);
    }
    return sum * dt;
}

/**
 * @brief Continuous von Karman structure function [rad^2].
 */
double val_sf_vonkarman(
    double r,
    double r0,
    double L0)
{
    if (L0 <= 0.0)
    {
        return 6.88 * pow(r / r0, 5.0 / 3.0);
    }
    double x = 2.0 * M_PI * r / L0;
    if (x < 1e-4)
    {
        return 6.88 * pow(r / r0, 5.0 / 3.0);
    }
    double gamma_5_6 = tgamma(5.0 / 6.0);
    double factor = pow(2.0, 1.0 / 6.0) / gamma_5_6;
    double kv = val_kv_5_6(x);
    double bracket = 1.0 - factor * pow(x, 5.0 / 6.0) * kv;
    return 0.17253 * pow(L0 / r0, 5.0 / 3.0) * bracket;
}

/**
 * @brief Discrete PSD ensemble structure function via FFTW backward transform.
 */
void val_sf_discrete(
    int        nmaster,
    double     dx_m,
    double     r0,
    double     L0,
    double     l0,
    int        nlags,
    const int *lags,
    double    *dref)
{
    long total = (long) nmaster * nmaster;
    fftwf_complex *psd = (fftwf_complex *) fftwf_alloc_complex(total);
    fftwf_complex *cov = (fftwf_complex *) fftwf_alloc_complex(total);
    fftwf_plan plan = fftwf_plan_dft_2d(nmaster, nmaster, psd, cov, FFTW_BACKWARD, FFTW_ESTIMATE);

    double k0 = (L0 > 0.0) ? ((double) nmaster * dx_m / L0) : 0.0;
    double km = (l0 > 0.0) ? ((5.92 / (2.0 * M_PI)) * (double) nmaster * dx_m / l0) : 0.0;
    double pref = 0.023 * pow(r0 / dx_m, -5.0 / 3.0) * pow((double) nmaster, 5.0 / 3.0);

    for (int iy = 0; iy < nmaster; iy++)
    {
        double ky = (iy < nmaster / 2) ? (double) iy : (double) (iy - nmaster);
        for (int ix = 0; ix < nmaster; ix++)
        {
            double kx = (ix < nmaster / 2) ? (double) ix : (double) (ix - nmaster);
            long idx = (long) iy * nmaster + ix;
            if (ix == 0 && iy == 0)
            {
                psd[idx][0] = 0.0f;
            }
            else
            {
                double k2 = kx * kx + ky * ky;
                double val = pref * pow(k2 + k0 * k0, -11.0 / 6.0);
                if (km > 0.0)
                {
                    val *= exp(-k2 / (km * km));
                }
                psd[idx][0] = (float) val;
            }
            psd[idx][1] = 0.0f;
        }
    }

    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    double c00 = cov[0][0];
    for (int i = 0; i < nlags; i++)
    {
        int lag = lags[i] % nmaster;
        dref[i] = 2.0 * (c00 - (double) cov[lag][0]);
    }

    fftwf_free(psd);
    fftwf_free(cov);
}

/**
 * @brief Measure empirical structure function along x and y for given lags.
 */
void val_measure_sf(
    const val_cube_t *cube,
    int               nlags,
    const int        *lags,
    double           *dmeas)
{
    long nx = cube->nx;
    long ny = cube->ny;
    long frame_size = nx * ny;

    for (int li = 0; li < nlags; li++)
    {
        int lag = lags[li];
        double sum_dx2 = 0.0, sum_dy2 = 0.0;
        long count_x = 0, count_y = 0;

        for (long f = 0; f < cube->nframes; f++)
        {
            const double *fr = &cube->data[f * frame_size];

            for (long y = 0; y < ny; y++)
            {
                for (long x = 0; x < nx - lag; x++)
                {
                    double d = fr[y * nx + (x + lag)] - fr[y * nx + x];
                    sum_dx2 += d * d;
                    count_x++;
                }
            }

            for (long y = 0; y < ny - lag; y++)
            {
                for (long x = 0; x < nx; x++)
                {
                    double d = fr[(y + lag) * nx + x] - fr[y * nx + x];
                    sum_dy2 += d * d;
                    count_y++;
                }
            }
        }

        double mx = (count_x > 0) ? (sum_dx2 / (double) count_x) : 0.0;
        double my = (count_y > 0) ? (sum_dy2 / (double) count_y) : 0.0;
        dmeas[li] = 0.5 * (mx + my);
    }
}

/**
 * @brief Theoretical two-axis tip+tilt variance over circular pupil.
 */
double val_tilt_variance_theory(
    double D,
    double r0,
    double L0)
{
    if (L0 <= 0.0)
    {
        return 0.896 * pow(D / r0, 5.0 / 3.0);
    }
    int nsteps = 20000;
    double log_fmin = log(1e-6 / D);
    double log_fmax = log(1e3 / D);
    double dlog_f = (log_fmax - log_fmin) / (double) nsteps;
    double sum = 0.0;

    for (int i = 0; i < nsteps; i++)
    {
        double f = exp(log_fmin + (i + 0.5) * dlog_f);
        double df = f * dlog_f;
        double psd = 0.023 * pow(r0, -5.0 / 3.0) * pow(f * f + 1.0 / (L0 * L0), -11.0 / 6.0);
        double x = M_PI * D * f;
        double j2 = jn(2, x);
        double filt = 4.0 * pow(2.0 * j2 / x, 2.0);
        sum += psd * filt * 2.0 * M_PI * f * df;
    }
    return sum;
}

/**
 * @brief Project one frame onto normalized Zernike tip and tilt modes.
 */
static void project_frame_tilt(
    const double *frame,
    long          nx,
    long          ny,
    double        rad,
    double        cx,
    double        cy,
    double        sum_z2sq,
    double        sum_z3sq,
    double       *out_tilt_sq)
{
    double sum_v = 0.0;
    long count = 0;

    for (long y = 0; y < ny; y++)
    {
        double yn = (y - cy) / rad;
        for (long x = 0; x < nx; x++)
        {
            double xn = (x - cx) / rad;
            if (xn * xn + yn * yn <= 1.0)
            {
                sum_v += frame[y * nx + x];
                count++;
            }
        }
    }
    double mean_v = (count > 0) ? (sum_v / (double) count) : 0.0;

    double sum_vz2 = 0.0, sum_vz3 = 0.0;
    for (long y = 0; y < ny; y++)
    {
        double yn = (y - cy) / rad;
        for (long x = 0; x < nx; x++)
        {
            double xn = (x - cx) / rad;
            if (xn * xn + yn * yn <= 1.0)
            {
                double v = frame[y * nx + x] - mean_v;
                sum_vz2 += v * 2.0 * xn;
                sum_vz3 += v * 2.0 * yn;
            }
        }
    }

    double a2 = (sum_z2sq > 0.0) ? (sum_vz2 / sum_z2sq) : 0.0;
    double a3 = (sum_z3sq > 0.0) ? (sum_vz3 / sum_z3sq) : 0.0;
    *out_tilt_sq = a2 * a2 + a3 * a3;
}

/**
 * @brief Measure tip+tilt variance across frames over circular pupil aperture.
 */
double val_measure_tilt_variance(
    const val_cube_t *cube)
{
    long nx = cube->nx, ny = cube->ny;
    double rad = (nx < ny ? nx : ny) / 2.0;
    double cx = (nx - 1.0) / 2.0, cy = (ny - 1.0) / 2.0;
    long frame_size = nx * ny;
    double sum_z2sq = 0.0, sum_z3sq = 0.0;

    for (long y = 0; y < ny; y++)
    {
        double yn = (y - cy) / rad;
        for (long x = 0; x < nx; x++)
        {
            double xn = (x - cx) / rad;
            if (xn * xn + yn * yn <= 1.0)
            {
                double z2 = 2.0 * xn, z3 = 2.0 * yn;
                sum_z2sq += z2 * z2;
                sum_z3sq += z3 * z3;
            }
        }
    }

    double total_var = 0.0;
    for (long f = 0; f < cube->nframes; f++)
    {
        double tilt_sq = 0.0;
        project_frame_tilt(&cube->data[f * frame_size], nx, ny, rad, cx, cy,
                           sum_z2sq, sum_z3sq, &tilt_sq);
        total_var += tilt_sq;
    }

    return (cube->nframes > 0) ? (total_var / (double) cube->nframes) : 0.0;
}

/**
 * @brief Compute normalized cross-correlation between two frames at integer shift.
 */
static double val_xcorr_single_lag(
    const double *f0,
    const double *f1,
    long          nx,
    long          ny,
    long          ix,
    long          iy)
{
    long y0_a = (iy < 0) ? -iy : 0;
    long y1_a = (iy < 0) ? ny : (ny - iy);
    long x0_a = (ix < 0) ? -ix : 0;
    long x1_a = (ix < 0) ? nx : (nx - ix);
    long count = (y1_a - y0_a) * (x1_a - x0_a);
    if (count <= 0)
    {
        return 0.0;
    }

    double sum_a = 0.0, sum_b = 0.0;
    for (long y = y0_a; y < y1_a; y++)
    {
        long row_a = y * nx;
        long row_b = (y + iy) * nx + ix;
        for (long x = x0_a; x < x1_a; x++)
        {
            sum_a += f0[row_a + x];
            sum_b += f1[row_b + x];
        }
    }
    double mean_a = sum_a / (double) count;
    double mean_b = sum_b / (double) count;

    double sum_ab = 0.0, sum_aa = 0.0, sum_bb = 0.0;
    for (long y = y0_a; y < y1_a; y++)
    {
        long row_a = y * nx;
        long row_b = (y + iy) * nx + ix;
        for (long x = x0_a; x < x1_a; x++)
        {
            double da = f0[row_a + x] - mean_a;
            double db = f1[row_b + x] - mean_b;
            sum_ab += da * db;
            sum_aa += da * da;
            sum_bb += db * db;
        }
    }
    double denom = sqrt(sum_aa * sum_bb);
    return (denom > 0.0) ? (sum_ab / denom) : 0.0;
}

/**
 * @brief Sub-pixel shift between two 2D frames via normalized spatial cross-correlation.
 */
void val_xcorr_shift(
    const double *f0,
    const double *f1,
    long          nx,
    long          ny,
    double       *out_dx,
    double       *out_dy)
{
    long r = (nx / 4 < 16) ? (nx / 4) : 16;
    if (r < 4)
    {
        r = 4;
    }

    long grid_dim = 2 * r + 1;
    double *grid = (double *) calloc((size_t) (grid_dim * grid_dim), sizeof(double));
    if (grid == NULL)
    {
        *out_dx = 0.0;
        *out_dy = 0.0;
        return;
    }

    double best_c = -1e30;
    long best_dx = 0, best_dy = 0;

    for (long iy = -r; iy <= r; iy++)
    {
        for (long ix = -r; ix <= r; ix++)
        {
            double c = val_xcorr_single_lag(f0, f1, nx, ny, ix, iy);
            grid[(iy + r) * grid_dim + (ix + r)] = c;
            if (c > best_c)
            {
                best_c = c;
                best_dx = ix;
                best_dy = iy;
            }
        }
    }

    double dx_sub = (double) best_dx;
    if (best_dx > -r && best_dx < r)
    {
        double cm = grid[(best_dy + r) * grid_dim + (best_dx - 1 + r)];
        double c0 = grid[(best_dy + r) * grid_dim + (best_dx + r)];
        double cp = grid[(best_dy + r) * grid_dim + (best_dx + 1 + r)];
        double denom = cm - 2.0 * c0 + cp;
        if (denom != 0.0)
        {
            dx_sub += 0.5 * (cm - cp) / denom;
        }
    }

    double dy_sub = (double) best_dy;
    if (best_dy > -r && best_dy < r)
    {
        double cm = grid[(best_dy - 1 + r) * grid_dim + (best_dx + r)];
        double c0 = grid[(best_dy + r) * grid_dim + (best_dx + r)];
        double cp = grid[(best_dy + 1 + r) * grid_dim + (best_dx + r)];
        double denom = cm - 2.0 * c0 + cp;
        if (denom != 0.0)
        {
            dy_sub += 0.5 * (cm - cp) / denom;
        }
    }

    free(grid);
    *out_dx = dx_sub;
    *out_dy = dy_sub;
}

/**
 * @brief Pearson correlation between two arrays.
 */
static double pearson_correlation(
    const double *a,
    const double *b,
    long          n)
{
    double sum_a = 0.0, sum_b = 0.0;
    for (long i = 0; i < n; i++)
    {
        sum_a += a[i];
        sum_b += b[i];
    }
    double ma = sum_a / (double) n, mb = sum_b / (double) n;

    double s_ab = 0.0, s_aa = 0.0, s_bb = 0.0;
    for (long i = 0; i < n; i++)
    {
        double da = a[i] - ma;
        double db = b[i] - mb;
        s_ab += da * db;
        s_aa += da * da;
        s_bb += db * db;
    }

    double denom = sqrt(s_aa * s_bb);
    return (denom > 0.0) ? (s_ab / denom) : 0.0;
}

/**
 * @brief Pearson correlation between two planes in a 3D cube.
 */
double val_cube_plane_corr(
    const val_cube_t *cube,
    int               plane_a,
    int               plane_b)
{
    long frame_size = cube->nx * cube->ny;
    const double *a = &cube->data[plane_a * frame_size];
    const double *b = &cube->data[plane_b * frame_size];
    return pearson_correlation(a, b, frame_size);
}

/**
 * @brief Worst-case Pearson correlation between frames separated by fixed lag.
 */
double val_lag_repeat_corr(
    const val_cube_t *cube,
    int               lag)
{
    long frame_size = cube->nx * cube->ny;
    long n = cube->nframes - lag;
    if (n <= 0)
    {
        return 0.0;
    }

    double worst_corr = 0.0;
    for (long f = 0; f < n; f++)
    {
        const double *f0 = &cube->data[f * frame_size];
        const double *f1 = &cube->data[(f + lag) * frame_size];
        double corr = fabs(pearson_correlation(f0, f1, frame_size));
        if (corr > worst_corr)
        {
            worst_corr = corr;
        }
    }
    return worst_corr;
}

/**
 * @brief Compute mean intensity and scintillation index of amplitude cube.
 */
void val_scintillation(
    const val_cube_t *cube,
    double           *out_mean,
    double           *out_scint_index)
{
    double sum_i = 0.0, sum_i2 = 0.0;

    for (long k = 0; k < cube->nelements; k++)
    {
        double a = cube->data[k];
        double intensity = a * a;
        sum_i += intensity;
        sum_i2 += intensity * intensity;
    }

    double n = (double) cube->nelements;
    double mean_i = (n > 0.0) ? (sum_i / n) : 0.0;
    double var_i = (n > 1.0) ? ((sum_i2 - (sum_i * sum_i) / n) / (n - 1.0)) : 0.0;

    *out_mean = mean_i;
    *out_scint_index = (mean_i > 0.0) ? (var_i / (mean_i * mean_i)) : 0.0;
}

/**
 * @brief Measure peak-to-peak high-frequency power variation across frames (breathing metric).
 */
double val_cube_breathing_ratio(
    const val_cube_t *cube)
{
    if (cube == NULL || cube->nframes < 2 || cube->nx < 2 || cube->ny < 2)
    {
        return 0.0;
    }

    long nx = cube->nx;
    long ny = cube->ny;
    long nframes = cube->nframes;
    long frame_pixels = nx * ny;

    double p_min = 1e30;
    double p_max = -1e30;
    double p_sum = 0.0;

    for (long t = 0; t < nframes; t++)
    {
        const double *f = &cube->data[t * frame_pixels];
        double diff_sum = 0.0;
        long count = 0;

        for (long y = 0; y < ny; y++)
        {
            long row = y * nx;
            for (long x = 0; x < nx - 1; x++)
            {
                double dx = f[row + x + 1] - f[row + x];
                diff_sum += dx * dx;
                count++;
            }
        }

        for (long y = 0; y < ny - 1; y++)
        {
            long row = y * nx;
            long row_next = (y + 1) * nx;
            for (long x = 0; x < nx; x++)
            {
                double dy = f[row_next + x] - f[row + x];
                diff_sum += dy * dy;
                count++;
            }
        }

        double pt = (count > 0) ? (diff_sum / (double) count) : 0.0;
        if (pt < p_min)
        {
            p_min = pt;
        }
        if (pt > p_max)
        {
            p_max = pt;
        }
        p_sum += pt;
    }

    double p_mean = p_sum / (double) nframes;
    if (p_mean <= 1e-30)
    {
        return 0.0;
    }
    return (p_max - p_min) / p_mean;
}

/**
 * @brief Compare scalar and vectorized SIMD extrusions across schemes and strides.
 */
int val_check_simd_parity(
    double  tol,
    char   *msg_buf,
    size_t  msg_size)
{
    long msize = 512, pup_size = 128;
    long mtot = msize * msize, ptot = pup_size * pup_size;
    float *master = (float *) malloc(sizeof(float) * mtot);
    float *out_s  = (float *) calloc((size_t) ptot, sizeof(float));
    float *out_v  = (float *) calloc((size_t) ptot, sizeof(float));

    if (!master || !out_s || !out_v)
    {
        free(master); free(out_s); free(out_v);
        snprintf(msg_buf, msg_size, "allocation failure in SIMD parity test");
        return 0;
    }

    for (long i = 0; i < mtot; i++)
    {
        master[i] = (float) sin((double) i * 0.1) + (float) cos((double) i * 0.03);
    }

    double max_diff = 0.0;
    for (int interp = 0; interp <= 1; interp++)
    {
        for (int os = 1; os <= 2; os++)
        {
            for (double sub = 0.0; sub < 1.0; sub += 0.25)
            {
                atmturb_extrude_params_t ep_s = {
                    .master   = master,
                    .msize    = msize,
                    .x0       = 15.3 + sub,
                    .y0       = 22.7 + sub,
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

                for (long i = 0; i < ptot; i++)
                {
                    double diff = fabs((double) out_s[i] - (double) out_v[i]);
                    if (diff > max_diff)
                    {
                        max_diff = diff;
                    }
                }
            }
        }
    }

    free(master); free(out_s); free(out_v);
    snprintf(msg_buf, msg_size, "SIMD parity max diff = %.3e (tol %.3e)", max_diff, tol);
    return (max_diff <= tol);
}
