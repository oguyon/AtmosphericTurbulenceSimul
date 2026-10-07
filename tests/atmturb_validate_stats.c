// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_validate_stats.c
 * @brief   CLI test runner for statistical turbulence validation in milkatmturb
 */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "atmturb_validate_math.h"

#define ARCSEC_RAD (3.14159265358979323846 / 180.0 / 3600.0)

/**
 * @brief Print formatted test verdict and return 0 (pass) or 1 (fail).
 */
static int verdict(
    int         ok,
    const char *msg)
{
    printf("%s: %s\n", ok ? "PASS" : "FAIL", msg);
    return ok ? 0 : 1;
}

/**
 * @brief Resolve r0 [m] taking airmass into account.
 */
static double resolve_r0(
    double r0_in,
    double seeing,
    double lam,
    double zenith)
{
    double r0 = (r0_in > 0.0) ? r0_in : (0.98 * lam / (seeing * ARCSEC_RAD));
    if (zenith != 0.0)
    {
        r0 *= pow(cos(zenith), 0.6);
    }
    return r0;
}

/**
 * @brief Subcommand: finite - verify all pixels are finite and non-flat.
 */
static int cmd_finite(
    int    argc,
    char **argv)
{
    if (argc < 2) return 1;
    const char *fname = argv[1];
    double min_std = 1e-3;

    for (int i = 2; i < argc; i++)
    {
        if (strcmp(argv[i], "--min-std") == 0 && i + 1 < argc) min_std = atof(argv[++i]);
    }

    val_cube_t cube;
    if (val_cube_load(fname, &cube) != 0)
    {
        return verdict(0, "failed to load FITS cube");
    }

    double std_val = 0.0;
    int is_fin = val_cube_check_finite(&cube, &std_val);
    int ok = is_fin && (std_val > min_std);

    char msg[256];
    snprintf(msg, sizeof(msg), "%s: finite=%d std=%.4g (min %.4g)",
             fname, is_fin, std_val, min_std);
    val_cube_free(&cube);
    return verdict(ok, msg);
}

/**
 * @brief Subcommand: same - verify two cubes match within relative tolerance.
 */
static int cmd_same(
    int    argc,
    char **argv)
{
    if (argc < 3)
    {
        fprintf(stderr, "Usage: atmturb-validate-stats same <fileA> <fileB> [--tol <val>]\n");
        return 1;
    }
    const char *fname_a = argv[1];
    const char *fname_b = argv[2];
    double tol = 1e-3;

    for (int i = 3; i < argc; i++)
    {
        if (strcmp(argv[i], "--tol") == 0 && i + 1 < argc)
        {
            tol = atof(argv[++i]);
        }
    }

    val_cube_t a, b;
    if (val_cube_load(fname_a, &a) != 0 || val_cube_load(fname_b, &b) != 0)
    {
        return verdict(0, "failed to load FITS cubes for comparison");
    }

    if (a.nx != b.nx || a.ny != b.ny || a.nframes != b.nframes)
    {
        val_cube_free(&a);
        val_cube_free(&b);
        return verdict(0, "shape mismatch between cubes");
    }

    double max_diff = 0.0;
    for (long i = 0; i < a.nelements; i++)
    {
        double diff = fabs(a.data[i] - b.data[i]);
        if (diff > max_diff)
        {
            max_diff = diff;
        }
    }

    double std_a = val_cube_std(&a);
    double scale = (std_a > 1e-30) ? std_a : 1e-30;
    double rel = max_diff / scale;

    char msg[256];
    snprintf(msg, sizeof(msg), "max |a-b| / std(a) = %.3g (tol %.3g)", rel, tol);
    val_cube_free(&a);
    val_cube_free(&b);
    return verdict(rel <= tol, msg);
}

/**
 * @brief Subcommand: ratio - compute RMS or variance ratio between two cubes.
 */
static int cmd_ratio(
    int    argc,
    char **argv)
{
    if (argc < 3)
    {
        fprintf(stderr, "Usage: atmturb-validate-stats ratio <fileA> <fileB> [options]\n");
        return 1;
    }
    const char *fa = argv[1], *fb = argv[2];
    double expect = 1.0, tol = 0.05;
    int is_var = 0, sub_piston = 0;

    for (int i = 3; i < argc; i++)
    {
        if (strcmp(argv[i], "--expect") == 0 && i + 1 < argc) expect = atof(argv[++i]);
        else if (strcmp(argv[i], "--tol") == 0 && i + 1 < argc) tol = atof(argv[++i]);
        else if (strcmp(argv[i], "--variance") == 0) is_var = 1;
        else if (strcmp(argv[i], "--subtract-piston") == 0) sub_piston = 1;
    }

    val_cube_t a, b;
    if (val_cube_load(fa, &a) != 0 || val_cube_load(fb, &b) != 0)
    {
        return verdict(0, "failed to load cubes for ratio");
    }

    if (sub_piston)
    {
        val_cube_subtract_piston(&a);
        val_cube_subtract_piston(&b);
    }

    double std_a = val_cube_std(&a), std_b = val_cube_std(&b);
    double meas = (std_b > 1e-30) ? (std_a / std_b) : 0.0;
    if (is_var) meas *= meas;

    double err = fabs(meas / expect - 1.0);
    char msg[256];
    snprintf(msg, sizeof(msg), "%s ratio %.5g, expected %.5g (err %.2f%%, tol %.1f%%)",
             is_var ? "variance" : "rms", meas, expect, 100.0 * err, 100.0 * tol);

    val_cube_free(&a);
    val_cube_free(&b);
    return verdict(err <= tol, msg);
}

/**
 * @brief Subcommand: corr - Pearson correlation between planes in a cube.
 */
static int cmd_corr(
    int    argc,
    char **argv)
{
    if (argc < 2)
    {
        fprintf(stderr, "Usage: atmturb-validate-stats corr <file> [options]\n");
        return 1;
    }
    const char *fname = argv[1];
    int pa = 0, pb = 1;
    double max_corr = 0.05;

    for (int i = 2; i < argc; i++)
    {
        if (strcmp(argv[i], "--plane-a") == 0 && i + 1 < argc) pa = atoi(argv[++i]);
        else if (strcmp(argv[i], "--plane-b") == 0 && i + 1 < argc) pb = atoi(argv[++i]);
        else if (strcmp(argv[i], "--max-corr") == 0 && i + 1 < argc) max_corr = atof(argv[++i]);
    }

    val_cube_t cube;
    if (val_cube_load(fname, &cube) != 0)
    {
        return verdict(0, "failed to load cube for corr");
    }

    if (pa >= cube.nframes || pb >= cube.nframes)
    {
        val_cube_free(&cube);
        return verdict(0, "plane index out of range");
    }

    double corr = val_cube_plane_corr(&cube, pa, pb);
    double abs_corr = fabs(corr);

    char msg[256];
    snprintf(msg, sizeof(msg), "correlation planes %d-%d: %.4f (max tol %.4f)",
             pa, pb, abs_corr, max_corr);
    val_cube_free(&cube);
    return verdict(abs_corr <= max_corr, msg);
}

/**
 * @brief Structure function evaluation and fitting helper.
 */
static int evaluate_sf_fit(
    int           nlags,
    const int    *lags,
    const double *dmeas,
    const double *dref,
    double        r0,
    double        tol,
    double        lag_tol,
    const char   *ref_name)
{
    double sum_ratio = 0.0, worst = 0.0;
    for (int i = 0; i < nlags; i++)
    {
        double ratio = (dref[i] > 1e-30) ? (dmeas[i] / dref[i]) : 1.0;
        sum_ratio += ratio;
        double dev = fabs(ratio - 1.0);
        if (dev > worst) worst = dev;
    }

    double mean_ratio = sum_ratio / (double) nlags;
    double r0_fit = r0 * pow(mean_ratio, -0.6);
    double err = fabs(r0_fit / r0 - 1.0);
    int ok = (err <= tol) && (lag_tol <= 0.0 || worst <= lag_tol);

    char msg[384];
    snprintf(msg, sizeof(msg),
             "sf[%s] r0 expected %.5g m, fitted %.5g m (err %.2f%%, tol %.1f%%), "
             "worst lag dev %.1f%% over lags %d..%d",
             ref_name, r0, r0_fit, 100.0 * err, 100.0 * tol,
             100.0 * worst, lags[0], lags[nlags - 1]);
    return verdict(ok, msg);
}

/**
 * @brief Subcommand: sf - check phase structure function against theory.
 */
static int cmd_sf(
    int    argc,
    char **argv)
{
    if (argc < 2) return 1;
    const char *fname = argv[1];
    char ref[32] = "theory";
    double r0_in = 0.0, seeing = 0.6, lam = 0.5e-6, zenith = 0.0;
    double L0 = 0.0, l0 = 0.0, pixscale = 0.01, tol = 0.10, lag_tol = 0.0;
    int msize = 4096, os = 1, lmin = 1, lmax = 0, sub_piston = 0;

    for (int i = 2; i < argc; i++)
    {
        if (strcmp(argv[i], "--reference") == 0 && i + 1 < argc) strncpy(ref, argv[++i], 31);
        else if (strcmp(argv[i], "--r0") == 0 && i + 1 < argc) r0_in = atof(argv[++i]);
        else if (strcmp(argv[i], "--seeing") == 0 && i + 1 < argc) seeing = atof(argv[++i]);
        else if (strcmp(argv[i], "--lam") == 0 && i + 1 < argc) lam = atof(argv[++i]);
        else if (strcmp(argv[i], "--zenith") == 0 && i + 1 < argc) zenith = atof(argv[++i]);
        else if (strcmp(argv[i], "--L0") == 0 && i + 1 < argc) L0 = atof(argv[++i]);
        else if (strcmp(argv[i], "--l0") == 0 && i + 1 < argc) l0 = atof(argv[++i]);
        else if (strcmp(argv[i], "--master-size") == 0 && i + 1 < argc) msize = atoi(argv[++i]);
        else if (strcmp(argv[i], "--oversample") == 0 && i + 1 < argc) os = atoi(argv[++i]);
        else if (strcmp(argv[i], "--pixscale") == 0 && i + 1 < argc) pixscale = atof(argv[++i]);
        else if (strcmp(argv[i], "--lag-min") == 0 && i + 1 < argc) lmin = atoi(argv[++i]);
        else if (strcmp(argv[i], "--lag-max") == 0 && i + 1 < argc) lmax = atoi(argv[++i]);
        else if (strcmp(argv[i], "--tol") == 0 && i + 1 < argc) tol = atof(argv[++i]);
        else if (strcmp(argv[i], "--lag-tol") == 0 && i + 1 < argc) lag_tol = atof(argv[++i]);
        else if (strcmp(argv[i], "--subtract-piston") == 0) sub_piston = 1;
    }

    val_cube_t cube;
    if (val_cube_load(fname, &cube) != 0) return verdict(0, "failed to load cube for sf");
    if (sub_piston) val_cube_subtract_piston(&cube);

    if (lmax <= 0) lmax = (cube.nx / 16 > lmin + 1) ? (int)(cube.nx / 16) : (lmin + 1);
    int nlags = lmax - lmin + 1;
    int *lags = (int *) malloc(nlags * sizeof(int));
    double *dmeas = (double *) malloc(nlags * sizeof(double));
    double *dref = (double *) malloc(nlags * sizeof(double));
    for (int i = 0; i < nlags; i++) lags[i] = lmin + i;

    val_measure_sf(&cube, nlags, lags, dmeas);
    double r0 = resolve_r0(r0_in, seeing, lam, zenith);

    if (strcmp(ref, "discrete") == 0)
    {
        int *m_lags = (int *) malloc(nlags * sizeof(int));
        for (int i = 0; i < nlags; i++) m_lags[i] = lags[i] * os;
        val_sf_discrete(msize, pixscale / (double) os, r0, L0, l0, nlags, m_lags, dref);
        free(m_lags);
    }
    else
    {
        for (int i = 0; i < nlags; i++)
        {
            dref[i] = val_sf_vonkarman(lags[i] * pixscale, r0, L0);
        }
    }

    int res = evaluate_sf_fit(nlags, lags, dmeas, dref, r0, tol, lag_tol, ref);
    free(lags); free(dmeas); free(dref); val_cube_free(&cube);
    return res;
}

/**
 * @brief Subcommand: shift - check frame displacement via cross-correlation.
 */
static int cmd_shift(
    int    argc,
    char **argv)
{
    if (argc < 2) return 1;
    const char *fa = argv[1], *fb = NULL;
    double edx = 0.0, edy = 0.0, tol = 0.10;
    int step = 1, argi = 2;

    if (argi < argc && argv[argi][0] != '-') fb = argv[argi++];
    for (; argi < argc; argi++)
    {
        if (strcmp(argv[argi], "--dx") == 0 && argi + 1 < argc) edx = atof(argv[++argi]);
        else if (strcmp(argv[argi], "--dy") == 0 && argi + 1 < argc) edy = atof(argv[++argi]);
        else if (strcmp(argv[argi], "--tol") == 0 && argi + 1 < argc) tol = atof(argv[++argi]);
        else if (strcmp(argv[argi], "--step") == 0 && argi + 1 < argc) step = atoi(argv[++argi]);
    }

    val_cube_t a, b;
    if (val_cube_load(fa, &a) != 0) return verdict(0, "failed to load primary shift cube");
    long npairs = (fb != NULL) ? a.nframes : (a.nframes - step);
    if (fb != NULL && val_cube_load(fb, &b) != 0)
    {
        val_cube_free(&a);
        return verdict(0, "failed to load secondary shift cube");
    }

    double sum_dx = 0.0, sum_dy = 0.0;
    long fsize = a.nx * a.ny;
    for (long f = 0; f < npairs; f++)
    {
        const double *f0 = (fb != NULL) ? &a.data[f * fsize] : &a.data[f * fsize];
        const double *f1 = (fb != NULL) ? &b.data[f * fsize] : &a.data[(f + step) * fsize];
        double dx = 0.0, dy = 0.0;
        val_xcorr_shift(f0, f1, a.nx, a.ny, &dx, &dy);
        sum_dx += dx; sum_dy += dy;
    }

    double m_dx = (npairs > 0) ? (sum_dx / (double) npairs) : 0.0;
    double m_dy = (npairs > 0) ? (sum_dy / (double) npairs) : 0.0;
    double err = hypot(m_dx - edx, m_dy - edy);

    char msg[256];
    snprintf(msg, sizeof(msg),
             "mean shift (%.3f, %.3f) px, expected (%.3f, %.3f) (err %.3f, tol %.3f)",
             m_dx, m_dy, edx, edy, err, tol);

    val_cube_free(&a);
    if (fb != NULL) val_cube_free(&b);
    return verdict(err <= tol, msg);
}

/**
 * @brief Subcommand: tilt - tip/tilt variance over circular pupil.
 */
static int cmd_tilt(
    int    argc,
    char **argv)
{
    if (argc < 2) return 1;
    const char *fname = argv[1];
    double r0_in = 0.0, seeing = 0.6, lam = 0.5e-6, zenith = 0.0;
    double L0 = 0.0, pixscale = 0.01, pupil_diam = 0.0, tol = 0.20;

    for (int i = 2; i < argc; i++)
    {
        if (strcmp(argv[i], "--r0") == 0 && i + 1 < argc) r0_in = atof(argv[++i]);
        else if (strcmp(argv[i], "--seeing") == 0 && i + 1 < argc) seeing = atof(argv[++i]);
        else if (strcmp(argv[i], "--lam") == 0 && i + 1 < argc) lam = atof(argv[++i]);
        else if (strcmp(argv[i], "--zenith") == 0 && i + 1 < argc) zenith = atof(argv[++i]);
        else if (strcmp(argv[i], "--L0") == 0 && i + 1 < argc) L0 = atof(argv[++i]);
        else if (strcmp(argv[i], "--pixscale") == 0 && i + 1 < argc) pixscale = atof(argv[++i]);
        else if (strcmp(argv[i], "--pupil-diam") == 0 && i + 1 < argc) pupil_diam = atof(argv[++i]);
        else if (strcmp(argv[i], "--tol") == 0 && i + 1 < argc) tol = atof(argv[++i]);
    }

    val_cube_t cube;
    if (val_cube_load(fname, &cube) != 0) return verdict(0, "failed to load cube for tilt");

    double var = val_measure_tilt_variance(&cube);
    double r0 = resolve_r0(r0_in, seeing, lam, zenith);
    long min_dim = (cube.nx < cube.ny) ? cube.nx : cube.ny;
    double diam = (pupil_diam > 0.0) ? pupil_diam : (double) min_dim * pixscale;
    double expect = val_tilt_variance_theory(diam, r0, L0);
    double err = (expect > 0.0) ? fabs(var / expect - 1.0) : 0.0;

    char msg[256];
    snprintf(msg, sizeof(msg),
             "tip/tilt variance %.4g rad^2, expected %.4g (err %.1f%%, tol %.0f%%)",
             var, expect, 100.0 * err, 100.0 * tol);

    val_cube_free(&cube);
    return verdict(err <= tol, msg);
}

/**
 * @brief Subcommand: repeat - temporal correlation at fixed lag.
 */
static int cmd_repeat(
    int    argc,
    char **argv)
{
    if (argc < 2) return 1;
    const char *fname = argv[1];
    int lag = 1;
    double tol = 0.05;

    for (int i = 2; i < argc; i++)
    {
        if (strcmp(argv[i], "--lag") == 0 && i + 1 < argc) lag = atoi(argv[++i]);
        else if (strcmp(argv[i], "--tol") == 0 && i + 1 < argc) tol = atof(argv[++i]);
    }

    val_cube_t cube;
    if (val_cube_load(fname, &cube) != 0) return verdict(0, "failed to load cube for repeat");
    val_cube_subtract_piston(&cube);

    double worst = val_lag_repeat_corr(&cube, lag);
    char msg[256];
    snprintf(msg, sizeof(msg), "lag %d max frame correlation = %.3f (tol %.3f)", lag, worst, tol);

    val_cube_free(&cube);
    return verdict(worst <= tol, msg);
}

/**
 * @brief Subcommand: scint - scintillation index of amplitude cube.
 */
static int cmd_scint(
    int    argc,
    char **argv)
{
    if (argc < 2) return 1;
    const char *fname = argv[1];
    double tol_mean = 0.05, tol_sigma = 0.10;

    for (int i = 2; i < argc; i++)
    {
        if (strcmp(argv[i], "--tol-mean") == 0 && i + 1 < argc) tol_mean = atof(argv[++i]);
        else if (strcmp(argv[i], "--tol-sigma") == 0 && i + 1 < argc) tol_sigma = atof(argv[++i]);
    }

    val_cube_t cube;
    if (val_cube_load(fname, &cube) != 0) return verdict(0, "failed to load cube for scint");

    double mean_i = 0.0, scint = 0.0;
    val_scintillation(&cube, &mean_i, &scint);

    int ok = (fabs(mean_i - 1.0) <= tol_mean) && (scint <= tol_sigma);
    char msg[256];
    snprintf(msg, sizeof(msg), "mean intensity %.4f (tol %.3f), scint index %.4f (tol %.3f)",
             mean_i, tol_mean, scint, tol_sigma);

    val_cube_free(&cube);
    return verdict(ok, msg);
}

/**
 * @brief Subcommand: breathing - peak-to-peak high-frequency power variation across frames.
 */
static int cmd_breathing(
    int    argc,
    char **argv)
{
    if (argc < 2) return 1;
    const char *fname = argv[1];
    double tol = 0.20;

    for (int i = 2; i < argc; i++)
    {
        if (strcmp(argv[i], "--tol") == 0 && i + 1 < argc) tol = atof(argv[++i]);
    }

    val_cube_t cube;
    if (val_cube_load(fname, &cube) != 0)
    {
        return verdict(0, "failed to load cube for breathing test");
    }

    double ratio = val_cube_breathing_ratio(&cube);
    int ok = (ratio <= tol);
    char msg[256];
    snprintf(msg, sizeof(msg), "breathing ratio = %.4f (tol %.4f)", ratio, tol);

    val_cube_free(&cube);
    return verdict(ok, msg);
}

/**
 * @brief Subcommand: simd-parity - test bitwise parity across CPU SIMD backends.
 */
static int cmd_simd_parity(
    int    argc,
    char **argv)
{
    double tol = 1e-4;
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "--tol") == 0 && i + 1 < argc) tol = atof(argv[++i]);
    }
    char msg[256];
    int ok = val_check_simd_parity(tol, msg, sizeof(msg));
    return verdict(ok, msg);
}

/**
 * @brief Print usage summary.
 */
static void print_usage(void)
{
    printf("Usage: atmturb-validate-stats <subcommand> [args...]\n\n"
           "Subcommands:\n"
           "  finite <file> [--min-std <val>]\n"
           "  sf <file> [options]\n"
           "  ratio <fileA> <fileB> [options]\n"
           "  same <fileA> <fileB> [--tol <val>]\n"
           "  shift <fileA> [<fileB>] [options]\n"
           "  tilt <file> [options]\n"
           "  repeat <file> [--lag <val>] [--tol <val>]\n"
           "  scint <file> [options]\n"
           "  corr <file> [options]\n"
           "  breathing <file> [--tol <val>]\n"
           "  simd-parity\n");
}

/**
 * @brief Main entry point dispatching to subcommands.
 */
int main(
    int   argc,
    char *argv[])
{
    if (argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
    {
        print_usage();
        return 0;
    }

    const char *cmd = argv[1];
    if (strcmp(cmd, "finite") == 0) return cmd_finite(argc - 1, &argv[1]);
    if (strcmp(cmd, "sf") == 0) return cmd_sf(argc - 1, &argv[1]);
    if (strcmp(cmd, "ratio") == 0) return cmd_ratio(argc - 1, &argv[1]);
    if (strcmp(cmd, "same") == 0) return cmd_same(argc - 1, &argv[1]);
    if (strcmp(cmd, "shift") == 0) return cmd_shift(argc - 1, &argv[1]);
    if (strcmp(cmd, "tilt") == 0) return cmd_tilt(argc - 1, &argv[1]);
    if (strcmp(cmd, "repeat") == 0) return cmd_repeat(argc - 1, &argv[1]);
    if (strcmp(cmd, "scint") == 0) return cmd_scint(argc - 1, &argv[1]);
    if (strcmp(cmd, "corr") == 0) return cmd_corr(argc - 1, &argv[1]);
    if (strcmp(cmd, "breathing") == 0) return cmd_breathing(argc - 1, &argv[1]);
    if (strcmp(cmd, "simd-parity") == 0) return cmd_simd_parity(argc - 1, &argv[1]);

    fprintf(stderr, "Unknown subcommand: %s\n", cmd);
    return 1;
}
