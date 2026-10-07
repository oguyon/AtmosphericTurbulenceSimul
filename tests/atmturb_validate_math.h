// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_validate_math.h
 * @brief   Core math, I/O, and statistical primitives for milkatmturb validation
 */

#ifndef ATMTURB_VALIDATE_MATH_H
#define ATMTURB_VALIDATE_MATH_H

#include <stddef.h>

/**
 * @brief In-memory representation of a 2D image or 3D cube.
 */
typedef struct
{
    long    naxis;
    long    nframes;
    long    ny;
    long    nx;
    long    nelements;
    double *data;
} val_cube_t;

/**
 * @brief Load a 2D or 3D FITS file into a float64 array.
 */
int val_cube_load(
    const char *filename,
    val_cube_t *cube);

/**
 * @brief Free heap memory allocated for a val_cube_t.
 */
void val_cube_free(
    val_cube_t *cube);

/**
 * @brief Subtract mean per frame (piston removal).
 */
void val_cube_subtract_piston(
    val_cube_t *cube);

/**
 * @brief Check if all elements are finite and compute sample standard deviation.
 */
int val_cube_check_finite(
    const val_cube_t *cube,
    double           *out_std);

/**
 * @brief Compute overall sample standard deviation.
 */
double val_cube_std(
    const val_cube_t *cube);

/**
 * @brief Continuous von Karman structure function [rad^2].
 */
double val_sf_vonkarman(
    double r,
    double r0,
    double L0);

/**
 * @brief Ensemble structure function of a discrete periodic FFT master screen.
 */
void val_sf_discrete(
    int        nmaster,
    double     dx_m,
    double     r0,
    double     L0,
    double     l0,
    int        nlags,
    const int *lags,
    double    *dref);

/**
 * @brief Measure empirical structure function along x and y.
 */
void val_measure_sf(
    const val_cube_t *cube,
    int               nlags,
    const int        *lags,
    double           *dmeas);

/**
 * @brief Theoretical two-axis Zernike tip+tilt variance over circular pupil.
 */
double val_tilt_variance_theory(
    double D,
    double r0,
    double L0);

/**
 * @brief Measure empirical tip+tilt variance over circular pupil.
 */
double val_measure_tilt_variance(
    const val_cube_t *cube);

/**
 * @brief Sub-pixel 2D shift via FFT cross-correlation peak.
 */
void val_xcorr_shift(
    const double *f0,
    const double *f1,
    long          nx,
    long          ny,
    double       *out_dx,
    double       *out_dy);

/**
 * @brief Pearson correlation between two planes in a 3D cube.
 */
double val_cube_plane_corr(
    const val_cube_t *cube,
    int               plane_a,
    int               plane_b);

/**
 * @brief Pearson correlation between frames separated by fixed lag.
 */
double val_lag_repeat_corr(
    const val_cube_t *cube,
    int               lag);

/**
 * @brief Compute mean intensity and scintillation index of amplitude cube.
 */
void val_scintillation(
    const val_cube_t *cube,
    double           *out_mean,
    double           *out_scint_index);

/**
 * @brief Measure peak-to-peak high-frequency power variation across frames (breathing metric).
 */
double val_cube_breathing_ratio(
    const val_cube_t *cube);

/**
 * @brief Compare scalar and vectorized SIMD extrusions across schemes and strides.
 */
int val_check_simd_parity(
    double  tol,
    char   *msg_buf,
    size_t  msg_size);

#endif /* ATMTURB_VALIDATE_MATH_H */
