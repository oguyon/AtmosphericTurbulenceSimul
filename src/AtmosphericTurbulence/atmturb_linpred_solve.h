// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_linpred_solve.h
 * @brief   Regularized normal-equation Cholesky solver for AR linear predictors (CPU & GPU)
 */

#ifndef ATMTURB_LINPRED_SOLVE_H
#define ATMTURB_LINPRED_SOLVE_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * struct atmturb_linpred_problem_t - Linear least-squares problem A * X = C
 * @mat_a:    Row-major history matrix [nbmvec x mvecsize].
 * @mat_c:    Row-major target matrix [nbmvec x nbpix].
 * @filt_out: Row-major output predictor matrix [nbpix x mvecsize] (X^T).
 * @nbmvec:   Number of training sample rows (N).
 * @mvecsize: Feature vector dimension (M).
 * @nbpix:    Target output dimension (P).
 * @reg_eps:  Relative singular-value regularization threshold.
 */
typedef struct
{
    const float *mat_a;
    const float *mat_c;
    float       *filt_out;
    long         nbmvec;
    long         mvecsize;
    long         nbpix;
    double       reg_eps;
} atmturb_linpred_problem_t;

/**
 * atmturb_linpred_solve_normal - Solve regularized normal equations (A^T A + mu I) X = A^T C
 * @prob: Problem descriptor and input/output host buffers.
 *
 * Return: 0 on success, -1 on invalid dimensions or solver failure.
 */
int atmturb_linpred_solve_normal(
    const atmturb_linpred_problem_t *prob);

#ifdef HAVE_CUDA

/**
 * atmturb_cuda_linpred_solve - Solve regularized normal equations on GPU with cuBLAS + cuSOLVER
 * @prob: Problem descriptor and input/output host buffers.
 *
 * Return: 0 on success, -1 on CUDA/cuSOLVER failure (caller may fall back to CPU).
 */
int atmturb_cuda_linpred_solve(
    const atmturb_linpred_problem_t *prob);

/**
 * atmturb_cuda_linpred_cleanup - Release cached cuBLAS/cuSOLVER handles and scratch buffers
 */
void atmturb_cuda_linpred_cleanup(void);

#endif /* HAVE_CUDA */

#ifdef __cplusplus
}
#endif

#endif /* ATMTURB_LINPRED_SOLVE_H */
