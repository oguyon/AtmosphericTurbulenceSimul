// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_linpred_solve.c
 * @brief   CPU (OpenBLAS/LAPACKE + pure-C fallback) and GPU normal-equation solver
 */

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef HAVE_OPENBLAS_LAPACKE
#include <cblas.h>
#include <lapacke.h>
#endif

#include "atmturb_linpred_solve.h"
#include "atmturb_simd.h"

#ifdef HAVE_OPENBLAS_LAPACKE

/**
 * linpred_estimate_lmax_blas - Estimate spectral radius of symmetric G via 8 power steps
 */
static double linpred_estimate_lmax_blas(
    const float *g_mat,
    int          m)
{
    float *vw = (float *) malloc(sizeof(float) * 2U * (size_t) m);
    if (vw == NULL)
    {
        return 1.0;
    }

    float *v          = vw;
    float *w          = vw + m;
    float  inv_sqrt_m = 1.0f / sqrtf((float) m);
    for (int i = 0; i < m; i++)
    {
        v[i] = inv_sqrt_m;
    }

    double rq = 1.0;
    for (int it = 0; it < 8; it++)
    {
        cblas_ssymv(CblasColMajor, CblasUpper, m, 1.0f, g_mat, m, v, 1, 0.0f, w, 1);
        double dot_vw = 0.0;
        double dot_ww = 0.0;
        for (int i = 0; i < m; i++)
        {
            double wi = (double) w[i];
            dot_vw += (double) v[i] * wi;
            dot_ww += wi * wi;
        }
        rq = dot_vw;
        if (dot_ww <= 0.0)
        {
            break;
        }
        float inv_n = (float) (1.0 / sqrt(dot_ww));
        for (int i = 0; i < m; i++)
        {
            v[i] = w[i] * inv_n;
        }
    }

    free(vw);
    return (rq > 0.0) ? rq : 1.0;
}

/**
 * linpred_eig_fallback_blas - FP64 symmetric eigensolver fallback for rank-deficient inputs
 */
static int linpred_eig_fallback_blas(
    const atmturb_linpred_problem_t *prob,
    double                           mu)
{
    int     m     = (int) prob->mvecsize;
    int     n     = (int) prob->nbmvec;
    int     p     = (int) prob->nbpix;
    size_t  mm    = (size_t) m * (size_t) m;
    size_t  mp    = (size_t) m * (size_t) p;
    double *gd    = (double *) calloc(mm, sizeof(double));
    double *evals = (double *) malloc(sizeof(double) * (size_t) m);
    double *rd    = (double *) malloc(sizeof(double) * mp);
    double *zd    = (double *) malloc(sizeof(double) * mp);
    if (gd == NULL || evals == NULL || rd == NULL || zd == NULL)
    {
        free(gd);
        free(evals);
        free(rd);
        free(zd);
        return -1;
    }

    float *gf = (float *) malloc(sizeof(float) * mm);
    if (gf == NULL)
    {
        free(gd);
        free(evals);
        free(rd);
        free(zd);
        return -1;
    }
    cblas_ssyrk(CblasColMajor, CblasUpper, CblasNoTrans, m, n, 1.0f, prob->mat_a, m, 0.0f, gf, m);
    for (size_t idx = 0; idx < mm; idx++)
    {
        gd[idx] = (double) gf[idx];
    }
    free(gf);
    for (size_t idx = 0; idx < mp; idx++)
    {
        rd[idx] = (double) prob->filt_out[idx];
    }

    int info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', m, gd, m, evals);
    if (info == 0)
    {
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    m, p, m, 1.0, gd, m, rd, m, 0.0, zd, m);
        for (int col = 0; col < p; col++)
        {
            for (int k = 0; k < m; k++)
            {
                double ev                                  = (evals[k] > 0.0) ? evals[k] : 0.0;
                zd[(size_t) col * (size_t) m + (size_t) k] /= (ev + mu);
            }
        }
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    m, p, m, 1.0, gd, m, zd, m, 0.0, rd, m);
        for (size_t idx = 0; idx < mp; idx++)
        {
            prob->filt_out[idx] = (float) rd[idx];
        }
    }

    free(gd);
    free(evals);
    free(rd);
    free(zd);
    return (info == 0) ? 0 : -1;
}

/**
 * linpred_solve_cpu_blas - Solve regularized normal equations via OpenBLAS + LAPACKE Cholesky
 */
static int linpred_solve_cpu_blas(
    const atmturb_linpred_problem_t *prob)
{
    int    m     = (int) prob->mvecsize;
    int    n     = (int) prob->nbmvec;
    int    p     = (int) prob->nbpix;
    float *g_mat = (float *) malloc(sizeof(float) * (size_t) m * (size_t) m);
    if (g_mat == NULL)
    {
        return -1;
    }

    cblas_ssyrk(CblasColMajor, CblasUpper, CblasNoTrans,
                m, n, 1.0f, prob->mat_a, m, 0.0f, g_mat, m);
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                m, p, n, 1.0f, prob->mat_a, m, prob->mat_c, p, 0.0f, prob->filt_out, m);

    double lmax   = linpred_estimate_lmax_blas(g_mat, m);
    double eps_sq = prob->reg_eps * prob->reg_eps;
    double rel    = (eps_sq > 1e-6) ? eps_sq : 1e-6;
    double mu     = rel * lmax;
    float  mu_f   = (float) mu;
    for (int i = 0; i < m; i++)
    {
        g_mat[(size_t) i * (size_t) m + (size_t) i] += mu_f;
    }

    int info = LAPACKE_spotrf(LAPACK_COL_MAJOR, 'U', m, g_mat, m);
    if (info == 0)
    {
        info = LAPACKE_spotrs(LAPACK_COL_MAJOR, 'U', m, p, g_mat, m, prob->filt_out, m);
    }
    free(g_mat);

    if (info != 0)
    {
        return linpred_eig_fallback_blas(prob, mu);
    }
    return 0;
}

#else /* !HAVE_OPENBLAS_LAPACKE */

/**
 * linpred_cholesky_solve_c - Pure-C FP64 Cholesky decomposition and multi-RHS triangular solve
 */
static int linpred_cholesky_solve_c(
    double *l_mat,
    double *rhs,
    long    m,
    long    p,
    float  *filt_out)
{
    for (long i = 0; i < m; i++)
    {
        for (long j = 0; j <= i; j++)
        {
            double sum = l_mat[i * m + j];
            for (long k = 0; k < j; k++)
            {
                sum -= l_mat[i * m + k] * l_mat[j * m + k];
            }
            if (i == j)
            {
                if (sum <= 0.0)
                {
                    return -1;
                }
                l_mat[i * m + i] = sqrt(sum);
            }
            else
            {
                l_mat[i * m + j] = sum / l_mat[j * m + j];
            }
        }
    }

    #pragma omp parallel for schedule(static)
    for (long col = 0; col < p; col++)
    {
        double *b = &rhs[col * m];
        for (long i = 0; i < m; i++)
        {
            double sum = b[i];
            for (long k = 0; k < i; k++)
            {
                sum -= l_mat[i * m + k] * b[k];
            }
            b[i] = sum / l_mat[i * m + i];
        }
        for (long i = m - 1; i >= 0; i--)
        {
            double sum = b[i];
            for (long k = i + 1; k < m; k++)
            {
                sum -= l_mat[k * m + i] * b[k];
            }
            b[i]                  = sum / l_mat[i * m + i];
            filt_out[col * m + i] = (float) b[i];
        }
    }
    return 0;
}

/**
 * linpred_solve_cpu_fallback - Pure-C + OpenMP regularized normal-equation solver
 */
static int linpred_solve_cpu_fallback(
    const atmturb_linpred_problem_t *prob)
{
    long    m   = prob->mvecsize;
    long    n   = prob->nbmvec;
    long    p   = prob->nbpix;
    double *g   = (double *) calloc((size_t) m * (size_t) m, sizeof(double));
    double *rhs = (double *) calloc((size_t) p * (size_t) m, sizeof(double));
    if (g == NULL || rhs == NULL)
    {
        free(g);
        free(rhs);
        return -1;
    }

    #pragma omp parallel for schedule(dynamic, 4)
    for (long i = 0; i < m; i++)
    {
        for (long s = 0; s < n; s++)
        {
            double       ai    = (double) prob->mat_a[s * m + i];
            const float *a_row = &prob->mat_a[s * m];
            for (long j = 0; j <= i; j++)
            {
                g[i * m + j] += ai * (double) a_row[j];
            }
        }
    }

    #pragma omp parallel for schedule(static)
    for (long col = 0; col < p; col++)
    {
        double *r_col = &rhs[col * m];
        for (long s = 0; s < n; s++)
        {
            double       cs    = (double) prob->mat_c[s * p + col];
            const float *a_row = &prob->mat_a[s * m];
            for (long i = 0; i < m; i++)
            {
                r_col[i] += (double) a_row[i] * cs;
            }
        }
    }

    double trace = 0.0;
    for (long i = 0; i < m; i++)
    {
        trace += g[i * m + i];
    }
    double eps_sq = prob->reg_eps * prob->reg_eps;
    double mu     = ((eps_sq > 1e-6) ? eps_sq : 1e-6) * ((trace > 0.0) ? trace : 1.0);
    for (long i = 0; i < m; i++)
    {
        g[i * m + i] += mu;
    }

    int ret = linpred_cholesky_solve_c(g, rhs, m, p, prob->filt_out);
    free(g);
    free(rhs);
    return ret;
}

#endif /* HAVE_OPENBLAS_LAPACKE */

/**
 * atmturb_linpred_solve_normal - Solve regularized normal equations (A^T A + mu I) X = A^T C
 * @prob: Problem descriptor and input/output host buffers.
 *
 * Return: 0 on success, -1 on invalid dimensions or solver failure.
 */
int atmturb_linpred_solve_normal(
    const atmturb_linpred_problem_t *prob)
{
    if (prob == NULL || prob->mat_a == NULL || prob->mat_c == NULL || prob->filt_out == NULL ||
        prob->nbmvec <= 0 || prob->mvecsize <= 0 || prob->nbpix <= 0)
    {
        return -1;
    }

#ifdef HAVE_CUDA
    if (atmturb_simd_is_gpu() && (prob->mvecsize >= 128 || atmturb_simd_is_gpu_forced()))
    {
        if (atmturb_cuda_linpred_solve(prob) == 0)
        {
            return 0;
        }
    }
#endif

#ifdef HAVE_OPENBLAS_LAPACKE
    return linpred_solve_cpu_blas(prob);
#else
    return linpred_solve_cpu_fallback(prob);
#endif
}
