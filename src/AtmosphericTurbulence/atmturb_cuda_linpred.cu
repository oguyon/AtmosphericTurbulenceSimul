// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_linpred.cu
 * @brief   GPU cuBLAS + cuSOLVER Cholesky normal-equation solver for AR linear predictors
 */

#ifdef HAVE_CUDA

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <math.h>
#include <stddef.h>

#include "atmturb_cuda.h"
#include "atmturb_linpred_solve.h"

namespace {

struct linpred_gpu_cache_t
{
    cublasHandle_t     cublas;
    cusolverDnHandle_t cusolver;
    cudaStream_t       stream;
    int                handles_ready;
    float             *d_a;
    size_t             cap_a;
    float             *d_c;
    size_t             cap_c;
    float             *d_g;
    size_t             cap_g;
    float             *d_r;
    size_t             cap_r;
    float             *d_vw;
    size_t             cap_vw;
    float             *d_work;
    size_t             cap_work;
    double            *d_lmax;
    int               *d_info;
};

static linpred_gpu_cache_t s_lp = {};

__global__ void linpred_init_vec_kernel(
    float *__restrict__ d_v,
    int                 m,
    float               inv_sqrt_m)
{
    int idx = (int) blockIdx.x * (int) blockDim.x + (int) threadIdx.x;
    if (idx < m)
    {
        d_v[idx] = inv_sqrt_m;
    }
}

__global__ void linpred_rayleigh_norm_kernel(
    float        *__restrict__ d_v,
    const float  *__restrict__ d_w,
    double       *__restrict__ d_lmax,
    int                        m)
{
    __shared__ double s_rq[256];
    __shared__ double s_wn[256];

    int    tid = (int) threadIdx.x;
    double rq  = 0.0;
    double wn  = 0.0;

    for (int i = tid; i < m; i += 256)
    {
        double vi = (double) d_v[i];
        double wi = (double) d_w[i];
        rq += vi * wi;
        wn += wi * wi;
    }
    s_rq[tid] = rq;
    s_wn[tid] = wn;
    __syncthreads();

    for (int stride = 128; stride > 0; stride >>= 1)
    {
        if (tid < stride)
        {
            s_rq[tid] += s_rq[tid + stride];
            s_wn[tid] += s_wn[tid + stride];
        }
        __syncthreads();
    }

    if (tid == 0)
    {
        d_lmax[0] = s_rq[0];
    }

    double inv_norm = (s_wn[0] > 0.0) ? (1.0 / sqrt(s_wn[0])) : 0.0;
    for (int i = tid; i < m; i += 256)
    {
        d_v[i] = (float) ((double) d_w[i] * inv_norm);
    }
}

__global__ void linpred_add_reg_kernel(
    float        *__restrict__ d_g,
    int                        m,
    const double *__restrict__ d_lmax,
    double                     eps_sq)
{
    int idx = (int) blockIdx.x * (int) blockDim.x + (int) threadIdx.x;
    if (idx < m)
    {
        double lmax = d_lmax[0];
        double base = (lmax > 0.0) ? lmax : 1.0;
        double rel  = (eps_sq > 1e-6) ? eps_sq : 1e-6;
        float  mu   = (float) (rel * base);
        d_g[(size_t) idx * (size_t) m + (size_t) idx] += mu;
    }
}

static int linpred_grow_buf(
    void  **ptr,
    size_t *cap,
    size_t  need_bytes)
{
    if (*ptr != NULL && *cap >= need_bytes)
    {
        return 0;
    }
    if (*ptr != NULL)
    {
        cudaFree(*ptr);
        *ptr = NULL;
        *cap = 0;
    }
    if (cudaMalloc(ptr, need_bytes) != cudaSuccess)
    {
        *ptr = NULL;
        return -1;
    }
    *cap = need_bytes;
    return 0;
}

static int linpred_ensure_handles(void)
{
    if (s_lp.handles_ready)
    {
        return 0;
    }
    if (cudaStreamCreateWithFlags(&s_lp.stream, cudaStreamNonBlocking) != cudaSuccess)
    {
        s_lp.stream = 0;
        return -1;
    }
    if (cublasCreate(&s_lp.cublas) != CUBLAS_STATUS_SUCCESS)
    {
        cudaStreamDestroy(s_lp.stream);
        s_lp.stream = 0;
        return -1;
    }
    if (cusolverDnCreate(&s_lp.cusolver) != CUSOLVER_STATUS_SUCCESS)
    {
        cublasDestroy(s_lp.cublas);
        cudaStreamDestroy(s_lp.stream);
        s_lp.cublas = NULL;
        s_lp.stream = 0;
        return -1;
    }
    cublasSetStream(s_lp.cublas, s_lp.stream);
    cusolverDnSetStream(s_lp.cusolver, s_lp.stream);

    if (cudaMalloc((void **) &s_lp.d_lmax, sizeof(double)) != cudaSuccess ||
        cudaMalloc((void **) &s_lp.d_info, sizeof(int)) != cudaSuccess)
    {
        atmturb_cuda_linpred_cleanup();
        return -1;
    }
    s_lp.handles_ready = 1;
    return 0;
}

static int linpred_prepare_workspace(
    int m,
    int n,
    int p,
    int *out_lwork)
{
    size_t bytes_a  = sizeof(float) * (size_t) n * (size_t) m;
    size_t bytes_c  = sizeof(float) * (size_t) n * (size_t) p;
    size_t bytes_g  = sizeof(float) * (size_t) m * (size_t) m;
    size_t bytes_r  = sizeof(float) * (size_t) m * (size_t) p;
    size_t bytes_vw = sizeof(float) * 2U * (size_t) m;

    if (linpred_grow_buf((void **) &s_lp.d_a, &s_lp.cap_a, bytes_a) != 0 ||
        linpred_grow_buf((void **) &s_lp.d_c, &s_lp.cap_c, bytes_c) != 0 ||
        linpred_grow_buf((void **) &s_lp.d_g, &s_lp.cap_g, bytes_g) != 0 ||
        linpred_grow_buf((void **) &s_lp.d_r, &s_lp.cap_r, bytes_r) != 0 ||
        linpred_grow_buf((void **) &s_lp.d_vw, &s_lp.cap_vw, bytes_vw) != 0)
    {
        return -1;
    }

    int lwork = 0;
    if (cusolverDnSpotrf_bufferSize(
            s_lp.cusolver, CUBLAS_FILL_MODE_UPPER, m, s_lp.d_g, m, &lwork) !=
        CUSOLVER_STATUS_SUCCESS)
    {
        return -1;
    }
    size_t work_bytes = sizeof(float) * (size_t) ((lwork > 1) ? lwork : 1);
    if (linpred_grow_buf((void **) &s_lp.d_work, &s_lp.cap_work, work_bytes) != 0)
    {
        return -1;
    }
    *out_lwork = lwork;
    return 0;
}

} // namespace

extern "C" int atmturb_cuda_linpred_solve(
    const atmturb_linpred_problem_t *prob)
{
    if (prob == NULL || prob->mat_a == NULL || prob->mat_c == NULL || prob->filt_out == NULL ||
        prob->nbmvec <= 0 || prob->mvecsize <= 0 || prob->nbpix <= 0 ||
        !atmturb_cuda_device_available())
    {
        return -1;
    }

    int m     = (int) prob->mvecsize;
    int n     = (int) prob->nbmvec;
    int p     = (int) prob->nbpix;
    int lwork = 0;
    if (linpred_ensure_handles() != 0 || linpred_prepare_workspace(m, n, p, &lwork) != 0)
    {
        return -1;
    }

    size_t bytes_a = sizeof(float) * (size_t) n * (size_t) m;
    size_t bytes_c = sizeof(float) * (size_t) n * (size_t) p;
    size_t bytes_r = sizeof(float) * (size_t) m * (size_t) p;
    cudaMemcpyAsync(s_lp.d_a, prob->mat_a, bytes_a, cudaMemcpyHostToDevice, s_lp.stream);
    cudaMemcpyAsync(s_lp.d_c, prob->mat_c, bytes_c, cudaMemcpyHostToDevice, s_lp.stream);

    const float one  = 1.0f;
    const float zero = 0.0f;
    if (cublasSsyrk(s_lp.cublas, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N,
                    m, n, &one, s_lp.d_a, m, &zero, s_lp.d_g, m) != CUBLAS_STATUS_SUCCESS ||
        cublasSgemm(s_lp.cublas, CUBLAS_OP_N, CUBLAS_OP_T,
                    m, p, n, &one, s_lp.d_a, m, s_lp.d_c, p, &zero, s_lp.d_r, m) !=
            CUBLAS_STATUS_SUCCESS)
    {
        return -1;
    }

    float *d_v  = s_lp.d_vw;
    float *d_w  = s_lp.d_vw + m;
    int    grid = (m + 255) / 256;
    linpred_init_vec_kernel<<<grid, 256, 0, s_lp.stream>>>(d_v, m, 1.0f / sqrtf((float) m));
    for (int it = 0; it < 8; it++)
    {
        cublasSsymv(s_lp.cublas, CUBLAS_FILL_MODE_UPPER, m, &one, s_lp.d_g, m,
                    d_v, 1, &zero, d_w, 1);
        linpred_rayleigh_norm_kernel<<<1, 256, 0, s_lp.stream>>>(d_v, d_w, s_lp.d_lmax, m);
    }
    linpred_add_reg_kernel<<<grid, 256, 0, s_lp.stream>>>(
        s_lp.d_g, m, s_lp.d_lmax, prob->reg_eps * prob->reg_eps);

    if (cusolverDnSpotrf(s_lp.cusolver, CUBLAS_FILL_MODE_UPPER, m, s_lp.d_g, m,
                         s_lp.d_work, lwork, s_lp.d_info) != CUSOLVER_STATUS_SUCCESS ||
        cusolverDnSpotrs(s_lp.cusolver, CUBLAS_FILL_MODE_UPPER, m, p, s_lp.d_g, m,
                         s_lp.d_r, m, s_lp.d_info) != CUSOLVER_STATUS_SUCCESS)
    {
        return -1;
    }

    int h_info = -1;
    cudaMemcpyAsync(prob->filt_out, s_lp.d_r, bytes_r, cudaMemcpyDeviceToHost, s_lp.stream);
    cudaMemcpyAsync(&h_info, s_lp.d_info, sizeof(int), cudaMemcpyDeviceToHost, s_lp.stream);
    if (cudaStreamSynchronize(s_lp.stream) != cudaSuccess || h_info != 0)
    {
        return -1;
    }
    return 0;
}

extern "C" void atmturb_cuda_linpred_cleanup(void)
{
    if (s_lp.stream != 0)
    {
        cudaStreamSynchronize(s_lp.stream);
    }
    cudaFree(s_lp.d_a);
    cudaFree(s_lp.d_c);
    cudaFree(s_lp.d_g);
    cudaFree(s_lp.d_r);
    cudaFree(s_lp.d_vw);
    cudaFree(s_lp.d_work);
    cudaFree(s_lp.d_lmax);
    cudaFree(s_lp.d_info);
    if (s_lp.cusolver != NULL)
    {
        cusolverDnDestroy(s_lp.cusolver);
    }
    if (s_lp.cublas != NULL)
    {
        cublasDestroy(s_lp.cublas);
    }
    if (s_lp.stream != 0)
    {
        cudaStreamDestroy(s_lp.stream);
    }
    s_lp = linpred_gpu_cache_t{};
}

#endif // HAVE_CUDA
