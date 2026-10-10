// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_rytov.cu
 * @brief   CUDA and cuFFT accelerated multi-layer Rytov propagation
 */

#include <cuda_runtime.h>
#include <cufft.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "atmturb_cuda.h"
#include "atmturb_cuda_rytov.h"
#include "atmturb_cuda_rytov_internal.h"

// Persistent cuFFT plans and configuration cache
static cufftHandle s_plan_1d          = 0;
static cufftHandle s_plan_2d_r2c      = 0;
static cufftHandle s_plan_2d_c2r_bat  = 0;
static long        s_cached_pad       = 0;
static int         s_cached_has_sec   = -1;

// Persistent CUDA stream and Graph execution instance
static cudaStream_t                   s_stream            = 0;
static cudaGraphExec_t                s_graph_exec        = NULL;
static long                           s_graph_pad         = 0;
static long                           s_graph_pup         = 0;
static int                            s_graph_nsuper      = 0;
static int                            s_graph_has_sec     = 0;
static int                            s_graph_sec_shared  = 0;
static int                            s_graph_interp      = 0;
static int                            s_graph_guard       = 0;
static int                            s_graph_moisan      = -1;
static int                            s_graph_os          = 0;
static const float                   *s_graph_masters     = NULL;
static atmturb_cuda_rytov_sublayer_t *s_d_active_sublayers = NULL;
static size_t                         s_cached_active_sub = 0;

// Persistent device precomputed filter buffers
static float        *s_d_laplace_inv  = NULL;
static cufftComplex *s_d_exp_x        = NULL;
static cufftComplex *s_d_exp_y        = NULL;
static float        *s_d_filt_a_pri   = NULL;
static float        *s_d_filt_b_pri   = NULL;
static float        *s_d_filt_a_sec   = NULL;
static float        *s_d_filt_b_sec   = NULL;
static cufftComplex *s_d_chrom_ramp   = NULL;
static size_t        s_cached_filters = 0;

// Persistent device frame working buffers
static float        *s_d_super_pha    = NULL;
static float        *s_d_super_spha   = NULL;
static float        *s_d_bounds       = NULL;
static cufftComplex *s_d_hat_bounds   = NULL;
static cufftComplex *s_d_spec         = NULL;
static cufftComplex *s_d_acc_batch    = NULL;
static float        *s_d_inv_out      = NULL;
static float        *s_d_frame_pha    = NULL;
static float        *s_d_frame_amp    = NULL;
static float        *s_d_frame_spha   = NULL;
static float        *s_d_frame_samp   = NULL;
static float        *s_d_sum_I        = NULL;
static size_t        s_cached_work    = 0;

/**
 * atmturb_cuda_forward_spectrum - Compute 2D spectrum with optional Moisan adjustment
 * @d_pha: Padded device phase grid.
 * @use_moisan: 1 if Moisan boundary decomposition active, 0 to bypass.
 * @pad: Grid linear dimension in pixels.
 * @n_half: Half-complex dimension along X (pad / 2 + 1).
 * @stream: Active CUDA execution stream.
 */
static void atmturb_cuda_forward_spectrum(
    const float  *d_pha,
    int           use_moisan,
    int           pad,
    int           n_half,
    cudaStream_t  stream)
{
    if (use_moisan)
    {
        atmturb_cuda_extract_bounds_launch(d_pha, s_d_bounds, pad, stream);
        cufftExecR2C(s_plan_1d, (cufftReal *) s_d_bounds, (cufftComplex *) s_d_hat_bounds);
        cufftExecR2C(s_plan_2d_r2c, (cufftReal *) d_pha, (cufftComplex *) s_d_spec);
        atmturb_cuda_moisan_adjust_launch(
            s_d_spec, s_d_hat_bounds, s_d_exp_x, s_d_exp_y, s_d_laplace_inv, pad, n_half, stream);
    }
    else
    {
        cufftExecR2C(s_plan_2d_r2c, (cufftReal *) d_pha, (cufftComplex *) s_d_spec);
    }
}

/**
 * atmturb_cuda_rytov_cleanup - Free cached GPU buffers, plans, and stream
 */
void atmturb_cuda_rytov_cleanup(void)
{
    if (s_graph_exec != NULL)
    {
        cudaGraphExecDestroy(s_graph_exec);
        s_graph_exec = NULL;
    }
    if (s_plan_1d != 0)
    {
        cufftDestroy(s_plan_1d);
        s_plan_1d = 0;
    }
    if (s_plan_2d_r2c != 0)
    {
        cufftDestroy(s_plan_2d_r2c);
        s_plan_2d_r2c = 0;
    }
    if (s_plan_2d_c2r_bat != 0)
    {
        cufftDestroy(s_plan_2d_c2r_bat);
        s_plan_2d_c2r_bat = 0;
    }
    s_cached_pad     = 0;
    s_cached_has_sec = -1;

    if (s_stream != 0)
    {
        cudaStreamDestroy(s_stream);
        s_stream = 0;
    }
    s_graph_pad         = 0;
    s_graph_pup         = 0;
    s_graph_nsuper      = 0;
    s_graph_has_sec     = 0;
    s_graph_sec_shared  = 0;
    s_graph_interp      = 0;
    s_graph_guard       = 0;
    s_graph_moisan      = -1;
    s_graph_os          = 0;
    s_graph_masters     = NULL;

    if (s_d_active_sublayers != NULL)
    {
        cudaFree(s_d_active_sublayers);
        s_d_active_sublayers = NULL;
    }
    s_cached_active_sub = 0;

    if (s_d_laplace_inv != NULL) { cudaFree(s_d_laplace_inv); s_d_laplace_inv = NULL; }
    if (s_d_exp_x != NULL)       { cudaFree(s_d_exp_x);       s_d_exp_x = NULL; }
    if (s_d_exp_y != NULL)       { cudaFree(s_d_exp_y);       s_d_exp_y = NULL; }
    if (s_d_filt_a_pri != NULL)  { cudaFree(s_d_filt_a_pri);  s_d_filt_a_pri = NULL; }
    if (s_d_filt_b_pri != NULL)  { cudaFree(s_d_filt_b_pri);  s_d_filt_b_pri = NULL; }
    if (s_d_filt_a_sec != NULL)  { cudaFree(s_d_filt_a_sec);  s_d_filt_a_sec = NULL; }
    if (s_d_filt_b_sec != NULL)  { cudaFree(s_d_filt_b_sec);  s_d_filt_b_sec = NULL; }
    if (s_d_chrom_ramp != NULL)  { cudaFree(s_d_chrom_ramp);  s_d_chrom_ramp = NULL; }
    s_cached_filters = 0;

    if (s_d_super_pha != NULL)  { cudaFree(s_d_super_pha);  s_d_super_pha = NULL; }
    if (s_d_super_spha != NULL) { cudaFree(s_d_super_spha); s_d_super_spha = NULL; }
    if (s_d_bounds != NULL)     { cudaFree(s_d_bounds);     s_d_bounds = NULL; }
    if (s_d_hat_bounds != NULL) { cudaFree(s_d_hat_bounds); s_d_hat_bounds = NULL; }
    if (s_d_spec != NULL)       { cudaFree(s_d_spec);       s_d_spec = NULL; }
    if (s_d_acc_batch != NULL)  { cudaFree(s_d_acc_batch);  s_d_acc_batch = NULL; }
    if (s_d_inv_out != NULL)    { cudaFree(s_d_inv_out);    s_d_inv_out = NULL; }
    if (s_d_frame_pha != NULL)  { cudaFree(s_d_frame_pha);  s_d_frame_pha = NULL; }
    if (s_d_frame_amp != NULL)  { cudaFree(s_d_frame_amp);  s_d_frame_amp = NULL; }
    if (s_d_frame_spha != NULL) { cudaFree(s_d_frame_spha); s_d_frame_spha = NULL; }
    if (s_d_frame_samp != NULL) { cudaFree(s_d_frame_samp); s_d_frame_samp = NULL; }
    if (s_d_sum_I != NULL)      { cudaFree(s_d_sum_I);      s_d_sum_I = NULL; }
    s_cached_work = 0;
}

/**
 * atmturb_cuda_rytov_init_plans - Initialize persistent cuFFT plans and execution stream
 * @pad_size: Linear dimension of padded compute grid.
 * @has_sec: 1 if secondary wavelength enabled, 0 otherwise.
 *
 * Return: 0 on success, -1 on failure.
 */
static int atmturb_cuda_rytov_init_plans(
    long pad_size,
    int  has_sec)
{
    if (s_stream == 0)
    {
        if (cudaStreamCreate(&s_stream) != cudaSuccess)
        {
            return -1;
        }
    }

    if (s_plan_1d != 0 && s_cached_pad == pad_size && s_cached_has_sec == has_sec)
    {
        return 0;
    }
    atmturb_cuda_rytov_cleanup();

    if (cudaStreamCreate(&s_stream) != cudaSuccess)
    {
        return -1;
    }

    int n_inv = has_sec ? 4 : 2;
    int n[2] = {(int) pad_size, (int) pad_size};
    long n_half = pad_size / 2 + 1;
    long n_spec = pad_size * n_half;
    long pad_pixels = pad_size * pad_size;

    if (cufftPlan1d(&s_plan_1d, (int) pad_size, CUFFT_R2C, 2) != CUFFT_SUCCESS ||
        cufftPlan2d(&s_plan_2d_r2c, (int) pad_size, (int) pad_size, CUFFT_R2C) != CUFFT_SUCCESS ||
        cufftPlanMany(&s_plan_2d_c2r_bat, 2, n,
                      NULL, 1, (int) n_spec,
                      NULL, 1, (int) pad_pixels,
                      CUFFT_C2R, n_inv) != CUFFT_SUCCESS)
    {
        atmturb_cuda_rytov_cleanup();
        return -1;
    }

    cufftSetStream(s_plan_1d, s_stream);
    cufftSetStream(s_plan_2d_r2c, s_stream);
    cufftSetStream(s_plan_2d_c2r_bat, s_stream);

    s_cached_pad     = pad_size;
    s_cached_has_sec = has_sec;
    return 0;
}

/**
 * atmturb_cuda_rytov_init_filters - Precompute and upload static Rytov filters to device
 * @params: Propagation configuration and host filter pointers.
 *
 * Return: 0 on success, -1 on failure.
 */
static int atmturb_cuda_rytov_init_filters(
    const atmturb_cuda_rytov_params_t *params)
{
    long pad = params->pad_size;
    long n_half = pad / 2 + 1;
    long n_spec = pad * n_half;
    int ns = params->nsuper;
    size_t f_bytes = (size_t) (ns * n_spec) * sizeof(float);

    if (s_d_filt_a_pri != NULL && s_cached_filters == f_bytes)
    {
        return 0;
    }

    cudaMalloc((void **) &s_d_laplace_inv, sizeof(float) * (size_t) n_spec);
    cudaMalloc((void **) &s_d_exp_x, sizeof(cufftComplex) * (size_t) n_half);
    cudaMalloc((void **) &s_d_exp_y, sizeof(cufftComplex) * (size_t) pad);
    cudaMalloc((void **) &s_d_filt_a_pri, f_bytes);
    cudaMalloc((void **) &s_d_filt_b_pri, f_bytes);
    cudaMalloc((void **) &s_d_filt_a_sec, f_bytes);
    cudaMalloc((void **) &s_d_filt_b_sec, f_bytes);
    cudaMalloc((void **) &s_d_chrom_ramp, (size_t) (ns * n_spec) * sizeof(cufftComplex));

    if (params->use_moisan)
    {
        cudaMemcpy(s_d_laplace_inv, params->h_laplace_inv, sizeof(float) * (size_t) n_spec,
                   cudaMemcpyHostToDevice);
        cudaMemcpy(s_d_exp_x, params->h_exp_x, sizeof(cufftComplex) * (size_t) n_half,
                   cudaMemcpyHostToDevice);
        cudaMemcpy(s_d_exp_y, params->h_exp_y, sizeof(cufftComplex) * (size_t) pad,
                   cudaMemcpyHostToDevice);
    }

    for (int m = 0; m < ns; m++)
    {
        size_t off = (size_t) (m * n_spec);
        cudaMemcpy(s_d_filt_a_pri + off, params->h_filt_a_pri[m],
                   sizeof(float) * (size_t) n_spec, cudaMemcpyHostToDevice);
        cudaMemcpy(s_d_filt_b_pri + off, params->h_filt_b_pri[m],
                   sizeof(float) * (size_t) n_spec, cudaMemcpyHostToDevice);
        if (params->has_sec)
        {
            cudaMemcpy(s_d_filt_a_sec + off, params->h_filt_a_sec[m],
                       sizeof(float) * (size_t) n_spec, cudaMemcpyHostToDevice);
            cudaMemcpy(s_d_filt_b_sec + off, params->h_filt_b_sec[m],
                       sizeof(float) * (size_t) n_spec, cudaMemcpyHostToDevice);
            if (params->sec_shared && params->h_chrom_ramp != NULL &&
                params->h_chrom_ramp[m] != NULL)
            {
                cudaMemcpy(s_d_chrom_ramp + off, params->h_chrom_ramp[m],
                           sizeof(cufftComplex) * (size_t) n_spec, cudaMemcpyHostToDevice);
            }
        }
    }

    s_cached_filters = f_bytes;
    return 0;
}

/**
 * atmturb_cuda_rytov_init_workspace - Allocate device scratchpads for one frame
 * @params: Active simulation geometry and dimensions.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
static int atmturb_cuda_rytov_init_workspace(
    const atmturb_cuda_rytov_params_t *params)
{
    long pad_size   = params->pad_size;
    long pup_size   = params->pup_size;
    long pad_pixels = pad_size * pad_size;
    long pup_pixels = pup_size * pup_size;
    long n_half     = pad_size / 2 + 1;
    long n_spec     = pad_size * n_half;
    int  n_inv      = params->has_sec ? 4 : 2;
    size_t sub_bytes = (size_t) (params->nsuper * params->max_sublayers) *
                       sizeof(atmturb_cuda_rytov_sublayer_t);

    if (s_d_active_sublayers != NULL && s_cached_active_sub < sub_bytes)
    {
        cudaFree(s_d_active_sublayers);
        s_d_active_sublayers = NULL;
        s_cached_active_sub = 0;
    }
    if (s_d_active_sublayers == NULL)
    {
        if (cudaMalloc((void **) &s_d_active_sublayers, sub_bytes) != cudaSuccess)
        {
            return -1;
        }
        s_cached_active_sub = sub_bytes;
    }

    if (s_d_super_pha != NULL && s_cached_work == (size_t) pad_pixels)
    {
        return 0;
    }

    cudaMalloc((void **) &s_d_super_pha, sizeof(float) * (size_t) pad_pixels);
    cudaMalloc((void **) &s_d_super_spha, sizeof(float) * (size_t) pad_pixels);
    cudaMalloc((void **) &s_d_bounds, sizeof(float) * (size_t) (2 * pad_size));
    cudaMalloc((void **) &s_d_hat_bounds, sizeof(cufftComplex) * (size_t) (2 * n_half));
    cudaMalloc((void **) &s_d_spec, sizeof(cufftComplex) * (size_t) n_spec);
    cudaMalloc((void **) &s_d_acc_batch, sizeof(cufftComplex) * (size_t) (n_inv * n_spec));
    cudaMalloc((void **) &s_d_inv_out, sizeof(float) * (size_t) (n_inv * pad_pixels));
    cudaMalloc((void **) &s_d_frame_pha, sizeof(float) * (size_t) pup_pixels);
    cudaMalloc((void **) &s_d_frame_amp, sizeof(float) * (size_t) pup_pixels);
    cudaMalloc((void **) &s_d_frame_spha, sizeof(float) * (size_t) pup_pixels);
    cudaMalloc((void **) &s_d_frame_samp, sizeof(float) * (size_t) pup_pixels);
    cudaMalloc((void **) &s_d_sum_I, sizeof(float));

    s_cached_work = (size_t) pad_pixels;
    return 0;
}

/**
 * atmturb_cuda_rytov_step_ops - Enqueue GPU operations for one Rytov propagation frame
 * @params: Propagation parameters.
 * @d_masters: Device master screens buffer.
 * @stream: Stream on which kernels and FFTs are enqueued.
 */
static void atmturb_cuda_rytov_step_ops(
    const atmturb_cuda_rytov_params_t *params,
    const float                       *d_masters,
    cudaStream_t                       stream)
{
    long pad        = params->pad_size;
    long guard      = params->guard_pix;
    long pup_size   = params->pup_size;
    long pup_pixels = pup_size * pup_size;
    long n_half     = pad / 2 + 1;
    long n_spec     = pad * n_half;
    long pad_pixels = pad * pad;
    int  has_sec    = params->has_sec;
    int  n_inv      = has_sec ? 4 : 2;
    int  os         = (params->os > 1) ? params->os : 1;

    cufftComplex *d_acc_dphi_pri = s_d_acc_batch + 0 * n_spec;
    cufftComplex *d_acc_chi_pri  = s_d_acc_batch + 1 * n_spec;
    cufftComplex *d_acc_dphi_sec = s_d_acc_batch + 2 * n_spec;
    cufftComplex *d_acc_chi_sec  = s_d_acc_batch + 3 * n_spec;

    const float *d_dphi_out_pri = s_d_inv_out + 0 * pad_pixels;
    const float *d_chi_out_pri  = s_d_inv_out + 1 * pad_pixels;
    const float *d_dphi_out_sec = s_d_inv_out + 2 * pad_pixels;
    const float *d_chi_out_sec  = s_d_inv_out + 3 * pad_pixels;

    cudaMemsetAsync(s_d_frame_pha, 0, sizeof(float) * (size_t) pup_pixels, stream);
    cudaMemsetAsync(s_d_frame_spha, 0, sizeof(float) * (size_t) pup_pixels, stream);
    cudaMemsetAsync(s_d_acc_batch, 0, sizeof(cufftComplex) * (size_t) (n_inv * n_spec), stream);

    for (int m = 0; m < params->nsuper; m++)
    {
        int n_sub = params->super_nlayers[m];
        const atmturb_cuda_rytov_sublayer_t *sl_ptr =
            s_d_active_sublayers + (size_t) (m * params->max_sublayers);

        atmturb_cuda_extrude_sl_launch(
            params->interp, d_masters, (int) params->msize, (int) pad, (int) pup_size,
            (int) guard, n_sub, has_sec, os, sl_ptr, s_d_super_pha, s_d_super_spha,
            s_d_frame_pha, s_d_frame_spha, stream);

        if (params->super_dist_m[m] > 0.0)
        {
            atmturb_cuda_forward_spectrum(s_d_super_pha, params->use_moisan,
                                          (int) pad, (int) n_half, stream);
            size_t off = (size_t) (m * n_spec);

            if (has_sec && params->sec_shared && params->h_chrom_ramp != NULL &&
                params->h_chrom_ramp[m] != NULL)
            {
                atmturb_cuda_filter_mac_dual_launch(
                    d_acc_dphi_pri, d_acc_chi_pri, d_acc_dphi_sec, d_acc_chi_sec,
                    s_d_spec, s_d_chrom_ramp + off,
                    s_d_filt_a_pri + off, s_d_filt_b_pri + off,
                    s_d_filt_a_sec + off, s_d_filt_b_sec + off,
                    (int) n_spec, stream);
            }
            else
            {
                atmturb_cuda_filter_mac_launch(
                    d_acc_dphi_pri, d_acc_chi_pri, s_d_spec,
                    s_d_filt_a_pri + off, s_d_filt_b_pri + off, (int) n_spec, stream);
                if (has_sec)
                {
                    atmturb_cuda_forward_spectrum(s_d_super_spha, params->use_moisan,
                                                  (int) pad, (int) n_half, stream);
                    atmturb_cuda_filter_mac_launch(
                        d_acc_dphi_sec, d_acc_chi_sec, s_d_spec,
                        s_d_filt_a_sec + off, s_d_filt_b_sec + off, (int) n_spec, stream);
                }
            }
        }
    }

    cufftExecC2R(s_plan_2d_c2r_bat, (cufftComplex *) s_d_acc_batch, (cufftReal *) s_d_inv_out);

    cudaMemsetAsync(s_d_sum_I, 0, sizeof(float), stream);
    atmturb_cuda_assemble_launch((int) pup_size, (int) guard, (int) pad,
                                 s_d_frame_pha, s_d_frame_amp,
                                 d_dphi_out_pri, d_chi_out_pri,
                                 s_d_sum_I, stream);
    atmturb_cuda_normalize_launch(s_d_frame_amp, (int) pup_pixels, s_d_sum_I, stream);

    if (has_sec)
    {
        cudaMemsetAsync(s_d_sum_I, 0, sizeof(float), stream);
        atmturb_cuda_assemble_launch((int) pup_size, (int) guard, (int) pad,
                                     s_d_frame_spha, s_d_frame_samp,
                                     d_dphi_out_sec, d_chi_out_sec,
                                     s_d_sum_I, stream);
        atmturb_cuda_normalize_launch(s_d_frame_samp, (int) pup_pixels, s_d_sum_I, stream);
    }
}

/**
 * atmturb_cuda_rytov_init_graph - Capture or reuse CUDA Graph for Rytov frame pipeline
 * @params: Propagation parameters.
 * @d_masters: Device master screens buffer.
 *
 * Return: 0 on success, -1 on failure.
 */
static int atmturb_cuda_rytov_init_graph(
    const atmturb_cuda_rytov_params_t *params,
    const float                       *d_masters)
{
    long pad      = params->pad_size;
    long pup_size = params->pup_size;
    int  ns       = params->nsuper;
    int  has_sec  = params->has_sec;
    int  sec_sh   = params->sec_shared;
    int  interp   = params->interp;
    int  guard    = (int) params->guard_pix;
    int  moisan   = params->use_moisan;
    int  os       = (params->os > 1) ? params->os : 1;

    if (s_graph_exec != NULL &&
        s_graph_pad == pad &&
        s_graph_pup == pup_size &&
        s_graph_nsuper == ns &&
        s_graph_has_sec == has_sec &&
        s_graph_sec_shared == sec_sh &&
        s_graph_interp == interp &&
        s_graph_guard == guard &&
        s_graph_moisan == moisan &&
        s_graph_os == os &&
        s_graph_masters == d_masters)
    {
        return 0;
    }

    if (s_graph_exec != NULL)
    {
        cudaGraphExecDestroy(s_graph_exec);
        s_graph_exec = NULL;
    }

    cudaGraph_t graph;
    if (cudaStreamBeginCapture(s_stream, cudaStreamCaptureModeGlobal) != cudaSuccess)
    {
        return -1;
    }

    atmturb_cuda_rytov_step_ops(params, d_masters, s_stream);

    if (cudaStreamEndCapture(s_stream, &graph) != cudaSuccess)
    {
        return -1;
    }

    if (cudaGraphInstantiate(&s_graph_exec, graph, NULL, NULL, 0) != cudaSuccess)
    {
        cudaGraphDestroy(graph);
        s_graph_exec = NULL;
        return -1;
    }

    cudaGraphDestroy(graph);
    s_graph_pad        = pad;
    s_graph_pup        = pup_size;
    s_graph_nsuper     = ns;
    s_graph_has_sec    = has_sec;
    s_graph_sec_shared = sec_sh;
    s_graph_interp     = interp;
    s_graph_guard      = guard;
    s_graph_moisan     = moisan;
    s_graph_os         = os;
    s_graph_masters    = d_masters;
    return 0;
}

/**
 * atmturb_cuda_rytov_render - GPU-accelerated Rytov diffractive propagation
 * @params: Parameters, geometry, master screens, and destination buffers.
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_cuda_rytov_render(
    const atmturb_cuda_rytov_params_t *params)
{
    if (params->nblayers <= 0 || params->pup_size <= 0 || params->nbframes <= 0)
    {
        return -1;
    }

    const float *d_masters = atmturb_cuda_sync_device_masters(
        params->msize, params->nblayers, params->h_masters);
    if (d_masters == NULL)
    {
        return -1;
    }

    long pad        = params->pad_size;
    long pup_size   = params->pup_size;
    long pup_pixels = pup_size * pup_size;

    if (atmturb_cuda_rytov_init_plans(pad, params->has_sec) != 0 ||
        atmturb_cuda_rytov_init_filters(params) != 0 ||
        atmturb_cuda_rytov_init_workspace(params) != 0)
    {
        return -1;
    }

    atmturb_cuda_rytov_init_graph(params, d_masters);

    size_t sub_step_bytes = (size_t) (params->nsuper * params->max_sublayers) *
                             sizeof(atmturb_cuda_rytov_sublayer_t);

    if (params->nbframes == 1)
    {
        cudaMemcpyAsync(s_d_active_sublayers, params->frame_sublayers, sub_step_bytes,
                        cudaMemcpyHostToDevice, s_stream);

        if (s_graph_exec != NULL)
        {
            cudaGraphLaunch(s_graph_exec, s_stream);
        }
        else
        {
            atmturb_cuda_rytov_step_ops(params, d_masters, s_stream);
        }

        cudaMemcpyAsync(params->pha, s_d_frame_pha, sizeof(float) * (size_t) pup_pixels,
                        cudaMemcpyDeviceToHost, s_stream);
        cudaMemcpyAsync(params->amp, s_d_frame_amp, sizeof(float) * (size_t) pup_pixels,
                        cudaMemcpyDeviceToHost, s_stream);
        if (params->has_sec)
        {
            cudaMemcpyAsync(params->spha, s_d_frame_spha, sizeof(float) * (size_t) pup_pixels,
                            cudaMemcpyDeviceToHost, s_stream);
            cudaMemcpyAsync(params->samp, s_d_frame_samp, sizeof(float) * (size_t) pup_pixels,
                            cudaMemcpyDeviceToHost, s_stream);
        }
        cudaStreamSynchronize(s_stream);
        return 0;
    }

    size_t n_sub_tot = (size_t) (params->nbframes * params->nsuper * params->max_sublayers);
    atmturb_cuda_rytov_sublayer_t *d_all_sublayers = NULL;
    cudaMalloc((void **) &d_all_sublayers, n_sub_tot * sizeof(atmturb_cuda_rytov_sublayer_t));
    cudaMemcpyAsync(d_all_sublayers, params->frame_sublayers,
                    n_sub_tot * sizeof(atmturb_cuda_rytov_sublayer_t),
                    cudaMemcpyHostToDevice, s_stream);

    for (long t = 0; t < params->nbframes; t++)
    {
        const atmturb_cuda_rytov_sublayer_t *src_t =
            d_all_sublayers + (size_t) t * (params->nsuper * params->max_sublayers);
        cudaMemcpyAsync(s_d_active_sublayers, src_t, sub_step_bytes,
                        cudaMemcpyDeviceToDevice, s_stream);

        if (s_graph_exec != NULL)
        {
            cudaGraphLaunch(s_graph_exec, s_stream);
        }
        else
        {
            atmturb_cuda_rytov_step_ops(params, d_masters, s_stream);
        }

        size_t frame_off = (size_t) (t * pup_pixels);
        cudaMemcpyAsync(params->pha + frame_off, s_d_frame_pha,
                        sizeof(float) * (size_t) pup_pixels,
                        cudaMemcpyDeviceToHost, s_stream);
        cudaMemcpyAsync(params->amp + frame_off, s_d_frame_amp,
                        sizeof(float) * (size_t) pup_pixels,
                        cudaMemcpyDeviceToHost, s_stream);

        if (params->has_sec)
        {
            cudaMemcpyAsync(params->spha + frame_off, s_d_frame_spha,
                            sizeof(float) * (size_t) pup_pixels,
                            cudaMemcpyDeviceToHost, s_stream);
            cudaMemcpyAsync(params->samp + frame_off, s_d_frame_samp,
                            sizeof(float) * (size_t) pup_pixels,
                            cudaMemcpyDeviceToHost, s_stream);
        }
    }

    cudaStreamSynchronize(s_stream);
    cudaFree(d_all_sublayers);
    return 0;
}
