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

// Persistent CUDA stream and Graph execution instances
static cudaStream_t                   s_stream                   = 0;
static cudaGraphExec_t                s_graph_exec[2]            = {NULL, NULL};
static long                           s_graph_pad                = 0;
static long                           s_graph_pup                = 0;
static int                            s_graph_nsuper             = 0;
static int                            s_graph_has_sec            = 0;
static int                            s_graph_sec_shared         = 0;
static int                            s_graph_interp             = 0;
static int                            s_graph_guard              = 0;
static int                            s_graph_moisan             = -1;
static int                            s_graph_os                 = 0;
static int                            s_graph_lowfreq            = -1;
static int                            s_graph_piston             = -1;
static const float                   *s_graph_masters            = NULL;
static atmturb_cuda_rytov_sublayer_t *s_d_active_sublayers[2]    = {NULL, NULL};
static size_t                         s_cached_active_sub        = 0;
static atmturb_cuda_lowfreq_mode_t   *s_d_active_lf_modes[2]     = {NULL, NULL};
static size_t                         s_cached_active_lf         = 0;

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
static float        *s_d_frame_pha[2]  = {NULL, NULL};
static float        *s_d_frame_amp[2]  = {NULL, NULL};
static float        *s_d_frame_spha[2] = {NULL, NULL};
static float        *s_d_frame_samp[2] = {NULL, NULL};
static float        *s_d_sum_I[2]      = {NULL, NULL};
static size_t        s_cached_work    = 0;

/**
 * atmturb_cuda_rytov_cleanup - Free cached GPU buffers, plans, and stream
 */
void atmturb_cuda_rytov_cleanup(void)
{
    for (int s = 0; s < 2; s++)
    {
        if (s_graph_exec[s] != NULL)
        {
            cudaGraphExecDestroy(s_graph_exec[s]);
            s_graph_exec[s] = NULL;
        }
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
    s_graph_lowfreq     = -1;
    s_graph_piston      = -1;
    s_graph_masters     = NULL;

    for (int s = 0; s < 2; s++)
    {
        if (s_d_active_sublayers[s] != NULL)
        {
            cudaFree(s_d_active_sublayers[s]);
            s_d_active_sublayers[s] = NULL;
        }
        if (s_d_active_lf_modes[s] != NULL)
        {
            cudaFree(s_d_active_lf_modes[s]);
            s_d_active_lf_modes[s] = NULL;
        }
    }
    s_cached_active_sub = 0;
    s_cached_active_lf  = 0;

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
    for (int s = 0; s < 2; s++)
    {
        if (s_d_frame_pha[s] != NULL)  { cudaFree(s_d_frame_pha[s]);  s_d_frame_pha[s] = NULL; }
        if (s_d_frame_amp[s] != NULL)  { cudaFree(s_d_frame_amp[s]);  s_d_frame_amp[s] = NULL; }
        if (s_d_frame_spha[s] != NULL) { cudaFree(s_d_frame_spha[s]); s_d_frame_spha[s] = NULL; }
        if (s_d_frame_samp[s] != NULL) { cudaFree(s_d_frame_samp[s]); s_d_frame_samp[s] = NULL; }
        if (s_d_sum_I[s] != NULL)      { cudaFree(s_d_sum_I[s]);      s_d_sum_I[s] = NULL; }
    }
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

    for (int s = 0; s < 2; s++)
    {
        if (s_d_active_sublayers[s] != NULL && s_cached_active_sub < sub_bytes)
        {
            cudaFree(s_d_active_sublayers[s]);
            s_d_active_sublayers[s] = NULL;
        }
        if (s_d_active_sublayers[s] == NULL)
        {
            if (cudaMalloc((void **) &s_d_active_sublayers[s], sub_bytes) != cudaSuccess)
            {
                return -1;
            }
        }
    }
    s_cached_active_sub = sub_bytes;

    if (params->lowfreq && params->max_lf_modes > 0)
    {
        size_t lf_bytes = (size_t) (params->nsuper * params->max_lf_modes) *
                          sizeof(atmturb_cuda_lowfreq_mode_t);
        for (int s = 0; s < 2; s++)
        {
            if (s_d_active_lf_modes[s] != NULL && s_cached_active_lf < lf_bytes)
            {
                cudaFree(s_d_active_lf_modes[s]);
                s_d_active_lf_modes[s] = NULL;
            }
            if (s_d_active_lf_modes[s] == NULL)
            {
                if (cudaMalloc((void **) &s_d_active_lf_modes[s], lf_bytes) != cudaSuccess)
                {
                    return -1;
                }
            }
        }
        s_cached_active_lf = lf_bytes;
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
    for (int s = 0; s < 2; s++)
    {
        cudaMalloc((void **) &s_d_frame_pha[s], sizeof(float) * (size_t) pup_pixels);
        cudaMalloc((void **) &s_d_frame_amp[s], sizeof(float) * (size_t) pup_pixels);
        cudaMalloc((void **) &s_d_frame_spha[s], sizeof(float) * (size_t) pup_pixels);
        cudaMalloc((void **) &s_d_frame_samp[s], sizeof(float) * (size_t) pup_pixels);
        cudaMalloc((void **) &s_d_sum_I[s], sizeof(float));
    }

    s_cached_work = (size_t) pad_pixels;
    return 0;
}

/**
 * atmturb_cuda_rytov_accum_spectrum - Accumulate optical field spectrum for one superlayer
 * @params: Propagation parameters.
 * @m: Superlayer index.
 * @stream: CUDA execution stream.
 */
static void atmturb_cuda_rytov_accum_spectrum(
    const atmturb_cuda_rytov_params_t *params,
    int                                m,
    cudaStream_t                       stream)
{
    long pad        = params->pad_size;
    long n_half     = pad / 2 + 1;
    long n_spec     = pad * n_half;
    int  has_sec    = params->has_sec;
    size_t off      = (size_t) (m * n_spec);

    cufftComplex *d_acc_dphi_pri = s_d_acc_batch + 0 * n_spec;
    cufftComplex *d_acc_chi_pri  = s_d_acc_batch + 1 * n_spec;
    cufftComplex *d_acc_dphi_sec = s_d_acc_batch + 2 * n_spec;
    cufftComplex *d_acc_chi_sec  = s_d_acc_batch + 3 * n_spec;

    int is_dual = (has_sec && params->sec_shared && params->h_chrom_ramp != NULL &&
                   params->h_chrom_ramp[m] != NULL);

    if (params->use_moisan)
    {
        atmturb_cuda_extract_bounds_launch(s_d_super_pha, s_d_bounds, (int) pad, stream);
        cufftExecR2C(s_plan_1d, (cufftReal *) s_d_bounds, (cufftComplex *) s_d_hat_bounds);
        cufftExecR2C(s_plan_2d_r2c, (cufftReal *) s_d_super_pha, (cufftComplex *) s_d_spec);

        if (is_dual)
        {
            atmturb_cuda_moisan_dual_mac_t args;
            args.d_acc_dphi_pri = d_acc_dphi_pri;
            args.d_acc_chi_pri  = d_acc_chi_pri;
            args.d_acc_dphi_sec = d_acc_dphi_sec;
            args.d_acc_chi_sec  = d_acc_chi_sec;
            args.d_spec         = s_d_spec;
            args.d_ramp         = s_d_chrom_ramp + off;
            args.d_hat_bounds   = s_d_hat_bounds;
            args.d_exp_x        = s_d_exp_x;
            args.d_exp_y        = s_d_exp_y;
            args.d_laplace_inv  = s_d_laplace_inv;
            args.d_filt_a_pri   = s_d_filt_a_pri + off;
            args.d_filt_b_pri   = s_d_filt_b_pri + off;
            args.d_filt_a_sec   = s_d_filt_a_sec + off;
            args.d_filt_b_sec   = s_d_filt_b_sec + off;
            args.pad_size       = (int) pad;
            args.n_half         = (int) n_half;
            atmturb_cuda_moisan_dual_mac_launch(&args, stream);
        }
        else
        {
            atmturb_cuda_moisan_mac_t args;
            args.d_acc_dphi    = d_acc_dphi_pri;
            args.d_acc_chi     = d_acc_chi_pri;
            args.d_spec        = s_d_spec;
            args.d_hat_bounds  = s_d_hat_bounds;
            args.d_exp_x       = s_d_exp_x;
            args.d_exp_y       = s_d_exp_y;
            args.d_laplace_inv = s_d_laplace_inv;
            args.d_filt_a      = s_d_filt_a_pri + off;
            args.d_filt_b      = s_d_filt_b_pri + off;
            args.pad_size      = (int) pad;
            args.n_half        = (int) n_half;
            atmturb_cuda_moisan_filter_mac_launch(&args, stream);

            if (has_sec)
            {
                atmturb_cuda_extract_bounds_launch(s_d_super_spha, s_d_bounds, (int) pad, stream);
                cufftExecR2C(s_plan_1d, (cufftReal *) s_d_bounds, (cufftComplex *) s_d_hat_bounds);
                cufftExecR2C(s_plan_2d_r2c, (cufftReal *) s_d_super_spha, (cufftComplex *) s_d_spec);
                args.d_acc_dphi = d_acc_dphi_sec;
                args.d_acc_chi  = d_acc_chi_sec;
                args.d_filt_a   = s_d_filt_a_sec + off;
                args.d_filt_b   = s_d_filt_b_sec + off;
                atmturb_cuda_moisan_filter_mac_launch(&args, stream);
            }
        }
    }
    else
    {
        cufftExecR2C(s_plan_2d_r2c, (cufftReal *) s_d_super_pha, (cufftComplex *) s_d_spec);
        if (is_dual)
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
                cufftExecR2C(s_plan_2d_r2c, (cufftReal *) s_d_super_spha,
                             (cufftComplex *) s_d_spec);
                atmturb_cuda_filter_mac_launch(
                    d_acc_dphi_sec, d_acc_chi_sec, s_d_spec,
                    s_d_filt_a_sec + off, s_d_filt_b_sec + off, (int) n_spec, stream);
            }
        }
    }
}

/**
 * atmturb_cuda_rytov_step_ops - Enqueue GPU operations for one Rytov propagation frame
 * @params: Propagation parameters.
 * @d_masters: Device master screens buffer.
 * @slot: Active double-buffer slot index (0 or 1).
 * @stream: Stream on which kernels and FFTs are enqueued.
 */
static void atmturb_cuda_rytov_step_ops(
    const atmturb_cuda_rytov_params_t *params,
    const float                       *d_masters,
    int                                slot,
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

    const float *d_dphi_out_pri = s_d_inv_out + 0 * pad_pixels;
    const float *d_chi_out_pri  = s_d_inv_out + 1 * pad_pixels;
    const float *d_dphi_out_sec = s_d_inv_out + 2 * pad_pixels;
    const float *d_chi_out_sec  = s_d_inv_out + 3 * pad_pixels;

    cudaMemsetAsync(s_d_frame_pha[slot], 0, sizeof(float) * (size_t) pup_pixels, stream);
    cudaMemsetAsync(s_d_frame_spha[slot], 0, sizeof(float) * (size_t) pup_pixels, stream);
    cudaMemsetAsync(s_d_acc_batch, 0, sizeof(cufftComplex) * (size_t) (n_inv * n_spec), stream);

    for (int m = 0; m < params->nsuper; m++)
    {
        int n_sub = params->super_nlayers[m];
        const atmturb_cuda_rytov_sublayer_t *sl_ptr =
            s_d_active_sublayers[slot] + (size_t) (m * params->max_sublayers);

        int n_lf = (params->lowfreq && params->super_nlf_modes != NULL &&
                    params->super_nlf_modes[m] > 0 && s_d_active_lf_modes[slot] != NULL)
                   ? params->super_nlf_modes[m] : 0;
        const atmturb_cuda_lowfreq_mode_t *lf_ptr = (n_lf > 0)
                   ? s_d_active_lf_modes[slot] + (size_t) (m * params->max_lf_modes) : NULL;

        atmturb_cuda_extrude_sl_launch(
            params->interp, d_masters, (int) params->msize, (int) pad, (int) pup_size,
            (int) guard, n_sub, has_sec, os, sl_ptr, n_lf, lf_ptr,
            s_d_super_pha, s_d_super_spha, s_d_frame_pha[slot], s_d_frame_spha[slot], stream);

        if (params->super_dist_m[m] > 0.0)
        {
            atmturb_cuda_rytov_accum_spectrum(params, m, stream);
        }
    }

    cufftExecC2R(s_plan_2d_c2r_bat, (cufftComplex *) s_d_acc_batch, (cufftReal *) s_d_inv_out);

    cudaMemsetAsync(s_d_sum_I[slot], 0, sizeof(float), stream);
    atmturb_cuda_assemble_launch((int) pup_size, (int) guard, (int) pad,
                                 s_d_frame_pha[slot], s_d_frame_amp[slot],
                                 d_dphi_out_pri, d_chi_out_pri,
                                 s_d_sum_I[slot], stream);
    atmturb_cuda_normalize_launch(s_d_frame_amp[slot], (int) pup_pixels, s_d_sum_I[slot], stream);

    if (has_sec)
    {
        cudaMemsetAsync(s_d_sum_I[slot], 0, sizeof(float), stream);
        atmturb_cuda_assemble_launch((int) pup_size, (int) guard, (int) pad,
                                     s_d_frame_spha[slot], s_d_frame_samp[slot],
                                     d_dphi_out_sec, d_chi_out_sec,
                                     s_d_sum_I[slot], stream);
        atmturb_cuda_normalize_launch(s_d_frame_samp[slot], (int) pup_pixels, s_d_sum_I[slot], stream);
    }

    if (params->remove_piston)
    {
        atmturb_cuda_remove_piston_launch(
            s_d_frame_pha[slot], s_d_sum_I[slot], (int) pup_pixels, stream);
        if (has_sec)
        {
            atmturb_cuda_remove_piston_launch(
                s_d_frame_spha[slot], s_d_sum_I[slot], (int) pup_pixels, stream);
        }
    }
}

/**
 * atmturb_cuda_rytov_init_graph - Capture or reuse CUDA Graph for Rytov frame pipeline
 * @params: Propagation parameters.
 * @d_masters: Device master screens buffer.
 * @slot: Active double-buffer slot index (0 or 1).
 *
 * Return: 0 on success, -1 on failure.
 */
static int atmturb_cuda_rytov_init_graph(
    const atmturb_cuda_rytov_params_t *params,
    const float                       *d_masters,
    int                                slot)
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

    if (s_graph_exec[slot] != NULL &&
        s_graph_pad == pad &&
        s_graph_pup == pup_size &&
        s_graph_nsuper == ns &&
        s_graph_has_sec == has_sec &&
        s_graph_sec_shared == sec_sh &&
        s_graph_interp == interp &&
        s_graph_guard == guard &&
        s_graph_moisan == moisan &&
        s_graph_os == os &&
        s_graph_lowfreq == params->lowfreq &&
        s_graph_piston == params->remove_piston &&
        s_graph_masters == d_masters)
    {
        return 0;
    }

    if (s_graph_exec[slot] != NULL)
    {
        cudaGraphExecDestroy(s_graph_exec[slot]);
        s_graph_exec[slot] = NULL;
    }

    cudaGraph_t graph;
    if (cudaStreamBeginCapture(s_stream, cudaStreamCaptureModeGlobal) != cudaSuccess)
    {
        return -1;
    }

    atmturb_cuda_rytov_step_ops(params, d_masters, slot, s_stream);

    if (cudaStreamEndCapture(s_stream, &graph) != cudaSuccess)
    {
        return -1;
    }

    if (cudaGraphInstantiate(&s_graph_exec[slot], graph, NULL, NULL, 0) != cudaSuccess)
    {
        cudaGraphDestroy(graph);
        s_graph_exec[slot] = NULL;
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
    s_graph_lowfreq    = params->lowfreq;
    s_graph_piston     = params->remove_piston;
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

    int slot = params->slot & 1;
    atmturb_cuda_rytov_init_graph(params, d_masters, slot);

    size_t sub_step_bytes = (size_t) (params->nsuper * params->max_sublayers) *
                             sizeof(atmturb_cuda_rytov_sublayer_t);

    if (params->nbframes == 1)
    {
        cudaMemcpyAsync(s_d_active_sublayers[slot], params->frame_sublayers, sub_step_bytes,
                        cudaMemcpyHostToDevice, s_stream);

        if (params->lowfreq && params->frame_lf_modes != NULL && s_d_active_lf_modes[slot] != NULL)
        {
            size_t lf_step_bytes = (size_t) (params->nsuper * params->max_lf_modes) *
                                   sizeof(atmturb_cuda_lowfreq_mode_t);
            cudaMemcpyAsync(s_d_active_lf_modes[slot], params->frame_lf_modes, lf_step_bytes,
                            cudaMemcpyHostToDevice, s_stream);
        }

        if (s_graph_exec[slot] != NULL)
        {
            cudaGraphLaunch(s_graph_exec[slot], s_stream);
        }
        else
        {
            atmturb_cuda_rytov_step_ops(params, d_masters, slot, s_stream);
        }

        if (params->pha != NULL)
        {
            cudaMemcpyAsync(params->pha, s_d_frame_pha[slot], sizeof(float) * (size_t) pup_pixels,
                            cudaMemcpyDeviceToHost, s_stream);
        }
        if (params->amp != NULL)
        {
            cudaMemcpyAsync(params->amp, s_d_frame_amp[slot], sizeof(float) * (size_t) pup_pixels,
                            cudaMemcpyDeviceToHost, s_stream);
        }
        if (params->has_sec && params->spha != NULL)
        {
            cudaMemcpyAsync(params->spha, s_d_frame_spha[slot], sizeof(float) * (size_t) pup_pixels,
                            cudaMemcpyDeviceToHost, s_stream);
        }
        if (params->has_sec && params->samp != NULL)
        {
            cudaMemcpyAsync(params->samp, s_d_frame_samp[slot], sizeof(float) * (size_t) pup_pixels,
                            cudaMemcpyDeviceToHost, s_stream);
        }
        if (params->async_stream)
        {
            return 0;
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

    atmturb_cuda_lowfreq_mode_t *d_all_lf_modes = NULL;
    size_t lf_step_bytes = 0;
    if (params->lowfreq && params->frame_lf_modes != NULL && s_d_active_lf_modes != NULL)
    {
        lf_step_bytes = (size_t) (params->nsuper * params->max_lf_modes) *
                        sizeof(atmturb_cuda_lowfreq_mode_t);
        size_t n_lf_tot = (size_t) (params->nbframes * params->nsuper * params->max_lf_modes);
        cudaMalloc((void **) &d_all_lf_modes, n_lf_tot * sizeof(atmturb_cuda_lowfreq_mode_t));
        cudaMemcpyAsync(d_all_lf_modes, params->frame_lf_modes,
                        n_lf_tot * sizeof(atmturb_cuda_lowfreq_mode_t),
                        cudaMemcpyHostToDevice, s_stream);
    }

    for (long t = 0; t < params->nbframes; t++)
    {
        int bslot = (int) (t & 1);
        const atmturb_cuda_rytov_sublayer_t *src_t =
            d_all_sublayers + (size_t) t * (params->nsuper * params->max_sublayers);
        cudaMemcpyAsync(s_d_active_sublayers[bslot], src_t, sub_step_bytes,
                        cudaMemcpyDeviceToDevice, s_stream);

        if (d_all_lf_modes != NULL)
        {
            const atmturb_cuda_lowfreq_mode_t *src_lf_t =
                d_all_lf_modes + (size_t) t * (params->nsuper * params->max_lf_modes);
            cudaMemcpyAsync(s_d_active_lf_modes[bslot], src_lf_t, lf_step_bytes,
                            cudaMemcpyDeviceToDevice, s_stream);
        }

        if (s_graph_exec[bslot] != NULL)
        {
            cudaGraphLaunch(s_graph_exec[bslot], s_stream);
        }
        else
        {
            atmturb_cuda_rytov_step_ops(params, d_masters, bslot, s_stream);
        }

        size_t frame_off = (size_t) (t * pup_pixels);
        if (params->pha != NULL)
        {
            cudaMemcpyAsync(params->pha + frame_off, s_d_frame_pha[bslot],
                            sizeof(float) * (size_t) pup_pixels,
                            cudaMemcpyDeviceToHost, s_stream);
        }
        if (params->amp != NULL)
        {
            cudaMemcpyAsync(params->amp + frame_off, s_d_frame_amp[bslot],
                            sizeof(float) * (size_t) pup_pixels,
                            cudaMemcpyDeviceToHost, s_stream);
        }

        if (params->has_sec && params->spha != NULL)
        {
            cudaMemcpyAsync(params->spha + frame_off, s_d_frame_spha[bslot],
                            sizeof(float) * (size_t) pup_pixels,
                            cudaMemcpyDeviceToHost, s_stream);
        }
        if (params->has_sec && params->samp != NULL)
        {
            cudaMemcpyAsync(params->samp + frame_off, s_d_frame_samp[bslot],
                            sizeof(float) * (size_t) pup_pixels,
                            cudaMemcpyDeviceToHost, s_stream);
        }
    }

    cudaStreamSynchronize(s_stream);
    cudaFree(d_all_sublayers);
    if (d_all_lf_modes != NULL)
    {
        cudaFree(d_all_lf_modes);
    }
    return 0;
}

int atmturb_cuda_rytov_create_event(
    void **ev)
{
    if (ev == NULL)
    {
        return -1;
    }
    cudaEvent_t e;
    if (cudaEventCreateWithFlags(&e, cudaEventDisableTiming) != cudaSuccess)
    {
        *ev = NULL;
        return -1;
    }
    *ev = (void *) e;
    return 0;
}

void atmturb_cuda_rytov_destroy_event(
    void *ev)
{
    if (ev != NULL)
    {
        cudaEventDestroy((cudaEvent_t) ev);
    }
}

int atmturb_cuda_rytov_record_event(
    void *ev)
{
    if (ev == NULL || s_stream == 0)
    {
        return -1;
    }
    return (cudaEventRecord((cudaEvent_t) ev, s_stream) == cudaSuccess) ? 0 : -1;
}

int atmturb_cuda_rytov_wait_event(
    void *ev)
{
    if (ev == NULL)
    {
        return -1;
    }
    return (cudaEventSynchronize((cudaEvent_t) ev) == cudaSuccess) ? 0 : -1;
}
