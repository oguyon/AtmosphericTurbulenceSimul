// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    benchmark_perf.c
 * @brief   Comprehensive performance benchmark comparing SIMD, OpenMP, and CUDA
 */

#define _GNU_SOURCE
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>

#include "atmturb_lowfreq.h"
#include "atmturb_simd.h"
#include "atmturb_cuda.h"
#include "wfprop_fresnel_cuda.h"

#ifndef PI
#    define PI 3.14159265358979323846264338328
#endif

typedef struct
{
    float re;
    float im;
} complex_float;

/**
 * get_time_sec - High-precision monotonic timestamp in seconds
 *
 * Return: Current monotonic time in seconds.
 */
static double get_time_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double) ts.tv_sec + (double) ts.tv_nsec * 1e-9;
}

/**
 * run_extrusion_scalar - Execute wavefront extrusion with scalar fallback
 * @masters: Pointers to master phase screens.
 * @msize: Master screen linear dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 * @nbframes: Number of simulation time steps.
 * @nblayers: Number of turbulence layers.
 * @vx: X-velocity per layer [pix/frame].
 * @vy: Y-velocity per layer [pix/frame].
 * @cn2: Relative refractive index structure constant per layer.
 * @out_pha: Destination phase buffer.
 */
static void run_extrusion_scalar(
    const float *const *masters,
    long                msize,
    long                pup_size,
    long                nbframes,
    long                nblayers,
    const double       *vx,
    const double       *vy,
    const double       *cn2,
    float              *out_pha)
{
    long frame_pixels = pup_size * pup_size;
    for (long t = 0; t < nbframes; t++)
    {
        float *pha_slice = &out_pha[t * frame_pixels];
        for (long i = 0; i < frame_pixels; i++)
        {
            pha_slice[i] = 0.0f;
        }

        for (long k = 0; k < nblayers; k++)
        {
            double cur_x = (double) (t + 1) * vx[k];
            double cur_y = (double) (t + 1) * vy[k];
            float weight = (float) sqrt(cn2[k]);
            atmturb_extrude_params_t ep = {
                .master   = masters[k],
                .msize    = msize,
                .x0       = cur_x,
                .y0       = cur_y,
                .pup_size = pup_size,
                .os       = 1,
                .interp   = ATMTURB_INTERP_BILINEAR,
                .weight   = weight,
                .out_pha  = pha_slice
            };
            atmturb_extrude_accumulate_scalar(&ep);
        }
    }
}

/**
 * run_extrusion_avx2_single - Execute wavefront extrusion with single-threaded AVX2
 * @masters: Pointers to master phase screens.
 * @msize: Master screen linear dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 * @nbframes: Number of simulation time steps.
 * @nblayers: Number of turbulence layers.
 * @vx: X-velocity per layer [pix/frame].
 * @vy: Y-velocity per layer [pix/frame].
 * @cn2: Relative refractive index structure constant per layer.
 * @out_pha: Destination phase buffer.
 */
static void run_extrusion_avx2_single(
    const float *const *masters,
    long                msize,
    long                pup_size,
    long                nbframes,
    long                nblayers,
    const double       *vx,
    const double       *vy,
    const double       *cn2,
    float              *out_pha)
{
    long frame_pixels = pup_size * pup_size;
    for (long t = 0; t < nbframes; t++)
    {
        float *pha_slice = &out_pha[t * frame_pixels];
        for (long i = 0; i < frame_pixels; i++)
        {
            pha_slice[i] = 0.0f;
        }

        for (long k = 0; k < nblayers; k++)
        {
            double cur_x = (double) (t + 1) * vx[k];
            double cur_y = (double) (t + 1) * vy[k];
            float weight = (float) sqrt(cn2[k]);
            atmturb_extrude_params_t ep = {
                .master   = masters[k],
                .msize    = msize,
                .x0       = cur_x,
                .y0       = cur_y,
                .pup_size = pup_size,
                .os       = 1,
                .interp   = ATMTURB_INTERP_BILINEAR,
                .weight   = weight,
                .out_pha  = pha_slice
            };
            atmturb_extrude_accumulate_avx2(&ep);
        }
    }
}

/**
 * run_extrusion_avx2_omp - Execute wavefront extrusion with AVX2 and OpenMP threading
 * @masters: Pointers to master phase screens.
 * @msize: Master screen linear dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 * @nbframes: Number of simulation time steps.
 * @nblayers: Number of turbulence layers.
 * @vx: X-velocity per layer [pix/frame].
 * @vy: Y-velocity per layer [pix/frame].
 * @cn2: Relative refractive index structure constant per layer.
 * @out_pha: Destination phase buffer.
 */
static void run_extrusion_avx2_omp(
    const float *const *masters,
    long                msize,
    long                pup_size,
    long                nbframes,
    long                nblayers,
    const double       *vx,
    const double       *vy,
    const double       *cn2,
    float              *out_pha)
{
    long frame_pixels = pup_size * pup_size;
    #pragma omp parallel for schedule(dynamic)
    for (long t = 0; t < nbframes; t++)
    {
        float *pha_slice = &out_pha[t * frame_pixels];
        for (long i = 0; i < frame_pixels; i++)
        {
            pha_slice[i] = 0.0f;
        }

        for (long k = 0; k < nblayers; k++)
        {
            double cur_x = (double) (t + 1) * vx[k];
            double cur_y = (double) (t + 1) * vy[k];
            float weight = (float) sqrt(cn2[k]);
            atmturb_extrude_params_t ep = {
                .master   = masters[k],
                .msize    = msize,
                .x0       = cur_x,
                .y0       = cur_y,
                .pup_size = pup_size,
                .os       = 1,
                .interp   = ATMTURB_INTERP_BILINEAR,
                .weight   = weight,
                .out_pha  = pha_slice
            };
            atmturb_extrude_accumulate_avx2(&ep);
        }
    }
}

/**
 * run_extrusion_avx512_single - Execute wavefront extrusion with single-threaded AVX-512
 * @masters: Pointers to master phase screens.
 * @msize: Master screen linear dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 * @nbframes: Number of simulation time steps.
 * @nblayers: Number of turbulence layers.
 * @vx: X-velocity per layer [pix/frame].
 * @vy: Y-velocity per layer [pix/frame].
 * @cn2: Relative refractive index structure constant per layer.
 * @out_pha: Destination phase buffer.
 */
static void run_extrusion_avx512_single(
    const float *const *masters,
    long                msize,
    long                pup_size,
    long                nbframes,
    long                nblayers,
    const double       *vx,
    const double       *vy,
    const double       *cn2,
    float              *out_pha)
{
    long frame_pixels = pup_size * pup_size;
    for (long t = 0; t < nbframes; t++)
    {
        float *pha_slice = &out_pha[t * frame_pixels];
        for (long i = 0; i < frame_pixels; i++)
        {
            pha_slice[i] = 0.0f;
        }

        for (long k = 0; k < nblayers; k++)
        {
            double cur_x = (double) (t + 1) * vx[k];
            double cur_y = (double) (t + 1) * vy[k];
            float weight = (float) sqrt(cn2[k]);
            atmturb_extrude_params_t ep = {
                .master   = masters[k],
                .msize    = msize,
                .x0       = cur_x,
                .y0       = cur_y,
                .pup_size = pup_size,
                .os       = 1,
                .interp   = ATMTURB_INTERP_BILINEAR,
                .weight   = weight,
                .out_pha  = pha_slice
            };
            atmturb_extrude_accumulate_avx512(&ep);
        }
    }
}

/**
 * run_extrusion_avx512_omp - Execute wavefront extrusion with AVX-512 and OpenMP threading
 * @masters: Pointers to master phase screens.
 * @msize: Master screen linear dimension in pixels.
 * @pup_size: Output pupil dimension in pixels.
 * @nbframes: Number of simulation time steps.
 * @nblayers: Number of turbulence layers.
 * @vx: X-velocity per layer [pix/frame].
 * @vy: Y-velocity per layer [pix/frame].
 * @cn2: Relative refractive index structure constant per layer.
 * @out_pha: Destination phase buffer.
 */
static void run_extrusion_avx512_omp(
    const float *const *masters,
    long                msize,
    long                pup_size,
    long                nbframes,
    long                nblayers,
    const double       *vx,
    const double       *vy,
    const double       *cn2,
    float              *out_pha)
{
    long frame_pixels = pup_size * pup_size;
    #pragma omp parallel for schedule(dynamic)
    for (long t = 0; t < nbframes; t++)
    {
        float *pha_slice = &out_pha[t * frame_pixels];
        for (long i = 0; i < frame_pixels; i++)
        {
            pha_slice[i] = 0.0f;
        }

        for (long k = 0; k < nblayers; k++)
        {
            double cur_x = (double) (t + 1) * vx[k];
            double cur_y = (double) (t + 1) * vy[k];
            float weight = (float) sqrt(cn2[k]);
            atmturb_extrude_params_t ep = {
                .master   = masters[k],
                .msize    = msize,
                .x0       = cur_x,
                .y0       = cur_y,
                .pup_size = pup_size,
                .os       = 1,
                .interp   = ATMTURB_INTERP_BILINEAR,
                .weight   = weight,
                .out_pha  = pha_slice
            };
            atmturb_extrude_accumulate_avx512(&ep);
        }
    }
}

/**
 * benchmark_extrusion - Benchmark wavefront extrusion across compute engines
 * @pup_size: Linear dimension of pupil in pixels.
 * @nbframes: Number of simulation time steps.
 * @nblayers: Number of turbulent layers.
 */
static void benchmark_extrusion(
    long pup_size,
    long nbframes,
    long nblayers)
{
    long msize = 2048;
    printf("\n========================================================================\n");
    printf(" Benchmark 1: Wavefront Extrusion (%ldx%ld pupil, %ld frames, %ld layers)\n",
           pup_size, pup_size, nbframes, nblayers);
    printf("========================================================================\n");

    float **masters = (float **) malloc(sizeof(float *) * (size_t) nblayers);
    for (long k = 0; k < nblayers; k++)
    {
        masters[k] = (float *) malloc(sizeof(float) * (size_t) (msize * msize));
        for (long i = 0; i < msize * msize; i++)
        {
            masters[k][i] = (float) rand() / (float) RAND_MAX;
        }
    }

    double *vx  = (double *) malloc(sizeof(double) * (size_t) nblayers);
    double *vy  = (double *) malloc(sizeof(double) * (size_t) nblayers);
    double *cn2 = (double *) malloc(sizeof(double) * (size_t) nblayers);
    for (long k = 0; k < nblayers; k++)
    {
        vx[k]  = 0.5 + 0.2 * (double) k;
        vy[k]  = 0.3 - 0.1 * (double) k;
        cn2[k] = 1.0 / (double) (k + 1);
    }

    long total_pixels = nbframes * pup_size * pup_size;
    float *out_pha  = (float *) malloc(sizeof(float) * (size_t) total_pixels);
    float *out_amp  = (float *) malloc(sizeof(float) * (size_t) total_pixels);
    float *out_spha = (float *) malloc(sizeof(float) * (size_t) total_pixels);
    float *out_samp = (float *) malloc(sizeof(float) * (size_t) total_pixels);

    // Warm up
    run_extrusion_scalar((const float *const *) masters, msize, pup_size, 2,
                         nblayers, vx, vy, cn2, out_pha);

    // 1. Scalar (1 thread)
    double t0 = get_time_sec();
    run_extrusion_scalar((const float *const *) masters, msize, pup_size,
                         nbframes, nblayers, vx, vy, cn2, out_pha);
    double t_scalar = get_time_sec() - t0;

    // 2. AVX2 (1 thread)
    t0 = get_time_sec();
    run_extrusion_avx2_single((const float *const *) masters, msize, pup_size,
                              nbframes, nblayers, vx, vy, cn2, out_pha);
    double t_avx2_single = get_time_sec() - t0;

    // 3. AVX2 + OpenMP (multi-threaded)
    t0 = get_time_sec();
    run_extrusion_avx2_omp((const float *const *) masters, msize, pup_size,
                           nbframes, nblayers, vx, vy, cn2, out_pha);
    double t_avx2_omp = get_time_sec() - t0;

    // 4. AVX-512 (if CPU supports it)
    int have_avx512 = 0;
#if defined(__x86_64__) || defined(_M_X64)
    __builtin_cpu_init();
    have_avx512 = __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq");
#endif
    double t_avx512_single = 0.0;
    double t_avx512_omp = 0.0;
    if (have_avx512)
    {
        t0 = get_time_sec();
        run_extrusion_avx512_single((const float *const *) masters, msize, pup_size,
                                    nbframes, nblayers, vx, vy, cn2, out_pha);
        t_avx512_single = get_time_sec() - t0;

        t0 = get_time_sec();
        run_extrusion_avx512_omp((const float *const *) masters, msize, pup_size,
                                 nbframes, nblayers, vx, vy, cn2, out_pha);
        t_avx512_omp = get_time_sec() - t0;
    }

    // 5. CUDA GPU
    double t_cuda_cold = 0.0;
    double t_cuda_warm = 0.0;
    int have_gpu = atmturb_cuda_device_available();
    if (have_gpu)
    {
        atmturb_cuda_sim_params_t params = {
            .nblayers  = nblayers,
            .msize     = msize,
            .pup_size  = pup_size,
            .nbframes  = nbframes,
            .Scoeff    = 1.05,
            .h_masters = (const float *const *) masters,
            .vxpix     = vx,
            .vypix     = vy,
            .cn2       = cn2
        };
        atmturb_cuda_sim_outputs_t outputs = {
            .pha  = out_pha,
            .amp  = out_amp,
            .spha = out_spha,
            .samp = out_samp
        };

        // Cold upload (includes PCIe upload of 80MB master screens + alloc)
        atmturb_cuda_cleanup();
        t0 = get_time_sec();
        atmturb_wfs_render_frames_cuda(&params, &outputs);
        t_cuda_cold = get_time_sec() - t0;

        // Steady-state cached (re-uses resident screens in GPU VRAM)
        t0 = get_time_sec();
        atmturb_wfs_render_frames_cuda(&params, &outputs);
        t_cuda_warm = get_time_sec() - t0;
    }

    printf("%-30s | %10s | %12s | %8s\n", "Engine", "Time (ms)", "Throughput", "Speedup");
    printf("-------------------------------+------------+--------------+---------\n");
    printf("%-30s | %10.2f | %8.2f MP/s | %7.2fx\n", "Scalar (1 core)",
           t_scalar * 1000.0, (double) total_pixels / t_scalar / 1e6, 1.0);
    printf("%-30s | %10.2f | %8.2f MP/s | %7.2fx\n", "AVX2 SIMD (1 core)",
           t_avx2_single * 1000.0, (double) total_pixels / t_avx2_single / 1e6,
           t_scalar / t_avx2_single);
    printf("%-30s | %10.2f | %8.2f MP/s | %7.2fx\n", "AVX2 + OpenMP (24 threads)",
           t_avx2_omp * 1000.0, (double) total_pixels / t_avx2_omp / 1e6,
           t_scalar / t_avx2_omp);
    if (have_avx512)
    {
        printf("%-30s | %10.2f | %8.2f MP/s | %7.2fx\n", "AVX-512 SIMD (1 core)",
               t_avx512_single * 1000.0, (double) total_pixels / t_avx512_single / 1e6,
               t_scalar / t_avx512_single);
        printf("%-30s | %10.2f | %8.2f MP/s | %7.2fx\n", "AVX-512 + OpenMP",
               t_avx512_omp * 1000.0, (double) total_pixels / t_avx512_omp / 1e6,
               t_scalar / t_avx512_omp);
    }
    if (have_gpu)
    {
        printf("%-30s | %10.2f | %8.2f MP/s | %7.2fx\n", "CUDA GPU (Cold + 80MB Upload)",
               t_cuda_cold * 1000.0, (double) total_pixels / t_cuda_cold / 1e6,
               t_scalar / t_cuda_cold);
        printf("%-30s | %10.2f | %8.2f MP/s | %7.2fx\n", "CUDA GPU (VRAM Cached)",
               t_cuda_warm * 1000.0, (double) total_pixels / t_cuda_warm / 1e6,
               t_scalar / t_cuda_warm);
    }

    for (long k = 0; k < nblayers; k++)
    {
        free(masters[k]);
    }
    free(masters);
    free(vx);
    free(vy);
    free(cn2);
    free(out_pha);
    free(out_amp);
    free(out_spha);
    free(out_samp);
}

/**
 * run_fresnel_cpu_naive - Naive un-cached Fresnel propagation (plans FFTW per step)
 * @buf: Complex optical field.
 * @size: Grid dimension in pixels.
 * @pupil_scale: Physical pixel scale [m/pixel].
 * @z: Propagation distance [m].
 * @lambda: Light wavelength [m].
 */
static void run_fresnel_cpu_naive(
    complex_float *buf,
    long           size,
    double         pupil_scale,
    double         z,
    double         lambda)
{
    fftwf_plan forward = fftwf_plan_dft_2d((int) size, (int) size,
                                           (fftwf_complex *) buf, (fftwf_complex *) buf,
                                           FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_plan inverse = fftwf_plan_dft_2d((int) size, (int) size,
                                           (fftwf_complex *) buf, (fftwf_complex *) buf,
                                           FFTW_BACKWARD, FFTW_ESTIMATE);

    fftwf_execute(forward);

    double coeff = PI * z * lambda / (pupil_scale * size) / (pupil_scale * size);
    float co1 = (float) (size * size);
    long n0h = size / 2;

    #pragma omp parallel for
    for (long jj = 0; jj < size; jj++)
    {
        long jj1 = size * jj;
        long jj2 = (jj - n0h) * (jj - n0h);
        for (long ii = 0; ii < size; ii++)
        {
            long ii1 = jj1 + ii;
            long ii2 = ii - n0h;
            double sqdist = (double) (ii2 * ii2 + jj2);
            float angle = (float) (-coeff * sqdist);
            float s, c;
            sincosf(angle, &s, &c);
            float re = buf[ii1].re / co1;
            float im = buf[ii1].im / co1;
            buf[ii1].re = re * c - im * s;
            buf[ii1].im = re * s + im * c;
        }
    }

    fftwf_execute(inverse);

    fftwf_destroy_plan(forward);
    fftwf_destroy_plan(inverse);
}

/**
 * build_fresnel_tf - Precalculate 2D Fresnel optical transfer function
 * @tf: Destination transfer function array.
 * @size: Grid dimension in pixels.
 * @pupil_scale: Physical pixel scale [m/pixel].
 * @z: Propagation distance [m].
 * @lambda: Light wavelength [m].
 */
static void build_fresnel_tf(
    complex_float *tf,
    long           size,
    double         pupil_scale,
    double         z,
    double         lambda)
{
    double coeff = PI * z * lambda / (pupil_scale * size) / (pupil_scale * size);
    float co1 = (float) (size * size);
    long n0h = size / 2;

    #pragma omp parallel for
    for (long jj = 0; jj < size; jj++)
    {
        long jj1 = size * jj;
        long jj2 = (jj - n0h) * (jj - n0h);
        for (long ii = 0; ii < size; ii++)
        {
            long ii1 = jj1 + ii;
            long ii2 = ii - n0h;
            double sqdist = (double) (ii2 * ii2 + jj2);
            float angle = (float) (-coeff * sqdist);
            float s, c;
            sincosf(angle, &s, &c);
            tf[ii1].re = c / co1;
            tf[ii1].im = s / co1;
        }
    }
}

/**
 * run_fresnel_cpu_optimized - Optimized Fresnel with cached FFTW plans and SIMD transfer
 * @buf: Optical complex field.
 * @tf: Precomputed transfer function.
 * @forward: Cached forward 2D FFTW plan.
 * @inverse: Cached inverse 2D FFTW plan.
 * @size: Grid dimension in pixels.
 */
static void run_fresnel_cpu_optimized(
    complex_float       *buf,
    const complex_float *tf,
    fftwf_plan           forward,
    fftwf_plan           inverse,
    long                 size)
{
    fftwf_execute(forward);
    atmturb_complex_mul_array((float *) buf, (const float *) buf,
                              (const float *) tf, size * size);
    fftwf_execute(inverse);
}

/**
 * benchmark_fresnel - Benchmark diffractive propagation across CPU and GPU engines
 * @size: Linear grid dimension in pixels.
 * @iters: Number of test iterations.
 */
static void benchmark_fresnel(
    long size,
    int  iters)
{
    printf("\n========================================================================\n");
    printf(" Benchmark 2: Fresnel Diffractive Propagation (%ldx%ld grid, %d iters)\n",
           size, size, iters);
    printf("========================================================================\n");

    long nbelem = size * size;
    complex_float *h_in  = (complex_float *) malloc(sizeof(complex_float) * (size_t) nbelem);
    complex_float *h_out = (complex_float *) malloc(sizeof(complex_float) * (size_t) nbelem);
    complex_float *h_cpu = (complex_float *) malloc(sizeof(complex_float) * (size_t) nbelem);
    complex_float *tf    = (complex_float *) malloc(sizeof(complex_float) * (size_t) nbelem);

    for (long i = 0; i < nbelem; i++)
    {
        h_in[i].re = (float) rand() / (float) RAND_MAX;
        h_in[i].im = 0.0f;
    }

    double pupil_scale = 0.01;
    double z           = 1000.0;
    double lambda      = 0.5e-6;

    // 1. CPU Naive (plans created per frame, sincos per pixel)
    memcpy(h_cpu, h_in, sizeof(complex_float) * (size_t) nbelem);
    run_fresnel_cpu_naive(h_cpu, size, pupil_scale, z, lambda);
    double t0 = get_time_sec();
    for (int it = 0; it < iters; it++)
    {
        memcpy(h_cpu, h_in, sizeof(complex_float) * (size_t) nbelem);
        run_fresnel_cpu_naive(h_cpu, size, pupil_scale, z, lambda);
    }
    double t_cpu_naive = (get_time_sec() - t0) / (double) iters;

    // 2. CPU Optimized (cached FFTW plans + precomputed TF + AVX2 complex multiplication)
    build_fresnel_tf(tf, size, pupil_scale, z, lambda);
    fftwf_plan fwd = fftwf_plan_dft_2d((int) size, (int) size,
                                       (fftwf_complex *) h_cpu, (fftwf_complex *) h_cpu,
                                       FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_plan inv = fftwf_plan_dft_2d((int) size, (int) size,
                                       (fftwf_complex *) h_cpu, (fftwf_complex *) h_cpu,
                                       FFTW_BACKWARD, FFTW_ESTIMATE);

    memcpy(h_cpu, h_in, sizeof(complex_float) * (size_t) nbelem);
    run_fresnel_cpu_optimized(h_cpu, tf, fwd, inv, size);
    t0 = get_time_sec();
    for (int it = 0; it < iters; it++)
    {
        memcpy(h_cpu, h_in, sizeof(complex_float) * (size_t) nbelem);
        run_fresnel_cpu_optimized(h_cpu, tf, fwd, inv, size);
    }
    double t_cpu_opt = (get_time_sec() - t0) / (double) iters;
    fftwf_destroy_plan(fwd);
    fftwf_destroy_plan(inv);

    // 3. CUDA GPU (cuFFT)
    double t_cuda = 0.0;
    int have_gpu = wfprop_fresnel_device_available();
    if (have_gpu)
    {
        wfprop_fresnel_propagate_cuda(h_in, h_out, size, size, pupil_scale, z, lambda, 0);
        t0 = get_time_sec();
        for (int it = 0; it < iters; it++)
        {
            wfprop_fresnel_propagate_cuda(h_in, h_out, size, size, pupil_scale, z, lambda, 0);
        }
        t_cuda = (get_time_sec() - t0) / (double) iters;
    }

    printf("%-32s | %10s | %14s | %8s\n", "Engine", "Time (ms)", "Props / sec", "Speedup");
    printf("---------------------------------+------------+----------------+---------\n");
    printf("%-32s | %10.3f | %11.1f /s | %7.2fx\n", "CPU Naive (Plan Per Frame)",
           t_cpu_naive * 1000.0, 1.0 / t_cpu_naive, 1.0);
    printf("%-32s | %10.3f | %11.1f /s | %7.2fx\n", "CPU Opt (Plan Cache + AVX2 TF)",
           t_cpu_opt * 1000.0, 1.0 / t_cpu_opt, t_cpu_naive / t_cpu_opt);
    if (have_gpu)
    {
        printf("%-32s | %10.3f | %11.1f /s | %7.2fx\n", "CUDA GPU (cuFFT A5500)",
               t_cuda * 1000.0, 1.0 / t_cuda, t_cpu_naive / t_cuda);
    }

    free(h_in);
    free(h_out);
    free(h_cpu);
    free(tf);
}

/**
 * benchmark_simd_kernels - Microbenchmark raw vectorized mathematical kernels
 * @n: Number of elements in test arrays.
 * @iters: Number of benchmark repetitions.
 */
static void benchmark_simd_kernels(
    long n,
    int  iters)
{
    printf("\n========================================================================\n");
    printf(" Benchmark 3: SIMD Kernel Microbenchmarks (%ld elements, %d iters)\n", n, iters);
    printf("========================================================================\n");

    float *a    = (float *) malloc(sizeof(float) * (size_t) (2 * n));
    float *b    = (float *) malloc(sizeof(float) * (size_t) (2 * n));
    float *dest = (float *) malloc(sizeof(float) * (size_t) (2 * n));

    for (long i = 0; i < 2 * n; i++)
    {
        a[i] = (float) rand() / (float) RAND_MAX;
        b[i] = (float) rand() / (float) RAND_MAX;
    }

    int have_avx512 = 0;
#if defined(__x86_64__) || defined(_M_X64)
    __builtin_cpu_init();
    have_avx512 = __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq");
#endif

    if (have_avx512)
    {
        printf("%-30s | %10s | %10s | %10s | %8s\n", "Kernel Operation",
               "Scalar ms", "AVX2 ms", "AVX-512 ms", "Speedup");
        printf("-------------------------------+------------+------------+------------+---------\n");
    }
    else
    {
        printf("%-30s | %10s | %10s | %8s\n", "Kernel Operation",
               "Scalar ms", "AVX2 ms", "Speedup");
        printf("-------------------------------+------------+------------+---------\n");
    }

    // 1. atmturb_complex_mul_array
    double t0 = get_time_sec();
    for (int it = 0; it < iters; it++)
    {
        atmturb_complex_mul_array_scalar(dest, a, b, n);
    }
    double t_cmul_scalar = get_time_sec() - t0;

    t0 = get_time_sec();
    for (int it = 0; it < iters; it++)
    {
        atmturb_complex_mul_array_avx2(dest, a, b, n);
    }
    double t_cmul_avx2 = get_time_sec() - t0;

    double t_cmul_512 = 0.0;
    if (have_avx512)
    {
        t0 = get_time_sec();
        for (int it = 0; it < iters; it++)
        {
            atmturb_complex_mul_array_avx512(dest, a, b, n);
        }
        t_cmul_512 = get_time_sec() - t0;
        printf("%-30s | %10.2f | %10.2f | %10.2f | %7.2fx\n", "Complex Array Multiply",
               t_cmul_scalar * 1000.0, t_cmul_avx2 * 1000.0, t_cmul_512 * 1000.0,
               t_cmul_scalar / t_cmul_512);
    }
    else
    {
        printf("%-30s | %10.2f | %10.2f | %7.2fx\n", "Complex Array Multiply",
               t_cmul_scalar * 1000.0, t_cmul_avx2 * 1000.0, t_cmul_scalar / t_cmul_avx2);
    }

    // 2. atmturb_scale_float_array
    t0 = get_time_sec();
    for (int it = 0; it < iters; it++)
    {
        atmturb_scale_float_array_scalar(dest, a, 1.25f, n);
    }
    double t_scale_scalar = get_time_sec() - t0;

    t0 = get_time_sec();
    for (int it = 0; it < iters; it++)
    {
        atmturb_scale_float_array_avx2(dest, a, 1.25f, n);
    }
    double t_scale_avx2 = get_time_sec() - t0;

    double t_scale_512 = 0.0;
    if (have_avx512)
    {
        t0 = get_time_sec();
        for (int it = 0; it < iters; it++)
        {
            atmturb_scale_float_array_avx512(dest, a, 1.25f, n);
        }
        t_scale_512 = get_time_sec() - t0;
        printf("%-30s | %10.2f | %10.2f | %10.2f | %7.2fx\n", "Float Array Scale (4x unroll)",
               t_scale_scalar * 1000.0, t_scale_avx2 * 1000.0, t_scale_512 * 1000.0,
               t_scale_scalar / t_scale_512);
    }
    else
    {
        printf("%-30s | %10.2f | %10.2f | %7.2fx\n", "Float Array Scale (4x unroll)",
               t_scale_scalar * 1000.0, t_scale_avx2 * 1000.0, t_scale_scalar / t_scale_avx2);
    }

    // 3. atmturb_add_float_array
    t0 = get_time_sec();
    for (int it = 0; it < iters; it++)
    {
        atmturb_add_float_array_scalar(dest, a, n);
    }
    double t_add_scalar = get_time_sec() - t0;

    t0 = get_time_sec();
    for (int it = 0; it < iters; it++)
    {
        atmturb_add_float_array_avx2(dest, a, n);
    }
    double t_add_avx2 = get_time_sec() - t0;

    double t_add_512 = 0.0;
    if (have_avx512)
    {
        t0 = get_time_sec();
        for (int it = 0; it < iters; it++)
        {
            atmturb_add_float_array_avx512(dest, a, n);
        }
        t_add_512 = get_time_sec() - t0;
        printf("%-30s | %10.2f | %10.2f | %10.2f | %7.2fx\n", "Float Array Accumulate",
               t_add_scalar * 1000.0, t_add_avx2 * 1000.0, t_add_512 * 1000.0,
               t_add_scalar / t_add_512);
    }
    else
    {
        printf("%-30s | %10.2f | %10.2f | %7.2fx\n", "Float Array Accumulate",
               t_add_scalar * 1000.0, t_add_avx2 * 1000.0, t_add_scalar / t_add_avx2);
    }

    // 4. atmturb_init_phase_amp
    t0 = get_time_sec();
    for (int it = 0; it < iters; it++)
    {
        atmturb_init_phase_amp_scalar(a, b, n);
    }
    double t_init_scalar = get_time_sec() - t0;

    t0 = get_time_sec();
    for (int it = 0; it < iters; it++)
    {
        atmturb_init_phase_amp_avx2(a, b, n);
    }
    double t_init_avx2 = get_time_sec() - t0;

    double t_init_512 = 0.0;
    if (have_avx512)
    {
        t0 = get_time_sec();
        for (int it = 0; it < iters; it++)
        {
            atmturb_init_phase_amp_avx512(a, b, n);
        }
        t_init_512 = get_time_sec() - t0;
        printf("%-30s | %10.2f | %10.2f | %10.2f | %7.2fx\n", "Phase & Amp Dual Init (2x unroll)",
               t_init_scalar * 1000.0, t_init_avx2 * 1000.0, t_init_512 * 1000.0,
               t_init_scalar / t_init_512);
    }
    else
    {
        printf("%-30s | %10.2f | %10.2f | %7.2fx\n", "Phase & Amp Dual Init (2x unroll)",
               t_init_scalar * 1000.0, t_init_avx2 * 1000.0, t_init_scalar / t_init_avx2);
        printf("\n  [AVX-512 kernels compiled & supported in library; host CPU lacks AVX-512]\n");
    }

    free(a);
    free(b);
    free(dest);
}

/**
 * benchmark_lowfreq - Benchmark low-frequency subharmonic mode extrusion
 * @pup_size: Linear dimension of square pupil.
 * @iters: Number of benchmark repetitions.
 */
static void benchmark_lowfreq(
    long pup_size,
    int  iters)
{
    printf("\n========================================================================\n");
    printf(" Benchmark 4: Low-Frequency Subharmonics (%ldx%ld, 24 modes, %d iters)\n",
           pup_size, pup_size, iters);
    printf("========================================================================\n");

    atmturb_lowfreq_t lf;
    atmturb_lowfreq_init(&lf, 2048, 25.0, 100.0, 0.1, 42);

    long ptot = pup_size * pup_size;
    float *out_s = (float *) calloc((size_t) ptot, sizeof(float));
    float *out_v = (float *) calloc((size_t) ptot, sizeof(float));

    atmturb_lowfreq_params_t lp_s = {
        .lf         = &lf,
        .screen_idx = 0,
        .x0         = 12.34,
        .y0         = 56.78,
        .pup_size   = pup_size,
        .os         = 1,
        .weight     = 1.0f,
        .out_pha    = out_s
    };
    atmturb_lowfreq_params_t lp_v = lp_s;
    lp_v.out_pha = out_v;

    double t0 = get_time_sec();
    for (int it = 0; it < iters; it++)
    {
        atmturb_extrude_lowfreq_scalar(&lp_s);
    }
    double t_scalar = get_time_sec() - t0;

    t0 = get_time_sec();
    for (int it = 0; it < iters; it++)
    {
        atmturb_extrude_lowfreq_avx2(&lp_v);
    }
    double t_avx2 = get_time_sec() - t0;

    printf("%-30s | %10s | %10s | %8s\n", "Operation", "Scalar ms", "AVX2 ms", "Speedup");
    printf("-------------------------------+------------+------------+---------\n");
    printf("%-30s | %10.2f | %10.2f | %7.2fx\n", "Subharmonics (24 modes)",
           t_scalar * 1000.0, t_avx2 * 1000.0, t_scalar / t_avx2);

    free(out_s);
    free(out_v);
}

/**
 * main - Run full performance benchmark suite
 *
 * Return: 0 on success.
 */
int main(void)
{
    printf("########################################################################\n");
    printf("#  milkatmturb Performance Benchmark: SIMD (AVX2) vs OpenMP vs CUDA   #\n");
    printf("########################################################################\n");
    printf(" Active SIMD Vectorization ISA: %s\n", atmturb_simd_active_isa());

    // SIMD Microbenchmarks
    benchmark_simd_kernels(1024 * 1024, 500);

    // Low-Frequency Subharmonic Extrusion
    benchmark_lowfreq(256, 100);

    // Extrusion at different scales
    benchmark_extrusion(64, 100, 5);
    benchmark_extrusion(256, 100, 5);
    benchmark_extrusion(512, 100, 7);

    // Fresnel propagation at different scales
    benchmark_fresnel(256, 20);
    benchmark_fresnel(512, 20);
    benchmark_fresnel(1024, 10);
    benchmark_fresnel(2048, 5);

    printf("\nAll benchmarks completed successfully!\n");
    return 0;
}
