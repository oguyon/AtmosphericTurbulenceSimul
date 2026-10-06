// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    benchmark_perf.c
 * @brief   Performance benchmark comparing Scalar, AVX2, AVX-512, OpenMP, and CUDA
 */

#define _GNU_SOURCE
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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

static double get_time_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static void run_extrusion_scalar(const float *const *masters, long msize, long pup_size,
                                 long nbframes, long nblayers, const double *vx,
                                 const double *vy, const double *cn2, float *out_pha)
{
    long frame_pixels = pup_size * pup_size;
    for (long t = 0; t < nbframes; t++)
    {
        float *pha_slice = &out_pha[t * frame_pixels];
        for (long i = 0; i < frame_pixels; i++) pha_slice[i] = 0.0f;

        for (long k = 0; k < nblayers; k++)
        {
            double cur_x = (double)(t + 1) * vx[k];
            double cur_y = (double)(t + 1) * vy[k];
            float weight = (float)sqrt(cn2[k]);
            atmturb_extrude_accumulate_scalar(masters[k], msize, cur_x, cur_y,
                                              pup_size, weight, pha_slice);
        }
    }
}

static void run_extrusion_avx2_single(const float *const *masters, long msize, long pup_size,
                                      long nbframes, long nblayers, const double *vx,
                                      const double *vy, const double *cn2, float *out_pha)
{
    long frame_pixels = pup_size * pup_size;
    for (long t = 0; t < nbframes; t++)
    {
        float *pha_slice = &out_pha[t * frame_pixels];
        for (long i = 0; i < frame_pixels; i++) pha_slice[i] = 0.0f;

        for (long k = 0; k < nblayers; k++)
        {
            double cur_x = (double)(t + 1) * vx[k];
            double cur_y = (double)(t + 1) * vy[k];
            float weight = (float)sqrt(cn2[k]);
            atmturb_extrude_accumulate_avx2(masters[k], msize, cur_x, cur_y,
                                            pup_size, weight, pha_slice);
        }
    }
}

static void run_extrusion_avx2_omp(const float *const *masters, long msize, long pup_size,
                                   long nbframes, long nblayers, const double *vx,
                                   const double *vy, const double *cn2, float *out_pha)
{
    long frame_pixels = pup_size * pup_size;
    #pragma omp parallel for schedule(dynamic)
    for (long t = 0; t < nbframes; t++)
    {
        float *pha_slice = &out_pha[t * frame_pixels];
        for (long i = 0; i < frame_pixels; i++) pha_slice[i] = 0.0f;

        for (long k = 0; k < nblayers; k++)
        {
            double cur_x = (double)(t + 1) * vx[k];
            double cur_y = (double)(t + 1) * vy[k];
            float weight = (float)sqrt(cn2[k]);
            atmturb_extrude_accumulate_avx2(masters[k], msize, cur_x, cur_y,
                                            pup_size, weight, pha_slice);
        }
    }
}

static void benchmark_extrusion(long pup_size, long nbframes, long nblayers)
{
    long msize = 2048;
    printf("\n========================================================================\n");
    printf(" Benchmark 1: Wavefront Extrusion (%ldx%ld pupil, %ld frames, %ld layers)\n",
           pup_size, pup_size, nbframes, nblayers);
    printf("========================================================================\n");

    float **masters = (float **)malloc(sizeof(float *) * nblayers);
    for (long k = 0; k < nblayers; k++)
    {
        masters[k] = (float *)malloc(sizeof(float) * msize * msize);
        for (long i = 0; i < msize * msize; i++)
        {
            masters[k][i] = (float)rand() / (float)RAND_MAX;
        }
    }

    double *vx = (double *)malloc(sizeof(double) * nblayers);
    double *vy = (double *)malloc(sizeof(double) * nblayers);
    double *cn2 = (double *)malloc(sizeof(double) * nblayers);
    for (long k = 0; k < nblayers; k++)
    {
        vx[k] = 0.5 + 0.2 * k;
        vy[k] = 0.3 - 0.1 * k;
        cn2[k] = 1.0 / (k + 1);
    }

    long total_pixels = nbframes * pup_size * pup_size;
    float *out_pha = (float *)malloc(sizeof(float) * total_pixels);
    float *out_amp = (float *)malloc(sizeof(float) * total_pixels);
    float *out_spha = (float *)malloc(sizeof(float) * total_pixels);
    float *out_samp = (float *)malloc(sizeof(float) * total_pixels);

    // Warm up
    run_extrusion_scalar((const float *const *)masters, msize, pup_size, 2, nblayers, vx, vy, cn2, out_pha);

    // 1. Scalar (1 thread)
    double t0 = get_time_sec();
    run_extrusion_scalar((const float *const *)masters, msize, pup_size, nbframes, nblayers, vx, vy, cn2, out_pha);
    double t_scalar = get_time_sec() - t0;

    // 2. AVX2 (1 thread)
    t0 = get_time_sec();
    run_extrusion_avx2_single((const float *const *)masters, msize, pup_size, nbframes, nblayers, vx, vy, cn2, out_pha);
    double t_avx2_single = get_time_sec() - t0;

    // 3. AVX2 + OpenMP (multi-threaded)
    t0 = get_time_sec();
    run_extrusion_avx2_omp((const float *const *)masters, msize, pup_size, nbframes, nblayers, vx, vy, cn2, out_pha);
    double t_avx2_omp = get_time_sec() - t0;

    // 4. CUDA GPU
    double t_cuda = 0.0;
    int have_gpu = atmturb_cuda_device_available();
    if (have_gpu)
    {
        atmturb_cuda_sim_params_t params = {
            .nblayers = nblayers,
            .msize = msize,
            .pup_size = pup_size,
            .nbframes = nbframes,
            .Scoeff = 1.05,
            .h_masters = (const float *const *)masters,
            .vxpix = vx,
            .vypix = vy,
            .cn2 = cn2
        };
        atmturb_cuda_sim_outputs_t outputs = {
            .pha = out_pha,
            .amp = out_amp,
            .spha = out_spha,
            .samp = out_samp
        };

        // Warm up CUDA
        atmturb_wfs_render_frames_cuda(&params, &outputs);

        t0 = get_time_sec();
        atmturb_wfs_render_frames_cuda(&params, &outputs);
        t_cuda = get_time_sec() - t0;
    }

    printf("%-26s | %10s | %12s | %8s\n", "Engine", "Time (ms)", "Throughput", "Speedup");
    printf("---------------------------+------------+--------------+---------\n");
    printf("%-26s | %10.2f | %8.2f MP/s | %7.2fx\n", "Scalar (1 core)",
           t_scalar * 1000.0, (double)total_pixels / t_scalar / 1e6, 1.0);
    printf("%-26s | %10.2f | %8.2f MP/s | %7.2fx\n", "AVX2 SIMD (1 core)",
           t_avx2_single * 1000.0, (double)total_pixels / t_avx2_single / 1e6, t_scalar / t_avx2_single);
    printf("%-26s | %10.2f | %8.2f MP/s | %7.2fx\n", "AVX2 + OpenMP (24 threads)",
           t_avx2_omp * 1000.0, (double)total_pixels / t_avx2_omp / 1e6, t_scalar / t_avx2_omp);
    if (have_gpu)
    {
        printf("%-26s | %10.2f | %8.2f MP/s | %7.2fx\n", "NVIDIA RTX A5500 CUDA GPU",
               t_cuda * 1000.0, (double)total_pixels / t_cuda / 1e6, t_scalar / t_cuda);
    }

    for (long k = 0; k < nblayers; k++) free(masters[k]);
    free(masters);
    free(vx);
    free(vy);
    free(cn2);
    free(out_pha);
    free(out_amp);
    free(out_spha);
    free(out_samp);
}

#include <fftw3.h>

static void run_fresnel_cpu(complex_float *buf, long size, double pupil_scale, double z, double lambda)
{
    fftwf_plan forward = fftwf_plan_dft_2d((int)size, (int)size,
                                           (fftwf_complex *)buf, (fftwf_complex *)buf,
                                           FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_plan inverse = fftwf_plan_dft_2d((int)size, (int)size,
                                           (fftwf_complex *)buf, (fftwf_complex *)buf,
                                           FFTW_BACKWARD, FFTW_ESTIMATE);

    fftwf_execute(forward);

    double coeff = PI * z * lambda / (pupil_scale * size) / (pupil_scale * size);
    float co1 = (float)(size * size);
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
            double sqdist = (double)(ii2 * ii2 + jj2);
            float angle = (float)(-coeff * sqdist);
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

static void benchmark_fresnel(long size, int iters)
{
    printf("\n========================================================================\n");
    printf(" Benchmark 2: Fresnel Diffractive Propagation (%ldx%ld grid, %d iters)\n",
           size, size, iters);
    printf("========================================================================\n");

    long nbelem = size * size;
    complex_float *h_in  = (complex_float *)malloc(sizeof(complex_float) * nbelem);
    complex_float *h_out = (complex_float *)malloc(sizeof(complex_float) * nbelem);
    complex_float *h_cpu = (complex_float *)malloc(sizeof(complex_float) * nbelem);

    for (long i = 0; i < nbelem; i++)
    {
        h_in[i].re = (float)rand() / (float)RAND_MAX;
        h_in[i].im = 0.0f;
    }

    double pupil_scale = 0.01;
    double z = 1000.0;
    double lambda = 0.5e-6;

    // 1. CPU (FFTW + OpenMP)
    memcpy(h_cpu, h_in, sizeof(complex_float) * nbelem);
    run_fresnel_cpu(h_cpu, size, pupil_scale, z, lambda); // warm up
    double t0 = get_time_sec();
    for (int it = 0; it < iters; it++)
    {
        memcpy(h_cpu, h_in, sizeof(complex_float) * nbelem);
        run_fresnel_cpu(h_cpu, size, pupil_scale, z, lambda);
    }
    double t_cpu = (get_time_sec() - t0) / iters;

    // 2. CUDA GPU (cuFFT)
    double t_cuda = 0.0;
    int have_gpu = wfprop_fresnel_device_available();
    if (have_gpu)
    {
        wfprop_fresnel_propagate_cuda(h_in, h_out, size, size, pupil_scale, z, lambda, 0); // warm up

        t0 = get_time_sec();
        for (int it = 0; it < iters; it++)
        {
            wfprop_fresnel_propagate_cuda(h_in, h_out, size, size, pupil_scale, z, lambda, 0);
        }
        t_cuda = (get_time_sec() - t0) / iters;
    }

    printf("%-26s | %10s | %14s | %8s\n", "Engine", "Time (ms)", "Props / sec", "Speedup");
    printf("---------------------------+------------+----------------+---------\n");
    printf("%-26s | %10.3f | %11.1f /s | %7.2fx\n", "CPU (FFTW + OpenMP 24t)",
           t_cpu * 1000.0, 1.0 / t_cpu, 1.0);
    if (have_gpu)
    {
        printf("%-26s | %10.3f | %11.1f /s | %7.2fx\n", "CUDA GPU (cuFFT A5500)",
               t_cuda * 1000.0, 1.0 / t_cuda, t_cpu / t_cuda);
    }

    free(h_in);
    free(h_out);
    free(h_cpu);
}

int main(void)
{
    printf("########################################################################\n");
    printf("#  milkatmturb Performance Benchmark: SIMD (AVX2) vs OpenMP vs CUDA   #\n");
    printf("########################################################################\n");

    // Extrusion at different scales
    benchmark_extrusion(64, 100, 5);    // Small pupil
    benchmark_extrusion(256, 100, 5);   // Medium pupil
    benchmark_extrusion(512, 100, 7);   // ELT scale pupil

    // Fresnel propagation at different scales
    benchmark_fresnel(256, 20);
    benchmark_fresnel(512, 20);
    benchmark_fresnel(1024, 10);
    benchmark_fresnel(2048, 5);

    printf("\nBenchmark completed successfully!\n");
    return 0;
}
