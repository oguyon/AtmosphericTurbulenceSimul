// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    benchmark_streaming.c
 * @brief   Benchmark testing for 2D wavefront streaming mode (CPU vs GPU)
 */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "atmturb_geometry.h"
#include "atmturb_profile.h"
#include "atmturb_rolling.h"
#include "atmturb_rytov.h"
#include "atmturb_simd.h"
#include "atmturb_cuda.h"
#include "atmturb_wfs_render.h"
#include "atmturb_wfs_render_cuda.h"

#define BENCH_MSIZE 1024
#define BENCH_NLAYERS 7

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
 * init_synthetic_screen - Generate synthetic phase screen
 * @screen: Destination buffer for phase screen.
 * @size: Linear screen dimension in pixels.
 * @seed: Random seed for synthesis.
 */
static void init_synthetic_screen(
    float *screen,
    long   size,
    int    seed)
{
    srand((unsigned int) seed);
    for (long i = 0; i < size * size; i++)
    {
        screen[i] = ((float) rand() / (float) RAND_MAX - 0.5f) * 10.0f;
    }
}

/**
 * run_cpu_geom_step_scalar - Render one geometric frame using scalar extrusion
 * @rsim: Rolling context.
 * @geom: Observing geometry.
 * @t: Frame index.
 * @pup_size: Pupil size.
 * @pha: Destination phase buffer.
 */
static void run_cpu_geom_step_scalar(
    const atmturb_rolling_t *rsim,
    const atmturb_geom_t    *geom,
    long                     t,
    long                     pup_size,
    float                   *pha)
{
    long npix = pup_size * pup_size;
    for (long i = 0; i < npix; i++)
    {
        pha[i] = 0.0f;
    }
    for (int k = 0; k < geom->nlayers; k++)
    {
        atmturb_extrude_params_t ep = {
            .master   = rsim->layers[k].screens[0].data,
            .msize    = rsim->msize,
            .x0       = geom->layers[k].x0 + (double) t * geom->layers[k].vx_pix,
            .y0       = geom->layers[k].y0 + (double) t * geom->layers[k].vy_pix,
            .pup_size = pup_size,
            .os       = 1,
            .interp   = ATMTURB_INTERP_BILINEAR,
            .weight   = (float) geom->layers[k].weight,
            .out_pha  = pha
        };
        atmturb_extrude_accumulate_scalar(&ep);
    }
}

/**
 * run_cpu_geom_step_avx2 - Render one geometric frame using AVX2 extrusion
 * @rsim: Rolling context.
 * @geom: Observing geometry.
 * @t: Frame index.
 * @pup_size: Pupil size.
 * @pha: Destination phase buffer.
 */
static void run_cpu_geom_step_avx2(
    const atmturb_rolling_t *rsim,
    const atmturb_geom_t    *geom,
    long                     t,
    long                     pup_size,
    float                   *pha)
{
    long npix = pup_size * pup_size;
    for (long i = 0; i < npix; i++)
    {
        pha[i] = 0.0f;
    }
    for (int k = 0; k < geom->nlayers; k++)
    {
        atmturb_extrude_params_t ep = {
            .master   = rsim->layers[k].screens[0].data,
            .msize    = rsim->msize,
            .x0       = geom->layers[k].x0 + (double) t * geom->layers[k].vx_pix,
            .y0       = geom->layers[k].y0 + (double) t * geom->layers[k].vy_pix,
            .pup_size = pup_size,
            .os       = 1,
            .interp   = ATMTURB_INTERP_BILINEAR,
            .weight   = (float) geom->layers[k].weight,
            .out_pha  = pha
        };
        atmturb_extrude_accumulate_avx2(&ep);
    }
}

/**
 * run_cpu_omp_geom_step - Render one geometric frame with multi-threaded row/layer reduction
 * @rsim: Rolling context.
 * @geom: Observing geometry.
 * @t: Frame index.
 * @pup_size: Pupil size.
 * @pha: Destination phase buffer.
 * @nthreads: Number of worker threads.
 * @thread_bufs: Pre-allocated thread scratch buffers.
 */
static void run_cpu_omp_geom_step(
    const atmturb_rolling_t *rsim,
    const atmturb_geom_t    *geom,
    long                     t,
    long                     pup_size,
    float                   *pha,
    int                      nthreads,
    float                  **thread_bufs)
{
    long npix = pup_size * pup_size;
    #pragma omp parallel num_threads(nthreads)
    {
        int tid = omp_get_thread_num();
        float *buf = (tid == 0) ? pha : thread_bufs[tid];
        for (long i = 0; i < npix; i++)
        {
            buf[i] = 0.0f;
        }

        #pragma omp for schedule(static)
        for (int k = 0; k < geom->nlayers; k++)
        {
            atmturb_extrude_params_t ep = {
                .master   = rsim->layers[k].screens[0].data,
                .msize    = rsim->msize,
                .x0       = geom->layers[k].x0 + (double) t * geom->layers[k].vx_pix,
                .y0       = geom->layers[k].y0 + (double) t * geom->layers[k].vy_pix,
                .pup_size = pup_size,
                .os       = 1,
                .interp   = ATMTURB_INTERP_BILINEAR,
                .weight   = (float) geom->layers[k].weight,
                .out_pha  = buf
            };
            atmturb_extrude_accumulate_avx2(&ep);
        }
    }

    for (int tid = 1; tid < nthreads; tid++)
    {
        atmturb_add_float_array(pha, thread_bufs[tid], npix);
    }
}

/**
 * bench_piston_removal - Benchmark piston removal and publish strategies
 * @pup_size: Linear dimension of pupil.
 * @nframes: Number of frames to test.
 */
static void bench_piston_removal(
    long pup_size,
    int  nframes)
{
    long npix = pup_size * pup_size;
    size_t nbytes = (size_t) npix * sizeof(float);
    float *in_pha = (float *) malloc(nbytes);
    float *out_shm = (float *) malloc(nbytes);
    if (!in_pha || !out_shm)
    {
        free(in_pha);
        free(out_shm);
        return;
    }
    init_synthetic_screen(in_pha, pup_size, 42);

    atmturb_remove_piston_stream_scalar(in_pha, out_shm, npix);
    atmturb_remove_piston_stream_avx2(in_pha, out_shm, npix);

    double t0 = get_time_sec();
    for (int i = 0; i < nframes; i++)
    {
        double sum = 0.0;
        for (long j = 0; j < npix; j++)
        {
            sum += (double) in_pha[j];
        }
        float mean = (float) (sum / (double) npix);
        for (long j = 0; j < npix; j++)
        {
            in_pha[j] -= mean;
        }
        memcpy(out_shm, in_pha, nbytes);
    }
    double t1 = get_time_sec();
    double dt_scalar = (t1 - t0) / (double) nframes;

    t0 = get_time_sec();
    for (int i = 0; i < nframes; i++)
    {
        atmturb_remove_piston_stream_avx2(in_pha, out_shm, npix);
    }
    t1 = get_time_sec();
    double dt_avx2 = (t1 - t0) / (double) nframes;

    printf("  Piston & Publish (%4ldx%4ld): Scalar: %6.2f us | AVX2: %6.2f us (%.1fx faster)\n",
           pup_size, pup_size, dt_scalar * 1e6, dt_avx2 * 1e6, dt_scalar / dt_avx2);

    free(in_pha);
    free(out_shm);
}

/**
 * bench_geom_streaming - Benchmark geometric 2D streaming across CPU and GPU
 * @rsim: Rolling context.
 * @geom: Observing geometry.
 * @pup_size: Linear dimension of pupil.
 * @nframes: Number of frames to test.
 */
static void bench_geom_streaming(
    const atmturb_rolling_t *rsim,
    const atmturb_geom_t    *geom,
    long                     pup_size,
    int                      nframes)
{
    long npix = pup_size * pup_size;
    float *pha = (float *) malloc((size_t) npix * sizeof(float));
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
    if (nthreads > 8)
    {
        nthreads = 8;
    }
#endif
    float **thread_bufs = (float **) malloc((size_t) nthreads * sizeof(float *));
    for (int i = 0; i < nthreads; i++)
    {
        thread_bufs[i] = (float *) malloc((size_t) npix * sizeof(float));
    }

    run_cpu_geom_step_avx2(rsim, geom, 0, pup_size, pha);

    double t0 = get_time_sec();
    for (int t = 0; t < nframes; t++)
    {
        run_cpu_geom_step_scalar(rsim, geom, t, pup_size, pha);
    }
    double dt_scalar = (get_time_sec() - t0) / (double) nframes;

    t0 = get_time_sec();
    for (int t = 0; t < nframes; t++)
    {
        run_cpu_geom_step_avx2(rsim, geom, t, pup_size, pha);
    }
    double dt_avx2 = (get_time_sec() - t0) / (double) nframes;

    t0 = get_time_sec();
    for (int t = 0; t < nframes; t++)
    {
        run_cpu_omp_geom_step(rsim, geom, t, pup_size, pha, nthreads, thread_bufs);
    }
    double dt_omp = (get_time_sec() - t0) / (double) nframes;

    double dt_gpu = 0.0;
    atmturb_cuda_geom_stream_t *gstream =
        atmturb_wfs_cuda_geom_stream_init(rsim, geom, BENCH_MSIZE, pup_size);
    if (gstream != NULL)
    {
        atmturb_wfs_cuda_geom_stream_render_step(gstream, 0, pha, NULL, NULL, NULL);
        t0 = get_time_sec();
        for (int t = 0; t < nframes; t++)
        {
            atmturb_wfs_cuda_geom_stream_render_step(gstream, t, pha, NULL, NULL, NULL);
        }
        dt_gpu = (get_time_sec() - t0) / (double) nframes;
        atmturb_wfs_cuda_geom_stream_free(gstream);
    }

    double mp = (double) npix * 1e-6;
    printf("  Geom Stream (%4ldx%4ld, %d layers):\n", pup_size, pup_size, geom->nlayers);
    printf("    - CPU Scalar (1T):    %8.2f us/frame | %7.1f fps | %6.1f MP/s | 1.0x\n",
           dt_scalar * 1e6, 1.0 / dt_scalar, mp / dt_scalar);
    printf("    - CPU AVX2   (1T):    %8.2f us/frame | %7.1f fps | %6.1f MP/s | %4.1fx\n",
           dt_avx2 * 1e6, 1.0 / dt_avx2, mp / dt_avx2, dt_scalar / dt_avx2);
    printf("    - CPU AVX2   (%dT):   %8.2f us/frame | %7.1f fps | %6.1f MP/s | %4.1fx\n",
           nthreads, dt_omp * 1e6, 1.0 / dt_omp, mp / dt_omp, dt_scalar / dt_omp);
    if (dt_gpu > 0.0)
    {
        printf("    - CUDA GPU Streaming: %8.2f us/frame | %7.1f fps | %6.1f MP/s | %4.1fx\n",
               dt_gpu * 1e6, 1.0 / dt_gpu, mp / dt_gpu, dt_scalar / dt_gpu);
    }

    for (int i = 0; i < nthreads; i++)
    {
        free(thread_bufs[i]);
    }
    free(thread_bufs);
    free(pha);
}

/**
 * bench_rytov_streaming - Benchmark Rytov diffractive streaming (CPU vs GPU)
 * @rsim: Rolling context.
 * @geom: Observing geometry.
 * @prof: Atmospheric profile.
 * @params: Observing parameters.
 * @pup_size: Linear dimension of pupil.
 * @nframes: Number of frames to test.
 */
static void bench_rytov_streaming(
    const atmturb_rolling_t    *rsim,
    const atmturb_geom_t       *geom,
    const atmturb_profile_t    *prof,
    const atmturb_obs_params_t *params,
    long                        pup_size,
    int                         nframes)
{
    long npix = pup_size * pup_size;
    float *pha  = (float *) malloc((size_t) npix * sizeof(float));
    float *amp  = (float *) malloc((size_t) npix * sizeof(float));
    float *spha = (float *) malloc((size_t) npix * sizeof(float));
    float *samp = (float *) malloc((size_t) npix * sizeof(float));

    double dt_cpu = 0.0;
    atmturb_rytov_plan_t rplan;
    atmturb_rytov_ctx_t rctx;
    double z_bin = 2000.0;
    if (atmturb_rytov_plan_init(&rplan, prof, geom, pup_size, params->pupil_scale_m,
                                params->lambda_ref_m, params->lambda_s_m, z_bin, 0) == 0)
    {
        atmturb_rytov_ctx_init(&rctx, rplan.pad_size);
        atmturb_rytov_render_step(&rctx, &rplan, rsim, geom, 0, 0.001,
                                  BENCH_MSIZE, pup_size, pha, amp, spha, samp);

        double t0 = get_time_sec();
        for (int t = 0; t < nframes; t++)
        {
            atmturb_rytov_render_step(&rctx, &rplan, rsim, geom, t, 0.001,
                                      BENCH_MSIZE, pup_size, pha, amp, spha, samp);
        }
        dt_cpu = (get_time_sec() - t0) / (double) nframes;
        atmturb_rytov_ctx_free(&rctx);
        atmturb_rytov_plan_free(&rplan);
    }

    double dt_gpu = 0.0;
    atmturb_cuda_rytov_stream_t *crstream =
        atmturb_cuda_rytov_stream_init(rsim, geom, prof, params, BENCH_MSIZE, pup_size);
    if (crstream != NULL)
    {
        atmturb_cuda_rytov_stream_render_step(crstream, rsim, geom, 0, 0.001,
                                             pha, amp, spha, samp);
        double t0 = get_time_sec();
        for (int t = 0; t < nframes; t++)
        {
            atmturb_cuda_rytov_stream_render_step(crstream, rsim, geom, t, 0.001,
                                                 pha, amp, spha, samp);
        }
        dt_gpu = (get_time_sec() - t0) / (double) nframes;
        atmturb_cuda_rytov_stream_free(crstream);
    }

    double mp = (double) npix * 1e-6;
    printf("  Rytov Diffractive Stream (%4ldx%4ld):\n", pup_size, pup_size);
    if (dt_cpu > 0.0)
    {
        printf("    - CPU AVX2 Rytov:     %8.2f us/frame | %7.1f fps | %6.1f MP/s | 1.0x\n",
               dt_cpu * 1e6, 1.0 / dt_cpu, mp / dt_cpu);
    }
    if (dt_gpu > 0.0 && dt_cpu > 0.0)
    {
        printf("    - CUDA GPU Rytov:     %8.2f us/frame | %7.1f fps | %6.1f MP/s | %4.1fx\n",
               dt_gpu * 1e6, 1.0 / dt_gpu, mp / dt_gpu, dt_cpu / dt_gpu);
    }

    free(pha);
    free(amp);
    free(spha);
    free(samp);
}

/**
 * main - Streaming benchmark harness entry point
 * @argc: Command argument count.
 * @argv: Command argument vector.
 *
 * Return: 0 on success, non-zero on error.
 */
int main(
    int   argc,
    char *argv[])
{
    (void) argc;
    (void) argv;

    printf("========================================================================\n");
    printf(" milkatmturb 2D Streaming Performance Benchmark (CPU vs GPU)\n");
    printf("========================================================================\n");
    printf("ISA active: %s | CUDA available: %s\n",
           atmturb_simd_active_isa(),
           atmturb_cuda_device_available() ? "YES" : "NO");

    atmturb_rolling_t rsim;
    memset(&rsim, 0, sizeof(rsim));
    rsim.nlayers = BENCH_NLAYERS;
    rsim.msize   = BENCH_MSIZE;
    rsim.layers  = (atmturb_rolling_layer_t *) calloc(
        (size_t) BENCH_NLAYERS, sizeof(atmturb_rolling_layer_t));

    for (int k = 0; k < BENCH_NLAYERS; k++)
    {
        rsim.layers[k].nscreens = 1;
        rsim.layers[k].screens = (atmturb_rolling_screen_t *) calloc(
            1, sizeof(atmturb_rolling_screen_t));
        rsim.layers[k].screens[0].data =
            (float *) malloc(sizeof(float) * BENCH_MSIZE * BENCH_MSIZE);
        init_synthetic_screen(rsim.layers[k].screens[0].data, BENCH_MSIZE, 100 + k);
    }

    atmturb_geom_t geom;
    memset(&geom, 0, sizeof(geom));
    geom.nlayers = BENCH_NLAYERS;
    geom.oversample = 1.0;
    geom.interp = ATMTURB_INTERP_BILINEAR;
    geom.layers = (atmturb_layer_geom_t *) calloc(
        (size_t) BENCH_NLAYERS, sizeof(atmturb_layer_geom_t));
    for (int k = 0; k < BENCH_NLAYERS; k++)
    {
        geom.layers[k].weight   = 0.35;
        geom.layers[k].weight_s = 0.35;
        geom.layers[k].vx_pix   = 1.5 + (double) k * 0.5;
        geom.layers[k].vy_pix   = -1.0 + (double) k * 0.3;
        geom.layers[k].x0       = 50.0 + (double) k * 10.0;
        geom.layers[k].y0       = 50.0 + (double) k * 10.0;
    }

    atmturb_profile_t prof;
    memset(&prof, 0, sizeof(prof));
    atmturb_profile_load(NULL, &prof);

    atmturb_obs_params_t params;
    memset(&params, 0, sizeof(params));
    params.pupil_scale_m = 0.05;
    params.lambda_ref_m  = 0.65e-6;
    params.lambda_s_m    = 0.85e-6;
    params.time_step_s   = 0.001;

    long sizes[] = {64, 128, 256, 512};
    int n_sizes = sizeof(sizes) / sizeof(sizes[0]);

    printf("\n--- Section 1: Piston Removal & Publish to Shared Memory ---\n");
    for (int i = 0; i < n_sizes; i++)
    {
        bench_piston_removal(sizes[i], 1000);
    }

    printf("\n--- Section 2: Geometric Wavefront Streaming (7 layers) ---\n");
    for (int i = 0; i < n_sizes; i++)
    {
        int nframes = (sizes[i] <= 128) ? 300 : (sizes[i] == 256 ? 150 : 50);
        bench_geom_streaming(&rsim, &geom, sizes[i], nframes);
    }

    printf("\n--- Section 3: Diffractive Rytov Wavefront Streaming ---\n");
    for (int i = 0; i < n_sizes; i++)
    {
        int nframes = (sizes[i] <= 128) ? 100 : (sizes[i] == 256 ? 50 : 20);
        bench_rytov_streaming(&rsim, &geom, &prof, &params, sizes[i], nframes);
    }

    printf("\n========================================================================\n");
    printf(" Benchmark Completed Successfully\n");
    printf("========================================================================\n");

    atmturb_profile_free(&prof);
    free(geom.layers);
    for (int k = 0; k < BENCH_NLAYERS; k++)
    {
        free(rsim.layers[k].screens[0].data);
        free(rsim.layers[k].screens);
    }
    free(rsim.layers);
    return 0;
}
