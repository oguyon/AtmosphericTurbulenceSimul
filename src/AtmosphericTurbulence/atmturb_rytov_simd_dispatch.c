// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_rytov_simd_dispatch.c
 * @brief   Dynamic runtime CPU ISA dispatcher for Rytov SIMD kernels
 */

#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "atmturb_rytov_simd.h"

typedef struct
{
    void (*accumulate_filters)(
        fftwf_complex *restrict       acc_dphi,
        fftwf_complex *restrict       acc_chi,
        const fftwf_complex *restrict spec,
        const float *restrict         filt_a,
        const float *restrict         filt_b,
        long                          ntot);
    void (*accumulate_filters_rotated)(
        fftwf_complex *restrict       acc_dphi,
        fftwf_complex *restrict       acc_chi,
        const fftwf_complex *restrict spec,
        const fftwf_complex *restrict ramp,
        const float *restrict         filt_a,
        const float *restrict         filt_b,
        long                          ntot);
    void (*decompose_row)(
        fftwf_complex *restrict       spec_row,
        const fftwf_complex *restrict hat_a,
        const fftwf_complex *restrict exp_x,
        const float *restrict         laplace_inv_row,
        float                         ey_re,
        float                         ey_im,
        float                         b_re,
        float                         b_im,
        long                          n_half);
    const char *isa_name;
} atmturb_rytov_simd_ops_t;

static atmturb_rytov_simd_ops_t g_rytov_simd_ops;
static int g_rytov_simd_initialized = 0;

/**
 * atmturb_rytov_simd_init_dispatch - Initialize function pointers based on CPU capabilities
 */
static void atmturb_rytov_simd_init_dispatch(void)
{
    if (g_rytov_simd_initialized)
    {
        return;
    }

#if defined(__x86_64__) || defined(_M_X64)
    __builtin_cpu_init();

    const char *env_simd = getenv("ATMTURB_SIMD");
    if (env_simd != NULL)
    {
        if (strcasecmp(env_simd, "SCALAR") == 0)
        {
            goto use_scalar;
        }
        if (strcasecmp(env_simd, "AVX2") == 0 && __builtin_cpu_supports("avx2"))
        {
            goto use_avx2;
        }
    }

    if (__builtin_cpu_supports("avx2"))
    {
        goto use_avx2;
    }
    goto use_scalar;

use_avx2:
    g_rytov_simd_ops.accumulate_filters         = atmturb_rytov_accumulate_filters_avx2;
    g_rytov_simd_ops.accumulate_filters_rotated = atmturb_rytov_accumulate_filters_rotated_avx2;
    g_rytov_simd_ops.decompose_row              = atmturb_rytov_decompose_row_avx2;
    g_rytov_simd_ops.isa_name                   = "AVX2";
    g_rytov_simd_initialized                    = 1;
    return;

use_scalar:
#endif

    g_rytov_simd_ops.accumulate_filters         = atmturb_rytov_accumulate_filters_scalar;
    g_rytov_simd_ops.accumulate_filters_rotated = atmturb_rytov_accumulate_filters_rotated_scalar;
    g_rytov_simd_ops.decompose_row              = atmturb_rytov_decompose_row_scalar;
    g_rytov_simd_ops.isa_name                   = "Scalar";
    g_rytov_simd_initialized                    = 1;
}

__attribute__((constructor)) static void atmturb_rytov_simd_constructor(void)
{
    atmturb_rytov_simd_init_dispatch();
}

/**
 * atmturb_rytov_simd_accumulate_filters - Dispatch filter accumulation to optimal CPU kernel
 * @acc_dphi: Phase correction accumulator.
 * @acc_chi: Log-amplitude accumulator.
 * @spec: Periodic spectrum input.
 * @filt_a: Phase filter array.
 * @filt_b: Amplitude filter array.
 * @ntot: Elements in half-spectrum.
 */
void atmturb_rytov_simd_accumulate_filters(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          ntot)
{
    if (__builtin_expect(!g_rytov_simd_initialized, 0))
    {
        atmturb_rytov_simd_init_dispatch();
    }
    g_rytov_simd_ops.accumulate_filters(acc_dphi, acc_chi, spec, filt_a, filt_b, ntot);
}

/**
 * atmturb_rytov_simd_accumulate_filters_rotated - Dispatch rotated accumulation to kernel
 * @acc_dphi: Phase correction accumulator.
 * @acc_chi: Log-amplitude accumulator.
 * @spec: Periodic spectrum input.
 * @ramp: Complex chromatic phase factor array.
 * @filt_a: Phase filter array.
 * @filt_b: Amplitude filter array.
 * @ntot: Elements in half-spectrum.
 */
void atmturb_rytov_simd_accumulate_filters_rotated(
    fftwf_complex *restrict       acc_dphi,
    fftwf_complex *restrict       acc_chi,
    const fftwf_complex *restrict spec,
    const fftwf_complex *restrict ramp,
    const float *restrict         filt_a,
    const float *restrict         filt_b,
    long                          ntot)
{
    if (__builtin_expect(!g_rytov_simd_initialized, 0))
    {
        atmturb_rytov_simd_init_dispatch();
    }
    g_rytov_simd_ops.accumulate_filters_rotated(acc_dphi, acc_chi, spec, ramp,
                                                filt_a, filt_b, ntot);
}

/**
 * atmturb_rytov_simd_decompose_row - Dispatch Moisan decomposition row to optimal kernel
 * @spec_row: Target spectrum row complex array.
 * @hat_a: 1D Fourier transform of boundary A.
 * @exp_x: Precomputed complex exponential array along X.
 * @laplace_inv_row: Inverse discrete Laplacian row weights.
 * @ey_re: Real part of exp_y for current row.
 * @ey_im: Imaginary part of exp_y for current row.
 * @b_re: Real part of hat_b for current row.
 * @b_im: Imaginary part of hat_b for current row.
 * @n_half: Number of complex elements in half-row.
 */
void atmturb_rytov_simd_decompose_row(
    fftwf_complex *restrict       spec_row,
    const fftwf_complex *restrict hat_a,
    const fftwf_complex *restrict exp_x,
    const float *restrict         laplace_inv_row,
    float                         ey_re,
    float                         ey_im,
    float                         b_re,
    float                         b_im,
    long                          n_half)
{
    if (__builtin_expect(!g_rytov_simd_initialized, 0))
    {
        atmturb_rytov_simd_init_dispatch();
    }
    g_rytov_simd_ops.decompose_row(spec_row, hat_a, exp_x, laplace_inv_row,
                                   ey_re, ey_im, b_re, b_im, n_half);
}
