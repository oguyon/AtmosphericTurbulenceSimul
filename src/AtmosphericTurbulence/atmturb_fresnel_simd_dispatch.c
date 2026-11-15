// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_fresnel_simd_dispatch.c
 * @brief   Dynamic runtime CPU ISA dispatcher for Fresnel SIMD kernels
 */

#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "atmturb_fresnel_simd.h"

typedef struct
{
    void (*modulate_field)(
        fftwf_complex *restrict field,
        const float   *restrict phase,
        long                    npix);
    void (*init_modulate_field)(
        fftwf_complex *restrict field,
        const float   *restrict phase,
        long                    npix);
    void (*extract_unwrapped)(
        const fftwf_complex *restrict field,
        float               *restrict pha,
        float               *restrict amp,
        long                          npix);
    void (*phasor_masked)(
        fftwf_complex *restrict field,
        const float   *restrict opd,
        const float   *restrict mask,
        float                   scale,
        long                    npix);
    void (*spectrum_strip)(
        fftwf_complex *restrict out,
        const float   *restrict u1,
        const float   *restrict u2,
        float                   fx0,
        const float   *restrict spec_p);
    const char *isa_name;
} atmturb_fresnel_simd_ops_t;

static atmturb_fresnel_simd_ops_t g_fresnel_simd_ops;
static int g_fresnel_simd_initialized = 0;

/**
 * atmturb_fresnel_simd_init_dispatch - Initialize function pointers based on CPU capabilities
 */
static void atmturb_fresnel_simd_init_dispatch(void)
{
    if (g_fresnel_simd_initialized)
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
        if ((strcasecmp(env_simd, "AVX512") == 0 || strcasecmp(env_simd, "AVX-512") == 0) &&
            __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq") &&
            __builtin_cpu_supports("avx512vl"))
        {
            goto use_avx512;
        }
        if (strcasecmp(env_simd, "AVX2") == 0 && __builtin_cpu_supports("avx2"))
        {
            goto use_avx2;
        }
    }

    if (__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq") &&
        __builtin_cpu_supports("avx512vl"))
    {
        goto use_avx512;
    }
    if (__builtin_cpu_supports("avx2"))
    {
        goto use_avx2;
    }
    goto use_scalar;

use_avx512:
    g_fresnel_simd_ops.modulate_field      = atmturb_fresnel_modulate_field_avx512;
    g_fresnel_simd_ops.init_modulate_field = atmturb_fresnel_init_modulate_field_avx512;
    g_fresnel_simd_ops.extract_unwrapped   = atmturb_fresnel_extract_unwrapped_avx512;
    g_fresnel_simd_ops.phasor_masked       = atmturb_fresnel_phasor_masked_avx512;
    g_fresnel_simd_ops.spectrum_strip      = atmturb_fresnel_spectrum_strip_avx512;
    g_fresnel_simd_ops.isa_name            = "AVX-512";
    g_fresnel_simd_initialized             = 1;
    return;

use_avx2:
    g_fresnel_simd_ops.modulate_field      = atmturb_fresnel_modulate_field_avx2;
    g_fresnel_simd_ops.init_modulate_field = atmturb_fresnel_init_modulate_field_avx2;
    g_fresnel_simd_ops.extract_unwrapped   = atmturb_fresnel_extract_unwrapped_avx2;
    g_fresnel_simd_ops.phasor_masked       = atmturb_fresnel_phasor_masked_avx2;
    g_fresnel_simd_ops.spectrum_strip      = atmturb_fresnel_spectrum_strip_avx2;
    g_fresnel_simd_ops.isa_name            = "AVX2";
    g_fresnel_simd_initialized             = 1;
    return;

use_scalar:
#endif

    g_fresnel_simd_ops.modulate_field      = atmturb_fresnel_modulate_field_scalar;
    g_fresnel_simd_ops.init_modulate_field = atmturb_fresnel_init_modulate_field_scalar;
    g_fresnel_simd_ops.extract_unwrapped   = atmturb_fresnel_extract_unwrapped_scalar;
    g_fresnel_simd_ops.phasor_masked       = atmturb_fresnel_phasor_masked_scalar;
    g_fresnel_simd_ops.spectrum_strip      = atmturb_fresnel_spectrum_strip_scalar;
    g_fresnel_simd_ops.isa_name            = "Scalar";
    g_fresnel_simd_initialized             = 1;
}

/**
 * atmturb_fresnel_simd_constructor - Shared-library load-time initializer for Fresnel SIMD table
 */
__attribute__((constructor)) static void atmturb_fresnel_simd_constructor(void)
{
    atmturb_fresnel_simd_init_dispatch();
}

/**
 * atmturb_fresnel_simd_modulate_field - Dispatch complex field modulation to optimal kernel
 * @field: In-out complex field array.
 * @phase: Real phase screen array in radians.
 * @npix: Total number of pixels.
 */
void atmturb_fresnel_simd_modulate_field(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    npix)
{
    if (__builtin_expect(!g_fresnel_simd_initialized, 0))
    {
        atmturb_fresnel_simd_init_dispatch();
    }
    g_fresnel_simd_ops.modulate_field(field, phase, npix);
}

/**
 * atmturb_fresnel_simd_init_modulate_field - Dispatch unit plane wave init to optimal kernel
 * @field: Output complex field array.
 * @phase: Real phase screen array in radians.
 * @npix: Total number of pixels.
 */
void atmturb_fresnel_simd_init_modulate_field(
    fftwf_complex *restrict field,
    const float   *restrict phase,
    long                    npix)
{
    if (__builtin_expect(!g_fresnel_simd_initialized, 0))
    {
        atmturb_fresnel_simd_init_dispatch();
    }
    g_fresnel_simd_ops.init_modulate_field(field, phase, npix);
}

/**
 * atmturb_fresnel_simd_extract_unwrapped - Dispatch unwrapped extraction to optimal kernel
 * @field: Complex optical field at pupil.
 * @pha: Phase array holding accumulated geometric phase (updated in-place).
 * @amp: Destination amplitude array (can be NULL).
 * @npix: Total number of pixels.
 */
void atmturb_fresnel_simd_extract_unwrapped(
    const fftwf_complex *restrict field,
    float               *restrict pha,
    float               *restrict amp,
    long                          npix)
{
    if (__builtin_expect(!g_fresnel_simd_initialized, 0))
    {
        atmturb_fresnel_simd_init_dispatch();
    }
    g_fresnel_simd_ops.extract_unwrapped(field, pha, amp, npix);
}

/**
 * atmturb_fresnel_simd_phasor_masked - Dispatch masked pupil complex field synthesis
 */
void atmturb_fresnel_simd_phasor_masked(
    fftwf_complex *restrict field,
    const float   *restrict opd,
    const float   *restrict mask,
    float                   scale,
    long                    npix)
{
    if (__builtin_expect(!g_fresnel_simd_initialized, 0))
    {
        atmturb_fresnel_simd_init_dispatch();
    }
    g_fresnel_simd_ops.phasor_masked(field, opd, mask, scale, npix);
}

/**
 * atmturb_fresnel_simd_spectrum_strip - Dispatch 8-lane von Karman Fourier mode synthesis
 */
void atmturb_fresnel_simd_spectrum_strip(
    fftwf_complex *restrict out,
    const float   *restrict u1,
    const float   *restrict u2,
    float                   fx0,
    const float   *restrict spec_p)
{
    if (__builtin_expect(!g_fresnel_simd_initialized, 0))
    {
        atmturb_fresnel_simd_init_dispatch();
    }
    g_fresnel_simd_ops.spectrum_strip(out, u1, u2, fx0, spec_p);
}
