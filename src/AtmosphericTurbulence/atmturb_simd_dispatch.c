// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_simd_dispatch.c
 * @brief   Dynamic runtime CPU ISA dispatcher for turbulence SIMD kernels
 */

#include "atmturb_simd.h"

typedef struct
{
    void (*extrude_accumulate)(const float *master, long msize, double x0, double y0,
                               long pup_size, float weight, float *out_pha);
    void (*scale_float_array)(float *dest, const float *src, float scale, long n);
    void (*init_phase_amp)(float *pha, float *amp, long n);
    const char *isa_name;
} atmturb_simd_ops_t;

static atmturb_simd_ops_t g_simd_ops;
static int g_simd_initialized = 0;

/**
 * atmturb_simd_init_dispatch - Initialize function pointers based on CPU capabilities
 */
static void atmturb_simd_init_dispatch(void)
{
    if (g_simd_initialized)
    {
        return;
    }

#if defined(__x86_64__) || defined(_M_X64)
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq"))
    {
        g_simd_ops.extrude_accumulate = atmturb_extrude_accumulate_avx512;
        g_simd_ops.scale_float_array  = atmturb_scale_float_array_avx512;
        g_simd_ops.init_phase_amp     = atmturb_init_phase_amp_avx512;
        g_simd_ops.isa_name           = "AVX-512";
        g_simd_initialized            = 1;
        return;
    }
    if (__builtin_cpu_supports("avx2"))
    {
        g_simd_ops.extrude_accumulate = atmturb_extrude_accumulate_avx2;
        g_simd_ops.scale_float_array  = atmturb_scale_float_array_avx2;
        g_simd_ops.init_phase_amp     = atmturb_init_phase_amp_avx2;
        g_simd_ops.isa_name           = "AVX2";
        g_simd_initialized            = 1;
        return;
    }
#endif

    g_simd_ops.extrude_accumulate = atmturb_extrude_accumulate_scalar;
    g_simd_ops.scale_float_array  = atmturb_scale_float_array_scalar;
    g_simd_ops.init_phase_amp     = atmturb_init_phase_amp_scalar;
    g_simd_ops.isa_name           = "Scalar";
    g_simd_initialized            = 1;
}

__attribute__((constructor)) static void atmturb_simd_constructor(void)
{
    atmturb_simd_init_dispatch();
}

/**
 * atmturb_simd_active_isa - Query name of active vectorized ISA implementation
 *
 * Return: String name of active ISA ("AVX-512", "AVX2", or "Scalar").
 */
const char *atmturb_simd_active_isa(void)
{
    if (__builtin_expect(!g_simd_initialized, 0))
    {
        atmturb_simd_init_dispatch();
    }
    return g_simd_ops.isa_name;
}

/**
 * atmturb_extrude_accumulate - Dispatch bilinear extrusion to optimal CPU kernel
 * @master: Pointer to master phase screen data.
 * @msize: Master screen dimension.
 * @x0: Sub-pixel X coordinate offset.
 * @y0: Sub-pixel Y coordinate offset.
 * @pup_size: Linear dimension of the extracted pupil.
 * @weight: Layer Cn2 amplitude weight factor.
 * @out_pha: Output pupil phase array to accumulate into.
 */
void atmturb_extrude_accumulate(const float *master, long msize, double x0, double y0,
                                long pup_size, float weight, float *out_pha)
{
    if (__builtin_expect(!g_simd_initialized, 0))
    {
        atmturb_simd_init_dispatch();
    }
    g_simd_ops.extrude_accumulate(master, msize, x0, y0, pup_size, weight, out_pha);
}

/**
 * atmturb_scale_float_array - Dispatch array scaling to optimal CPU kernel
 * @dest: Output float array.
 * @src: Input float array.
 * @scale: Multiplicative scale factor.
 * @n: Number of elements.
 */
void atmturb_scale_float_array(float *dest, const float *src, float scale, long n)
{
    if (__builtin_expect(!g_simd_initialized, 0))
    {
        atmturb_simd_init_dispatch();
    }
    g_simd_ops.scale_float_array(dest, src, scale, n);
}

/**
 * atmturb_init_phase_amp - Dispatch array initialization to optimal CPU kernel
 * @pha: Output phase array (zeroed).
 * @amp: Output amplitude array (set to 1.0).
 * @n: Number of elements.
 */
void atmturb_init_phase_amp(float *pha, float *amp, long n)
{
    if (__builtin_expect(!g_simd_initialized, 0))
    {
        atmturb_simd_init_dispatch();
    }
    g_simd_ops.init_phase_amp(pha, amp, n);
}
