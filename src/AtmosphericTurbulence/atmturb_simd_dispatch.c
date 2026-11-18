// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_simd_dispatch.c
 * @brief   Dynamic runtime CPU ISA dispatcher for turbulence SIMD kernels
 */

#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "atmturb_simd.h"
#ifdef HAVE_CUDA
#include "atmturb_cuda.h"
#endif

typedef struct
{
    void (*extrude_bilinear)(const atmturb_extrude_params_t *params);
    void (*extrude_bicubic)(const atmturb_extrude_params_t *params);
    void (*scale_float_array)(float *dest, const float *src, float scale, long n);
    void (*init_phase_amp)(float *pha, float *amp, long n);
    void (*extrude_lowfreq)(const atmturb_lowfreq_params_t *params);
    void (*add_float_array)(float *dest, const float *src, long n);
    void (*complex_mul_array)(float *dest, const float *src1, const float *src2, long n);
    void (*remove_piston_stream)(float *dst, const float *src, long npix);
    const char *isa_name;
} atmturb_simd_ops_t;

static atmturb_simd_ops_t g_simd_ops;
static const char        *g_cpu_isa_name     = "Scalar";
static int                g_simd_gpu_forced  = 0;
static int                g_simd_initialized = 0;

/**
 * atmturb_simd_bind_scalar - Populate dispatch table with portable scalar kernels
 */
static void atmturb_simd_bind_scalar(void)
{
    g_simd_ops.extrude_bilinear     = atmturb_extrude_accumulate_bilinear_scalar;
    g_simd_ops.extrude_bicubic      = atmturb_extrude_accumulate_bicubic_scalar;
    g_simd_ops.scale_float_array    = atmturb_scale_float_array_scalar;
    g_simd_ops.init_phase_amp       = atmturb_init_phase_amp_scalar;
    g_simd_ops.extrude_lowfreq      = atmturb_extrude_lowfreq_scalar;
    g_simd_ops.add_float_array      = atmturb_add_float_array_scalar;
    g_simd_ops.complex_mul_array    = atmturb_complex_mul_array_scalar;
    g_simd_ops.remove_piston_stream = atmturb_remove_piston_stream_scalar;
    g_simd_ops.isa_name             = "Scalar";
    g_cpu_isa_name                  = "Scalar";
}

#if defined(__x86_64__) || defined(_M_X64)
/**
 * atmturb_simd_bind_avx2 - Populate dispatch table with AVX2 kernels
 */
static void atmturb_simd_bind_avx2(void)
{
    g_simd_ops.extrude_bilinear     = atmturb_extrude_accumulate_bilinear_avx2;
    g_simd_ops.extrude_bicubic      = atmturb_extrude_accumulate_bicubic_avx2;
    g_simd_ops.scale_float_array    = atmturb_scale_float_array_avx2;
    g_simd_ops.init_phase_amp       = atmturb_init_phase_amp_avx2;
    g_simd_ops.extrude_lowfreq      = atmturb_extrude_lowfreq_avx2;
    g_simd_ops.add_float_array      = atmturb_add_float_array_avx2;
    g_simd_ops.complex_mul_array    = atmturb_complex_mul_array_avx2;
    g_simd_ops.remove_piston_stream = atmturb_remove_piston_stream_avx2;
    g_simd_ops.isa_name             = "AVX2";
    g_cpu_isa_name                  = "AVX2";
}

/**
 * atmturb_simd_bind_avx512 - Populate dispatch table with AVX-512 kernels
 */
static void atmturb_simd_bind_avx512(void)
{
    g_simd_ops.extrude_bilinear     = atmturb_extrude_accumulate_bilinear_avx512;
    g_simd_ops.extrude_bicubic      = atmturb_extrude_accumulate_bicubic_avx512;
    g_simd_ops.scale_float_array    = atmturb_scale_float_array_avx512;
    g_simd_ops.init_phase_amp       = atmturb_init_phase_amp_avx512;
    g_simd_ops.extrude_lowfreq      = atmturb_extrude_lowfreq_avx512;
    g_simd_ops.add_float_array      = atmturb_add_float_array_avx512;
    g_simd_ops.complex_mul_array    = atmturb_complex_mul_array_avx512;
    g_simd_ops.remove_piston_stream = atmturb_remove_piston_stream_avx512;
    g_simd_ops.isa_name             = "AVX-512";
    g_cpu_isa_name                  = "AVX-512";
}

/**
 * atmturb_simd_bind_best_x86 - Bind highest supported host x86 SIMD kernel table
 */
static void atmturb_simd_bind_best_x86(void)
{
    if (__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq"))
    {
        atmturb_simd_bind_avx512();
    }
    else if (__builtin_cpu_supports("avx2"))
    {
        atmturb_simd_bind_avx2();
    }
    else
    {
        atmturb_simd_bind_scalar();
    }
}
#endif

/**
 * atmturb_simd_init_dispatch - Initialize function pointers based on CPU/GPU capabilities
 */
static void atmturb_simd_init_dispatch(void)
{
    if (g_simd_initialized)
    {
        return;
    }

#if defined(__x86_64__) || defined(_M_X64)
    __builtin_cpu_init();
    atmturb_simd_bind_best_x86();

    const char *env_simd = getenv("ATMTURB_SIMD");
    if (env_simd != NULL)
    {
#ifdef HAVE_CUDA
        if ((strcasecmp(env_simd, "CUDA") == 0 || strcasecmp(env_simd, "GPU") == 0) &&
            atmturb_cuda_device_available())
        {
            g_simd_ops.isa_name = "CUDA GPU";
            g_simd_gpu_forced   = 1;
            g_simd_initialized  = 1;
            return;
        }
#endif
        if (strcasecmp(env_simd, "SCALAR") == 0)
        {
            atmturb_simd_bind_scalar();
            g_simd_initialized = 1;
            return;
        }
        if ((strcasecmp(env_simd, "AVX512") == 0 || strcasecmp(env_simd, "AVX-512") == 0) &&
            __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq"))
        {
            atmturb_simd_bind_avx512();
            g_simd_initialized = 1;
            return;
        }
        if (strcasecmp(env_simd, "AVX2") == 0 && __builtin_cpu_supports("avx2"))
        {
            atmturb_simd_bind_avx2();
            g_simd_initialized = 1;
            return;
        }
    }

#ifdef HAVE_CUDA
    if (atmturb_cuda_device_available())
    {
        g_simd_ops.isa_name = "CUDA GPU";
    }
#endif
#else
    atmturb_simd_bind_scalar();
#endif

    g_simd_initialized = 1;
}

/**
 * atmturb_simd_constructor - Shared-library load-time initializer for extrusion SIMD dispatch table
 */
__attribute__((constructor)) static void atmturb_simd_constructor(void)
{
    atmturb_simd_init_dispatch();
}

/**
 * atmturb_simd_active_isa - Query name of active vectorized ISA implementation
 *
 * Return: String name of active ISA ("AVX-512", "AVX2", "CUDA GPU", or "Scalar").
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
 * atmturb_simd_cpu_isa - Query name of bound host CPU SIMD ISA implementation
 */
const char *atmturb_simd_cpu_isa(void)
{
    if (__builtin_expect(!g_simd_initialized, 0))
    {
        atmturb_simd_init_dispatch();
    }
    return g_cpu_isa_name;
}

/**
 * atmturb_simd_is_gpu - Query if CUDA GPU execution is active
 *
 * Return: 1 if GPU mode is selected, 0 otherwise.
 */
int atmturb_simd_is_gpu(void)
{
    if (__builtin_expect(!g_simd_initialized, 0))
    {
        atmturb_simd_init_dispatch();
    }
    return (g_simd_ops.isa_name != NULL &&
            strcmp(g_simd_ops.isa_name, "CUDA GPU") == 0) ? 1 : 0;
}

/**
 * atmturb_simd_is_gpu_forced - Query if CUDA GPU was explicitly requested via ATMTURB_SIMD
 */
int atmturb_simd_is_gpu_forced(void)
{
    if (__builtin_expect(!g_simd_initialized, 0))
    {
        atmturb_simd_init_dispatch();
    }
    return g_simd_gpu_forced;
}

/**
 * atmturb_extrude_accumulate - Dispatch phase screen extrusion to optimal CPU kernel
 * @params: Extrusion configuration and data pointers.
 */
void atmturb_extrude_accumulate(
    const atmturb_extrude_params_t *params)
{
    if (params == NULL || params->master == NULL || params->out_pha == NULL ||
        params->pup_size <= 0 || params->weight == 0.0f)
    {
        return;
    }
    if (__builtin_expect(!g_simd_initialized, 0))
    {
        atmturb_simd_init_dispatch();
    }
    if (params->interp == ATMTURB_INTERP_BICUBIC)
    {
        g_simd_ops.extrude_bicubic(params);
    }
    else
    {
        g_simd_ops.extrude_bilinear(params);
    }
}

/**
 * atmturb_scale_float_array - Dispatch array scaling to optimal CPU kernel
 * @dest: Output float array.
 * @src: Input float array.
 * @scale: Multiplicative scale factor.
 * @n: Number of elements.
 */
void atmturb_scale_float_array(
    float       *dest,
    const float *src,
    float        scale,
    long         n)
{
    if (dest == NULL || src == NULL || n <= 0)
    {
        return;
    }
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
void atmturb_init_phase_amp(
    float *pha,
    float *amp,
    long   n)
{
    if (n <= 0 || (pha == NULL && amp == NULL))
    {
        return;
    }
    if (__builtin_expect(!g_simd_initialized, 0))
    {
        atmturb_simd_init_dispatch();
    }
    g_simd_ops.init_phase_amp(pha, amp, n);
}

/**
 * atmturb_extrude_lowfreq - Dispatch low-frequency mode accumulation to optimal CPU kernel
 * @params: Low-frequency configuration and data pointers.
 */
void atmturb_extrude_lowfreq(
    const atmturb_lowfreq_params_t *params)
{
    if (params == NULL || params->lf == NULL || params->out_pha == NULL ||
        params->pup_size <= 0 || params->weight == 0.0f)
    {
        return;
    }
    if (__builtin_expect(!g_simd_initialized, 0))
    {
        atmturb_simd_init_dispatch();
    }
    g_simd_ops.extrude_lowfreq(params);
}

/**
 * atmturb_add_float_array - Dispatch array addition to optimal CPU kernel
 * @dest: Output/accumulator float array.
 * @src: Input float array.
 * @n: Number of elements.
 */
void atmturb_add_float_array(
    float       *dest,
    const float *src,
    long         n)
{
    if (dest == NULL || src == NULL || n <= 0)
    {
        return;
    }
    if (__builtin_expect(!g_simd_initialized, 0))
    {
        atmturb_simd_init_dispatch();
    }
    g_simd_ops.add_float_array(dest, src, n);
}

/**
 * atmturb_complex_mul_array - Dispatch complex multiplication to optimal CPU kernel
 * @dest: Output complex float array (length 2 * n_complex).
 * @src1: First input complex float array.
 * @src2: Second input complex float array.
 * @n_complex: Number of complex elements.
 */
void atmturb_complex_mul_array(
    float       *dest,
    const float *src1,
    const float *src2,
    long         n_complex)
{
    if (dest == NULL || src1 == NULL || src2 == NULL || n_complex <= 0)
    {
        return;
    }
    if (__builtin_expect(!g_simd_initialized, 0))
    {
        atmturb_simd_init_dispatch();
    }
    g_simd_ops.complex_mul_array(dest, src1, src2, n_complex);
}

/**
 * atmturb_remove_piston_stream - Dispatch piston removal to optimal CPU kernel
 * @dst: Destination phase array (can be equal to src for in-place).
 * @src: Source phase array.
 * @npix: Total number of pixels.
 */
void atmturb_remove_piston_stream(
    float       *dst,
    const float *src,
    long         npix)
{
    if (dst == NULL || src == NULL || npix <= 0)
    {
        return;
    }
    if (__builtin_expect(!g_simd_initialized, 0))
    {
        atmturb_simd_init_dispatch();
    }
    g_simd_ops.remove_piston_stream(dst, src, npix);
}
