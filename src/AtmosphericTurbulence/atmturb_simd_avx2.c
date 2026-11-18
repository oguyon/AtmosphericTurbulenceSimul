// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_simd_avx2.c
 * @brief   AVX2 target instantiation of compiler-vectorized 2D phase extrusion kernels
 *
 * Compiled with `-mavx2 -mfma -fno-math-errno` so the compiler lowers seam-split bilinear
 * and Keys bicubic inner loops in `atmturb_simd_extrude_impl.h` to contiguous 256-bit AVX2
 * loads, shuffles, and FMA instructions without slow `vgatherdps` instructions.
 */

#include "atmturb_simd.h"

#if defined(__AVX2__)

#include "atmturb_simd_extrude_impl.h"

/**
 * atmturb_extrude_accumulate_bilinear_avx2 - AVX2 accelerated bilinear extrusion
 * @params: Extrusion configuration and data pointers.
 */
void atmturb_extrude_accumulate_bilinear_avx2(
    const atmturb_extrude_params_t *params)
{
    atmturb_extrude_bilinear_impl(params);
}

/**
 * atmturb_extrude_accumulate_bicubic_avx2 - AVX2 accelerated Keys bicubic extrusion
 * @params: Extrusion configuration and data pointers.
 */
void atmturb_extrude_accumulate_bicubic_avx2(
    const atmturb_extrude_params_t *params)
{
    atmturb_extrude_bicubic_impl(params);
}

/**
 * atmturb_extrude_accumulate_avx2 - AVX2 extrusion dispatcher (bilinear/bicubic)
 * @params: Extrusion configuration and data pointers.
 */
void atmturb_extrude_accumulate_avx2(
    const atmturb_extrude_params_t *params)
{
    if (params->interp == ATMTURB_INTERP_BICUBIC)
    {
        atmturb_extrude_bicubic_impl(params);
    }
    else
    {
        atmturb_extrude_bilinear_impl(params);
    }
}

#else /* !defined(__AVX2__) */

void atmturb_extrude_accumulate_avx2(
    const atmturb_extrude_params_t *params)
{
    atmturb_extrude_accumulate_scalar(params);
}

void atmturb_extrude_accumulate_bilinear_avx2(
    const atmturb_extrude_params_t *params)
{
    atmturb_extrude_accumulate_bilinear_scalar(params);
}

void atmturb_extrude_accumulate_bicubic_avx2(
    const atmturb_extrude_params_t *params)
{
    atmturb_extrude_accumulate_bicubic_scalar(params);
}

#endif /* defined(__AVX2__) */
