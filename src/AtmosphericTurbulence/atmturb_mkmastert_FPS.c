// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_mkmastert_FPS.c
 * @brief   MILK_RUNNABLE FPS compute unit for master phase screens generation
 *
 * Exposes `atmturb.mkmastert` and `milk-fpsexec-atmturb-mkmastert`. Synthesizes an
 * independent pair of 2D Kolmogorov / von Karman master phase screens via complex
 * Hermitian FFTW filtering normalized to r0 = 6.88^(3/5) pixels.
 */

#include <stdio.h>

#include "milk_runnable.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "AtmosphericTurbulence/AtmosphericTurbulence.h"
#include "atmturb_mkmastert_FPS.h"

static int32_t param_size                                    = 2048;
static float   param_outerscale                              = 100.0f;
static float   param_innerscale                              = 1.0f;
static int32_t param_precision                               = 0;
static char    param_screen0[FUNCTION_PARAMETER_STRMAXLEN]   = "turbm00_p0";
static char    param_screen1[FUNCTION_PARAMETER_STRMAXLEN]   = "turbm00_p1";
static int64_t param_seed                                    = 1;
static char    param_fitsout0[FUNCTION_PARAMETER_STRMAXLEN]  = "";
static char    param_fitsout1[FUNCTION_PARAMETER_STRMAXLEN]  = "";

#define ATMTURB_MKMASTERT_PARAMS(X)                                                \
    X(".size", &param_size, FPTYPE_INT32, 1, FPFLAG_DEFAULT_INPUT,                 \
      "Screen grid dimension [pixels]")                                            \
    X(".outerscale", &param_outerscale, FPTYPE_FLOAT32, 1, FPFLAG_DEFAULT_INPUT,   \
      "Outer scale in grid units [pixels]")                                        \
    X(".innerscale", &param_innerscale, FPTYPE_FLOAT32, 1, FPFLAG_DEFAULT_INPUT,   \
      "Inner scale in grid units [pixels]")                                        \
    X(".precision", &param_precision, FPTYPE_INT32, 1,                             \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL,                                  \
      "Precision (0=single, 1=double)")                                            \
    X(".screen0", param_screen0, FPTYPE_STREAMNAME, 1,                             \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL, "Output screen 0 image name")    \
    X(".screen1", param_screen1, FPTYPE_STREAMNAME, 1,                             \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL, "Output screen 1 image name")    \
    X(".seed", &param_seed, FPTYPE_INT64, 1,                                       \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL, "RNG seed (0 = time-based)")     \
    X(".fitsout0", param_fitsout0, FPTYPE_STRING, 1,                               \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL,                                  \
      "Optional FITS output file for screen 0 (empty = do not save)")              \
    X(".fitsout1", param_fitsout1, FPTYPE_STRING, 1,                               \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL,                                  \
      "Optional FITS output file for screen 1 (empty = do not save)")

MILK_PARAMS_DEFINE(atmturb_mkmastert, ATMTURB_MKMASTERT_PARAMS);

/**
 * atmturb_mkmastert_exec - Execute master phase screen generation
 *
 * Return: RETURN_SUCCESS on success, RETURN_FAILURE otherwise.
 */
static MILK_HOT errno_t atmturb_mkmastert_exec(void)
{
    printf("[milkatmturb] Outputs:\n");
    printf("[milkatmturb]   - stream: \"%s\" (2D SHM, %ldx%ld float32)\n",
           param_screen0, (long) param_size, (long) param_size);
    printf("[milkatmturb]   - stream: \"%s\" (2D SHM, %ldx%ld float32)\n",
           param_screen1, (long) param_size, (long) param_size);
    if (param_fitsout0[0] != '\0')
    {
        printf("[milkatmturb]   - file:   \"%s\" (FITS)\n", param_fitsout0);
    }
    if (param_fitsout1[0] != '\0')
    {
        printf("[milkatmturb]   - file:   \"%s\" (FITS)\n", param_fitsout1);
    }
    fflush(stdout);

    int ret = make_master_turbulence_screen_seeded(param_screen0, param_screen1,
                                                   (long) param_size,
                                                   param_outerscale, param_innerscale,
                                                   (long) param_precision,
                                                   (uint64_t) param_seed);
    if (ret != 0)
    {
        return RETURN_FAILURE;
    }

    if (param_fitsout0[0] != '\0')
    {
        printf("[milkatmturb] Writing file: \"%s\" (FITS)...\n", param_fitsout0);
        save_fits(param_screen0, param_fitsout0);
    }
    if (param_fitsout1[0] != '\0')
    {
        printf("[milkatmturb] Writing file: \"%s\" (FITS)...\n", param_fitsout1);
        save_fits(param_screen1, param_fitsout1);
    }

    return RETURN_SUCCESS;
}

/**
 * atmturb_mkmastert_compute - Execute master screen synthesis within procinfo loop
 * @self: Pointer to this MILK_RUNNABLE descriptor.
 *
 * Return: RETURN_SUCCESS on completion, or error code on failure.
 */
static MILK_HOT errno_t atmturb_mkmastert_compute(
    const MILK_RUNNABLE *self)
{
    DEBUG_TRACE_FSTART();
    errno_t ret = RETURN_SUCCESS;

    MILK_PROCINFO_INIT(self)
    MILK_PROCINFO_LOOPSTART(self)
    ret = atmturb_mkmastert_exec();
    MILK_PROCINFO_END(self)

    DEBUG_TRACE_FEXIT();
    return ret;
}

const MILK_RUNNABLE runnable_atmturb_mkmastert = {
    .kind = MILK_RUNNABLE_FPS,
    .app  = {
        .fps_name         = "atmturb_mkmastert",
        .cmdkey           = "mkmastert",
        .description      = "Generate master Kolmogorov/von Karman turbulence screens",
        .description_long = "Generates a pair of normalized master phase screens in Fourier "
                            "domain with outer and inner scale cutoff filtering."
    },
    .srcfile = __FILE__,
    MILK_PARAMS_REF(atmturb_mkmastert),
    .compute = atmturb_mkmastert_compute,
};

EXPORT_FPSEXEC(milk, "atmturb-mkmastert", runnable_atmturb_mkmastert);
