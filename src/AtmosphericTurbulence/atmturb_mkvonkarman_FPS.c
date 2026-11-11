// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_mkvonkarman_FPS.c
 * @brief   MILK_RUNNABLE FPS compute unit for von Karman wind velocity simulation
 *
 * Exposes `atmturb.mkvonkarman` and `milk-fpsexec-atmturb-mkvonkarman`. Synthesizes
 * 3-component [u, v, w] 1D turbulent wind velocity fluctuations into a 3D shared-memory
 * image (`vksize x 1 x 3`) and optional FITS file.
 */

#include <stdio.h>

#include "milk_runnable.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "AtmosphericTurbulence/AtmosphericTurbulence.h"
#include "atmturb_mkvonkarman_FPS.h"

static int32_t param_vksize                                = 8192;
static float   param_pixscale                              = 0.1f;
static float   param_sigmawind                             = 20.0f;
static float   param_lwind                                 = 50.0f;
static char    param_outname[FUNCTION_PARAMETER_STRMAXLEN] = "vKwind";
static int64_t param_seed                                  = 1;
static char    param_fitsout[FUNCTION_PARAMETER_STRMAXLEN] = "";

#define ATMTURB_MKVONKARMAN_PARAMS(X)                                              \
    X(".vksize", &param_vksize, FPTYPE_INT32, 1, FPFLAG_DEFAULT_INPUT,             \
      "Sample count of 1D series")                                                 \
    X(".pixscale", &param_pixscale, FPTYPE_FLOAT32, 1, FPFLAG_DEFAULT_INPUT,       \
      "Physical sampling step [m]")                                                \
    X(".sigmawind", &param_sigmawind, FPTYPE_FLOAT32, 1, FPFLAG_DEFAULT_INPUT,     \
      "Velocity standard deviation [m/s]")                                         \
    X(".lwind", &param_lwind, FPTYPE_FLOAT32, 1, FPFLAG_DEFAULT_INPUT,             \
      "Turbulence outer scale [m]")                                                \
    X(".outname", param_outname, FPTYPE_STREAMNAME, 1,                             \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL,                                  \
      "Output 3D image name (vksize x 1 x 3)")                                     \
    X(".seed", &param_seed, FPTYPE_INT64, 1,                                       \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL, "RNG seed (0 = time-based)")     \
    X(".fitsout", param_fitsout, FPTYPE_STRING, 1,                                 \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL,                                  \
      "Optional FITS output file (empty = do not save)")

MILK_PARAMS_DEFINE(atmturb_mkvonkarman, ATMTURB_MKVONKARMAN_PARAMS);

/**
 * atmturb_mkvonkarman_exec - Execute von Karman wind series synthesis
 *
 * Return: RETURN_SUCCESS on success, RETURN_FAILURE otherwise.
 */
static MILK_HOT errno_t atmturb_mkvonkarman_exec(void)
{
    printf("[milkatmturb] Outputs:\n");
    printf("[milkatmturb]   - stream: \"%s\" (3D SHM, %ldx1x3 float32)\n",
           param_outname, (long) param_vksize);
    if (param_fitsout[0] != '\0')
    {
        printf("[milkatmturb]   - file:   \"%s\" (FITS)\n", param_fitsout);
    }
    fflush(stdout);

    long id = make_AtmosphericTurbulence_vonKarmanWind((long) param_vksize,
                                                       param_pixscale,
                                                       param_sigmawind,
                                                       param_lwind,
                                                       (long) param_seed,
                                                       param_outname);
    if (id < 0)
    {
        return RETURN_FAILURE;
    }

    if (param_fitsout[0] != '\0')
    {
        printf("[milkatmturb] Writing file: \"%s\" (FITS)...\n", param_fitsout);
        save_fits(param_outname, param_fitsout);
    }

    return RETURN_SUCCESS;
}

/**
 * atmturb_mkvonkarman_compute - Execute von Karman synthesis within procinfo loop
 * @self: Pointer to this MILK_RUNNABLE descriptor.
 *
 * Return: RETURN_SUCCESS on completion, or error code on failure.
 */
static MILK_HOT errno_t atmturb_mkvonkarman_compute(
    const MILK_RUNNABLE *self)
{
    DEBUG_TRACE_FSTART();
    errno_t ret = RETURN_SUCCESS;

    MILK_PROCINFO_INIT(self)
    MILK_PROCINFO_LOOPSTART(self)
    ret = atmturb_mkvonkarman_exec();
    MILK_PROCINFO_END(self)

    DEBUG_TRACE_FEXIT();
    return ret;
}

const MILK_RUNNABLE runnable_atmturb_mkvonkarman = {
    .kind = MILK_RUNNABLE_FPS,
    .app  = {
        .fps_name         = "atmturb_mkvonkarman",
        .cmdkey           = "mkvonkarman",
        .description      = "Generate 1D von Karman wind velocity series [u, v, w]",
        .description_long = "Synthesizes a 3-channel 1D time series representing longitudinal, "
                            "transverse, and vertical turbulent wind fluctuations according to "
                            "the von Karman power spectrum."
    },
    .srcfile = __FILE__,
    MILK_PARAMS_REF(atmturb_mkvonkarman),
    .compute = atmturb_mkvonkarman_compute,
};

EXPORT_FPSEXEC(milk, "atmturb-mkvonkarman", runnable_atmturb_mkvonkarman);
