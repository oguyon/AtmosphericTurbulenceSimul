// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_mkhvturb_FPS.c
 * @brief   MILK_RUNNABLE FPS compute unit for Hufnagel-Valley turbulence profile generation
 *
 * Exposes `atmturb.mkhvturb` and `milk-fpsexec-atmturb-mkhvturb`. Discretizes a
 * continuous Hufnagel-Valley Cn2(h) and Bufton/legacy wind profile into a multi-layer
 * ASCII table ready for `atmturb.mkwfs`.
 */

#include <stdio.h>

#include "milk_runnable.h"
#include "AtmosphericTurbulence/AtmosphericTurbulence.h"
#include "atmturb_mkhvturb_FPS.h"

static double   param_wspeed                                = 20.0;
static double   param_r0                                    = 0.15;
static double   param_sitealt                               = 4200.0;
static int32_t  param_nblayer                               = 20;
static int32_t  param_wind_model                            = 1;
static uint64_t param_seed                                  = 1;
static char     param_outfile[FUNCTION_PARAMETER_STRMAXLEN] = "turbHV.prof";

#define ATMTURB_MKHVTURB_PARAMS(X)                                                 \
    X(".wspeed", &param_wspeed, FPTYPE_FLOAT64, 1, FPFLAG_DEFAULT_INPUT,           \
      "High altitude wind speed [m/s]")                                            \
    X(".r0", &param_r0, FPTYPE_FLOAT64, 1, FPFLAG_DEFAULT_INPUT,                   \
      "Fried parameter r0 [m] at 0.55um")                                          \
    X(".sitealt", &param_sitealt, FPTYPE_FLOAT64, 1, FPFLAG_DEFAULT_INPUT,         \
      "Site altitude [m]")                                                         \
    X(".nblayers", &param_nblayer, FPTYPE_INT32, 1, FPFLAG_DEFAULT_INPUT,          \
      "Number of layers")                                                          \
    X(".outfile", param_outfile, FPTYPE_FILENAME, 1,                               \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL,                                  \
      "Output profile filename (default: turbHV.prof)")                            \
    X(".wind_model", &param_wind_model, FPTYPE_INT32, 1,                           \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL,                                  \
      "Wind model (0 = legacy, 1 = Bufton)")                                       \
    X(".seed", &param_seed, FPTYPE_UINT64, 1,                                      \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL,                                  \
      "RNG seed (default: 1)")

MILK_PARAMS_DEFINE(atmturb_mkhvturb, ATMTURB_MKHVTURB_PARAMS);

/**
 * atmturb_mkhvturb_exec - Execute Hufnagel-Valley Cn2 profile synthesis
 *
 * Return: RETURN_SUCCESS on success, RETURN_FAILURE on error.
 */
static MILK_HOT errno_t atmturb_mkhvturb_exec(void)
{
    printf("[milkatmturb] Outputs:\n");
    printf("[milkatmturb]   - file:   \"%s\" (profile text)\n", param_outfile);
    fflush(stdout);

    if (AtmosphericTurbulence_makeHV_CN2prof_opt(param_wspeed, param_r0, param_sitealt,
                                                 (long) param_nblayer, param_outfile,
                                                 (int) param_wind_model, param_seed) != 0)
    {
        return RETURN_FAILURE;
    }

    return RETURN_SUCCESS;
}

/**
 * atmturb_mkhvturb_compute - Execute HV profile generation in procinfo loop
 * @self: Pointer to this MILK_RUNNABLE descriptor.
 *
 * Return: RETURN_SUCCESS on completion, RETURN_FAILURE on error.
 */
static MILK_HOT errno_t atmturb_mkhvturb_compute(
    const MILK_RUNNABLE *self)
{
    DEBUG_TRACE_FSTART();
    errno_t ret = RETURN_SUCCESS;

    MILK_PROCINFO_INIT(self)
    MILK_PROCINFO_LOOPSTART(self)
    ret = atmturb_mkhvturb_exec();
    MILK_PROCINFO_END(self)

    DEBUG_TRACE_FEXIT();
    return ret;
}

const MILK_RUNNABLE runnable_atmturb_mkhvturb = {
    .kind = MILK_RUNNABLE_FPS,
    .app  = {
        .fps_name         = "atmturb_mkhvturb",
        .cmdkey           = "mkhvturb",
        .description      = "Generate Hufnagel-Valley turbulence profile",
        .description_long = "Computes a discretized Hufnagel-Valley vertical turbulence "
                            "profile (Cn2, wind speed, direction, outer and inner scales) "
                            "based on high-altitude wind speed, Fried parameter r0, "
                            "site altitude, and layer count."
    },
    .srcfile = __FILE__,
    MILK_PARAMS_REF(atmturb_mkhvturb),
    .compute = atmturb_mkhvturb_compute,
};

EXPORT_FPSEXEC(milk, "atmturb-mkhvturb", runnable_atmturb_mkhvturb);
