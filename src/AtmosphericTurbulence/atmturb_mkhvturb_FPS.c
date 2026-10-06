// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_mkhvturb_FPS.c
 * @brief   FPS V2 compute unit for Hufnagel-Valley turbulence profile generation
 */

#ifdef MILK_NO_CLI
#    include "CLIcore_standalone.h"
#else
#    include "CLIcore.h"
#endif
#include "fps.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "AtmosphericTurbulence/AtmosphericTurbulence.h"

/* ================================================================
 * 1.  FPS COMPONENT IDENTITY
 * ============================================================= */

static FPS_APP_INFO FPS_app_info = {
    .fps_name         = "atmturb_mkhvturb",
    .cmdkey           = "mkhvturb_fps",
    .description      = "Generate Hufnagel-Valley turbulence profile",
    .description_long = "Computes a discretized Hufnagel-Valley vertical turbulence profile (Cn2, "
                        "wind speed, direction, outer and inner scales) based on high-altitude wind "
                        "speed, Fried parameter r0, site altitude, and layer count."
};

/* ================================================================
 * 2.  LOCAL PARAMETER VARIABLES
 * ============================================================= */

static double  param_wspeed  = 20.0;
static double  param_r0      = 0.15;
static double  param_sitealt = 4200.0;
static int32_t param_nblayer = 20;
static char    param_outfile[FUNCTION_PARAMETER_STRMAXLEN] = "turbHV.prof";

/* ================================================================
 * 3.  UNIFIED PARAMETER TABLE (X-Macro)
 * ============================================================= */

#define FPS_PARAMS(X)                                                                            \
    X(".wspeed", &param_wspeed, FPTYPE_FLOAT64, 1, FPFLAG_DEFAULT_INPUT,                        \
      "High altitude wind speed [m/s]")                                                          \
    X(".r0", &param_r0, FPTYPE_FLOAT64, 1, FPFLAG_DEFAULT_INPUT,                                 \
      "Fried parameter r0 [m] at 0.55um")                                                        \
    X(".sitealt", &param_sitealt, FPTYPE_FLOAT64, 1, FPFLAG_DEFAULT_INPUT, "Site altitude [m]")  \
    X(".nblayers", &param_nblayer, FPTYPE_INT32, 1, FPFLAG_DEFAULT_INPUT, "Number of layers")    \
    X(".outfile", &param_outfile, FPTYPE_FILENAME, 1, FPFLAG_DEFAULT_INPUT,                      \
      "Output profile filename (default: turbHV.prof)")

/* ================================================================
 * 4.  COMPUTATION LOGIC
 * ============================================================= */

static MILK_HOT errno_t fpsexec(void)
{
    AtmosphericTurbulence_makeHV_CN2prof(param_wspeed, param_r0, param_sitealt,
                                        (long) param_nblayer, param_outfile);

    return RETURN_SUCCESS;
}

/* ================================================================
 * 5.  BINDINGS, FARG, AND CLI DATA
 * ============================================================= */

FPS_V2_SECTION5(FPS_PARAMS)

/* ================================================================
 * 6.  COMPUTE WRAPPER
 * ============================================================= */

static MILK_HOT errno_t __attribute__((unused)) compute_function(void)
{
    INSERT_STD_PROCINFO_COMPUTEFUNC_START
    fpsexec();
    INSERT_STD_PROCINFO_COMPUTEFUNC_END

    return RETURN_SUCCESS;
}

/* ================================================================
 * 7.  MILK MODULE REGISTRATION
 * ============================================================= */

#if !defined(FPS_STANDALONE) && !defined(MILK_NO_CLI)
static errno_t CLIfunction(void)
{
    return safe_fps_generic_CLIfunction(&FPS_app_info, farg, &CLIcmddata, my_bindings,
                                        nb_bindings, compute_function);
}

errno_t CLIADDCMD_milkatmturb__atmturb_mkhvturb_FPS(void)
{
    safe_fps_fill_farg_examples(farg, my_bindings, nb_bindings);
    INSERT_STD_CLIREGISTERFUNC
    return RETURN_SUCCESS;
}
#endif

/* ================================================================
 * 8.  STANDALONE ENTRY POINT
 * ============================================================= */

#ifdef FPS_STANDALONE
FPS_MAIN_STANDALONE_V2(FPS_app_info, FPS_PARAMS, compute_function)
#endif
