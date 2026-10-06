// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_mkvonkarman_FPS.c
 * @brief   FPS V2 compute unit for von Karman wind velocity simulation
 */

#ifdef MILK_NO_CLI
#    include "CLIcore_standalone.h"
#else
#    include "CLIcore.h"
#endif
#include "fps.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "AtmosphericTurbulence/AtmosphericTurbulence.h"
#include "atmturb_mkvonkarman_FPS.h"

/* ================================================================
 * 1.  FPS COMPONENT IDENTITY
 * ============================================================= */

static FPS_APP_INFO FPS_app_info = {
    .fps_name         = "atmturb_mkvonkarman",
    .cmdkey           = "mkvonkarman",
    .description      = "Generate 1D von Karman wind velocity series [u, v, w]",
    .description_long = "Synthesizes a 3-channel 1D time series representing longitudinal, "
                        "transverse, and vertical turbulent wind fluctuations according to "
                        "the von Karman power spectrum."
};

/* ================================================================
 * 2.  LOCAL PARAMETER VARIABLES
 * ============================================================= */

static int32_t param_vksize    = 8192;
static float   param_pixscale  = 0.1f;
static float   param_sigmawind = 20.0f;
static float   param_lwind     = 50.0f;
static char    param_outname[FUNCTION_PARAMETER_STRMAXLEN] = "vKwind";

/* ================================================================
 * 3.  UNIFIED PARAMETER TABLE (X-Macro)
 * ============================================================= */

#define FPS_PARAMS(X)                                                                                  \
    X(".vksize", &param_vksize, FPTYPE_INT32, 1, FPFLAG_DEFAULT_INPUT, "Sample count of 1D series")     \
    X(".pixscale", &param_pixscale, FPTYPE_FLOAT32, 1, FPFLAG_DEFAULT_INPUT,                            \
      "Physical sampling step [m]")                                                                     \
    X(".sigmawind", &param_sigmawind, FPTYPE_FLOAT32, 1, FPFLAG_DEFAULT_INPUT,                          \
      "Velocity standard deviation [m/s]")                                                              \
    X(".lwind", &param_lwind, FPTYPE_FLOAT32, 1, FPFLAG_DEFAULT_INPUT,                                  \
      "Turbulence outer scale [m]")                                                                     \
    X(".outname", &param_outname, FPTYPE_STREAMNAME, 1, FPFLAG_DEFAULT_INPUT,                           \
      "Output 3D image name (vksize x 1 x 3)")

/* ================================================================
 * 4.  COMPUTATION LOGIC
 * ============================================================= */

/**
 * fpsexec - Execute von Karman wind series synthesis
 *
 * Return: RETURN_SUCCESS on success, error code otherwise.
 */
static MILK_HOT errno_t fpsexec(void)
{
    make_AtmosphericTurbulence_vonKarmanWind((long) param_vksize, param_pixscale,
                                            param_sigmawind, param_lwind, 0, param_outname);

    return RETURN_SUCCESS;
}

/* ================================================================
 * 5.  BINDINGS, FARG, AND CLI DATA
 * ============================================================= */

FPS_V2_SECTION5(FPS_PARAMS)

/* ================================================================
 * 6.  COMPUTE WRAPPER
 * ============================================================= */

/**
 * compute_function - Orchestrator wrapper for processinfo tracing
 *
 * Return: RETURN_SUCCESS on completion.
 */
static MILK_HOT errno_t __attribute__((unused)) compute_function(void)
{
    DEBUG_TRACE_FSTART();

    INSERT_STD_PROCINFO_COMPUTEFUNC_START
    fpsexec();
    INSERT_STD_PROCINFO_COMPUTEFUNC_END

    DEBUG_TRACE_FEXIT();
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

/**
 * CLIADDCMD_milkatmturb__atmturb_mkvonkarman_FPS - Register command with CLI
 *
 * Return: RETURN_SUCCESS on success.
 */
errno_t CLIADDCMD_milkatmturb__atmturb_mkvonkarman_FPS(void)
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
