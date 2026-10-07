// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    wfprop_fresnel_FPS.c
 * @brief   FPS V2 compute unit for Fresnel wavefront propagation
 */

#ifdef MILK_NO_CLI
#    include "CLIcore_standalone.h"
#else
#    include "CLIcore.h"
#endif
#include "fps.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "WFpropagate/WFpropagate.h"
#include "wfprop_fresnel_FPS.h"

/* ================================================================
 * 1.  FPS COMPONENT IDENTITY
 * ============================================================= */

static FPS_APP_INFO FPS_app_info = {
    .fps_name         = "wfprop_fresnel",
    .cmdkey           = "fresnelprop",
    .description      = "Fresnel diffractive propagation of optical wavefront",
    .description_long = "Propagates a 2D complex optical field over a distance z using "
                        "the Fourier-domain Fresnel quadratic phase transfer function."
};

/* ================================================================
 * 2.  LOCAL PARAMETER VARIABLES
 * ============================================================= */

static double param_pupilscale = 0.01;
static double param_z          = 1000.0;
static double param_lambda     = 0.5e-6;
static char   param_inname[FUNCTION_PARAMETER_STRMAXLEN]  = "wfin";
static char   param_outname[FUNCTION_PARAMETER_STRMAXLEN] = "wfout";

/* ================================================================
 * 3.  UNIFIED PARAMETER TABLE (X-Macro)
 * ============================================================= */

#define FPS_PARAMS(X)                                                           \
    X(".inname", &param_inname, FPTYPE_STREAMNAME, 1, FPFLAG_DEFAULT_INPUT,     \
      "Input complex image")                                                    \
    X(".outname", &param_outname, FPTYPE_STREAMNAME, 1, FPFLAG_DEFAULT_INPUT,    \
      "Output complex image")                                                   \
    X(".pupilscale", &param_pupilscale, FPTYPE_FLOAT64, 1,                      \
      FPFLAG_DEFAULT_INPUT, "Pupil sampling scale [m/pixel]")                    \
    X(".distance", &param_z, FPTYPE_FLOAT64, 1, FPFLAG_DEFAULT_INPUT,           \
      "Propagation distance [m]")                                               \
    X(".lambda", &param_lambda, FPTYPE_FLOAT64, 1, FPFLAG_DEFAULT_INPUT,         \
      "Optical wavelength [m]")

/* ================================================================
 * 4.  COMPUTATION LOGIC
 * ============================================================= */

/**
 * fpsexec - Execute Fresnel diffractive propagation
 *
 * Return: RETURN_SUCCESS on success, error code otherwise.
 */
static MILK_HOT errno_t fpsexec(void)
{
    Fresnel_propagate_wavefront(param_inname, param_outname, param_pupilscale,
                                param_z, param_lambda);

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
 * CLIADDCMD_milkatmturb__wfprop_fresnel_FPS - Register command with CLI
 *
 * Return: RETURN_SUCCESS on success.
 */
errno_t CLIADDCMD_milkatmturb__wfprop_fresnel_FPS(void)
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
