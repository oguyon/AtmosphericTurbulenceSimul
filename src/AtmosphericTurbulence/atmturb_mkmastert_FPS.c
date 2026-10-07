// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_mkmastert_FPS.c
 * @brief   FPS V2 compute unit for master phase screens generation
 */

#ifdef MILK_NO_CLI
#    include "CLIcore_standalone.h"
#else
#    include "CLIcore.h"
#endif
#include "fps.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "AtmosphericTurbulence/AtmosphericTurbulence.h"
#include "atmturb_mkmastert_FPS.h"

/* ================================================================
 * 1.  FPS COMPONENT IDENTITY
 * ============================================================= */

static FPS_APP_INFO FPS_app_info = {
    .fps_name         = "atmturb_mkmastert",
    .cmdkey           = "mkmastert",
    .description      = "Generate master Kolmogorov/von Karman turbulence screens",
    .description_long = "Generates a pair of normalized master phase screens in Fourier domain "
                        "with outer and inner scale cutoff filtering."
};

/* ================================================================
 * 2.  LOCAL PARAMETER VARIABLES
 * ============================================================= */

static int32_t param_size       = 2048;
static float   param_outerscale = 100.0f;
static float   param_innerscale = 1.0f;
static int32_t param_precision  = 0;
static char    param_screen0[FUNCTION_PARAMETER_STRMAXLEN] = "turbm00_p0";
static char    param_screen1[FUNCTION_PARAMETER_STRMAXLEN] = "turbm00_p1";
static int64_t param_seed       = 1;
static char    param_fitsout0[FUNCTION_PARAMETER_STRMAXLEN] = "";
static char    param_fitsout1[FUNCTION_PARAMETER_STRMAXLEN] = "";

/* ================================================================
 * 3.  UNIFIED PARAMETER TABLE (X-Macro)
 * ============================================================= */

#define FPS_PARAMS(X)                                                          \
    X(".size", &param_size, FPTYPE_INT32, 1, FPFLAG_DEFAULT_INPUT,             \
      "Screen grid dimension [pixels]")                                        \
    X(".outerscale", &param_outerscale, FPTYPE_FLOAT32, 1,                     \
      FPFLAG_DEFAULT_INPUT, "Outer scale in grid units [pixels]")              \
    X(".innerscale", &param_innerscale, FPTYPE_FLOAT32, 1,                     \
      FPFLAG_DEFAULT_INPUT, "Inner scale in grid units [pixels]")              \
    X(".precision", &param_precision, FPTYPE_INT32, 1, FPFLAG_DEFAULT_INPUT,   \
      "Precision (0=single, 1=double)")                                        \
    X(".screen0", &param_screen0, FPTYPE_STREAMNAME, 1, FPFLAG_DEFAULT_INPUT,  \
      "Output screen 0 image name")                                            \
    X(".screen1", &param_screen1, FPTYPE_STREAMNAME, 1, FPFLAG_DEFAULT_INPUT,  \
      "Output screen 1 image name")                                            \
    X(".seed", &param_seed, FPTYPE_INT64, 1, FPFLAG_DEFAULT_INPUT,             \
      "RNG seed (0 = time-based)")                                             \
    X(".fitsout0", &param_fitsout0, FPTYPE_STRING, 1, FPFLAG_DEFAULT_INPUT,    \
      "Optional FITS output file for screen 0 (empty = do not save)")          \
    X(".fitsout1", &param_fitsout1, FPTYPE_STRING, 1, FPFLAG_DEFAULT_INPUT,    \
      "Optional FITS output file for screen 1 (empty = do not save)")

/* ================================================================
 * 4.  COMPUTATION LOGIC
 * ============================================================= */

/**
 * fpsexec - Execute master phase screen generation
 *
 * Return: RETURN_SUCCESS on success, error code otherwise.
 */
static MILK_HOT errno_t fpsexec(void)
{
    make_master_turbulence_screen_seeded(param_screen0, param_screen1, (long) param_size,
                                         param_outerscale, param_innerscale,
                                         (long) param_precision, (uint64_t) param_seed);

    if (param_fitsout0[0] != '\0')
    {
        save_fits(param_screen0, param_fitsout0);
    }
    if (param_fitsout1[0] != '\0')
    {
        save_fits(param_screen1, param_fitsout1);
    }

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
 * CLIADDCMD_milkatmturb__atmturb_mkmastert_FPS - Register command with CLI
 *
 * Return: RETURN_SUCCESS on success.
 */
errno_t CLIADDCMD_milkatmturb__atmturb_mkmastert_FPS(void)
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
