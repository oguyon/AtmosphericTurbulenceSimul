// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_mkwfs_FPS.c
 * @brief   FPS V2 compute unit for atmospheric turbulence wavefront series generation
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
    .fps_name         = "atmturb_mkwfs",
    .cmdkey           = "mkwfs_fps",
    .description      = "Generate atmospheric turbulence wavefront series",
    .description_long = "Generates multi-layer atmospheric turbulence wavefront series using "
                        "precomputed or generated phase screens and parameters from configuration file."
};

/* ================================================================
 * 2.  LOCAL PARAMETER VARIABLES
 * ============================================================= */

static float   param_slambda   = 1650.0f;
static int32_t param_precision = 0;
static char    param_conffile[FUNCTION_PARAMETER_STRMAXLEN] = "WFsim.conf";

/* ================================================================
 * 3.  UNIFIED PARAMETER TABLE (X-Macro)
 * ============================================================= */

#define FPS_PARAMS(X)                                                                            \
    X(".slambda", &param_slambda, FPTYPE_FLOAT32, 1, FPFLAG_DEFAULT_INPUT, "Wavelength [um or nm]") \
    X(".precision", &param_precision, FPTYPE_INT32, 1, FPFLAG_DEFAULT_INPUT,                    \
      "Precision (0=single, 1=double)")                                                          \
    X(".conffile", &param_conffile, FPTYPE_FILENAME, 0, FPFLAG_DEFAULT_INPUT,                    \
      "Configuration file (default: WFsim.conf)")

/* ================================================================
 * 4.  COMPUTATION LOGIC
 * ============================================================= */

static MILK_HOT errno_t fpsexec(void)
{
    if (param_conffile[0] != '\0')
    {
        AtmosphericTurbulence_change_configuration_file(param_conffile);
    }

    float wavel = param_slambda;
    if (wavel > 20.0f)
    {
        // Convert from nm to um if > 20
        wavel *= 1e-3f;
    }

    make_AtmosphericTurbulence_wavefront_series(wavel, (long) param_precision);

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

errno_t CLIADDCMD_milkatmturb__atmturb_mkwfs_FPS(void)
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
