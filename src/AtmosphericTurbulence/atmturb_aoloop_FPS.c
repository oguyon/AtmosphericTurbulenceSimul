// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_aoloop_FPS.c
 * @brief   FPS V2 compute unit for closed-loop AO simulation and science PSF synthesis
 */

#ifdef MILK_NO_CLI
#    include "CLIcore_standalone.h"
#else
#    include "CLIcore.h"
#endif
#include "fps.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "AtmosphericTurbulence/AtmosphericTurbulence.h"
#include "atmturb_aoloop_FPS.h"

/* ================================================================
 * 1.  FPS COMPONENT IDENTITY
 * ============================================================= */

static FPS_APP_INFO FPS_app_info = {
    .fps_name         = "atmturb_aoloop",
    .cmdkey           = "aoloop",
    .description      = "Run adaptive optics closed loop and synthesize science PSF",
    .description_long = "Applies leaky integrator or PID wavefront correction with latency "
                        "to a simulated wavefront series and synthesizes Nyquist-sampled "
                        "long-exposure science PSFs."
};

/* ================================================================
 * 2.  LOCAL PARAMETER VARIABLES
 * ============================================================= */

static char    param_wfname[FUNCTION_PARAMETER_STRMAXLEN]  = "outarraypha";
static char    param_psfname[FUNCTION_PARAMETER_STRMAXLEN] = "PSFcumul";
static char    param_fitsout[FUNCTION_PARAMETER_STRMAXLEN] = "";
static int32_t param_loop_mode = 1;
static double  param_gain      = 0.5;
static double  param_leak      = 0.001;
static int32_t param_delay     = 1;
static double  param_scilambda = 1.65;
static double  param_teldiam   = 8.0;

/* ================================================================
 * 3.  UNIFIED PARAMETER TABLE (X-Macro)
 * ============================================================= */

#define FPS_PARAMS(X)                                                                             \
    X(".wfname", &param_wfname, FPTYPE_FILENAME, 1, FPFLAG_DEFAULT_INPUT,                        \
      "Input wavefront phase cube filename or stream name (default: outarraypha)")               \
    X(".psfname", &param_psfname, FPTYPE_STREAMNAME, 1, FPFLAG_DEFAULT_INPUT,                     \
      "Output PSF stream name (default: PSFcumul)")                                              \
    X(".fitsout", &param_fitsout, FPTYPE_STRING, 1, FPFLAG_DEFAULT_INPUT,                         \
      "Optional FITS output file (empty = do not save)")                                         \
    X(".loop_mode", &param_loop_mode, FPTYPE_INT32, 1, FPFLAG_DEFAULT_INPUT,                      \
      "Loop mode: 0 = Open Loop, 1 = Leaky Integrator, 2 = PID (default: 1)")                   \
    X(".gain", &param_gain, FPTYPE_FLOAT64, 1, FPFLAG_DEFAULT_INPUT,                              \
      "Loop gain (default: 0.5)")                                                                \
    X(".leak", &param_leak, FPTYPE_FLOAT64, 1, FPFLAG_DEFAULT_INPUT,                              \
      "Integrator leak factor (default: 0.001)")                                                 \
    X(".delay", &param_delay, FPTYPE_INT32, 1, FPFLAG_DEFAULT_INPUT,                              \
      "Loop latency delay in frames (default: 1)")                                               \
    X(".scilambda", &param_scilambda, FPTYPE_FLOAT64, 1, FPFLAG_DEFAULT_INPUT,                    \
      "Science wavelength in microns (default: 1.65)")                                           \
    X(".teldiam", &param_teldiam, FPTYPE_FLOAT64, 1, FPFLAG_DEFAULT_INPUT,                        \
      "Telescope primary mirror diameter in meters (default: 8.0)")

/* ================================================================
 * 4.  COMPUTATION LOGIC
 * ============================================================= */

/**
 * fpsexec - Execute closed-loop AO simulation
 *
 * Return: RETURN_SUCCESS on success, RETURN_FAILURE otherwise.
 */
static MILK_HOT errno_t fpsexec(void)
{
    atmturb_ao_params_t p;
    atmturb_ao_init_params(&p);

    p.in_wfname     = param_wfname;
    p.out_psfname   = param_psfname;
    p.out_fitsname  = (param_fitsout[0] != '\0') ? param_fitsout : NULL;
    p.loop_mode     = (int) param_loop_mode;
    p.gain          = param_gain;
    p.leak          = param_leak;
    p.loop_delay    = (int) param_delay;
    p.lambda_sci_m  = param_scilambda * 1e-6;
    p.tel_diam_m    = param_teldiam;

    atmturb_ao_results_t res;
    if (atmturb_ao_sim_run(&p, &res) != 0)
    {
        return RETURN_FAILURE;
    }

    printf("[milkatmturb] AO simulation complete: Strehl = %.4f (first = %.4f, last = %.4f)\n",
           res.strehl_cumul, res.strehl_first, res.strehl_last);
    fflush(stdout);

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
 * CLIADDCMD_milkatmturb__atmturb_aoloop_FPS - Register command with CLI
 *
 * Return: RETURN_SUCCESS on success.
 */
errno_t CLIADDCMD_milkatmturb__atmturb_aoloop_FPS(void)
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
