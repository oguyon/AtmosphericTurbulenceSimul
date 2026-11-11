// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_aoloop_FPS.c
 * @brief   MILK_RUNNABLE FPS compute unit for closed-loop AO simulation and PSF synthesis
 *
 * Exposes `atmturb.aoloop` and `milk-fpsexec-atmturb-aoloop`. Simulates open-loop,
 * leaky-integrator, or PID wavefront control with integer frame latency, synthesizes
 * Nyquist-sampled science PSFs, and exports cumulative Strehl in calculator variable
 * `strehl`.
 */

#include <stdio.h>

#include "milk_runnable.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "AtmosphericTurbulence/AtmosphericTurbulence.h"
#include "atmturb_aoloop_FPS.h"

static char    param_wfname[FUNCTION_PARAMETER_STRMAXLEN]  = "outarraypha";
static char    param_psfname[FUNCTION_PARAMETER_STRMAXLEN] = "PSFcumul";
static char    param_fitsout[FUNCTION_PARAMETER_STRMAXLEN] = "";
static int32_t param_loop_mode                             = 1;
static double  param_gain                                  = 0.5;
static double  param_leak                                  = 0.001;
static int32_t param_delay                                 = 1;
static double  param_scilambda                             = 1.65;
static double  param_teldiam                               = 8.0;
static double  param_reflambda                             = 0.50;
static double  param_pupil_scale                           = 0.0;
static double  param_timestep                              = 0.001;

#define ATMTURB_AOLOOP_PARAMS(X)                                                   \
    X(".wfname", param_wfname, FPTYPE_FILENAME, 1,                                 \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL,                                  \
      "Input wavefront phase cube filename or stream name (default: outarraypha)") \
    X(".psfname", param_psfname, FPTYPE_STREAMNAME, 1,                             \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL,                                  \
      "Output PSF stream name (default: PSFcumul)")                                \
    X(".fitsout", param_fitsout, FPTYPE_STRING, 1,                                 \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL,                                  \
      "Optional FITS output file (empty = do not save)")                           \
    X(".loop_mode", &param_loop_mode, FPTYPE_INT32, 1,                             \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL,                                  \
      "Loop mode: 0 = Open Loop, 1 = Leaky Integrator, 2 = PID (default: 1)")      \
    X(".gain", &param_gain, FPTYPE_FLOAT64, 1,                                     \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL, "Loop gain (default: 0.5)")      \
    X(".leak", &param_leak, FPTYPE_FLOAT64, 1,                                     \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL,                                  \
      "Integrator leak factor (default: 0.001)")                                   \
    X(".delay", &param_delay, FPTYPE_INT32, 1,                                     \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL,                                  \
      "Loop latency delay in frames (default: 1)")                                 \
    X(".scilambda", &param_scilambda, FPTYPE_FLOAT64, 1,                           \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL,                                  \
      "Science wavelength in microns (default: 1.65)")                             \
    X(".teldiam", &param_teldiam, FPTYPE_FLOAT64, 1,                               \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL,                                  \
      "Telescope primary mirror diameter in meters (default: 8.0)")                \
    X(".reflambda", &param_reflambda, FPTYPE_FLOAT64, 1,                           \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL,                                  \
      "Reference wavefront wavelength in microns (default: 0.50)")                 \
    X(".pupil_scale", &param_pupil_scale, FPTYPE_FLOAT64, 1,                       \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL,                                  \
      "Pupil pixel scale in meters/pixel (0 = auto from teldiam/grid)")            \
    X(".time_step", &param_timestep, FPTYPE_FLOAT64, 1,                            \
      FPFLAG_DEFAULT_INPUT | FPFLAG_CLI_OPTIONAL,                                  \
      "Simulation time step in seconds (default: 0.001)")

MILK_PARAMS_DEFINE(atmturb_aoloop, ATMTURB_AOLOOP_PARAMS);

/**
 * atmturb_aoloop_exec - Execute closed-loop AO simulation
 *
 * Return: RETURN_SUCCESS on success, RETURN_FAILURE otherwise.
 */
static MILK_HOT errno_t atmturb_aoloop_exec(void)
{
    printf("[milkatmturb] Outputs:\n");
    printf("[milkatmturb]   - stream: \"%s\" (2D/3D SHM PSF)\n", param_psfname);
    if (param_fitsout[0] != '\0')
    {
        printf("[milkatmturb]   - file:   \"%s\" (FITS)\n", param_fitsout);
    }
    fflush(stdout);

    atmturb_ao_params_t p;
    atmturb_ao_init_params(&p);

    p.in_wfname     = param_wfname;
    p.out_psfname   = param_psfname;
    p.out_fitsname  = (param_fitsout[0] != '\0') ? param_fitsout : NULL;
    p.loop_mode     = (int) param_loop_mode;
    p.gain          = param_gain;
    p.leak          = param_leak;
    p.loop_delay    = (int) param_delay;
    p.lambda_ref_m  = (param_reflambda > 0.0) ? (param_reflambda * 1e-6) : 0.50e-6;
    p.lambda_sci_m  = param_scilambda * 1e-6;
    p.tel_diam_m    = param_teldiam;
    p.pupil_scale_m = param_pupil_scale;
    p.time_step_s   = param_timestep;

    atmturb_ao_results_t res;
    if (atmturb_ao_sim_run(&p, &res) != 0)
    {
        return RETURN_FAILURE;
    }

    create_variable_ID("strehl", res.strehl_cumul);

    printf("[milkatmturb] AO simulation complete: Strehl = %.4f (first = %.4f, last = %.4f)\n",
           res.strehl_cumul, res.strehl_first, res.strehl_last);
    fflush(stdout);

    return RETURN_SUCCESS;
}

/**
 * atmturb_aoloop_compute - Execute AO simulation within procinfo loop
 * @self: Pointer to this MILK_RUNNABLE descriptor.
 *
 * Return: RETURN_SUCCESS on completion, or error code on failure.
 */
static MILK_HOT errno_t atmturb_aoloop_compute(
    const MILK_RUNNABLE *self)
{
    DEBUG_TRACE_FSTART();
    errno_t ret = RETURN_SUCCESS;

    MILK_PROCINFO_INIT(self)
    MILK_PROCINFO_LOOPSTART(self)
    ret = atmturb_aoloop_exec();
    MILK_PROCINFO_END(self)

    DEBUG_TRACE_FEXIT();
    return ret;
}

const MILK_RUNNABLE runnable_atmturb_aoloop = {
    .kind = MILK_RUNNABLE_FPS,
    .app  = {
        .fps_name         = "atmturb_aoloop",
        .cmdkey           = "aoloop",
        .description      = "Run adaptive optics closed loop and synthesize science PSF",
        .description_long = "Applies open-loop (.loop_mode 0), leaky integrator (.loop_mode 1), "
                            "or PID (.loop_mode 2) wavefront correction with frame latency to a "
                            "phase cube, synthesizes Nyquist-sampled science PSFs, and writes "
                            "long-exposure Strehl ratio to calculator variable `strehl`.",
        .outvars          = "strehl"
    },
    .srcfile = __FILE__,
    MILK_PARAMS_REF(atmturb_aoloop),
    .compute = atmturb_aoloop_compute,
};

EXPORT_FPSEXEC(milk, "atmturb-aoloop", runnable_atmturb_aoloop);
