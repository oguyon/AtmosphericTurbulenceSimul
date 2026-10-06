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
#include "atmturb_types.h"

#include <unistd.h>

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

static float   param_slambda     = 1650.0f;
static int32_t param_precision   = 0;
static int32_t param_wfsize      = 256;
static float   param_pupil_scale = 0.04f;
static float   param_seeing      = 0.6f;
static float   param_time_step   = 0.001f;
static float   param_time_span   = 0.05f;
static char    param_prof_file[FUNCTION_PARAMETER_STRMAXLEN] = "turbul.prof";
static int32_t param_master_size = 2048;
static int32_t param_save_fits   = 1;
static int32_t param_amplitude   = 0;
static int32_t param_fresnel     = 0;
static char    param_out_phase[FUNCTION_PARAMETER_STRMAXLEN] = "outarraypha";
static char    param_out_ampl[FUNCTION_PARAMETER_STRMAXLEN]  = "outarrayamp";
static char    param_conffile[FUNCTION_PARAMETER_STRMAXLEN]  = "WFsim.conf";

/* ================================================================
 * 3.  UNIFIED PARAMETER TABLE (X-Macro)
 * ============================================================= */

#define FPS_PARAMS(X)                                                                                  \
    X(".slambda", &param_slambda, FPTYPE_FLOAT32, 1, FPFLAG_DEFAULT_INPUT, "Wavelength [um or nm]")       \
    X(".precision", &param_precision, FPTYPE_INT32, 1, FPFLAG_DEFAULT_INPUT,                            \
      "Precision (0=single, 1=double)")                                                                \
    X(".wfsize", &param_wfsize, FPTYPE_INT32, 0, FPFLAG_DEFAULT_INPUT, "Output grid dimension [pix]")    \
    X(".pupil_scale", &param_pupil_scale, FPTYPE_FLOAT32, 0, FPFLAG_DEFAULT_INPUT,                    \
      "Pupil sampling scale [m/pix]")                                                                  \
    X(".seeing", &param_seeing, FPTYPE_FLOAT32, 0, FPFLAG_DEFAULT_INPUT,                               \
      "Zenith seeing at reference wavelength [arcsec]")                                                \
    X(".time_step", &param_time_step, FPTYPE_FLOAT32, 0, FPFLAG_DEFAULT_INPUT,                        \
      "Simulation time step between frames [s]")                                                       \
    X(".time_span", &param_time_span, FPTYPE_FLOAT32, 0, FPFLAG_DEFAULT_INPUT,                        \
      "Wavefront cube duration [s]")                                                                   \
    X(".prof_file", &param_prof_file, FPTYPE_FILENAME, 0, FPFLAG_DEFAULT_INPUT,                        \
      "Turbulence profile file (default: turbul.prof)")                                                \
    X(".master_size", &param_master_size, FPTYPE_INT32, 0, FPFLAG_DEFAULT_INPUT,                      \
      "Master phase screen dimension [pix]")                                                           \
    X(".save_fits", &param_save_fits, FPTYPE_INT32, 0, FPFLAG_DEFAULT_INPUT,                          \
      "Save FITS data cubes to disk (0=no, 1=yes)")                                                    \
    X(".amplitude", &param_amplitude, FPTYPE_INT32, 0, FPFLAG_DEFAULT_INPUT,                          \
      "Compute amplitude in addition to phase (0/1)")                                                  \
    X(".fresnel", &param_fresnel, FPTYPE_INT32, 0, FPFLAG_DEFAULT_INPUT,                              \
      "Diffractive Fresnel inter-layer propagation (0/1)")                                             \
    X(".out_phase", &param_out_phase, FPTYPE_STREAMNAME, 0, FPFLAG_DEFAULT_INPUT,                      \
      "Output phase stream name")                                                                      \
    X(".out_ampl", &param_out_ampl, FPTYPE_STREAMNAME, 0, FPFLAG_DEFAULT_INPUT,                        \
      "Output amplitude stream name")                                                                  \
    X(".conffile", &param_conffile, FPTYPE_FILENAME, 0, FPFLAG_DEFAULT_INPUT,                          \
      "Optional legacy WFsim.conf file")

/* ================================================================
 * 4.  COMPUTATION LOGIC
 * ============================================================= */

/**
 * atmturb_mkwfs_sync_to_conf - Copy FPS parameters to simulation state
 */
static void atmturb_mkwfs_sync_to_conf(void)
{
    CONF_WFsize              = (long)param_wfsize;
    CONF_PUPIL_SCALE         = param_pupil_scale;
    CONF_SEEING              = param_seeing;
    CONF_WFTIME_STEP         = param_time_step;
    CONF_TIME_SPAN           = param_time_span;
    CONF_WFOUTPUT            = (int)param_save_fits;
    CONF_WAVEFRONT_AMPLITUDE = (int)param_amplitude;
    CONF_FRESNEL_PROPAGATION = (int)param_fresnel;
    CONF_MASTER_SIZE         = (long)param_master_size;

    if (param_prof_file[0] != '\0')
    {
        strncpy(CONF_TURBULENCE_PROF_FILE, param_prof_file,
                sizeof(CONF_TURBULENCE_PROF_FILE) - 1);
        CONF_TURBULENCE_PROF_FILE[sizeof(CONF_TURBULENCE_PROF_FILE) - 1] = '\0';
    }

    if (param_out_phase[0] != '\0')
    {
        strncpy(CONF_WF_PHASE_NAME, param_out_phase, sizeof(CONF_WF_PHASE_NAME) - 1);
        CONF_WF_PHASE_NAME[sizeof(CONF_WF_PHASE_NAME) - 1] = '\0';
    }

    if (param_out_ampl[0] != '\0')
    {
        strncpy(CONF_WF_AMPL_NAME, param_out_ampl, sizeof(CONF_WF_AMPL_NAME) - 1);
        CONF_WF_AMPL_NAME[sizeof(CONF_WF_AMPL_NAME) - 1] = '\0';
    }
}

/**
 * atmturb_mkwfs_sync_from_conf - Update FPS parameters from simulation state
 */
static void atmturb_mkwfs_sync_from_conf(void)
{
    param_wfsize      = (int32_t)CONF_WFsize;
    param_pupil_scale = CONF_PUPIL_SCALE;
    param_seeing      = CONF_SEEING;
    param_time_step   = CONF_WFTIME_STEP;
    param_time_span   = CONF_TIME_SPAN;
    param_save_fits   = (int32_t)CONF_WFOUTPUT;
    param_amplitude   = (int32_t)CONF_WAVEFRONT_AMPLITUDE;
    param_fresnel     = (int32_t)CONF_FRESNEL_PROPAGATION;
    param_master_size = (int32_t)CONF_MASTER_SIZE;

    if (CONF_TURBULENCE_PROF_FILE[0] != '\0')
    {
        strncpy(param_prof_file, CONF_TURBULENCE_PROF_FILE, sizeof(param_prof_file) - 1);
        param_prof_file[sizeof(param_prof_file) - 1] = '\0';
    }
}

static MILK_HOT errno_t fpsexec(void)
{
    if (param_conffile[0] != '\0' && access(param_conffile, R_OK) == 0)
    {
        AtmosphericTurbulence_change_configuration_file(param_conffile);
        if (AtmosphericTurbulence_ReadConf() != 0)
        {
            return RETURN_FAILURE;
        }
        atmturb_mkwfs_sync_from_conf();
    }
    else if (param_conffile[0] != '\0' && strcmp(param_conffile, "WFsim.conf") != 0)
    {
        printf("ERROR: Configuration file \"%s\" not found.\n", param_conffile);
        return RETURN_FAILURE;
    }
    else
    {
        CONFFILE[0] = '\0';
        atmturb_mkwfs_sync_to_conf();
    }

    float wavel = param_slambda;
    if (wavel > 20.0f)
    {
        // Convert from nm to um if > 20
        wavel *= 1e-3f;
    }

    if (make_AtmosphericTurbulence_wavefront_series(wavel, (long) param_precision) != 0)
    {
        return RETURN_FAILURE;
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
