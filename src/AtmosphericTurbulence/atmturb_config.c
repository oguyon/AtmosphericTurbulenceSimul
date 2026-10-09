// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_config.c
 * @brief   Atmospheric turbulence configuration file reading and air state equations
 */

#include <unistd.h>
#include "AtmosphericTurbulence.h"
#include "atmturb_types.h"

// Physical constants
double C_me = 9.10938291e-31;     // electron mass [kg]
double C_e0 = 8.854187817620e-12; // vacuum permittivity [F.m-1]
double C_Na = 6.0221413e23;        // Avogadro constant [mol-1]
double C_e = 1.60217657e-19;       // electron charge [C]
double C_ls = 2.686777447e25;      // Loschmidt constant [m-3]
double rhocoeff = 1.0;

// Configuration defaults and state
char CONFFILE[200] = "";

float CONF_LAMBDA = 0.5e-6f;
float CONF_SEEING = 0.8f;
char CONF_TURBULENCE_PROF_FILE[200] = "turbul.prof";
float CONF_ZANGLE = 0.0f;
float CONF_PARALLACTIC_ANGLE = 0.0f;
float CONF_SITE_ALT = -1.0f;
uint64_t CONF_SEED = 1ULL;
float CONF_SOURCE_Xpos = 0.0f;
float CONF_SOURCE_Ypos = 0.0f;

int CONF_WFOUTPUT = 1;
char CONF_WF_FILE_PREFIX[200] = "wf";
char CONF_WF_PHASE_NAME[100] = "outarraypha";
char CONF_WF_AMPL_NAME[100]  = "outarrayamp";
int CONF_SHM_OUTPUT = 0;
int CONF_STREAM_MODE = 0;

int CONF_MAKE_SWAVEFRONT = 0;
int CONF_SWF_WRITE2DISK = 0;
char CONF_SWF_FILE_PREFIX[200] = "swf";
int CONF_SHM_SOUTPUT = 0;
char CONF_SHM_SPREFIX[100] = "shmswf";
int CONF_SHM_SOUTPUTM = 0;

long CONF_WFsize = 512;
float CONF_PUPIL_SCALE = 0.01f;

int CONF_ATMWF_REALTIME = 0;
double CONF_ATMWF_REALTIMEFACTOR = 1.0;
float CONF_WFTIME_STEP = 0.001f;
float CONF_TIME_SPAN = 1.0f;
long CONF_NB_TSPAN = 1000;
long CONF_SIMTDELAY = 0;
int CONF_WAITFORSEM = 0;
char CONF_WAITSEMIMNAME[100] = "";

int CONF_SKIP_EXISTING = 0;
long CONF_WF_RAW_SIZE = 512;
long CONF_MASTER_SIZE = 4096;
int CONF_OVERSAMPLE = 2;
int CONF_INTERP = 1;
int CONF_LOWFREQ = 1;
int CONF_ROLLING = 1;
float CONF_BOIL_TIME = 0.0f;

int CONF_FRESNEL_PROPAGATION = 0;
int CONF_WAVEFRONT_AMPLITUDE = 0;
float CONF_FRESNEL_PROPAGATION_BIN = 100.0f;
int CONF_FRESNEL_RYTOV_SEC_EXACT = 0;
int CONF_FRESNEL_GUARD_PIX = 0;

/**
 * AtmosphericTurbulence_change_configuration_file - Set active configuration file path
 * @fname: Path to configuration file.
 *
 * Return: 0 on success.
 */
int AtmosphericTurbulence_change_configuration_file(const char *fname)
{
    snprintf(CONFFILE, sizeof(CONFFILE), "%s", fname);
    return 0;
}

/**
 * Z_Air - Compute compressibility factor for moist air
 * @P: Pressure in Pascals.
 * @T: Temperature in Celsius.
 * @RH: Relative humidity in percent.
 *
 * Return: Compressibility factor Z.
 */
double Z_Air(double P, double T, double RH)
{
    const double P0 = 101325.0;
    const double A = 1.2378847e-5, B = -1.9121316e-2, C = 33.93711047, D = -6.3431645e3;
    const double alpha = 1.00062, beta = 3.14e-8, gamma = 5.6e-7;
    const double a0 = 1.58123e-6, a1 = -2.9331e-8, a2 = 1.1043e-10, b1 = -2.051e-8;
    const double c0 = 1.9898e-4, c1 = -2.376e-6, d = 1.83e-11, e = -0.765e-8;

    double TK = T + 273.15;
    double Psv = exp(A * TK * TK + B * TK + C + D / TK);
    double f = alpha + beta * P + gamma * T * T;
    double xv = RH * 0.01 * f * Psv / P;

    double Z = 1.0 - P / TK * (a0 + a1 * T + a2 * T * T + (c0 + b1 * T) * xv +
                              (c0 + c1 * T) * xv * xv) +
               P * P / TK / TK * (d + e * xv * xv);

    double TK0 = 273.15;
    double xv0 = 0.0;
    double Z0 = 1.0 - P0 / TK0 * (a0 + a1 * 0.0 + a2 * 0.0 + (c0) * xv0 +
                                 (c0) * xv0 * xv0) +
                P0 * P0 / TK0 / TK0 * (d + e * xv0 * xv0);

    double rhocoeff1 = (101325.0 / P) * (TK / 273.15) * (Z / Z0);
    rhocoeff = rhocoeff1;

    return Z;
}

/**
 * Z_N2 - Compute compressibility factor for pure nitrogen
 * @P: Pressure in Pascals.
 * @T: Temperature in Kelvin.
 *
 * Return: Compressibility factor Z.
 */
double Z_N2(double P, double T)
{
    const double An = 4.446e-6;
    const double Bn = 6.4e-13;
    const double Cn = -1.07e-16;

    double tc = T - 273.15;
    double Z = 1.0 - 101325.0 * (P / 101325.0)
               * (0.449805 - 0.01177 * tc + 0.00006 * tc * tc) * 1e-8;
    double Z0 = 1.0 - 101325.0 * (0.449805) * 1e-8;

    double rhocoeff1 = (101325.0 / P) * (T / 273.15) * (Z / Z0);
    double rho0 = C_ls / C_Na;
    double rho = rho0 / rhocoeff1;

    rhocoeff = rhocoeff1 * (1.0 + Bn / An * rho0 + Cn / An * rho0 * rho0) /
               (1.0 + Bn / An * rho + Cn / An * rho * rho) /
               (1.0 + Bn / An * rho + Cn / An * rho * rho);

    return Z;
}

/**
 * atmturb_read_conf_turbulence - Read turbulence and atmosphere parameters
 */
static void atmturb_read_conf_turbulence(void)
{
    char keyword[200], content[200];

    snprintf(keyword, sizeof(keyword), "TURBULENCE_REF_WAVEL");
    read_config_parameter(CONFFILE, keyword, content);
    CONF_LAMBDA = (float)(atof(content) * 1e-6);

    snprintf(keyword, sizeof(keyword), "TURBULENCE_SEEING");
    read_config_parameter(CONFFILE, keyword, content);
    CONF_SEEING = (float)atof(content);

    snprintf(keyword, sizeof(keyword), "TURBULENCE_PROF_FILE");
    read_config_parameter(CONFFILE, keyword, CONF_TURBULENCE_PROF_FILE);

    snprintf(keyword, sizeof(keyword), "ZENITH_ANGLE");
    read_config_parameter(CONFFILE, keyword, content);
    CONF_ZANGLE = (float)atof(content);

    snprintf(keyword, sizeof(keyword), "PARALLACTIC_ANGLE");
    if (read_config_parameter_exists(CONFFILE, keyword) == 1)
    {
        read_config_parameter(CONFFILE, keyword, content);
        CONF_PARALLACTIC_ANGLE = (float)atof(content);
    }

    snprintf(keyword, sizeof(keyword), "SITE_ALT");
    if (read_config_parameter_exists(CONFFILE, keyword) == 1)
    {
        read_config_parameter(CONFFILE, keyword, content);
        CONF_SITE_ALT = (float)atof(content);
    }

    snprintf(keyword, sizeof(keyword), "SEED");
    if (read_config_parameter_exists(CONFFILE, keyword) == 1)
    {
        read_config_parameter(CONFFILE, keyword, content);
        CONF_SEED = (uint64_t)strtoull(content, NULL, 10);
    }

    snprintf(keyword, sizeof(keyword), "SOURCE_XPOS");
    read_config_parameter(CONFFILE, keyword, content);
    CONF_SOURCE_Xpos = (float)atof(content);

    snprintf(keyword, sizeof(keyword), "SOURCE_YPOS");
    read_config_parameter(CONFFILE, keyword, content);
    CONF_SOURCE_Ypos = (float)atof(content);
}

/**
 * atmturb_read_conf_output - Read output stream and file formatting settings
 */
static void atmturb_read_conf_output(void)
{
    char keyword[200], content[200];

    snprintf(keyword, sizeof(keyword), "WFOUTPUT");
    if (read_config_parameter_exists(CONFFILE, keyword) == 1)
    {
        read_config_parameter(CONFFILE, keyword, content);
        CONF_WFOUTPUT = atoi(content);
    }

    snprintf(keyword, sizeof(keyword), "WF_FILE_PREFIX");
    read_config_parameter(CONFFILE, keyword, CONF_WF_FILE_PREFIX);

    snprintf(keyword, sizeof(keyword), "SHM_OUTPUT");
    if (read_config_parameter_exists(CONFFILE, keyword) == 1)
    {
        read_config_parameter(CONFFILE, keyword, content);
        CONF_SHM_OUTPUT = atoi(content);
    }

    snprintf(keyword, sizeof(keyword), "STREAM_MODE");
    if (read_config_parameter_exists(CONFFILE, keyword) == 1)
    {
        read_config_parameter(CONFFILE, keyword, content);
        CONF_STREAM_MODE = atoi(content);
    }

    snprintf(keyword, sizeof(keyword), "MAKE_SWAVEFRONT");
    read_config_parameter(CONFFILE, keyword, content);
    CONF_MAKE_SWAVEFRONT = atoi(content);

    snprintf(keyword, sizeof(keyword), "SWF_WRITE2DISK");
    if (read_config_parameter_exists(CONFFILE, keyword) == 1)
    {
        read_config_parameter(CONFFILE, keyword, content);
        CONF_SWF_WRITE2DISK = atoi(content);
        if (CONF_WFOUTPUT == 2)
        {
            CONF_WFOUTPUT = 0;
        }
        if (CONF_SWF_WRITE2DISK == 1)
        {
            CONF_WFOUTPUT |= 2;
        }
    }

    snprintf(keyword, sizeof(keyword), "SWF_FILE_PREFIX");
    read_config_parameter(CONFFILE, keyword, CONF_SWF_FILE_PREFIX);

    snprintf(keyword, sizeof(keyword), "SHM_SOUTPUT");
    if (read_config_parameter_exists(CONFFILE, keyword) == 1)
    {
        read_config_parameter(CONFFILE, keyword, content);
        CONF_SHM_SOUTPUT = atoi(content);
    }

    snprintf(keyword, sizeof(keyword), "SHM_SPREFIX");
    read_config_parameter(CONFFILE, keyword, CONF_SHM_SPREFIX);

    snprintf(keyword, sizeof(keyword), "SHM_SOUTPUTM");
    read_config_parameter(CONFFILE, keyword, content);
    CONF_SHM_SOUTPUTM = atoi(content);

    snprintf(keyword, sizeof(keyword), "PUPIL_SCALE");
    read_config_parameter(CONFFILE, keyword, content);
    CONF_PUPIL_SCALE = (float)atof(content);

    snprintf(keyword, sizeof(keyword), "WFsize");
    read_config_parameter(CONFFILE, keyword, content);
    CONF_WFsize = atol(content);
}

/**
 * atmturb_read_conf_timing - Read simulation timing and synchronization settings
 */
static void atmturb_read_conf_timing(void)
{
    char keyword[200], content[200];

    snprintf(keyword, sizeof(keyword), "REALTIME");
    read_config_parameter(CONFFILE, keyword, content);
    CONF_ATMWF_REALTIME = atoi(content);

    snprintf(keyword, sizeof(keyword), "REALTIMEFACTOR");
    read_config_parameter(CONFFILE, keyword, content);
    CONF_ATMWF_REALTIMEFACTOR = atof(content);

    snprintf(keyword, sizeof(keyword), "WFTIME_STEP");
    read_config_parameter(CONFFILE, keyword, content);
    CONF_WFTIME_STEP = (float)atof(content);

    snprintf(keyword, sizeof(keyword), "TIME_SPAN");
    read_config_parameter(CONFFILE, keyword, content);
    CONF_TIME_SPAN = (float)atof(content);

    snprintf(keyword, sizeof(keyword), "NB_TSPAN");
    read_config_parameter(CONFFILE, keyword, content);
    CONF_NB_TSPAN = atol(content);

    snprintf(keyword, sizeof(keyword), "SIMTDELAY");
    read_config_parameter(CONFFILE, keyword, content);
    CONF_SIMTDELAY = atol(content);

    snprintf(keyword, sizeof(keyword), "WAITFORSEM");
    read_config_parameter(CONFFILE, keyword, content);
    CONF_WAITFORSEM = atoi(content);

    snprintf(keyword, sizeof(keyword), "WAITSEMIMNAME");
    read_config_parameter(CONFFILE, keyword, CONF_WAITSEMIMNAME);
}

/**
 * atmturb_read_conf_screen - Read phase screen sizing and interpolation settings
 */
static void atmturb_read_conf_screen(void)
{
    char keyword[200], content[200];

    snprintf(keyword, sizeof(keyword), "SKIP_EXISTING");
    CONF_SKIP_EXISTING = (read_config_parameter_exists(CONFFILE, keyword) == 1) ? 1 : 0;

    snprintf(keyword, sizeof(keyword), "WF_RAW_SIZE");
    read_config_parameter(CONFFILE, keyword, content);
    CONF_WF_RAW_SIZE = atol(content);

    snprintf(keyword, sizeof(keyword), "MASTER_SIZE");
    read_config_parameter(CONFFILE, keyword, content);
    CONF_MASTER_SIZE = atol(content);

    snprintf(keyword, sizeof(keyword), "MASTER_OVERSAMPLE");
    if (read_config_parameter_exists(CONFFILE, keyword) == 1)
    {
        read_config_parameter(CONFFILE, keyword, content);
        CONF_OVERSAMPLE = atoi(content);
    }
    else
    {
        snprintf(keyword, sizeof(keyword), "OVERSAMPLE");
        if (read_config_parameter_exists(CONFFILE, keyword) == 1)
        {
            read_config_parameter(CONFFILE, keyword, content);
            CONF_OVERSAMPLE = atoi(content);
        }
    }

    snprintf(keyword, sizeof(keyword), "INTERP");
    if (read_config_parameter_exists(CONFFILE, keyword) == 1)
    {
        read_config_parameter(CONFFILE, keyword, content);
        CONF_INTERP = atoi(content);
    }

    snprintf(keyword, sizeof(keyword), "LOWFREQ");
    if (read_config_parameter_exists(CONFFILE, keyword) == 1)
    {
        read_config_parameter(CONFFILE, keyword, content);
        CONF_LOWFREQ = atoi(content);
    }

    snprintf(keyword, sizeof(keyword), "ROLLING");
    if (read_config_parameter_exists(CONFFILE, keyword) == 1)
    {
        read_config_parameter(CONFFILE, keyword, content);
        CONF_ROLLING = atoi(content);
    }

    snprintf(keyword, sizeof(keyword), "BOIL_TIME");
    if (read_config_parameter_exists(CONFFILE, keyword) == 1)
    {
        read_config_parameter(CONFFILE, keyword, content);
        CONF_BOIL_TIME = (float) atof(content);
    }
}

/**
 * atmturb_read_conf_pupil_files - Read optional pupil mask FITS files
 */
static void atmturb_read_conf_pupil_files(void)
{
    char keyword[200], content[200];

    snprintf(keyword, sizeof(keyword), "PUPIL_AMPL_FILE");
    if (read_config_parameter_exists(CONFFILE, keyword) == 1)
    {
        read_config_parameter(CONFFILE, keyword, content);
        load_fits(content, "ST_pa", 1);
    }

    snprintf(keyword, sizeof(keyword), "PUPIL_PHA_FILE");
    if (read_config_parameter_exists(CONFFILE, keyword) == 1)
    {
        read_config_parameter(CONFFILE, keyword, content);
        load_fits(content, "ST_pp", 1);
    }
}

/**
 * atmturb_read_conf_modes - Read compute modes and pupil masks
 */
static void atmturb_read_conf_modes(void)
{
    char keyword[200], content[200];

    atmturb_read_conf_screen();

    snprintf(keyword, sizeof(keyword), "WAVEFRONT_AMPLITUDE");
    read_config_parameter(CONFFILE, keyword, content);
    CONF_WAVEFRONT_AMPLITUDE = atoi(content);

    if (CONF_WAVEFRONT_AMPLITUDE == 1)
    {
        snprintf(keyword, sizeof(keyword), "FRESNEL_PROPAGATION");
        read_config_parameter(CONFFILE, keyword, content);
        CONF_FRESNEL_PROPAGATION = atoi(content);
        if (CONF_FRESNEL_PROPAGATION == 3)
        {
            CONF_FRESNEL_PROPAGATION     = 2;
            CONF_FRESNEL_RYTOV_SEC_EXACT = 1;
        }
    }
    else
    {
        CONF_FRESNEL_PROPAGATION = 0;
    }

    snprintf(keyword, sizeof(keyword), "FRESNEL_PROPAGATION_BIN");
    read_config_parameter(CONFFILE, keyword, content);
    CONF_FRESNEL_PROPAGATION_BIN = (float)atof(content);

    snprintf(keyword, sizeof(keyword), "FRESNEL_RYTOV_EXACT");
    if (read_config_parameter_exists(CONFFILE, keyword) == 1 ||
        read_config_parameter_exists(CONFFILE, "FRESNEL_RYTOV_SEC_EX") == 1)
    {
        read_config_parameter(CONFFILE, "FRESNEL_RYTOV_SEC_EXACT", content);
        if (strcmp(content, "-") == 0)
        {
            read_config_parameter(CONFFILE, "FRESNEL_RYTOV_EXACT", content);
        }
        CONF_FRESNEL_RYTOV_SEC_EXACT = atoi(content);
    }
    else if (CONF_FRESNEL_RYTOV_SEC_EXACT != 1)
    {
        CONF_FRESNEL_RYTOV_SEC_EXACT = 0;
    }

    snprintf(keyword, sizeof(keyword), "FRESNEL_GUARD_PIX");
    if (read_config_parameter_exists(CONFFILE, keyword) == 1)
    {
        read_config_parameter(CONFFILE, keyword, content);
        CONF_FRESNEL_GUARD_PIX = atoi(content);
        if (CONF_FRESNEL_GUARD_PIX < 0)
        {
            CONF_FRESNEL_GUARD_PIX = 0;
        }
    }

    atmturb_read_conf_pupil_files();
}

/**
 * atmturb_write_default_profile - Create standard 7-layer turbulence profile
 * @fname: Path to output profile file.
 *
 * Return: 0 on success, -1 on failure.
 */
static int atmturb_write_default_profile(const char *fname)
{
    FILE *fp = fopen(fname, "w");
    if (fp == NULL)
    {
        return -1;
    }

    fprintf(fp, "# altitude(m)   relativeCN2     speed(m/s)      direction(rad)\n\n");
    fprintf(fp, " 4215     5.32        6.5     1.47\n");
    fprintf(fp, " 4230     1.47        6.55    1.57\n");
    fprintf(fp, " 4349     1.08        6.6     1.67\n");
    fprintf(fp, " 5007     2.11        6.7     1.77\n");
    fprintf(fp, "12000     1.83       22.0     3.10\n");
    fprintf(fp, "16200     1.48        9.5     3.20\n");
    fprintf(fp, "23701     0.697       5.6     3.30\n");
    fclose(fp);

    printf("[milkatmturb] Created default atmospheric profile \"%s\"\n", fname);
    return 0;
}

/**
 * atmturb_write_default_config - Create standard default simulation configuration
 * @fname: Path to output configuration file.
 *
 * Return: 0 on success, -1 on failure.
 */
static int atmturb_write_default_config(const char *fname)
{
    FILE *fp = fopen(fname, "w");
    if (fp == NULL)
    {
        return -1;
    }

    fprintf(fp, "# Atmospheric Turbulence Simulation Default Configuration\n\n");
    fprintf(fp, "TURBULENCE_REF_WAVEL      0.500000\n");
    fprintf(fp, "TURBULENCE_SEEING         0.600000\n");
    fprintf(fp, "TURBULENCE_PROF_FILE      turbul.prof\n");
    fprintf(fp, "ZENITH_ANGLE              0.0\n");
    fprintf(fp, "SOURCE_XPOS               0.0\n");
    fprintf(fp, "SOURCE_YPOS               0.0\n");
    fprintf(fp, "WFOUTPUT                  1\n");
    fprintf(fp, "WF_FILE_PREFIX            wf\n");
    fprintf(fp, "SHM_OUTPUT                1\n");
    fprintf(fp, "MAKE_SWAVEFRONT           0\n");
    fprintf(fp, "SLAMBDA                   1.650000\n");
    fprintf(fp, "SWF_WRITE2DISK            0\n");
    fprintf(fp, "SWF_FILE_PREFIX           swf\n");
    fprintf(fp, "SHM_SOUTPUT               0\n");
    fprintf(fp, "SHM_SPREFIX               shmswf\n");
    fprintf(fp, "SHM_SOUTPUTM              0\n");
    fprintf(fp, "WFsize                    256\n");
    fprintf(fp, "PUPIL_SCALE               0.040000\n");
    fprintf(fp, "REALTIME                  0\n");
    fprintf(fp, "REALTIMEFACTOR            1.0\n");
    fprintf(fp, "WFTIME_STEP               0.001000\n");
    fprintf(fp, "TIME_SPAN                 0.050000\n");
    fprintf(fp, "NB_TSPAN                  1\n");
    fprintf(fp, "SIMTDELAY                 0\n");
    fprintf(fp, "WAITFORSEM                0\n");
    fprintf(fp, "WAITSEMIMNAME             wfsimwait\n");
    fprintf(fp, "SKIP_EXISTING             0\n");
    fprintf(fp, "WF_RAW_SIZE               256\n");
    fprintf(fp, "MASTER_SIZE               4096\n");
    fprintf(fp, "MASTER_OVERSAMPLE         2\n");
    fprintf(fp, "INTERP                    1\n");
    fprintf(fp, "LOWFREQ                   1\n");
    fprintf(fp, "ROLLING                   1\n");
    fprintf(fp, "BOIL_TIME                 0.0\n");
    fprintf(fp, "WAVEFRONT_AMPLITUDE       0\n");
    fprintf(fp, "FRESNEL_PROPAGATION       0\n");
    fprintf(fp, "FRESNEL_PROPAGATION_BIN   100.0\n");
    fprintf(fp, "FRESNEL_GUARD_PIX         0\n");
    fclose(fp);

    printf("[milkatmturb] Created default simulation configuration \"%s\"\n", fname);
    if (access("turbul.prof", R_OK) != 0)
    {
        atmturb_write_default_profile("turbul.prof");
    }
    return 0;
}

/**
 * AtmosphericTurbulence_ReadConf - Read full simulation configuration from CONFFILE
 *
 * Return: 0 on success, -1 on failure.
 */
int AtmosphericTurbulence_ReadConf(void)
{
    if (CONFFILE[0] == '\0')
    {
        return 0;
    }

    if (access(CONFFILE, R_OK) != 0)
    {
        if (strcmp(CONFFILE, "WFsim.conf") == 0)
        {
            printf("[milkatmturb] Notice: \"%s\" not found in working directory.\n", CONFFILE);
            if (atmturb_write_default_config(CONFFILE) != 0)
            {
                printf("[milkatmturb] Warning: could not write \"%s\", using built-in defaults.\n",
                       CONFFILE);
                return 0;
            }
        }
        else
        {
            printf("ERROR: Configuration file \"%s\" not found.\n", CONFFILE);
            return -1;
        }
    }

    atmturb_read_conf_turbulence();
    atmturb_read_conf_output();
    atmturb_read_conf_timing();
    atmturb_read_conf_modes();
    return 0;
}
