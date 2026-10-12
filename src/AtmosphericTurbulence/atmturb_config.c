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
int CONF_STREAM_CBSIZE = 0;

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
int CONF_FRESNEL_RYTOV_MOISAN = -1;
int CONF_FRESNEL_RYTOV_ZINTERP = 0;
int CONF_FRESNEL_REFRACT_PATH = 0;
int CONF_FRESNEL_SCINT_WEIGHT = 0;

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
    const double a0 = 1.58123e-6, a1 = -2.9331e-8, a2 = 1.1043e-10;
    const double b0 = 5.707e-6, b1 = -2.051e-8;
    const double c0 = 1.9898e-4, c1 = -2.376e-6, d = 1.83e-11, e = -0.765e-8;

    double TK = T + 273.15;
    double Psv = exp(A * TK * TK + B * TK + C + D / TK);
    double f = alpha + beta * P + gamma * T * T;
    double xv = RH * 0.01 * f * Psv / P;

    double Z = 1.0 - P / TK * (a0 + a1 * T + a2 * T * T + (b0 + b1 * T) * xv +
                              (c0 + c1 * T) * xv * xv) +
               P * P / TK / TK * (d + e * xv * xv);

    double TK0 = 273.15;
    double xv0 = 0.0;
    double Z0 = 1.0 - P0 / TK0 * (a0 + a1 * 0.0 + a2 * 0.0 + (b0) * xv0 +
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
 * atmturb_read_conf_str - Read string parameter from configuration file
 * @conffile: Path to configuration file.
 * @key: Parameter keyword.
 * @out: Destination buffer for parameter string.
 * @out_sz: Size of destination buffer.
 * @default_val: Default value if key is not found.
 *
 * Return: 1 if key was found and valid, 0 if default_val was used.
 */
static int atmturb_read_conf_str(
    const char *conffile,
    const char *key,
    char       *out,
    size_t      out_sz,
    const char *default_val)
{
    if (out == NULL || out_sz == 0)
    {
        return 0;
    }
    if (default_val != NULL && out != default_val)
    {
        snprintf(out, out_sz, "%.*s", (int) (out_sz - 1), default_val);
    }
    if (conffile == NULL || conffile[0] == '\0' || key == NULL || key[0] == '\0')
    {
        return 0;
    }

    FILE *fp = fopen(conffile, "r");
    if (fp == NULL)
    {
        return 0;
    }

    char line[512];
    char keyw[128];
    char val[256];
    int found = 0;

    while (fgets(line, sizeof(line), fp) != NULL)
    {
        char *p = line;
        while (*p == ' ' || *p == '\t')
        {
            p++;
        }
        if (*p == '#' || *p == '\0' || *p == '\n')
        {
            continue;
        }

        if (sscanf(p, "%127s %255s", keyw, val) == 2)
        {
            if (strcmp(keyw, key) == 0)
            {
                if (strcmp(val, "-") != 0 && strcmp(val, "---") != 0)
                {
                    snprintf(out, out_sz, "%.*s", (int) (out_sz - 1), val);
                    found = 1;
                }
                break;
            }
        }
    }

    fclose(fp);
    return found;
}

static long atmturb_read_conf_long(
    const char *conffile,
    const char *key,
    long        default_val)
{
    char buf[64];
    if (atmturb_read_conf_str(conffile, key, buf, sizeof(buf), NULL))
    {
        char *end = NULL;
        long v = strtol(buf, &end, 10);
        if (end != buf)
        {
            return v;
        }
    }
    return default_val;
}

static int atmturb_read_conf_int(
    const char *conffile,
    const char *key,
    int         default_val)
{
    return (int) atmturb_read_conf_long(conffile, key, (long) default_val);
}

static double atmturb_read_conf_double(
    const char *conffile,
    const char *key,
    double      default_val)
{
    char buf[64];
    if (atmturb_read_conf_str(conffile, key, buf, sizeof(buf), NULL))
    {
        char *end = NULL;
        double v = strtod(buf, &end);
        if (end != buf)
        {
            return v;
        }
    }
    return default_val;
}

static float atmturb_read_conf_float(
    const char *conffile,
    const char *key,
    float       default_val)
{
    return (float) atmturb_read_conf_double(conffile, key, (double) default_val);
}

static uint64_t atmturb_read_conf_uint64(
    const char *conffile,
    const char *key,
    uint64_t    default_val)
{
    char buf[64];
    if (atmturb_read_conf_str(conffile, key, buf, sizeof(buf), NULL))
    {
        char *end = NULL;
        unsigned long long v = strtoull(buf, &end, 10);
        if (end != buf)
        {
            return (uint64_t) v;
        }
    }
    return default_val;
}

/**
 * atmturb_read_conf_turbulence - Read turbulence and atmosphere parameters
 */
static void atmturb_read_conf_turbulence(void)
{
    float wavel_um = atmturb_read_conf_float(CONFFILE, "TURBULENCE_REF_WAVEL", -1.0f);
    if (wavel_um > 0.0f)
    {
        CONF_LAMBDA = wavel_um * 1e-6f;
    }
    CONF_SEEING = atmturb_read_conf_float(CONFFILE, "TURBULENCE_SEEING", CONF_SEEING);
    atmturb_read_conf_str(CONFFILE, "TURBULENCE_PROF_FILE", CONF_TURBULENCE_PROF_FILE,
                          sizeof(CONF_TURBULENCE_PROF_FILE), CONF_TURBULENCE_PROF_FILE);
    CONF_ZANGLE = atmturb_read_conf_float(CONFFILE, "ZENITH_ANGLE", CONF_ZANGLE);
    CONF_PARALLACTIC_ANGLE = atmturb_read_conf_float(CONFFILE, "PARALLACTIC_ANGLE",
                                                     CONF_PARALLACTIC_ANGLE);
    CONF_SITE_ALT = atmturb_read_conf_float(CONFFILE, "SITE_ALT", CONF_SITE_ALT);
    CONF_SEED = atmturb_read_conf_uint64(CONFFILE, "SEED", CONF_SEED);
    CONF_SOURCE_Xpos = atmturb_read_conf_float(CONFFILE, "SOURCE_XPOS", CONF_SOURCE_Xpos);
    CONF_SOURCE_Ypos = atmturb_read_conf_float(CONFFILE, "SOURCE_YPOS", CONF_SOURCE_Ypos);
}

/**
 * atmturb_read_conf_output - Read output stream and file formatting settings
 */
static void atmturb_read_conf_output(void)
{
    CONF_WFOUTPUT       = atmturb_read_conf_int(CONFFILE, "WFOUTPUT", CONF_WFOUTPUT);
    atmturb_read_conf_str(CONFFILE, "WF_FILE_PREFIX", CONF_WF_FILE_PREFIX,
                          sizeof(CONF_WF_FILE_PREFIX), CONF_WF_FILE_PREFIX);
    CONF_SHM_OUTPUT     = atmturb_read_conf_int(CONFFILE, "SHM_OUTPUT", CONF_SHM_OUTPUT);
    CONF_STREAM_MODE    = atmturb_read_conf_int(CONFFILE, "STREAM_MODE", CONF_STREAM_MODE);
    CONF_STREAM_CBSIZE  = atmturb_read_conf_int(CONFFILE, "STREAM_CBSIZE", CONF_STREAM_CBSIZE);
    CONF_MAKE_SWAVEFRONT = atmturb_read_conf_int(CONFFILE, "MAKE_SWAVEFRONT", CONF_MAKE_SWAVEFRONT);

    int swf_disk = atmturb_read_conf_int(CONFFILE, "SWF_WRITE2DISK", -1);
    if (swf_disk >= 0)
    {
        CONF_SWF_WRITE2DISK = swf_disk;
        if (CONF_WFOUTPUT == 2)
        {
            CONF_WFOUTPUT = 0;
        }
        if (CONF_SWF_WRITE2DISK == 1)
        {
            CONF_WFOUTPUT |= 2;
        }
    }
    atmturb_read_conf_str(CONFFILE, "SWF_FILE_PREFIX", CONF_SWF_FILE_PREFIX,
                          sizeof(CONF_SWF_FILE_PREFIX), CONF_SWF_FILE_PREFIX);
    CONF_SHM_SOUTPUT    = atmturb_read_conf_int(CONFFILE, "SHM_SOUTPUT", CONF_SHM_SOUTPUT);
    atmturb_read_conf_str(CONFFILE, "SHM_SPREFIX", CONF_SHM_SPREFIX,
                          sizeof(CONF_SHM_SPREFIX), CONF_SHM_SPREFIX);
    CONF_SHM_SOUTPUTM   = atmturb_read_conf_int(CONFFILE, "SHM_SOUTPUTM", CONF_SHM_SOUTPUTM);
    CONF_PUPIL_SCALE    = atmturb_read_conf_float(CONFFILE, "PUPIL_SCALE", CONF_PUPIL_SCALE);
    CONF_WFsize         = atmturb_read_conf_long(CONFFILE, "WFsize", CONF_WFsize);
}

/**
 * atmturb_read_conf_timing - Read simulation timing and synchronization settings
 */
static void atmturb_read_conf_timing(void)
{
    CONF_ATMWF_REALTIME       = atmturb_read_conf_int(CONFFILE, "REALTIME", CONF_ATMWF_REALTIME);
    CONF_ATMWF_REALTIMEFACTOR = atmturb_read_conf_double(CONFFILE, "REALTIMEFACTOR",
                                                         CONF_ATMWF_REALTIMEFACTOR);
    CONF_WFTIME_STEP          = atmturb_read_conf_float(CONFFILE, "WFTIME_STEP", CONF_WFTIME_STEP);
    CONF_TIME_SPAN            = atmturb_read_conf_float(CONFFILE, "TIME_SPAN", CONF_TIME_SPAN);
    CONF_NB_TSPAN             = atmturb_read_conf_long(CONFFILE, "NB_TSPAN", CONF_NB_TSPAN);
    CONF_SIMTDELAY            = atmturb_read_conf_long(CONFFILE, "SIMTDELAY", CONF_SIMTDELAY);
    CONF_WAITFORSEM           = atmturb_read_conf_int(CONFFILE, "WAITFORSEM", CONF_WAITFORSEM);
    atmturb_read_conf_str(CONFFILE, "WAITSEMIMNAME", CONF_WAITSEMIMNAME,
                          sizeof(CONF_WAITSEMIMNAME), CONF_WAITSEMIMNAME);
}

/**
 * atmturb_read_conf_screen - Read phase screen sizing and interpolation settings
 */
static void atmturb_read_conf_screen(void)
{
    char tmp[64];
    CONF_SKIP_EXISTING = atmturb_read_conf_str(CONFFILE, "SKIP_EXISTING", tmp,
                                              sizeof(tmp), NULL) ? 1 : 0;
    CONF_WF_RAW_SIZE   = atmturb_read_conf_long(CONFFILE, "WF_RAW_SIZE", CONF_WF_RAW_SIZE);
    CONF_MASTER_SIZE   = atmturb_read_conf_long(CONFFILE, "MASTER_SIZE", CONF_MASTER_SIZE);
    CONF_OVERSAMPLE    = atmturb_read_conf_int(
        CONFFILE, "OVERSAMPLE",
        atmturb_read_conf_int(CONFFILE, "MASTER_OVERSAMPLE", CONF_OVERSAMPLE));
    CONF_INTERP        = atmturb_read_conf_int(CONFFILE, "INTERP", CONF_INTERP);
    CONF_LOWFREQ       = atmturb_read_conf_int(CONFFILE, "LOWFREQ", CONF_LOWFREQ);
    CONF_ROLLING       = atmturb_read_conf_int(CONFFILE, "ROLLING", CONF_ROLLING);
    CONF_BOIL_TIME     = atmturb_read_conf_float(CONFFILE, "BOIL_TIME", CONF_BOIL_TIME);
}

/**
 * atmturb_read_conf_pupil_files - Read optional pupil mask FITS files
 */
static void atmturb_read_conf_pupil_files(void)
{
    char path[256];
    if (atmturb_read_conf_str(CONFFILE, "PUPIL_AMPL_FILE", path, sizeof(path), NULL))
    {
        load_fits(path, "ST_pa", 1);
    }
    if (atmturb_read_conf_str(CONFFILE, "PUPIL_PHA_FILE", path, sizeof(path), NULL))
    {
        load_fits(path, "ST_pp", 1);
    }
}

/**
 * atmturb_read_conf_modes - Read compute modes and pupil masks
 */
static void atmturb_read_conf_modes(void)
{
    atmturb_read_conf_screen();

    CONF_WAVEFRONT_AMPLITUDE = atmturb_read_conf_int(CONFFILE, "WAVEFRONT_AMPLITUDE",
                                                     CONF_WAVEFRONT_AMPLITUDE);
    if (CONF_WAVEFRONT_AMPLITUDE == 1)
    {
        CONF_FRESNEL_PROPAGATION = atmturb_read_conf_int(CONFFILE, "FRESNEL_PROPAGATION",
                                                         CONF_FRESNEL_PROPAGATION);
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

    CONF_FRESNEL_PROPAGATION_BIN = atmturb_read_conf_float(CONFFILE, "FRESNEL_PROPAGATION_BIN",
                                                           CONF_FRESNEL_PROPAGATION_BIN);
    CONF_FRESNEL_RYTOV_SEC_EXACT = atmturb_read_conf_int(
        CONFFILE, "FRESNEL_RYTOV_SEC_EXACT",
        atmturb_read_conf_int(CONFFILE, "FRESNEL_RYTOV_EXACT",
        atmturb_read_conf_int(CONFFILE, "FRESNEL_RYTOV_SEC_EX", CONF_FRESNEL_RYTOV_SEC_EXACT)));

    CONF_FRESNEL_GUARD_PIX = atmturb_read_conf_int(CONFFILE, "FRESNEL_GUARD_PIX",
                                                   CONF_FRESNEL_GUARD_PIX);
    if (CONF_FRESNEL_GUARD_PIX < 0)
    {
        CONF_FRESNEL_GUARD_PIX = 0;
    }
    CONF_FRESNEL_RYTOV_MOISAN = atmturb_read_conf_int(CONFFILE, "FRESNEL_RYTOV_MOISAN", -1);

    CONF_FRESNEL_RYTOV_ZINTERP = atmturb_read_conf_int(
        CONFFILE, "FRESNEL_RYTOV_ZINTERP",
        atmturb_read_conf_int(CONFFILE, "FRESNEL_RYTOV_ZINT", CONF_FRESNEL_RYTOV_ZINTERP));
    CONF_FRESNEL_REFRACT_PATH = atmturb_read_conf_int(CONFFILE, "FRESNEL_REFRACT_PATH",
                                                      CONF_FRESNEL_REFRACT_PATH);
    CONF_FRESNEL_SCINT_WEIGHT = atmturb_read_conf_int(CONFFILE, "FRESNEL_SCINT_WEIGHT",
                                                      CONF_FRESNEL_SCINT_WEIGHT);

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
    fprintf(fp, "FRESNEL_RYTOV_MOISAN     -1\n");
    fprintf(fp, "FRESNEL_RYTOV_ZINT        0\n");
    fprintf(fp, "FRESNEL_REFRACT_PATH      0\n");
    fprintf(fp, "FRESNEL_SCINT_WEIGHT      0\n");
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
