// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmmod_config_parser.c
 * @brief   Atmospheric model configuration file parser and builder
 */

#include <unistd.h>
#include "CLIcore.h"
#include "COREMOD_tools/COREMOD_tools.h"
#include "AtmosphereModel.h"
#include "atmmod_types.h"

/**
 * atmmod_read_conf_str - Read string parameter from configuration file
 * @conffile: Path to configuration file.
 * @key: Parameter keyword.
 * @out: Destination buffer for parameter string.
 * @out_sz: Size of destination buffer.
 *
 * Return: 1 if key was found and valid, 0 otherwise.
 */
static int atmmod_read_conf_str(
    const char *conffile,
    const char *key,
    char       *out,
    size_t      out_sz)
{
    if (conffile == NULL || conffile[0] == '\0' || key == NULL || out == NULL || out_sz == 0)
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

/**
 * atmmod_read_conf_float - Read floating-point parameter with fallback default
 * @conffile: Path to configuration file.
 * @key: Parameter keyword.
 * @default_val: Default value if key is not found.
 *
 * Return: Parsed float or default_val.
 */
static float atmmod_read_conf_float(
    const char *conffile,
    const char *key,
    float       default_val)
{
    char buf[64];
    if (atmmod_read_conf_str(conffile, key, buf, sizeof(buf)))
    {
        char *end = NULL;
        float v = strtof(buf, &end);
        if (end != buf)
        {
            return v;
        }
    }
    return default_val;
}

/**
 * atmmod_read_conf_int - Read integer parameter with fallback default
 * @conffile: Path to configuration file.
 * @key: Parameter keyword.
 * @default_val: Default value if key is not found.
 *
 * Return: Parsed integer or default_val.
 */
static int atmmod_read_conf_int(
    const char *conffile,
    const char *key,
    int         default_val)
{
    char buf[64];
    if (atmmod_read_conf_str(conffile, key, buf, sizeof(buf)))
    {
        char *end = NULL;
        long v = strtol(buf, &end, 10);
        if (end != buf)
        {
            return (int) v;
        }
    }
    return default_val;
}

/**
 * atmmod_read_site_config - Read atmospheric and site parameters from file
 * @conffile: Path to configuration file.
 */
static void atmmod_read_site_config(const char *conffile)
{
    ZenithAngle        = atmmod_read_conf_float(conffile, "ZENITH_ANGLE", ZenithAngle);
    TimeDayOfYear      = atmmod_read_conf_int(conffile, "TIME_DAY_OF_YEAR", TimeDayOfYear);
    TimeLocalSolarTime = atmmod_read_conf_float(conffile, "TIME_LOCAL_SOLAR_TIME",
                                                TimeLocalSolarTime);
    SiteLat            = atmmod_read_conf_float(conffile, "SITE_LATITUDE", SiteLat);
    SiteLong           = atmmod_read_conf_float(conffile, "SITE_LONGITUDE", SiteLong);
    SiteAlt            = atmmod_read_conf_float(conffile, "SITE_ALT", SiteAlt);
    SiteTPauto         = atmmod_read_conf_int(conffile, "SITE_TP_AUTO", SiteTPauto);
    SiteTemp           = atmmod_read_conf_float(conffile, "SITE_TEMP", SiteTemp);
    SitePress          = atmmod_read_conf_float(conffile, "SITE_PRESS", SitePress);
    CO2_ppm            = atmmod_read_conf_float(conffile, "CO2_PPM", CO2_ppm);
    SiteH2OMethod      = atmmod_read_conf_float(conffile, "SITE_H2O_METHOD", SiteH2OMethod);
    SiteTPW            = atmmod_read_conf_float(conffile, "SITE_TPW", SiteTPW);
    SiteRH             = atmmod_read_conf_float(conffile, "SITE_RH", SiteRH);
    SitePWSH           = atmmod_read_conf_float(conffile, "SITE_PW_SCALEH", SitePWSH);
}

/**
 * atmmod_load_all_species_tables - Load optical dispersion tables for all species
 */
static void atmmod_load_all_species_tables(void)
{
    atmmod_load_species_ria("N2", "./RefractiveIndices/RIA_N2.dat", &initRIA_N2,
                            &RIA_N2_NBpts, &RIA_N2_lambda, &RIA_N2_ri, &RIA_N2_abs);
    atmmod_load_species_ria("O2", "./RefractiveIndices/RIA_O2.dat", &initRIA_O2,
                            &RIA_O2_NBpts, &RIA_O2_lambda, &RIA_O2_ri, &RIA_O2_abs);
    atmmod_load_species_ria("Ar", "./RefractiveIndices/RIA_Ar.dat", &initRIA_Ar,
                            &RIA_Ar_NBpts, &RIA_Ar_lambda, &RIA_Ar_ri, &RIA_Ar_abs);
    atmmod_load_species_ria("H2O", "./RefractiveIndices/RIA_H2O.dat", &initRIA_H2O,
                            &RIA_H2O_NBpts, &RIA_H2O_lambda, &RIA_H2O_ri, &RIA_H2O_abs);
    atmmod_load_species_ria("CO2", "./RefractiveIndices/RIA_CO2.dat", &initRIA_CO2,
                            &RIA_CO2_NBpts, &RIA_CO2_lambda, &RIA_CO2_ri, &RIA_CO2_abs);
    atmmod_load_species_ria("Ne", "./RefractiveIndices/RIA_Ne.dat", &initRIA_Ne,
                            &RIA_Ne_NBpts, &RIA_Ne_lambda, &RIA_Ne_ri, &RIA_Ne_abs);
    atmmod_load_species_ria("He", "./RefractiveIndices/RIA_He.dat", &initRIA_He,
                            &RIA_He_NBpts, &RIA_He_lambda, &RIA_He_ri, &RIA_He_abs);
    atmmod_load_species_ria("CH4", "./RefractiveIndices/RIA_CH4.dat", &initRIA_CH4,
                            &RIA_CH4_NBpts, &RIA_CH4_lambda, &RIA_CH4_ri, &RIA_CH4_abs);
    atmmod_load_species_ria("Kr", "./RefractiveIndices/RIA_Kr.dat", &initRIA_Kr,
                            &RIA_Kr_NBpts, &RIA_Kr_lambda, &RIA_Kr_ri, &RIA_Kr_abs);
    atmmod_load_species_ria("H2", "./RefractiveIndices/RIA_H2.dat", &initRIA_H2,
                            &RIA_H2_NBpts, &RIA_H2_lambda, &RIA_H2_ri, &RIA_H2_abs);
    atmmod_load_species_ria("O3", "./RefractiveIndices/RIA_O3.dat", &initRIA_O3,
                            &RIA_O3_NBpts, &RIA_O3_lambda, &RIA_O3_ri, &RIA_O3_abs);
    atmmod_load_species_ria("N", "./RefractiveIndices/RIA_N.dat", &initRIA_N,
                            &RIA_N_NBpts, &RIA_N_lambda, &RIA_N_ri, &RIA_N_abs);
    atmmod_load_species_ria("O", "./RefractiveIndices/RIA_O.dat", &initRIA_O,
                            &RIA_O_NBpts, &RIA_O_lambda, &RIA_O_ri, &RIA_O_abs);
    atmmod_load_species_ria("H", "./RefractiveIndices/RIA_H.dat", &initRIA_H,
                            &RIA_H_NBpts, &RIA_H_lambda, &RIA_H_ri, &RIA_H_abs);
}

/**
 * atmmod_write_site_dispersion - Write wavelength dispersion curve at site altitude
 * @fname: Output filename.
 */
static void atmmod_write_site_dispersion(const char *fname)
{
    FILE *fp = fopen(fname, "w");
    if (fp == NULL)
    {
        return;
    }

    for (double lam = 0.2e-6; lam < 2.0e-6; lam *= 1.0 + 1e-3)
    {
        double n = 1.0 + AtmosphereModel_stdAtmModel_N(SiteAlt, (float)lam, 0);
        fprintf(fp, "%.8g %.14f %.14f\n", lam, n, v_ABSCOEFF);
    }
    fclose(fp);
}

/**
 * AtmosphereModel_Create_from_CONF - Build model from configuration file
 * @CONFFILE: Path to configuration file.
 * @slambda: Secondary test wavelength in meters.
 *
 * Return: 0 on success.
 */
int AtmosphereModel_Create_from_CONF(
    const char *CONFFILE,
    float       slambda)
{
    if (access(CONFFILE, R_OK) != 0)
    {
        printf("ERROR: Configuration file \"%s\" not found.\n", CONFFILE);
        return -1;
    }

    atmmod_read_site_config(CONFFILE);
    atmmod_load_all_species_tables();

    printf("Building reference atmosphere model ...\n");
    AtmosphereModel_build_stdAtmModel("atm.txt");

    atmmod_write_site_dispersion("RindexSite.txt");

    FILE *fp = fopen("Refract.txt", "w");
    if (fp != NULL)
    {
        fclose(fp);
    }

    atmmod_reset_precomp_indices();

    AtmosphereModel_RefractionPath(0.55e-6, ZenithAngle, 1);
    int r = system("mv refractpath.txt refractpath_0550.txt");
    (void)r;

    AtmosphereModel_RefractionPath(slambda, ZenithAngle, 1);
    char command[200];
    snprintf(command, sizeof(command), "mv refractpath.txt refractpath_%04ld.txt",
             (long)(1e9 * slambda + 0.5));
    r = system(command);
    (void)r;

    return 0;
}
