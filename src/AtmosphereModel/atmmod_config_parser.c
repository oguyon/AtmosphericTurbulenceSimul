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
 * atmmod_read_site_config - Read atmospheric and site parameters from file
 * @conffile: Path to configuration file.
 */
static void atmmod_read_site_config(char *conffile)
{
    char keyword[200];
    char content[200];

    snprintf(keyword, sizeof(keyword), "ZENITH_ANGLE");
    read_config_parameter(conffile, keyword, content);
    ZenithAngle = (float)atof(content);

    snprintf(keyword, sizeof(keyword), "TIME_DAY_OF_YEAR");
    read_config_parameter(conffile, keyword, content);
    TimeDayOfYear = atoi(content);

    snprintf(keyword, sizeof(keyword), "TIME_LOCAL_SOLAR_TIME");
    read_config_parameter(conffile, keyword, content);
    TimeLocalSolarTime = (float)atof(content);

    snprintf(keyword, sizeof(keyword), "SITE_LATITUDE");
    read_config_parameter(conffile, keyword, content);
    SiteLat = (float)atof(content);

    snprintf(keyword, sizeof(keyword), "SITE_LONGITUDE");
    read_config_parameter(conffile, keyword, content);
    SiteLong = (float)atof(content);

    snprintf(keyword, sizeof(keyword), "SITE_ALT");
    read_config_parameter(conffile, keyword, content);
    SiteAlt = (float)atof(content);

    snprintf(keyword, sizeof(keyword), "SITE_TP_AUTO");
    read_config_parameter(conffile, keyword, content);
    SiteTPauto = atoi(content);

    snprintf(keyword, sizeof(keyword), "SITE_TEMP");
    read_config_parameter(conffile, keyword, content);
    SiteTemp = (float)atof(content);

    snprintf(keyword, sizeof(keyword), "SITE_PRESS");
    read_config_parameter(conffile, keyword, content);
    SitePress = (float)atof(content);

    snprintf(keyword, sizeof(keyword), "CO2_PPM");
    read_config_parameter(conffile, keyword, content);
    CO2_ppm = (float)atof(content);

    snprintf(keyword, sizeof(keyword), "SITE_H2O_METHOD");
    read_config_parameter(conffile, keyword, content);
    SiteH2OMethod = (float)atof(content);

    snprintf(keyword, sizeof(keyword), "SITE_TPW");
    read_config_parameter(conffile, keyword, content);
    SiteTPW = (float)atof(content);

    snprintf(keyword, sizeof(keyword), "SITE_RH");
    read_config_parameter(conffile, keyword, content);
    SiteRH = (float)atof(content);

    snprintf(keyword, sizeof(keyword), "SITE_PW_SCALEH");
    read_config_parameter(conffile, keyword, content);
    SitePWSH = (float)atof(content);
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

    for (double lam = 0.2e-6; lam < 2.0e-6; lam *= 1.0 + 1e-6)
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
int AtmosphereModel_Create_from_CONF(char *CONFFILE, float slambda)
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
