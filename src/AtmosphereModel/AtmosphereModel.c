// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    AtmosphereModel.c
 * @brief   Atmospheric model module lifecycle, state allocation, and CLI registration
 */

#include "CLIcore.h"
#include "AtmosphereModel.h"
#include "atmmod_types.h"

extern DATA data;

// Site and atmospheric configuration parameters
float ZenithAngle = 0.0f;
int TimeDayOfYear = 0;
float TimeLocalSolarTime = 0.0f;
float SiteLat = 0.0f;
float SiteLong = 0.0f;
float SiteAlt = 0.0f;
float CO2_ppm = 390.0f;

int SiteTPauto = 1;
float SiteTemp = 273.15f;
float SitePress = 1.0f;

float SiteH2OMethod = 1.0f;
float SiteTPW = 1.0f;
float SiteRH = 20.0f;
float SitePWSH = 2000.0f;
float alpha1H2O = 0.0f;

// Profile arrays (densities in cm^-3, 10000 vertical bins of 10m each)
int initAtmosphereModel = 0;
float *densN2 = NULL;
float *densO2 = NULL;
float *densAr = NULL;
float *densH2O = NULL;
float *densCO2 = NULL;
float *densNe = NULL;
float *densHe = NULL;
float *densCH4 = NULL;
float *densKr = NULL;
float *densH2 = NULL;
float *densO3 = NULL;
float *densN = NULL;
float *densO = NULL;
float *densH = NULL;
float *denstot = NULL;

float *density = NULL;
float *temperature = NULL;
float *pressure = NULL;
float *RH = NULL;
float dens0 = 0.0f;

double v_ABSCOEFF = 0.0;
double v_TRANSM = 1.0;

/**
 * init_AtmosphereModel - Initialize module and allocate profile arrays
 *
 * Return: 0 on success.
 */
int init_AtmosphereModel(void)
{
    if (initAtmosphereModel == 0)
    {
        denstot = (float *)malloc(sizeof(float) * ATMMOD_NB_BINS);
        densN2 = (float *)malloc(sizeof(float) * ATMMOD_NB_BINS);
        densO2 = (float *)malloc(sizeof(float) * ATMMOD_NB_BINS);
        densAr = (float *)malloc(sizeof(float) * ATMMOD_NB_BINS);
        densH2O = (float *)malloc(sizeof(float) * ATMMOD_NB_BINS);
        densCO2 = (float *)malloc(sizeof(float) * ATMMOD_NB_BINS);
        densNe = (float *)malloc(sizeof(float) * ATMMOD_NB_BINS);
        densHe = (float *)malloc(sizeof(float) * ATMMOD_NB_BINS);
        densCH4 = (float *)malloc(sizeof(float) * ATMMOD_NB_BINS);
        densKr = (float *)malloc(sizeof(float) * ATMMOD_NB_BINS);
        densH2 = (float *)malloc(sizeof(float) * ATMMOD_NB_BINS);
        densO3 = (float *)malloc(sizeof(float) * ATMMOD_NB_BINS);
        densO = (float *)malloc(sizeof(float) * ATMMOD_NB_BINS);
        densN = (float *)malloc(sizeof(float) * ATMMOD_NB_BINS);
        densH = (float *)malloc(sizeof(float) * ATMMOD_NB_BINS);

        density = (float *)malloc(sizeof(float) * ATMMOD_NB_BINS);
        temperature = (float *)malloc(sizeof(float) * ATMMOD_NB_BINS);
        pressure = (float *)malloc(sizeof(float) * ATMMOD_NB_BINS);
        RH = (float *)malloc(sizeof(float) * ATMMOD_NB_BINS);

        initAtmosphereModel = 1;
    }

    return 0;
}
