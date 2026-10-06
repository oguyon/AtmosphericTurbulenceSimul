// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmmod_types.h
 * @brief   Internal types, state declarations, and prototypes for AtmosphereModel
 */

#ifndef ATMMOD_TYPES_H
#define ATMMOD_TYPES_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ATMMOD_NB_BINS 10000
#define ATMMOD_LOSCHMIDT 2.6867805e25 // [m-3] at STP (101.325 kPa, 273.15 K)

// Shared site configuration variables
extern float ZenithAngle;
extern int TimeDayOfYear;
extern float TimeLocalSolarTime;
extern float SiteLat;
extern float SiteLong;
extern float SiteAlt;
extern float CO2_ppm;
extern int SiteTPauto;
extern float SiteTemp;
extern float SitePress;
extern float SiteH2OMethod;
extern float SiteTPW;
extern float SiteRH;
extern float SitePWSH;
extern float alpha1H2O;

// Atmospheric profile state
extern int initAtmosphereModel;
extern float *densN2;
extern float *densO2;
extern float *densAr;
extern float *densH2O;
extern float *densCO2;
extern float *densNe;
extern float *densHe;
extern float *densCH4;
extern float *densKr;
extern float *densH2;
extern float *densO3;
extern float *densN;
extern float *densO;
extern float *densH;
extern float *denstot;
extern float *density;
extern float *temperature;
extern float *pressure;
extern float *RH;
extern float dens0;
extern double v_ABSCOEFF;
extern double v_TRANSM;

// Refractive indices and absorption coefficient table pointers
extern long RIA_NBpts;
extern double RIA_lambda_min;
extern double RIA_lambda_max;

extern int initRIA_N2;
extern long RIA_N2_NBpts;
extern double *RIA_N2_lambda;
extern double *RIA_N2_ri;
extern double *RIA_N2_abs;

extern int initRIA_O2;
extern long RIA_O2_NBpts;
extern double *RIA_O2_lambda;
extern double *RIA_O2_ri;
extern double *RIA_O2_abs;

extern int initRIA_Ar;
extern long RIA_Ar_NBpts;
extern double *RIA_Ar_lambda;
extern double *RIA_Ar_ri;
extern double *RIA_Ar_abs;

extern int initRIA_H2O;
extern long RIA_H2O_NBpts;
extern double *RIA_H2O_lambda;
extern double *RIA_H2O_ri;
extern double *RIA_H2O_abs;

extern int initRIA_CO2;
extern long RIA_CO2_NBpts;
extern double *RIA_CO2_lambda;
extern double *RIA_CO2_ri;
extern double *RIA_CO2_abs;

extern int initRIA_Ne;
extern long RIA_Ne_NBpts;
extern double *RIA_Ne_lambda;
extern double *RIA_Ne_ri;
extern double *RIA_Ne_abs;

extern int initRIA_He;
extern long RIA_He_NBpts;
extern double *RIA_He_lambda;
extern double *RIA_He_ri;
extern double *RIA_He_abs;

extern int initRIA_CH4;
extern long RIA_CH4_NBpts;
extern double *RIA_CH4_lambda;
extern double *RIA_CH4_ri;
extern double *RIA_CH4_abs;

extern int initRIA_Kr;
extern long RIA_Kr_NBpts;
extern double *RIA_Kr_lambda;
extern double *RIA_Kr_ri;
extern double *RIA_Kr_abs;

extern int initRIA_H2;
extern long RIA_H2_NBpts;
extern double *RIA_H2_lambda;
extern double *RIA_H2_ri;
extern double *RIA_H2_abs;

extern int initRIA_O3;
extern long RIA_O3_NBpts;
extern double *RIA_O3_lambda;
extern double *RIA_O3_ri;
extern double *RIA_O3_abs;

extern int initRIA_N;
extern long RIA_N_NBpts;
extern double *RIA_N_lambda;
extern double *RIA_N_ri;
extern double *RIA_N_abs;

extern int initRIA_O;
extern long RIA_O_NBpts;
extern double *RIA_O_lambda;
extern double *RIA_O_ri;
extern double *RIA_O_abs;

extern int initRIA_H;
extern long RIA_H_NBpts;
extern double *RIA_H_lambda;
extern double *RIA_H_ri;
extern double *RIA_H_abs;

extern int lliprecompN2;
extern int lliprecompO2;
extern int lliprecompAr;
extern int lliprecompH2O;
extern int lliprecompCO2;
extern int lliprecompNe;
extern int lliprecompHe;
extern int lliprecompCH4;
extern int lliprecompKr;
extern int lliprecompH2;
extern int lliprecompO3;
extern int lliprecompO;
extern int lliprecompN;
extern int lliprecompH;

// Internal helpers
void atmmod_reset_precomp_indices(void);
void atmmod_load_species_ria(const char *name, const char *path, int *init_flag,
                             long *nbpts_out, double **lambda_out, double **ri_out,
                             double **abs_out);

#endif // ATMMOD_TYPES_H
