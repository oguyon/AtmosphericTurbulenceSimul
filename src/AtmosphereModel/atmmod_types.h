// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmmod_types.h
 * @brief   Internal constants, global state declarations, and helpers for AtmosphereModel
 *
 * Declares the 10000-bin (0 to 100 km at 10 m vertical resolution) chemical
 * density and thermodynamic arrays, per-species Refractive Index and Absorption
 * (RIA) lookup tables, and internal helper prototypes shared across Level-1
 * AtmosphereModel translation units.
 */

#ifndef ATMMOD_TYPES_H
#define ATMMOD_TYPES_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ATMMOD_NB_BINS   10000        // Vertical altitude bins (10 m step, 0..100 km)
#define ATMMOD_LOSCHMIDT 2.6867805e25 // Number density [m^-3] at STP (101.325 kPa, 273.15 K)

/* Shared site and meteorological configuration parameters */
extern float ZenithAngle;        // Observation zenith angle [rad]
extern int   TimeDayOfYear;      // Day of year [1..366]
extern float TimeLocalSolarTime; // Local apparent solar time [h]
extern float SiteLat;            // Observatory latitude [deg]
extern float SiteLong;           // Observatory longitude [deg]
extern float SiteAlt;            // Observatory altitude above sea level [m]
extern float CO2_ppm;            // Tropospheric CO2 mixing ratio [ppm]
extern int   SiteTPauto;         // 1 = use NRLMSISE-00 site T/P, 0 = match SiteTemp/SitePress
extern float SiteTemp;           // Manual site temperature [K]
extern float SitePress;          // Manual site pressure [atm]
extern float SiteH2OMethod;      // H2O profile mode (1: TPW+scaleH, 2: RH+scaleH, 3: RH+TPW)
extern float SiteTPW;            // Total precipitable water vapor [mm]
extern float SiteRH;             // Ground relative humidity [%]
extern float SitePWSH;           // Precipitable water exponential scale height [m]
extern float alpha1H2O;          // Fitted H2O scaling factor [dimensionless]

/* Vertical atmospheric profile state (ATMMOD_NB_BINS elements each) */
extern int    initAtmosphereModel; // 1 once profile arrays are allocated
extern int    stdAtmModel_loaded;  // 1 once a profile has been built or loaded
extern float *densN2;              // N2 number density profile [cm^-3]
extern float *densO2;              // O2 number density profile [cm^-3]
extern float *densAr;              // Ar number density profile [cm^-3]
extern float *densH2O;             // H2O number density profile [cm^-3]
extern float *densCO2;             // CO2 number density profile [cm^-3]
extern float *densNe;              // Ne number density profile [cm^-3]
extern float *densHe;              // He number density profile [cm^-3]
extern float *densCH4;             // CH4 number density profile [cm^-3]
extern float *densKr;              // Kr number density profile [cm^-3]
extern float *densH2;              // H2 number density profile [cm^-3]
extern float *densO3;              // O3 number density profile [cm^-3]
extern float *densN;               // Atomic N number density profile [cm^-3]
extern float *densO;               // Atomic O number density profile [cm^-3]
extern float *densH;               // Atomic H number density profile [cm^-3]
extern float *denstot;             // Total gas number density profile [cm^-3]
extern float *density;             // Mass density profile [g/cm^3]
extern float *temperature;         // Kinetic temperature profile [K]
extern float *pressure;            // Atmospheric pressure profile [atm]
extern float *RH;                  // Fractional relative humidity profile [dimensionless]
extern float  dens0;               // Sea-level total number density [cm^-3]
extern double v_ABSCOEFF;          // Last evaluated linear absorption coefficient [m^-1]
extern double v_TRANSM;            // Last integrated ray-path transmission [0..1]

/* Refractive index and absorption (RIA) wavelength table state */
extern long   RIA_NBpts;      // Sample count from most recent RIA header read
extern double RIA_lambda_min; // Minimum wavelength of most recent RIA table [m]
extern double RIA_lambda_max; // Maximum wavelength of most recent RIA table [m]

extern int     initRIA_N2;
extern long    RIA_N2_NBpts;
extern double *RIA_N2_lambda;
extern double *RIA_N2_ri;
extern double *RIA_N2_abs;

extern int     initRIA_O2;
extern long    RIA_O2_NBpts;
extern double *RIA_O2_lambda;
extern double *RIA_O2_ri;
extern double *RIA_O2_abs;

extern int     initRIA_Ar;
extern long    RIA_Ar_NBpts;
extern double *RIA_Ar_lambda;
extern double *RIA_Ar_ri;
extern double *RIA_Ar_abs;

extern int     initRIA_H2O;
extern long    RIA_H2O_NBpts;
extern double *RIA_H2O_lambda;
extern double *RIA_H2O_ri;
extern double *RIA_H2O_abs;

extern int     initRIA_CO2;
extern long    RIA_CO2_NBpts;
extern double *RIA_CO2_lambda;
extern double *RIA_CO2_ri;
extern double *RIA_CO2_abs;

extern int     initRIA_Ne;
extern long    RIA_Ne_NBpts;
extern double *RIA_Ne_lambda;
extern double *RIA_Ne_ri;
extern double *RIA_Ne_abs;

extern int     initRIA_He;
extern long    RIA_He_NBpts;
extern double *RIA_He_lambda;
extern double *RIA_He_ri;
extern double *RIA_He_abs;

extern int     initRIA_CH4;
extern long    RIA_CH4_NBpts;
extern double *RIA_CH4_lambda;
extern double *RIA_CH4_ri;
extern double *RIA_CH4_abs;

extern int     initRIA_Kr;
extern long    RIA_Kr_NBpts;
extern double *RIA_Kr_lambda;
extern double *RIA_Kr_ri;
extern double *RIA_Kr_abs;

extern int     initRIA_H2;
extern long    RIA_H2_NBpts;
extern double *RIA_H2_lambda;
extern double *RIA_H2_ri;
extern double *RIA_H2_abs;

extern int     initRIA_O3;
extern long    RIA_O3_NBpts;
extern double *RIA_O3_lambda;
extern double *RIA_O3_ri;
extern double *RIA_O3_abs;

extern int     initRIA_N;
extern long    RIA_N_NBpts;
extern double *RIA_N_lambda;
extern double *RIA_N_ri;
extern double *RIA_N_abs;

extern int     initRIA_O;
extern long    RIA_O_NBpts;
extern double *RIA_O_lambda;
extern double *RIA_O_ri;
extern double *RIA_O_abs;

extern int     initRIA_H;
extern long    RIA_H_NBpts;
extern double *RIA_H_lambda;
extern double *RIA_H_ri;
extern double *RIA_H_abs;

/* Cached wavelength search indices per species (-1 when uninitialized) */
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

/**
 * atmmod_reset_precomp_indices - Reset cached wavelength lookup indices to -1
 */
void atmmod_reset_precomp_indices(void);

/**
 * atmmod_load_species_ria - Allocate and load RIA dispersion arrays for a species
 * @name:       Species identifier string (e.g. "N2", "H2O").
 * @path:       Path to RIA ASCII data file.
 * @init_flag:  Output flag set to 1 on successful load, 0 otherwise.
 * @nbpts_out:  Output number of wavelength samples loaded.
 * @lambda_out: Output pointer to allocated wavelength array [m].
 * @ri_out:     Output pointer to allocated refractive index array.
 * @abs_out:    Output pointer to allocated absorption coefficient array [m^-1].
 */
void atmmod_load_species_ria(
    const char *name,
    const char *path,
    int        *init_flag,
    long       *nbpts_out,
    double    **lambda_out,
    double    **ri_out,
    double    **abs_out);

#endif // ATMMOD_TYPES_H
