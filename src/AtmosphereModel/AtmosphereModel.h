// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    AtmosphereModel.h
 * @brief   Atmospheric model public interface (Level 1)
 *
 * Computes vertical atmospheric density, temperature, pressure, and humidity
 * profiles from sea level to 100 km via NRLMSISE-00, evaluates polychromatic
 * 14-species Lorentz-Lorenz refractivity (n - 1) and absorption, and traces
 * curved optical ray paths through spherical atmospheric shells.
 */

#ifndef ATMOSPHEREMODEL_H
#define ATMOSPHEREMODEL_H

#include "atmturb_c23_compat.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Observatory site coordinates shared across AtmosphereModel routines */
extern float SiteLat;  // Site geodetic latitude [deg]
extern float SiteLong; // Site geodetic longitude [deg]
extern float SiteAlt;  // Site altitude above sea level [m]

/**
 * init_AtmosphereModel - Initialize module and allocate profile arrays
 *
 * Return: 0 on success.
 */
ATMTURB_NODISCARD int init_AtmosphereModel(void);

/**
 * AirMixture_N - Compute refractive index minus 1 for a 14-species gas mixture
 * @lambda:   Optical wavelength [m].
 * @dens_N2:  N2 number density [cm^-3].
 * @dens_O2:  O2 number density [cm^-3].
 * @dens_Ar:  Ar number density [cm^-3].
 * @dens_H2O: H2O number density [cm^-3].
 * @dens_CO2: CO2 number density [cm^-3].
 * @dens_Ne:  Ne number density [cm^-3].
 * @dens_He:  He number density [cm^-3].
 * @dens_CH4: CH4 number density [cm^-3].
 * @dens_Kr:  Kr number density [cm^-3].
 * @dens_H2:  H2 number density [cm^-3].
 * @dens_O3:  O3 number density [cm^-3].
 * @dens_N:   Atomic N number density [cm^-3].
 * @dens_O:   Atomic O number density [cm^-3].
 * @dens_H:   Atomic H number density [cm^-3].
 *
 * Evaluates (n - 1) via the Lorentz-Lorenz molar refractivity relation using
 * loaded RIA tables (falling back to OpticsMaterials dispersion) and updates
 * the global absorption coefficient v_ABSCOEFF [m^-1].
 *
 * Return: Mixture refractivity (n - 1) [dimensionless].
 */
double AirMixture_N(
    double lambda,
    double dens_N2,
    double dens_O2,
    double dens_Ar,
    double dens_H2O,
    double dens_CO2,
    double dens_Ne,
    double dens_He,
    double dens_CH4,
    double dens_Kr,
    double dens_H2,
    double dens_O3,
    double dens_N,
    double dens_O,
    double dens_H);

/**
 * AtmosphereModel_stdAtmModel_N - Refractivity (n - 1) at specified altitude
 * @alt:    Altitude above sea level [m].
 * @lambda: Optical wavelength [m].
 * @mode:   Verbosity flag (1 for debug output, 0 for silent).
 *
 * Interpolates loaded vertical species profiles at @alt (or falls back to an
 * exponential Edlen approximation when no profile is loaded).
 *
 * Return: Air refractivity (n - 1) [dimensionless].
 */
float AtmosphereModel_stdAtmModel_N(
    float alt,
    float lambda,
    int   mode);

/**
 * AtmosphereModel_H2O_Saturation - Water vapor saturation pressure
 * @T: Temperature [K].
 *
 * Return: Saturation vapor pressure [Pa] using the IAPWS-95 formulation.
 */
double AtmosphereModel_H2O_Saturation(
    double T);

/**
 * AtmosphereModel_save_stdAtmModel - Save standard atmosphere profile to ASCII file
 * @fname: Output file path.
 *
 * Return: 0 on success, -1 on failure.
 */
int AtmosphereModel_save_stdAtmModel(
    const char *fname);

/**
 * AtmosphereModel_build_stdAtmModel - Generate atmosphere profile via NRLMSISE-00
 * @fname: Output file path to save the generated profile.
 *
 * Return: 0 on success, -1 on failure.
 */
int AtmosphereModel_build_stdAtmModel(
    const char *fname);

/**
 * AtmosphereModel_load_stdAtmModel - Load atmosphere profile from ASCII file
 * @fname: Input file path.
 *
 * Return: 0 on success, -1 on failure.
 */
int AtmosphereModel_load_stdAtmModel(
    const char *fname);

/**
 * AtmosphereModel_Create_from_CONF - Build site atmosphere model from config file
 * @CONFFILE: Path to key-value site configuration file.
 * @slambda:  Secondary reference wavelength for refraction ray trace [m].
 *
 * Return: 0 on success, -1 if @CONFFILE cannot be read.
 */
int AtmosphereModel_Create_from_CONF(
    const char *CONFFILE,
    float       slambda);

/**
 * AtmosphereModel_RefractionPath - Trace refracted ray path through atmosphere
 * @lambda:    Optical wavelength [m].
 * @Zangle:    Apparent zenith angle at ground [rad].
 * @WritePath: 1 to write trajectory to refractpath.txt, 0 otherwise.
 *
 * Iteratively solves for the upward launch angle matching @Zangle and records
 * total atmospheric transmission in v_TRANSM.
 *
 * Return: Atmospheric refraction deflection angle [arcsec].
 */
double AtmosphereModel_RefractionPath(
    double lambda,
    double Zangle,
    int    WritePath);

/**
 * ATMOSPHEREMODEL_loadRIA_readsize - Read header metadata of an RIA table file
 * @fname: Path to RIA table file.
 *
 * Populates globals RIA_NBpts, RIA_lambda_min [m], and RIA_lambda_max [m].
 *
 * Return: 1 on success, 0 on file open failure.
 */
int ATMOSPHEREMODEL_loadRIA_readsize(
    const char *fname);

/**
 * ATMOSPHEREMODEL_loadRIA - Load wavelength, refractive index, and absorption arrays
 * @fname:  Path to RIA table file.
 * @lptr:   Output wavelength array [m].
 * @RIptr:  Output refractive index array [dimensionless].
 * @absptr: Output absorption coefficient array [m^-1 at STP].
 *
 * Return: 0 on success, -1 on failure.
 */
int ATMOSPHEREMODEL_loadRIA(
    const char *fname,
    double     *lptr,
    double     *RIptr,
    double     *absptr);

#ifdef __cplusplus
}
#endif

#endif // ATMOSPHEREMODEL_H
