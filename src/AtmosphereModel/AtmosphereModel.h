// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    AtmosphereModel.h
 * @brief   Atmospheric model public interface
 */

#ifndef ATMOSPHEREMODEL_H
#define ATMOSPHEREMODEL_H

#ifdef __cplusplus
extern "C" {
#endif

// Shared site coordinate parameters
extern float SiteLat;
extern float SiteLong;
extern float SiteAlt;

/**
 * init_AtmosphereModel - Initialize module and allocate profile arrays
 *
 * Return: 0 on success.
 */
int init_AtmosphereModel(void);

/**
 * AirMixture_N - Compute refractive index minus 1 for gas mixture
 *
 * Computes (n - 1) for a 14-species atmospheric gas mixture using the
 * Lorentz-Lorenz relation and tabulated / Sellmeier dispersion values.
 * Absorption coefficient is accumulated into v_ABSCOEFF.
 */
double AirMixture_N(double lambda, double dens_N2, double dens_O2, double dens_Ar,
                    double dens_H2O, double dens_CO2, double dens_Ne, double dens_He,
                    double dens_CH4, double dens_Kr, double dens_H2, double dens_O3,
                    double dens_N, double dens_O, double dens_H);

/**
 * AtmosphereModel_stdAtmModel_N - Refractive index minus 1 at specified altitude
 * @alt: Altitude above sea level in meters.
 * @lambdaum: Wavelength in meters.
 * @mode: Verbosity / test mode flag (1 for verbose output, 0 for silent).
 *
 * Return: Refractivity (n - 1).
 */
float AtmosphereModel_stdAtmModel_N(float alt, float lambdaum, int mode);

/**
 * AtmosphereModel_H2O_Saturation - Water vapor saturation pressure
 * @T: Temperature in Kelvin.
 *
 * Return: Saturation vapor pressure in Pascals using IAPWS-95 formulation.
 */
double AtmosphereModel_H2O_Saturation(double T);

/**
 * AtmosphereModel_save_stdAtmModel - Save standard atmosphere profile to file
 * @fname: Output filename.
 *
 * Return: 0 on success.
 */
int AtmosphereModel_save_stdAtmModel(char *fname);

/**
 * AtmosphereModel_build_stdAtmModel - Generate profile from NRLMSISE-00
 * @fname: Output filename to save generated model.
 *
 * Return: 0 on success.
 */
int AtmosphereModel_build_stdAtmModel(char *fname);

/**
 * AtmosphereModel_load_stdAtmModel - Load atmosphere profile from file
 * @fname: Input filename.
 *
 * Return: 0 on success.
 */
int AtmosphereModel_load_stdAtmModel(char *fname);

/**
 * AtmosphereModel_Create_from_CONF - Build model from configuration file
 * @CONFFILE: Path to configuration file.
 * @slambda: Secondary test wavelength in meters.
 *
 * Return: 0 on success.
 */
int AtmosphereModel_Create_from_CONF(char *CONFFILE, float slambda);

/**
 * AtmosphereModel_RefractionPath - Compute ray trajectory through atmospheric layers
 * @lambda: Optical wavelength in meters.
 * @Zangle: Zenith angle at ground in radians.
 * @WritePath: Flag (1 to write refractpath.txt, 0 otherwise).
 *
 * Return: Atmospheric refraction deflection in arcseconds.
 */
double AtmosphereModel_RefractionPath(double lambda, double Zangle, int WritePath);

/**
 * ATMOSPHEREMODEL_loadRIA_readsize - Read header of RIA data file
 * @fname: RIA file path.
 *
 * Return: 1 on success, 0 on file open failure.
 */
int ATMOSPHEREMODEL_loadRIA_readsize(char *fname);

/**
 * ATMOSPHEREMODEL_loadRIA - Load wavelength, index, and absorption arrays
 * @fname: RIA file path.
 * @lptr: Pointer to wavelength array.
 * @RIptr: Pointer to refractive index array.
 * @absptr: Pointer to absorption coefficient array.
 *
 * Return: 0 on success.
 */
int ATMOSPHEREMODEL_loadRIA(char *fname, double *lptr, double *RIptr, double *absptr);

#ifdef __cplusplus
}
#endif

#endif // ATMOSPHEREMODEL_H
