// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    optmat_types.h
 * @brief   Internal declarations and prototypes for OpticsMaterials module
 */

#ifndef OPTMAT_TYPES_H
#define OPTMAT_TYPES_H

#include <stddef.h>

struct MaterialIndex {
    char *name;
    int   code;
};

extern const struct MaterialIndex MatCode[];

/**
 * @brief Calculate refractive index for gaseous materials
 *
 * @param material Material code (Vacuum, Air, N2, O2, Ar, He, H2, H2Og, CO2, Ne, O)
 * @param lambdaum Wavelength in microns
 * @param lambdaa Wavelength in Angstroms
 * @return double Refractive index n
 */
double optmat_calc_refractive_index_gas(int material, double lambdaum, double lambdaa);

/**
 * @brief Calculate refractive index for solid optical materials
 *
 * @param material Material code (Mirror, SiO2, Si, CaF2)
 * @param lambdaum Wavelength in microns
 * @return double Refractive index n
 */
double optmat_calc_refractive_index_solid(int material, double lambdaum);

/**
 * @brief Calculate refractive index for PMGI resist from tabulated data
 *
 * @param lambdanm Wavelength in nanometers
 * @return double Refractive index n
 */
double optmat_calc_refractive_index_pmgi(double lambdanm);

/**
 * @brief Calculate refractive index for PMMA resist from tabulated data
 *
 * @param lambdanm Wavelength in nanometers
 * @return double Refractive index n
 */
double optmat_calc_refractive_index_pmma(double lambdanm);

#endif // OPTMAT_TYPES_H
