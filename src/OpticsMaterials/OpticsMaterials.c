// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    OpticsMaterials.c
 * @brief   Optical material dispersion and refractive index lookups
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "OpticsMaterials.h"
#include "optmat_types.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846264338328
#endif

const struct MaterialIndex MatCode[] = {
    { "Mirror",   0 },
    {   "SiO2",   1 },
    {     "Si",   2 },
    {   "PMGI",   3 },
    {   "PMMA",   4 },
    {     "N2",   5 },
    {     "O2",   6 },
    {     "Ar",   7 },
    {     "He",   8 },
    {     "H2",   9 },
    {   "H2Og",  10 },
    {    "H2O",  10 },
    {    "CO2",  11 },
    {     "Ne",  12 },
    {      "O",  13 },
    {   "CaF2",  14 },
    { "Vacuum", 100 },
    {    "Air", 101 },
    {     NULL,   0 } /* end marker */
};

/**
 * init_OpticsMaterials - Initialize OpticsMaterials module
 *
 * Return: 0 on success.
 */
int init_OpticsMaterials(void)
{
    return 0;
}

/**
 * OPTICSMATERIALS_code - Lookup numeric material identifier by name
 * @name: Name of the optical material
 *
 * Return: Integer material code, or -1 if not recognized.
 */
int OPTICSMATERIALS_code(char *name)
{
    if (name == NULL)
    {
        return -1;
    }
    for (int i = 0; MatCode[i].name != NULL; i++)
    {
        if (strcmp(name, MatCode[i].name) == 0)
        {
            return MatCode[i].code;
        }
    }

    return -1;
}

/**
 * OPTICSMATERIALS_name - Lookup material name string by numeric code
 * @code: Integer material code
 *
 * Return: String name of the optical material, or NULL if not recognized.
 */
char *OPTICSMATERIALS_name(int code)
{
    for (int i = 0; MatCode[i].name != NULL; i++)
    {
        if (code == MatCode[i].code)
        {
            return MatCode[i].name;
        }
    }

    return NULL;
}

/**
 * OPTICSMATERIALS_n - Calculate refractive index of an optical material
 * @material: Material identifier code
 * @lambda: Wavelength in meters
 *
 * Orchestrates lookup across solids, gases, and resists based on material code.
 *
 * Return: Refractive index n at specified wavelength.
 */
double OPTICSMATERIALS_n(int material, double lambda)
{
    if (material < 0)
    {
        return 1.0;
    }

    double lambdaum = lambda * 1.0e6;
    double lambdanm = lambda * 1.0e9;
    double lambdaa  = lambdanm * 10.0;

    switch (material)
    {
    case 0:  // Mirror
    case 1:  // SiO2
    case 2:  // Si
    case 14: // CaF2
        return optmat_calc_refractive_index_solid(material, lambdaum);

    case 3: // PMGI resist
        return optmat_calc_refractive_index_pmgi(lambdanm);

    case 4: // PMMA resist
        return optmat_calc_refractive_index_pmma(lambdanm);

    default: // Gases & vacuum (100, 101, 5..13)
        return optmat_calc_refractive_index_gas(material, lambdaum, lambdaa);
    }
}

/**
 * OPTICSMATERIALS_pha_lambda - Compute phase offset across material thickness
 * @material: Material identifier code
 * @z: Physical thickness in meters
 * @lambda: Wavelength in meters
 *
 * Computes optical phase delay relative to vacuum across thickness z.
 *
 * Return: Phase delay in radians.
 */
double OPTICSMATERIALS_pha_lambda(int material, double z, double lambda)
{
    double n    = OPTICSMATERIALS_n(material, lambda);
    double nair = 1.0; // Vacuum reference

    return 2.0 * M_PI * (n - nair) * z / lambda;
}
