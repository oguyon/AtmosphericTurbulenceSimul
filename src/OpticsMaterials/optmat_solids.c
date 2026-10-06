// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    optmat_solids.c
 * @brief   Sellmeier dispersion equations for solid optical materials
 */

#include <math.h>
#include <stdio.h>
#include "optmat_types.h"

/**
 * optmat_calc_refractive_index_solid - Calculate refractive index for solid materials
 * @material: Material identifier (0: Mirror, 1: SiO2, 2: Si, 14: CaF2)
 * @lambdaum: Wavelength in microns
 *
 * Evaluates Sellmeier dispersion relations for solid glasses and semiconductors.
 *
 * Return: Refractive index n, or 0.0 if material is not recognized.
 */
double optmat_calc_refractive_index_solid(int material, double lambdaum)
{
    const double SiO2_n0 = 1.28604141;
    const double SiO2_B1 = 1.07044083;
    const double SiO2_C1 = 1.00585997e-2;
    const double SiO2_B2 = 1.10202242;
    const double SiO2_C2 = 100.0;

    const double Si_B1 = 10.6684293;
    const double Si_C1 = 0.301516485;
    const double Si_B2 = 0.003043475;
    const double Si_C2 = 1.13475115;
    const double Si_B3 = 1.54133408;
    const double Si_C3 = 1104.0;

    const double CaF2_B1 = 0.69913;
    const double CaF2_C1 = 0.09374;
    const double CaF2_B2 = 0.11994;
    const double CaF2_C2 = 21.18;
    const double CaF2_B3 = 4.35181;
    const double CaF2_C3 = 38.46;

    double l2 = lambdaum * lambdaum;

    switch (material)
    {
    case 0: // Mirror
        return 3.0;

    case 1: // SiO2
        return sqrt(SiO2_n0 + (SiO2_B1 * l2) / (l2 - SiO2_C1) + (SiO2_B2 * l2) / (l2 - SiO2_C2));

    case 2: // Si
        return sqrt(1.0 + (Si_B1 * l2) / (l2 - Si_C1 * Si_C1) +
                    (Si_B2 * l2) / (l2 - Si_C2 * Si_C2) +
                    (Si_B3 * l2) / (l2 - Si_C3 * Si_C3));

    case 14: // CaF2
        return sqrt(1.33973 + (CaF2_B1 * l2) / (l2 - CaF2_C1 * CaF2_C1) +
                    (CaF2_B2 * l2) / (l2 - CaF2_C2 * CaF2_C2) +
                    (CaF2_B3 * l2) / (l2 - CaF2_C3 * CaF2_C3));

    default:
        return 0.0;
    }
}
