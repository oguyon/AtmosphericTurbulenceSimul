// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    optmat_gases.c
 * @brief   Dispersion formulas for atmospheric and noble gases
 */

#include <math.h>
#include <stdio.h>
#include "optmat_types.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846264338328
#endif

/**
 * optmat_calc_refractive_index_gas - Calculate refractive index for gaseous materials
 * @material: Material code (100: Vacuum, 101: Air, 5: N2, 6: O2, 7: Ar, 8: He, 9: H2,
 *            10: H2Og, 11: CO2, 12: Ne, 13: O)
 * @lambdaum: Wavelength in microns
 * @lambdaa: Wavelength in Angstroms
 *
 * Computes the refractive index of atmospheric gases at standard conditions.
 *
 * Return: Refractive index n, or 0.0 if not recognized.
 */
double optmat_calc_refractive_index_gas(int material, double lambdaum, double lambdaa)
{
    const double LoschmidtConstant = 2.6867805e25;
    double       l2                = lambdaum * lambdaum;
    double       LL;

    switch (material)
    {
    case 100: // Vacuum
        return 1.0;

    case 101: // Air (1 atm)
        return 1.0 + (5792105.0e-8 / (238.0185 - 1.0 / l2)) +
               (167917.0e-8 / (57.362 - 1.0 / l2));

    case 5: { // N2 (Guyon 2016)
        const double a1 = 9.30052e-5, b1 = 0.0951983;
        const double a2 = 5.74577e-05, b2 = 0.0264892;
        const double a3 = 0.000143257, b3 = 0.0687966;
        return 1.0 + (a1 * l2) / (l2 - b1 * b1) + (a2 * l2) / (l2 - b2 * b2) +
               (a3 * l2) / (l2 - b3 * b3);
    }

    case 6: { // O2 (Zhang et al. 2008)
        double n0 = 1.0 + 1e-8 * (15532.45 + 456402.97 / (50.0 - 1.0 / l2));
        LL        = (n0 * n0 - 1.0) / (n0 * n0 + 2.0);
        LL *= 293.15 / 273.15;
        return sqrt((2.0 * LL + 1.0) / (1.0 - LL));
    }

    case 7: // Ar (Peck and Fisher)
        return 1.0 + 1.0e-8 * (6322.05 + 2811641.7 / (144.0 - 1.0 / l2));

    case 8: // He (Weber 2003)
        return 1.0 + 0.01470091 / (423.98 - 1.0 / l2);

    case 9: // H2 (Leonard 1974)
        return 1.0 + 0.0175329 / (128.905 - 1.0 / l2);

    case 10: { // H2O vapor (Ciddor 1996 converted to STP)
        const double w0 = 295.235, w1 = 2.6422, w2 = 20.032380, w3 = 0.004028;
        double       l4 = l2 * l2;
        double       l6 = l4 * l2;
        double       n0 = 1.0 + 1.0e-8 * 1.022 * (w0 + w1 / l2 + w2 / l4 + w3 / l6);
        LL              = (n0 * n0 - 1.0) / (n0 * n0 + 2.0);
        LL *= (293.15 / 273.15) * (101325.0 / 1333.0);
        return sqrt((2.0 * LL + 1.0) / (1.0 - LL));
    }

    case 11: { // CO2 (Bideau-Mehu et al. 1973)
        return 1.0 + 0.06991 / (166.175 - 1.0 / l2) + 0.00144720 / (79.609 - 1.0 / l2) +
               0.0000642941 / (56.3064 - 1.0 / l2) + 0.0000521306 / (46.0196 - 1.0 / l2) +
               0.00000146847 / (0.0584738 - 1.0 / l2);
    }

    case 12: { // Ne (Weber 2003 / Bideau-Mehu et al. 1981)
        return 1.0 + 0.012055 * ((0.1063 * l2) / (184.661 * l2 - 1.0) +
                                 (1.8290 * l2) / (376.840 * l2 - 1.0));
    }

    case 13: { // O (atomic oxygen, Ivanova and Kologrivov 1970)
        double pol = 4.5e-30 * (1.0 + 6.2e5 / (lambdaa * lambdaa));
        double A   = (4.0 * M_PI / 3.0) * pol * LoschmidtConstant;
        return sqrt((2.0 * A + 1.0) / (1.0 - A));
    }

    default:
        return 0.0;
    }
}
