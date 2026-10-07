// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmmod_gas_mixture.c
 * @brief   Atmospheric gas mixture refractive index calculation using Lorentz-Lorenz
 */

#include "AtmosphereModel.h"
#include "OpticsMaterials.h"
#include "atmmod_types.h"

/**
 * atmmod_eval_species_lorentz_lorenz - Compute Lorentz-Lorenz term for a gas species
 * @lambda: Wavelength in meters.
 * @dens: Number density in cm^-3.
 * @Z: Compressibility factor.
 * @initRIA: 1 if RIA table is loaded, 0 otherwise.
 * @nbpts: Number of points in RIA table.
 * @ria_lambda: RIA wavelength array.
 * @ria_ri: RIA refractive index array.
 * @ria_abs: RIA absorption coefficient array.
 * @lliprecomp: Pointer to cached table lookup index.
 * @mat_name: Material name string for OpticsMaterials fallback.
 * @abscoeff_accum: Pointer to accumulate absorption coefficient into.
 *
 * Return: Lorentz-Lorenz molar refractivity contribution.
 */
static double atmmod_eval_species_lorentz_lorenz(
    double lambda, double dens, double Z, int initRIA, long nbpts,
    const double *ria_lambda, const double *ria_ri, const double *ria_abs,
    int *lliprecomp, const char *mat_name, double *abscoeff_accum)
{
    double n = 1.0;
    double abscoeff = 0.0;

    if (initRIA == 1 && nbpts > 1 && ria_lambda != NULL)
    {
        long lli = *lliprecomp;
        if (lli < 0 || lli >= nbpts - 1)
        {
            lli = nbpts / 2;
        }

        long llistep = 100;
        int llidir = 1;
        while (llistep != 1)
        {
            llistep = (long)(0.3 * llistep);
            if (llistep == 0)
            {
                llistep = 1;
            }
            while ((ria_lambda[lli] * llidir < lambda * llidir) &&
                   (lli < nbpts - llistep) && (lli > llistep))
            {
                lli += llidir * llistep;
            }
            llidir = -llidir;
        }

        double denom = ria_lambda[lli + 1] - ria_lambda[lli];
        double alpha = (denom != 0.0) ? (lambda - ria_lambda[lli]) / denom : 0.0;
        *lliprecomp = (int)lli;
        n = (1.0 - alpha) * ria_ri[lli] + alpha * ria_ri[lli + 1];
        abscoeff = ria_abs[lli];
    }
    else
    {
        n = OPTICSMATERIALS_n(OPTICSMATERIALS_code((char *)mat_name), lambda);
        abscoeff = 0.0;
    }

    double LL = (n * n - 1.0) / (n * n + 2.0);
    double tmpc = dens / (ATMMOD_LOSCHMIDT / 1e6 / Z);
    LL *= tmpc;
    *abscoeff_accum += tmpc * abscoeff;

    return LL;
}

/**
 * AirMixture_N - Compute refractive index minus 1 for gas mixture
 * @lambda: Optical wavelength in meters.
 * @dens_N2: N2 density in cm^-3.
 * @dens_O2: O2 density in cm^-3.
 * @dens_Ar: Ar density in cm^-3.
 * @dens_H2O: H2O density in cm^-3.
 * @dens_CO2: CO2 density in cm^-3.
 * @dens_Ne: Ne density in cm^-3.
 * @dens_He: He density in cm^-3.
 * @dens_CH4: CH4 density in cm^-3.
 * @dens_Kr: Kr density in cm^-3.
 * @dens_H2: H2 density in cm^-3.
 * @dens_O3: O3 density in cm^-3.
 * @dens_N: Atomic N density in cm^-3.
 * @dens_O: Atomic O density in cm^-3.
 * @dens_H: Atomic H density in cm^-3.
 *
 * Return: Refractivity (n - 1).
 */
double AirMixture_N(double lambda, double dens_N2, double dens_O2, double dens_Ar,
                    double dens_H2O, double dens_CO2, double dens_Ne, double dens_He,
                    double dens_CH4, double dens_Kr, double dens_H2, double dens_O3,
                    double dens_N, double dens_O, double dens_H)
{
    double LL_total = 0.0;
    double abscoeff_accum = 0.0;

    // Compressibility factors at STP (1.013 bar, 15 deg C; ref: Air Liquide Encyclopedia)
    LL_total += atmmod_eval_species_lorentz_lorenz(
        lambda, dens_N2, 0.99971, initRIA_N2, RIA_N2_NBpts, RIA_N2_lambda,
        RIA_N2_ri, RIA_N2_abs, &lliprecompN2, "N2", &abscoeff_accum);
    LL_total += atmmod_eval_species_lorentz_lorenz(
        lambda, dens_O2, 0.99924, initRIA_O2, RIA_O2_NBpts, RIA_O2_lambda,
        RIA_O2_ri, RIA_O2_abs, &lliprecompO2, "O2", &abscoeff_accum);
    LL_total += atmmod_eval_species_lorentz_lorenz(
        lambda, dens_Ar, 0.99925, initRIA_Ar, RIA_Ar_NBpts, RIA_Ar_lambda,
        RIA_Ar_ri, RIA_Ar_abs, &lliprecompAr, "Ar", &abscoeff_accum);
    LL_total += atmmod_eval_species_lorentz_lorenz(
        lambda, dens_H2O, 1.0, initRIA_H2O, RIA_H2O_NBpts, RIA_H2O_lambda,
        RIA_H2O_ri, RIA_H2O_abs, &lliprecompH2O, "H2O", &abscoeff_accum);
    LL_total += atmmod_eval_species_lorentz_lorenz(
        lambda, dens_CO2, 0.99435, initRIA_CO2, RIA_CO2_NBpts, RIA_CO2_lambda,
        RIA_CO2_ri, RIA_CO2_abs, &lliprecompCO2, "CO2", &abscoeff_accum);
    LL_total += atmmod_eval_species_lorentz_lorenz(
        lambda, dens_Ne, 1.0005, initRIA_Ne, RIA_Ne_NBpts, RIA_Ne_lambda,
        RIA_Ne_ri, RIA_Ne_abs, &lliprecompNe, "Ne", &abscoeff_accum);
    LL_total += atmmod_eval_species_lorentz_lorenz(
        lambda, dens_He, 1.0005, initRIA_He, RIA_He_NBpts, RIA_He_lambda,
        RIA_He_ri, RIA_He_abs, &lliprecompHe, "He", &abscoeff_accum);
    LL_total += atmmod_eval_species_lorentz_lorenz(
        lambda, dens_CH4, 0.99802, initRIA_CH4, RIA_CH4_NBpts, RIA_CH4_lambda,
        RIA_CH4_ri, RIA_CH4_abs, &lliprecompCH4, "CH4", &abscoeff_accum);
    LL_total += atmmod_eval_species_lorentz_lorenz(
        lambda, dens_Kr, 0.99768, initRIA_Kr, RIA_Kr_NBpts, RIA_Kr_lambda,
        RIA_Kr_ri, RIA_Kr_abs, &lliprecompKr, "Kr", &abscoeff_accum);
    LL_total += atmmod_eval_species_lorentz_lorenz(
        lambda, dens_H2, 1.0006, initRIA_H2, RIA_H2_NBpts, RIA_H2_lambda,
        RIA_H2_ri, RIA_H2_abs, &lliprecompH2, "H2", &abscoeff_accum);
    LL_total += atmmod_eval_species_lorentz_lorenz(
        lambda, dens_O3, 1.0, initRIA_O3, RIA_O3_NBpts, RIA_O3_lambda,
        RIA_O3_ri, RIA_O3_abs, &lliprecompO3, "O3", &abscoeff_accum);
    LL_total += atmmod_eval_species_lorentz_lorenz(
        lambda, dens_N, 1.0, initRIA_N, RIA_N_NBpts, RIA_N_lambda,
        RIA_N_ri, RIA_N_abs, &lliprecompN, "N", &abscoeff_accum);
    LL_total += atmmod_eval_species_lorentz_lorenz(
        lambda, dens_O, 1.0, initRIA_O, RIA_O_NBpts, RIA_O_lambda,
        RIA_O_ri, RIA_O_abs, &lliprecompO, "O", &abscoeff_accum);
    LL_total += atmmod_eval_species_lorentz_lorenz(
        lambda, dens_H, 1.0, initRIA_H, RIA_H_NBpts, RIA_H_lambda,
        RIA_H_ri, RIA_H_abs, &lliprecompH, "H", &abscoeff_accum);

    v_ABSCOEFF = abscoeff_accum;
    double n = sqrt((2.0 * LL_total + 1.0) / (1.0 - LL_total));
    return n - 1.0;
}
