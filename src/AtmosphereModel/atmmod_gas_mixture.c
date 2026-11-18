// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmmod_gas_mixture.c
 * @brief   Atmospheric gas mixture refractive index calculation using Lorentz-Lorenz
 */

#include <math.h>
#include <stdatomic.h>
#include <stdint.h>

#include "AtmosphereModel.h"
#include "OpticsMaterials/OpticsMaterials.h"
#include "atmmod_types.h"

#define ATMMOD_NUM_SPECIES 14

/**
 * struct atmmod_species_desc_t - Static descriptor for an atmospheric gas species
 * @init_ria:   Pointer to species RIA table loaded flag.
 * @nbpts:      Pointer to RIA sample count.
 * @ria_lambda: Pointer to RIA wavelength array pointer.
 * @ria_ri:     Pointer to RIA refractive index array pointer.
 * @ria_abs:    Pointer to RIA absorption coefficient array pointer.
 * @lliprecomp: Pointer to legacy cached lookup index.
 * @inv_c_stp:  Reciprocal STP Loschmidt scaling factor Z / (N_L * 1e-6).
 * @mat_name:   Material identifier string for OpticsMaterials fallback.
 */
typedef struct
{
    const int     *init_ria;
    const long    *nbpts;
    double *const *ria_lambda;
    double *const *ria_ri;
    double *const *ria_abs;
    int           *lliprecomp;
    double         inv_c_stp;
    const char    *mat_name;
} atmmod_species_desc_t;

// Compressibility factors at STP (1.013 bar, 15 deg C; ref: Air Liquide Encyclopedia)
static const atmmod_species_desc_t g_species_desc[ATMMOD_NUM_SPECIES] = {
    {&initRIA_N2,  &RIA_N2_NBpts,  &RIA_N2_lambda,  &RIA_N2_ri,  &RIA_N2_abs,
     &lliprecompN2,  0.99971 / (ATMMOD_LOSCHMIDT / 1e6), "N2"},
    {&initRIA_O2,  &RIA_O2_NBpts,  &RIA_O2_lambda,  &RIA_O2_ri,  &RIA_O2_abs,
     &lliprecompO2,  0.99924 / (ATMMOD_LOSCHMIDT / 1e6), "O2"},
    {&initRIA_Ar,  &RIA_Ar_NBpts,  &RIA_Ar_lambda,  &RIA_Ar_ri,  &RIA_Ar_abs,
     &lliprecompAr,  0.99925 / (ATMMOD_LOSCHMIDT / 1e6), "Ar"},
    {&initRIA_H2O, &RIA_H2O_NBpts, &RIA_H2O_lambda, &RIA_H2O_ri, &RIA_H2O_abs,
     &lliprecompH2O, 1.00000 / (ATMMOD_LOSCHMIDT / 1e6), "H2O"},
    {&initRIA_CO2, &RIA_CO2_NBpts, &RIA_CO2_lambda, &RIA_CO2_ri, &RIA_CO2_abs,
     &lliprecompCO2, 0.99435 / (ATMMOD_LOSCHMIDT / 1e6), "CO2"},
    {&initRIA_Ne,  &RIA_Ne_NBpts,  &RIA_Ne_lambda,  &RIA_Ne_ri,  &RIA_Ne_abs,
     &lliprecompNe,  1.00050 / (ATMMOD_LOSCHMIDT / 1e6), "Ne"},
    {&initRIA_He,  &RIA_He_NBpts,  &RIA_He_lambda,  &RIA_He_ri,  &RIA_He_abs,
     &lliprecompHe,  1.00050 / (ATMMOD_LOSCHMIDT / 1e6), "He"},
    {&initRIA_CH4, &RIA_CH4_NBpts, &RIA_CH4_lambda, &RIA_CH4_ri, &RIA_CH4_abs,
     &lliprecompCH4, 0.99802 / (ATMMOD_LOSCHMIDT / 1e6), "CH4"},
    {&initRIA_Kr,  &RIA_Kr_NBpts,  &RIA_Kr_lambda,  &RIA_Kr_ri,  &RIA_Kr_abs,
     &lliprecompKr,  0.99768 / (ATMMOD_LOSCHMIDT / 1e6), "Kr"},
    {&initRIA_H2,  &RIA_H2_NBpts,  &RIA_H2_lambda,  &RIA_H2_ri,  &RIA_H2_abs,
     &lliprecompH2,  1.00060 / (ATMMOD_LOSCHMIDT / 1e6), "H2"},
    {&initRIA_O3,  &RIA_O3_NBpts,  &RIA_O3_lambda,  &RIA_O3_ri,  &RIA_O3_abs,
     &lliprecompO3,  1.00000 / (ATMMOD_LOSCHMIDT / 1e6), "O3"},
    {&initRIA_N,   &RIA_N_NBpts,   &RIA_N_lambda,   &RIA_N_ri,   &RIA_N_abs,
     &lliprecompN,   1.00000 / (ATMMOD_LOSCHMIDT / 1e6), "N"},
    {&initRIA_O,   &RIA_O_NBpts,   &RIA_O_lambda,   &RIA_O_ri,   &RIA_O_abs,
     &lliprecompO,   1.00000 / (ATMMOD_LOSCHMIDT / 1e6), "O"},
    {&initRIA_H,   &RIA_H_NBpts,   &RIA_H_lambda,   &RIA_H_ri,   &RIA_H_abs,
     &lliprecompH,   1.00000 / (ATMMOD_LOSCHMIDT / 1e6), "H"}
};

static _Atomic uint64_t g_mixture_cache_gen = 1;

static _Thread_local double       s_cached_lambda   = -1.0;
static _Thread_local uint64_t     s_cached_gen      = 0;
static _Thread_local unsigned int s_cached_ria_mask = 0;
static _Thread_local int          s_mat_codes_ready = 0;
static _Thread_local int          s_mat_codes[ATMMOD_NUM_SPECIES];
static _Thread_local int          s_lli_idx[ATMMOD_NUM_SPECIES];
static _Thread_local double       s_ll_unit[ATMMOD_NUM_SPECIES];
static _Thread_local double       s_abs_unit[ATMMOD_NUM_SPECIES];

/**
 * atmmod_invalidate_mixture_cache - Invalidate cached per-wavelength gas mixture coefficients
 */
void atmmod_invalidate_mixture_cache(void)
{
    atomic_fetch_add_explicit(&g_mixture_cache_gen, 1, memory_order_relaxed);
}

/**
 * atmmod_eval_species_unit_coeffs - Evaluate unit Lorentz-Lorenz and absorption for one species
 * @k:      Species index in g_species_desc [0..13].
 * @lambda: Optical wavelength in meters.
 */
static void atmmod_eval_species_unit_coeffs(
    int    k,
    double lambda)
{
    const atmmod_species_desc_t *d       = &g_species_desc[k];
    double                       n       = 1.0;
    double                       absc    = 0.0;
    long                         nbpts   = *d->nbpts;
    const double                *ria_lam = *d->ria_lambda;

    if (*d->init_ria == 1 && nbpts > 1 && ria_lam != NULL)
    {
        const double *ria_ri  = *d->ria_ri;
        const double *ria_abs = *d->ria_abs;
        long          lli     = s_lli_idx[k];
        if (lli < 0 || lli >= nbpts - 1)
        {
            lli = (*d->lliprecomp >= 0 && *d->lliprecomp < nbpts - 1)
                      ? (long) (*d->lliprecomp)
                      : (nbpts / 2);
        }

        long llistep = 100;
        int  llidir  = 1;
        while (llistep != 1)
        {
            llistep = (long) (0.3 * (double) llistep);
            if (llistep == 0)
            {
                llistep = 1;
            }
            while ((ria_lam[lli] * (double) llidir < lambda * (double) llidir) &&
                   (lli < nbpts - llistep) && (lli > llistep))
            {
                lli += (long) llidir * llistep;
            }
            llidir = -llidir;
        }

        double denom     = ria_lam[lli + 1] - ria_lam[lli];
        double alpha     = (denom != 0.0) ? (lambda - ria_lam[lli]) / denom : 0.0;
        s_lli_idx[k]     = (int) lli;
        *d->lliprecomp   = (int) lli;
        n                = (1.0 - alpha) * ria_ri[lli] + alpha * ria_ri[lli + 1];
        absc             = ria_abs[lli];
    }
    else
    {
        n    = OPTICSMATERIALS_n(s_mat_codes[k], lambda);
        absc = 0.0;
    }

    double n2     = n * n;
    s_ll_unit[k]  = ((n2 - 1.0) / (n2 + 2.0)) * d->inv_c_stp;
    s_abs_unit[k] = absc * d->inv_c_stp;
}

/**
 * atmmod_refresh_wavelength_cache - Update thread-local species unit coefficients on cache miss
 * @lambda: Optical wavelength in meters.
 */
static void atmmod_refresh_wavelength_cache(
    double lambda)
{
    uint64_t     gen      = atomic_load_explicit(&g_mixture_cache_gen, memory_order_relaxed);
    unsigned int ria_mask = 0;
    for (int k = 0; k < ATMMOD_NUM_SPECIES; k++)
    {
        ria_mask |= (unsigned int) ((*g_species_desc[k].init_ria != 0) ? (1u << k) : 0u);
    }

    if (lambda == s_cached_lambda && gen == s_cached_gen && ria_mask == s_cached_ria_mask)
    {
        return;
    }

    if (!s_mat_codes_ready)
    {
        for (int k = 0; k < ATMMOD_NUM_SPECIES; k++)
        {
            s_mat_codes[k] = OPTICSMATERIALS_code((char *) g_species_desc[k].mat_name);
            s_lli_idx[k]   = -1;
        }
        s_mat_codes_ready = 1;
    }
    else if (gen != s_cached_gen)
    {
        for (int k = 0; k < ATMMOD_NUM_SPECIES; k++)
        {
            s_lli_idx[k] = -1;
        }
    }

    for (int k = 0; k < ATMMOD_NUM_SPECIES; k++)
    {
        atmmod_eval_species_unit_coeffs(k, lambda);
    }

    s_cached_lambda   = lambda;
    s_cached_gen      = gen;
    s_cached_ria_mask = ria_mask;
}

/**
 * AirMixture_N - Compute refractive index minus 1 for gas mixture
 * @lambda:   Optical wavelength in meters.
 * @dens_N2:  N2 density in cm^-3.
 * @dens_O2:  O2 density in cm^-3.
 * @dens_Ar:  Ar density in cm^-3.
 * @dens_H2O: H2O density in cm^-3.
 * @dens_CO2: CO2 density in cm^-3.
 * @dens_Ne:  Ne density in cm^-3.
 * @dens_He:  He density in cm^-3.
 * @dens_CH4: CH4 density in cm^-3.
 * @dens_Kr:  Kr density in cm^-3.
 * @dens_H2:  H2 density in cm^-3.
 * @dens_O3:  O3 density in cm^-3.
 * @dens_N:   Atomic N density in cm^-3.
 * @dens_O:   Atomic O density in cm^-3.
 * @dens_H:   Atomic H density in cm^-3.
 *
 * Return: Refractivity (n - 1).
 */
/* size-exempt: public API signature in AtmosphereModel.h requires 15 species arguments */
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
    double dens_H)
{
    atmmod_refresh_wavelength_cache(lambda);

    const double dens[ATMMOD_NUM_SPECIES] = {
        dens_N2,  dens_O2, dens_Ar, dens_H2O, dens_CO2, dens_Ne, dens_He,
        dens_CH4, dens_Kr, dens_H2, dens_O3,  dens_N,   dens_O,  dens_H};

    double ll_total       = 0.0;
    double abscoeff_accum = 0.0;
    for (int k = 0; k < ATMMOD_NUM_SPECIES; k++)
    {
        ll_total       += dens[k] * s_ll_unit[k];
        abscoeff_accum += dens[k] * s_abs_unit[k];
    }

    v_ABSCOEFF = abscoeff_accum;
    return sqrt((2.0 * ll_total + 1.0) / (1.0 - ll_total)) - 1.0;
}
