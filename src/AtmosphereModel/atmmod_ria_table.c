// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmmod_ria_table.c
 * @brief   Atmospheric gas refractive index and absorption table loader
 */

#include "AtmosphereModel.h"
#include "atmmod_types.h"

// Definition of shared table parameters
long RIA_NBpts = 0;
double RIA_lambda_min = 0.0;
double RIA_lambda_max = 0.0;

int initRIA_N2 = 0;
long RIA_N2_NBpts = 0;
double *RIA_N2_lambda = NULL;
double *RIA_N2_ri = NULL;
double *RIA_N2_abs = NULL;

int initRIA_O2 = 0;
long RIA_O2_NBpts = 0;
double *RIA_O2_lambda = NULL;
double *RIA_O2_ri = NULL;
double *RIA_O2_abs = NULL;

int initRIA_Ar = 0;
long RIA_Ar_NBpts = 0;
double *RIA_Ar_lambda = NULL;
double *RIA_Ar_ri = NULL;
double *RIA_Ar_abs = NULL;

int initRIA_H2O = 0;
long RIA_H2O_NBpts = 0;
double *RIA_H2O_lambda = NULL;
double *RIA_H2O_ri = NULL;
double *RIA_H2O_abs = NULL;

int initRIA_CO2 = 0;
long RIA_CO2_NBpts = 0;
double *RIA_CO2_lambda = NULL;
double *RIA_CO2_ri = NULL;
double *RIA_CO2_abs = NULL;

int initRIA_Ne = 0;
long RIA_Ne_NBpts = 0;
double *RIA_Ne_lambda = NULL;
double *RIA_Ne_ri = NULL;
double *RIA_Ne_abs = NULL;

int initRIA_He = 0;
long RIA_He_NBpts = 0;
double *RIA_He_lambda = NULL;
double *RIA_He_ri = NULL;
double *RIA_He_abs = NULL;

int initRIA_CH4 = 0;
long RIA_CH4_NBpts = 0;
double *RIA_CH4_lambda = NULL;
double *RIA_CH4_ri = NULL;
double *RIA_CH4_abs = NULL;

int initRIA_Kr = 0;
long RIA_Kr_NBpts = 0;
double *RIA_Kr_lambda = NULL;
double *RIA_Kr_ri = NULL;
double *RIA_Kr_abs = NULL;

int initRIA_H2 = 0;
long RIA_H2_NBpts = 0;
double *RIA_H2_lambda = NULL;
double *RIA_H2_ri = NULL;
double *RIA_H2_abs = NULL;

int initRIA_O3 = 0;
long RIA_O3_NBpts = 0;
double *RIA_O3_lambda = NULL;
double *RIA_O3_ri = NULL;
double *RIA_O3_abs = NULL;

int initRIA_N = 0;
long RIA_N_NBpts = 0;
double *RIA_N_lambda = NULL;
double *RIA_N_ri = NULL;
double *RIA_N_abs = NULL;

int initRIA_O = 0;
long RIA_O_NBpts = 0;
double *RIA_O_lambda = NULL;
double *RIA_O_ri = NULL;
double *RIA_O_abs = NULL;

int initRIA_H = 0;
long RIA_H_NBpts = 0;
double *RIA_H_lambda = NULL;
double *RIA_H_ri = NULL;
double *RIA_H_abs = NULL;

int lliprecompN2 = -1;
int lliprecompO2 = -1;
int lliprecompAr = -1;
int lliprecompH2O = -1;
int lliprecompCO2 = -1;
int lliprecompNe = -1;
int lliprecompHe = -1;
int lliprecompCH4 = -1;
int lliprecompKr = -1;
int lliprecompH2 = -1;
int lliprecompO3 = -1;
int lliprecompO = -1;
int lliprecompN = -1;
int lliprecompH = -1;

/**
 * atmmod_reset_precomp_indices - Reset table lookup indices
 */
void atmmod_reset_precomp_indices(void)
{
    lliprecompN2 = -1;
    lliprecompO2 = -1;
    lliprecompAr = -1;
    lliprecompH2O = -1;
    lliprecompCO2 = -1;
    lliprecompNe = -1;
    lliprecompHe = -1;
    lliprecompCH4 = -1;
    lliprecompKr = -1;
    lliprecompH2 = -1;
    lliprecompO3 = -1;
    lliprecompO = -1;
    lliprecompN = -1;
    lliprecompH = -1;
}

/**
 * ATMOSPHEREMODEL_loadRIA_readsize - Read header metadata of RIA table file
 * @fname: Path to RIA table file.
 *
 * Return: 1 on success, 0 on failure.
 */
int ATMOSPHEREMODEL_loadRIA_readsize(
    const char *fname)
{
    FILE *fp = fopen(fname, "r");
    if (fp == NULL)
    {
        printf("ERROR: cannot open file \"%s\"\n", fname);
        return 0;
    }

    long nbpt = 0;
    double lmin = 0.0;
    double lmax = 0.0;
    if (fscanf(fp, "# %ld %lf %lf\n", &nbpt, &lmin, &lmax) == 3)
    {
        printf("%ld points from %g m to %g m\n", nbpt, lmin, lmax);
        RIA_NBpts = nbpt;
        RIA_lambda_min = lmin;
        RIA_lambda_max = lmax;
    }
    fclose(fp);
    return 1;
}

/**
 * ATMOSPHEREMODEL_loadRIA - Load RIA table values into allocated buffers
 * @fname: Path to RIA table file.
 * @lptr: Output wavelength array buffer.
 * @RIptr: Output refractive index array buffer.
 * @absptr: Output absorption coefficient array buffer.
 *
 * Return: 0 on success, -1 on failure.
 */
int ATMOSPHEREMODEL_loadRIA(
    const char *fname,
    double     *lptr,
    double     *RIptr,
    double     *absptr)
{
    FILE *fp = fopen(fname, "r");
    if (fp == NULL)
    {
        printf("ERROR: cannot open file \"%s\"\n", fname);
        return -1;
    }

    long nbpt = 0;
    double lmin = 0.0;
    double lmax = 0.0;
    if (fscanf(fp, "# %ld %lf %lf\n", &nbpt, &lmin, &lmax) != 3)
    {
        fclose(fp);
        return -1;
    }
    printf("%ld points from %g m to %g m\n", nbpt, lmin, lmax);

    for (long i = 0; i < nbpt; i++)
    {
        double v0 = 0.0;
        double v1 = 0.0;
        double v2 = 0.0;
        if (fscanf(fp, "%lf %lf %lf\n", &v0, &v1, &v2) == 3)
        {
            lptr[i] = v0;
            RIptr[i] = v1;
            absptr[i] = v2;
        }
    }
    fclose(fp);
    return 0;
}

/**
 * atmmod_load_species_ria - Load RIA table for a single chemical species
 * @name: Species name (e.g. "N2", "O2").
 * @path: File path to table.
 * @init_flag: Pointer to species initialization flag.
 * @nbpts_out: Pointer to store number of points.
 * @lambda_out: Pointer to store allocated wavelength array.
 * @ri_out: Pointer to store allocated refractive index array.
 * @abs_out: Pointer to store allocated absorption array.
 */
void atmmod_load_species_ria(const char *name, const char *path, int *init_flag,
                             long *nbpts_out, double **lambda_out, double **ri_out,
                             double **abs_out)
{
    printf("Reading Refractive Index and Abs for %s\n", name);
    char pathbuf[256];
    snprintf(pathbuf, sizeof(pathbuf), "%s", path);
    if (ATMOSPHEREMODEL_loadRIA_readsize(pathbuf) == 1)
    {
        long nb = RIA_NBpts;
        *nbpts_out = nb;
        *lambda_out = (double *)malloc(sizeof(double) * nb);
        *ri_out = (double *)malloc(sizeof(double) * nb);
        *abs_out = (double *)malloc(sizeof(double) * nb);
        ATMOSPHEREMODEL_loadRIA(pathbuf, *lambda_out, *ri_out, *abs_out);
        *init_flag = 1;
    }
}
