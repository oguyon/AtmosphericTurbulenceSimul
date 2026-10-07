// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_profile.c
 * @brief   Atmospheric turbulence layer profile loading, parsing, and normalization
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "atmturb_profile.h"

/**
 * atmturb_profile_free - Deallocate turbulence layer array
 * @prof: Pointer to profile container.
 */
void atmturb_profile_free(
    atmturb_profile_t *prof)
{
    if (prof == NULL)
    {
        return;
    }
    free(prof->layers);
    prof->layers = NULL;
    prof->nlayers = 0;
    prof->total_cn2_raw = 0.0;
}

/**
 * atmturb_profile_init_default - Populate standard 7-layer atmospheric profile
 * @prof: Output profile container to initialize.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_profile_init_default(
    atmturb_profile_t *prof)
{
    static const double def_alt[7] = {4215.0, 4230.0, 4349.0, 5007.0, 12000.0, 16200.0, 23701.0};
    static const double def_cn2[7] = {5.32,   1.47,   1.08,   2.11,   1.83,    1.48,    0.697};
    static const double def_spd[7] = {6.5,    6.55,   6.6,    6.7,   22.0,     9.5,     5.6};
    static const double def_dir[7] = {1.47,   1.57,   1.67,   1.77,   3.10,    3.20,    3.30};
    static const double def_l0[7]  = {25.0,   25.0,   25.0,   25.0,  100.0,   100.0,   100.0};

    atmturb_profile_free(prof);
    prof->layers = (atmturb_layer_t *) calloc(7, sizeof(atmturb_layer_t));
    if (prof->layers == NULL)
    {
        return -1;
    }

    prof->nlayers = 7;
    double raw_sum = 0.0;
    for (int k = 0; k < 7; k++)
    {
        raw_sum += def_cn2[k];
    }
    prof->total_cn2_raw = raw_sum;

    for (int k = 0; k < 7; k++)
    {
        prof->layers[k].alt_m = def_alt[k];
        prof->layers[k].cn2_frac = def_cn2[k] / raw_sum;
        prof->layers[k].speed_mps = def_spd[k];
        prof->layers[k].dir_rad = def_dir[k];
        prof->layers[k].L0_m = def_l0[k];
        prof->layers[k].l0_m = 0.0;
        prof->layers[k].sigma_wind_mps = 0.0;
        prof->layers[k].L_wind_m = 500.0;
    }

    return 0;
}

/**
 * atmturb_profile_is_data_line - Test whether a profile line carries layer data
 * @line: NUL-terminated text line.
 *
 * Return: 1 if the line is neither blank nor a '#' comment, 0 otherwise.
 */
static int atmturb_profile_is_data_line(
    const char *line)
{
    while (*line != '\0' && isspace((unsigned char) *line))
    {
        line++;
    }
    return (*line != '\0' && *line != '#') ? 1 : 0;
}

/**
 * atmturb_profile_parse_layer - Parse one profile text line into layer struct
 * @line: Text line from profile file.
 * @lineno: Line number for error logging.
 * @layer: Destination layer structure.
 *
 * Return: 0 on success, -1 on invalid format or negative Cn2.
 */
static int atmturb_profile_parse_layer(
    const char      *line,
    long             lineno,
    atmturb_layer_t *layer)
{
    double alt = 0.0, cn2 = 0.0, spd = 0.0, dir = 0.0;
    double L0 = -1.0, l0 = 0.0, sigma_w = 0.0, L_w = 500.0;

    int nf = sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf",
                    &alt, &cn2, &spd, &dir, &L0, &l0, &sigma_w, &L_w);
    if (nf < 4)
    {
        printf("ERROR: profile line %ld: expected at least 4 fields, got %d\n", lineno, nf);
        return -1;
    }
    if (cn2 < 0.0)
    {
        printf("ERROR: profile line %ld: negative Cn2 (%g)\n", lineno, cn2);
        return -1;
    }

    layer->alt_m = alt;
    layer->cn2_frac = cn2; // Temporary raw value, normalized later
    layer->speed_mps = spd;
    layer->dir_rad = dir;
    layer->L0_m = L0;
    layer->l0_m = l0;
    layer->sigma_wind_mps = sigma_w;
    layer->L_wind_m = L_w;

    return 0;
}

/**
 * atmturb_profile_load - Load and normalize turbulence profile from disk
 * @fname: Path to text profile file (or NULL to load built-in default).
 * @prof: Output profile container to initialize.
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_profile_load(
    const char        *fname,
    atmturb_profile_t *prof)
{
    atmturb_profile_free(prof);
    if (fname == NULL || fname[0] == '\0')
    {
        return atmturb_profile_init_default(prof);
    }

    FILE *fp = fopen(fname, "r");
    if (fp == NULL)
    {
        if (strcmp(fname, "turbul.prof") == 0)
        {
            printf("[milkatmturb] Notice: Profile \"%s\" not found, using default 7 layers.\n",
                   fname);
            return atmturb_profile_init_default(prof);
        }
        printf("ERROR: cannot open profile \"%s\"\n", fname);
        return -1;
    }

    char line[2000];
    int count = 0;
    while (fgets(line, sizeof(line), fp) != NULL)
    {
        count += atmturb_profile_is_data_line(line);
    }
    if (count == 0)
    {
        printf("ERROR: profile \"%s\" contains no data lines\n", fname);
        fclose(fp);
        return -1;
    }

    prof->layers = (atmturb_layer_t *) calloc((size_t) count, sizeof(atmturb_layer_t));
    if (prof->layers == NULL)
    {
        fclose(fp);
        return -1;
    }
    prof->nlayers = count;

    rewind(fp);
    int k = 0;
    long lineno = 0;
    double raw_sum = 0.0;
    while (fgets(line, sizeof(line), fp) != NULL && k < count)
    {
        lineno++;
        if (!atmturb_profile_is_data_line(line))
        {
            continue;
        }
        if (atmturb_profile_parse_layer(line, lineno, &prof->layers[k]) != 0)
        {
            atmturb_profile_free(prof);
            fclose(fp);
            return -1;
        }
        raw_sum += prof->layers[k].cn2_frac;
        k++;
    }
    fclose(fp);

    if (raw_sum <= 0.0)
    {
        printf("ERROR: total Cn2 sum is zero or non-positive (%g)\n", raw_sum);
        atmturb_profile_free(prof);
        return -1;
    }

    prof->total_cn2_raw = raw_sum;
    for (int i = 0; i < count; i++)
    {
        prof->layers[i].cn2_frac /= raw_sum;
    }

    printf("[milkatmturb] Profile \"%s\": loaded %d layers (raw sum Cn2 = %g)\n",
           fname, count, raw_sum);
    fflush(stdout);

    return 0;
}
