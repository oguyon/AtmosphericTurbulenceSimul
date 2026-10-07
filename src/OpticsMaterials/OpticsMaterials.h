// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    OpticsMaterials.h
 * @brief   Optical material dispersion and refractive index lookups
 */

#ifndef _OPTICSMATERIALS_H
#define _OPTICSMATERIALS_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * init_OpticsMaterials - Initialize OpticsMaterials module
 *
 * Return: 0 on success.
 */
int init_OpticsMaterials(void);

/**
 * OPTICSMATERIALS_code - Lookup numeric material identifier by name
 * @name: Name of the optical material
 *
 * Return: Integer material code, or -1 if not recognized.
 */
int OPTICSMATERIALS_code(
    const char *name);

/**
 * OPTICSMATERIALS_name - Lookup material name string by numeric code
 * @code: Integer material code
 *
 * Return: String name of the optical material, or NULL if not recognized.
 */
const char *OPTICSMATERIALS_name(
    int         code);

/**
 * OPTICSMATERIALS_n - Calculate refractive index of an optical material
 * @material: Material identifier code
 * @lambda: Wavelength in meters
 *
 * Return: Refractive index n, or 0.0 if not recognized.
 */
double OPTICSMATERIALS_n(
    int         material,
    double      lambda);

/**
 * OPTICSMATERIALS_pha_lambda - Phase offset as function of mask thickness and lambda
 * @material: Material identifier code
 * @z: Physical thickness in meters
 * @lambda: Wavelength in meters
 *
 * Return: Phase shift in radians.
 */
double OPTICSMATERIALS_pha_lambda(
    int         material,
    double      z,
    double      lambda);

#ifdef __cplusplus
}
#endif

#endif // _OPTICSMATERIALS_H
