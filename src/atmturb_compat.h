// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_compat.h
 * @brief   Compatibility layer for milk framework-dev image APIs
 */

#ifndef _ATMTURB_COMPAT_H
#define _ATMTURB_COMPAT_H

#include <string.h>

#include "CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "linopt_imtools/compute_SVDpseudoInverse.h"

#ifndef dcimg
#define dcimg (milk_data.image)
#endif
#ifndef dcnimg
#define dcnimg (milk_data.NB_MAX_IMAGE)
#endif
#ifndef dcvar
#define dcvar (milk_data.variable)
#endif

// Compatibility aliases for legacy datatype fields and defines
#define atype datatype

#ifndef FLOAT
#define FLOAT _DATATYPE_FLOAT
#endif
#ifndef DOUBLE
#define DOUBLE _DATATYPE_DOUBLE
#endif
#ifndef COMPLEX_FLOAT
#define COMPLEX_FLOAT _DATATYPE_COMPLEX_FLOAT
#endif
#ifndef COMPLEX_DOUBLE
#define COMPLEX_DOUBLE _DATATYPE_COMPLEX_DOUBLE
#endif

#define save_db_fits(in, out) save_fits(in, out)

// image_ID wrapper
static inline imageID atmturb_image_ID(
    const char *name)
{
    imageID ID = image_ID(name, dcimg, dcnimg);
    if (ID < 0)
    {
        ID = read_sharedmem_image(name, dcimg, dcnimg);
    }
    return ID;
}

// delete_image_ID wrapper
static inline errno_t atmturb_delete_image_ID(
    const char *name)
{
    return delete_image_ID(name, DELETE_IMAGE_ERRMODE_IGNORE);
}

// arith_set_pixel compatibility helper
static inline imageID arith_set_pixel(
    const char *ID_name,
    double      value,
    long        x,
    long        y)
{
    imageID ID = image_ID(ID_name, dcimg, dcnimg);
    if (ID < 0)
    {
        return ID;
    }
    long naxes[2];
    naxes[0] = dcimg[ID].md[0].size[0];
    naxes[1] = dcimg[ID].md[0].size[1];
    if (x < 0 || x >= naxes[0] || y < 0 || y >= naxes[1])
    {
        return ID;
    }
    dcimg[ID].md[0].write = 1;
    if (dcimg[ID].md[0].datatype == _DATATYPE_FLOAT)
    {
        dcimg[ID].array.F[y * naxes[0] + x] = (float)value;
    }
    else if (dcimg[ID].md[0].datatype == _DATATYPE_DOUBLE)
    {
        dcimg[ID].array.D[y * naxes[0] + x] = value;
    }
    dcimg[ID].md[0].write = 0;
    dcimg[ID].md[0].cnt0++;
    COREMOD_MEMORY_image_set_sempost_byID(ID, -1);
    return ID;
}

// create_2Dimage_ID wrappers
static inline imageID atmturb_create_2Dimage_ID(
    const char *name,
    uint32_t    xsize,
    uint32_t    ysize)
{
    imageID ID = -1;
    create_2Dimage_ID(name, xsize, ysize, &ID);
    return ID;
}

static inline imageID atmturb_create_2Dimage_ID_double(
    const char *name,
    uint32_t    xsize,
    uint32_t    ysize)
{
    imageID ID = -1;
    create_2Dimage_ID_double(name, xsize, ysize, &ID);
    return ID;
}

static inline imageID atmturb_create_2DCimage_ID(
    const char *name,
    uint32_t    xsize,
    uint32_t    ysize)
{
    imageID ID = -1;
    create_2DCimage_ID(name, xsize, ysize, &ID);
    return ID;
}

static inline imageID atmturb_create_2DCimage_ID_double(
    const char *name,
    uint32_t    xsize,
    uint32_t    ysize)
{
    imageID ID = -1;
    create_2DCimage_ID_double(name, xsize, ysize, &ID);
    return ID;
}

// create_3Dimage_ID wrappers
static inline imageID atmturb_create_3Dimage_ID(
    const char *name,
    uint32_t    xsize,
    uint32_t    ysize,
    uint32_t    zsize)
{
    imageID ID = -1;
    uint32_t sz[3] = {xsize, ysize, zsize};
    create_image_ID(name, 3, sz, _DATATYPE_FLOAT, 1, 10, 0, &ID);
    return ID;
}

static inline imageID atmturb_create_3Dimage_ID_double(
    const char *name,
    uint32_t    xsize,
    uint32_t    ysize,
    uint32_t    zsize)
{
    imageID ID = -1;
    create_3Dimage_ID_double(name, xsize, ysize, zsize, &ID);
    return ID;
}

// create_image_ID wrapper (accepts both uint32_t* and long*)
static inline imageID atmturb_create_image_ID(
    const char *name,
    long        naxis,
    const void *size,
    uint8_t     datatype,
    int         shared,
    int         nbkw)
{
    uint32_t sz[8];
    const long *lsz = (const long *)size;
    for (int i = 0; i < naxis && i < 8; i++)
    {
        sz[i] = (uint32_t)lsz[i];
    }
    imageID ID = -1;
    create_image_ID(name, naxis, sz, datatype, shared, nbkw, 0, &ID);
    return ID;
}

// load_fits wrappers
static inline imageID atmturb_load_fits_2(
    const char *file,
    const char *name)
{
    imageID ID = -1;
    errno_t ret = load_fits(file, name, LOADFITS_ERRMODE_IGNORE, &ID);
    if (ID < 0 && ret == RETURN_SUCCESS)
    {
        ID = read_sharedmem_image(name, dcimg, dcnimg);
    }
    return ID;
}

static inline imageID atmturb_load_fits_3(
    const char *file,
    const char *name,
    int         errmode)
{
    imageID ID = -1;
    errno_t ret = load_fits(file, name, errmode, &ID);
    if (ID < 0 && ret == RETURN_SUCCESS)
    {
        ID = read_sharedmem_image(name, dcimg, dcnimg);
    }
    return ID;
}

// linopt_compute_SVDpseudoInverse wrapper
static inline errno_t atmturb_linopt_compute_SVDpseudoInverse(
    const char *r,
    const char *c,
    double      eps,
    long        maxnb,
    const char *vt)
{
    imageID outID = -1;
    return linopt_compute_SVDpseudoInverse(r, c, eps, maxnb, vt, &outID);
}

// arith_image_zero helper
static inline errno_t arith_image_zero(
    const char *name)
{
    imageID ID = image_ID(name, dcimg, dcnimg);
    if (ID != -1)
    {
        memset(dcimg[ID].array.raw, 0,
               ImageStreamIO_typesize(dcimg[ID].md[0].datatype) * dcimg[ID].md[0].nelement);
        return 0;
    }
    return -1;
}

// Overload macros using __VA_ARGS__
#define GET_MACRO_LOADFITS(_1, _2, _3, _4, NAME, ...) NAME
#define load_fits(...) \
    GET_MACRO_LOADFITS(__VA_ARGS__, load_fits_4, atmturb_load_fits_3, \
                       atmturb_load_fits_2)(__VA_ARGS__)
#define load_fits_4(f, n, e, id) load_fits(f, n, e, id)

#define GET_MACRO_IMAGE_ID(_1, _2, _3, NAME, ...) NAME
#define image_ID(...) \
    GET_MACRO_IMAGE_ID(__VA_ARGS__, image_ID_3, image_ID_2, atmturb_image_ID)(__VA_ARGS__)
#define image_ID_3(n, a, nb) image_ID(n, a, nb)

#define GET_MACRO_DELETE_IMG(_1, _2, NAME, ...) NAME
#define delete_image_ID(...) \
    GET_MACRO_DELETE_IMG(__VA_ARGS__, delete_image_ID_2, atmturb_delete_image_ID)(__VA_ARGS__)
#define delete_image_ID_2(n, e) delete_image_ID(n, e)

#define GET_MACRO_CREATE_2D(_1, _2, _3, _4, NAME, ...) NAME
#define create_2Dimage_ID(...) \
    GET_MACRO_CREATE_2D(__VA_ARGS__, create_2Dimage_ID_4, atmturb_create_2Dimage_ID)(__VA_ARGS__)
#define create_2Dimage_ID_4(n, x, y, id) create_2Dimage_ID(n, x, y, id)

#define GET_MACRO_CREATE_2D_DBL(_1, _2, _3, _4, NAME, ...) NAME
#define create_2Dimage_ID_double(...) \
    GET_MACRO_CREATE_2D_DBL(__VA_ARGS__, create_2Dimage_ID_double_4, \
                            atmturb_create_2Dimage_ID_double)(__VA_ARGS__)
#define create_2Dimage_ID_double_4(n, x, y, id) create_2Dimage_ID_double(n, x, y, id)

#define GET_MACRO_CREATE_2DC(_1, _2, _3, _4, NAME, ...) NAME
#define create_2DCimage_ID(...) \
    GET_MACRO_CREATE_2DC(__VA_ARGS__, create_2DCimage_ID_4, atmturb_create_2DCimage_ID)(__VA_ARGS__)
#define create_2DCimage_ID_4(n, x, y, id) create_2DCimage_ID(n, x, y, id)

#define GET_MACRO_CREATE_2DC_DBL(_1, _2, _3, _4, NAME, ...) NAME
#define create_2DCimage_ID_double(...) \
    GET_MACRO_CREATE_2DC_DBL(__VA_ARGS__, create_2DCimage_ID_double_4, \
                             atmturb_create_2DCimage_ID_double)(__VA_ARGS__)
#define create_2DCimage_ID_double_4(n, x, y, id) create_2DCimage_ID_double(n, x, y, id)

#define GET_MACRO_CREATE_3D(_1, _2, _3, _4, _5, NAME, ...) NAME
#define create_3Dimage_ID(...) \
    GET_MACRO_CREATE_3D(__VA_ARGS__, create_3Dimage_ID_5, atmturb_create_3Dimage_ID)(__VA_ARGS__)
#define create_3Dimage_ID_5(n, x, y, z, id) create_3Dimage_ID(n, x, y, z, id)

#define GET_MACRO_CREATE_3D_DBL(_1, _2, _3, _4, _5, NAME, ...) NAME
#define create_3Dimage_ID_double(...) \
    GET_MACRO_CREATE_3D_DBL(__VA_ARGS__, create_3Dimage_ID_double_5, \
                            atmturb_create_3Dimage_ID_double)(__VA_ARGS__)
#define create_3Dimage_ID_double_5(n, x, y, z, id) create_3Dimage_ID_double(n, x, y, z, id)

#define GET_MACRO_CREATE_IMG(_1, _2, _3, _4, _5, _6, _7, _8, NAME, ...) NAME
#define create_image_ID(...) \
    GET_MACRO_CREATE_IMG(__VA_ARGS__, create_image_ID_8, create_image_ID_7, \
                         atmturb_create_image_ID)(__VA_ARGS__)
#define create_image_ID_8(n, ax, s, dt, sh, nkw, cb, id) \
    create_image_ID(n, ax, s, dt, sh, nkw, cb, id)

#define GET_MACRO_SVD(_1, _2, _3, _4, _5, _6, NAME, ...) NAME
#define linopt_compute_SVDpseudoInverse(...) \
    GET_MACRO_SVD(__VA_ARGS__, linopt_compute_SVDpseudoInverse_6, \
                  atmturb_linopt_compute_SVDpseudoInverse)(__VA_ARGS__)
#define linopt_compute_SVDpseudoInverse_6(r, c, e, m, v, o) \
    linopt_compute_SVDpseudoInverse(r, c, e, m, v, o)

#endif // _ATMTURB_COMPAT_H
