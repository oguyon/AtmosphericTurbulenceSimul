// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_compat.h
 * @brief   Internal convenience shims over milk_data and core image APIs
 *
 * Bridges concise return-by-imageID helper signatures used across milkatmturb
 * onto the explicit out-pointer and milk_data.image APIs of milk
 * (feat/bash-engine). Also provides variadic argument-count dispatch macros so
 * both legacy short-form calls and full-argument core calls compile cleanly.
 */

#ifndef _ATMTURB_COMPAT_H
#define _ATMTURB_COMPAT_H

#include <string.h>

#include "atmturb_c23_compat.h"
#include "milk_runnable.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "linopt_imtools/compute_SVDpseudoInverse.h"

/* Global workspace table shortcuts into milk_data */
#ifndef dcimg
#define dcimg (milk_data.image)
#endif
#ifndef dcnimg
#define dcnimg (milk_data.NB_MAX_IMAGE)
#endif
#ifndef dcvar
#define dcvar (milk_data.variable)
#endif

/* Compatibility aliases for ImageStreamIO datatype fields and constants */
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

/**
 * atmturb_image_ID - Resolve an image by name in memory or shared memory
 * @name: Image or shared-memory stream name.
 *
 * Searches the process-local milk_data.image table first and falls back to
 * attaching an existing ImageStreamIO shared-memory stream of the same name.
 *
 * Return: Valid imageID (>= 0) on success, or -1 if not found.
 */
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

/**
 * atmturb_delete_image_ID - Delete a named image ignoring missing-image errors
 * @name: Image name in milk_data.image.
 *
 * Return: RETURN_SUCCESS or error code from delete_image_ID().
 */
static inline errno_t atmturb_delete_image_ID(
    const char *name)
{
    return delete_image_ID(name, DELETE_IMAGE_ERRMODE_IGNORE);
}

/**
 * arith_set_pixel - Write a single 2D pixel and post stream semaphores
 * @ID_name: Target 2D float or double image name.
 * @value:   Scalar pixel value to store.
 * @x:       Zero-based column index [pixels].
 * @y:       Zero-based row index [pixels].
 *
 * Return: Resolved imageID (>= 0), or -1 if the image does not exist.
 */
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

/**
 * atmturb_create_2Dimage_ID - Create a 2D single-precision float image
 * @name:  Image name to allocate in milk_data.image.
 * @xsize: Width along axis 0 [pixels].
 * @ysize: Height along axis 1 [pixels].
 *
 * Return: Newly allocated imageID, or -1 on failure.
 */
static inline imageID atmturb_create_2Dimage_ID(
    const char *name,
    uint32_t    xsize,
    uint32_t    ysize)
{
    imageID ID = -1;
    create_2Dimage_ID(name, xsize, ysize, &ID);
    return ID;
}

/**
 * atmturb_create_2Dimage_ID_double - Create a 2D double-precision float image
 * @name:  Image name to allocate in milk_data.image.
 * @xsize: Width along axis 0 [pixels].
 * @ysize: Height along axis 1 [pixels].
 *
 * Return: Newly allocated imageID, or -1 on failure.
 */
static inline imageID atmturb_create_2Dimage_ID_double(
    const char *name,
    uint32_t    xsize,
    uint32_t    ysize)
{
    imageID ID = -1;
    create_2Dimage_ID_double(name, xsize, ysize, &ID);
    return ID;
}

/**
 * atmturb_create_2DCimage_ID - Create a 2D complex single-precision image
 * @name:  Image name to allocate in milk_data.image.
 * @xsize: Width along axis 0 [pixels].
 * @ysize: Height along axis 1 [pixels].
 *
 * Return: Newly allocated imageID, or -1 on failure.
 */
static inline imageID atmturb_create_2DCimage_ID(
    const char *name,
    uint32_t    xsize,
    uint32_t    ysize)
{
    imageID ID = -1;
    create_2DCimage_ID(name, xsize, ysize, &ID);
    return ID;
}

/**
 * atmturb_create_2DCimage_ID_double - Create a 2D complex double image
 * @name:  Image name to allocate in milk_data.image.
 * @xsize: Width along axis 0 [pixels].
 * @ysize: Height along axis 1 [pixels].
 *
 * Return: Newly allocated imageID, or -1 on failure.
 */
static inline imageID atmturb_create_2DCimage_ID_double(
    const char *name,
    uint32_t    xsize,
    uint32_t    ysize)
{
    imageID ID = -1;
    create_2DCimage_ID_double(name, xsize, ysize, &ID);
    return ID;
}

/**
 * atmturb_create_3Dimage_ID - Create a shared 3D float image cube
 * @name:  Stream/image name to allocate.
 * @xsize: Width along axis 0 [pixels].
 * @ysize: Height along axis 1 [pixels].
 * @zsize: Depth along axis 2 [frames/slices].
 *
 * Allocates a 3D _DATATYPE_FLOAT cube with shared=1 and 10 keyword slots so
 * phase screens and wavefront cubes are accessible across processes.
 *
 * Return: Newly allocated imageID, or -1 on failure.
 */
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

/**
 * atmturb_create_3Dimage_ID_double - Create a 3D double-precision image cube
 * @name:  Image name to allocate in milk_data.image.
 * @xsize: Width along axis 0 [pixels].
 * @ysize: Height along axis 1 [pixels].
 * @zsize: Depth along axis 2 [frames/slices].
 *
 * Return: Newly allocated imageID, or -1 on failure.
 */
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

/**
 * atmturb_create_image_ID - Create an N-dimensional image from a long* size array
 * @name:     Image or shared stream name.
 * @naxis:    Number of dimensions (1..8).
 * @size:     Pointer to array of @naxis dimension sizes (stored as long).
 * @datatype: ImageStreamIO _DATATYPE_* code.
 * @shared:   1 to back with POSIX shared memory, 0 for local heap memory.
 * @nbkw:     Number of metadata keyword slots to allocate.
 *
 * Return: Newly allocated imageID, or -1 on failure.
 */
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

/**
 * atmturb_load_fits_2 - Load a FITS file into memory ignoring missing-file errors
 * @file: Path to input FITS file on disk.
 * @name: Destination image name in milk_data.image.
 *
 * Return: Loaded imageID (>= 0), or -1 on failure.
 */
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

/**
 * atmturb_load_fits_3 - Load a FITS file with explicit error-handling mode
 * @file:    Path to input FITS file on disk.
 * @name:    Destination image name in milk_data.image.
 * @errmode: LOADFITS_ERRMODE_* flag controlling failure reporting.
 *
 * Return: Loaded imageID (>= 0), or -1 on failure.
 */
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

/**
 * atmturb_linopt_compute_SVDpseudoInverse - 5-argument SVD pseudoinverse shim
 * @r:     Input matrix image name.
 * @c:     Output pseudoinverse matrix image name.
 * @eps:   Singular-value truncation threshold.
 * @maxnb: Maximum number of singular modes to retain.
 * @vt:    Output V^T matrix image name.
 *
 * Return: RETURN_SUCCESS on completion, or error code from milklinoptimtools.
 */
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

/**
 * arith_image_zero - Zero all bytes of an in-memory image buffer
 * @name: Target image name in milk_data.image.
 *
 * Return: 0 on success, -1 if @name is not found.
 */
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

/*
 * Argument-count dispatch macros routing short-form calls to the helpers above
 * while forwarding full-argument calls directly to the core milk functions.
 */
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
