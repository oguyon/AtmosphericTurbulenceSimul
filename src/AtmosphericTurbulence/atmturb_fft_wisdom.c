// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_fft_wisdom.c
 * @brief   Thread-safe and process-atomic FFTW3 single-precision wisdom caching
 */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <fftw3.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "atmturb_fft_wisdom.h"

static int      s_wisdom_loaded = 0;
static uint64_t s_wisdom_hash   = 0;
static size_t   s_wisdom_len    = 0;

/**
 * atmturb_wisdom_fnv1a - Compute 64-bit FNV-1a digest of an ASCII wisdom string
 * @str: Null-terminated wisdom string.
 * @len: Output byte length of string.
 *
 * Return: 64-bit FNV-1a hash.
 */
static uint64_t atmturb_wisdom_fnv1a(
    const char *str,
    size_t     *len)
{
    uint64_t h = 14695981039346656037ULL;
    size_t   n = 0;
    if (str != NULL)
    {
        while (str[n] != '\0')
        {
            h ^= (uint64_t) (unsigned char) str[n];
            h *= 1099511628211ULL;
            n++;
        }
    }
    *len = n;
    return h;
}

/**
 * atmturb_get_wisdom_path - Resolve persistent FFTW wisdom file path
 * @path:   Destination string buffer.
 * @maxlen: Buffer capacity in bytes.
 */
static void atmturb_get_wisdom_path(
    char  *path,
    size_t maxlen)
{
    const char *env_path = getenv("ATMTURB_FFTW_WISDOM");
    if (env_path != NULL && env_path[0] != '\0')
    {
        snprintf(path, maxlen, "%s", env_path);
        return;
    }

    const char *home = getenv("HOME");
    if (home != NULL && home[0] != '\0')
    {
        char milk_dir[512];
        snprintf(milk_dir, sizeof(milk_dir), "%s/.milk", home);
        (void) mkdir(milk_dir, 0755);
        snprintf(path, maxlen, "%s/.milk/atmturb_fftw_wisdom.txt", home);
        return;
    }

    snprintf(path, maxlen, "/tmp/atmturb_fftw_wisdom_%lu.txt", (unsigned long) getuid());
}

/**
 * atmturb_fftwf_plan_flags - Query FFTW planner flags for persistent 2D/1D transforms
 *
 * Return: FFTW_MEASURE by default, or FFTW_ESTIMATE if ATMTURB_FFTW_MEASURE=0.
 */
unsigned int atmturb_fftwf_plan_flags(void)
{
    const char *env_meas = getenv("ATMTURB_FFTW_MEASURE");
    if (env_meas != NULL && env_meas[0] != '\0' && atoi(env_meas) == 0)
    {
        return FFTW_ESTIMATE;
    }
    return FFTW_MEASURE;
}

/**
 * atmturb_fftwf_wisdom_load - Import persistent single-precision FFTW wisdom from disk
 */
void atmturb_fftwf_wisdom_load(void)
{
#ifdef _OPENMP
    #pragma omp critical(milkatmturb_fftw_wisdom)
#endif
    {
        if (!s_wisdom_loaded)
        {
            char path[512];
            atmturb_get_wisdom_path(path, sizeof(path));
            (void) fftwf_import_wisdom_from_filename(path);

            char *wstr = fftwf_export_wisdom_to_string();
            if (wstr != NULL)
            {
                s_wisdom_hash = atmturb_wisdom_fnv1a(wstr, &s_wisdom_len);
                fftwf_free(wstr);
            }
            s_wisdom_loaded = 1;
        }
    }
}

/**
 * atmturb_fftwf_wisdom_save - Atomically export updated FFTW wisdom if new plans were measured
 */
void atmturb_fftwf_wisdom_save(void)
{
#ifdef _OPENMP
    #pragma omp critical(milkatmturb_fftw_wisdom)
#endif
    {
        char *wstr = fftwf_export_wisdom_to_string();
        if (wstr != NULL)
        {
            size_t   len = 0;
            uint64_t h   = atmturb_wisdom_fnv1a(wstr, &len);
            if (len > 0 && (len != s_wisdom_len || h != s_wisdom_hash))
            {
                char path[512];
                char tmp_path[544];
                atmturb_get_wisdom_path(path, sizeof(path));
                snprintf(tmp_path, sizeof(tmp_path), "%s.tmp.XXXXXX", path);

                int fd = mkstemp(tmp_path);
                if (fd >= 0)
                {
                    ssize_t written = write(fd, wstr, len);
                    close(fd);
                    if (written == (ssize_t) len && rename(tmp_path, path) == 0)
                    {
                        s_wisdom_hash = h;
                        s_wisdom_len  = len;
                    }
                    else
                    {
                        unlink(tmp_path);
                    }
                }
            }
            fftwf_free(wstr);
        }
    }
}
