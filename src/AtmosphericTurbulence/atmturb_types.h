// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_types.h
 * @brief   Internal types, state, and prototypes for AtmosphericTurbulence module
 */

#ifndef ATMTURB_TYPES_H
#define ATMTURB_TYPES_H

#include <assert.h>
#include <ctype.h>
#include <malloc.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "CLIcore.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_tools/COREMOD_tools.h"
#include "OpticsMaterials/OpticsMaterials.h"
#include "WFpropagate/WFpropagate.h"
#include "atmturb_compat.h"
#include "fft/fft.h"
#include "image_basic/image_basic.h"
#include "image_filter/image_filter.h"
#include "image_gen/image_gen.h"
#include "info/info.h"
#include "linalgebra/linalgebra.h"
#include "linopt_imtools/linopt_imtools.h"
#include "psf/psf.h"
#include "statistic/statistic.h"
#include "atmturb_profile.h"
#include "atmturb_geometry.h"

#ifndef PI
#define PI 3.14159265358979323846264338328
#endif

// Physical constants
extern double C_me;
extern double C_e0;
extern double C_Na;
extern double C_e;
extern double C_ls;
extern double rhocoeff;

// Configuration variables
extern char CONFFILE[200];
extern float CONF_LAMBDA;
extern float CONF_SEEING;
extern char CONF_TURBULENCE_PROF_FILE[200];
extern float CONF_ZANGLE;
extern float CONF_PARALLACTIC_ANGLE;
extern float CONF_SITE_ALT;
extern uint64_t CONF_SEED;
extern float CONF_SOURCE_Xpos;
extern float CONF_SOURCE_Ypos;

extern int CONF_WFOUTPUT;
extern char CONF_WF_FILE_PREFIX[200];
extern char CONF_WF_PHASE_NAME[100];
extern char CONF_WF_AMPL_NAME[100];
extern int CONF_SHM_OUTPUT;

extern int CONF_MAKE_SWAVEFRONT;
extern int CONF_SWF_WRITE2DISK;
extern char CONF_SWF_FILE_PREFIX[200];
extern int CONF_SHM_SOUTPUT;
extern char CONF_SHM_SPREFIX[100];
extern int CONF_SHM_SOUTPUTM;

extern long CONF_WFsize;
extern float CONF_PUPIL_SCALE;

extern int CONF_ATMWF_REALTIME;
extern double CONF_ATMWF_REALTIMEFACTOR;
extern float CONF_WFTIME_STEP;
extern float CONF_TIME_SPAN;
extern long CONF_NB_TSPAN;
extern long CONF_SIMTDELAY;
extern int CONF_WAITFORSEM;
extern char CONF_WAITSEMIMNAME[100];

extern int CONF_SKIP_EXISTING;
extern long CONF_WF_RAW_SIZE;
extern long CONF_MASTER_SIZE;
extern int CONF_OVERSAMPLE;
extern int CONF_INTERP;
extern int CONF_LOWFREQ;
extern int CONF_ROLLING;
extern float CONF_BOIL_TIME;

extern int CONF_FRESNEL_PROPAGATION;
extern int CONF_WAVEFRONT_AMPLITUDE;
extern float CONF_FRESNEL_PROPAGATION_BIN;

// Layer configuration structure
typedef struct
{
    double alt;
    double cn2;
    double spd;
    double dir;
    double outerscale;
    double innerscale;
    double sigmawspeed;
    double lwind;
} atmturb_layer_info_t;

// Module internal prototypes
int AtmosphericTurbulence_ReadConf(void);
double Z_Air(double P, double T, double RH);
double Z_N2(double P, double T);

int atmturb_screen_create_dist(const char *name, long size, double outer_f0,
                              double rlim, int rlim_mode, long precision);

int atmturb_bin_2d_pixels(const float *src, float *dst, long in_x, long in_y, int factor);

/**
 * atmturb_rng_splitmix64 - Advance 64-bit SplitMix PRNG state
 * @state: Pointer to 64-bit state variable.
 *
 * Return: Pseudorandom 64-bit unsigned integer.
 */
static inline uint64_t atmturb_rng_splitmix64(uint64_t *state)
{
    uint64_t z = (*state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

/**
 * atmturb_rng_gaussian_pair - Sample standard normal pair via Box-Muller transform
 * @state: Pointer to PRNG state.
 * @g0: Output pointer for first standard Gaussian sample.
 * @g1: Output pointer for second standard Gaussian sample.
 */
static inline void atmturb_rng_gaussian_pair(uint64_t *state, double *g0, double *g1)
{
    double u1 = ((atmturb_rng_splitmix64(state) >> 11) + 1.0) * (1.0 / 9007199254740992.0);
    double u2 = (atmturb_rng_splitmix64(state) >> 11) * (1.0 / 9007199254740992.0);
    double r = sqrt(-2.0 * log(u1));
    double theta = 2.0 * M_PI * u2;
    *g0 = r * cos(theta);
    *g1 = r * sin(theta);
}

/**
 * atmturb_resolve_seed - Map a user seed to a non-zero RNG seed
 * @seed: User seed; 0 requests a time-based seed.
 *
 * Return: @seed if non-zero, otherwise a seed derived from the real-time clock.
 */
static inline uint64_t atmturb_resolve_seed(uint64_t seed)
{
    if (seed != 0)
    {
        return seed;
    }
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    uint64_t s = (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
    return atmturb_rng_splitmix64(&s) | 1ULL;
}

/**
 * atmturb_rng_stream_seed - Derive a decorrelated RNG state for sub-stream @stream
 * @seed: Base seed.
 * @stream: Sub-stream index (component, layer, row, ...).
 *
 * Return: Initial PRNG state for the sub-stream.
 */
static inline uint64_t atmturb_rng_stream_seed(uint64_t seed, uint64_t stream)
{
    uint64_t s = seed ^ ((stream + 1ULL) * 0xd1b54a32d192ed03ULL);
    atmturb_rng_splitmix64(&s);
    return atmturb_rng_splitmix64(&s);
}

#endif // ATMTURB_TYPES_H
