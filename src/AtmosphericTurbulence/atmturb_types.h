// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_types.h
 * @brief   Internal types, global configuration state, and PRNG helpers
 *
 * Declares shared simulation configuration globals (CONF_*), physical constants,
 * legacy layer records, internal helper prototypes, and deterministic SplitMix64
 * / Box-Muller pseudorandom number generators used across AtmosphericTurbulence.
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

#include "milk_runnable.h"
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

/* Fundamental physical constants (SI units) */
extern double C_me;     // Electron rest mass [kg]
extern double C_e0;     // Vacuum permittivity [F/m]
extern double C_Na;     // Avogadro constant [mol^-1]
extern double C_e;      // Elementary charge [C]
extern double C_ls;     // Speed of light in vacuum [m/s]
extern double rhocoeff; // Classical electron radius / refractivity constant

/* Observing & site configuration globals */
extern char     CONFFILE[200];                  // Active WFsim.conf file path
extern float    CONF_LAMBDA;                    // Reference seeing wavelength [m]
extern float    CONF_SEEING;                    // Zenith seeing FWHM at CONF_LAMBDA [arcsec]
extern char     CONF_TURBULENCE_PROF_FILE[200]; // Path to multi-layer Cn2 profile file
extern float    CONF_ZANGLE;                    // Observation zenith angle [rad]
extern float    CONF_PARALLACTIC_ANGLE;         // Parallactic angle [rad]
extern float    CONF_SITE_ALT;                  // Observatory elevation ASL [m] (-1 = auto)
extern uint64_t CONF_SEED;                      // Master PRNG seed (0 = clock-seeded)
extern float    CONF_SOURCE_Xpos;               // Off-axis source X angle [rad]
extern float    CONF_SOURCE_Ypos;               // Off-axis source Y angle [rad]

/* Primary wavefront output & 2D SHM streaming configuration */
extern int  CONF_WFOUTPUT;           // FITS output bitmask (0=none, 1=ref, 2=sci, 3=both)
extern char CONF_WF_FILE_PREFIX[200];// Output FITS file prefix for reference cubes
extern char CONF_WF_PHASE_NAME[100]; // Output reference phase stream/image name
extern char CONF_WF_AMPL_NAME[100];  // Output reference amplitude stream/image name
extern int  CONF_SHM_OUTPUT;         // Publish 3D cubes to shared memory (0/1)
extern int  CONF_STREAM_MODE;        // 0=3D cube, 1=continuous 2D, 2=finite 2D, 3=unpaced 2D
extern int  CONF_STREAM_CBSIZE;      // ImageStreamIO circular history depth [frames]

/* Secondary chromatic wavefront output configuration */
extern int  CONF_MAKE_SWAVEFRONT;     // Synthesize secondary wavelength channel (0/1)
extern int  CONF_SWF_WRITE2DISK;      // Save secondary FITS cubes to disk (0/1)
extern char CONF_SWF_FILE_PREFIX[200];// Output FITS file prefix for secondary cubes
extern int  CONF_SHM_SOUTPUT;         // Publish secondary streams to shared memory (0/1)
extern char CONF_SHM_SPREFIX[100];    // Stream name prefix for secondary channel
extern int  CONF_SHM_SOUTPUTM;        // Secondary SHM format mode

/* Grid sampling & temporal pacing configuration */
extern long   CONF_WFsize;               // Synthesized pupil grid linear size [pixels]
extern float  CONF_PUPIL_SCALE;          // Pupil spatial sampling scale [m/pixel]
extern int    CONF_ATMWF_REALTIME;       // Enable wall-clock frame pacing (0/1)
extern double CONF_ATMWF_REALTIMEFACTOR; // Real-time pacing speed multiplier
extern float  CONF_WFTIME_STEP;          // Physics time step per frame [s]
extern float  CONF_TIME_SPAN;            // Duration per wavefront cube block [s]
extern long   CONF_NB_TSPAN;             // Number of consecutive time-span blocks
extern long   CONF_SIMTDELAY;            // Artificial per-frame delay [us]
extern int    CONF_WAITFORSEM;           // Synchronize frame emission on external semaphore (0/1)
extern char   CONF_WAITSEMIMNAME[100];   // External trigger stream name for CONF_WAITFORSEM

/* Master screen synthesis & extrusion fidelity configuration */
extern int   CONF_SKIP_EXISTING; // Reuse existing FITS cubes on disk when present (0/1)
extern long  CONF_WF_RAW_SIZE;   // Unbinned raw wavefront grid size [pixels]
extern long  CONF_MASTER_SIZE;   // Master Fourier phase screen linear size [pixels]
extern int   CONF_OVERSAMPLE;    // Master screen spatial oversampling factor (1 or 2)
extern int   CONF_INTERP;        // Sub-pixel kernel (0=bilinear, 1=Keys bicubic)
extern int   CONF_LOWFREQ;       // Enable 24-mode Johansson-Gavel subharmonics (0/1)
extern int   CONF_ROLLING;       // Enable multi-epoch rolling screen cross-fade (0/1)
extern float CONF_BOIL_TIME;     // Maximum rolling screen epoch duration [s] (0=auto)

/* Diffractive propagation (Fresnel & Rytov) configuration */
extern int   CONF_FRESNEL_PROPAGATION;     // 0=geometric, 1=Fresnel split-step, 2=Rytov
extern int   CONF_WAVEFRONT_AMPLITUDE;     // Compute scintillation amplitude screens (0/1)
extern float CONF_FRESNEL_PROPAGATION_BIN; // Super-layer distance binning threshold [m]
extern int   CONF_FRESNEL_RYTOV_SEC_EXACT; // Force separate secondary Rytov FFT (0/1)
extern int   CONF_FRESNEL_GUARD_PIX;       // Guard-band margin around pupil [pixels]
extern int   CONF_FRESNEL_RYTOV_MOISAN;    // Moisan periodic+smooth (-1=auto, 0=off, 1=on)
extern int   CONF_FRESNEL_RYTOV_ZINTERP;   // Linear distance-node interpolation in Rytov (0/1)
extern int   CONF_FRESNEL_REFRACT_PATH;    // Use curved refracted optical path length (0/1)
extern int   CONF_FRESNEL_SCINT_WEIGHT;    // z^(5/6) scintillation-weighted centroids (0/1)

/**
 * struct atmturb_layer_info_t - Parsed atmospheric layer record from config/profile
 * @alt:         Layer altitude above sea level [m].
 * @cn2:         Layer Cn2 integral weight or fractional strength.
 * @spd:         Mean horizontal wind speed [m/s].
 * @dir:         Wind direction azimuth [rad].
 * @outerscale:  Turbulence outer scale L0 [m].
 * @innerscale:  Turbulence inner scale l0 [m].
 * @sigmawspeed: Turbulent wind speed standard deviation [m/s].
 * @lwind:       Wind velocity turbulence outer scale [m].
 */
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

/**
 * AtmosphericTurbulence_ReadConf - Parse active CONFFILE into CONF_* globals
 *
 * Return: 0 on success, -1 if CONFFILE cannot be opened.
 */
int AtmosphericTurbulence_ReadConf(void);

/**
 * Z_Air - Compressibility factor of moist air (CIPM-81 / Picard et al. 2008)
 * @P:  Total atmospheric pressure [Pa].
 * @T:  Air temperature [deg C].
 * @RH: Relative humidity [%].
 *
 * Return: Dimensionless compressibility factor Z (also updates global rhocoeff).
 */
double Z_Air(
    double P,
    double T,
    double RH);

/**
 * Z_N2 - Virial compressibility factor of pure molecular nitrogen
 * @P: Gas pressure [Pa].
 * @T: Thermodynamic temperature [K].
 *
 * Return: Dimensionless compressibility factor Z for N2.
 */
double Z_N2(
    double P,
    double T);

/**
 * atmturb_screen_create_dist - Populate 2D radial spatial-frequency distance map
 * @name:      Destination image name in milk_data.image.
 * @size:      Square 2D grid dimension [pixels].
 * @outer_f0:  Outer-scale spatial frequency 1/L0 [cycles/pixel].
 * @rlim:      Radial cutoff distance [pixels].
 * @rlim_mode: Cutoff behavior flag (0=none, 1=zero inside/outside limit).
 * @precision: Output datatype selector (0=float, 1=double).
 *
 * Return: Allocated imageID (>= 0) on success, or -1 on failure.
 */
int atmturb_screen_create_dist(
    const char *name,
    long        size,
    double      outer_f0,
    double      rlim,
    int         rlim_mode,
    long        precision);

/**
 * atmturb_bin_2d_pixels - Spatially average-bin a 2D float frame by 2^factor
 * @src:    Input 2D float buffer (in_x * in_y elements).
 * @dst:    Output 2D float buffer ((in_x >> factor) * (in_y >> factor) elements).
 * @in_x:   Input width [pixels].
 * @in_y:   Input height [pixels].
 * @factor: Power-of-two decimation exponent (>= 0).
 *
 * Return: 0 on success, -1 on invalid dimensions or factor.
 */
int atmturb_bin_2d_pixels(
    const float *src,
    float       *dst,
    long         in_x,
    long         in_y,
    int          factor);

/**
 * atmturb_rng_splitmix64 - Advance 64-bit SplitMix PRNG state
 * @state: Pointer to 64-bit state variable.
 *
 * Return: Pseudorandom 64-bit unsigned integer.
 */
static inline uint64_t atmturb_rng_splitmix64(
    uint64_t *state)
{
    uint64_t z = (*state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

/**
 * atmturb_rng_gaussian_pair - Sample standard normal pair via Box-Muller transform
 * @state: Pointer to PRNG state.
 * @g0:    Output pointer for first standard Gaussian sample.
 * @g1:    Output pointer for second standard Gaussian sample.
 */
static inline void atmturb_rng_gaussian_pair(
    uint64_t *state,
    double   *g0,
    double   *g1)
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
static inline uint64_t atmturb_resolve_seed(
    uint64_t seed)
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
 * @seed:   Base seed.
 * @stream: Sub-stream index (component, layer, row, ...).
 *
 * Return: Initial PRNG state for the sub-stream.
 */
static inline uint64_t atmturb_rng_stream_seed(
    uint64_t seed,
    uint64_t stream)
{
    uint64_t s = seed ^ ((stream + 1ULL) * 0xd1b54a32d192ed03ULL);
    atmturb_rng_splitmix64(&s);
    return atmturb_rng_splitmix64(&s);
}

#endif // ATMTURB_TYPES_H
