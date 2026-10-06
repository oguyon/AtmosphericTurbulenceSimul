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
extern float CONF_SOURCE_Xpos;
extern float CONF_SOURCE_Ypos;

extern int CONF_WFOUTPUT;
extern char CONF_WF_FILE_PREFIX[200];
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

#endif // ATMTURB_TYPES_H
