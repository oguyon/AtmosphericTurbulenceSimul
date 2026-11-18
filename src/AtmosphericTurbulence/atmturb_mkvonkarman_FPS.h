// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_mkvonkarman_FPS.h
 * @brief   MILK_RUNNABLE FPS compute unit header for von Karman wind velocity simulation
 */

#ifndef ATMTURB_MKVONKARMAN_FPS_H
#define ATMTURB_MKVONKARMAN_FPS_H

#include "milk_runnable.h"

/**
 * runnable_atmturb_mkvonkarman - FPS runnable descriptor for von Karman wind cubes
 *
 * Exposes CLI command `atmturb.mkvonkarman` and standalone executable
 * `milk-fpsexec-atmturb-mkvonkarman` to synthesize 3-component [u, v, w] series.
 */
extern const MILK_RUNNABLE runnable_atmturb_mkvonkarman;

#endif // ATMTURB_MKVONKARMAN_FPS_H
