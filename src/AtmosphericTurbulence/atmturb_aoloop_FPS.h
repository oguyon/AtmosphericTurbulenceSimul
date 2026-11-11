// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_aoloop_FPS.h
 * @brief   MILK_RUNNABLE FPS compute unit header for closed-loop AO PSF simulation
 */

#ifndef ATMTURB_AOLOOP_FPS_H
#define ATMTURB_AOLOOP_FPS_H

#include "milk_runnable.h"

/**
 * runnable_atmturb_aoloop - FPS runnable descriptor for closed-loop AO PSF synthesis
 *
 * Exposes CLI command `atmturb.aoloop` and standalone executable
 * `milk-fpsexec-atmturb-aoloop`; publishes cumulative Strehl to outvar `strehl`.
 */
extern const MILK_RUNNABLE runnable_atmturb_aoloop;

#endif // ATMTURB_AOLOOP_FPS_H
