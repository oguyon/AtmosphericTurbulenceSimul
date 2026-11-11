// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_mkhvturb_FPS.h
 * @brief   MILK_RUNNABLE FPS compute unit header for Hufnagel-Valley profile generation
 */

#ifndef ATMTURB_MKHVTURB_FPS_H
#define ATMTURB_MKHVTURB_FPS_H

#include "milk_runnable.h"

/**
 * runnable_atmturb_mkhvturb - FPS runnable descriptor for Hufnagel-Valley profiles
 *
 * Exposes CLI command `atmturb.mkhvturb` and standalone executable
 * `milk-fpsexec-atmturb-mkhvturb` to write discretized Cn2 and wind layer files.
 */
extern const MILK_RUNNABLE runnable_atmturb_mkhvturb;

#endif // ATMTURB_MKHVTURB_FPS_H
