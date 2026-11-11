// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_mkmastert_FPS.h
 * @brief   MILK_RUNNABLE FPS compute unit header for master phase screens generation
 */

#ifndef ATMTURB_MKMASTERT_FPS_H
#define ATMTURB_MKMASTERT_FPS_H

#include "milk_runnable.h"

/**
 * runnable_atmturb_mkmastert - FPS runnable descriptor for master phase screens
 *
 * Exposes CLI command `atmturb.mkmastert` and standalone executable
 * `milk-fpsexec-atmturb-mkmastert` for Fourier synthesis of von Karman screens.
 */
extern const MILK_RUNNABLE runnable_atmturb_mkmastert;

#endif // ATMTURB_MKMASTERT_FPS_H
