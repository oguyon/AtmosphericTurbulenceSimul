// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_mkwfs_FPS.h
 * @brief   MILK_RUNNABLE FPS compute unit header for wavefront series generation
 */

#ifndef ATMTURB_MKWFS_FPS_H
#define ATMTURB_MKWFS_FPS_H

#include "milk_runnable.h"

/**
 * runnable_atmturb_mkwfs - FPS runnable descriptor for wavefront generation
 *
 * Exposes CLI command `atmturb.mkwfs` and standalone executable
 * `milk-fpsexec-atmturb-mkwfs` for 3D cube and 2D SHM wavefront simulation.
 */
extern const MILK_RUNNABLE runnable_atmturb_mkwfs;

#endif // ATMTURB_MKWFS_FPS_H
