// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cli.h
 * @brief   MILK_RUNNABLE CLI command descriptors for milkatmturb
 */

#ifndef ATMTURB_CLI_H
#define ATMTURB_CLI_H

#include "milk_runnable.h"

/** CLI runnable `atmturb.mkatm`: build NRLMSISE-00 atmosphere & refraction files. */
extern const MILK_RUNNABLE runnable_atmturb_mkatm;

/** CLI runnable `atmturb.measr0`: measure Fried parameter r0 [pix] (outvar `r0pix`). */
extern const MILK_RUNNABLE runnable_atmturb_measr0;

/** CLI runnable `atmturb.linpred_build`: train full-aperture AR linear predictor. */
extern const MILK_RUNNABLE runnable_atmturb_linpred_build;

/** CLI runnable `atmturb.linpred_apply`: apply trained linear predictor to cube. */
extern const MILK_RUNNABLE runnable_atmturb_linpred_apply;

/** CLI runnable `atmturb.contract`: spatially decimate wavefront FITS cube(s) by 2^k. */
extern const MILK_RUNNABLE runnable_atmturb_contract;

#endif // ATMTURB_CLI_H
