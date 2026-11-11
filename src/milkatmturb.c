// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    milkatmturb.c
 * @brief   MILK_MODULE_RUNNABLES registry for the milkatmturb plugin
 *
 * Registers all FPS and CLI runnable components under the "atmturb" namespace
 * (mkwfs, mkhvturb, mkmastert, mkvonkarman, aoloop, mkatm, measr0,
 * linpred.build, linpred.apply, and contract) and declares runtime module
 * dependencies for milk-cli and milk.bash dynamic loading.
 */

#define MODULE_SHORTNAME_DEFAULT "atmturb"
#define MODULE_DESCRIPTION       "Atmospheric Turbulence Simulation"

#include "milk_runnable.h"
#include "milkatmturb.h"
#include "AtmosphericTurbulence/atmturb_aoloop_FPS.h"
#include "AtmosphericTurbulence/atmturb_cli.h"
#include "AtmosphericTurbulence/atmturb_mkhvturb_FPS.h"
#include "AtmosphericTurbulence/atmturb_mkmastert_FPS.h"
#include "AtmosphericTurbulence/atmturb_mkvonkarman_FPS.h"
#include "AtmosphericTurbulence/atmturb_mkwfs_FPS.h"

MODULE_DEPS("milkCOREMODmemory", "milkCOREMODarith", "milkCOREMODiofits",
            "milkCOREMODtools", "milkfft", "milkOpticsMaterials",
            "milkWFpropagate");

static const MILK_RUNNABLE *const runnables[] = {
    &runnable_atmturb_mkwfs,
    &runnable_atmturb_mkhvturb,
    &runnable_atmturb_mkmastert,
    &runnable_atmturb_mkvonkarman,
    &runnable_atmturb_aoloop,
    &runnable_atmturb_mkatm,
    &runnable_atmturb_measr0,
    &runnable_atmturb_linpred_build,
    &runnable_atmturb_linpred_apply,
    &runnable_atmturb_contract,
};

MILK_MODULE_RUNNABLES(milkatmturb, _module_deps, runnables);
