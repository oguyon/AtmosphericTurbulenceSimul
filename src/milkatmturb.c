// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    milkatmturb.c
 * @brief   Atmospheric Turbulence Simulation module for milk
 */

#define MODULE_SHORTNAME_DEFAULT "atmturb"
#define MODULE_DESCRIPTION       "Atmospheric Turbulence Simulation"

#include "CLIcore.h"
#include "milkatmturb.h"

errno_t CLIADDCMD_milkatmturb__atmturb_mkwfs_FPS(void);
errno_t CLIADDCMD_milkatmturb__atmturb_mkhvturb_FPS(void);

static errno_t init_module_CLI()
{
    init_AtmosphereModel();
    init_OpticsMaterials();
    init_WFpropagate();
    init_AtmosphericTurbulence();

    CLIADDCMD_milkatmturb__atmturb_mkwfs_FPS();
    CLIADDCMD_milkatmturb__atmturb_mkhvturb_FPS();

    return RETURN_SUCCESS;
}

MILK_MODULE(milkatmturb, init_module_CLI, NULL);
