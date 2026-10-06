// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    WFpropagate.c
 * @brief   Wavefront propagation module initialization and CLI registration
 */

#include "CLIcore.h"
#include "WFpropagate.h"

extern DATA data;

/**
 * Fresnel_propagate_wavefront_cli - CLI wrapper for Fresnel propagation
 *
 * Return: 0 on success, 1 on argument check failure.
 */
static int Fresnel_propagate_wavefront_cli(void)
{
    if (CLI_checkarg(1, 4) + CLI_checkarg(2, 3) + CLI_checkarg(3, 1) + CLI_checkarg(4, 1) +
            CLI_checkarg(5, 1) ==
        0)
    {
        Fresnel_propagate_wavefront(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string,
                                    data.cmdargtoken[3].val.numf, data.cmdargtoken[4].val.numf,
                                    data.cmdargtoken[5].val.numf);
        return 0;
    }

    return 1;
}

/**
 * init_WFpropagate - Initialize module and register CLI commands
 *
 * Return: 0 on success.
 */
int init_WFpropagate(void)
{
    RegisterCLIcommand("fresnelpw", __FILE__, Fresnel_propagate_wavefront_cli,
                       "Fresnel propagate wavefront",
                       "<input image> <output image> <pupil scale m/s> <prop dist> <lambda>",
                       "fresnelpw in out 0.01 1000 0.0000005",
                       "int Fresnel_propagate_wavefront(char *in, char *out, double PUPIL_SCALE, "
                       "double z, double lambda)");

    return 0;
}
