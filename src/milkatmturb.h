// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    milkatmturb.h
 * @brief   Umbrella public header for the milkatmturb plugin
 *
 * Exposes C23/C11 portability macros together with the public APIs of
 * AtmosphereModel (refractivity and ray tracing), AtmosphericTurbulence
 * (phase synthesis, propagation, streaming, and AO metrics), OpticsMaterials,
 * and WFpropagate.
 */

#ifndef _MILK_ATMTURB_H
#define _MILK_ATMTURB_H

#include "atmturb_c23_compat.h"
#include "AtmosphereModel/AtmosphereModel.h"
#include "AtmosphericTurbulence/AtmosphericTurbulence.h"
#include "OpticsMaterials/OpticsMaterials.h"
#include "WFpropagate/WFpropagate.h"

#endif // _MILK_ATMTURB_H
