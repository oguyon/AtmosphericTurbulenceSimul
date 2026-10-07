// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_rolling.c
 * @brief   Multi-epoch rolling screen management and cross-fading
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "AtmosphericTurbulence.h"
#include "atmturb_geometry.h"
#include "atmturb_lowfreq.h"
#include "atmturb_profile.h"
#include "atmturb_rolling.h"
#include "atmturb_types.h"

/**
 * atmturb_rolling_ensure_screen - Validate or create master screen image in FP32
 * @name: Master screen image name.
 * @msize: Screen linear dimension in pixels.
 *
 * Return: Valid Image ID, or -1 on failure.
 */
static imageID atmturb_rolling_ensure_screen(
    const char *name,
    long        msize)
{
    imageID id = image_ID(name);
    if (id < 0)
    {
        return -1;
    }
    if ((long) dcimg[id].md[0].size[0] != msize || (long) dcimg[id].md[0].size[1] != msize)
    {
        return -1;
    }
    if (dcimg[id].md[0].datatype == _DATATYPE_FLOAT)
    {
        return id;
    }
    if (dcimg[id].md[0].datatype != _DATATYPE_DOUBLE)
    {
        return -1;
    }

    long ntot = msize * msize;
    float *tmp = (float *) malloc(sizeof(float) * (size_t) ntot);
    if (tmp == NULL)
    {
        return -1;
    }
    for (long ii = 0; ii < ntot; ii++)
    {
        tmp[ii] = (float) dcimg[id].array.D[ii];
    }
    delete_image_ID(name);
    id = create_2Dimage_ID(name, msize, msize);
    memcpy(dcimg[id].array.F, tmp, sizeof(float) * (size_t) ntot);
    free(tmp);
    return id;
}

/**
 * atmturb_rolling_init_screen_pair - Synthesize or load a pair of master screens
 * @name0: Name of screen 0.
 * @name1: Name of screen 1.
 * @spec: Screen generation specifications.
 * @id0: Output image ID for screen 0.
 * @id1: Output image ID for screen 1.
 *
 * Return: 0 on success, -1 on failure.
 */
static int atmturb_rolling_init_screen_pair(
    const char                  *name0,
    const char                  *name1,
    const atmturb_screen_spec_t *spec,
    imageID                     *id0,
    imageID                     *id1)
{
    long msize = spec->size;
    if (!CONF_SKIP_EXISTING || image_ID(name0) < 0 || image_ID(name1) < 0)
    {
        delete_image_ID(name0);
        delete_image_ID(name1);
        *id0 = create_2Dimage_ID(name0, msize, msize);
        *id1 = create_2Dimage_ID(name1, msize, msize);
        if (*id0 < 0 || *id1 < 0)
        {
            return -1;
        }
        if (atmturb_generate_screen_pair(spec, dcimg[*id0].array.F, dcimg[*id1].array.F) != 0)
        {
            return -1;
        }
    }
    *id0 = atmturb_rolling_ensure_screen(name0, msize);
    *id1 = atmturb_rolling_ensure_screen(name1, msize);
    return (*id0 >= 0 && *id1 >= 0) ? 0 : -1;
}

/**
 * atmturb_rolling_init_layer_screens - Allocate and populate screens for a single layer
 * @rl: Target layer rolling container.
 * @k: Layer index.
 * @prof: Active atmospheric profile.
 * @geom: Computed observing geometry.
 * @msize: Master screen linear dimension.
 * @precision: FFT precision mode.
 * @seed: Master PRNG seed.
 *
 * Return: 0 on success, -1 on failure.
 */
static int atmturb_rolling_init_layer_screens(
    atmturb_rolling_layer_t *rl,
    int                      k,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    long                     msize,
    long                     precision,
    uint64_t                 seed)
{
    double L0_pix = (prof->layers[k].L0_m > 0.0)
                        ? (prof->layers[k].L0_m / geom->dx_master_m) : 0.0;
    double l0_pix = (prof->layers[k].l0_m > 0.0)
                        ? (prof->layers[k].l0_m / geom->dx_master_m) : 0.0;

    atmturb_lowfreq_init(&rl->lf_base, msize, geom->r0_ref_pix, L0_pix, l0_pix, seed);

    for (int s = 0; s < rl->nscreens; s += 2)
    {
        char name0[200], name1[200];
        snprintf(name0, sizeof(name0), "turbm%02d_p%d", k, s);
        snprintf(name1, sizeof(name1), "turbm%02d_p%d", k, s + 1);

        atmturb_screen_spec_t spec;
        memset(&spec, 0, sizeof(spec));
        spec.size      = msize;
        spec.r0_pix    = geom->r0_ref_pix;
        spec.L0_pix    = L0_pix;
        spec.l0_pix    = l0_pix;
        spec.seed      = atmturb_rng_stream_seed(seed, (uint64_t) (k + s * 1000));
        spec.precision = (int) precision;

        imageID id0 = -1, id1 = -1;
        if (atmturb_rolling_init_screen_pair(name0, name1, &spec, &id0, &id1) != 0)
        {
            return -1;
        }

        rl->screens[s].image_id = id0;
        rl->screens[s].data     = dcimg[id0].array.F;
        uint64_t seed0 = atmturb_rng_stream_seed(seed, (uint64_t) (10000 + k + s * 1000));
        atmturb_lowfreq_draw_amplitudes(rl->screens[s].are, rl->screens[s].aim,
                                        msize, geom->r0_ref_pix, L0_pix, l0_pix, seed0);

        if (s + 1 < rl->nscreens)
        {
            rl->screens[s + 1].image_id = id1;
            rl->screens[s + 1].data     = dcimg[id1].array.F;
            uint64_t seed1 = atmturb_rng_stream_seed(seed, (uint64_t) (10000 + k + (s + 1) * 1000));
            atmturb_lowfreq_draw_amplitudes(rl->screens[s + 1].are, rl->screens[s + 1].aim,
                                            msize, geom->r0_ref_pix, L0_pix, l0_pix, seed1);
        }
    }
    return 0;
}

/**
 * atmturb_rolling_init - Allocate and generate rolling screens and low-order modes
 */
int atmturb_rolling_init(
    atmturb_rolling_t       *r,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    long                     msize,
    long                     nbframes,
    double                   time_step_s,
    long                     precision,
    uint64_t                 seed)
{
    memset(r, 0, sizeof(*r));
    r->nlayers = prof->nlayers;
    r->msize   = msize;
    r->rolling = geom->rolling;
    r->lowfreq = geom->lowfreq;

    r->layers = (atmturb_rolling_layer_t *) calloc((size_t) r->nlayers,
                                                   sizeof(atmturb_rolling_layer_t));
    if (r->layers == NULL)
    {
        return -1;
    }

    for (int k = 0; k < r->nlayers; k++)
    {
        atmturb_rolling_layer_t *rl = &r->layers[k];
        double t_dec = geom->layers[k].t_dec_s;
        rl->t_dec_s  = (t_dec > 0.0) ? t_dec : 10.0;

        int nscreens = 1;
        if (r->rolling)
        {
            double sim_duration = (double) (nbframes - 1) * time_step_s;
            long max_epoch = (long) (sim_duration / rl->t_dec_s);
            nscreens = (int) (max_epoch + 2);
        }
        rl->nscreens = nscreens;
        rl->screens  = (atmturb_rolling_screen_t *) calloc((size_t) (nscreens + 1),
                                                           sizeof(atmturb_rolling_screen_t));
        if (rl->screens == NULL)
        {
            atmturb_rolling_free(r);
            return -1;
        }

        if (atmturb_rolling_init_layer_screens(rl, k, prof, geom, msize, precision, seed) != 0)
        {
            atmturb_rolling_free(r);
            return -1;
        }
    }
    return 0;
}

/**
 * atmturb_rolling_free - Release all resources associated with rolling context
 */
void atmturb_rolling_free(
    atmturb_rolling_t *r)
{
    if (r == NULL || r->layers == NULL)
    {
        return;
    }
    for (int k = 0; k < r->nlayers; k++)
    {
        if (r->layers[k].screens != NULL)
        {
            free(r->layers[k].screens);
        }
    }
    free(r->layers);
    r->layers = NULL;
}

/**
 * atmturb_rolling_get_frame - Evaluate screen pointers, weights, and mode blending for frame t
 */
void atmturb_rolling_get_frame(
    const atmturb_rolling_t *r,
    int                      layer_idx,
    long                     t,
    double                   time_step_s,
    atmturb_rolling_eval_t  *out)
{
    const atmturb_rolling_layer_t *rl = &r->layers[layer_idx];
    if (!r->rolling || rl->nscreens <= 1)
    {
        out->scrA = rl->screens[0].data;
        out->scrB = NULL;
        out->wA   = 1.0f;
        out->wB   = 0.0f;
        if (r->lowfreq)
        {
            memcpy(out->are_eff, rl->screens[0].are, sizeof(out->are_eff));
            memcpy(out->aim_eff, rl->screens[0].aim, sizeof(out->aim_eff));
        }
        return;
    }

    double t_time = (double) t * time_step_s;
    long epoch = (long) (t_time / rl->t_dec_s);
    if (epoch >= rl->nscreens - 1)
    {
        epoch = rl->nscreens - 2;
    }

    double t_intra = t_time - (double) epoch * rl->t_dec_s;
    double theta = (0.5 * M_PI) * (t_intra / rl->t_dec_s);
    float wA = (float) cos(theta);
    float wB = (float) sin(theta);

    out->scrA = rl->screens[epoch].data;
    out->scrB = rl->screens[epoch + 1].data;
    out->wA   = wA;
    out->wB   = wB;

    if (r->lowfreq)
    {
        const float *are0 = rl->screens[epoch].are;
        const float *aim0 = rl->screens[epoch].aim;
        const float *are1 = rl->screens[epoch + 1].are;
        const float *aim1 = rl->screens[epoch + 1].aim;

        for (int m = 0; m < ATMTURB_LOWFREQ_NMODES; m++)
        {
            out->are_eff[m] = wA * are0[m] + wB * are1[m];
            out->aim_eff[m] = wA * aim0[m] + wB * aim1[m];
        }
    }
}
