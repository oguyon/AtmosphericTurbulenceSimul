// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_rolling.c
 * @brief   Multi-epoch rolling screen management and variance-preserving cross-fading
 *
 * Allocates pools of independent master phase screens and subharmonic amplitudes
 * per turbulent layer and blends consecutive epochs via quadrature weights
 * w_A = cos(theta), w_B = sin(theta) (w_A^2 + w_B^2 = 1) to prevent wrap-around
 * periodicity and simulate temporal boiling.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "AtmosphericTurbulence.h"
#include "atmturb_boiling.h"
#include "atmturb_geometry.h"
#include "atmturb_lowfreq.h"
#include "atmturb_profile.h"
#include "atmturb_rolling.h"
#include "atmturb_screens_batch.h"
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
    if (id < 0)
    {
        free(tmp);
        return -1;
    }
    memcpy(dcimg[id].array.F, tmp, sizeof(float) * (size_t) ntot);
    free(tmp);
    return id;
}

typedef struct
{
    atmturb_screen_spec_t *specs;
    float                **out_a;
    float                **out_b;
    int                   *slot_a;
    int                   *slot_b;
    const float          **host_tags;
    int                    npairs;
    int                    all_fresh;
} atmturb_rolling_batch_state_t;

/**
 * atmturb_rolling_calc_nscreens - Compute required master screen pool depth for one layer
 */
static int atmturb_rolling_calc_nscreens(
    int    rolling,
    double t_dec_s,
    long   nbframes,
    double time_step_s)
{
    if (!rolling)
    {
        return 1;
    }

    double sim_duration = (nbframes > 1) ? ((double) (nbframes - 1) * time_step_s) : 20.0;
    long   max_epoch    = (long) (sim_duration / t_dec_s);
    int    nscreens     = (int) (max_epoch + 2);
    if (nscreens < 4)
    {
        nscreens = 4;
    }
    if (CONF_STREAM_MODE > 0 && nscreens > 6)
    {
        nscreens = 6;
    }
    else if (nscreens > 32)
    {
        nscreens = 32;
    }
    return nscreens;
}

/**
 * atmturb_rolling_prepare_pair_ids - Reuse or allocate dcimg master screen buffers for a pair
 */
static int atmturb_rolling_prepare_pair_ids(
    int      k,
    int      s,
    int      has_second,
    long     msize,
    imageID *id0,
    imageID *id1,
    int     *reused)
{
    char name0[200];
    char name1[200];
    snprintf(name0, sizeof(name0), "turbm%02d_p%d", k, s);
    snprintf(name1, sizeof(name1), "turbm%02d_p%d", k, s + 1);

    if (CONF_SKIP_EXISTING)
    {
        *id0 = atmturb_rolling_ensure_screen(name0, msize);
        *id1 = has_second ? atmturb_rolling_ensure_screen(name1, msize) : 0;
        if (*id0 >= 0 && (!has_second || *id1 >= 0))
        {
            *reused = 1;
            return 0;
        }
    }

    *reused = 0;
    delete_image_ID(name0);
    *id0 = create_2Dimage_ID(name0, msize, msize);
    if (*id0 < 0)
    {
        return -1;
    }
    if (has_second)
    {
        delete_image_ID(name1);
        *id1 = create_2Dimage_ID(name1, msize, msize);
        if (*id1 < 0)
        {
            return -1;
        }
    }
    return 0;
}

/**
 * atmturb_rolling_queue_layer - Setup dcimg buffers, lowfreq modes, and batch queue for layer k
 */
static int atmturb_rolling_queue_layer(
    atmturb_rolling_layer_t       *rl,
    int                            k,
    int                            vram_base,
    const atmturb_profile_t       *prof,
    const atmturb_geom_t          *geom,
    long                           msize,
    long                           precision,
    uint64_t                       seed,
    atmturb_rolling_batch_state_t *bs)
{
    double L0_pix = (prof->layers[k].L0_m > 0.0) ? (prof->layers[k].L0_m / geom->dx_master_m) : 0.0;
    double l0_pix = (prof->layers[k].l0_m > 0.0) ? (prof->layers[k].l0_m / geom->dx_master_m) : 0.0;
    atmturb_lowfreq_init(&rl->lf_base, msize, geom->r0_ref_pix, L0_pix, l0_pix, seed);

    for (int s = 0; s < rl->nscreens; s += 2)
    {
        int     has_second = (s + 1 < rl->nscreens) ? 1 : 0;
        imageID id0        = -1;
        imageID id1        = -1;
        int     reused     = 0;
        if (atmturb_rolling_prepare_pair_ids(k, s, has_second, msize, &id0, &id1, &reused) != 0)
        {
            return -1;
        }

        rl->screens[s].image_id      = id0;
        rl->screens[s].data          = dcimg[id0].array.F;
        bs->host_tags[vram_base + s] = rl->screens[s].data;
        uint64_t seed0 = atmturb_rng_stream_seed(seed, (uint64_t) (10000 + k + s * 1000));
        atmturb_lowfreq_draw_amplitudes(rl->screens[s].are, rl->screens[s].aim,
                                        msize, geom->r0_ref_pix, L0_pix, l0_pix, seed0);

        if (has_second)
        {
            rl->screens[s + 1].image_id          = id1;
            rl->screens[s + 1].data              = dcimg[id1].array.F;
            bs->host_tags[vram_base + s + 1]     = rl->screens[s + 1].data;
            uint64_t seed1 = atmturb_rng_stream_seed(seed, (uint64_t) (10000 + k + (s + 1) * 1000));
            atmturb_lowfreq_draw_amplitudes(rl->screens[s + 1].are, rl->screens[s + 1].aim,
                                            msize, geom->r0_ref_pix, L0_pix, l0_pix, seed1);
        }

        if (reused)
        {
            bs->all_fresh = 0;
            continue;
        }

        int p = bs->npairs++;
        memset(&bs->specs[p], 0, sizeof(bs->specs[p]));
        bs->specs[p].size      = msize;
        bs->specs[p].r0_pix    = geom->r0_ref_pix;
        bs->specs[p].L0_pix    = L0_pix;
        bs->specs[p].l0_pix    = l0_pix;
        bs->specs[p].seed      = atmturb_rng_stream_seed(seed, (uint64_t) (k + s * 1000));
        bs->specs[p].precision = (int) precision;
        bs->out_a[p]           = rl->screens[s].data;
        bs->out_b[p]           = has_second ? rl->screens[s + 1].data : NULL;
        bs->slot_a[p]          = vram_base + s;
        bs->slot_b[p]          = has_second ? (vram_base + s + 1) : -1;
    }
    return 0;
}

/**
 * atmturb_rolling_free_batch_state - Release temporary batch queue buffers
 */
static void atmturb_rolling_free_batch_state(
    atmturb_rolling_batch_state_t *bs)
{
    free(bs->specs);
    free(bs->out_a);
    free(bs->out_b);
    free(bs->slot_a);
    free(bs->slot_b);
    free(bs->host_tags);
}

/**
 * atmturb_rolling_init - Allocate and generate rolling screens and low-order modes
 * @r:           Output rolling context to initialize.
 * @prof:        Active multi-layer atmospheric turbulence profile.
 * @geom:        Computed observing geometry and layer decorrelation times.
 * @msize:       Master phase screen linear dimension [pixels].
 * @nbframes:    Total simulation frame count (used to size epoch screen pool).
 * @time_step_s: Physics time step per frame [s].
 * @precision:   Master screen FFT precision (0=single, 1=double).
 * @seed:        Master PRNG seed (0 selects clock-derived seed).
 *
 * Return: 0 on success, -1 on allocation or screen synthesis failure.
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
    seed       = atmturb_resolve_seed(seed);
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

    int abs_boil = (geom->boil_mode < 0) ? -geom->boil_mode : geom->boil_mode;
    if (abs_boil == 1 || abs_boil == 2)
    {
        int rc = atmturb_boiling_init(r, prof, geom, msize, nbframes,
                                      time_step_s, precision, seed);
        if (rc != 0)
        {
            atmturb_rolling_free(r);
        }
        return rc;
    }

    int max_pairs     = 0;
    int total_screens = 0;
    for (int k = 0; k < r->nlayers; k++)
    {
        atmturb_rolling_layer_t *rl = &r->layers[k];
        double t_dec  = geom->layers[k].t_dec_s;
        rl->t_dec_s   = (t_dec > 0.0) ? t_dec : 10.0;
        rl->nscreens  = atmturb_rolling_calc_nscreens(r->rolling, rl->t_dec_s,
                                                      nbframes, time_step_s);
        rl->screens   = (atmturb_rolling_screen_t *) calloc((size_t) (rl->nscreens + 1),
                                                            sizeof(atmturb_rolling_screen_t));
        if (rl->screens == NULL)
        {
            atmturb_rolling_free(r);
            return -1;
        }
        max_pairs     += (rl->nscreens + 1) / 2;
        total_screens += rl->nscreens;
    }

    atmturb_rolling_batch_state_t bs;
    memset(&bs, 0, sizeof(bs));
    bs.all_fresh = 1;
    bs.specs     = (atmturb_screen_spec_t *) calloc((size_t) (max_pairs + 1), sizeof(*bs.specs));
    bs.out_a     = (float **) calloc((size_t) (max_pairs + 1), sizeof(float *));
    bs.out_b     = (float **) calloc((size_t) (max_pairs + 1), sizeof(float *));
    bs.slot_a    = (int *) calloc((size_t) (max_pairs + 1), sizeof(int));
    bs.slot_b    = (int *) calloc((size_t) (max_pairs + 1), sizeof(int));
    bs.host_tags = (const float **) calloc((size_t) (total_screens + 1), sizeof(float *));
    if (!bs.specs || !bs.out_a || !bs.out_b || !bs.slot_a || !bs.slot_b || !bs.host_tags)
    {
        atmturb_rolling_free_batch_state(&bs);
        atmturb_rolling_free(r);
        return -1;
    }

    int vram_off = 0;
    int ret      = 0;
    for (int k = 0; k < r->nlayers; k++)
    {
        if (atmturb_rolling_queue_layer(&r->layers[k], k, vram_off, prof, geom,
                                        msize, precision, seed, &bs) != 0)
        {
            ret = -1;
            break;
        }
        vram_off += r->layers[k].nscreens;
    }

    if (ret == 0 && bs.npairs > 0)
    {
        atmturb_screen_batch_req_t req;
        memset(&req, 0, sizeof(req));
        req.specs              = bs.specs;
        req.out_a              = bs.out_a;
        req.out_b              = bs.out_b;
        req.npairs             = bs.npairs;
        req.vram_slot_a        = bs.all_fresh ? bs.slot_a : NULL;
        req.vram_slot_b        = bs.all_fresh ? bs.slot_b : NULL;
        req.total_vram_screens = bs.all_fresh ? total_screens : 0;
        req.vram_host_tags     = bs.all_fresh ? bs.host_tags : NULL;
        ret = atmturb_generate_screen_pairs_batch(&req);
    }

    atmturb_rolling_free_batch_state(&bs);
    if (ret != 0)
    {
        atmturb_rolling_free(r);
    }
    return ret;
}

/**
 * atmturb_rolling_free - Release all resources associated with rolling context
 * @r: Rolling context whose per-layer screen descriptors are freed (safe on NULL).
 */
void atmturb_rolling_free(
    atmturb_rolling_t *r)
{
    if (r == NULL)
    {
        return;
    }
    if (r->boil != NULL)
    {
        atmturb_boiling_free(r->boil);
        r->boil = NULL;
    }
    if (r->layers == NULL)
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
 * @r:           Initialized rolling simulation context.
 * @layer_idx:   Turbulence layer index (0..r->nlayers - 1).
 * @t:           Zero-based simulation frame index.
 * @time_step_s: Physics time step per frame [s].
 * @out:         Output descriptor receiving active screen pointers, quadrature weights,
 *               and blended subharmonic mode amplitudes.
 */
void atmturb_rolling_get_frame(
    const atmturb_rolling_t *r,
    int                      layer_idx,
    long                     t,
    double                   time_step_s,
    atmturb_rolling_eval_t  *out)
{
    if (r->boil != NULL)
    {
        atmturb_boiling_get_frame(r, layer_idx, t, out);
        return;
    }

    const atmturb_rolling_layer_t *rl = &r->layers[layer_idx];
    if (!r->rolling || rl->nscreens <= 1)
    {
        out->scrA = rl->screens[0].data;
        out->scrB = NULL;
        out->wA   = 1.0f;
        out->wB   = 0.0f;
        out->idxA = 0;
        out->idxB = 0;
        if (r->lowfreq)
        {
            memcpy(out->are_eff, rl->screens[0].are, sizeof(out->are_eff));
            memcpy(out->aim_eff, rl->screens[0].aim, sizeof(out->aim_eff));
        }
        return;
    }

    double t_time = (double) t * time_step_s;
    long epoch = (long) (t_time / rl->t_dec_s);
    long idxA, idxB;
    if (epoch >= rl->nscreens - 1)
    {
        idxA = epoch % rl->nscreens;
        idxB = (epoch + 1) % rl->nscreens;
    }
    else
    {
        idxA = epoch;
        idxB = epoch + 1;
    }

    double t_intra = t_time - (double) epoch * rl->t_dec_s;
    double theta = (0.5 * M_PI) * (t_intra / rl->t_dec_s);
    float wA = (float) cos(theta);
    float wB = (float) sin(theta);

    out->scrA = rl->screens[idxA].data;
    out->scrB = rl->screens[idxB].data;
    out->wA   = wA;
    out->wB   = wB;
    out->idxA = (int) idxA;
    out->idxB = (int) idxB;

    if (r->lowfreq)
    {
        const float *are0 = rl->screens[idxA].are;
        const float *aim0 = rl->screens[idxA].aim;
        const float *are1 = rl->screens[idxB].are;
        const float *aim1 = rl->screens[idxB].aim;

        for (int m = 0; m < ATMTURB_LOWFREQ_NMODES; m++)
        {
            out->are_eff[m] = wA * are0[m] + wB * are1[m];
            out->aim_eff[m] = wA * aim0[m] + wB * aim1[m];
        }
    }
}

/**
 * atmturb_rolling_advance - Advance active Fourier boiling state up to frame target_t
 */
int atmturb_rolling_advance(
    atmturb_rolling_t *r,
    long               target_t,
    int                gpu_only)
{
    if (r == NULL || r->boil == NULL)
    {
        return 0;
    }
    return atmturb_boiling_advance(r, target_t, gpu_only);
}

/**
 * atmturb_rolling_next_window - Advance and compute next collision-free 3D cube frame window
 */
int atmturb_rolling_next_window(
    atmturb_rolling_t *r,
    long               t0,
    long               nbframes,
    long              *out_t1)
{
    return atmturb_boiling_next_window(r, t0, nbframes, out_t1);
}

/**
 * atmturb_rolling_is_windowed - Query whether 3D cube rendering requires multi-window stepping
 */
int atmturb_rolling_is_windowed(
    const atmturb_rolling_t *r)
{
    return (r != NULL && r->boil != NULL && r->boil->windowed) ? 1 : 0;
}
