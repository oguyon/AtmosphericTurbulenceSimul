// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_boiling.c
 * @brief   Lifecycle, keyframe ring/window management, and blending for Fourier boiling
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>

#include "AtmosphericTurbulence.h"
#include "atmturb_boiling.h"
#include "atmturb_fft_wisdom.h"
#include "atmturb_screens_batch.h"
#include "atmturb_simd.h"
#include "atmturb_types.h"

/**
 * atmturb_boil_keyframe_seeds - Derive deterministic PRNG seeds for layer k at keyframe kf
 */
static void atmturb_boil_keyframe_seeds(
    uint64_t base_seed,
    int      k,
    long     kf,
    uint64_t out_seeds[2])
{
    if (kf <= 0)
    {
        out_seeds[0] = atmturb_rng_stream_seed(base_seed, (uint64_t) (200000 + k * 2));
        out_seeds[1] = atmturb_rng_stream_seed(base_seed, (uint64_t) (200000 + k * 2 + 1));
    }
    else
    {
        uint64_t tag = 1000000ULL + (uint64_t) k * 1000003ULL + (uint64_t) kf;
        out_seeds[0] = atmturb_rng_stream_seed(base_seed, tag);
        out_seeds[1] = 0ULL;
    }
}

/**
 * atmturb_boil_step_keyframe - Materialize initial or next AR keyframe into layer slot
 */
static int atmturb_boil_step_keyframe(
    atmturb_rolling_t *r,
    int                k,
    long               kf,
    int                gpu_only,
    float             *c2r_out)
{
    atmturb_boiling_t       *b    = r->boil;
    atmturb_boil_layer_t    *bl   = &b->layers[k];
    atmturb_rolling_layer_t *rl   = &r->layers[k];
    int                      slot = (int) (kf % (long) rl->nscreens);
    float                   *hdst = gpu_only ? NULL : rl->screens[slot].data;
    uint64_t                 seeds[2];

    atmturb_boil_keyframe_seeds(b->seed, k, kf, seeds);

    if (b->gpu != NULL)
    {
        if (atmturb_cuda_boiling_step_layer(b->gpu, k, slot, (kf == 0) ? 1 : 0,
                                            seeds, hdst) != 0)
        {
            return -1;
        }
        bl->max_kf_ready = kf;
        return 0;
    }

    if (kf == 0)
    {
        atmturb_boil_init_spectrum_cpu(bl, b->c2r_work, b->msize, b->ar_order, seeds);
    }
    else
    {
        atmturb_boil_step_spectrum_cpu(bl, b->c2r_work, b->msize, b->ar_order, seeds[0]);
    }

    fftwf_execute_dft_c2r(b->c2r_plan, b->c2r_work, c2r_out);
    memcpy(rl->screens[slot].data, c2r_out,
           sizeof(float) * (size_t) (b->msize * b->msize));
    bl->max_kf_ready = kf;
    return 0;
}

/**
 * atmturb_boil_config_layer - Compute cadence, spectral constants, and pool size for layer k
 */
static int atmturb_boil_config_layer(
    atmturb_boiling_t       *b,
    atmturb_rolling_layer_t *rl,
    atmturb_boil_layer_t    *bl,
    const atmturb_layer_t   *prof_l,
    const atmturb_geom_t    *geom,
    int                      k,
    long                     nbframes,
    double                   time_step_s)
{
    double sb     = geom->layers[k].sigma_boil_mps;
    double L0_pix = (prof_l->L0_m > 0.0) ? (prof_l->L0_m / geom->dx_master_m) : 0.0;
    double l0_pix = (prof_l->l0_m > 0.0) ? (prof_l->l0_m / geom->dx_master_m) : 0.0;

    atmturb_lowfreq_init(&rl->lf_base, b->msize, geom->r0_ref_pix, L0_pix, l0_pix, b->seed);
    bl->max_kf_ready = -1;
    bl->max_lf_ready = -1;
    bl->active       = (sb > 0.0) ? 1 : 0;

    if (!bl->active)
    {
        bl->k_sub    = 1;
        bl->lf_cap   = 1;
        bl->rho_bar  = 1.0;
        rl->nscreens = 1;
        return 0;
    }

    double eta_f   = (sqrt(2.0) * M_PI * sb * time_step_s) / geom->dx_master_m;
    double eta_lim = (b->ar_order == 2) ? 0.20 : 0.10;
    int    ks      = (geom->boil_mode < 0 || eta_f >= eta_lim || eta_f <= 0.0)
                         ? 1 : (int) floor(eta_lim / eta_f);
    bl->k_sub = (ks < 1) ? 1 : ((ks > 16) ? 16 : ks);

    if (CONF_STREAM_MODE > 0)
    {
        rl->nscreens = 4;
        bl->lf_cap   = 16;
    }
    else
    {
        long last_t   = (nbframes > 1) ? (nbframes - 1) : 0;
        long total_kf = (last_t % bl->k_sub > 0) ? (last_t / bl->k_sub + 2)
                                                 : (last_t / bl->k_sub + 1);
        int  pcap     = (b->msize <= 512) ? 64 : ((b->msize <= 1024) ? 32
                        : ((b->msize <= 2048) ? 16 : 8));
        rl->nscreens  = (total_kf <= (long) pcap) ? (int) total_kf : pcap;
        if (total_kf > (long) pcap)
        {
            b->windowed = 1;
        }
        bl->lf_cap = (total_kf <= (long) pcap) ? ((nbframes > 1) ? nbframes : 1)
                                               : (long) (rl->nscreens * bl->k_sub + 4);
    }

    float r0   = (geom->r0_ref_pix > 0.0) ? (float) geom->r0_ref_pix : powf(6.88f, 0.6f);
    float k0   = (L0_pix > 0.0) ? ((float) b->msize / (float) L0_pix) : 0.0f;
    float km   = (l0_pix > 0.0) ? ((5.92f / (2.0f * (float) M_PI)) * (float) b->msize
                                   / (float) l0_pix) : 0.0f;
    float pref = 0.023f * powf(r0, -5.0f / 3.0f) * powf((float) b->msize, 5.0f / 3.0f);
    double dt_kf       = (double) bl->k_sub * time_step_s;
    bl->k0_sq          = k0 * k0;
    bl->sqrt_pref_real = sqrtf(pref);
    bl->sqrt_pref_pair = sqrtf(0.5f * pref);
    bl->inv_two_km2    = (km > 0.0f) ? (0.5f / (km * km)) : 0.0f;
    bl->alpha_dt       = (float) ((2.0 * M_PI * sb * dt_kf)
                                  / ((double) b->msize * geom->dx_master_m));
    bl->rho_bar        = (bl->k_sub > 1)
                             ? atmturb_boil_compute_rho1(bl, b->msize, b->ar_order) : 1.0;
    return 0;
}

/**
 * atmturb_boil_alloc_screens_and_lf - Allocate dcimg buffers and subharmonic AR rings
 */
static int atmturb_boil_alloc_screens_and_lf(
    atmturb_rolling_t       *r,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    double                   time_step_s,
    long                     precision,
    const float            **vram_tags)
{
    atmturb_boiling_t *b = r->boil;
    for (int k = 0; k < r->nlayers; k++)
    {
        atmturb_rolling_layer_t *rl = &r->layers[k];
        atmturb_boil_layer_t    *bl = &b->layers[k];
        rl->screens = (atmturb_rolling_screen_t *) calloc((size_t) rl->nscreens,
                                                          sizeof(atmturb_rolling_screen_t));
        bl->lf_ring_re = (float *) calloc((size_t) (bl->lf_cap * ATMTURB_LOWFREQ_NMODES),
                                          sizeof(float));
        bl->lf_ring_im = (float *) calloc((size_t) (bl->lf_cap * ATMTURB_LOWFREQ_NMODES),
                                          sizeof(float));
        if (!rl->screens || !bl->lf_ring_re || !bl->lf_ring_im)
        {
            return -1;
        }

        for (int s = 0; s < rl->nscreens; s++)
        {
            char name[200];
            snprintf(name, sizeof(name), "turbm%02d_p%d", k, s);
            delete_image_ID(name);
            imageID id = create_2Dimage_ID(name, b->msize, b->msize);
            if (id < 0)
            {
                return -1;
            }
            rl->screens[s].image_id           = id;
            rl->screens[s].data               = dcimg[id].array.F;
            vram_tags[bl->vram_offset + s]    = rl->screens[s].data;
        }

        double L0_p = (prof->layers[k].L0_m > 0.0)
                          ? (prof->layers[k].L0_m / geom->dx_master_m) : 0.0;
        double l0_p = (prof->layers[k].l0_m > 0.0)
                          ? (prof->layers[k].l0_m / geom->dx_master_m) : 0.0;
        if (!bl->active)
        {
            atmturb_screen_spec_t sp = {
                .size      = b->msize,
                .r0_pix    = geom->r0_ref_pix,
                .L0_pix    = L0_p,
                .l0_pix    = l0_p,
                .seed      = atmturb_rng_stream_seed(b->seed, (uint64_t) k),
                .precision = (int) precision
            };
            if (atmturb_generate_screen_pair(&sp, rl->screens[0].data, NULL) != 0)
            {
                return -1;
            }
            uint64_t s0 = atmturb_rng_stream_seed(b->seed, (uint64_t) (10000 + k));
            atmturb_lowfreq_draw_amplitudes(rl->screens[0].are, rl->screens[0].aim,
                                            b->msize, geom->r0_ref_pix, L0_p, l0_p, s0);
        }
        else if (r->lowfreq)
        {
            atmturb_boil_lf_spec_t lfs = {
                .msize       = b->msize,
                .r0_pix      = geom->r0_ref_pix,
                .L0_pix      = L0_p,
                .l0_pix      = l0_p,
                .gamma_scale = (geom->layers[k].sigma_boil_mps * time_step_s)
                               / geom->dx_master_m,
                .ar_order    = b->ar_order,
                .seed0       = atmturb_rng_stream_seed(b->seed, (uint64_t) (300000 + k * 2)),
                .seed1       = atmturb_rng_stream_seed(b->seed, (uint64_t) (300000 + k * 2 + 1))
            };
            atmturb_boil_init_lowfreq_layer(bl, &rl->lf_base, &lfs);
        }
    }
    return 0;
}

/**
 * atmturb_boil_alloc_cpu_backend - Allocate CPU half-spectra and cached FFTW C2R plan
 */
static int atmturb_boil_alloc_cpu_backend(
    atmturb_boiling_t *b,
    float             *c2r_buf)
{
    size_t ncomp = (size_t) (b->msize * b->hcols);
    b->c2r_work  = (fftwf_complex *) fftwf_alloc_complex(ncomp);
    if (!b->c2r_work || !c2r_buf)
    {
        return -1;
    }

    for (int k = 0; k < b->nlayers; k++)
    {
        atmturb_boil_layer_t *bl = &b->layers[k];
        if (!bl->active)
        {
            continue;
        }
        bl->spec_cur = (fftwf_complex *) fftwf_alloc_complex(ncomp);
        if (b->ar_order == 2)
        {
            bl->spec_aux = (fftwf_complex *) fftwf_alloc_complex(ncomp);
        }
        if (!bl->spec_cur || (b->ar_order == 2 && !bl->spec_aux))
        {
            return -1;
        }
    }

    atmturb_fftwf_wisdom_load();
    atmturb_fftwf_plan_threads(0);
    b->c2r_plan = fftwf_plan_dft_c2r_2d((int) b->msize, (int) b->msize,
                                        b->c2r_work, c2r_buf,
                                        atmturb_fftwf_plan_flags());
    atmturb_fftwf_wisdom_save();
    return (b->c2r_plan != NULL) ? 0 : -1;
}

/**
 * atmturb_boil_try_gpu_backend - Attempt GPU cuFFT boiling backend initialization
 */
static void atmturb_boil_try_gpu_backend(
    atmturb_boiling_t *b,
    int                total_vram,
    long               precision)
{
    if (precision != 0 || !atmturb_simd_is_gpu() ||
        (!atmturb_simd_is_gpu_forced() && b->msize < 512 && CONF_WFsize < 128))
    {
        return;
    }

    atmturb_cuda_boil_layer_spec_t *gs = (atmturb_cuda_boil_layer_spec_t *) calloc(
        (size_t) b->nlayers, sizeof(*gs));
    if (gs == NULL)
    {
        return;
    }
    for (int k = 0; k < b->nlayers; k++)
    {
        gs[k].k0_sq          = b->layers[k].k0_sq;
        gs[k].sqrt_pref_pair = b->layers[k].sqrt_pref_pair;
        gs[k].sqrt_pref_real = b->layers[k].sqrt_pref_real;
        gs[k].inv_two_km2    = b->layers[k].inv_two_km2;
        gs[k].alpha_dt       = b->layers[k].alpha_dt;
        gs[k].vram_offset    = b->layers[k].vram_offset;
        gs[k].active         = b->layers[k].active;
    }
    b->gpu = atmturb_cuda_boiling_init(b->nlayers, b->msize, b->ar_order, total_vram, gs);
    free(gs);
}

/* Scratch aligned output buffer stored alongside c2r_work for FFTW C2R execution */
static float *s_boil_c2r_out = NULL;
static long   s_boil_c2r_sz  = 0;

static float *atmturb_boil_get_c2r_out(
    long msize)
{
    long ntot = msize * msize;
    if (s_boil_c2r_out == NULL || s_boil_c2r_sz < ntot)
    {
        if (s_boil_c2r_out != NULL)
        {
            fftwf_free(s_boil_c2r_out);
        }
        s_boil_c2r_out = (float *) fftwf_alloc_real((size_t) ntot);
        s_boil_c2r_sz  = (s_boil_c2r_out != NULL) ? ntot : 0;
    }
    return s_boil_c2r_out;
}

__attribute__((destructor)) static void atmturb_boil_cleanup_static(void)
{
    if (s_boil_c2r_out != NULL)
    {
        fftwf_free(s_boil_c2r_out);
        s_boil_c2r_out = NULL;
        s_boil_c2r_sz  = 0;
    }
}

/**
 * atmturb_boiling_init - Initialize multi-layer co-moving Fourier AR(1)/AR(2) context
 */
int atmturb_boiling_init(
    atmturb_rolling_t       *r,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    long                     msize,
    long                     nbframes,
    double                   time_step_s,
    long                     precision,
    uint64_t                 seed)
{
    int abs_mode = (geom->boil_mode < 0) ? -geom->boil_mode : geom->boil_mode;
    atmturb_boiling_t *b = (atmturb_boiling_t *) calloc(1, sizeof(*b));
    if (b == NULL)
    {
        return -1;
    }
    r->boil     = b;
    b->ar_order = (abs_mode == 2) ? 2 : 1;
    b->msize    = msize;
    b->hcols    = msize / 2 + 1;
    b->nlayers  = prof->nlayers;
    b->seed     = seed;
    b->layers   = (atmturb_boil_layer_t *) calloc((size_t) b->nlayers, sizeof(*b->layers));
    if (b->layers == NULL)
    {
        return -1;
    }

    int total_vram = 0;
    for (int k = 0; k < b->nlayers; k++)
    {
        b->layers[k].vram_offset = total_vram;
        atmturb_boil_config_layer(b, &r->layers[k], &b->layers[k], &prof->layers[k],
                                  geom, k, nbframes, time_step_s);
        total_vram += r->layers[k].nscreens;
    }

    const float **vram_tags = (const float **) calloc((size_t) (total_vram + 1),
                                                      sizeof(float *));
    if (vram_tags == NULL ||
        atmturb_boil_alloc_screens_and_lf(r, prof, geom, time_step_s, precision, vram_tags) != 0)
    {
        free(vram_tags);
        return -1;
    }

    atmturb_boil_try_gpu_backend(b, total_vram, precision);
    if (b->gpu == NULL)
    {
        float *c2r_buf = atmturb_boil_get_c2r_out(msize);
        if (atmturb_boil_alloc_cpu_backend(b, c2r_buf) != 0)
        {
            free(vram_tags);
            return -1;
        }
    }

    for (int k = 0; k < b->nlayers && b->gpu != NULL; k++)
    {
        if (!b->layers[k].active)
        {
            atmturb_cuda_boiling_upload_static(b->gpu, b->layers[k].vram_offset,
                                               r->layers[k].screens[0].data);
        }
    }

    long init_target = (CONF_STREAM_MODE > 0) ? 1 : (b->windowed ? 0 : (nbframes - 1));
    if (init_target < 0)
    {
        init_target = 0;
    }
    int rc = atmturb_boiling_advance(r, init_target, 0);
    if (rc == 0 && b->gpu != NULL)
    {
        rc = atmturb_cuda_boiling_commit_cache(b->gpu, vram_tags);
    }
    free(vram_tags);
    return rc;
}

/**
 * atmturb_boiling_free - Release all FFTW, CUDA, and ring resources for Fourier boiling
 */
void atmturb_boiling_free(
    atmturb_boiling_t *b)
{
    if (b == NULL)
    {
        return;
    }
    atmturb_cuda_boiling_free(b->gpu);
    if (b->c2r_plan != NULL)
    {
        fftwf_destroy_plan(b->c2r_plan);
    }
    if (b->c2r_work != NULL)
    {
        fftwf_free(b->c2r_work);
    }
    if (b->layers != NULL)
    {
        for (int k = 0; k < b->nlayers; k++)
        {
            fftwf_free(b->layers[k].spec_cur);
            fftwf_free(b->layers[k].spec_aux);
            free(b->layers[k].lf_ring_re);
            free(b->layers[k].lf_ring_im);
        }
        free(b->layers);
    }
    free(b);
}

/**
 * atmturb_boiling_advance - Step AR keyframes and subharmonics up to frame target_t
 */
int atmturb_boiling_advance(
    atmturb_rolling_t *r,
    long               target_t,
    int                gpu_only)
{
    if (r == NULL || r->boil == NULL || target_t < 0)
    {
        return 0;
    }
    atmturb_boiling_t *b       = r->boil;
    float             *c2r_out = (b->gpu == NULL) ? atmturb_boil_get_c2r_out(b->msize) : NULL;

    for (int k = 0; k < b->nlayers; k++)
    {
        atmturb_boil_layer_t *bl = &b->layers[k];
        if (!bl->active)
        {
            continue;
        }

        if (r->lowfreq)
        {
            for (long f = bl->max_lf_ready + 1; f <= target_t; f++)
            {
                uint64_t tag = 5000000ULL + (uint64_t) k * 1000003ULL + (uint64_t) f;
                uint64_t ls  = atmturb_rng_stream_seed(b->seed, tag);
                atmturb_boil_step_lowfreq_layer(bl, b->ar_order, f, ls);
            }
        }

        long e      = target_t / (long) bl->k_sub;
        long rem    = target_t % (long) bl->k_sub;
        long req_kf = (rem > 0) ? (e + 1) : e;
        for (long kf = bl->max_kf_ready + 1; kf <= req_kf; kf++)
        {
            if (atmturb_boil_step_keyframe(r, k, kf, gpu_only, c2r_out) != 0)
            {
                return -1;
            }
        }
    }
    return 0;
}

/**
 * atmturb_boiling_next_window - Advance and return next conflict-free 3D cube frame window
 */
int atmturb_boiling_next_window(
    atmturb_rolling_t *r,
    long               t0,
    long               nbframes,
    long              *out_t1)
{
    if (t0 >= nbframes)
    {
        *out_t1 = nbframes;
        return 0;
    }
    if (r == NULL || r->boil == NULL || !r->boil->windowed)
    {
        *out_t1 = nbframes;
        return 1;
    }

    atmturb_boiling_t *b  = r->boil;
    long               t1 = nbframes;
    for (int k = 0; k < b->nlayers; k++)
    {
        const atmturb_boil_layer_t *bl = &b->layers[k];
        if (!bl->active)
        {
            continue;
        }
        long ks     = (long) bl->k_sub;
        long ns     = (long) r->layers[k].nscreens;
        long t1_lim = (t0 / ks + ns - 1) * ks + 1;
        if (t1_lim < t1)
        {
            t1 = t1_lim;
        }
    }
    if (t1 <= t0)
    {
        t1 = t0 + 1;
    }
    (void) atmturb_boiling_advance(r, t1 - 1, 0);
    *out_t1 = t1;
    return 1;
}

/**
 * atmturb_boiling_get_frame - Evaluate correlated keyframe weights and frame-rate subharmonics
 */
void atmturb_boiling_get_frame(
    const atmturb_rolling_t *r,
    int                      layer_idx,
    long                     t,
    atmturb_rolling_eval_t  *out)
{
    const atmturb_rolling_layer_t *rl = &r->layers[layer_idx];
    const atmturb_boil_layer_t    *bl = &r->boil->layers[layer_idx];

    if (!bl->active || rl->nscreens <= 1)
    {
        out->scrA = rl->screens[0].data;
        out->scrB = NULL;
        out->wA   = 1.0f;
        out->wB   = 0.0f;
        out->idxA = 0;
        out->idxB = 0;
        if (r->lowfreq)
        {
            const float *src_re = bl->active
                                      ? &bl->lf_ring_re[(t % bl->lf_cap) * ATMTURB_LOWFREQ_NMODES]
                                      : rl->screens[0].are;
            const float *src_im = bl->active
                                      ? &bl->lf_ring_im[(t % bl->lf_cap) * ATMTURB_LOWFREQ_NMODES]
                                      : rl->screens[0].aim;
            memcpy(out->are_eff, src_re, sizeof(out->are_eff));
            memcpy(out->aim_eff, src_im, sizeof(out->aim_eff));
        }
        return;
    }

    long e    = t / (long) bl->k_sub;
    long rem  = t % (long) bl->k_sub;
    int  idxA = (int) (e % (long) rl->nscreens);
    out->scrA = rl->screens[idxA].data;
    out->idxA = idxA;

    if (rem == 0)
    {
        out->scrB = NULL;
        out->wA   = 1.0f;
        out->wB   = 0.0f;
        out->idxB = idxA;
    }
    else
    {
        int    idxB = (int) ((e + 1) % (long) rl->nscreens);
        double u    = (double) rem / (double) bl->k_sub;
        double nu   = sqrt(1.0 - 2.0 * u * (1.0 - u) * (1.0 - bl->rho_bar));
        out->scrB   = rl->screens[idxB].data;
        out->idxB   = idxB;
        out->wA     = (float) ((1.0 - u) / nu);
        out->wB     = (float) (u / nu);
    }

    if (r->lowfreq)
    {
        long lf_off = (t % bl->lf_cap) * ATMTURB_LOWFREQ_NMODES;
        memcpy(out->are_eff, &bl->lf_ring_re[lf_off], sizeof(out->are_eff));
        memcpy(out->aim_eff, &bl->lf_ring_im[lf_off], sizeof(out->aim_eff));
    }
}
