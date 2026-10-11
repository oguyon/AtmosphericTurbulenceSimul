// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_wfs_cuda_common.c
 * @brief   Shared helper implementations for CUDA wavefront rendering and streaming
 */

#include <stdlib.h>
#include <string.h>

#include "atmturb_wfs_cuda_common.h"

/**
 * atmturb_wfs_cuda_init_screen_pool - Assemble pointers to master phase screens
 * @r: Rolling simulation context.
 * @nlayers: Number of simulation layers.
 * @out_masters: Pointer to store allocated screen pointer array.
 * @out_offsets: Pointer to store allocated layer screen offsets.
 * @out_total: Pointer to store total screen count.
 *
 * Flattens master phase screen arrays across all turbulent layers into a contiguous pointer list
 * for upload and caching in GPU device memory.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_wfs_cuda_init_screen_pool(
    const atmturb_rolling_t  *r,
    int                       nlayers,
    const float            ***out_masters,
    int                     **out_offsets,
    int                      *out_total)
{
    int total = 0;
    for (int k = 0; k < nlayers; k++)
    {
        total += r->layers[k].nscreens;
    }

    const float **masters = (const float **) malloc(sizeof(float *) * (size_t) total);
    int *offsets = (int *) malloc(sizeof(int) * (size_t) nlayers);
    if (!masters || !offsets)
    {
        free(masters);
        free(offsets);
        return -1;
    }

    int off = 0;
    for (int k = 0; k < nlayers; k++)
    {
        offsets[k] = off;
        for (int s = 0; s < r->layers[k].nscreens; s++)
        {
            masters[off + s] = r->layers[k].screens[s].data;
        }
        off += r->layers[k].nscreens;
    }

    *out_masters = masters;
    *out_offsets = offsets;
    *out_total   = total;
    return 0;
}

/**
 * atmturb_wfs_cuda_build_frame_sublayers - Populate sublayer parameters for a single frame
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @nsuper: Number of superlayers.
 * @supers: Array of superlayer definitions.
 * @screen_offsets: Screen index offsets per layer.
 * @t: Frame index.
 * @time_step_s: Time step in seconds.
 * @offset_os: Offset in oversampled master pixels for guard band alignment.
 * @max_sublayers: Maximum allocated sublayers per superlayer.
 * @sublayers: Output array of sublayer interpolation parameters.
 *
 * Computes wind-displaced pupil coordinates, chromatic secondary shifts, and rolling screen
 * interpolation weights for each layer comprising the active superlayers at frame t.
 */
void atmturb_wfs_cuda_build_frame_sublayers(
    const atmturb_rolling_t       *r,
    const atmturb_geom_t          *geom,
    int                            nsuper,
    const atmturb_superlayer_t    *supers,
    const int                     *screen_offsets,
    long                           t,
    double                         time_step_s,
    double                         offset_os,
    int                            max_sublayers,
    atmturb_cuda_rytov_sublayer_t *sublayers)
{
    for (int m = 0; m < nsuper; m++)
    {
        const atmturb_superlayer_t *sl = &supers[m];
        size_t base = (size_t) (m * max_sublayers);

        for (int j = 0; j < sl->nlayers; j++)
        {
            int k = sl->layer_indices[j];
            float ws = (sl->layer_weights != NULL) ? sl->layer_weights[j] : 1.0f;
            atmturb_rolling_eval_t rev;
            atmturb_rolling_get_frame(r, k, t, time_step_s, &rev);

            double dx = atmturb_geom_get_layer_dx(&geom->layers[k], t);
            double dy = atmturb_geom_get_layer_dy(&geom->layers[k], t);

            float cur_x  = (float) (geom->layers[k].x0 + dx - offset_os);
            float cur_y  = (float) (geom->layers[k].y0 + dy - offset_os);
            float cur_xs = (float) (geom->layers[k].xs0 + dx - offset_os);
            float cur_ys = (float) (geom->layers[k].ys0 + dy - offset_os);

            double w_base   = geom->layers[k].weight * (double) ws;
            double ws_base  = geom->layers[k].weight_s * (double) ws;
            int    k_offset = screen_offsets ? screen_offsets[k] : k;

            size_t idx_a = base + (size_t) (2 * j);
            sublayers[idx_a].k     = k_offset + rev.idxA;
            sublayers[idx_a].w_pri = (float) (w_base * (double) rev.wA);
            sublayers[idx_a].w_sec = (float) (ws_base * (double) rev.wA);
            sublayers[idx_a].x     = cur_x;
            sublayers[idx_a].y     = cur_y;
            sublayers[idx_a].xs    = cur_xs;
            sublayers[idx_a].ys    = cur_ys;

            size_t idx_b = base + (size_t) (2 * j + 1);
            if (rev.wB > 0.0f && rev.idxB >= 0)
            {
                sublayers[idx_b].k     = k_offset + rev.idxB;
                sublayers[idx_b].w_pri = (float) (w_base * (double) rev.wB);
                sublayers[idx_b].w_sec = (float) (ws_base * (double) rev.wB);
                sublayers[idx_b].x     = cur_x;
                sublayers[idx_b].y     = cur_y;
                sublayers[idx_b].xs    = cur_xs;
                sublayers[idx_b].ys    = cur_ys;
            }
            else
            {
                memset(&sublayers[idx_b], 0, sizeof(atmturb_cuda_rytov_sublayer_t));
            }
        }
    }
}

/**
 * atmturb_wfs_cuda_build_frame_lf_modes - Populate subharmonic low-frequency modes for a frame
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @nsuper: Number of superlayers.
 * @supers: Array of superlayer definitions.
 * @t: Frame index.
 * @time_step_s: Time step in seconds.
 * @offset_os: Offset in oversampled master pixels for guard band alignment.
 * @max_lf_modes: Maximum allocated modes per superlayer.
 * @lf_modes: Output array of subharmonic Fourier mode parameters.
 *
 * Precomputes spatial wavevectors, primary/secondary initial phases, and temporal amplitudes
 * for all active low-frequency subharmonics across superlayers.
 */
void atmturb_wfs_cuda_build_frame_lf_modes(
    const atmturb_rolling_t     *r,
    const atmturb_geom_t        *geom,
    int                          nsuper,
    const atmturb_superlayer_t  *supers,
    long                         t,
    double                       time_step_s,
    double                       offset_os,
    int                          max_lf_modes,
    atmturb_cuda_lowfreq_mode_t *lf_modes)
{
    if (r->lowfreq == 0 || lf_modes == NULL)
    {
        return;
    }

    for (int m = 0; m < nsuper; m++)
    {
        const atmturb_superlayer_t *sl = &supers[m];
        size_t base = (size_t) (m * max_lf_modes);
        int mode_idx = 0;

        for (int j = 0; j < sl->nlayers; j++)
        {
            int k = sl->layer_indices[j];
            float ws = (sl->layer_weights != NULL) ? sl->layer_weights[j] : 1.0f;
            atmturb_rolling_eval_t rev;
            atmturb_rolling_get_frame(r, k, t, time_step_s, &rev);

            const atmturb_layer_geom_t *lg = &geom->layers[k];
            double dx = atmturb_geom_get_layer_dx(lg, t);
            double dy = atmturb_geom_get_layer_dy(lg, t);

            double x  = lg->x0 + dx - offset_os;
            double y  = lg->y0 + dy - offset_os;
            double xs = lg->xs0 + dx - offset_os;
            double ys = lg->ys0 + dy - offset_os;

            double w_pri = lg->weight * (double) ws;
            double w_sec = lg->weight_s * (double) ws;
            const atmturb_lowfreq_t *lf_base = &r->layers[k].lf_base;

            for (int lm = 0; lm < lf_base->nmodes; lm++)
            {
                float kx = lf_base->kx[lm];
                float ky = lf_base->ky[lm];
                float are = rev.are_eff[lm];
                float aim = rev.aim_eff[lm];

                atmturb_cuda_lowfreq_mode_t *out = &lf_modes[base + mode_idx++];
                out->are_pri = (float) (w_pri * (double) are);
                out->aim_pri = (float) (w_pri * (double) aim);
                out->are_sec = (float) (w_sec * (double) are);
                out->aim_sec = (float) (w_sec * (double) aim);
                out->th0_pri = (float) ((double) kx * x + (double) ky * y);
                out->th0_sec = (float) ((double) kx * xs + (double) ky * ys);
                out->dkx     = (float) ((double) kx * (double) geom->oversample);
                out->dky     = (float) ((double) ky * (double) geom->oversample);
            }
        }
    }
}
