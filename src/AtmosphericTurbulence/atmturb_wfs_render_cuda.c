// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_wfs_render_cuda.c
 * @brief   CUDA GPU dispatch wrappers for wavefront series rendering
 */

#ifdef HAVE_CUDA

#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "CommandLineInterface/CLIcore.h"
#include "atmturb_cuda.h"
#include "atmturb_cuda_rytov.h"
#include "atmturb_geometry.h"
#include "atmturb_profile.h"
#include "atmturb_rolling.h"
#include "atmturb_rytov.h"
#include "atmturb_types.h"
#include "atmturb_wfs_render.h"
#include "atmturb_wfs_render_cuda.h"

int atmturb_wfs_render_cuda(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    long                        master_size,
    long                        pup_size,
    long                        nbframes,
    const atmturb_wfs_images_t *imgs)
{
    int nl = geom->nlayers;
    size_t m_sz = (size_t) nl * sizeof(float *);
    size_t d_sz = (size_t) nl * sizeof(double);
    void *buf = malloc(m_sz + 5 * d_sz);
    if (!buf)
    {
        return -1;
    }
    const float **masters = (const float **) buf;
    double *vx  = (double *) ((char *) buf + m_sz);
    double *vy  = vx + nl;
    double *cn2 = vy + nl;
    double *x0  = cn2 + nl;
    double *y0  = x0 + nl;

    for (int k = 0; k < nl; k++)
    {
        masters[k] = r->layers[k].screens[0].data;
        vx[k]      = geom->layers[k].vx_pix;
        vy[k]      = geom->layers[k].vy_pix;
        cn2[k]     = geom->layers[k].weight * geom->layers[k].weight;
        x0[k]      = geom->layers[k].x0;
        y0[k]      = geom->layers[k].y0;
    }

    double w0 = geom->layers[0].weight;
    double scoeff = (w0 > 0.0) ? (geom->layers[0].weight_s / w0) : 1.0;
    atmturb_cuda_sim_params_t cparams = {
        .nblayers  = nl,
        .msize     = master_size,
        .pup_size  = pup_size,
        .nbframes  = nbframes,
        .Scoeff    = scoeff,
        .h_masters = (const float *const *) masters,
        .vxpix     = vx,
        .vypix     = vy,
        .cn2       = cn2,
        .x0        = x0,
        .y0        = y0
    };
    long total_pixels = nbframes * pup_size * pup_size;
    float *amp_ptr  = &dcimg[imgs->ID_amp].array.F[0];
    float *samp_ptr = &dcimg[imgs->ID_samp].array.F[0];

    #pragma omp parallel for
    for (long i = 0; i < total_pixels; i++)
    {
        amp_ptr[i]  = 1.0f;
        samp_ptr[i] = 1.0f;
    }

    atmturb_cuda_sim_outputs_t outputs = {
        .pha  = &dcimg[imgs->ID_pha].array.F[0],
        .amp  = NULL,
        .spha = &dcimg[imgs->ID_spha].array.F[0],
        .samp = NULL
    };

    int ret = atmturb_wfs_render_frames_cuda(&cparams, &outputs);
    free(buf);
    return ret;
}

static void atmturb_wfs_cuda_build_sublayers(
    const atmturb_rolling_t       *r,
    const atmturb_geom_t          *geom,
    const atmturb_rytov_plan_t    *plan,
    long                           nbframes,
    double                         time_step_s,
    int                            max_sublayers,
    atmturb_cuda_rytov_sublayer_t *sublayers)
{
    double offset_os = (double) plan->guard_pix * (double) geom->oversample;

    #pragma omp parallel for collapse(2)
    for (long t = 0; t < nbframes; t++)
    {
        for (int m = 0; m < plan->nsuper; m++)
        {
            const atmturb_superlayer_t *sl = &plan->supers[m];
            size_t base = (size_t) ((t * plan->nsuper + m) * max_sublayers);

            for (int j = 0; j < sl->nlayers; j++)
            {
                int k = sl->layer_indices[j];
                float ws = (sl->layer_weights != NULL) ? sl->layer_weights[j] : 1.0f;
                atmturb_rolling_eval_t rev;
                atmturb_rolling_get_frame(r, k, t, time_step_s, &rev);

                double dx = (geom->layers[k].traj_x != NULL)
                            ? geom->layers[k].traj_x[t]
                            : ((double) t * geom->layers[k].vx_pix);
                double dy = (geom->layers[k].traj_y != NULL)
                            ? geom->layers[k].traj_y[t]
                            : ((double) t * geom->layers[k].vy_pix);

                sublayers[base + j].k     = k;
                sublayers[base + j].w_pri = (float) (geom->layers[k].weight * (double) ws *
                                                     (double) rev.wA);
                sublayers[base + j].w_sec = (float) (geom->layers[k].weight_s * (double) ws *
                                                     (double) rev.wA);
                sublayers[base + j].x     = (float) (geom->layers[k].x0 + dx - offset_os);
                sublayers[base + j].y     = (float) (geom->layers[k].y0 + dy - offset_os);
                sublayers[base + j].xs    = (float) (geom->layers[k].xs0 + dx - offset_os);
                sublayers[base + j].ys    = (float) (geom->layers[k].ys0 + dy - offset_os);
            }
        }
    }
}

int atmturb_wfs_render_rytov_cuda(
    const atmturb_rolling_t    *r,
    const atmturb_geom_t       *geom,
    const atmturb_profile_t    *prof,
    const atmturb_obs_params_t *params,
    long                        master_size,
    long                        pup_size,
    long                        nbframes,
    const atmturb_wfs_images_t *imgs)
{
    if (r->lowfreq != 0 || !atmturb_cuda_device_available())
    {
        return -1;
    }

    atmturb_rytov_plan_t plan;
    double z_bin = (double) CONF_FRESNEL_PROPAGATION_BIN;
    if (atmturb_rytov_plan_init(&plan, prof, geom, pup_size, params->pupil_scale_m,
                                params->lambda_ref_m, params->lambda_s_m, z_bin,
                                CONF_FRESNEL_RYTOV_SEC_EXACT) != 0)
    {
        return -1;
    }

    int max_sublayers = 1;
    for (int m = 0; m < plan.nsuper; m++)
    {
        if (plan.supers[m].nlayers > max_sublayers)
        {
            max_sublayers = plan.supers[m].nlayers;
        }
    }
    if (max_sublayers > 64)
    {
        atmturb_rytov_plan_free(&plan);
        return -1;
    }

    size_t n_sub_tot = (size_t) (nbframes * plan.nsuper * max_sublayers);
    atmturb_cuda_rytov_sublayer_t *sublayers =
        (atmturb_cuda_rytov_sublayer_t *) calloc(n_sub_tot, sizeof(atmturb_cuda_rytov_sublayer_t));
    int *super_nlayers   = (int *) malloc(sizeof(int) * (size_t) plan.nsuper);
    double *super_dist_m = (double *) malloc(sizeof(double) * (size_t) plan.nsuper);
    const float **masters = (const float **) malloc(sizeof(float *) * (size_t) geom->nlayers);

    if (!sublayers || !super_nlayers || !super_dist_m || !masters)
    {
        free(sublayers);
        free(super_nlayers);
        free(super_dist_m);
        free(masters);
        atmturb_rytov_plan_free(&plan);
        return -1;
    }

    for (int k = 0; k < geom->nlayers; k++)
    {
        masters[k] = r->layers[k].screens[0].data;
    }
    for (int m = 0; m < plan.nsuper; m++)
    {
        super_nlayers[m] = plan.supers[m].nlayers;
        super_dist_m[m]  = plan.supers[m].dist_m;
    }

    atmturb_wfs_cuda_build_sublayers(r, geom, &plan, nbframes, params->time_step_s,
                                     max_sublayers, sublayers);

    atmturb_cuda_rytov_params_t cparams = {
        .msize           = master_size,
        .pup_size        = pup_size,
        .pad_size        = plan.pad_size,
        .guard_pix       = plan.guard_pix,
        .nbframes        = nbframes,
        .nblayers        = geom->nlayers,
        .nsuper          = plan.nsuper,
        .has_sec         = (plan.lambda_s_m > 0.0) ? 1 : 0,
        .sec_shared      = plan.sec_shared,
        .os              = (int) geom->oversample,
        .interp          = geom->interp,
        .h_masters       = masters,
        .h_laplace_inv   = plan.laplace_inv,
        .h_exp_x         = (const void *) plan.exp_x,
        .h_exp_y         = (const void *) plan.exp_y,
        .h_filt_a_pri    = (const float *const *) plan.filter_a_pri,
        .h_filt_b_pri    = (const float *const *) plan.filter_b_pri,
        .h_filt_a_sec    = (const float *const *) plan.filter_a_sec,
        .h_filt_b_sec    = (const float *const *) plan.filter_b_sec,
        .h_chrom_ramp    = (const void *const *) plan.chrom_ramp,
        .frame_sublayers = sublayers,
        .super_nlayers   = super_nlayers,
        .super_dist_m    = super_dist_m,
        .max_sublayers   = max_sublayers,
        .pha             = &dcimg[imgs->ID_pha].array.F[0],
        .amp             = &dcimg[imgs->ID_amp].array.F[0],
        .spha            = (plan.lambda_s_m > 0.0) ? &dcimg[imgs->ID_spha].array.F[0] : NULL,
        .samp            = (plan.lambda_s_m > 0.0) ? &dcimg[imgs->ID_samp].array.F[0] : NULL
    };

    int ret = atmturb_cuda_rytov_render(&cparams);

    free(sublayers);
    free(super_nlayers);
    free(super_dist_m);
    free(masters);
    atmturb_rytov_plan_free(&plan);

    return ret;
}

#endif // HAVE_CUDA
