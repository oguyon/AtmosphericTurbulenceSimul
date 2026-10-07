// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_fresnel.h
 * @brief   Super-layer binning and diffractive propagation plan
 */

#ifndef ATMTURB_FRESNEL_H
#define ATMTURB_FRESNEL_H

#include "atmturb_geometry.h"
#include "atmturb_profile.h"
#include "atmturb_rolling.h"
#include "WFpropagate/wfprop_fresnel_engine.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * struct atmturb_superlayer_t - Grouped atmospheric layers for Fresnel diffraction
 * @nlayers: Count of profile layers binned into this super-layer.
 * @layer_indices: Array of indices in geom->layers / rsim->layers.
 * @dist_m: Line-of-sight centroid distance from telescope pupil [meters].
 * @step_dist_m: Propagation distance to the next super-layer or ground [meters].
 */
typedef struct
{
    int     nlayers;
    int    *layer_indices;
    double  dist_m;
    double  step_dist_m;
} atmturb_superlayer_t;

/**
 * struct atmturb_fresnel_plan_t - Precomputed diffractive propagation plan
 * @nsuper: Number of super-layers.
 * @supers: Array of super-layer definitions.
 * @grid_size: Linear dimension of 2D grid in pixels.
 * @pixscale_m: Grid physical pixel size [m/pixel].
 * @lambda_ref_m: Primary reference wavelength [m].
 * @lambda_s_m: Secondary observing wavelength [m].
 * @tf_pri: Array of nsuper precomputed transfer functions for primary wavelength.
 * @tf_sec: Array of nsuper precomputed transfer functions for secondary wavelength.
 */
typedef struct
{
    int                   nsuper;
    atmturb_superlayer_t *supers;
    long                  grid_size;
    double                pixscale_m;
    double                lambda_ref_m;
    double                lambda_s_m;
    fftwf_complex       **tf_pri;
    fftwf_complex       **tf_sec;
} atmturb_fresnel_plan_t;

/**
 * struct atmturb_fresnel_ctx_t - Thread-local scratchpad and engine for diffractive render
 * @eng: Initialized FFTW Fresnel propagation engine.
 * @field_pri: Working complex array for primary optical field.
 * @field_sec: Working complex array for secondary optical field.
 * @super_pha: Scratchpad for super-layer primary phase screen.
 * @super_spha: Scratchpad for super-layer secondary phase screen.
 */
typedef struct
{
    wfprop_fresnel_engine_t eng;
    fftwf_complex          *field_pri;
    fftwf_complex          *field_sec;
    float                  *super_pha;
    float                  *super_spha;
} atmturb_fresnel_ctx_t;

/**
 * atmturb_fresnel_plan_init - Initialize super-layer binning and transfer functions
 * @plan: Propagation plan container to populate.
 * @prof: Active turbulence profile.
 * @geom: Computed observing geometry.
 * @pup_size: Linear dimension of wavefront in pixels.
 * @pixscale_m: Physical pixel scale [m/pixel].
 * @lambda_ref_m: Primary reference wavelength [m].
 * @lambda_s_m: Secondary observing wavelength [m].
 * @z_bin_m: Altitude binning distance threshold [m].
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_fresnel_plan_init(
    atmturb_fresnel_plan_t  *plan,
    const atmturb_profile_t *prof,
    const atmturb_geom_t    *geom,
    long                     pup_size,
    double                   pixscale_m,
    double                   lambda_ref_m,
    double                   lambda_s_m,
    double                   z_bin_m);

/**
 * atmturb_fresnel_plan_free - Release resources held by Fresnel propagation plan
 * @plan: Propagation plan container to tear down.
 */
void atmturb_fresnel_plan_free(
    atmturb_fresnel_plan_t *plan);

/**
 * atmturb_fresnel_ctx_init - Initialize thread-local Fresnel render context
 * @ctx: Thread context structure to initialize.
 * @pup_size: Linear dimension of wavefront in pixels.
 *
 * Return: 0 on success, -1 on allocation failure.
 */
int atmturb_fresnel_ctx_init(
    atmturb_fresnel_ctx_t *ctx,
    long                   pup_size);

/**
 * atmturb_fresnel_ctx_free - Release resources held by thread context
 * @ctx: Thread context structure to tear down.
 */
void atmturb_fresnel_ctx_free(
    atmturb_fresnel_ctx_t *ctx);

/**
 * atmturb_fresnel_render_step - Render one frame using multi-layer diffractive propagation
 * @ctx: Thread-local execution context.
 * @plan: Precomputed diffractive propagation plan.
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @t: Simulation frame index.
 * @time_step_s: Time step between frames in seconds.
 * @master_size: Master screen dimension in pixels.
 * @pup_size: Linear dimension of wavefront in pixels.
 * @pha_slice: Destination primary phase slice (accumulates unwrapped diffractive phase).
 * @amp_slice: Destination primary amplitude slice.
 * @spha_slice: Destination secondary phase slice (accumulates unwrapped diffractive phase).
 * @samp_slice: Destination secondary amplitude slice.
 */
void atmturb_fresnel_render_step(
    atmturb_fresnel_ctx_t        *ctx,
    const atmturb_fresnel_plan_t *plan,
    const atmturb_rolling_t      *r,
    const atmturb_geom_t         *geom,
    long                          t,
    double                        time_step_s,
    long                          master_size,
    long                          pup_size,
    float                        *pha_slice,
    float                        *amp_slice,
    float                        *spha_slice,
    float                        *samp_slice);

/**
 * atmturb_wfs_render_layer - Render one turbulence layer into phase slices for frame t
 * @r: Rolling simulation context.
 * @geom: Computed observing geometry.
 * @k: Layer index.
 * @t: Frame index.
 * @time_step_s: Time step in seconds.
 * @master_size: Master screen dimension in pixels.
 * @pup_size: Pupil dimension in pixels.
 * @pha_slice: Primary phase frame slice.
 * @spha_slice: Secondary phase frame slice.
 */
void atmturb_wfs_render_layer(
    const atmturb_rolling_t *r,
    const atmturb_geom_t    *geom,
    int                      k,
    long                     t,
    double                   time_step_s,
    long                     master_size,
    long                     pup_size,
    float                   *pha_slice,
    float                   *spha_slice);

#ifdef __cplusplus
}
#endif

#endif // ATMTURB_FRESNEL_H
