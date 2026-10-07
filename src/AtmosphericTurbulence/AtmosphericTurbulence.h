// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    AtmosphericTurbulence.h
 * @brief   Atmospheric turbulence simulation public interface
 */

#ifndef ATMOSPHERETURBULENCE_H
#define ATMOSPHERETURBULENCE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * init_AtmosphericTurbulence - Initialize AtmosphericTurbulence module state
 *
 * Return: 0 on success.
 */
int init_AtmosphericTurbulence(void);

/**
 * AtmosphericTurbulence_change_configuration_file - Set active configuration file path
 * @fname: Path to configuration file.
 *
 * Return: 0 on success.
 */
int AtmosphericTurbulence_change_configuration_file(
    const char *fname);

/**
 * make_AtmosphericTurbulence_vonKarmanWind - Generate 3D velocity cube [u, v, w]
 * @vKsize: Sample length of 1D series.
 * @pixscale: Physical sampling step in meters.
 * @sigmawind: Velocity standard deviation in m/s.
 * @Lwind: Wind velocity turbulence outer scale in meters.
 * @seed: RNG seed (0 = time-based); u, v, w use independent streams.
 * @IDout_name: Output 3D image name (vKsize x 1 x 3).
 *
 * Return: Output image ID on success.
 */
long make_AtmosphericTurbulence_vonKarmanWind(
    long        vKsize,
    float       pixscale,
    float       sigmawind,
    float       Lwind,
    long        seed,
    const char *IDout_name);

/**
 * struct atmturb_wind_traj_params_t - Trajectory generation parameters
 * @nbframes: Number of frames in simulation series.
 * @dt_s: Time step per frame in seconds.
 * @vx_pix: Mean pupil-plane x velocity [master pixels / frame].
 * @vy_pix: Mean pupil-plane y velocity [master pixels / frame].
 * @dx_master_m: Physical master pixel scale in meters.
 * @sigma_wind_mps: Wind velocity standard deviation [m/s].
 * @L_wind_m: Wind turbulence outer scale [m].
 * @seed: Resolved PRNG seed.
 */
typedef struct
{
    long     nbframes;
    double   dt_s;
    double   vx_pix;
    double   vy_pix;
    double   dx_master_m;
    double   sigma_wind_mps;
    double   L_wind_m;
    uint64_t seed;
} atmturb_wind_traj_params_t;

/**
 * atmturb_wind_synthesize_trajectory - Synthesize cumulative 2D trajectory with turbulent wind
 * @params: Trajectory generation input parameters.
 * @traj_x: Output buffer of size nbframes for cumulative x offset [pixels].
 * @traj_y: Output buffer of size nbframes for cumulative y offset [pixels].
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_wind_synthesize_trajectory(
    const atmturb_wind_traj_params_t *params,
    double                           *traj_x,
    double                           *traj_y);


/**
 * struct atmturb_screen_spec_t - Specification for phase screen generation
 * @size: Grid dimension in pixels (N).
 * @r0_pix: Fried parameter in pixel units (>0, or <= 0 for default 6.88^0.6).
 * @L0_pix: Outer scale in pixel units (<= 0 for infinite).
 * @l0_pix: Inner scale in pixel units (<= 0 for none).
 * @seed: PRNG seed value (0 selects time-based seed).
 * @precision: Internal FFT precision (0 = single, 1 = double).
 */
typedef struct
{
    long     size;
    double   r0_pix;
    double   L0_pix;
    double   l0_pix;
    uint64_t seed;
    int      precision;
} atmturb_screen_spec_t;

/**
 * atmturb_generate_screen_pair - Generate two independent phase screens from spectrum
 * @spec: Generation parameters and scales.
 * @screen_a: Output buffer of size * size floats (or NULL to skip).
 * @screen_b: Output buffer of size * size floats (or NULL to skip).
 *
 * Return: 0 on success, non-zero error code otherwise.
 */
int atmturb_generate_screen_pair(
    const atmturb_screen_spec_t *spec,
    float                       *screen_a,
    float                       *screen_b);

/**
 * make_master_turbulence_screen_seeded - Generate von Karman master screens with seed
 * @ID_name1: Output name for screen 1.
 * @ID_name2: Output name for screen 2.
 * @size: Grid dimension in pixels.
 * @outerscale: Outer scale in pixels.
 * @innerscale: Inner scale in pixels.
 * @WFprecision: Precision flag (0=single, 1=double).
 * @seed: PRNG seed value (0 = time-based).
 *
 * Return: 0 on success.
 */
int make_master_turbulence_screen_seeded(
    const char *ID_name1,
    const char *ID_name2,
    long        size,
    float       outerscale,
    float       innerscale,
    long        WFprecision,
    uint64_t    seed);

/**
 * atmturb_measure_r0_pix - Estimate effective r0 in pixels from structure function at lag 1
 * @data: Pointer to 2D float screen array.
 * @size: Grid dimension in pixels.
 *
 * Return: Estimated r0 in pixel units.
 */
double atmturb_measure_r0_pix(
    const float *data,
    long         size);

/**
 * make_master_turbulence_screen - Generate von Karman master turbulence screens
 * @ID_name1: Output name for screen 1.
 * @ID_name2: Output name for screen 2.
 * @size: Grid dimension in pixels.
 * @outerscale: Outer scale in pixels.
 * @innerscale: Inner scale in pixels.
 * @WFprecision: Precision flag (0=single, 1=double).
 *
 * Return: 0 on success.
 */
int make_master_turbulence_screen(
    const char *ID_name1,
    const char *ID_name2,
    long        size,
    float       outerscale,
    float       innerscale,
    long        WFprecision);

/**
 * make_master_turbulence_screen_pow - Generate phase screen with arbitrary power-law PSD
 * @ID_name1: Output name for screen 1.
 * @ID_name2: Output name for screen 2.
 * @size: Grid dimension in pixels.
 * @power: Spatial frequency power-law exponent.
 *
 * Return: 0 on success.
 */
int make_master_turbulence_screen_pow(
    const char *ID_name1,
    const char *ID_name2,
    long        size,
    float       power);

/**
 * contract_wavefront_series - Decimate series of complex wavefront cubes by factor of 2
 * @in_prefix: Input filename prefix.
 * @out_prefix: Output filename prefix.
 * @NB_files: Number of files in sequence.
 *
 * Return: 0 on success.
 */
int contract_wavefront_series(
    const char *in_prefix,
    const char *out_prefix,
    long        NB_files);

/**
 * contract_wavefront_cube - Decimate complex wavefront cube by 2^factor
 * @ina_file: Input amplitude FITS file.
 * @inp_file: Input phase FITS file.
 * @outa_file: Output amplitude FITS file.
 * @outp_file: Output phase FITS file.
 * @factor: Power of 2 decimation exponent.
 *
 * Return: 0 on success.
 */
int contract_wavefront_cube(
    const char *ina_file,
    const char *inp_file,
    const char *outa_file,
    const char *outp_file,
    int         factor);

/**
 * contract_wavefront_cube_phaseonly - Decimate phase-only wavefront cube by 2^factor
 * @inp_file: Input phase FITS file.
 * @outp_file: Output phase FITS file.
 * @factor: Power of 2 decimation exponent.
 *
 * Return: 0 on success.
 */
int contract_wavefront_cube_phaseonly(
    const char *inp_file,
    const char *outp_file,
    int         factor);

/**
 * make_AtmosphericTurbulence_wavefront_series - Run full wavefront simulation series
 * @slambdaum: Secondary observing wavelength in um.
 * @WFprecision: Precision mode flag.
 *
 * Return: 0 on success, -1 on failure.
 */
int make_AtmosphericTurbulence_wavefront_series(
    float slambdaum,
    long  WFprecision);

/**
 * measure_wavefront_series - Process wavefront series and extract PSF metrics
 * @factor: Decimation/binning factor.
 *
 * Return: 0 on success.
 */
int measure_wavefront_series(
    float factor);

/**
 * AtmosphericTurbulence_mkTestTTseq - Generate synthetic tip-tilt test sequence
 * @dt: Time step in seconds.
 * @NBpts: Number of points per block.
 * @NBblocks: Number of blocks.
 * @measnoise: Measurement noise RMS.
 * @ACCnmode: Accelerometer mode flag.
 * @ACCnoise: Accelerometer noise RMS.
 * @MODE: Generation mode flag.
 *
 * Return: 0 on success.
 */
int AtmosphericTurbulence_mkTestTTseq(
    double dt,
    long   NBpts,
    long   NBblocks,
    double measnoise,
    int    ACCnmode,
    double ACCnoise,
    int    MODE);

/**
 * AtmosphericTurbulence_Build_LinPredictor_Full - Build full-aperture AR linear prediction matrix
 * @WFin_name: Input wavefront cube name.
 * @WFmask_name: Pupil mask image name.
 * @PForder: Autoregressive predictor order (history steps).
 * @PFlag: Prediction horizon lag.
 * @SVDeps: Singular value cutoff tolerance.
 * @RegLambda: Tikhonov regularization parameter.
 *
 * Return: 0 on success, -1 on failure.
 */
int AtmosphericTurbulence_Build_LinPredictor_Full(
    const char *WFin_name,
    const char *WFmask_name,
    int         PForder,
    float       PFlag,
    double      SVDeps,
    double      RegLambda);

/**
 * AtmosphericTurbulence_Apply_LinPredictor_Full - Apply full-aperture linear predictor
 * @MODE: Application mode flag.
 * @WFin_name: Input wavefront cube name.
 * @WFmask_name: Pupil mask image name.
 * @PForder: Filter AR order.
 * @PFlag: Prediction horizon lag.
 * @WFoutp_name: Output predicted wavefront cube name.
 * @WFoutf_name: Output residual error cube name.
 *
 * Return: 0 on success.
 */
int AtmosphericTurbulence_Apply_LinPredictor_Full(
    int         MODE,
    const char *WFin_name,
    const char *WFmask_name,
    int         PForder,
    float       PFlag,
    const char *WFoutp_name,
    const char *WFoutf_name);

/**
 * AtmosphericTurbulence_LinPredictor_filt_2DKernelExtract - Extract shift-invariant 2D kernels
 * @IDfilt_name: Input filter matrix image name.
 * @IDmask_name: Active pupil mask image name.
 * @krad: Extraction neighborhood radius in pixels.
 * @IDkern_name: Output 3D kernel image name.
 *
 * Return: Output image ID on success.
 */
long AtmosphericTurbulence_LinPredictor_filt_2DKernelExtract(
    const char *IDfilt_name,
    const char *IDmask_name,
    long        krad,
    const char *IDkern_name);

/**
 * AtmosphericTurbulence_LinPredictor_filt_Expand - Expand 2D kernel across full aperture
 * @IDfilt_name: Input 2D/3D shift-invariant kernel.
 * @IDmask_name: Active pupil mask.
 *
 * Return: Output expanded matrix image ID.
 */
long AtmosphericTurbulence_LinPredictor_filt_Expand(
    const char *IDfilt_name,
    const char *IDmask_name);

/**
 * AtmosphericTurbulence_Build_LinPredictor - Train linear predictor at single reference pixel
 * @NB_WFstep: Number of input wavefront frames.
 * @WFphaNoise: Measurement phase noise level.
 * @WFPlag: Prediction lag in frames.
 * @WFP_NBstep: History steps (AR order).
 * @WFP_xyrad: Spatial footprint radius around target pixel.
 * @WFPiipix: Target pixel X coordinate.
 * @WFPjjpix: Target pixel Y coordinate.
 * @slambdaum: Observing wavelength in um.
 *
 * Return: 0 on success.
 */
int AtmosphericTurbulence_Build_LinPredictor(
    long   NB_WFstep,
    double WFphaNoise,
    long   WFPlag,
    long   WFP_NBstep,
    long   WFP_xyrad,
    long   WFPiipix,
    long   WFPjjpix,
    float  slambdaum);

/**
 * AtmosphericTurbulence_psfCubeContrast - Measure dark-hole contrast profile from PSF cube
 * @IDwfc_name: Wavefront / PSF cube name.
 * @IDmask_name: Mask image name (1 inside dark hole, 0 elsewhere).
 * @IDpsfc_name: Output contrast profile image name.
 *
 * Return: Number of processed slices.
 */
long AtmosphericTurbulence_psfCubeContrast(
    const char *IDwfc_name,
    const char *IDmask_name,
    const char *IDpsfc_name);

/**
 * AtmosphericTurbulence_Test_LinPredictor - Verify linear predictor on wavefront series
 * @NB_WFstep: Number of test steps.
 * @WFphaNoise: Measurement phase noise level.
 * @IDWFPfilt_name: Filter kernel image name.
 * @WFPlag: Prediction lag.
 * @WFPiipix: Target pixel X.
 * @WFPjjpix: Target pixel Y.
 * @slambdaum: Wavelength in um.
 *
 * Return: 0 on success.
 */
int AtmosphericTurbulence_Test_LinPredictor(
    long        NB_WFstep,
    double      WFphaNoise,
    const char *IDWFPfilt_name,
    long        WFPlag,
    long        WFPiipix,
    long        WFPjjpix,
    float       slambdaum);

/**
 * measure_wavefront_series_expoframes - Integrate frames into specified exposure duration
 * @etime: Exposure duration in seconds.
 * @outfile: Destination output text file.
 *
 * Return: 0 on success.
 */
int measure_wavefront_series_expoframes(
    float       etime,
    const char *outfile);

/**
 * frame_select_PSF - Lucky imaging frame selection from logged PSF series
 * @logfile: Text log containing frame index and metric columns.
 * @NBfiles: Number of log rows / frames.
 * @frac: Selection fraction (e.g. 0.1 for top 10%).
 *
 * Return: 0 on success.
 */
int frame_select_PSF(
    const char *logfile,
    long        NBfiles,
    float       frac);

/**
 * AtmosphericTurbulence_WFprocess - Standalone test routine for automated AO grid evaluation
 *
 * Return: 0 on success.
 */
int AtmosphericTurbulence_WFprocess(void);

/**
 * struct atmturb_ao_params_t - Configuration for AO closed-loop simulation
 * @in_wfname: Input phase cube name or FITS file path (default: "outarraypha").
 * @in_ampname: Optional input amplitude cube name or FITS file path (NULL = 1.0).
 * @out_psfname: Output PSF cube or cumulative image name (default: "PSFcumul").
 * @out_fitsname: Optional FITS file to write cumulative PSF (NULL = do not save).
 * @pup_size: Linear dimension of pupil array [pixels] (0 = auto from cube).
 * @nbframes: Number of frames to process (0 = all frames in cube).
 * @lambda_ref_m: Reference wavelength of phase cube in meters (default: 0.55e-6).
 * @lambda_sci_m: Science wavelength in meters (default: 1.65e-6).
 * @tel_diam_m: Telescope primary mirror diameter in meters (default: 8.0).
 * @pupil_scale_m: Pupil pixel scale in meters/pixel (<= 0 for auto).
 * @loop_mode: 0 = Open loop, 1 = Leaky integrator, 2 = PID.
 * @gain: Loop gain (default: 0.5).
 * @leak: Integrator leak factor (default: 0.001).
 * @loop_delay: Hardware latency delay in frames (default: 1).
 * @Kp: Proportional gain for PID (default: 0.5).
 * @Ki: Integral gain for PID (default: 0.0).
 * @Kd: Derivative gain for PID (default: 0.0).
 * @save_psfcube: 1 to save full 3D PSF cube, 0 for cumulative 2D PSF only.
 */
typedef struct
{
    const char *in_wfname;
    const char *in_ampname;
    const char *out_psfname;
    const char *out_fitsname;
    long        pup_size;
    long        nbframes;
    double      lambda_ref_m;
    double      lambda_sci_m;
    double      tel_diam_m;
    double      pupil_scale_m;
    int         loop_mode;
    double      gain;
    double      leak;
    int         loop_delay;
    double      Kp;
    double      Ki;
    double      Kd;
    int         save_psfcube;
} atmturb_ao_params_t;

/**
 * struct atmturb_ao_results_t - Simulation summary metrics
 * @strehl_cumul: Long-exposure cumulative Strehl ratio.
 * @strehl_first: First frame Strehl ratio.
 * @strehl_last: Last frame Strehl ratio.
 * @dl_peak: Peak intensity of diffraction-limited reference PSF.
 * @open_loop_wfe_rms: Open-loop RMS wavefront error in meters.
 * @closed_loop_wfe_rms: Closed-loop RMS residual wavefront error in meters.
 */
typedef struct
{
    double strehl_cumul;
    double strehl_first;
    double strehl_last;
    double dl_peak;
    double open_loop_wfe_rms;
    double closed_loop_wfe_rms;
} atmturb_ao_results_t;

/**
 * atmturb_ao_init_params - Populate default parameters for AO simulation
 * @params: Structure to initialize.
 */
void atmturb_ao_init_params(
    atmturb_ao_params_t *params);

/**
 * atmturb_ao_sim_run - Execute closed-loop AO simulation and generate science PSF
 * @params: Simulation parameters.
 * @results: Output results container (optional, can be NULL).
 *
 * Return: 0 on success, -1 on failure.
 */
int atmturb_ao_sim_run(
    const atmturb_ao_params_t *params,
    atmturb_ao_results_t      *results);

/**
 * AtmosphericTurbulence_makePSF - Run closed-loop AO simulation and generate PSF
 * @Kp: Proportional feedback gain.
 * @Ki: Integral feedback gain.
 * @Kd: Derivative feedback gain.
 * @Kdgain: Dynamic gain multiplier.
 *
 * Return: Peak intensity or Strehl estimate.
 */
double AtmosphericTurbulence_makePSF(
    double Kp,
    double Ki,
    double Kd,
    double Kdgain);

/**
 * AtmosphericTurbulence_makeHV_CN2prof - Generate Hufnagel-Valley Cn2 profile file
 * @wspeed: Upper atmospheric wind speed [m/s].
 * @r0: Target Fried parameter in meters.
 * @sitealt: Observatory elevation in meters.
 * @NBlayer: Number of vertical discrete layers.
 * @outfile: Destination profile file path.
 *
 * Return: 0 on success, -1 on failure.
 */
int AtmosphericTurbulence_makeHV_CN2prof(
    double      wspeed,
    double      r0,
    double      sitealt,
    long        NBlayer,
    const char *outfile);

/**
 * AtmosphericTurbulence_makeHV_CN2prof_opt - Generate HV profile with model options
 * @wspeed: Upper atmospheric wind speed [m/s].
 * @r0: Target Fried parameter in meters.
 * @sitealt: Observatory elevation in meters.
 * @NBlayer: Number of vertical discrete layers.
 * @outfile: Destination profile file path.
 * @wind_model: Wind velocity model (0 = legacy, 1 = Bufton).
 * @seed: Master RNG seed value.
 *
 * Return: 0 on success, -1 on failure.
 */
int AtmosphericTurbulence_makeHV_CN2prof_opt(
    double      wspeed,
    double      r0,
    double      sitealt,
    long        NBlayer,
    const char *outfile,
    int         wind_model,
    uint64_t    seed);

#ifdef __cplusplus
}
#endif

#endif // ATMOSPHERETURBULENCE_H
