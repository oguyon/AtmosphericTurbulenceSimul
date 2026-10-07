// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    AtmosphericTurbulence.h
 * @brief   Atmospheric turbulence simulation public interface
 */

#ifndef ATMOSPHERETURBULENCE_H
#define ATMOSPHERETURBULENCE_H

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

#ifdef __cplusplus
}
#endif

#endif // ATMOSPHERETURBULENCE_H
