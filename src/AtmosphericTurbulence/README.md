# AtmosphericTurbulence

**Architectural Layer**: Level 2 (Turbulence Simulation, Wavefront Generation & Predictive Filtering)

## Purpose

The `AtmosphericTurbulence` module provides multi-layer atmospheric turbulence simulation using the frozen-flow (Taylor) hypothesis and von Karman / Kolmogorov phase statistics. It generates master phase screens, synthesizes 3D turbulent wind velocity fields, computes Hufnagel-Valley $C_n^2$ profiles, extrudes multi-layer phase screens, calculates Fresnel or geometric diffractive wavefront cubes, evaluates point spread function (PSF) metrics, and constructs autoregressive linear predictors for predictive adaptive optics control.

## Components

- `AtmosphericTurbulence.c` / `AtmosphericTurbulence.h`: Module initialization and milk CLI command dispatchers.
- `atmturb_types.h`: Internal types, simulation configurations (`CONF_*`), physical constants, and cross-file prototypes.
- `atmturb_config.c`: Configuration file reader (`AtmosphericTurbulence_ReadConf`) and air compressibility equations of state (`Z_Air`, `Z_N2`).
- `atmturb_screens.c`: von Karman and power-law master phase screen generator (`make_master_turbulence_screen`).
- `atmturb_wind.c`: 1D von Karman turbulent wind velocity synthesis (`make_AtmosphericTurbulence_vonKarmanWind`).
- `atmturb_hvturb.c`: Hufnagel-Valley $C_n^2$ vertical profile generator with Fried parameter fitting (`AtmosphericTurbulence_makeHV_CN2prof`).
- `atmturb_cube_ops.c`: Complex and phase-only wavefront cube decimation and spatial binning (`contract_wavefront_cube`).
- `atmturb_psf_metrics.c`: Encircled energy aperture photometry, contrast extraction, and lucky imaging frame selection (`measure_wavefront_series`, `frame_select_PSF`).
- `atmturb_psf_sim.c`: End-to-end closed-loop adaptive optics simulation with PID feedback and PSF accumulation (`AtmosphericTurbulence_makePSF`).
- `atmturb_linpred_full.c`: Full-aperture autoregressive linear prediction matrix building and evaluation (`AtmosphericTurbulence_Build_LinPredictor_Full`).
- `atmturb_linpred_pixel.c`: Pixel-level predictor training and shift-invariant 2D kernel extraction (`AtmosphericTurbulence_LinPredictor_filt_2DKernelExtract`).
- `atmturb_linpred_test.c`: Synthetic tip-tilt vibration test sequences and predictor verification (`AtmosphericTurbulence_mkTestTTseq`).
- `atmturb_wfs_sim.c`: Multi-layer extruded wavefront time-series simulation engine (`make_AtmosphericTurbulence_wavefront_series`).
- `atmturb_mkwfs_FPS.c`: FPS V2 compute unit for CLI and standalone execution (`milk-fpsexec-atmturb-mkwfs`).
- `atmturb_mkhvturb_FPS.c`: FPS V2 compute unit for Cn2 profile generation (`milk-fpsexec-atmturb-mkhvturb`).

## Public Headers

- `AtmosphericTurbulence.h`:
  - `int init_AtmosphericTurbulence(void)`
  - `int AtmosphericTurbulence_change_configuration_file(char *fname)`
  - `long make_AtmosphericTurbulence_vonKarmanWind(...)`
  - `int make_master_turbulence_screen(...)`
  - `int make_master_turbulence_screen_pow(...)`
  - `int contract_wavefront_series(...)`
  - `int contract_wavefront_cube(...)`
  - `int contract_wavefront_cube_phaseonly(...)`
  - `int make_AtmosphericTurbulence_wavefront_series(...)`
  - `int measure_wavefront_series(...)`
  - `int AtmosphericTurbulence_mkTestTTseq(...)`
  - `int AtmosphericTurbulence_Build_LinPredictor_Full(...)`
  - `int AtmosphericTurbulence_Apply_LinPredictor_Full(...)`
  - `long AtmosphericTurbulence_LinPredictor_filt_2DKernelExtract(...)`
  - `long AtmosphericTurbulence_LinPredictor_filt_Expand(...)`
  - `int AtmosphericTurbulence_Build_LinPredictor(...)`
  - `long AtmosphericTurbulence_psfCubeContrast(...)`
  - `int AtmosphericTurbulence_Test_LinPredictor(...)`
  - `int measure_wavefront_series_expoframes(...)`
  - `int frame_select_PSF(...)`
  - `int AtmosphericTurbulence_WFprocess(void)`
  - `int AtmosphericTurbulence_makeHV_CN2prof(...)`

## Dependencies

- `CLIcore` (CLI registration, argument parsing, configuration readers)
- `AtmosphereModel` (refractive index and chromatic dispersion scaling)
- `OpticsMaterials` (material dispersion)
- `WFpropagate` (Fresnel diffractive propagation)
- `milkCOREMODarith`, `milkCOREMODmemory`, `milkCOREMODtools`, `milkCOREMODiofits`
- `milkfft`, `milkpsf`, `milklinalgebra`, `milkimagegen`, `milkimagefilter`, `milkstatistic`
- OpenMP (multi-threaded simulation acceleration)
