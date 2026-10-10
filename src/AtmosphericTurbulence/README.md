# AtmosphericTurbulence

**Architectural Layer**: Level 2 (Turbulence Simulation, Wavefront Generation & Filtering)

## Purpose

The `AtmosphericTurbulence` module provides multi-layer atmospheric turbulence simulation
using the frozen-flow (Taylor) hypothesis and von Karman / Kolmogorov phase statistics.
It generates master phase screens, synthesizes 3D turbulent wind velocity fields, computes
Hufnagel-Valley $C_n^2$ profiles, extrudes multi-layer phase screens, calculates Fresnel or
geometric diffractive wavefront cubes, evaluates point spread function (PSF) metrics, and
constructs autoregressive linear predictors for predictive adaptive optics control.

## Components

- `AtmosphericTurbulence.c` / `AtmosphericTurbulence.h`: Module initialization and milk CLI command
  dispatchers.
- `atmturb_types.h`: Internal types, simulation configurations (`CONF_*`), physical constants,
  and cross-file prototypes.
- `atmturb_config.c`: Configuration file reader (`AtmosphericTurbulence_ReadConf`) and air
  compressibility equations of state (`Z_Air`, `Z_N2`).
- `atmturb_screens.c`: von Karman and power-law master phase screen generator
  (`atmturb_generate_screen_pair`, `make_master_turbulence_screen_seeded`).
- `atmturb_wind.c`: 1D von Karman turbulent wind velocity synthesis
  (`make_AtmosphericTurbulence_vonKarmanWind`).
- `atmturb_hvturb.c`: Hufnagel-Valley $C_n^2$ vertical profile generator with Fried parameter
  fitting (`AtmosphericTurbulence_makeHV_CN2prof`).
- `atmturb_cube_ops.c`: Complex and phase-only wavefront cube decimation and spatial binning
  (`contract_wavefront_cube`).
- `atmturb_psf_metrics.c`: Encircled energy aperture photometry, contrast extraction, and lucky
  imaging frame selection (`measure_wavefront_series`, `frame_select_PSF`).
- `atmturb_psf_sim.c`: End-to-end closed-loop adaptive optics simulation with PID feedback and PSF
  accumulation (`AtmosphericTurbulence_makePSF`).
- `atmturb_linpred_full.c`: Full-aperture autoregressive linear prediction matrix building and
  evaluation (`AtmosphericTurbulence_Build_LinPredictor_Full`).
- `atmturb_linpred_pixel.c`: Pixel-level predictor training and shift-invariant 2D kernel
  extraction (`AtmosphericTurbulence_LinPredictor_filt_2DKernelExtract`).
- `atmturb_linpred_test.c`: Synthetic tip-tilt vibration test sequences and predictor verification
  (`AtmosphericTurbulence_mkTestTTseq`).
- `atmturb_wfs_sim.c`: Multi-layer extruded wavefront time-series simulation engine
  (`make_AtmosphericTurbulence_wavefront_series`).
- `atmturb_wfs_render.c` / `atmturb_wfs_render.h`: High-performance frame rendering orchestrator
  dispatching between geometric, split-step Fresnel, CUDA GPU, and Rytov Fourier propagation.
- `atmturb_wfs_stream.c` / `atmturb_wfs_stream.h`: Real-time streaming simulation engine with
  ImageStreamIO shared memory synchronization and persistent propagation engines.
- `atmturb_geometry.c` / `atmturb_geometry.h`: Chromatic atmospheric refraction ray bending,
  dispersion offsets, and curved ray path length integration.
- `atmturb_superlayer.c` / `atmturb_superlayer.h`: Vertical altitude binning, scintillation-weighted
  centroids, and super-layer reduction for multi-layer diffractive propagation.
- `atmturb_superlayer_nodes.c`: Piecewise-linear hat interpolation node generation for
  2nd-order $O(\Delta z^2)$ distance convergence.
- `atmturb_fresnel.c` / `atmturb_fresnel.h`: Multi-layer split-step Fresnel diffractive propagation
  engine (`FRESNEL_PROPAGATION=1`).
- `atmturb_rytov.c`, `atmturb_rytov_kernels.c`, `atmturb_rytov.h`: First-order Rytov Fourier-space
  propagation engine (`FRESNEL_PROPAGATION=2`) with Moisan periodic-plus-smooth decomposition.
- `atmturb_cuda_rytov.cu`, `atmturb_cuda_rytov_kernels.cu`, `atmturb_cuda_rytov.h`: CUDA and cuFFT
  GPU acceleration for Rytov wavefront synthesis, batched transforms, and CUDA Graph execution.
- `atmturb_cuda.cu` / `atmturb_cuda.h`: CUDA GPU accelerated multi-layer wavefront extrusion kernel.
- `atmturb_simd.h`: Declarations for SIMD-accelerated extrusion, scaling, initialization, and
  ISA queries.
- `atmturb_simd_scalar.c`: Portable scalar reference implementation for extrusion and array ops.
- `atmturb_simd_avx2.c`: AVX2 and FMA 256-bit vectorized compute kernels (8 floats/vector).
- `atmturb_simd_avx512.c`: AVX-512 512-bit vectorized compute kernels (16 floats/vector).
- `atmturb_simd_dispatch.c`: Dynamic CPU capability detector and function pointer dispatcher
  (`atmturb_simd_active_isa`).
- `atmturb_mkwfs_FPS.c`: FPS V2 compute unit for CLI and standalone execution
  (`milk-fpsexec-atmturb-mkwfs`).
- `atmturb_mkhvturb_FPS.c`: FPS V2 compute unit for Cn2 profile generation
  (`milk-fpsexec-atmturb-mkhvturb`).
- `atmturb_mkmastert_FPS.c`: FPS V2 compute unit for master turbulence screens
  (`milk-fpsexec-atmturb-mkmastert`).
- `atmturb_mkvonkarman_FPS.c`: FPS V2 compute unit for von Karman wind synthesis
  (`milk-fpsexec-atmturb-mkvonkarman`).

## Public Headers

- `AtmosphericTurbulence.h`:
  - `int init_AtmosphericTurbulence(void)`
  - `int AtmosphericTurbulence_change_configuration_file(const char *fname)`
  - `long make_AtmosphericTurbulence_vonKarmanWind(...)`
  - `int atmturb_generate_screen_pair(...)`
  - `int make_master_turbulence_screen_seeded(...)`
  - `double atmturb_measure_r0_pix(...)`
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
