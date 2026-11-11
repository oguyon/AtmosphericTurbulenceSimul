# AtmosphericTurbulence

**Architectural Layer**: Level 2 (Turbulence Simulation, Wavefront Generation & Filtering)

## Purpose

The `AtmosphericTurbulence` module provides multi-layer atmospheric turbulence simulation
using the frozen-flow (Taylor) hypothesis and von Karman / Kolmogorov phase statistics.
It generates master phase screens, synthesizes 3D turbulent wind velocity fields, computes
Hufnagel-Valley $C_n^2$ profiles, extrudes multi-layer phase screens, calculates Fresnel or
geometric diffractive wavefront cubes, evaluates point spread function (PSF) metrics,
streams 2D wavefront frames in real time to shared memory, and performs closed-loop adaptive
optics simulation with PSF synthesis.

## Components

- `AtmosphericTurbulence.c` / `AtmosphericTurbulence.h`: Module initialization and public API.
- `atmturb_types.h`: Internal types, simulation configurations (`CONF_*`), physical constants,
  and cross-file prototypes.
- `atmturb_config.c`: Configuration file reader (`AtmosphericTurbulence_ReadConf`) and air
  compressibility equations of state (`Z_Air`, `Z_N2`).
- `atmturb_profile.c` / `atmturb_profile.h`: Discrete altitude turbulence profile parser,
  validation, and sorting.
- `atmturb_screens.c`: von Karman and power-law master phase screen generator
  (`atmturb_generate_screen_pair`, `make_master_turbulence_screen_seeded`).
- `atmturb_wind.c`: 1D von Karman turbulent wind velocity synthesis and trajectory computation
  (`make_AtmosphericTurbulence_vonKarmanWind`, `atmturb_wind_synthesize_trajectory`).
- `atmturb_hvturb.c`: Hufnagel-Valley $C_n^2$ vertical profile generator with Fried parameter
  fitting (`AtmosphericTurbulence_makeHV_CN2prof`, `AtmosphericTurbulence_makeHV_CN2prof_opt`).
- `atmturb_cube_ops.c`: Complex and phase-only wavefront cube decimation and spatial binning
  (`contract_wavefront_cube`).
- `atmturb_psf_metrics.c`: Encircled energy aperture photometry, contrast extraction, and lucky
  imaging frame selection (`measure_wavefront_series`, `frame_select_PSF`).
- `atmturb_psf_sim.c`: End-to-end closed-loop adaptive optics simulation with PID feedback,
  latency modeling, and PSF accumulation (`atmturb_ao_sim_run`, `AtmosphericTurbulence_makePSF`).
- `atmturb_linpred_full.c`: Full-aperture autoregressive linear prediction matrix building and
  evaluation (`AtmosphericTurbulence_Build_LinPredictor_Full`).
- `atmturb_linpred_pixel.c`: Pixel-level predictor training and shift-invariant 2D kernel
  extraction (`AtmosphericTurbulence_LinPredictor_filt_2DKernelExtract`).
- `atmturb_linpred_test.c`: Synthetic tip-tilt vibration test sequences and predictor verification
  (`AtmosphericTurbulence_mkTestTTseq`).
- `atmturb_lowfreq.c` / `atmturb_lowfreq.h`: Analytic low-frequency subharmonic modes
  synthesis (Lane et al.) correcting Fourier tilt under-representation.
- `atmturb_rolling.c` / `atmturb_rolling.h`: Rolling cross-faded buffer memory management
  for seamless continuous phase screen evolution.
- `atmturb_wfs_sim.c`: Multi-layer extruded wavefront time-series simulation engine
  (`make_AtmosphericTurbulence_wavefront_series`).
- `atmturb_wfs_render.c` / `atmturb_wfs_render.h`: High-performance frame rendering orchestrator
  dispatching between geometric, split-step Fresnel, CUDA GPU, and Rytov Fourier propagation.
- `atmturb_wfs_render_cuda.c` / `atmturb_wfs_render_cuda.h`: CUDA rendering pipeline orchestrator
  for batch wavefront cube rendering.
- `atmturb_wfs_cuda_common.c` / `atmturb_wfs_cuda_common.h`: Shared CUDA memory pool allocation,
  device capability queries, and host-device memory transfers.
- `atmturb_wfs_stream.c`, `atmturb_wfs_stream.h`, `atmturb_wfs_stream_internal.h`: Real-time 2D
  streaming simulation engine with ImageStreamIO shared memory synchronization, persistent
  propagation engines, and internal stream-state definitions.
- `atmturb_wfs_stream_cuda.c`: CUDA accelerated real-time streaming engine with zero-copy and
  direct shared memory DMA transfers.
- `atmturb_geometry.c` / `atmturb_geometry.h`: Chromatic atmospheric refraction ray bending,
  dispersion offsets, and curved ray path length integration.
- `atmturb_superlayer.c`, `atmturb_superlayer.h`, `atmturb_superlayer_internal.h`: Vertical
  altitude binning, scintillation-weighted centroids, and super-layer reduction for multi-layer
  diffractive propagation.
- `atmturb_superlayer_nodes.c`: Piecewise-linear hat interpolation node generation for
  2nd-order $O(\Delta z^2)$ distance convergence.
- `atmturb_fresnel.c` / `atmturb_fresnel.h`: Multi-layer split-step Fresnel diffractive propagation
  engine (`FRESNEL_PROPAGATION=1`).
- `atmturb_fresnel_simd_scalar.c`, `atmturb_fresnel_simd_avx2.c`,
  `atmturb_fresnel_simd_dispatch.c`, `atmturb_fresnel_simd.h`: Vectorized split-step Fresnel
  phase addition and normalization kernels.
- `atmturb_rytov.c`, `atmturb_rytov_kernels.c`, `atmturb_rytov.h`, `atmturb_rytov_internal.h`:
  First-order Rytov Fourier-space propagation engine (`FRESNEL_PROPAGATION=2`) with Moisan
  periodic-plus-smooth decomposition.
- `atmturb_rytov_simd_scalar.c`, `atmturb_rytov_simd_avx2.c`,
  `atmturb_rytov_simd_dispatch.c`, `atmturb_rytov_simd.h`: Vectorized Rytov diffractive quadratic
  phase filters and complex multiply-accumulate kernels.
- `atmturb_cuda.cu` / `atmturb_cuda.h`: CUDA GPU accelerated multi-layer wavefront extrusion and
  spatial interpolation kernels.
- `atmturb_cuda_piston.cu` / `atmturb_cuda_piston.h`: CUDA GPU piston subtraction kernel.
- `atmturb_cuda_fresnel.cu`, `atmturb_cuda_fresnel_kernels.cu`, `atmturb_cuda_fresnel.h`,
  `atmturb_cuda_fresnel_kernels.h`: CUDA GPU split-step Fresnel propagation engine with CUDA Graph
  capture.
- `atmturb_cuda_rytov.cu`, `atmturb_cuda_rytov_extrude.cu`, `atmturb_cuda_rytov_filter.cu`,
  `atmturb_cuda_rytov.h`, `atmturb_cuda_rytov_internal.h`: CUDA and cuFFT GPU acceleration for
  Rytov wavefront synthesis, batched transforms, and CUDA Graph execution.
- `atmturb_simd.h`: Declarations for SIMD-accelerated extrusion, scaling, initialization, and
  ISA queries.
- `atmturb_simd_scalar.c`: Portable scalar reference implementation for extrusion and array ops.
- `atmturb_simd_avx2.c`: AVX2 and FMA 256-bit 2D phase-screen extrusion kernels (bilinear/bicubic).
- `atmturb_simd_avx2_ops.c`: AVX2 and FMA 256-bit 1D array, low-frequency subharmonic, and piston
  subtraction kernels.
- `atmturb_simd_avx512.c`: AVX-512 512-bit vectorized compute kernels (16 floats/vector).
- `atmturb_simd_dispatch.c`: Dynamic CPU capability detector and function pointer dispatcher
  (`atmturb_simd_active_isa`).
- `atmturb_mkwfs_FPS.c` / `atmturb_mkwfs_FPS.h`: `MILK_RUNNABLE_FPS` compute unit for CLI and
  standalone execution (`milk-fpsexec-atmturb-mkwfs`, `atmturb.mkwfs`).
- `atmturb_mkhvturb_FPS.c` / `atmturb_mkhvturb_FPS.h`: `MILK_RUNNABLE_FPS` compute unit for
  $C_n^2$ profile generation (`milk-fpsexec-atmturb-mkhvturb`, `atmturb.mkhvturb`).
- `atmturb_mkmastert_FPS.c` / `atmturb_mkmastert_FPS.h`: `MILK_RUNNABLE_FPS` compute unit for
  master turbulence screens (`milk-fpsexec-atmturb-mkmastert`, `atmturb.mkmastert`).
- `atmturb_mkvonkarman_FPS.c` / `atmturb_mkvonkarman_FPS.h`: `MILK_RUNNABLE_FPS` compute unit
  for von Karman wind synthesis (`milk-fpsexec-atmturb-mkvonkarman`, `atmturb.mkvonkarman`).
- `atmturb_aoloop_FPS.c` / `atmturb_aoloop_FPS.h`: `MILK_RUNNABLE_FPS` compute unit for
  closed-loop AO simulation and PSF synthesis (`milk-fpsexec-atmturb-aoloop`, `atmturb.aoloop`).
- `atmturb_cli.c` / `atmturb_cli.h`: `MILK_RUNNABLE_CLI` command wrappers (`atmturb.mkatm`,
  `atmturb.measr0`, `atmturb.linpred_build`, `atmturb.linpred_apply`, `atmturb.contract`).

## Public Headers

- `AtmosphericTurbulence.h`:
  - `int init_AtmosphericTurbulence(void)`
  - `int AtmosphericTurbulence_change_configuration_file(const char *fname)`
  - `long make_AtmosphericTurbulence_vonKarmanWind(...)`
  - `int atmturb_wind_synthesize_trajectory(...)`
  - `int atmturb_generate_screen_pair(...)`
  - `int make_master_turbulence_screen_seeded(...)`
  - `double atmturb_measure_r0_pix(...)`
  - `int make_master_turbulence_screen(...)`
  - `int make_master_turbulence_screen_pow(...)`
  - `int contract_wavefront_series(...)`
  - `int contract_wavefront_cube(...)`
  - `int contract_wavefront_cube_phaseonly(...)`
  - `int make_AtmosphericTurbulence_wavefront_series(...)`
  - `int make_AtmosphericTurbulence_wavefront_stream(...)`
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
  - `void atmturb_ao_init_params(atmturb_ao_params_t *params)`
  - `int atmturb_ao_sim_run(...)`
  - `double AtmosphericTurbulence_makePSF(...)`
  - `int AtmosphericTurbulence_makeHV_CN2prof(...)`
  - `int AtmosphericTurbulence_makeHV_CN2prof_opt(...)`

## Dependencies

- `milkcomputeunit`, `milkfps`, `milkprocessinfo`, `milkdata`, `ImageStreamIO` (`milk_runnable.h`)
- `AtmosphereModel` (refractive index and chromatic dispersion scaling)
- `OpticsMaterials` (material dispersion)
- `WFpropagate` (Fresnel diffractive propagation)
- `milkCOREMODarith`, `milkCOREMODmemory`, `milkCOREMODtools`, `milkCOREMODiofits`
- `milkfft`, `milkpsf`, `milklinalgebra`, `milkimagegen`, `milkimagefilter`, `milkstatistic`
- OpenMP (multi-threaded simulation acceleration)
- FFTW3 / FFTW3F (Fourier domain phase screen generation and Rytov propagation)
- NVIDIA CUDA Toolkit & cuFFT (optional GPU acceleration)
