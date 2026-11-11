# milkatmturb — Atmospheric Turbulence Simulation Plugin for milk

`milkatmturb` is a high-performance atmospheric turbulence and wavefront propagation plugin
for the [milk](https://github.com/milk-org/milk) framework (`feat/bash-engine` branch). It simulates
dynamic optical effects of Earth's atmosphere on astronomical wavefronts across visible and
near-infrared wavelengths, supporting multi-layer Fresnel diffraction, Fourier-space Rytov
diffractive approximation, chromatic atmospheric refraction, real-time shared memory streaming,
closed-loop adaptive optics (AO) simulation, and linear predictive filtering.

Key features include:
- **Hardware Acceleration**: Highly optimized AVX2 and AVX-512 SIMD vectorization, multi-threaded
  OpenMP, and NVIDIA CUDA / cuFFT GPU acceleration with CUDA Graph stream capture.
- **Physical Fidelity**: Kolmogorov and von Karman statistics, Hufnagel-Valley $C_n^2$ profiles,
  3D turbulent wind velocity fields, analytic low-frequency subharmonics (Lane et al.), and
  Moisan periodic-plus-smooth decomposition eliminating boundary edge artifacts.
- **Real-Time Streaming**: Frame-by-frame 2D shared memory streaming via ImageStreamIO with
  sub-millisecond latency, zero-copy GPU DMA, and decoupled simulation physics timing.
- **Adaptive Optics & Metrics**: Integrated closed-loop simulation (leaky integrator, PID) with
  hardware latency modeling, long-exposure science PSF accumulation, and Strehl ratio tracking.

---

## Quick-Start Examples

### 1. From the Shell (Command Line)

Run an end-to-end turbulence simulation and closed-loop AO analysis in 3 simple commands:

```bash
# Step 1: Generate a 20-layer Hufnagel-Valley turbulence profile
# (20 m/s wind, r0=0.15m at 4200m altitude)
milk-fpsexec-atmturb-mkhvturb exec 20.0 0.15 4200 20 turbHV.prof

# Step 2: Run wavefront simulation at 1650 nm (H-band) in single precision
milk-fpsexec-atmturb-mkwfs exec 1650.0 0

# Step 3: Run closed-loop AO simulation and synthesize long-exposure science PSF
milk-fpsexec-atmturb-aoloop exec
```

Or run directly via `milk-cli -c` or `milk.bash`:
```bash
milk-cli -c "atmturb.mkhvturb 20.0 0.15 4200 20 turbHV.prof"
milk-cli -c "atmturb.mkwfs 1650.0 0"
milk-cli -c "atmturb.aoloop"

# Or source milk.bash directly in any GNU Bash script:
source milk.bash
atmturb.mkhvturb 20.0 0.15 4200 20 turbHV.prof
atmturb.mkwfs 1650.0 0
atmturb.aoloop
@calc -v strehl_val "strehl"
echo "Cumulative Strehl = $strehl_val"
```

---

### 2. From Inside the Interactive `milk` CLI

Start the milk interactive shell:
```bash
milk-cli
```

Inside `milk-cli`, run the following sequence:

```text
# 1. Create the vertical turbulence profile
milk-cli > atmturb.mkhvturb 20.0 0.15 4200 20 turbHV.prof

# 2. Generate two 2048x2048 master phase screens (outer scale = 50m, inner scale = 0.01m)
milk-cli > atmturb.mkmastert 2048 50.0 1.0 0 scr0 scr1

# 3. Extrude and propagate wavefronts at 1650 nm
milk-cli > atmturb.mkwfs 1650.0 0

# 4. Run closed-loop AO simulation (leaky integrator, gain 0.5, 1 frame delay)
milk-cli > atmturb.aoloop

# 5. Check generated images in memory
milk-cli > mem.listim

# 6. Save cumulative PSF output to FITS file
milk-cli > iofits.savefits PSFcumul "psfcumul.fits"

# Exit
milk-cli > exit
```

---

## Build and Installation

`milkatmturb` supports both standalone out-of-tree builds against an installed `milk`
(`feat/bash-engine` branch) and in-tree builds inside `milk/plugins/`:

### Option A: Out-of-Tree Build (Recommended)
```bash
cmake -B _build -DCMAKE_BUILD_TYPE=Release -DUSE_CUDA=ON \
    -DCMAKE_INSTALL_PREFIX=/path/to/milk-install
cmake --build _build -j$(nproc)
cmake --install _build   # installs library, stubs, headers, and refreshes milk-cli --catalog
```

### Option B: In-Tree Build Inside `milk/plugins/`
```bash
ln -s /path/to/AtmosphericTurbulenceSimul /path/to/milk/plugins/milkatmturb
cd /path/to/milk
cmake -B _build -DCMAKE_BUILD_TYPE=Release
cmake --build _build --target milkatmturb \
    milk-fpsexec-atmturb-mkwfs milk-fpsexec-atmturb-mkhvturb \
    milk-fpsexec-atmturb-mkmastert milk-fpsexec-atmturb-mkvonkarman \
    milk-fpsexec-atmturb-aoloop -j$(nproc)
cmake --install _build
```

### Build Requirements & Options
- **C Standard**: Builds with C23 (`-std=c23` / `-std=gnu23`) by default, with automatic fallback
  to C11 via `atmturb_c23_compat.h`.
- **OpenMP & FFTW3**: OpenMP multi-threading and single/double precision FFTW (`fftw3`, `fftw3f`)
  are required.
- **SIMD Vectorization**: AVX2 and AVX-512 kernels are compiled on x86_64 architectures; runtime
  CPU capability detection automatically selects the fastest available vector instruction set.
- **CUDA GPU Acceleration**: NVIDIA CUDA Toolkit and cuFFT support is enabled automatically when
  detected, or explicitly via `-DUSE_CUDA=ON`.
- **Dependencies**: Links Level-1 shared modules `milkOpticsMaterials` (chromatic dispersion) and
  `milkWFpropagate` (split-step Fresnel diffraction).

---

## Command Line Usage

### A. Standalone FPS Executables

`milkatmturb` provides dedicated standalone executables using the `milk` FPS V2
(Function Parameter Structure) framework.

#### `milk-fpsexec-atmturb-mkhvturb` — Generate Hufnagel-Valley Profile
Computes a discretized vertical turbulence profile ($C_n^2$, wind velocity,
direction, inner and outer scales).

```bash
# Syntax: milk-fpsexec-atmturb-mkhvturb exec <wspeed> <r0> <sitealt> <nblayers> <outfile>
milk-fpsexec-atmturb-mkhvturb exec 20.0 0.15 4200 20 turbHV.prof
```

| Positional Arg | Keyword | Type | Default | Description |
| :---: | :--- | :--- | :--- | :--- |
| 0 | `.wspeed` | `FLOAT64` | `20.0` | High altitude wind speed [m/s] |
| 1 | `.r0` | `FLOAT64` | `0.15` | Fried parameter $r_0$ [m] at $0.55\,\mu\text{m}$ |
| 2 | `.sitealt` | `FLOAT64` | `4200.0` | Site altitude above sea level [m] |
| 3 | `.nblayers`| `INT32` | `20` | Number of discretized atmospheric layers |
| 4 | `.outfile` | `FILENAME`| `turbHV.prof` | Output profile filename |

#### `milk-fpsexec-atmturb-mkwfs` — Generate Wavefront Series / Stream
Simulates dynamic multi-layer turbulence, propagating light through layers
and outputting FITS cubes or shared memory streams.

```bash
# Syntax: milk-fpsexec-atmturb-mkwfs exec <slambda> <precision>
milk-fpsexec-atmturb-mkwfs exec 1650.0 0
```

| Positional Arg | Keyword | Type | Default | Description |
| :---: | :--- | :--- | :--- | :--- |
| 0 | `.slambda` | `FLOAT32` | `1650.0` | Science wavelength [$\text{nm}$ or $\mu\text{m}$] |
| 1 | `.precision`| `INT32` | `0` | Precision mode (`0` = single `float`, `1` = `double`) |

Key runtime parameters (configurable via `milk-fps-set` or command line options):

| Keyword | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `.stream_mode` | `INT32` | `0` | Mode: `0` = 3D cube, `1` = continuous stream, `2` = finite |
| `.stream_cbsize`| `INT32` | `0` | Circular history depth (`0` = bypass CPU memcpy) |
| `.wfsize` | `INT32` | `512` | Wavefront grid dimension [pixels] |
| `.pupil_scale` | `FLOAT32` | `0.02` | Sampling scale [m/pixel] |
| `.seeing` | `FLOAT32` | `0.65` | Zenith seeing at reference wavelength [arcsec] |
| `.time_step` | `FLOAT32` | `0.001` | Physics time between time slices [seconds] |
| `.time_span` | `FLOAT32` | `0.1` | Duration per file/cube [seconds] |
| `.prof_file` | `FILENAME`| `turbul.prof` | Vertical $C_n^2$ and wind profile file |
| `.master_size` | `INT32` | `4096` | Master phase screen dimension [pixels] |
| `.save_fits` | `INT32` | `0` | Save FITS cubes (`0` = no, `1` = ref, `2` = sci, `3` = all) |
| `.amplitude` | `INT32` | `0` | Compute amplitude scintillation in addition to phase |
| `.fresnel` | `INT32` | `0` | Diffraction mode (`0` = geom, `1` = split-step, `2` = Rytov) |
| `.fresnel_guard`| `INT32` | `0` | Rytov guard band margin [pixels] (`0` = none) |
| `.fresnel_moisan`| `INT32`| `-1` | Moisan decomposition (`-1` = auto, `0` = bypass, `1` = on) |
| `.fresnel_zinterp`| `INT32`| `0` | Piecewise-linear $z$-interpolation between nodes |
| `.fresnel_refract_path`| `INT32`| `0` | Curved ray refracted path length integration |
| `.fresnel_scint_weight`| `INT32`| `0` | Scintillation-weighted superlayer centroid |
| `.interp` | `INT32` | `1` | Extrusion interpolation (`0` = bilinear, `1` = Keys bicubic)|
| `.lowfreq` | `INT32` | `1` | Analytic low-frequency subharmonics (`0` = off, `1` = on) |
| `.rolling` | `INT32` | `1` | Rolling cross-faded screens (`0` = off, `1` = on) |
| `.boil_time` | `FLOAT32` | `0.0` | Maximum epoch duration for cross-fade [s] (`0` = auto) |
| `.out_phase` | `STREAMNAME`| `outarraypha`| Output phase stream name |
| `.out_ampl` | `STREAMNAME`| `outarrayamp`| Output amplitude stream name |
| `.conffile` | `FILENAME`| `WFsim.conf` | Optional simulation configuration file |

#### `milk-fpsexec-atmturb-mkmastert` — Generate Master Turbulence Screens
Generates a pair of normalized master phase screens in the Fourier domain with
outer and inner scale filtering.

```bash
# Syntax: milk-fpsexec-atmturb-mkmastert exec <size> <outerscale> <innerscale> \
#                                            <precision> <screen0> <screen1>
milk-fpsexec-atmturb-mkmastert exec 2048 50.0 1.0 0 scr0 scr1
```

| Positional Arg | Keyword | Type | Default | Description |
| :---: | :--- | :--- | :--- | :--- |
| 0 | `.size` | `INT32` | `2048` | Screen grid dimension [pixels] |
| 1 | `.outerscale` | `FLOAT32` | `100.0` | Outer scale in grid units [pixels] |
| 2 | `.innerscale` | `FLOAT32` | `1.0` | Inner scale in grid units [pixels] |
| 3 | `.precision` | `INT32` | `0` | Precision mode (`0` = single, `1` = double) |
| 4 | `.screen0` | `STREAMNAME` | `turbm00_p0` | Output screen 0 image name |
| 5 | `.screen1` | `STREAMNAME` | `turbm00_p1` | Output screen 1 image name |

#### `milk-fpsexec-atmturb-mkvonkarman` — Generate von Karman Wind Velocity Series
Synthesizes a 3-channel 1D time series $[u, v, w]$ representing longitudinal,
transverse, and vertical turbulent wind fluctuations.

```bash
# Syntax: milk-fpsexec-atmturb-mkvonkarman exec <vksize> <pixscale> <sigmawind> <lwind> \
#                                             <outname> [<seed>] [<fitsout>]
milk-fpsexec-atmturb-mkvonkarman exec 8192 0.1 20.0 50.0 vkwind 1 vkwind.fits
```

| Positional Arg | Keyword | Type | Default | Description |
| :---: | :--- | :--- | :--- | :--- |
| 0 | `.vksize` | `INT32` | `8192` | Sample count of 1D series |
| 1 | `.pixscale` | `FLOAT32` | `0.1` | Physical sampling step [m] |
| 2 | `.sigmawind` | `FLOAT32` | `20.0` | Velocity standard deviation [m/s] |
| 3 | `.lwind` | `FLOAT32` | `50.0` | Turbulence outer scale [m] |
| 4 | `.outname` | `STREAMNAME` | `vKwind` | Output image name ($v_{\text{sz}} \times 1 \times 3$) |
| 5 | `.seed` | `INT64` | `1` | RNG seed (`0` = time-based) |
| 6 | `.fitsout` | `STRING` | `""` | Optional FITS output file (empty = do not save) |

#### `milk-fpsexec-atmturb-aoloop` — Run Adaptive Optics Closed Loop & PSF Synthesis
Applies leaky integrator or PID wavefront correction with latency to simulated
wavefronts and synthesizes Nyquist-sampled long-exposure science PSFs.

```bash
# Syntax: milk-fpsexec-atmturb-aoloop exec [<wfname>] [<psfname>] [<fitsout>] \
#                                         [<loop_mode>] [<gain>] [<leak>] [<delay>] \
#                                         [<scilambda>] [<teldiam>] [<reflambda>] \
#                                         [<pupil_scale>] [<time_step>]
milk-fpsexec-atmturb-aoloop exec outarraypha PSFcumul psf.fits 1 0.5 0.001 1 1.65 8.0 0.50 0.0 0.001
```

| Positional Arg | Keyword | Type | Default | Description |
| :---: | :--- | :--- | :--- | :--- |
| 0 | `.wfname` | `FILENAME` | `outarraypha` | Input wavefront phase cube or stream |
| 1 | `.psfname` | `STREAMNAME`| `PSFcumul` | Output cumulative PSF stream name |
| 2 | `.fitsout` | `STRING` | `""` | Optional FITS output file |
| 3 | `.loop_mode` | `INT32` | `1` | Mode: `0` = Open Loop, `1` = Leaky, `2` = PID |
| 4 | `.gain` | `FLOAT64` | `0.5` | Loop gain |
| 5 | `.leak` | `FLOAT64` | `0.001` | Integrator leak factor |
| 6 | `.delay` | `INT32` | `1` | Hardware latency delay [frames] |
| 7 | `.scilambda` | `FLOAT64` | `1.65` | Science wavelength [$\mu\text{m}$] |
| 8 | `.teldiam` | `FLOAT64` | `8.0` | Telescope primary mirror diameter [m] |
| 9 | `.reflambda` | `FLOAT64` | `0.50` | Reference wavefront wavelength [$\mu\text{m}$] |
| 10 | `.pupil_scale`| `FLOAT64` | `0.0` | Pupil scale [m/pix] (`0` = auto) |
| 11 | `.time_step` | `FLOAT64` | `0.001` | Simulation time step [seconds] |

---

### B. Real-Time 2D Shared Memory Streaming

In streaming mode (`stream_mode = 1` or `2`), `milkatmturb` writes 2D wavefront phase and
amplitude frames directly into ImageStreamIO shared memory (`/milk/shm/outarraypha.im.shm`).
Downstream processes (such as wavefront sensors, AO loops, or telemetry loggers) can synchronize
via POSIX semaphores with sub-millisecond latency.

```bash
# Initialize FPS and enable continuous streaming mode
milk-fpsexec-atmturb-mkwfs fpsinit
milk-fps-set atmturb_mkwfs.stream_mode 1

# Launch detached background streaming compute unit (use `run` for foreground finite streams)
milk-fpsexec-atmturb-mkwfs runstart

# Inspect active shared memory streams
milk-stream-info outarraypha
milk-stream-info outarrayamp

# Dynamically adjust science wavelength or time step on the fly
milk-fps-set atmturb_mkwfs.slambda 850.0
milk-fps-set atmturb_mkwfs.time_step 0.0005

# Stop background streaming execution cleanly
milk-fps-runstop atmturb_mkwfs
```

---

### C. Shell Helper Scripts

Helper scripts are located in `src/AtmosphericTurbulence/scripts/`:

- **`mkHVturbprof`**:
  ```bash
  # Creates turbHV.prof automatically
  ./src/AtmosphericTurbulence/scripts/mkHVturbprof 20.0 0.15 4200 20
  ```

- **`runturb`**:
  Automated driver script that generates `WFsim.conf`, builds the turbulence profile,
  and launches simulation:
  ```bash
  # Run simulation for 1.65 um with median seeing profile (0.65")
  ./src/AtmosphericTurbulence/scripts/runturb -T HVmed 1.65

  # Run simulation for 0.85 um with good seeing profile (0.40") in double precision (-d)
  ./src/AtmosphericTurbulence/scripts/runturb -T HVgood -d 0.85
  ```

---

## Interactive `milk` CLI Usage

Start `milk-cli`:
```bash
milk-cli
```

### Command Discovery
List all commands provided by the `milkatmturb` plugin:
```text
milk-cli > @help milkatmturb
```

Get detailed syntax and documentation for a specific command:
```text
milk-cli > atmturb.mkwfs -h
milk-cli > atmturb.mkhvturb -h
milk-cli > atmturb.mkmastert -h
milk-cli > atmturb.aoloop -h
```

### Registered Commands Reference

| Command | Syntax | Description |
| :--- | :--- | :--- |
| `atmturb.mkhvturb` | `<wspeed> <r0> <sitealt> <nb> [out] [wmod] [seed]` | Make HV profile |
| `atmturb.mkmastert` | `<sz> <L0> <l0> [prec] [scr0] [scr1] [seed] ...` | Master phase screens |
| `atmturb.mkwfs` | `[slambda] [precision] [conffile]` | Wavefront series / 2D stream |
| `atmturb.mkvonkarman` | `<sz> <scale> <sig> <L0> [out] [seed] [fits]` | von Karman wind series |
| `atmturb.aoloop` | `[wf] [psf] [fits] [mode] [gain] ...` | AO closed loop (`strehl` outvar) |
| `atmturb.mkatm` | `[conffile] [slambda]` | Build NRLMSISE-00 atmosphere & refraction |
| `atmturb.measr0` | `<screen>` | Measure $r_0$ [pix] (`r0pix` outvar) |
| `atmturb.linpred_build`| `<inwf> <order> <delay> <svdlim> <outpf>` | Train linear predictor |
| `atmturb.linpred_apply`| `<inwf> <delay> <pf> <outwf>` | Apply linear wavefront predictor |
| `atmturb.contract` | `<incube> <factor> <outcube> [phase_only]` | Decimate wavefront FITS cube |

---

## Configuration File Reference (`WFsim.conf`)

`atmturb.mkwfs` reads default parameters from `WFsim.conf`:

```ini
# --- ATMOSPHERE & SEEING ---
TURBULENCE_REF_WAVEL     0.550000   # Reference wavelength [micron]
TURBULENCE_SEEING        0.650000   # Zenith seeing at reference wavelength [arcsec]
TURBULENCE_PROF_FILE  turbHV.prof   # Vertical Cn2 and wind profile file
ZENITH_ANGLE             0.000000   # Zenith angle [rad]

# --- GRID & SAMPLING ---
WFsize                        512   # Output wavefront pixel dimension (512x512)
PUPIL_SCALE              0.020000   # Spatial resolution [m/pixel]
MASTER_SIZE                  4096   # Master phase screen size [pixels]
MASTER_OVERSAMPLE               2   # Master screen oversampling factor (1 or 2)
INTERP                          1   # Interpolation: 0 = bilinear, 1 = Keys bicubic
LOWFREQ                         1   # 1: Enable analytic low-frequency subharmonics
ROLLING                         1   # 1: Enable rolling cross-faded screens
BOIL_TIME                0.000000   # Maximum epoch duration for cross-fade [s]

# --- TEMPORAL EVOLUTION ---
WFTIME_STEP              0.001000   # Simulation time step [seconds] (1 ms = 1 kHz)
TIME_SPAN                0.100000   # Continuous file time span [seconds]
NB_TSPAN                      100   # Total consecutive time spans to simulate

# --- OUTPUT & PROPAGATION MODES ---
SHM_SOUTPUT                     1   # 1: Stream output to shared memory
SHM_SPREFIX             wfsim_out   # Shared memory stream prefix
SWF_WRITE2DISK                  0   # 1: Save wavefront cubes to FITS files
WAVEFRONT_AMPLITUDE             1   # 1: Compute amplitude (scintillation)
FRESNEL_PROPAGATION             2   # 0: geometric, 1: split-step Fresnel, 2: Rytov
FRESNEL_PROPAGATION_BIN    1000.0   # Vertical superlayer binning distance [meters]
FRESNEL_GUARD_PIX              16   # Continuous guard band margin [pixels]
FRESNEL_RYTOV_MOISAN           -1   # Moisan decomposition: -1 auto, 0 off, 1 on
FRESNEL_RYTOV_ZINT              1   # 1: 2nd-order piecewise-linear z-interpolation
FRESNEL_REFRACT_PATH            0   # 1: Integrate curved ray path length
FRESNEL_SCINT_WEIGHT            0   # 1: Scintillation-weighted centroid
```

---

## Wavefront Propagation & Scintillation Engines

`milkatmturb` supports three propagation regimes selected by `FRESNEL_PROPAGATION`:

- **Mode 0: Geometric Ray Tracing (`FRESNEL_PROPAGATION=0`)**:
  Pure geometric phase accumulation along lines of sight with zero diffractive amplitude
  scintillation. Accelerated via AVX2 and AVX-512 SIMD vectorization and CUDA GPU kernels.
  Supports Keys bicubic interpolation, rolling buffer memory management, and Lane et al.
  low-frequency subharmonics.
- **Mode 1: Multi-Layer Split-Step Fresnel (`FRESNEL_PROPAGATION=1`)**:
  Full optical field diffraction propagated iteratively between altitude superlayers via
  `milkWFpropagate`. Computes amplitude scintillation and phase distortion across thick layers.
  Accelerated with AVX2 SIMD phase update kernels and a full CUDA / cuFFT GPU engine with
  CUDA Graph stream capture and zero-copy direct SHM DMA.
- **Mode 2: Fourier-Space Rytov Diffractive Approximation (`FRESNEL_PROPAGATION=2`)**:
  High-speed analytical multi-layer diffractive propagation in Fourier space:
  - **Moisan Periodic-plus-Smooth Decomposition**: Decomposes extruded phases into periodic and
    smooth harmonic components, eliminating boundary edge artifacts without artificial apodization.
  - **Shared Primary/Secondary Spectrum**: Evaluates a single forward FFT per superlayer shared
    between primary and secondary wavelengths, applying a static chromatic phase ramp.
  - **Continuous Guard Band Padding (`FRESNEL_GUARD_PIX`)**: Eliminates boundary wrap-around.
  - **Piecewise-Linear $z$-Interpolation (`FRESNEL_RYTOV_ZINT`)**: Second-order $O(\Delta z^2)$
    convergence between propagation nodes.
  - **Refracted Ray Path Length (`FRESNEL_REFRACT_PATH`)**: Integrates curved ray trajectories.
  - **Scintillation-Weighted Centroid (`FRESNEL_SCINT_WEIGHT`)**: Concave $z^{5/6}$ weighting.
  - **Hardware Acceleration**: Vectorized AVX2 CPU kernels, batched cuFFT, and CUDA Graph capture.

---

## Performance & Benchmarks

End-to-end runtime benchmarks measured with `tests/benchmark_wavefront_series.sh` and
`tests/benchmark_streaming.sh` on an Intel Xeon / NVIDIA RTX workstation:

### Propagation Mode Comparison (Multi-Layer Turbulence, 5 Layers)

| Grid | Split-Step | Geometric | Rytov AVX2 | Rytov CUDA |
|---|---|---|---|---|
| **256x256** (50 fr) | 116.8 fps (1.0x) | 183.8 fps (1.6x) | 152.0 fps (1.3x) | 124.4 fps (1.1x) |
| **512x512** (20 fr) | 18.7 fps (1.0x) | 26.7 fps (1.4x) | 24.2 fps (1.3x) | 22.2 fps (1.2x) |

### SIMD Vectorization Scaling (Geometric Extrusion)

| Grid Size | SCALAR | AVX2 | CUDA (Render) |
|---|---|---|---|
| **64x64** (200 fr) | 1,242 fps (5.1 MP/s) | 1,258 fps (5.2 MP/s) | 2,131 fps |
| **256x256** (100 fr) | 237 fps (15.5 MP/s) | 265 fps (17.4 MP/s) | 939 fps |
| **512x512** (50 fr) | 52 fps (13.6 MP/s) | 58 fps (15.2 MP/s) | 424 fps |

### Dedicated Benchmark Binaries
- `atmturb-benchmark-perf`: Evaluates CPU and GPU kernel execution times across grid sizes.
- `atmturb-benchmark-stream`: Measures streaming latency and sustained FPS in shared memory.

---

## Testing & Verification

The plugin includes comprehensive test suites covering components, physics, and streaming:

```bash
# Verify standalone FPS executables
bash tests/test_fps_components.sh

# Verify wavefront cube generation and FITS output
bash tests/test_wavefront_series.sh

# Verify real-time 2D shared memory streaming
bash tests/test_wavefront_stream.sh

# Quantitative physics regression tests (Kolmogorov scaling, structure functions)
bash tests/test_turbulence_physics.sh

# Check code size and function scope ratchet
./scripts/check_code_size.sh
```

---

## Directory Organization

```
milkatmturb/
├── CMakeLists.txt              # CMake build definitions (library + standalones)
├── README.md                   # Plugin documentation and usage guide
├── scripts/                    # Code size ratchet & validation scripts
│   ├── check_code_size.sh      # Enforcement script
│   └── code_size_baseline.txt  # Ratchet baseline
├── src/
│   ├── milkatmturb.h           # Master plugin header
│   ├── milkatmturb.c           # MILK_MODULE_RUNNABLES registration & entry point
│   ├── atmturb_compat.h        # milk feat/bash-engine compatibility shims
│   ├── atmturb_c23_compat.h    # C23/C11 compatibility layer
│   ├── AtmosphereModel/        # NRLMSISE-00, vertical density, refractive index
│   │   ├── AtmosphereModel.c/.h
│   │   ├── atmmod_*.c
│   │   └── README.md
│   └── AtmosphericTurbulence/  # Turbulence screens, wavefront simulation, FPS & CLI
│       ├── atmturb_mkwfs_FPS.c      # Wavefront series standalone & CLI runnable
│       ├── atmturb_mkhvturb_FPS.c   # Cn2 profile standalone & CLI runnable
│       ├── atmturb_mkmastert_FPS.c  # Master screen standalone & CLI runnable
│       ├── atmturb_mkvonkarman_FPS.c# Wind velocity standalone & CLI runnable
│       ├── atmturb_aoloop_FPS.c     # AO closed loop standalone & CLI runnable
│       ├── atmturb_cli.c            # CLI utility runnables (mkatm, measr0, linpred, ...)
│       ├── atmturb_simd_*.c         # SIMD scalar, AVX2, AVX-512 kernels
│       ├── atmturb_cuda*.cu         # CUDA GPU kernels (extrusion, Fresnel, Rytov)
│       ├── scripts/                 # runturb and mkHVturbprof shell drivers
│       └── README.md
└── tests/                      # Integration, physics, and benchmark suites
    ├── test_fps_components.sh
    ├── test_wavefront_series.sh
    ├── test_wavefront_stream.sh
    ├── test_turbulence_physics.sh
    ├── benchmark_wavefront_series.sh
    └── benchmark_streaming.sh
```

---

## License & Credits

- Original Authors: Olivier Guyon et al.
- Developed with support from the National Science Foundation (award #1006063).
- License: LGPL-3.0-or-later.
