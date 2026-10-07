# milkatmturb — Atmospheric Turbulence Simulation Plugin for milk

`milkatmturb` is a high-performance atmospheric turbulence and wavefront propagation plugin
for the [milk](https://github.com/milk-org/milk) framework (`framework-dev` branch). It simulates
the dynamic optical effects of Earth's atmosphere on astronomical wavefronts across visible and
near-infrared wavelengths, supporting multi-layer Fresnel diffraction, chromatic dispersion,
real-time shared memory streaming, and linear predictive filtering.

---

## Quick-Start Examples

### 1. From the Shell (Command Line)

You can run an end-to-end simulation in 2 simple commands:

```bash
# Step 1: Generate a 20-layer Hufnagel-Valley turbulence profile
# (20 m/s wind, r0=0.15m at 4200m altitude)
milk-fpsexec-atmturb-mkhvturb exec 20.0 0.15 4200 20 turbHV.prof

# Step 2: Run simulation at 1650 nm (H-band) in single precision
milk-fpsexec-atmturb-mkwfs exec 1650.0 0
```

Or run directly via `milk-cli -c`:
```bash
milk-cli -c "atmturb.mkHVturbprof 20.0 0.15 4200 20 turbHV.prof"
milk-cli -c "atmturb.mkwfs 1650.0 0"
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
milk-cli > atmturb.mkHVturbprof 20.0 0.15 4200 20 turbHV.prof

# 2. Generate two 2048x2048 master phase screens (outer scale = 50m, inner scale = 0.01m)
milk-cli > atmturb.mkmastert scr0 scr1 2048 50.0 0.01

# 3. Extrude and propagate wavefronts at 1650 nm
milk-cli > atmturb.mkwfs 1650.0 0

# 4. Check generated images in memory
milk-cli > listim

# 5. Save output to FITS file
milk-cli > savefits wfphase "wfphase.fits"

# Exit
milk-cli > exitCLI
```

---

## Build and Installation

`milkatmturb` is built directly as a plugin inside `milk` (`framework-dev` branch):

1. **Clone or symlink this repository into `milk/plugins/`**:
   ```bash
   ln -s /path/to/AtmosphericTurbulenceSimul /path/to/milk/plugins/milkatmturb
   ```

2. **Configure and build**:
   ```bash
   cd /path/to/milk
   cmake -B _build -DCMAKE_BUILD_TYPE=Release
   cmake --build _build --target milkatmturb \
       milk-fpsexec-atmturb-mkwfs milk-fpsexec-atmturb-mkhvturb \
       milk-fpsexec-atmturb-mkmastert milk-fpsexec-atmturb-mkvonkarman \
       milk-fpsexec-wfprop-fresnel -j
   ```

3. **Install to system**:
   ```bash
   sudo cmake --install _build
   ```

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

#### `milk-fpsexec-atmturb-mkwfs` — Generate Wavefront Series
Simulates dynamic multi-layer turbulence, propagating light through layers
and outputting FITS cubes or shared memory streams.

```bash
# Syntax: milk-fpsexec-atmturb-mkwfs exec <slambda> <precision> [conffile]
milk-fpsexec-atmturb-mkwfs exec 1650.0 0 WFsim.conf
```

| Positional Arg | Keyword | Type | Default | Description |
| :---: | :--- | :--- | :--- | :--- |
| 0 | `.slambda` | `FLOAT32` | `1650.0` | Science wavelength [$\text{nm}$ or $\mu\text{m}$] |
| 1 | `.precision`| `INT32` | `0` | Precision mode (`0` = single precision `float`, `1` = `double`) |
| - | `.conffile` | `FILENAME`| `WFsim.conf` | Simulation configuration file |

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
# Syntax: milk-fpsexec-atmturb-mkvonkarman exec <vksize> <pixscale> <sigmawind> <lwind> <outname>
milk-fpsexec-atmturb-mkvonkarman exec 8192 0.1 20.0 50.0 vkwind
```

| Positional Arg | Keyword | Type | Default | Description |
| :---: | :--- | :--- | :--- | :--- |
| 0 | `.vksize` | `INT32` | `8192` | Sample count of 1D series |
| 1 | `.pixscale` | `FLOAT32` | `0.1` | Physical sampling step [m] |
| 2 | `.sigmawind` | `FLOAT32` | `20.0` | Velocity standard deviation [m/s] |
| 3 | `.lwind` | `FLOAT32` | `50.0` | Turbulence outer scale [m] |
| 4 | `.outname` | `STREAMNAME` | `vKwind` | Output 3D image name ($v_{\text{size}} \times 1 \times 3$) |

#### `milk-fpsexec-wfprop-fresnel` — Fresnel Wavefront Propagation
Propagates a 2D complex optical field across distance $z$ using the Fourier
Fresnel quadratic phase transfer function.

```bash
# Syntax: milk-fpsexec-wfprop-fresnel exec <inname> <outname> <pupilscale> <distance> <lambda>
milk-fpsexec-wfprop-fresnel exec wfin wfout 0.01 1000.0 0.5e-6
```

| Positional Arg | Keyword | Type | Default | Description |
| :---: | :--- | :--- | :--- | :--- |
| 0 | `.inname` | `STREAMNAME` | `wfin` | Input complex image name |
| 1 | `.outname` | `STREAMNAME` | `wfout` | Output complex image name |
| 2 | `.pupilscale` | `FLOAT64` | `0.01` | Pupil sampling scale [m/pixel] |
| 3 | `.distance` | `FLOAT64` | `1000.0` | Propagation distance [m] |
| 4 | `.lambda` | `FLOAT64` | `0.5e-6` | Optical wavelength [m] |

#### Background Daemon Mode (tmux)
Run continuous simulation in the background and control it in real-time:
```bash
# Launch in detached tmux session
milk-fpsexec-atmturb-mkwfs -tmux runstart

# Check current status and FPS parameters
milk-fps-info atmturb_mkwfs

# Dynamically change wavelength to 850 nm while running
milk-fps-set atmturb_mkwfs.slambda 850.0

# Stop background execution
milk-fpsexec-atmturb-mkwfs runstop
```

---

### B. Shell Helper Scripts

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

## Interactive milk CLI Usage

Start `milk-cli`:
```bash
milk-cli
```

### Command Discovery
List all commands provided by the `milkatmturb` plugin:
```text
milk-cli > m? milkatmturb
```

Get detailed syntax and documentation for a specific command:
```text
milk-cli > cmd? atmturb.mkwfs
milk-cli > cmd? atmturb.mkHVturbprof
milk-cli > cmd? atmturb.mkmastert
```

### Complete Command Reference

| Command | Syntax | Description |
| :--- | :--- | :--- |
| `atmturb.mkHVturbprof` | `<wspeed> <r0> <sitealt> <nblayers> <outfile>` | Make Hufnagel-Valley turbulence profile |
| `atmturb.mkmastert` | `<scr0> <scr1> <size> <outerscale> <innerscale>` | Generate 2 master phase screens |
| `atmturb.mkwfs` | `<wavelength_nm> <precision>` | Generate wavefront series (`0`: float, `1`: double) |
| `atmturb.mkwfs_fps` | `<wavelength_nm> <precision>` | Generate wavefront series with FPS process tracking |
| `atmturb.mkhvturb_fps` | `<wspeed> <r0> <sitealt> <nblayers> [outfile]` | Generate HV profile via FPS |
| `atmturb.mkatmospheremodel` | `<conffile>` | Create vertical atmospheric composition model |
| `atmturb.mkvonKarmanWind` | `<vKsize> <pixscale> <sigmawind> <Lwind> <size> <out>` | Generate von Kármán wind model screen |
| `atmturb.fresnelpw` | `<in_re> <in_im> <z> <lambda> <out_re> <out_im>` | Fresnel propagate optical field |
| `atmturb.atmturbmeasexpo` | `<etime_s> <out_name>` | Measure long-exposure PSF from wavefront series |
| `atmturb.atmturbwfpredictf` | `<in> <mask> <order> <lag> <svdeps> <reglambda>` | Build linear predictor from wavefront series |
| `atmturb.atmturbwfpapply` | `<mode> <in> <mask> <order> <lag> <outp> <outf>` | Apply linear predictor filter |

---

## Configuration File Reference (`WFsim.conf`)

`atmturb.mkwfs` reads parameters from `WFsim.conf`. A typical configuration includes:

```ini
# --- ATMOSPHERE & SEEING ---
TURBULENCE_REF_WAVEL     0.550000   # Reference wavelength [micron]
TURBULENCE_SEEING        0.650000   # Zenith seeing at reference wavelength [arcsec]
TURBULENCE_PROF_FILE  turbHV.prof   # Vertical Cn2 and wind profile file
ZENITH_ANGLE             0.000000   # Zenith angle [rad]

# --- GRID & SAMPLING ---
WFsize                        512   # Output wavefront pixel dimension (512x512)
PUPIL_SCALE              0.020000   # Spatial resolution [m/pixel] (e.g. 512 * 0.02m = 10.24m)
MASTER_SIZE                  4096   # Master phase screen size [pixels]

# --- TEMPORAL EVOLUTION ---
WFTIME_STEP              0.001000   # Simulation time step [seconds] (1 ms = 1 kHz)
TIME_SPAN                0.100000   # Continuous file time span [seconds]
NB_TSPAN                      100   # Total consecutive time spans to simulate

# --- OUTPUT MODES ---
SHM_SOUTPUT                     1   # 1: Stream output to shared memory
SHM_SPREFIX             wfsim_out   # Shared memory stream prefix
SWF_WRITE2DISK                  0   # 1: Save wavefront cubes to FITS files
WAVEFRONT_AMPLITUDE             1   # 1: Compute wavefront amplitude (scintillation)
FRESNEL_PROPAGATION             1   # 1: Enable chromatic diffractive layer propagation
```

---

## Directory Organization

```
milkatmturb/
├── CMakeLists.txt              # CMake build definitions (library + standalones)
├── README.md                   # Plugin documentation and usage guide
└── src/
    ├── milkatmturb.h           # Master plugin header
    ├── milkatmturb.c           # MILK_MODULE registration & entry point
    ├── atmturb_compat.h        # milk framework-dev compatibility shims
    ├── AtmosphereModel/        # NRLMSISE-00, vertical density, refractive index
    ├── AtmosphericTurbulence/  # Turbulence screens, wavefront simulation, FPS
    │   ├── atmturb_mkwfs_FPS.c     # Standalone binary + CLI FPS unit for mkwfs
    │   ├── atmturb_mkhvturb_FPS.c  # Standalone binary + CLI FPS unit for mkHVturbprof
    │   └── scripts/                # runturb and mkHVturbprof shell drivers
    ├── OpticsMaterials/        # Sellmeier dispersion relations
    └── WFpropagate/            # Fresnel diffractive optical propagation
```

---

## License & Credits

- Original Authors: Olivier Guyon et al.
- Developed with support from the National Science Foundation (award #1006063).
- License: LGPL-3.0-or-later.
