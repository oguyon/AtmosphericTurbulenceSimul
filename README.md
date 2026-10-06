# milkatmturb — Atmospheric Turbulence Simulation Plugin for milk

`milkatmturb` is a high-performance atmospheric turbulence and wavefront propagation plugin for the [milk](https://github.com/milk-org/milk) framework (`framework-dev` branch). It simulates the dynamic optical effects of Earth's atmosphere on astronomical wavefronts across visible and near-infrared wavelengths, supporting multi-layer Fresnel propagation, chromatic dispersion, real-time shared memory streaming, and linear predictive filtering.

---

## Quick Start: Build and Installation

`milkatmturb` builds directly as a plugin inside `milk`:

1. **Symlink or copy into `milk/plugins/`**:
   ```bash
   ln -s /path/to/AtmosphericTurbulenceSimul /path/to/milk/plugins/milkatmturb
   ```

2. **Build with milk**:
   ```bash
   cd /path/to/milk
   cmake -B _build -DCMAKE_BUILD_TYPE=Release
   cmake --build _build --target milkatmturb milk-fpsexec-atmturb-mkwfs milk-fpsexec-atmturb-mkhvturb -j
   ```

3. **Install**:
   ```bash
   sudo cmake --install _build
   ```

---

## 1. Using from the Command Line (Shell)

You can run atmospheric turbulence simulations directly from bash/zsh without entering an interactive prompt.

### A. Standalone FPS Executables

The plugin provides dedicated standalone executables following the `milk` FPS V2 (Function Parameter Structure) pattern.

#### `milk-fpsexec-atmturb-mkhvturb` — Generate Hufnagel-Valley Profile
Computes a discretized vertical turbulence profile ($C_n^2$, wind velocity, direction, inner and outer scales).

- **Inspect parameters and usage**:
  ```bash
  milk-fpsexec-atmturb-mkhvturb --help
  ```

- **One-shot execution (`exec`)**:
  ```bash
  # Syntax: milk-fpsexec-atmturb-mkhvturb exec <wspeed> <r0> <sitealt> <nblayers> <outfile>
  milk-fpsexec-atmturb-mkhvturb exec 20.0 0.15 4200 20 turbHV.prof
  ```
  | Parameter | Keyword | Type | Default | Description |
  | :--- | :--- | :--- | :--- | :--- |
  | 0 | `.wspeed` | `FLOAT64` | `20.0` | High altitude wind speed [m/s] |
  | 1 | `.r0` | `FLOAT64` | `0.15` | Fried parameter $r_0$ [m] at $0.55\,\mu\text{m}$ |
  | 2 | `.sitealt` | `FLOAT64` | `4200.0` | Site altitude above sea level [m] |
  | 3 | `.nblayers`| `INT32`   | `20` | Number of discretized atmospheric layers |
  | 4 | `.outfile` | `FILENAME`| `turbHV.prof` | Output text profile filename |

#### `milk-fpsexec-atmturb-mkwfs` — Generate Wavefront Series
Generates multi-layer turbulent wavefront sequences, propagating light through layers and writing to FITS cubes or shared memory streams.

- **Inspect parameters and usage**:
  ```bash
  milk-fpsexec-atmturb-mkwfs --help
  ```

- **One-shot execution (`exec`)**:
  ```bash
  # Syntax: milk-fpsexec-atmturb-mkwfs exec <slambda> <precision> [conffile]
  milk-fpsexec-atmturb-mkwfs exec 1650.0 0
  ```
  | Parameter | Keyword | Type | Default | Description |
  | :--- | :--- | :--- | :--- | :--- |
  | 0 | `.slambda` | `FLOAT32` | `1650.0` | Science wavelength [$\text{nm}$ or $\mu\text{m}$] |
  | 1 | `.precision`| `INT32` | `0` | `0`: single precision (`float`), `1`: double precision |
  | - | `.conffile` | `FILENAME`| `WFsim.conf` | Simulation configuration file |

- **Continuous Background Run (tmux / daemon mode)**:
  ```bash
  # Start in background tmux session
  milk-fpsexec-atmturb-mkwfs -tmux runstart

  # Inspect live status
  milk-fps-info atmturb_mkwfs

  # Stop execution
  milk-fpsexec-atmturb-mkwfs runstop
  ```

- **Tune parameters in real-time**:
  ```bash
  milk-fps-set atmturb_mkwfs.slambda 850.0
  ```

---

### B. One-Liner CLI Commands via `milk-cli -c`

Any plugin command registered in `milk` can be executed non-interactively using `milk-cli -c`:

```bash
# Generate turbulence profile
milk-cli -c "atmturb.mkHVturbprof 20.0 0.15 4200 20 turbHV.prof"

# Create master turbulence phase screens
milk-cli -c "atmturb.mkmastert scr0 scr1 2048 50.0 2.0"

# Generate wavefront series
milk-cli -c "atmturb.mkwfs 1650.0 0"
```

---

### C. Driver Shell Scripts

Convenience scripts are provided under `src/AtmosphericTurbulence/scripts/`:

- **`mkHVturbprof <wspeed> <r0> <sitealt> <nblayers>`**:
  Creates `turbHV.prof` automatically using the standalone FPS executable or `milk-cli`.
  ```bash
  ./mkHVturbprof 20.0 0.15 4200 20
  ```

- **`runturb [-hNdT] <slambdaum>`**:
  Configures `WFsim.conf`, generates turbulence profiles, and executes wavefront propagation.
  ```bash
  # Run for 1.65 um with median seeing profile
  ./runturb -T HVmed 1.65
  ```

---

## 2. Using from the Interactive milk CLI

Launch the interactive shell:
```bash
milk-cli
```

### A. Discovering Commands

Once inside `milk-cli`, `milkatmturb` registers all routines under the `atmturb` prefix:

```text
milk-cli > m? milkatmturb
```
This prints module metadata and all available commands:
```text
---- MODULE milkatmturb COMMANDS ---------
   atmturb.mkatmospheremodel  make Earth atmosphere model
   atmturb.fresnelpw          Fresnel propagate wavefront
   atmturb.mkwfs              make wavefront series
   atmturb.mkvonKarmanWind    make vonKarman wind model
   atmturb.mkmastert          make 2 master phase screens
   atmturb.mkHVturbprof       make Hufnager-Valley turbulence profile
   atmturb.atmturbmeasexpo    Measure long exposure time PSF from wavefront series
   atmturb.atmturbmktestTTs   make test TT sequence
   atmturb.atmturbwfpredictf  build full linear predictor from wavefront series
   atmturb.atmturbwfpapply    Apply full linear predictor from wavefront series
   atmturb.atmturbwfp2Dkern   collapse WF predictor into 2D kernel
   atmturb.atmturbwfpexp      Expand 3D filter cube into pixel-based filter
   atmturb.atmturbwfpredict   build linear predictor from wavefront series
   atmturb.atmturbmkpsfcc     measure contrast performance of WF cube
   atmturb.atmturbwfptest     Test linear predictor on wavefront series
   atmturb.mkwfs_fps          Generate atmospheric turbulence wavefront series (FPS)
   atmturb.mkhvturb_fps       Generate Hufnagel-Valley turbulence profile (FPS)
```

To see detailed help and syntax for any command:
```text
milk-cli > cmd? atmturb.mkwfs
milk-cli > cmd? atmturb.mkHVturbprof
```

---

### B. Interactive Simulation Workflow

#### 1. Generate an Atmosphere & Turbulence Profile
```text
milk-cli > atmturb.mkHVturbprof 20.0 0.15 4200 20 turbHV.prof
```

#### 2. Generate Master Phase Screens
Generate two uncorrelated master Kolmogorov / von Kármán phase screens (e.g. $2048 \times 2048$ pixels with $L_0 = 50\,\text{m}$, $l_0 = 2\,\text{m}$):
```text
milk-cli > atmturb.mkmastert scr0 scr1 2048 50.0 2.0
```

#### 3. Run Wavefront Propagation
Run the multi-layer wavefront extrusion and Fresnel propagation at $\lambda = 1650\,\text{nm}$:
```text
milk-cli > atmturb.mkwfs 1650.0 0
```
Or run through the unified FPS engine with live telemetry tracking:
```text
milk-cli > atmturb.mkwfs_fps 1650.0 0
```

#### 4. Inspect Output Images and Shared Memory Streams
Inspect created images and shared memory buffers:
```text
milk-cli > listim
```
Save output streams to FITS:
```text
milk-cli > savefits wfphase "wfphase.fits"
```

Exit the shell:
```text
milk-cli > exitCLI
```

---

## Configuration File (`WFsim.conf`)

`atmturb.mkwfs` reads parameters from `WFsim.conf`. Key settings include:

```ini
TURBULENCE_REF_WAVEL     0.550000   # Reference wavelength [micron]
TURBULENCE_SEEING        0.650000   # Zenith seeing at reference wavelength [arcsec]
TURBULENCE_PROF_FILE  turbHV.prof   # Vertical profile file (Cn2, wind, etc.)
ZENITH_ANGLE             0.000000   # Zenith angle [rad]

WFsize                        512   # Output grid pixel dimension
PUPIL_SCALE              0.020000   # Pupil scale [m/pixel]
WFTIME_STEP              0.001000   # Time step between frames [s]
TIME_SPAN                0.100000   # Time span per output FITS cube [s]

SHM_SOUTPUT                     1   # 1: Output directly to shared memory stream
SHM_SPREFIX             wfsim_out   # Shared memory stream name prefix
FRESNEL_PROPAGATION             1   # 1: Enable chromatic Fresnel propagation between layers
```

---

## Directory Organization

```
milkatmturb/
├── CMakeLists.txt              # CMake plugin & standalone build configuration
├── README.md                   # Plugin documentation
└── src/
    ├── milkatmturb.h           # Master plugin header
    ├── milkatmturb.c           # MILK_MODULE registration & entry point
    ├── atmturb_compat.h        # milk framework-dev compatibility shims
    ├── AtmosphereModel/        # NRLMSISE-00, atmospheric composition, refraction
    ├── AtmosphericTurbulence/  # Turbulence phase screens, wavefront simulation, FPS
    │   ├── atmturb_mkwfs_FPS.c     # Standalone binary + CLI FPS unit for mkwfs
    │   ├── atmturb_mkhvturb_FPS.c  # Standalone binary + CLI FPS unit for mkHVturbprof
    │   └── scripts/                # runturb and mkHVturbprof shell drivers
    ├── OpticsMaterials/        # Sellmeier dispersion relations
    └── WFpropagate/            # Fresnel diffractive propagation
```

---

## License & Credits

- Original Authors: Olivier Guyon et al.
- Supported by the National Science Foundation (award #1006063).
- License: LGPL-3.0-or-later.
