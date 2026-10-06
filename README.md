# milkatmturb - Atmospheric Turbulence Simulation Plugin for milk

`milkatmturb` is a high-performance atmospheric turbulence and wavefront propagation plugin for the [milk](https://github.com/milk-org/milk) framework (`framework-dev` branch). It simulates the dynamic optical effects of Earth's atmosphere on astronomical wavefronts across visible and near-infrared wavelengths.

---

## Key Features

- **Multi-layer Atmospheric Modeling**:
  - Chromatic diffractive Fresnel propagation between atmospheric layers.
  - Realistic atmospheric vertical composition and dispersion profiles (Sellmeier equations and NRLMSISE-00 model).
  - Hufnagel-Valley (HV) and customized $C_n^2$ profiles with altitude-dependent wind speed and direction.
- **Turbulence Realizations**:
  - von Kármán and Kolmogorov master phase screens with configurable inner ($l_0$) and outer ($L_0$) scales.
  - Multi-layer phase screen extrusion with geometric shear and atmospheric refraction.
- **Wavefront Prediction & Analysis**:
  - Linear predictor computation and testing on turbulent wavefront sequences.
  - Contrast evaluation and PSF computation over exposure intervals.
- **Real-Time Integration with milk**:
  - Low-latency shared memory stream outputs compatible with `ImageStreamIO`.
  - Clock synchronization, frame rates, and time scaling factors for hardware-in-the-loop and real-time AO simulations.
  - Standalone FPS V2 compute units with processinfo telemetry, CLI commands, and standalone executables.

---

## Directory & Module Structure

The plugin code is organized cleanly into modular subdirectories:

```
├── CMakeLists.txt              # CMake plugin build & standalone executable targets
├── README.md                   # This documentation
└── src
    ├── milkatmturb.h           # Master plugin header
    ├── milkatmturb.c           # MILK_MODULE registration & entry point
    ├── atmturb_compat.h        # milk framework-dev compatibility definitions
    ├── AtmosphereModel/        # Atmospheric density, composition, and refraction models
    │   ├── AtmosphereModel.c
    │   └── AtmosphereModel.h
    ├── AtmosphericTurbulence/  # Turbulence screens, wavefront series, FPS compute units
    │   ├── AtmosphericTurbulence.c
    │   ├── AtmosphericTurbulence.h
    │   ├── atmturb_mkwfs_FPS.c     # Standalone / CLI FPS compute unit for wavefront series
    │   ├── atmturb_mkhvturb_FPS.c  # Standalone / CLI FPS compute unit for HV profiles
    │   └── scripts/
    │       ├── runturb             # End-to-end turbulence simulation runner
    │       └── mkHVturbprof        # Script to create Hufnagel-Valley turbulence profile
    ├── OpticsMaterials/        # Sellmeier optical dispersion equations
    │   ├── OpticsMaterials.c
    │   └── OpticsMaterials.h
    └── WFpropagate/            # Fresnel diffractive propagation routines
        ├── WFpropagate.c
        └── WFpropagate.h
```

---

## Building and Installation

`milkatmturb` is built as an out-of-tree or in-tree plugin for `milk` (`framework-dev` branch).

### Option A: In-Tree Build via milk plugins directory (Recommended)

1. Clone or symlink this repository into `milk/plugins/milkatmturb`:
   ```bash
   ln -s /path/to/AtmosphericTurbulenceSimul /path/to/milk/plugins/milkatmturb
   ```

2. Configure and build `milk`:
   ```bash
   cd /path/to/milk
   cmake -B _build -DCMAKE_BUILD_TYPE=Release
   cmake --build _build --target milkatmturb milk-fpsexec-atmturb-mkwfs milk-fpsexec-atmturb-mkhvturb -j
   ```

3. Install:
   ```bash
   cmake --install _build
   ```

### Dependencies
- `milk` (`framework-dev` branch, including `CLIcore`, `ImageStreamIO`, `COREMOD_*`, `milkfft`, `milklinalgebra`)
- `OpenMP` (for multi-threaded screen extrusion and Fresnel propagation)
- `GSL` (GNU Scientific Library, for atmospheric modeling)
- `FFTW3` / `FFTW3F`
- `CFITSIO`

---

## Usage

### 1. Standalone FPS Executables (FPS V2 Architecture)

`milkatmturb` provides standalone executables following the milk FPS V2 pattern:

- **Generate Hufnagel-Valley Profile**:
  ```bash
  milk-fpsexec-atmturb-mkhvturb exec 20.0 0.15 4200 20 turbHV.prof
  ```
  Options and parameters:
  - `wspeed`: High-altitude wind speed [m/s] (default: `20.0`)
  - `r0`: Fried parameter at $0.55\,\mu\text{m}$ [m] (default: `0.15`)
  - `sitealt`: Site altitude [m] (default: `4200.0`)
  - `nblayers`: Number of discretized layers (default: `20`)
  - `outfile`: Output file path (default: `turbHV.prof`)

- **Generate Wavefront Series**:
  ```bash
  milk-fpsexec-atmturb-mkwfs exec 1650.0 0
  ```
  Options and parameters:
  - `slambda`: Science wavelength in $\mu\text{m}$ or $\text{nm}$ (e.g. `1650.0` or `1.65`)
  - `precision`: `0` for single precision, `1` for double precision
  - `conffile`: Configuration file (default: `WFsim.conf`)

Both executables support full FPS lifecycle control (`fpsinit`, `confstart`, `runstart`, `runstop`, `-tmux`, `-procinfo`, `--help`).

### 2. Interactive milk CLI

When launching `milk-cli`, `milkatmturb` commands are accessible under the `atmturb` module namespace:

```
milk-cli > m? milkatmturb
---- MODULE milkatmturb COMMANDS ---------
   atmturb.mkatmospheremodel  make Earth atmosphere model
   atmturb.fresnelpw          Fresnel propagate wavefront
   atmturb.mkwfs              make wavefront series
   atmturb.mkvonKarmanWind    make vonKarman wind model
   atmturb.mkmastert          make 2 master phase screens
   atmturb.mkHVturbprof       make Hufnager-Valley turbulence profile
   atmturb.mkwfs_fps          Generate atmospheric turbulence wavefront series (FPS)
   atmturb.mkhvturb_fps       Generate Hufnagel-Valley turbulence profile (FPS)
   ...
```

Run directly:
```bash
milk-cli -c "atmturb.mkHVturbprof 20.0 0.15 4200 20 turbHV.prof"
milk-cli -c "atmturb.mkwfs 1650.0 0"
```

### 3. Simulation Scripts

- `mkHVturbprof <windspeed> <r0> <sitealt> <nblayers>`: Generates `turbHV.prof` using the standalone FPS executable or `milk-cli`.
- `runturb [-hNdT] <slambdaum>`: Automated driver script that configures `WFsim.conf`, generates master turbulence screens, and executes the simulation.

---

## Credits & License

- Original implementation: Olivier Guyon et al.
- Developed with support from the National Science Foundation (award #1006063).
- License: LGPL-3.0-or-later.
