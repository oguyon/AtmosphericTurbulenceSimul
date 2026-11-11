# Atmospheric Turbulence Simulation (`milkatmturb`)

## Overview

`milkatmturb` provides multi-layer atmospheric turbulence synthesis, chromatic ray tracing,
diffractive wavefront propagation, real-time shared memory streaming, and closed-loop adaptive
optics (AO) evaluation for the `milk` (`feat/bash-engine`) framework.

Key physics and compute capabilities include:

- **Phase Screen Synthesis**: Power-law Kolmogorov and von Karman Fourier phase screens with outer
  scale $L_0$, inner scale $l_0$, analytic Lane/Johansson-Gavel low-frequency subharmonics, and
  variance-preserving quadrature ($\cos\theta / \sin\theta$) rolling screen evolution (boiling).
- **Chromatic Atmospheric Refraction**: Curved-ray path integration through 10000-layer standard or
  site-measured $(P, T, \mathrm{RH})$ vertical profiles using Edlen/Sellmeier equations or HITRAN
  Refractive Index & Absorption (RIA) tables.
- **Wavefront Propagation Modes (`FRESNEL_PROPAGATION`)**:
  - `0` (*Geometric*): Line-of-sight optical path difference integration (unit amplitude).
  - `1` (*Split-Step Fresnel*): Angular-spectrum diffractive propagation across discrete layers or
    scintillation-weighted super-layers with $O(\Delta z^2)$ piecewise-linear hat nodes.
  - `2` (*First-Order Rytov*): Weak-scintillation Fourier quadratic-phase filtering ($\cos/\sin$
    phase and log-amplitude kernels) with Moisan periodic-plus-smooth boundary decomposition.
- **Hardware Acceleration**: Runtime-dispatched x86 SIMD (`AVX-512`, `AVX2+FMA`, scalar fallback),
  OpenMP multi-threading, and optional NVIDIA CUDA + cuFFT GPU pipelines with CUDA Graph capture
  and pinned-memory ImageStreamIO DMA.

---

## Interactive `milk` Commands (`atmturb.*`)

When `libmilkatmturb.so` is loaded in `milk-cli` or `milk.bash`, commands are registered under the
`atmturb` module namespace:

### FPS Compute Units (`MILK_RUNNABLE_FPS`)

| Command | Standalone Binary | Description |
|---|---|---|
| `atmturb.mkwfs` | `milk-fpsexec-atmturb-mkwfs` | 3D wavefront series or 2D SHM stream |
| `atmturb.mkhvturb` | `milk-fpsexec-atmturb-mkhvturb` | Hufnagel-Valley $C_n^2(h)$ profile |
| `atmturb.mkmastert` | `milk-fpsexec-atmturb-mkmastert` | Seeded von Karman master screens |
| `atmturb.mkvonkarman` | `milk-fpsexec-atmturb-mkvonkarman` | 1D von Karman wind time series |
| `atmturb.aoloop` | `milk-fpsexec-atmturb-aoloop` | Closed-loop AO integrator & PSF |

Inspect parameters and output variables for any runnable with `-h`:

```bash
milk-cli -c "atmturb.mkwfs -h"
milk-fpsexec-atmturb-mkwfs --help
```

### Immediate CLI Utilities (`MILK_RUNNABLE_CLI`)

| Command | Description |
|---|---|
| `atmturb.mkatm` | Compute standard atmospheric profile and write `atmosphere.dat` |
| `atmturb.measr0` | Estimate Fried parameter $r_0$ `[px]` from 2D structure function |
| `atmturb.linpred_build` | Train autoregressive full-aperture linear predictor via SVD |
| `atmturb.linpred_apply` | Evaluate linear predictor on phase cube (`outwfc` / `outres`) |
| `atmturb.contract` | Spatially bin complex `(amp, pha)` wavefront cubes |

---

## Real-Time Streaming vs. Batch Cube Synthesis

`atmturb.mkwfs` (`milk-fpsexec-atmturb-mkwfs`) supports two output modes controlled by `stream2d`:

1. **Batch 3D Cube Mode (`stream2d=OFF`)**: Synthesizes `nbstep` frames into 3D `(size, size,
   nbstep)` amplitude (`<prefix>amp`) and phase (`<prefix>pha`) images and writes FITS cubes when
   `write_fits=ON`.
2. **Real-Time 2D Stream Mode (`stream2d=ON`)**: Continuously renders 2D `(size, size)` wavefront
   frames into persistent ImageStreamIO shared memory streams (`<prefix>pha`, `<prefix>amp`) with
   semaphore posting (`ImageStreamIO_sempost`) and optional frame-rate pacing (`dt_us`).

---

## Helper Scripts (`src/AtmosphericTurbulence/scripts/`)

Legacy and convenience workflow scripts are located in `src/AtmosphericTurbulence/scripts/` and
installed to `${CMAKE_INSTALL_PREFIX}/bin`:

- `mkHVturbprof`: Generates a multi-layer Hufnagel-Valley turbulence profile (`turb.prof`).
- `runturb`: End-to-end shell driver configuring site atmosphere, master screens, and wavefront
  propagation at a target wavelength `[um]` (run `runturb -h` for options).
- `syncscripts`: Symlinks helper scripts into a local working directory.

---

## Refractive Index Models

By default, air refractivity $(n - 1)$ is evaluated analytically from pressure, temperature, water
vapor partial pressure, and $\mathrm{CO}_2$ concentration via Edlen/Sellmeier dispersion formulas
(`AtmosphereModel_stdAtmModel_N`). High-resolution molecular line absorption and anomalous
dispersion tables (HITRAN RIA files) can also be loaded via `AtmosphereModel_Load_RIA_table`.
