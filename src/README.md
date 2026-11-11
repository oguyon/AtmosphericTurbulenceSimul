# src (milkatmturb Top-Level Plugin Layer)

**Architectural Layer**: Level 3 (Plugin Registration, Umbrella Header & Compatibility Shims)

## Purpose

The top-level `src/` directory provides the `milk` plugin registration table (`libmilkatmturb.so`),
the umbrella public header (`milkatmturb.h`), C23/C11 language feature shims, and internal
convenience wrappers over `milk` (`feat/bash-engine`) shared-memory image APIs.

## Submodules by Layer

- **`AtmosphereModel/`** (*Level 1*): Empirical atmospheric density/temperature profiles
  (NRLMSISE-00), 14-species Lorentz-Lorenz optical dispersion, water-vapor saturation (IAPWS-95),
  and spherical-shell refraction ray tracing. See `src/AtmosphereModel/README.md`.
- **`AtmosphericTurbulence/`** (*Level 2*): Kolmogorov / von Karman Fourier phase screen synthesis,
  Hufnagel-Valley profiles, sub-pixel multi-layer wind extrusion, Fresnel & Rytov wave propagation,
  CPU SIMD (AVX2 / AVX-512) and CUDA GPU acceleration, closed-loop AO simulation, linear wavefront
  prediction, and `MILK_RUNNABLE` CLI/FPS components. See `src/AtmosphericTurbulence/README.md`.

## Top-Level Files

- `milkatmturb.c`: Declares `MODULE_DEPS(...)` and exports the `MILK_MODULE_RUNNABLES` registry
  exposing all 10 `atmturb.*` runnable units (`mkwfs`, `mkhvturb`, `mkmastert`, `mkvonkarman`,
  `aoloop`, `mkatm`, `measr0`, `linpred.build`, `linpred.apply`, `contract`).
- `milkatmturb.h`: Public umbrella header re-exporting `atmturb_c23_compat.h`,
  `AtmosphereModel/AtmosphereModel.h`, `AtmosphericTurbulence/AtmosphericTurbulence.h`,
  `OpticsMaterials/OpticsMaterials.h`, and `WFpropagate/WFpropagate.h`.
- `atmturb_c23_compat.h`: Installed C23 portability header providing `bool`, `nullptr`,
  `ATMTURB_CONSTEXPR`, `ATMTURB_NODISCARD`, `ATMTURB_UNUSED`, `ATMTURB_FALLTHROUGH`,
  `ATMTURB_STATIC_ASSERT`, and `<stdbit.h>` bit-manipulation fallbacks for C11 compilers.
- `atmturb_compat.h`: Internal (uninstalled) convenience header mapping concise return-by-value
  image creation, lookup, FITS I/O, and SVD pseudoinverse helpers onto `milk_data` and core module
  functions.

## Standalone Executables

Each FPS-enabled `MILK_RUNNABLE` also builds as a standalone `FPS_MAIN_STANDALONE` executable:
- `milk-fpsexec-atmturb-mkwfs`
- `milk-fpsexec-atmturb-mkhvturb`
- `milk-fpsexec-atmturb-mkmastert`
- `milk-fpsexec-atmturb-mkvonkarman`
- `milk-fpsexec-atmturb-aoloop`

## Permitted Dependencies

- Submodules `AtmosphereModel` (Level 1) and `AtmosphericTurbulence` (Level 2)
- `milk` core and Level-1 shared modules (`milk_runnable.h`, `milkCOREMODmemory`,
  `milkCOREMODarith`, `milkCOREMODiofits`, `milkCOREMODtools`, `milkfft`, `milkOpticsMaterials`,
  `milkWFpropagate`, `milklinoptimtools`)
