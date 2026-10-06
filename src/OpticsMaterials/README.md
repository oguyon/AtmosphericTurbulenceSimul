# OpticsMaterials

**Architectural Layer**: Level 0 (Foundation / Physical Property Models)

## Purpose

The `OpticsMaterials` module provides refractive index evaluations and optical phase delay calculations for gases, optical glasses, semiconductors, and photoresists across visible and near-infrared wavelengths.

## Components

- `OpticsMaterials.c` / `OpticsMaterials.h`: Public API orchestrator for material name/code resolution and refractive index lookups.
- `optmat_types.h`: Internal data structures and function prototypes.
- `optmat_solids.c`: Sellmeier dispersion equations for solids (Fused silica $\text{SiO}_2$, silicon $\text{Si}$, calcium fluoride $\text{CaF}_2$, mirror).
- `optmat_gases.c`: Atmospheric gas dispersion relations (air, $\text{N}_2$, $\text{O}_2$, $\text{Ar}$, $\text{He}$, $\text{H}_2$, $\text{H}_2\text{O}$ vapor, $\text{CO}_2$, $\text{Ne}$, atomic oxygen $\text{O}$, vacuum).
- `optmat_pmgi.c`: High-resolution tabulated dispersion data and interpolation for PMGI resist ($400\,\text{nm} - 1000\,\text{nm}$).
- `optmat_pmma.c`: High-resolution tabulated dispersion data and interpolation for PMMA resist ($400\,\text{nm} - 1000\,\text{nm}$).

## Public Headers

- `OpticsMaterials.h`:
  - `int init_OpticsMaterials(void)`
  - `int OPTICSMATERIALS_code(char *name)`
  - `char* OPTICSMATERIALS_name(int code)`
  - `double OPTICSMATERIALS_n(int material, double lambda)`
  - `double OPTICSMATERIALS_pha_lambda(int material, double z, double lambda)`

## Dependencies

- C Standard Library (`math.h`, `stdio.h`, `stdlib.h`, `string.h`)
- No dependencies on higher-level simulation or CLI modules.
