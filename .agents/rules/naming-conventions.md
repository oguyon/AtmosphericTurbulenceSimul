---
description: Naming conventions for files, variables, functions, and structures.
---

# Naming Conventions

Maintain consistency with the existing `milkatmturb` codebase naming patterns.

## 1. Files & Folders
- Use lowercase `snake_case` with module prefix for C source and header files:
  - `atmturb_*.c` in `src/AtmosphericTurbulence/`
  - `atmmod_*.c` in `src/AtmosphereModel/`
  - `optmat_*.c` in `src/OpticsMaterials/`
  - `wfprop_*.c` in `src/WFpropagate/`
- Keep file names descriptive and under 600 lines.

## 2. Functions
- Public API functions declared in headers: `<module>_<verb>_<object>`:
  - Examples: `atmturb_extrude_accumulate`, `wfprop_fresnel_step`, `atmmod_compute_refraction`.
- Static helper functions: lowercase `snake_case` descriptive names.

## 3. Variables
- Local variables: short, lowercase `snake_case`.
- Loop indices:
  - Inner grid/pixel loops: use doubled letters `ii`, `jj`, `kk`.
  - Outer or layer loops: use descriptive names like `layer_idx`, `step_idx`, `frame_idx`.
  - Avoid single-character indices (`i`, `j`, `k`) in non-trivial loops to remain searchable.
- Dimension variables:
  - `xsize`, `ysize`: 2D image dimensions
  - `zsize`: 3D cube depth / number of frames or layers
  - `npupil`: pupil diameter / dimension in pixels
- Pointers: descriptive name, optionally suffixed with `_ptr` or typed arrays (`*pha`, `*amp`).

## 4. Structs & Types
- Struct names: use `UPPER_CASE` or `CamelCase` matching `milk` conventions:
  - `ATMTURB_PROFILE`, `WFPROP_FRESNEL_CONFIG`, etc.
