---
description: Inspect compiler assembly outputs.
---

# Inspect Machine Code

To inspect the generated assembly for hot loop vectorization in `milkatmturb`:
1. Compile the target source file to assembly:
   ```bash
   gcc -S -O3 -mavx2 -mfma -Isrc -Isrc/AtmosphericTurbulence \
       src/AtmosphericTurbulence/atmturb_simd_avx2.c -o atmturb_simd_avx2.s
   ```
2. Inspect the output `.s` file for vector instructions (e.g. `ymm` or `zmm` registers,
   `vfmadd213ps`, `vmulps`, `vaddps`).
3. Confirm that auto-vectorization did not fall back to scalar instructions in hot regions.
