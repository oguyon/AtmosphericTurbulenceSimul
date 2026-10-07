---
name: optimize-compute-function
description: Optimization checklist for hot simulation and wavefront propagation loops.
---

# Optimize Compute Functions

Use this methodology to audit and optimize performance-critical simulation and propagation loops.

## Optimization Checklist
- [ ] **No dynamic allocations:** Verify that no `malloc`, `calloc`, or `free` calls occur inside
  hot per-frame extrusion or Fresnel propagation loops.
- [ ] **Pointer restrict:** Add `restrict` (or `__restrict__`) to pointers in the hot path to
  assist vectorization.
- [ ] **Float promotions:** Avoid double promotions; use correct float suffixes (`f`) and
  functions (`sqrtf`, `sinf`, `cosf` vs `sqrt`, `sin`, `cos`).
- [ ] **Check loop indexes:** Ensure loop variables and bounds use matching types.
- [ ] **Avoid transcendentals:** Inline small calculations and replace `pow(x, 2)` with `x * x`.
- [ ] **FFTW plan reuse:** Ensure `fftwf_plan` is created once and reused across all simulation
  steps.
- [ ] **Cache friendliness:** Process arrays in row-major order (consecutive pixel memory access).
