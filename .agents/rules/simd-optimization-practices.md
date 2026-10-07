---
description: Mandatory rules and constraints for SIMD intrinsics and compiler vectorization.
---

# SIMD Optimization Practices

Follow these mandatory rules when writing, auditing, or refactoring SIMD compute kernels
in `milkatmturb`.

## 1. When to Use Explicit SIMD vs. Compiler Auto-Vectorization
* **Clean C by Default:** Write idiomatic C with `restrict` pointers and `#pragma omp simd`
  for element-wise array operations.
* **Explicit SIMD Kernels for Bottlenecks:**
  - Bilinear phase screen extrusion (`atmturb_extrude_accumulate`)
  - Float array scaling (`atmturb_scale_float_array`)
  - Phase and amplitude array initialization (`atmturb_init_phase_amp`)
  - Complex Fresnel propagation kernel multiplications

## 2. Requirements for Explicit SIMD Kernels
* **ISA Separation:** Keep architecture-specific kernels in separate files
  (`_scalar.c`, `_avx2.c`, `_avx512.c`) as in `src/AtmosphericTurbulence/atmturb_simd_*`.
* **Runtime Dynamic Dispatch:** Never rely on static `-march=native`. Dispatch via CPUID
  detection (`atmturb_simd_dispatch.c`) at startup.
* **Mandatory Scalar Fallback:** Every handwritten SIMD function MUST have a verified,
  portable scalar C fallback to guarantee execution on non-x86 or non-AVX machines.
* **Latency Hiding:** Use multiple independent vector accumulators in reduction or sum loops
  to saturate processor execution ports.
* **Alignment Safety:** Use `_mm256_loadu_ps` / `_mm512_loadu_ps` unless memory alignment
  is explicitly guaranteed.
