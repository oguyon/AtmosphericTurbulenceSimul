---
name: simd-optimization
description: Methodology for writing, auditing, and benchmarking SIMD compute kernels in
  milkatmturb.
---

# SIMD Optimization Guide

`milkatmturb` performs high-resolution atmospheric turbulence synthesis and wavefront propagation.
Hot compute loops must saturate CPU vector execution pipelines and cache bandwidth.

This skill documents the SIMD methodology in `milkatmturb`: **rely on compiler auto-vectorization
(`#pragma omp simd`) first, and introduce explicit SIMD intrinsics ONLY when proven necessary**.

---

## 1. Core Rule: Auto-Vectorization First, Intrinsics Last

1. **Pure-C `#pragma omp simd` First:**
   - Express hot loops in branchless C using `restrict` pointers, matching index/bound integer
     types, unit-stride access, and inline branchless polynomial approximations where `<math.h>`
     calls would inhibit vectorization (see `atmturb_rytov_fast_exp_clamped` in
     `src/AtmosphericTurbulence/atmturb_rytov_simd_avx2.c`).
2. **Per-ISA Target Instantiation (`*_impl.h`):**
   - Place shared `#pragma omp simd` C implementations in an internal `*_impl.h` header (such as
     `atmturb_simd_extrude_impl.h` and `atmturb_simd_ops_impl.h`) and include it from thin
     `_avx2.c` and `_avx512.c` translation units compiled with target-specific flags.
3. **Strict Gate for Explicit Intrinsics (`<immintrin.h>`):**
   - **Never** write `_mm256_*` or `_mm512_*` intrinsics if GCC/Clang already vectorizes the C
     loop cleanly.
   - Use explicit intrinsics **only** after verifying via `-fopt-info-vec-missed` (or Clang
     `-Rpass-missed=loop-vectorize`), assembly inspection (`objdump -d`), and `benchmark_perf`
     that the compiler cannot generate efficient vector instructions and intrinsics provide a
     measurable speedup.

---

## 2. Architecture & File Layout

Multi-ISA dispatch in `milkatmturb` follows the ISA Strategy Separation pattern:
- **`atmturb_*_simd.h`**: Public API declarations and function prototypes.
- **`atmturb_*_simd_impl.h`**: Shared pure-C `#pragma omp simd` kernel templates (preferred).
- **`atmturb_*_simd_scalar.c`**: Portable C scalar fallback (always compiled without AVX flags).
- **`atmturb_*_simd_avx2.c`**: Target TU compiled with `-mavx2 -mfma -fno-math-errno`.
- **`atmturb_*_simd_avx512.c`**: Target TU compiled with `-mavx512f -mavx512dq -mavx512bw
  -mavx512vl -mfma -fno-math-errno`.
- **`atmturb_*_simd_dispatch.c`**: Runtime CPU feature detection via `__builtin_cpu_supports`
  and function pointer binding (honoring `ATMTURB_SIMD=SCALAR|AVX2|AVX512`).

---

## 3. Verification & Kernel Guidelines

- **Verify Compiler Output First:** Inspect `-fopt-info-vec-optimized` / `-fopt-info-vec-missed`
  before reaching for `<immintrin.h>`. Often removing a branch, splitting a toroidal wrap seam,
  or adding `restrict` unlocks clean compiler `vfmaddps` codegen without intrinsics.
- **Latency Hiding:** Expose independent accumulators or `#pragma omp simd reduction(...)` in
  reduction loops to saturate FMA execution ports.
- **Unaligned Memory Loads (if intrinsics are required):** Use `_mm256_loadu_ps` and
  `_mm512_loadu_ps` unless allocations are explicitly aligned via `aligned_alloc()`.
- **Mandatory Scalar Fallback:** Every dispatched kernel must have a matching scalar
  implementation in `_scalar.c` to guarantee portability across all architectures.
- **Benchmark Verification:** Validate vectorized speedups against `benchmark_perf` in `tests/`.
