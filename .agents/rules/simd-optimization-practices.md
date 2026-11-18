---
description: Mandatory rules and constraints for SIMD intrinsics and compiler vectorization.
---

# SIMD Optimization Practices

Follow these mandatory rules when writing, auditing, or refactoring SIMD compute kernels
in `milkatmturb`.

## 1. Compiler Auto-Vectorization First; Explicit Intrinsics Only When Proven Necessary
* **Clean, Vectorizable C by Default:** Write branchless idiomatic C with `restrict` pointers,
  matched loop index types, inline polynomial math helpers, and `#pragma omp simd`.
* **Target-Compiled ISA Instantiation (`_impl.h` Pattern):** For multi-ISA dispatch (`AVX2`,
  `AVX-512`), place pure-C `#pragma omp simd` kernels in a shared `*_impl.h` header and
  instantiate thin `_avx2.c` and `_avx512.c` translation units compiled with target flags
  (e.g., `-mavx2 -mfma -fno-math-errno` and `-mavx512f ...`), as in `atmturb_simd_*` and
  `atmturb_rytov_simd_*`.
* **Strict Gate for Explicit Intrinsics (`<immintrin.h>`):** NEVER write architecture-specific
  SIMD intrinsics (`_mm256_*`, `_mm512_*`) for loops the compiler can vectorize. Introduce
  explicit intrinsics ONLY when compiler vectorization reports (`-fopt-info-vec-missed`),
  generated assembly inspection, and `benchmark_perf` measurements prove the compiler cannot
  emit efficient vector code (e.g., complex shuffle/permute patterns or bit-level operations)
  and intrinsics deliver a verified speedup.

## 2. Requirements for Multi-ISA & Intrinsics Kernels
* **ISA Separation:** Keep target-specific translation units separate (`_scalar.c`, `_avx2.c`,
  `_avx512.c`) as in `src/AtmosphericTurbulence/atmturb_simd_*`.
* **Runtime Dynamic Dispatch:** Never rely on static `-march=native` in portable builds.
  Dispatch via runtime CPU feature detection (`__builtin_cpu_supports`) at startup.
* **Mandatory Scalar Fallback:** Every multi-ISA or intrinsic kernel MUST have a verified,
  portable scalar C fallback (`_scalar.c`) to guarantee execution on non-x86 or non-AVX hosts.
* **Latency Hiding:** Ensure reduction loops expose independent accumulators or `#pragma omp simd
  reduction(...)` so the compiler or intrinsic loop saturates FMA execution ports.
* **Alignment Safety:** When explicit intrinsics are proven necessary, use unaligned loads/stores
  (`_mm256_loadu_ps` / `_mm512_loadu_ps`) unless alignment is guaranteed via `aligned_alloc()`.
