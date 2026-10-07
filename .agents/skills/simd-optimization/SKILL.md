---
name: simd-optimization
description: Methodology for writing, auditing, and benchmarking SIMD compute kernels in
  milkatmturb.
---

# SIMD Optimization Guide

`milkatmturb` performs high-resolution atmospheric turbulence synthesis and wavefront propagation.
Hot compute loops must saturate CPU vector execution pipelines and cache bandwidth.

This skill documents the SIMD architecture in `milkatmturb`, identifying where explicit intrinsics
are essential for maximum throughput and how kernels are organized across ISAs.

---

## 1. High-Value SIMD Opportunities in milkatmturb

1. **Bilinear Phase Screen Extrusion (`atmturb_extrude_accumulate`):**
   * *Mechanism:* Bilinear interpolation weights ($w_{00}, w_{01}, w_{10}, w_{11}$) broadcast
     across vector registers, computing 8 or 16 pupil pixels concurrently from master screens.
   * *Files:* `src/AtmosphericTurbulence/atmturb_simd_avx2.c`, `atmturb_simd_avx512.c`.
2. **Float Array Scaling (`atmturb_scale_float_array`):**
   * *Mechanism:* Scalar multiplier broadcast into vector register, executing 8 or 16 multiplies
     per cycle.
   * *Files:* `atmturb_simd_avx2.c`, `atmturb_simd_avx512.c`.
3. **Phase & Amplitude Initialization (`atmturb_init_phase_amp`):**
   * *Mechanism:* Vector zeroing for phase and vector unit fill for amplitude arrays.
   * *Files:* `atmturb_simd_avx2.c`, `atmturb_simd_avx512.c`.
4. **Fresnel Propagation Kernel Products:**
   * *Mechanism:* Complex multiplication of pupil fields by quadratic Fresnel phase factors.
   * *Files:* `src/WFpropagate/wfprop_fresnel.c`.

---

## 2. Architecture & File Layout

Explicit SIMD in `milkatmturb` follows the ISA Strategy Separation pattern:
- **`atmturb_simd.h`**: Public API declarations and function prototypes.
- **`atmturb_simd_scalar.c`**: Portable C scalar fallback (always compiled).
- **`atmturb_simd_avx2.c`**: AVX2 + FMA vectorized implementation (compiled with `-mavx2 -mfma`).
- **`atmturb_simd_avx512.c`**: AVX-512 implementation (compiled with `-mavx512f -mavx512dq -mfma`).
- **`atmturb_simd_dispatch.c`**: Runtime CPU feature detection via `__builtin_cpu_supports` and
  function pointer binding.

---

## 3. Kernel Guidelines

- **Latency Hiding:** Modern x86 FMA units require multiple independent accumulators to prevent
  pipeline stalls. Unroll reduction loops with 4 to 8 accumulators.
- **Unaligned Memory Loads:** Use `_mm256_loadu_ps` and `_mm512_loadu_ps` unless allocations
  are explicitly aligned to 32 or 64 bytes with `posix_memalign()`.
- **Mandatory Scalar Fallback:** Every handwritten SIMD function must have a matching scalar
  implementation to guarantee portability on non-x86 hardware.
- **Benchmark Verification:** Validate vectorized speedups against `benchmark_perf.c` in `tests/`.
