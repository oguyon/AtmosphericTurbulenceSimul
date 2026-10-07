---
description: Guidelines for high-performance loops, vectorization, and thread safety.
---

# Performance Practices

Simulation of multi-layer turbulence and Fresnel diffraction requires high computational throughput.
Follow these principles to optimize hot-path functions.

## 1. Compiler Attributes
- Annotate performance-critical helper functions (like bilinear interpolation or phase accumulation)
  with compiler qualifiers:
  - Use `__attribute__((hot))` on hot-path functions.
  - Use `restrict` (or `__restrict__`) on non-aliasing pointers to allow the compiler to vectorize.

## 2. Zero-Allocation on Hot Path
- Never call memory allocators (`malloc`, `calloc`, `realloc`, `free`) in per-frame simulation
  loops.
- Allocate master screens, pupil buffers, FFT working arrays, and scratch memory during setup
  and reuse them.

## 3. FFTW3 Performance & Plan Reuse
- Create FFTW plans (`fftwf_plan_dft_2d`, etc.) once during initialization with `FFTW_MEASURE`
  or `FFTW_PATIENT` if setup time permits, or `FFTW_ESTIMATE` for fast starts.
- Never destroy and re-create plans in the frame loop.
- Use FFTW single-precision (`fftwf_*`) for float arrays.

## 4. Loop Optimizations & OpenMP
- Use `#pragma omp parallel for` across independent atmospheric layers or 2D image rows.
- Minimize transcendental functions (`sin`, `cos`, `sqrt`, `pow`) in inner loops.
  - For float arrays, use `sinf`, `cosf`, `sqrtf`.
  - For integer powers, write explicit multiplications (`diff * diff`).
- Ensure unit-stride memory access (row-major order: outer loop `row` / `yy`, inner loop `col`
  / `xx`).
