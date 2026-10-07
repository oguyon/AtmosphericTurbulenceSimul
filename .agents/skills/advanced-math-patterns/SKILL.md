---
name: advanced-math-patterns
description: Reference for implementing high-performance mathematical operations in milkatmturb.
---

# Advanced Math Patterns

Numerical calculations in atmospheric turbulence, phase unwrapping, and Fresnel propagation
can easily become computational bottlenecks. Follow these standard patterns for high-
  performance math.

## 1. Array and Grid Vectorization (SIMD)
When writing operations across 2D pupil arrays or 3D cubes, guide the compiler to vectorize:
1. **Restrict Pointers:** Always use `restrict` on local pointers to indicate that memory
   regions do not overlap.
2. **Type Matching:** Loop index types MUST match bound variables (e.g. `long ii = 0; ii < npupil`).
3. **OpenMP + SIMD:** Use OpenMP pragmas to parallelize loops:
   ```c
   #pragma omp parallel for
   for (long ii = 0; ii < size; ii++)
   {
       out[ii] = in[ii] * factor;
   }
   ```

## 2. Avoiding Transcendentals in Hot Loops
Function calls like `sqrt()`, `sin()`, and `pow()` inside loops degrade performance:
- **Integer Powers:** Never use `pow(x, 2)`. Use `x * x` instead.
- **Float vs Double:** If processing floats, use `sqrtf()`, `sinf()`, `cosf()`. If doubles,
  use `sqrt()`, `sin()`, `cos()`. Do not mix types to avoid costly promotions.
- **Trigonometric Evaluation:** When both sine and cosine of an angle are needed (e.g. complex phase
  factor $e^{i\phi} = \cos\phi + i\sin\phi$), use `sincosf()` (GNU extension) or explicit
  pre-computed phasor lookups.

## 3. FFTW3 / FFTW3F Optimization
- Pre-allocate FFTW plans once during setup: `fftwf_plan_dft_2d(...)`.
- Use `fftwf_*` for single-precision arrays and `fftw_*` for double-precision arrays.
- Remember that FFTW does not normalize unnormalized DFTs: multiply inverse DFTs by $1.0f / N$.
