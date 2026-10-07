---
description: Thread safety and concurrency patterns for milkatmturb.
---

# Concurrency Practices

`milkatmturb` supports multi-threading via OpenMP and real-time streaming via `ImageStreamIO`
shared memory. Follow these patterns to avoid race conditions.

## OpenMP Parallelization
- Ensure loops parallelized via `#pragma omp parallel for` do not have race conditions on
  shared state.
- Keep loop index variables private to each thread (declaring them inside the `for` statement
  achieves this automatically in C99/C11).
- When accumulating statistics or metrics across threads, use `#pragma omp atomic` or reduction
  clauses:
  ```c
  #ifdef _OPENMP
  #pragma omp atomic
  #endif
  total_variance += local_var;
  ```
- Guard multi-layer phase screen extrusion and Fresnel propagation so each layer or independent
  pixel operates on private memory or disjoint memory slices.

## Shared Memory Semaphore Protocol
- When reading/writing streams via `ImageStreamIO`, use the library's semaphore functions
  (`ImageStreamIO_semwait`, `ImageStreamIO_sempost`) for synchronization.
- Post to all listening processes when updating a stream:
  ```c
  ImageStreamIO_sempost(&img, -1);
  ```
- Always check return values of `ImageStreamIO` calls for failure.

## Volatile Keyword
- Use `volatile sig_atomic_t` for signal handling flags (e.g., `stop_requested`):
  ```c
  extern volatile sig_atomic_t stop_requested;
  ```
- Do not use `volatile` on data arrays or pointers, as it disables compiler optimization and
  SIMD vectorization entirely.
