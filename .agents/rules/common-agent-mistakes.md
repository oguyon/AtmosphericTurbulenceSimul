---
trigger: always_on
---

# Common Agent Mistakes

Checklist of pitfalls that AI agents frequently hit when generating code for `milkatmturb`.
Check this list before finalizing any generated code.

## CMake & Build System
1. **Forgetting to add `.c` to `SOURCEFILES`.**
   If you add a source file used by the plugin library, it must be added to `SOURCEFILES` in
   the root `CMakeLists.txt`. If it is a public header, add it to `INCLUDEFILES`.
2. **Not linking required libraries.**
   When adding mathematical or framework capabilities, ensure required libraries are linked:
   `OpenMP::OpenMP_C`, `GSL::gsl`, `${FFTW_LIBRARIES}`, `${FFTWF_LIBRARIES}`, `ImageStreamIO`,
   `CLIcore`, `m`.

## Simulation Loops & Performance
3. **Allocation inside simulation loops.**
   Never call `malloc()`, `calloc()`, or `realloc()` inside hot wavefront extrusion, screen
   generation, or Fresnel propagation loops. Pre-allocate working arrays in setup functions.
4. **Recreating FFTW plans per frame.**
   Creating an FFTW plan (`fftwf_plan_dft_2d`) is computationally expensive and introduces
   locking. Always create plans during initialization and reuse them across frames.
5. **Leaving test SHM files.**
   When creating test `ImageStreamIO` streams in tests/debug, clean them up from `/dev/shm/`
   immediately afterwards.

## Compiler & Math Correctness
6. **Implicit double promotions.**
   Ensure that float math uses float variants (`sqrtf`, `sinf`, `cosf`) and float literals (`0.5f`)
   when working with single-precision floats, or standard double variants when using double.
7. **Type mismatches in loop indices.**
   Always match the loop index type to the bound variable type (e.g., `for (long ii = 0; ii <
     n; ii++)`
   when `n` is `long`) to prevent breaking compiler SIMD auto-vectorization.
8. **Implicit header includes.**
   Every `.c` file must include exactly the headers it uses. Do not rely on header side-effects.
9. **Lines > 100 characters.**
   Limit line length in C source code, scripts, and documentation files to 100 characters.

## Memory Alignment & Lifecycle
10. **`aligned_alloc(alignment, size)` sizing rule.**
    POSIX and ASan mandate that `size` MUST be an integral multiple of `alignment`
    (`size % alignment == 0`). Non-multiples abort under ASan with
    `invalid-aligned-alloc-alignment`. Pad `size` to a multiple of `alignment`.
11. **Dangling pointers to stack buffers.**
    Never assign local stack buffers (`char buf[...]`) to struct pointer fields or return them.
    Under `-O2`/`-O3` optimization, stack frames are popped and overwritten, causing memory
      corruption.

## Pre-Merge Verification
12. **Merging before running ratchet and tests.**
    Always run `./scripts/check_code_size.sh` and the automated test scripts before considering
    work complete.
