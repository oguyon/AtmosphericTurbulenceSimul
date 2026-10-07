---
name: diagnose-build-failure
description: Troubleshooting CMake and compilation failures in milkatmturb.
---

# Diagnose Build Failures

A guide for troubleshooting build and link errors in `milkatmturb`.

## 1. CMake Configuration Failures
- **Missing packages:** Look for pkg-config or CMake errors:
  - "fftw3 or fftw3f not found" -> Install `libfftw3-dev`.
  - "OpenMP not found" -> Install `libomp-dev`.
  - "CLIcore not found" -> Ensure `milk` is installed, or set `MILK_SOURCE_DIR` / `MILK_ROOT`.

## 2. Compiler Errors
- **Undefined references / symbols:** Check if headers are missing or if the source file was not
  added to `SOURCEFILES` in `CMakeLists.txt`.
- **Implicit declaration warnings:** Occur when a function is called without a matching header
  include. Fix by explicitly adding the `#include` of the header defining the function.
- **Size ratchet failures:** If compilation passes but `./scripts/check_code_size.sh` fails,
  decompose oversized functions (<60 lines) or files (<600 lines).

## 3. Linker Errors
- **Missing library linkage:** Ensure that the target library and executables link against
  `${FFTW_LIBRARIES}`, `${FFTWF_LIBRARIES}`, `OpenMP::OpenMP_C`, and `m`.
