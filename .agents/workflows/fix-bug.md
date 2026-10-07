---
description: Reproduce, debug, and fix a bug in milkatmturb.
---

# Fix Bug

Follow this sequence to fix bug reports:
1. **Reproduce:** Construct a test case, CLI command, or script that consistently reproduces
   the issue.
2. **Trace:** Check error outputs and logs. Use AddressSanitizer (ASan) if necessary by
   configuring CMake with `-DCMAKE_BUILD_TYPE=Debug` and `-fsanitize=address`.
3. **Fix:** Implement the bugfix using proper conventions (defensive coding, single-label cleanup,
   Allman braces, parameter alignment).
4. **Code Size Ratchet:** Verify `./scripts/check_code_size.sh` still passes with 0 violations.
5. **Verify:** Run compile and test validations:
   `bash tests/test_fps_components.sh` and `bash tests/test_wavefront_series.sh`.
