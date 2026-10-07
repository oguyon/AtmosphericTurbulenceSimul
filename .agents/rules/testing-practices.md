---
description: Testing guidelines for verification of milkatmturb simulations.
---

# Testing Practices

Ensure all changes are validated for regressions and correctness before completing tasks.

## 1. Local Compile & Test
After making changes to source or CMake configuration:
1. Rebuild the project:
   ```bash
   make -C _build -j$(nproc)
   ```
2. Run code size ratchet check:
   ```bash
   ./scripts/check_code_size.sh
   ```
3. Run component tests:
   ```bash
   bash tests/test_fps_components.sh
   ```
4. Run end-to-end wavefront simulation test:
   ```bash
   bash tests/test_wavefront_series.sh
   ```

## 2. Regression Testing
When fixing a bug:
1. Create or identify a test case that triggers the bug.
2. Confirm the bug reproduces prior to code changes.
3. Confirm the bug is fully resolved after your changes, without introducing build warnings
   or breaking other tests.
