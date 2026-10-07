---
description: Perform compile and verification testing on milkatmturb.
---

# Compile Test

To verify the build:
1. Navigate to the build directory (create if it does not exist):
   ```bash
   mkdir -p _build && cd _build
   ```
2. Configure using CMake (if not already configured):
   ```bash
   cmake .. -DCMAKE_BUILD_TYPE=Release
   ```
3. Compile with multiple jobs:
   ```bash
   make -j$(nproc)
   ```
4. Verify there are no warnings or errors.
5. Check code size constraints:
   ```bash
   ./scripts/check_code_size.sh
   ```
6. Run automated test scripts:
   ```bash
   bash tests/test_fps_components.sh
   bash tests/test_wavefront_series.sh
   ```
