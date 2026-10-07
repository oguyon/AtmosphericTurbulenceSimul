## Summary of Changes
<!-- Provide a clear, concise overview of what changed and why. -->

## Changes Checklist
- [ ] Code compiles cleanly with zero warnings (`-Wall -Wextra`)
- [ ] Code size ratchet passes (`./scripts/check_code_size.sh`)
- [ ] Automated tests pass (`bash tests/test_fps_components.sh`,
      `bash tests/test_wavefront_series.sh`)
- [ ] Adheres to C code style guide (Allman braces, <= 100 character lines,
      column-aligned parameters)
- [ ] No allocations (`malloc`, `calloc`) inside hot simulation loops (phase extrusion, Fresnel)
- [ ] No resource/memory leaks (verified with ASan/UBSan or Valgrind)
- [ ] New source files reflected in `CMakeLists.txt` and corresponding module `README.md`

## Testing Performed
<!-- Describe how these changes were tested (commands executed, scripts run, output verified). -->
```bash
# Example:
make -C _build -j$(nproc)
./scripts/check_code_size.sh
bash tests/test_fps_components.sh
bash tests/test_wavefront_series.sh
```

---
<!--
If agentic AI tools were used, disclose below per project guidelines:
Implemented by <model name>. Reviewed and signed off by O. Guyon.
Followed by a concise technical summary of what task the model performed.
-->
Implemented by . Reviewed and signed off by O. Guyon.
