---
name: refactor-c-source
description: Safely split and reorganize large C source files and long functions into clean modules.
---

# Refactoring C Source Files & Long Functions

This skill provides step-by-step guidance for safely decomposing oversized C source files
(> 600 lines) and long functions (> 60 lines) into modular, testable, and maintainable units
without introducing regressions or performance penalties.

---

## 1. Targeted Metrics & Code Hygiene

When refactoring, aim for the following targets:
- **Source Files**: $\le 600$ lines (hard limit $1000$ lines).
- **Function Bodies**: $\le 60$ lines (hard limit $150$ lines).
- **`main()` Functions**: $\le 40$ lines (hard limit $80$ lines).
- **Brace Style**: Strictly Allman (opening brace `{` on its own line).
- **Documentation**: Kernel-Doc comment above every non-trivial function with a clear,
  single-sentence summary on the first line.
- **Enforcement**: Run `./scripts/check_code_size.sh` after any change.

---

## 2. Function Decomposition Patterns

### Pattern A: Pipeline / Phase Extraction (Context Struct)
Long computational pipelines that interweave preparation, multi-layer propagation, and output
should be partitioned into sequential phase functions sharing a context structure:

```c
/* Instead of one 400-line monolithic runner: */
int atmturb_wfs_run(struct atmturb_sim_config *cfg)
{
    struct atmturb_sim_ctx ctx;

    if (atmturb_sim_init(&ctx, cfg) != 0)
    {
        return -1;
    }

    if (atmturb_sim_extrude_screens(&ctx) != 0 ||
        atmturb_sim_propagate_layers(&ctx) != 0 ||
        atmturb_sim_save_telemetry(&ctx) != 0)
    {
        atmturb_sim_cleanup(&ctx);
        return -1;
    }

    atmturb_sim_cleanup(&ctx);
    return 0;
}
```

### Pattern B: ISA Strategy Separation (SIMD Kernels)
When a numerical algorithm contains separate AVX2, AVX-512, and scalar variants in a single file:
1. Split into distinct files: `<algo>_scalar.c`, `<algo>_avx2.c`, `<algo>_avx512.c`.
2. Keep the dispatch logic and public entry point in `<algo>_dispatch.c` and `<algo>.h`.
3. Use runtime CPUID detection to bind function pointers at initialization.
4. Allows CMake to attach target architecture flags (`-mavx2`, `-mavx512f`) per-file.

### Pattern C: Predicate & Inline Extraction
Extract deeply nested boolean conditions into `static inline bool is_<condition>(...)` functions
in internal headers. This flattens control flow and documents the domain intent.

---

## 3. Safe Refactoring Workflow

Follow this ordered protocol when executing a refactoring:

1. **Establish Baseline**:
   - Verify code compiles and passes tests: `make -C _build -j$(nproc)`
   - Run code size ratchet check: `./scripts/check_code_size.sh`
2. **Move or Split with `git mv`**:
   - Use `git mv` when renaming or relocating files to preserve commit history.
3. **Extract Functions Incrementally**:
   - Extract one helper at a time.
   - Maintain column-aligned parameters and Allman braces.
   - Ensure header hygiene: each new file includes only its required headers.
4. **Update Build Files**:
   - Add new source files to `SOURCEFILES` in `CMakeLists.txt`.
5. **Verify Correctness & Invariants**:
   - Recompile: `make -C _build -j$(nproc)`.
   - Run ratchet check: `./scripts/check_code_size.sh`.
   - Run test scripts: `bash tests/test_fps_components.sh`.
6. **Update Module Documentation**:
   - Update the submodule `README.md` source file table per `.agents/rules/readme-update.md`.
