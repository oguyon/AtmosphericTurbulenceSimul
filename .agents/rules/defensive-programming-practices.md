---
description: Defensive coding practices, buffer safety, pointer discipline, and resource limits.
---

# Defensive Programming Practices

`milkatmturb` structures its code defensively to guarantee robustness and stability in optical
simulations and real-time streaming.

## 1. Buffer and String Safety
- **Ban unbounded functions:** `strcpy()`, `sprintf()`, and `strcat()` are strictly forbidden.
- **Use bounded alternatives:** Always use `strncpy()`, `snprintf()`, and `strncat()`. Ensure the
  destination size is explicitly passed and that strings are null-terminated.
- Check format truncation warnings: ensure buffers (`char path[512]`) are large enough for paths
  and file prefixes.

## 2. Pointer Discipline
- **Initialization:** Always initialize pointers to `NULL` or valid memory immediately upon
  declaration.
- **Dereference safety:** Validate pointers (especially incoming arguments) against `NULL` before
  dereferencing.
- **Dangling pointers:** Immediately after calling `free()` on a pointer, set it to `NULL`.

## 3. Input Validation
- **Untrusted input:** Validate command-line arguments, FPS parameter inputs, and config file
  values.
- **Bounds check:** Validate configuration options (e.g. wavelength $\lambda > 0$, grid size $> 0$,
  screen dimensions $\ge$ pupil size, layer altitudes $\ge 0$) during initialization, not inside
  hot loops.

## 4. Integer Arithmetic and Bounds Checking
- **Safe arithmetic:** Prevent integer overflow. Be mindful of signed versus unsigned conversions.
- **Array access:** Validate array indices before accessing memory.
- **Hoist checks:** Hoist bounds checking and size validation outside of simulation loops.

## 5. State Initialization
- **Zero-initialization:** Prefer `calloc()` over `malloc()` for allocating structs to avoid
  uninitialized memory bugs. If `malloc()` is used, initialize all fields immediately.
- **Pre-allocation:** Allocate state structures and working image arrays during initialization.
  Never allocate inside per-frame wavefront generation or Fresnel loops.

## 6. Format String Safety
- **Safe formatting:** Never pass user input directly as the format string to `printf()` or
  `fprintf()`. Always use a literal format string: `fprintf(stderr, "%s", msg)`.

## 7. Safe Signal Handling
- Keep signal handlers minimal. Set a `volatile sig_atomic_t` flag (e.g., `stop_requested = 1`)
  that the simulation loop safely checks.
