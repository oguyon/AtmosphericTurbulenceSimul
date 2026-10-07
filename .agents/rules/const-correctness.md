---
description: Use const on pointer parameters that are read-only, with milk-specific exceptions.
---

# Const-Correctness Policy

Use `const` to document intent and let the compiler catch accidental mutations. This rule
applies to **new code and code being actively refactored** — do not bulk-migrate existing functions.

## Required: String Parameters

All `char *` parameters that are not modified by the function MUST be `const char *`.

```c
/* WRONG — allows accidental mutation */
errno_t atmturb_load_profile(char *filename);

/* RIGHT */
errno_t atmturb_load_profile(const char *filename);
```

## Required: Input Pixel & Phase Data Pointers

Array pointers in compute functions that are read-only MUST be `const restrict`:

```c
/* WRONG — no const on read-only data */
void atmturb_extrude_screen(
    float       *restrict master,
    float       *restrict pupil,
    long                  npupil);

/* RIGHT */
void atmturb_extrude_screen(
    const float *restrict master,
    float       *restrict pupil,
    long                  npupil);
```

This applies to all typed array pointers (`float *`, `double *`, `complex float *`, etc.)
used for input data in compute-heavy functions.

## Encouraged: Public API Input Pointers

Input-only struct pointers in public API declarations (`.h` files) SHOULD be `const`:

```c
/* In .h — encouraged */
errno_t atmturb_print_profile_summary(
    const ATMTURB_PROFILE *prof);
```

Migrate opportunistically — when a function signature is already being changed for another reason,
add `const` to input pointers at the same time.

## Not Required: IMGID Parameters

`IMGID` is passed by value and contains a pointer to shared memory. Marking it `const`
only prevents reassigning the struct fields, not the underlying pixel data:

```c
/* Acceptable — const adds little value */
void func(IMGID img);
```

## Not Required: FPS Pointers

`FUNCTION_PARAMETER_STRUCT *fps` points to shared memory that other processes may modify
at any time. Using `const` is misleading because the data is inherently mutable.
It is acceptable to use `const` when the function genuinely does not write to FPS, but
it is not required:

```c
/* Acceptable — fps is SHM, always mutable */
errno_t fps_print_status(
    FUNCTION_PARAMETER_STRUCT *fps);
```

## Transition Strategy

- **New code:** follow this rule strictly.
- **Existing code:** add `const` when editing a function for other reasons.
- **Headers:** when updating a `.h` declaration, update the matching `.c` definition to match.
