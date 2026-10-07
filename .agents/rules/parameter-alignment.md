---
trigger: always_on
---

# Function Parameter Alignment

Multi-line function prototypes and definitions must use **column-aligned** parameter names for
readability. Pad each type (including any `*` prefix on the name) so that all parameter names
start at the same column.

## Rule

1. Each parameter on its own line, indented 4 spaces.
2. Pad the type with spaces so all parameter names begin at the same column within that prototype.
3. The alignment column is determined by the longest `type [*]` prefix among all parameters.
4. Use exactly **one space** between the last type token and the `*` or name when there is no need
   for padding (i.e., it is already the longest).
5. Qualifiers like `const` and `restrict` are part of the type for alignment purposes.

## Examples

**BAD** — unaligned:

```c
int atmturb_generate_screen(
    float *screen,
    long size,
    float outerscale,
    float innerscale)
```

**GOOD** — column-aligned:

```c
int atmturb_generate_screen(
    float       *screen,
    long         size,
    float        outerscale,
    float        innerscale)
```

**GOOD** — const and restrict:

```c
void atmturb_extrude_accumulate(
    const float *restrict master,
    long                  msize,
    double                x0,
    double                y0,
    long                  pup_size,
    float                 weight,
    float       *restrict out_pha)
```

## Scope

Applies to all `.c` and `.h` files in `src/`. Does not apply to single-parameter functions
(no alignment needed) or to macro arguments.
