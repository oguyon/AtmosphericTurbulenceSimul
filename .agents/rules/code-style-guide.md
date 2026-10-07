---
trigger: always_on
---

# Code Style Guide

- Use code blocks (`{ }`) to reduce the scope of variables as much as possible. Use this as an
  opportunity to add a comment immediately preceding the block explaining what the block is doing.
- Keep lines short, no more than 100 characters. This limit applies to C source and agent `.md`
  files.
- **Do not artificially split lines**: Avoid unnecessary vertical fragmentation. If an expression,
  variable assignment, or function call fits entirely within the 100-character limit, it must be
  written on a single line. Do not split short expressions across multiple lines, unless mandated
  by function prototype rules.
- Factorize code as much as possible. Extract shared logic into helper functions or macros instead
  of duplicating it.
- Loop indices should always be declared within the `for` statement (e.g., `for (int ii = 0; ...)`),
  unless they are meant to be used outside the loop. If used outside, add an explanatory comment.
- Function prototypes with arguments should be multi-line, with one line per argument.
- Document functions using Kernel-Doc style.
- Put a brief API-facing description in the `.h` file for declarations, and place the full
  Kernel-Doc
  comment above the function definition in the `.c` file.
- Follow the Linux kernel's C coding conventions for naming, scope, and control flow philosophy
  (early exit, `goto cleanup` pattern).
- **Brace style: Allman.** Place the opening brace on its own line for both function definitions
  and control-flow statements (`if`, `for`, `while`, `do`, `switch`).

  ```c
  /* Function definition */
  double compute_wavefront_phase(
      const double *restrict in,
      double       *restrict out,
      uint64_t               nelement)
  {
      /* Control flow — brace on next line */
      if (in == NULL)
      {
          return -1.0;
      }

      for (uint64_t ii = 0; ii < nelement; ii++)
      {
          out[ii] = in[ii] * 2.0;
      }

      return 0.0;
  }
  ```

- Enable and enforce compiler warnings (`-Wall`, `-Wextra`) during development to catch missing
  declarations early.
- Make sure every `.c` file strictly includes the exact headers it relies on, rather than
  implicitly relying on another header to include them.
- Conversely, do not include headers that the file does not use. Remove redundant `#include`s.
- Add a closing comment to any scope longer than ~10 lines:

  ```c
  #ifdef USE_IMAGESTREAMIO
  // ... many lines ...
  #endif // USE_IMAGESTREAMIO

  if (condition_met)
  {
      // ... many lines ...
  } // if (condition_met)
  ```

- Prefer early exit over deeply nested braces. The preferred control flow minimizes indentation:

  ```c
  /* GOOD — early exit */
  if (ptr == NULL)
  {
      return -1;
  }
  // main logic at low indentation ...

  /* BAD — unnecessary nesting */
  if (ptr != NULL)
  {
      // main logic indented ...
  }
  ```

- Use `restrict` on pointer parameters to array data in compute-heavy functions where pointers
  are guaranteed non-aliasing.
