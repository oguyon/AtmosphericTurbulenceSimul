---
description: Add a new function to an existing C file
---

# Add a Function to an Existing C File

Follow this workflow when adding a function to `milkatmturb`.

## Workflow
1. **Define prototype:** Add the prototype to the corresponding header file. Ensure parameters
   are column-aligned per `.agents/rules/parameter-alignment.md`.
2. **Implement logic:** Implement the function in the `.c` file using Allman brace style.
3. **Document:** Write Kernel-Doc style comments above the function definition.
4. **Header hygiene:** Ensure all headers required by the new logic are explicitly included in
   the `.c` file.
5. **Check code size:** Verify the function body remains under 60 lines (hard limit 150) and file
   under 600 lines (hard limit 1000). Run `./scripts/check_code_size.sh`.
6. **Compile and verify:** Compile the project to confirm there are no compiler warnings:
   ```bash
   make -C _build -j$(nproc)
   ```
