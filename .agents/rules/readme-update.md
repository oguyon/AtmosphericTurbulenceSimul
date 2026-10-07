---
description: Update module README when source files are added or removed.
---

# README Update Standard

When a task adds, removes, or renames a source file (`.c` or `.h`) within any module directory:
- `src/AtmosphericTurbulence/`
- `src/AtmosphereModel/`
- `src/OpticsMaterials/`
- `src/WFpropagate/`

You MUST update that module's `README.md`.

## Trigger Conditions

1. A new `.c` or `.h` file is added to a module.
2. An existing source file is renamed or deleted.
3. A source file's purpose changes significantly (warrants a description update).

## Required Actions

1. Open the module's `README.md` (e.g., `src/AtmosphericTurbulence/README.md`).
2. Update the **source file table**:
   ```markdown
   | File | Description |
   |------|-------------|
   | foo.c | Brief description |
   ```
3. If the file is new, add a row with a concise description of what the file does.
4. If the file was removed, delete its row.
5. If the file was renamed, update the filename in the table.
6. Summarize documentation updates in the task summary.
