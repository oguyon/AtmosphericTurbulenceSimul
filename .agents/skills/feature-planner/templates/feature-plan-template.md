# Feature Plan: [Title]

## Overview

**Goal**: [One-sentence description of the feature or refactoring]
**Scope**: [Submodule / Plugin / Standalone Executable] · [New Feature / Refactor / Optimization]
**Template**: [FPS V2 standalone / compute unit / module function / CLI command]

## Architecture Impact

### Affected Submodules

- **Affected submodules**: [`AtmosphericTurbulence` / `AtmosphereModel` / `OpticsMaterials` /
  `WFpropagate`]
- **New dependencies**: [list or "none"]

### Shared Memory & Image Streams

| Object | Type | Name | Details |
| ------ | ---- | ---- | ------- |
|        | Stream / FPS / Procinfo | | dtype, dims, purpose |

### CLI & FPS Surface

- **New commands**: [list or "none"]
- **New standalone executables**: [list or "none"]
- **Modified CLI commands**: [list or "none"]

## File Changes

| Action | File | Purpose |
| ------ | ---- | ------- |
| NEW    | `src/...` | description |
| MODIFY | `src/...` | what changes |
| MODIFY | `CMakeLists.txt` | build target updates |

## Documentation Updates

- [ ] Module README (`src/<Module>/README.md`)
- [ ] Top-level README (`README.md`)
- [ ] Code size ratchet verification (`./scripts/check_code_size.sh`)

## Implementation Phases

### Phase 1: [Foundation]
**Changes**: [list]
**Verify**: [compile-test, unit checks]

### Phase 2: [Core Logic]
**Changes**: [list]
**Verify**: [compile-test, test scripts]

### Phase 3: [Integration & Tests]
**Changes**: [CLI wiring, docs, tests]
**Verify**: [end-to-end tests, code size ratchet]

## Risks & Open Questions

1. [Any trade-offs or decisions requiring input]
