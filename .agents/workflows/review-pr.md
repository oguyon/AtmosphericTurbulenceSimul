---
description: Review a pull request for coding standards compliance
---

# Review a Pull Request

Use this workflow to systematically review a PR for compliance with `milkatmturb` conventions.

## 1. Read the PR

Get the PR diff and file list:
```bash
git diff main..HEAD --stat
git diff main..HEAD
```

## 2. Coding Style Check

For each modified `.c` and `.h` file, verify:
- [ ] Lines $\le 100$ characters
- [ ] Allman brace style
- [ ] Kernel-Doc comments above functions
- [ ] Explicit `#include` for every header used (header hygiene)
- [ ] Code blocks `{ }` minimize variable scope
- [ ] Multi-line function prototypes with column-aligned parameters

## 3. Code Size & Architecture Check

- [ ] `./scripts/check_code_size.sh` passes with 0 violations
- [ ] Files $\le 600$ lines, function bodies $\le 60$ lines
- [ ] No new cyclic dependencies between submodules
- [ ] Submodule `README.md` updated if source files were added/renamed/removed

## 4. Performance Check

For simulation and compute changes:
- [ ] `restrict` on pixel pointer parameters
- [ ] Float math uses `f` variants (`sqrtf`, `sinf`, `cosf`)
- [ ] No `malloc`/`calloc`/`free` in per-frame simulation loops
- [ ] No per-frame FFTW plan creation/destruction
- [ ] OpenMP loop indices properly scoped

## 5. Build and Test

```bash
make -C _build -j$(nproc)
./scripts/check_code_size.sh
bash tests/test_fps_components.sh
bash tests/test_wavefront_series.sh
```

## 6. Agentic Tool Disclosure Check

If agentic tools were used, verify the PR description contains:
- `Implemented by <model name>. Reviewed and signed off by O. Guyon.`
- Concise technical summary of the model's work.
