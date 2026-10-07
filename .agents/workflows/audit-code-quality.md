---
description: Code quality auditing of C sources.
---

# Audit Code Quality

Audit C source files and headers in `milkatmturb` for:
- [ ] Correct brace style (Allman: opening brace `{` on its own line).
- [ ] No lines exceeding 100 characters in `.c`, `.h`, `.md`, or scripts.
- [ ] Zero implicit header dependencies (include what is used).
- [ ] Column-aligned multi-line function prototypes.
- [ ] Code size compliance: verify `./scripts/check_code_size.sh` passes with 0 violations.
- [ ] Type consistency: loop index types match bound types, float functions for float arrays.
- [ ] Proper error handling, return value checks, and `goto cleanup;` resource release.
- [ ] Module `README.md` updated if source files were added, renamed, or deleted.
