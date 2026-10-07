---
description: Pull request standards, verification requirements, and agentic tool disclosure.
---

# Pull Request Standards

When preparing or submitting Pull Requests (PRs):

## 1. Agentic Tool Disclosure
Any PR implementing changes with the assistance of agentic tools must explicitly disclose this
in the PR description:
- **Attribution and Sign-off**: Include a statement matching the format:
  `Implemented by <model name>. Reviewed and signed off by O. Guyon.`
  Example:
  `Implemented by Gemini 3.8 flash. Reviewed and signed off by O. Guyon.`
- **Prompt Summary**: Include a short concise technical summary of what task the model performed
  (avoiding arbitrary option/phase numbers that mean nothing to the reader).

## 2. Verification Checklist
- **Compilation**: Code builds cleanly with zero warnings (`-Wall -Wextra`).
- **Code Size Ratchet**: `./scripts/check_code_size.sh` passes with zero violations.
- **Style**: Adheres to project C style (Allman braces, <= 100 chars, parameter alignment).
- **Testing**: Automated test scripts (`test_fps_components.sh`, `test_wavefront_series.sh`) pass.
- **Git History**: Clean, focused commit messages following conventional commit prefixes.
- **No Leaks**: No dangling shared memory files in `/dev/shm/` or leaked pointers/descriptors.
