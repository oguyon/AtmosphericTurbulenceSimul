---
name: pr-preparation
description: Standard check before preparing a Pull Request in milkatmturb.
---

# PR Preparation

Before opening a PR or merging code:

## Checklist
1. **Clean compilation:** Check that the code builds with zero warnings (`-Wall -Wextra`):
   ```bash
   make -C _build -j$(nproc)
   ```
2. **Code size ratchet:** Ensure `./scripts/check_code_size.sh` passes with zero violations.
3. **Format compliance:** Ensure code adheres to project C style (Allman braces, lines $\le 100$,
   column-aligned parameters).
4. **Tests passing:** Run test scripts:
   ```bash
   bash tests/test_fps_components.sh
   bash tests/test_wavefront_series.sh
   ```
5. **Clean Git tree:** Keep commit history structured and focus each commit on a single change.
6. **Agentic tool disclosure:** If agentic tools were used, disclose in the PR description with:
   - Notice: `Implemented by <model name>. Reviewed and signed off by O. Guyon.`
     (e.g., `Implemented by Gemini 3.8 flash. Reviewed and signed off by O. Guyon.`).
   - A short concise technical summary of what task the model performed.
7. **Clean SHM:** Ensure no lingering test streams in `/dev/shm/`.
