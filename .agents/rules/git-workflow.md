---
description: Git branching, committing, and worktree workflow.
---

# Git Workflow

The `AtmosphericTurbulenceSimul` repository uses a standard git workflow centered on the `main`
  branch.

## 1. Branching Strategy
- All development should occur on feature branches branched from `main`.
- Merge feature branches into `main` using Pull Requests after compiling, passing the code size
  ratchet, and verifying test scripts.

## 2. Commit Guidelines
- Write clear, descriptive commit messages using conventional prefixes:
  - `feat(<component>): ...`
  - `fix(<component>): ...`
  - `refactor(<component>): ...`
  - `perf(<component>): ...`
  - `docs(<component>): ...`
  - `test(<component>): ...`
- Example: `fix(atmturb): replace fragile CLI image ops in master screens with direct FFTW3`

## 3. Worktree Usage
For maintaining multiple branches or reviewing PRs, use Git worktrees rather than cloning
  multiple times:
```bash
git worktree add ../atmturb-review feature-branch
```

## 4. Pull Requests & Agentic Tool Disclosure
When opening a Pull Request:
- Ensure all code compiles cleanly without warnings, passes `./scripts/check_code_size.sh`, and
  passes test scripts.
- Disclose the use of agentic tools in the PR description with:
  - An attribution and sign-off message matching the format:
    `Implemented by <model name>. Reviewed and signed off by O. Guyon.`
    (e.g., `Implemented by Gemini 3.8 flash. Reviewed and signed off by O. Guyon.`).
  - A short concise technical summary of what task the model performed.
