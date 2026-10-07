---
description: Git worktree management for milkatmturb.
---

# Sync Worktree

To review PRs or work on isolated features without altering your active branch workspace:
1. Add a worktree:
   ```bash
   git worktree add ../atmturb-worktree main
   ```
2. Perform work and test builds in `../atmturb-worktree`.
3. Prune and remove when finished:
   ```bash
   git worktree remove ../atmturb-worktree
   ```
