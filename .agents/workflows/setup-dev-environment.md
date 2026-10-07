---
description: Setup development dependencies for milkatmturb.
---

# Setup Dev Environment

Install mandatory and optional development dependencies on Debian/Ubuntu:

```bash
# Core build tools
sudo apt update
sudo apt install build-essential cmake pkg-config

# Mathematical and simulation libraries
sudo apt install libfftw3-dev libgsl-dev libomp-dev

# Astronomy and image I/O
sudo apt install libcfitsio-dev
```

If building as a plugin inside `milk`, ensure `milk` is compiled and installed, or set
`MILK_ROOT` / `MILK_SOURCE_DIR` to the milk source path.
