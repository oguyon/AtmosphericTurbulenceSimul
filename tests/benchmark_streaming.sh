#!/usr/bin/env bash
#
# tests/benchmark_streaming.sh
# Automated benchmark testing runner for milkatmturb 2D streaming mode
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

BENCH_EXEC="${REPO_ROOT}/_build/atmturb-benchmark-stream"

if [[ ! -x "$BENCH_EXEC" ]]; then
    echo "Executable atmturb-benchmark-stream not found. Building..."
    make -C "${REPO_ROOT}/_build" atmturb-benchmark-stream -j$(nproc)
fi

echo "Running 2D wavefront streaming benchmark..."
"$BENCH_EXEC"
