#!/usr/bin/env bash
#
# tests/test_wavefront_stream.sh
# End-to-end integration test for milkatmturb 2D frame-by-frame SHM streaming.
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Locate milk-fpsexec-atmturb-mkwfs
MKWFS_EXEC="${MILK_MKWFS_EXEC:-}"
if [[ -z "$MKWFS_EXEC" ]]; then
    CANDIDATES=(
        "$REPO_ROOT/_build/milk-fpsexec-atmturb-mkwfs"
        "/home/oguyon/src/milk-framework-dev/_build/plugins/milkatmturb/milk-fpsexec-atmturb-mkwfs"
        "$(command -v milk-fpsexec-atmturb-mkwfs 2>/dev/null || true)"
    )
    for c in "${CANDIDATES[@]}"; do
        if [[ -x "$c" ]]; then
            MKWFS_EXEC="$c"
            break
        fi
    done
fi

if [[ -z "$MKWFS_EXEC" || ! -x "$MKWFS_EXEC" ]]; then
    echo "ERROR: milk-fpsexec-atmturb-mkwfs executable not found." >&2
    exit 1
fi

export LD_LIBRARY_PATH="$REPO_ROOT/_build:/usr/local/milk/lib:${LD_LIBRARY_PATH:-}"

echo "=== Running milkatmturb 2D wavefront streaming test ==="
echo "Executable: $MKWFS_EXEC"

TEST_TMPDIR="$(mktemp -d /tmp/test_milkatmturb_stream_XXXXXX)"
MILK_SHM_DIR="${MILK_SHM_DIR:-/milk/shm}"
FPS_NAME="tstr$$_mkwfs"

cleanup() {
    tmux kill-session -t "$FPS_NAME" 2>/dev/null || true
    milk-fps-rm "$FPS_NAME" 2>/dev/null || true
    rm -rf "$TEST_TMPDIR"
    rm -f "$MILK_SHM_DIR/${FPS_NAME}.fps.shm"
}
trap cleanup EXIT

cd "$TEST_TMPDIR"

# -------------------------------------------------------------
# Test 1: Finite 2D Stream Mode (stream_mode = 2)
# -------------------------------------------------------------
echo "--- Test 1: Finite 2D streaming (stream_mode = 2) ---"
"$MKWFS_EXEC" -n "$FPS_NAME" fpsinit
milk-fps-set "$FPS_NAME.wfsize" 128
milk-fps-set "$FPS_NAME.time_step" 0.01
milk-fps-set "$FPS_NAME.time_span" 0.05
milk-fps-set "$FPS_NAME.stream_mode" 2
milk-fps-set "$FPS_NAME.amplitude" 1
milk-fps-set "$FPS_NAME.save_fits" 0

"$MKWFS_EXEC" -n "$FPS_NAME" runstart

# Verify that outarraypha is 2D 128x128
PHAINF="$(milk-stream-info outarraypha)"
echo "$PHAINF" | grep -q "2D  128 x 128" || {
    echo "FAILED: outarraypha is not 2D 128x128!" >&2
    echo "$PHAINF"
    exit 1
}

# Verify that outarrayamp is 2D 128x128
AMPINF="$(milk-stream-info outarrayamp)"
echo "$AMPINF" | grep -q "2D  128 x 128" || {
    echo "FAILED: outarrayamp is not 2D 128x128!" >&2
    echo "$AMPINF"
    exit 1
}

CNT0_PHA=$(echo "$PHAINF" | awk '/cnt0/ {print $3}' | tr -dc '0-9')
if [[ "$CNT0_PHA" -lt 5 ]]; then
    echo "FAILED: outarraypha cnt0 ($CNT0_PHA) < expected 5 frames" >&2
    exit 1
fi
echo "  [OK] Finite 2D stream verified: 128x128, cnt0=$CNT0_PHA"

# -------------------------------------------------------------
# Test 2: Continuous 2D Stream Mode with runstop (stream_mode = 1)
# -------------------------------------------------------------
echo "--- Test 2: Continuous 2D streaming with runstop (stream_mode = 1) ---"
milk-fps-set "$FPS_NAME.stream_mode" 1
milk-fps-set "$FPS_NAME.time_step" 0.005

# Start inside tmux session
"$MKWFS_EXEC" -n "$FPS_NAME" -tmux runstart
sleep 1.5

# Sample counter
CNT_A=$(milk-stream-info outarraypha | awk '/cnt0/ {print $3}' | tr -dc '0-9')
sleep 0.8
CNT_B=$(milk-stream-info outarraypha | awk '/cnt0/ {print $3}' | tr -dc '0-9')

if [[ "$CNT_B" -le "$CNT_A" ]]; then
    echo "FAILED: outarraypha cnt0 did not advance ($CNT_A -> $CNT_B)" >&2
    exit 1
fi
echo "  [OK] Continuous streaming active: cnt0 advanced from $CNT_A to $CNT_B"

# Graceful stop via milk-fps-runstop
milk-fps-runstop "$FPS_NAME"
sleep 1.0

CNT_STOP=$(milk-stream-info outarraypha | awk '/cnt0/ {print $3}' | tr -dc '0-9')
sleep 0.5
CNT_AFTER=$(milk-stream-info outarraypha | awk '/cnt0/ {print $3}' | tr -dc '0-9')

if [[ "$CNT_STOP" -ne "$CNT_AFTER" ]]; then
    echo "FAILED: stream did not halt after runstop ($CNT_STOP -> $CNT_AFTER)" >&2
    exit 1
fi
echo "  [OK] Stream cleanly halted after runstop at frame $CNT_STOP"

echo "=== All 2D streaming tests passed successfully! ==="
