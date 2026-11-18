#!/usr/bin/env bash
#
# tests/test_fps_components.sh
# End-to-end test for milkatmturb standalone FPS executables, milk-cli, and milk.bash:
#   - milk-fpsexec-atmturb-mkhvturb
#   - milk-fpsexec-atmturb-mkmastert
#   - milk-fpsexec-atmturb-mkvonkarman
#   - milk-fpsexec-atmturb-aoloop
#   - wfprop.fresnelprop & atmturb.* via milk-cli and milk.bash
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
MILK_INSTALL_DIR="/home/oguyon/src-frameworkdev/milk/_install/milk-1.03.00"

if [[ -d "$MILK_INSTALL_DIR" ]]; then
    export PATH="$MILK_INSTALL_DIR/bin:${PATH}"
    export LD_LIBRARY_PATH="$REPO_ROOT/_build:$MILK_INSTALL_DIR/lib:${LD_LIBRARY_PATH:-}"
else
    export LD_LIBRARY_PATH="$REPO_ROOT/_build:${LD_LIBRARY_PATH:-}"
fi

find_executable() {
    local name="$1"
    local candidates=(
        "$REPO_ROOT/_build/$name"
        "$MILK_INSTALL_DIR/bin/$name"
        "$REPO_ROOT/../WFpropagate/_build/$name"
        "$(command -v "$name" 2>/dev/null || true)"
    )
    for c in "${candidates[@]}"; do
        if [[ -n "$c" && -x "$c" ]]; then
            echo "$c"
            return 0
        fi
    done
    return 1
}

MKHVTURB_EXEC=$(find_executable "milk-fpsexec-atmturb-mkhvturb") || {
    echo "ERROR: milk-fpsexec-atmturb-mkhvturb not found." >&2
    exit 1
}

MKMASTERT_EXEC=$(find_executable "milk-fpsexec-atmturb-mkmastert") || {
    echo "ERROR: milk-fpsexec-atmturb-mkmastert not found." >&2
    exit 1
}

MKVONKARMAN_EXEC=$(find_executable "milk-fpsexec-atmturb-mkvonkarman") || {
    echo "ERROR: milk-fpsexec-atmturb-mkvonkarman not found." >&2
    exit 1
}

MKAOLOOP_EXEC=$(find_executable "milk-fpsexec-atmturb-aoloop") || {
    echo "ERROR: milk-fpsexec-atmturb-aoloop not found." >&2
    exit 1
}

MILK_CLI_EXEC=$(find_executable "milk-cli" || true)
WFPROP_EXEC=$(find_executable "milk-fpsexec-wfprop-fresnel" || true)

echo "=== Testing milkatmturb FPS Components ==="
echo "mkhvturb:    $MKHVTURB_EXEC"
echo "mkmastert:   $MKMASTERT_EXEC"
echo "mkvonkarman: $MKVONKARMAN_EXEC"
echo "aoloop:      $MKAOLOOP_EXEC"
echo "milk-cli:    ${MILK_CLI_EXEC:-none}"

TEST_TMPDIR="$(mktemp -d /tmp/test_fps_components_XXXXXX)"
if [[ -d /milk/shm-fdev ]]; then
    export MILK_SHM_DIR="/milk/shm-fdev"
else
    export MILK_SHM_DIR="${MILK_SHM_DIR:-/milk/shm}"
fi
FPS_PREFIX="tfpsc$$"

cleanup() {
    rm -rf "$TEST_TMPDIR"
    rm -rf "$MILK_SHM_DIR/${FPS_PREFIX}"_*.fps.shm "$MILK_SHM_DIR/${FPS_PREFIX}"_*.fps.d
    rm -rf "$MILK_SHM_DIR"/atmturb_*.fps.shm "$MILK_SHM_DIR"/atmturb_*.fps.d
    rm -rf "$MILK_SHM_DIR"/wfprop_*.fps.shm "$MILK_SHM_DIR"/wfprop_*.fps.d
}
trap cleanup EXIT

cd "$TEST_TMPDIR"

# 1. Test Hufnagel-Valley profile generation
echo "Running atmturb-mkhvturb..."
"$MKHVTURB_EXEC" -n "${FPS_PREFIX}_hv" exec 21.0 0.15 4200.0 5 hv_test.prof
if [[ ! -f "hv_test.prof" ]]; then
    echo "FAILED: hv_test.prof was not created." >&2
    exit 1
fi
line_count=$(grep -v '^#' hv_test.prof | grep -v '^[[:space:]]*$' | wc -l)
if [[ "$line_count" -ne 5 ]]; then
    echo "FAILED: Expected 5 layer lines in hv_test.prof, got $line_count" >&2
    exit 1
fi
echo "  [OK] Hufnagel-Valley profile generated (5 layers)"

# 2. Test master turbulence screens generation
echo "Running atmturb-mkmastert..."
"$MKMASTERT_EXEC" -n "${FPS_PREFIX}_mt" exec 128 50.0 1.0 0 scr0 scr1
echo "  [OK] Master turbulence screens generated"

# 3. Test von Karman wind velocity series generation
echo "Running atmturb-mkvonkarman..."
"$MKVONKARMAN_EXEC" -n "${FPS_PREFIX}_vk" exec 512 0.1 20.0 50.0 vkwind 1 vkwind.fits
if [[ ! -f "vkwind.fits" ]]; then
    echo "FAILED: vkwind.fits was not created." >&2
    exit 1
fi
echo "  [OK] von Karman wind series synthesized"

# 4. Test AO closed-loop simulator standalone executable
echo "Running atmturb-aoloop..."
"$MKAOLOOP_EXEC" -n "${FPS_PREFIX}_ao" exec none PSFflat psf_flat.fits 0 0.0 0.0 0 1.65 8.0
if [[ ! -f "psf_flat.fits" ]]; then
    echo "FAILED: psf_flat.fits was not created." >&2
    exit 1
fi
echo "  [OK] AO closed-loop standalone executable verified"

# 5. Test Fresnel propagation component (standalone or milk-cli catalog)
echo "Running wfprop-fresnel..."
if [[ -n "$WFPROP_EXEC" && -x "$WFPROP_EXEC" ]]; then
    "$WFPROP_EXEC" -n "${FPS_PREFIX}_fr" exec wfin wfout 0.01 1000.0 0.5e-6
elif [[ -n "$MILK_CLI_EXEC" && -x "$MILK_CLI_EXEC" ]]; then
    WFPROP_CMD="mem.mk2Dim wfre 64 64; mem.mk2Dim wfim 64 64; @calc 'wfre = 1.0';"
    WFPROP_CMD+=" mem.ri2c wfre wfim wfin; wfprop.fresnelprop wfin wfout 0.01 1000.0 0.5e-6"
    "$MILK_CLI_EXEC" -c "$WFPROP_CMD" >/dev/null
else
    echo "ERROR: neither milk-fpsexec-wfprop-fresnel nor milk-cli found." >&2
    exit 1
fi
echo "  [OK] Fresnel diffractive propagation verified"

# 6. Test milk.bash & milk-cli catalog commands + calculator outvars
if [[ -n "$MILK_CLI_EXEC" && -x "$MILK_CLI_EXEC" && -f "$MILK_INSTALL_DIR/bin/milk.bash" ]]; then
    echo "Running milk.bash & outvars integration check..."
    STREHL_VAL=$(bash -c "
        MILK_QUIET=1 source '$MILK_INSTALL_DIR/bin/milk.bash'
        atmturb.mkmastert 64 25.0 1.0 0 bscr0 bscr1 7 >/dev/null
        atmturb.measr0 bscr0 >/dev/null
        atmturb.aoloop none PSFbash psf_bash.fits 0 0.0 0.0 0 1.65 8.0 >/dev/null
        @calc -v sval 'strehl'
        @unload >/dev/null 2>&1 || true
        printf '%s\n' \"\$sval\"
    " | tail -n 1 | tr -d '[:space:]')
    python3 -c "
import sys
val = float('${STREHL_VAL}')
if abs(val - 1.0) >= 1e-3:
    print(f'FAILED: expected strehl ~ 1.0 from outvars, got {val}', file=sys.stderr)
    sys.exit(1)
"
    echo "  [OK] milk.bash catalog commands & outvars verified (strehl=$STREHL_VAL)"
fi

echo "=== All FPS component tests passed successfully! ==="
