#!/usr/bin/env bash
#
# tests/test_fps_components.sh
# End-to-end test for milkatmturb standalone FPS executables:
#   - milk-fpsexec-atmturb-mkhvturb
#   - milk-fpsexec-atmturb-mkmastert
#   - milk-fpsexec-atmturb-mkvonkarman
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

find_executable() {
    local name="$1"
    local candidates=(
        "$REPO_ROOT/../../_build/plugins/milkatmturb/$name"
        "$REPO_ROOT/../_build/plugins/milkatmturb/$name"
        "/home/oguyon/src/milk-framework-dev/_build/plugins/milkatmturb/$name"
        "$(command -v "$name" 2>/dev/null || true)"
    )
    for c in "${candidates[@]}"; do
        if [[ -x "$c" ]]; then
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

WFPROP_EXEC=$(find_executable "milk-fpsexec-wfprop-fresnel") || {
    echo "ERROR: milk-fpsexec-wfprop-fresnel not found." >&2
    exit 1
}

echo "=== Testing milkatmturb FPS Components ==="
echo "mkhvturb:    $MKHVTURB_EXEC"
echo "mkmastert:   $MKMASTERT_EXEC"
echo "mkvonkarman: $MKVONKARMAN_EXEC"
echo "wfprop:      $WFPROP_EXEC"

TEST_TMPDIR="$(mktemp -d /tmp/test_fps_components_XXXXXX)"
trap 'rm -rf "$TEST_TMPDIR"' EXIT

cd "$TEST_TMPDIR"

# 1. Test Hufnagel-Valley profile generation
echo "Running atmturb-mkhvturb..."
"$MKHVTURB_EXEC" exec 21.0 0.15 4200.0 5 hv_test.prof
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
"$MKMASTERT_EXEC" exec 128 50.0 1.0 0 scr0 scr1
echo "  [OK] Master turbulence screens generated"

# 3. Test von Karman wind velocity series generation
echo "Running atmturb-mkvonkarman..."
"$MKVONKARMAN_EXEC" exec 512 0.1 20.0 50.0 vkwind
echo "  [OK] von Karman wind series synthesized"

# 4. Test Fresnel propagation component
echo "Running wfprop-fresnel..."
"$WFPROP_EXEC" exec wfin wfout 0.01 1000.0 0.5e-6
echo "  [OK] Fresnel diffractive propagation verified"

echo "=== All FPS component tests passed successfully! ==="
