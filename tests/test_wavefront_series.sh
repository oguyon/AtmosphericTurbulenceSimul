#!/usr/bin/env bash
#
# tests/test_wavefront_series.sh
# End-to-end integration test for milkatmturb wavefront series simulation.
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Locate milk-fpsexec-atmturb-mkwfs
MKWFS_EXEC="${MILK_MKWFS_EXEC:-}"

if [[ -z "$MKWFS_EXEC" ]]; then
    CANDIDATES=(
        "$REPO_ROOT/../../_build/plugins/milkatmturb/milk-fpsexec-atmturb-mkwfs"
        "$REPO_ROOT/../_build/plugins/milkatmturb/milk-fpsexec-atmturb-mkwfs"
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
    echo "Please set MILK_MKWFS_EXEC or build milk with milkatmturb plugin enabled." >&2
    exit 1
fi

VALIDATE_EXEC="${MILK_VALIDATE_EXEC:-}"
if [[ -z "$VALIDATE_EXEC" ]]; then
    CANDIDATES=(
        "$REPO_ROOT/../../_build/plugins/milkatmturb/atmturb-validate-stats"
        "$REPO_ROOT/../_build/plugins/milkatmturb/atmturb-validate-stats"
        "/home/oguyon/src/milk-framework-dev/_build/plugins/milkatmturb/atmturb-validate-stats"
        "$(command -v atmturb-validate-stats 2>/dev/null || true)"
    )
    for c in "${CANDIDATES[@]}"; do
        if [[ -x "$c" ]]; then
            VALIDATE_EXEC="$c"
            break
        fi
    done
fi

if [[ -z "$VALIDATE_EXEC" || ! -x "$VALIDATE_EXEC" ]]; then
    echo "ERROR: atmturb-validate-stats executable not found." >&2
    exit 1
fi

echo "=== Running milkatmturb wavefront series test ==="
echo "Executable: $MKWFS_EXEC"

# Create a temporary scratch workspace
TEST_TMPDIR="$(mktemp -d /tmp/test_milkatmturb_XXXXXX)"
MILK_SHM_DIR="${MILK_SHM_DIR:-/milk/shm}"
FPS_NAME="twfs$$_mkwfs"
# unique FPS name avoids reusing a stale parameter set; remove it on exit
trap 'rm -rf "$TEST_TMPDIR"; rm -f "$MILK_SHM_DIR/${FPS_NAME}.fps.shm"' EXIT

cd "$TEST_TMPDIR"

# 1. Generate minimal 3-layer turbulence profile
cat > turbul.prof << 'EOF'
# altitude(m)   relativeCN2     speed(m/s)      direction(rad)
500.0           5.0             10.0            0.5
5000.0          3.0             20.0            1.2
12000.0         1.0             30.0            2.5
EOF

# 2. Generate minimal WFsim.conf
cat > WFsim.conf << 'EOF'
TURBULENCE_REF_WAVEL       0.500000
TURBULENCE_SEEING          0.60000
TURBULENCE_PROF_FILE       turbul.prof
ZENITH_ANGLE               0.0
SOURCE_XPOS                0.0
SOURCE_YPOS                0.0

WFOUTPUT                   1
WF_FILE_PREFIX             ./wf_
SHM_OUTPUT                 0

MAKE_SWAVEFRONT            1
SLAMBDA                    0.700
SWF_WRITE2DISK             1
SWF_FILE_PREFIX            ./swf_
SHM_SOUTPUT                0
SHM_SPREFIX                shmswf
SHM_SOUTPUTM               0

WFsize                     64
PUPIL_SCALE                0.180000

REALTIME                   0
REALTIMEFACTOR             1.0
WFTIME_STEP                0.01
TIME_SPAN                  0.05
NB_TSPAN                   1
SIMTDELAY                  0
WAITFORSEM                 0
WAITSEMIMNAME              

SKIP_EXISTING              0
WF_RAW_SIZE                64
MASTER_SIZE                256
WAVEFRONT_AMPLITUDE        0
FRESNEL_PROPAGATION        0
FRESNEL_PROPAGATION_BIN    1000.0
EOF

# 3. Execute wavefront simulation
echo "Starting simulation in $TEST_TMPDIR..."
"$MKWFS_EXEC" -n "$FPS_NAME" exec 0.7 0

# 4. Verify expected output files exist
EXPECTED_FILES=(
    "outarraypha.fits"
    "outarrayamp.fits"
    "outsarraypha.fits"
    "outsarrayamp.fits"
)

for f in "${EXPECTED_FILES[@]}"; do
    if [[ ! -f "$f" ]]; then
        echo "FAILED: Expected output file '$f' was not created." >&2
        exit 1
    fi
    size=$(stat -c%s "$f")
    if [[ "$size" -le 2880 ]]; then
        echo "FAILED: File '$f' size ($size bytes) is too small to contain data." >&2
        exit 1
    fi
    echo "  [OK] Found $f ($size bytes)"
done

# 5. Check data integrity and statistical sanity via atmturb-validate-stats
echo "Validating phase and amplitude cubes..."
for f in outarraypha.fits outsarraypha.fits; do
    "$VALIDATE_EXEC" finite "$f" --min-std 0.01 || {
        echo "FAILED: statistical sanity check on $f" >&2
        exit 1
    }
done

for f in outarrayamp.fits outsarrayamp.fits; do
    "$VALIDATE_EXEC" scint "$f" --tol-mean 1e-4 --tol-sigma 1e-4 || {
        echo "FAILED: amplitude check on $f" >&2
        exit 1
    }
done

echo "=== All tests passed successfully! ==="
