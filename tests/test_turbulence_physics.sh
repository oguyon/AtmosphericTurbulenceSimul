#!/usr/bin/env bash
#
# tests/test_turbulence_physics.sh
# Quantitative physics regression tests for milkatmturb.
#
# Each scenario runs a standalone executable in a scratch directory and checks the output with
# tests/validate_turbulence_stats.py. Scenarios marked "xfail" document known physics defects
# that are scheduled for a later fix: they are reported but do not fail the suite (an
# unexpected pass is reported as XPASS so the marker can be removed).
#
# Usage: bash tests/test_turbulence_physics.sh [scenario_name ...]
#

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
VALIDATE=(python3 "$SCRIPT_DIR/validate_turbulence_stats.py")
MILK_SHM_DIR="${MILK_SHM_DIR:-/milk/shm}"
FPS_PREFIX="tphys$$"

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

MKWFS=$(find_executable milk-fpsexec-atmturb-mkwfs) || { echo "mkwfs not found" >&2; exit 1; }
MKVK=$(find_executable milk-fpsexec-atmturb-mkvonkarman) \
    || { echo "mkvonkarman not found" >&2; exit 1; }
MKHV=$(find_executable milk-fpsexec-atmturb-mkhvturb) || { echo "mkhvturb not found" >&2; exit 1; }

"${VALIDATE[@]}" --help >/dev/null 2>&1
if [[ $? -eq 77 ]]; then
    echo "SKIP: numpy/astropy/scipy not available"
    exit 0
fi

TEST_TMPDIR="$(mktemp -d /tmp/test_turb_physics_XXXXXX)"
cleanup() {
    rm -rf "$TEST_TMPDIR"
    rm -f "$MILK_SHM_DIR/${FPS_PREFIX}"_*.fps.shm
}
trap cleanup EXIT

# ----------------------------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------------------------

# write_conf <file> [KEY=VALUE ...] : write WFsim.conf with defaults overridden by arguments
write_conf() {
    local file="$1"
    shift
    declare -A kv=(
        [TURBULENCE_REF_WAVEL]=0.5 [TURBULENCE_SEEING]=0.6 [TURBULENCE_PROF_FILE]=turbul.prof
        [ZENITH_ANGLE]=0.0 [SOURCE_XPOS]=0.0 [SOURCE_YPOS]=0.0
        [WFOUTPUT]=1 [WF_FILE_PREFIX]=./wf_ [SHM_OUTPUT]=0
        [MAKE_SWAVEFRONT]=1 [SLAMBDA]=1.65 [SWF_WRITE2DISK]=1 [SWF_FILE_PREFIX]=./swf_
        [SHM_SOUTPUT]=0 [SHM_SPREFIX]=shmswf [SHM_SOUTPUTM]=0
        [WFsize]=128 [PUPIL_SCALE]=0.02
        [REALTIME]=0 [REALTIMEFACTOR]=1.0 [WFTIME_STEP]=0.01 [TIME_SPAN]=0.2 [NB_TSPAN]=1
        [SIMTDELAY]=0 [WAITFORSEM]=0 [WAITSEMIMNAME]=none
        [SKIP_EXISTING]=0 [WF_RAW_SIZE]=128 [MASTER_SIZE]=1024
        [WAVEFRONT_AMPLITUDE]=0 [FRESNEL_PROPAGATION]=0 [FRESNEL_PROPAGATION_BIN]=1000.0
    )
    local arg
    for arg in "$@"; do
        kv["${arg%%=*}"]="${arg#*=}"
    done
    : > "$file"
    local key
    for key in "${!kv[@]}"; do
        printf '%-26s %s\n' "$key" "${kv[$key]}" >> "$file"
    done
}

# run_mkwfs <fpsname> <slambda_um> <precision> : run wavefront series in the current directory
run_mkwfs() {
    "$MKWFS" -n "${FPS_PREFIX}_$1" exec "$2" "$3" > mkwfs.log 2>&1
}

# one_layer_profile <alt> <cn2> <speed> <dir> <L0> <l0>
one_layer_profile() {
    printf '# alt cn2 speed dir L0 l0\n%s %s %s %s %s %s\n' "$@" > turbul.prof
}

# ----------------------------------------------------------------------------------------------
# Scenarios (each runs in its own directory; return 0 = pass)
# ----------------------------------------------------------------------------------------------

# T2: single-layer phase screen has r0 = 0.98 lambda / seeing (scheduled fix: Phase 3)
scenario_T2_r0_single_layer() {
    one_layer_profile 4200 1.0 10.0 0.3 10000 0.0
    write_conf WFsim.conf
    run_mkwfs t2 1.65 0 || return 1
    "${VALIDATE[@]}" sf outarraypha.fits --seeing 0.6 --lam 0.5e-6 --L0 10000 \
        --pixscale 0.02 --lag-min 2 --lag-max 8 --tol 0.05 --subtract-piston
}

# T3: phase variance scales as 1/cos(z) (scheduled fix: Phase 3)
scenario_T3_airmass() {
    one_layer_profile 4200 1.0 10.0 0.3 10000 0.0
    mkdir -p z0 z60
    (cd z0 && cp ../turbul.prof . && write_conf WFsim.conf ZENITH_ANGLE=0.0 \
        && run_mkwfs t3a 1.65 0) || return 1
    (cd z60 && cp ../turbul.prof . && write_conf WFsim.conf ZENITH_ANGLE=1.0471976 \
        && run_mkwfs t3b 1.65 0) || return 1
    "${VALIDATE[@]}" ratio z60/outarraypha.fits z0/outarraypha.fits --expect 2.0 --variance \
        --tol 0.03 --subtract-piston
}

# T6: double-precision screen synthesis gives the same wavefronts as single precision
scenario_T6_precision() {
    one_layer_profile 4200 1.0 10.0 0.3 10000 0.0
    mkdir -p p0 p1
    (cd p0 && cp ../turbul.prof . && write_conf WFsim.conf && run_mkwfs t6a 1.65 0) || return 1
    (cd p1 && cp ../turbul.prof . && write_conf WFsim.conf && run_mkwfs t6b 1.65 1) || return 1
    "${VALIDATE[@]}" finite p1/outarraypha.fits --min-std 1e-6 || return 1
    "${VALIDATE[@]}" same p0/outarraypha.fits p1/outarraypha.fits --tol 1e-3
}

# T7: von Karman wind components are mutually uncorrelated
# (series spans ~5e4 outer scales so the sample correlation noise is ~0.01)
scenario_T7_wind_components() {
    "$MKVK" -n "${FPS_PREFIX}_t7" exec 262144 1.0 2.0 5.0 vkw 7 vkw.fits > mkvk.log 2>&1 \
        || return 1
    "${VALIDATE[@]}" corr vkw.fits --plane-a 0 --plane-b 1 --max-corr 0.05 || return 1
    "${VALIDATE[@]}" corr vkw.fits --plane-a 1 --plane-b 2 --max-corr 0.05 || return 1
    "${VALIDATE[@]}" corr vkw.fits --plane-a 0 --plane-b 2 --max-corr 0.05
}

# T8a: single-layer Hufnagel-Valley profile is well formed (no NaN, Cn2 fraction ~ 1)
scenario_T8a_hv_single_layer() {
    "$MKHV" -n "${FPS_PREFIX}_t8a" exec 21.0 0.15 4200.0 1 hv1.prof > mkhv.log 2>&1 \
        || return 1
    python3 - <<'EOF'
import math, sys
rows = [l.split() for l in open("hv1.prof") if l.strip() and not l.startswith("#")]
ok = len(rows) == 1 and all(math.isfinite(float(v)) for v in rows[0])
ok = ok and abs(float(rows[0][1]) - 1.0) < 0.01 and 4200.0 < float(rows[0][0]) < 30000.0
print(("PASS" if ok else "FAIL") + f": single-layer HV profile row = {rows}")
sys.exit(0 if ok else 1)
EOF
}

# T15: invalid geometry (PUPIL_SCALE = 0) is rejected without producing output
scenario_T15_invalid_pupil_scale() {
    one_layer_profile 4200 1.0 10.0 0.3 10000 0.0
    write_conf WFsim.conf PUPIL_SCALE=0.0
    run_mkwfs t15 1.65 0
    if [[ -f outarraypha.fits ]]; then
        echo "FAIL: output produced with PUPIL_SCALE = 0"
        return 1
    fi
    grep -q "PUPIL_SCALE must be > 0" mkwfs.log && echo "PASS: PUPIL_SCALE = 0 rejected"
}

# ----------------------------------------------------------------------------------------------
# Runner
# ----------------------------------------------------------------------------------------------

# name:expectation  (expectation = pass | xfail)
SCENARIOS=(
    "T2_r0_single_layer:xfail"
    "T3_airmass:xfail"
    "T6_precision:pass"
    "T7_wind_components:pass"
    "T8a_hv_single_layer:pass"
    "T15_invalid_pupil_scale:pass"
)

npass=0; nfail=0; nxfail=0; nxpass=0
for entry in "${SCENARIOS[@]}"; do
    name="${entry%%:*}"
    expect="${entry##*:}"
    if [[ $# -gt 0 && ! " $* " =~ " $name " ]]; then
        continue
    fi
    dir="$TEST_TMPDIR/$name"
    mkdir -p "$dir"
    out=$(cd "$dir" && "scenario_$name" 2>&1)
    status=$?
    if [[ $status -eq 0 && $expect == pass ]]; then
        result="PASS"; npass=$((npass + 1))
    elif [[ $status -ne 0 && $expect == xfail ]]; then
        result="XFAIL"; nxfail=$((nxfail + 1))
    elif [[ $status -eq 0 && $expect == xfail ]]; then
        result="XPASS"; nxpass=$((nxpass + 1))
    else
        result="FAIL"; nfail=$((nfail + 1))
    fi
    printf '[%-5s] %s\n' "$result" "$name"
    echo "$out" | grep -E "^(PASS|FAIL|SKIP)" | sed 's/^/          /'
    if [[ $result == FAIL ]]; then
        for log in "$dir"/*.log "$dir"/*/*.log; do
            [[ -f "$log" ]] && { echo "  --- $log (tail)"; tail -n 15 "$log"; }
        done
    fi
done

echo "=== physics tests: $npass passed, $nfail failed, $nxfail xfail, $nxpass xpass ==="
[[ $nfail -eq 0 ]]
