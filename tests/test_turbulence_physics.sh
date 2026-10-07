#!/usr/bin/env bash
#
# tests/test_turbulence_physics.sh
# Quantitative physics regression tests for milkatmturb.
#
# Each scenario runs a standalone executable in a scratch directory and checks the output with
# tests/atmturb_validate_stats. Scenarios marked "xfail" document known physics defects
# that are scheduled for a later fix: they are reported but do not fail the suite (an
# unexpected pass is reported as XPASS so the marker can be removed).
#
# Usage: bash tests/test_turbulence_physics.sh [scenario_name ...]
#

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
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

VALIDATE_EXEC=$(find_executable atmturb-validate-stats) || {
    echo "ERROR: atmturb-validate-stats not found. Please build the plugin." >&2
    exit 1
}
VALIDATE=("$VALIDATE_EXEC")

MKWFS=$(find_executable milk-fpsexec-atmturb-mkwfs) || { echo "mkwfs not found" >&2; exit 1; }
MKVK=$(find_executable milk-fpsexec-atmturb-mkvonkarman) \
    || { echo "mkvonkarman not found" >&2; exit 1; }
MKHV=$(find_executable milk-fpsexec-atmturb-mkhvturb) || { echo "mkhvturb not found" >&2; exit 1; }
MKMT=$(find_executable milk-fpsexec-atmturb-mkmastert) \
    || { echo "mkmastert not found" >&2; exit 1; }

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
        [ZENITH_ANGLE]=0.0 [PARALLACTIC_ANGLE]=0.0 [SITE_ALT]=-1.0 [SEED]=1
        [SOURCE_XPOS]=0.0 [SOURCE_YPOS]=0.0
        [WFOUTPUT]=1 [WF_FILE_PREFIX]=./wf_ [SHM_OUTPUT]=0
        [MAKE_SWAVEFRONT]=1 [SLAMBDA]=1.65 [SWF_WRITE2DISK]=1 [SWF_FILE_PREFIX]=./swf_
        [SHM_SOUTPUT]=0 [SHM_SPREFIX]=shmswf [SHM_SOUTPUTM]=0
        [WFsize]=128 [PUPIL_SCALE]=0.02
        [REALTIME]=0 [REALTIMEFACTOR]=1.0 [WFTIME_STEP]=0.01 [TIME_SPAN]=0.2 [NB_TSPAN]=1
        [SIMTDELAY]=0 [WAITFORSEM]=0 [WAITSEMIMNAME]=none
        [SKIP_EXISTING]=0 [WF_RAW_SIZE]=128 [MASTER_SIZE]=1024
        [MASTER_OVERSAMPLE]=2 [INTERP]=1
        [WAVEFRONT_AMPLITUDE]=0 [FRESNEL_PROPAGATION]=0 [FRESNEL_PROPAGATION_BIN]=1000.0
        [LOWFREQ]=0 [ROLLING]=0
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

# T1: master turbulence screen structure function matches discrete PSD expectation
# and generation is thread-count invariant (1 vs 8 threads)
scenario_T1_master_screen_sf() {
    "$MKMT" -n "${FPS_PREFIX}_t1" exec 1024 100.0 0.0 0 scr0 scr1 42 scr0.fits scr1.fits \
        > mkmt.log 2>&1 || return 1
    "${VALIDATE[@]}" sf scr0.fits --reference discrete --pixscale 1.0 --r0 3.1788734 \
        --L0 100.0 --master-size 1024 --oversample 1 --lag-min 2 --lag-max 16 --tol 0.03 || return 1
    "${VALIDATE[@]}" sf scr1.fits --reference discrete --pixscale 1.0 --r0 3.1788734 \
        --L0 100.0 --master-size 1024 --oversample 1 --lag-min 2 --lag-max 16 --tol 0.03 || return 1

    OMP_NUM_THREADS=1 "$MKMT" -n "${FPS_PREFIX}_t1_th1" exec 512 100.0 0.0 0 s1_0 s1_1 12345 \
        s1_0.fits s1_1.fits > mkmt_th1.log 2>&1 || return 1
    OMP_NUM_THREADS=8 "$MKMT" -n "${FPS_PREFIX}_t1_th8" exec 512 100.0 0.0 0 s8_0 s8_1 12345 \
        s8_0.fits s8_1.fits > mkmt_th8.log 2>&1 || return 1
    "${VALIDATE[@]}" same s1_0.fits s8_0.fits --tol 1e-6 || return 1
    "${VALIDATE[@]}" same s1_1.fits s8_1.fits --tol 1e-6
}

# T2: single-layer phase screen has r0 = 0.98 lambda / seeing
scenario_T2_r0_single_layer() {
    one_layer_profile 4200 1.0 0.0 0.0 10000 0.0
    write_conf WFsim.conf
    run_mkwfs t2 1.65 0 || return 1
    "${VALIDATE[@]}" sf outarraypha.fits --reference discrete --seeing 0.6 --lam 0.5e-6 \
        --L0 10000 --master-size 1024 --pixscale 0.02 --lag-min 2 --lag-max 8 --tol 0.08 \
        --subtract-piston
}

# T2b: multi-layer profile yields aggregate r0 matching target seeing
scenario_T2b_r0_multi_layer() {
    write_conf WFsim.conf
    run_mkwfs t2b 1.65 0 || return 1
    "${VALIDATE[@]}" sf outarraypha.fits --reference discrete --seeing 0.6 \
        --lam 0.5e-6 --L0 50 --master-size 1024 --pixscale 0.02 --lag-min 2 \
        --lag-max 8 --tol 0.10 --subtract-piston
}

# T3: phase variance scales as 1/cos(z)
scenario_T3_airmass() {
    one_layer_profile 4200 1.0 0.0 0.0 10000 0.0
    mkdir -p z0 z60
    (cd z0 && cp ../turbul.prof . && write_conf WFsim.conf ZENITH_ANGLE=0.0 \
        && run_mkwfs t3a 1.65 0) || return 1
    (cd z60 && cp ../turbul.prof . && write_conf WFsim.conf ZENITH_ANGLE=1.0471976 \
        && run_mkwfs t3b 1.65 0) || return 1
    "${VALIDATE[@]}" ratio z60/outarraypha.fits z0/outarraypha.fits --expect 2.0 --variance \
        --tol 0.03 --subtract-piston
}

# T4: chromatic phase scaling matches wavelength ratio and air dispersion
scenario_T4_chromatic_ratio() {
    one_layer_profile 4200 1.0 0.0 0.0 10000 0.0
    write_conf WFsim.conf MAKE_SWAVEFRONT=1 SLAMBDA=1.65
    run_mkwfs t4 1.65 0 || return 1
    "${VALIDATE[@]}" ratio outsarraypha.fits outarraypha.fits --expect 0.2967 \
        --tol 0.01 --subtract-piston
}

# T4b: elevated layer exhibits chromatic shear from differential refraction
scenario_T4b_differential_refraction() {
    one_layer_profile 15000 1.0 0.0 0.0 10000 0.0
    write_conf WFsim.conf ZENITH_ANGLE=0.785398 PARALLACTIC_ANGLE=0.0 SITE_ALT=4200.0 \
        MAKE_SWAVEFRONT=1 SLAMBDA=1.65 TIME_SPAN=0.01
    run_mkwfs t4b 1.65 0 || return 1
    "${VALIDATE[@]}" shift outarraypha.fits outsarraypha.fits --dx 0.0 --dy -2.26 --tol 0.05
}

# T5: wind advection produces correct frame-to-frame displacement
scenario_T5_wind_advection_shift() {
    one_layer_profile 4200 1.0 10.0 0.0 10000 0.0
    write_conf WFsim.conf ZENITH_ANGLE=0.0 PARALLACTIC_ANGLE=0.0 WFTIME_STEP=0.01 TIME_SPAN=0.2
    run_mkwfs t5 1.65 0 || return 1
    "${VALIDATE[@]}" shift outarraypha.fits --dx -5.0 --dy 0.0 --step 1 --tol 0.05
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

# T10: sub-pixel advection breathing is suppressed under bicubic + oversampled extrusion
scenario_T10_breathing() {
    one_layer_profile 4200 1.0 10.0 0.0 10000 0.0
    write_conf WFsim.conf PUPIL_SCALE=0.04 WFTIME_STEP=0.001 TIME_SPAN=0.02 MASTER_SIZE=2048 \
        MASTER_OVERSAMPLE=2 INTERP=1
    run_mkwfs t10 1.65 0 || return 1
    "${VALIDATE[@]}" breathing outarraypha.fits --tol 0.20
}

# T10b: SIMD parity between scalar, AVX2, and AVX-512 backends
scenario_T10b_simd_parity() {
    "${VALIDATE[@]}" simd-parity --tol 1e-4
}

# T11: analytic subharmonic modes restore tip/tilt variance within 15% of Noll theory
scenario_T11_tilt_variance() {
    one_layer_profile 4200 1.0 10.0 0.0 0.0 0.0
    write_conf WFsim.conf PUPIL_SCALE=0.01 WFTIME_STEP=0.01 TIME_SPAN=2.0 MASTER_SIZE=512 \
        MASTER_OVERSAMPLE=1 INTERP=1 LOWFREQ=1 ROLLING=0
    run_mkwfs t11 1.65 0 || return 1
    "${VALIDATE[@]}" tilt outarraypha.fits --seeing 0.6 --lam 0.5e-6 --pixscale 0.01 --tol 0.15
}

# T12: rolling cross-faded screens decorrelate at multiples of wrap time (corr < 0.10)
scenario_T12_wrap_decorrelation() {
    one_layer_profile 4200 1.0 10.0 0.0 0.0 0.0
    write_conf WFsim.conf PUPIL_SCALE=0.01 WFTIME_STEP=0.01 TIME_SPAN=1.5 MASTER_SIZE=500 \
        MASTER_OVERSAMPLE=1 INTERP=1 LOWFREQ=0 ROLLING=1 SEED=42
    run_mkwfs t12 1.65 0 || return 1
    python3 - <<'EOF'
import sys
import numpy as np
from astropy.io import fits
cube = fits.getdata("outarraypha.fits")
for f in range(len(cube)):
    cube[f] -= cube[f].mean()
c50 = np.mean([np.corrcoef(cube[f].flat, cube[f+50].flat)[0, 1] for f in range(50)])
c100 = np.mean([np.corrcoef(cube[f].flat, cube[f+100].flat)[0, 1] for f in range(50)])
ok = abs(c50) < 0.10 and abs(c100) < 0.10
msg = f": wrap correlation c50 = {c50:+.4f}, c100 = {c100:+.4f}"
print(("PASS" if ok else "FAIL") + msg)
sys.exit(0 if ok else 1)
EOF
}

# T13: smooth cross-fading across epoch boundaries (correlation >= 0.998, variance flat < 5%)
scenario_T13_epoch_stability() {
    one_layer_profile 4200 1.0 0.0 0.0 0.0 0.0
    write_conf WFsim.conf PUPIL_SCALE=0.01 WFTIME_STEP=0.01 TIME_SPAN=1.0 MASTER_SIZE=500 \
        MASTER_OVERSAMPLE=1 INTERP=1 LOWFREQ=0 ROLLING=1 BOIL_TIME=0.2 SEED=1
    run_mkwfs t13 1.65 0 || return 1
    python3 - <<'EOF'
import sys
import numpy as np
from astropy.io import fits
cube = fits.getdata("outarraypha.fits")
c_b = np.corrcoef(cube[19].flat, cube[20].flat)[0, 1]
var19 = float(np.var(cube[19]))
var20 = float(np.var(cube[20]))
var_jump = abs(var20 - var19) / var19
ok = c_b >= 0.998 and var_jump < 0.05
msg = f": boundary corr = {c_b:.5f} (>= 0.998), var jump = {var_jump*100:.2f}% (< 5%)"
print(("PASS" if ok else "FAIL") + msg)
sys.exit(0 if ok else 1)
EOF
}

# T14: maximum 1-pixel spatial and temporal gradients are within 6 sigma of Gaussian
scenario_T14_seam_absence() {
    one_layer_profile 4200 1.0 10.0 0.0 0.0 0.0
    write_conf WFsim.conf PUPIL_SCALE=0.01 WFTIME_STEP=0.01 TIME_SPAN=1.5 MASTER_SIZE=500 \
        MASTER_OVERSAMPLE=1 INTERP=1 LOWFREQ=1 ROLLING=1 SEED=1
    run_mkwfs t14 1.65 0 || return 1
    python3 - <<'EOF'
import sys
import numpy as np
from astropy.io import fits
cube = fits.getdata("outarraypha.fits")
dx = np.diff(cube, axis=2)
dy = np.diff(cube, axis=1)
dt = np.diff(cube, axis=0)
max_dx = float(np.max(np.abs(dx)) / np.std(dx))
max_dy = float(np.max(np.abs(dy)) / np.std(dy))
max_dt = float(np.max(np.abs(dt)) / np.std(dt))
ok = max_dx <= 6.0 and max_dy <= 6.0 and max_dt <= 6.0
msg = f": max gradient: dx = {max_dx:.2f}s, dy = {max_dy:.2f}s, dt = {max_dt:.2f}s"
print(("PASS" if ok else "FAIL") + msg)
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
    "T1_master_screen_sf:pass"
    "T2_r0_single_layer:pass"
    "T2b_r0_multi_layer:pass"
    "T3_airmass:pass"
    "T4_chromatic_ratio:pass"
    "T4b_differential_refraction:pass"
    "T5_wind_advection_shift:pass"
    "T6_precision:pass"
    "T7_wind_components:pass"
    "T8a_hv_single_layer:pass"
    "T10_breathing:pass"
    "T10b_simd_parity:pass"
    "T11_tilt_variance:pass"
    "T12_wrap_decorrelation:pass"
    "T13_epoch_stability:pass"
    "T14_seam_absence:pass"
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
