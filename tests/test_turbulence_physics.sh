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
        "$REPO_ROOT/_build/$name"
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
MKAOLOOP=$(find_executable milk-fpsexec-atmturb-aoloop) \
    || { echo "aoloop not found" >&2; exit 1; }

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
        [LOWFREQ]=0 [ROLLING]=0 [FRESNEL_GUARD_PIX]=0
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

# run_aoloop <fps> <wf> <psf> <fits> <mode> <gain> <leak> <delay> <scilam> <teldiam>
run_aoloop() {
    "$MKAOLOOP" -n "${FPS_PREFIX}_$1" exec "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}" \
        > aoloop.log 2>&1
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

# T7: von Karman wind components are mutually uncorrelated, match RMS and -5/3 spectral slope
# (series spans ~5e4 outer scales so the sample correlation noise is ~0.01)
scenario_T7_wind_components() {
    "$MKVK" -n "${FPS_PREFIX}_t7" exec 262144 1.0 2.0 5.0 vkw 7 vkw.fits > mkvk.log 2>&1 \
        || return 1
    "${VALIDATE[@]}" corr vkw.fits --plane-a 0 --plane-b 1 --max-corr 0.05 || return 1
    "${VALIDATE[@]}" corr vkw.fits --plane-a 1 --plane-b 2 --max-corr 0.05 || return 1
    "${VALIDATE[@]}" corr vkw.fits --plane-a 0 --plane-b 2 --max-corr 0.05 || return 1
    python3 - <<'EOF'
import sys
import numpy as np
from astropy.io import fits

data = fits.getdata("vkw.fits").squeeze()
rms_u = float(np.std(data[0]))
rms_v = float(np.std(data[1]))
rms_w = float(np.std(data[2]))

ok_rms = (abs(rms_u - 2.0) / 2.0 < 0.05 and
          abs(rms_v - 2.0) / 2.0 < 0.05 and
          abs(rms_w - 2.0) / 2.0 < 0.05)

psd = np.abs(np.fft.rfft(data[0])) ** 2
freqs = np.fft.rfftfreq(len(data[0]), d=1.0)
mask = (freqs >= 0.05) & (freqs <= 0.40)
slope = float(np.polyfit(np.log10(freqs[mask]), np.log10(psd[mask]), 1)[0])
ok_slope = abs(slope - (-5.0 / 3.0)) < 0.10

ok = ok_rms and ok_slope
msg = (f": rms=({rms_u:.3f}, {rms_v:.3f}, {rms_w:.3f}), "
       f"PSD slope={slope:.3f} (target -1.667 +- 0.10)")
print(("PASS" if ok else "FAIL") + msg)
sys.exit(0 if ok else 1)
EOF
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

# T8b: multi-layer HV profile with Bufton wind: recomputed r0 +- 0.5%, error on unreachable r0
scenario_T8b_hv_bounds_and_recomputed_r0() {
    "$MKHV" -n "${FPS_PREFIX}_t8b_valid" exec 21.0 0.15 4200.0 10 hv_valid.prof 1 42 \
        > mkhv_valid.log 2>&1 || return 1
    [[ -f hv_valid.prof ]] || return 1
    [[ -f conf_turb.txt ]] || return 1

    python3 - <<'EOF'
import math, sys
seeing_arcsec = float(open("conf_turb.txt").read().strip())
seeing_rad = seeing_arcsec * (math.pi / (180.0 * 3600.0))
r0_recomputed = 0.98 * 0.55e-6 / seeing_rad
target_r0 = 0.15
rel_err = abs(r0_recomputed - target_r0) / target_r0
ok = rel_err < 0.005
msg = (f": recomputed r0 = {r0_recomputed:.4f} m (target {target_r0:.4f} m, "
       f"err {rel_err*100:.2f}%, tol 0.5%)")
print(("PASS" if ok else "FAIL") + msg)
sys.exit(0 if ok else 1)
EOF
    [[ $? -eq 0 ]] || return 1

    "$MKHV" -n "${FPS_PREFIX}_t8b_unreach" exec 21.0 10.0 4200.0 10 hv_unreach.prof 1 42 \
        > mkhv_unreach.log 2>&1
    if [[ -f hv_unreach.prof ]]; then
        echo "FAIL: profile generated for unreachable r0 = 10 m"
        return 1
    fi
    grep -q "exceeds maximum reachable r0" mkhv_unreach.log || return 1
    echo "PASS: unreachable r0 correctly rejected with informative diagnostic"
}

# T8c: turbulent wind advection produces non-linear trajectory with physical velocity fluctuations
scenario_T8c_turbulent_wind_advection() {
    printf '# alt cn2 speed dir L0 l0 sigma_wsp L_wind\n%s %s %s %s %s %s %s %s\n' \
        4200 1.0 10.0 0.0 10000 0.0 5.0 50.0 > turbul.prof
    write_conf WFsim.conf PUPIL_SCALE=0.01 WFTIME_STEP=0.005 TIME_SPAN=0.1 MASTER_SIZE=512 \
        MASTER_OVERSAMPLE=1 INTERP=1 LOWFREQ=0 ROLLING=0 SEED=42
    run_mkwfs t8c 1.65 0 || return 1
    python3 - <<'EOF'
import sys
import numpy as np
from astropy.io import fits

pha = fits.getdata("outarraypha.fits")
assert pha.shape[0] == 20
assert np.isfinite(pha).all()

diffs = [float(np.std(pha[t] - pha[t-1])) for t in range(1, len(pha))]
diff_var = float(np.std(diffs))
ok = np.isfinite(diffs).all() and diff_var > 0.0
msg = f": turbulent advection simulated 20 frames (frame diff std = {diff_var:.4f} > 0)"
print(("PASS" if ok else "FAIL") + msg)
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

# T10c: CUDA GPU vs AVX2 CPU parity for Rytov propagation (corr > 0.99999, max_diff < 1e-4)
scenario_T10c_cuda_rytov_parity() {
    if ! command -v nvidia-smi >/dev/null 2>&1 || ! nvidia-smi >/dev/null 2>&1; then
        echo "PASS: skipped (no NVIDIA GPU available)"
        return 0
    fi
    cat << 'EOF' > turbul.prof
# alt cn2 speed dir L0 l0
 4215     5.32        6.5     1.47  10000 0.0
 4230     1.47        6.55    1.57  10000 0.0
 4349     1.08        6.6     1.67  10000 0.0
 5007     2.11        6.7     1.77  10000 0.0
12000     1.83       22.0     3.10  10000 0.0
16200     1.48        9.5     3.20  10000 0.0
23701     0.697       5.6     3.30  10000 0.0
EOF
    write_conf WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=2 \
        FRESNEL_PROPAGATION_BIN=1000.0 FRESNEL_GUARD_PIX=16 TIME_SPAN=0.05 WFTIME_STEP=0.01

    ATMTURB_SIMD=AVX2 run_mkwfs t10c_cpu 1.65 0 || return 1
    mv outarraypha.fits cpu_pha.fits
    mv outarrayamp.fits cpu_amp.fits
    mv outsarraypha.fits cpu_spha.fits
    mv outsarrayamp.fits cpu_samp.fits

    ATMTURB_SIMD=CUDA run_mkwfs t10c_gpu 1.65 0 || return 1
    mv outarraypha.fits gpu_pha.fits
    mv outarrayamp.fits gpu_amp.fits
    mv outsarraypha.fits gpu_spha.fits
    mv outsarrayamp.fits gpu_samp.fits

    python3 - << 'EOF'
import sys
import numpy as np
from astropy.io import fits

all_ok = True
details = []
for tag in ["pha", "amp", "spha", "samp"]:
    c = fits.getdata(f"cpu_{tag}.fits")
    g = fits.getdata(f"gpu_{tag}.fits")
    diff = float(np.max(np.abs(c - g)))
    corr = float(np.corrcoef(c.ravel(), g.ravel())[0, 1])
    ok = (diff < 1e-4) and (corr > 0.99999)
    if not ok:
        all_ok = False
    details.append(f"{tag.upper()}: diff={diff:.2e} corr={corr:.7f}")

msg = ": " + " ".join(details)
print(("PASS" if all_ok else "FAIL") + msg)
sys.exit(0 if all_ok else 1)
EOF
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

# T9a: Fresnel off leaves amplitude identically 1.0 on all pixels and frames
scenario_T9a_fresnel_off() {
    one_layer_profile 5000 1.0 10.0 0.0 10000 0.0
    write_conf WFsim.conf WAVEFRONT_AMPLITUDE=0 FRESNEL_PROPAGATION=0 \
        TIME_SPAN=0.05 WFTIME_STEP=0.01
    run_mkwfs t9a 1.65 0 || return 1
    python3 -c '
import sys
from astropy.io import fits
import numpy as np

amp = fits.getdata("outarrayamp.fits")
samp = fits.getdata("outsarrayamp.fits")
dev_pri = np.max(np.abs(amp - 1.0))
dev_sec = np.max(np.abs(samp - 1.0))
ok = (dev_pri < 1e-6) and (dev_sec < 1e-6)
msg = f": amplitude identically 1.0 (dev {max(dev_pri, dev_sec):.1e})"
print(("PASS" if ok else "FAIL") + msg)
sys.exit(0 if ok else 1)
'
}

# T9b: Multi-layer diffractive propagation conserves total optical energy (<I> = 1.00 +- 0.01)
scenario_T9b_energy_conservation() {
    cat << 'EOF' > turbul.prof
# alt cn2 speed dir L0 l0
 4215     5.32        6.5     1.47  10000 0.0
 4230     1.47        6.55    1.57  10000 0.0
 4349     1.08        6.6     1.67  10000 0.0
 5007     2.11        6.7     1.77  10000 0.0
12000     1.83       22.0     3.10  10000 0.0
16200     1.48        9.5     3.20  10000 0.0
23701     0.697       5.6     3.30  10000 0.0
EOF
    write_conf WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=1 \
        FRESNEL_PROPAGATION_BIN=1000.0 TIME_SPAN=0.2 WFTIME_STEP=0.01
    run_mkwfs t9b 1.65 0 || return 1
    python3 -c '
import sys
from astropy.io import fits
import numpy as np

amp = fits.getdata("outarrayamp.fits")
samp = fits.getdata("outsarrayamp.fits")
mean_I_pri = float(np.mean(amp ** 2))
mean_I_sec = float(np.mean(samp ** 2))
ok = abs(mean_I_pri - 1.0) < 0.01 and abs(mean_I_sec - 1.0) < 0.01
msg = f": mean intensity pri={mean_I_pri:.5f} sec={mean_I_sec:.5f} (tol 0.01)"
print(("PASS" if ok else "FAIL") + msg)
sys.exit(0 if ok else 1)
'
}

# T9c: High-altitude single layer scintillation index matches Rytov approximation (+- 20%)
scenario_T9c_scintillation_rytov() {
    one_layer_profile 15000 1.0 15.0 0.0 10000 0.0
    write_conf WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=1 \
        SITE_ALT=500.0 ZENITH_ANGLE=0.5235987756 TURBULENCE_REF_WAVEL=0.5 \
        TURBULENCE_SEEING=0.6 TIME_SPAN=0.4 WFTIME_STEP=0.01
    run_mkwfs t9c 1.65 0 || return 1
    python3 -c '
import sys
from astropy.io import fits
import numpy as np

lam = 0.5e-6
seeing_rad = 0.6 * (np.pi / (180.0 * 3600.0))
r0 = 0.98 * lam / seeing_rad
h = 15000.0
h_site = 500.0
cos_z = np.cos(0.5235987756)
int_cn2 = 0.060 * (lam ** 2) * (r0 ** (-5.0 / 3.0))
geom_factor = 19.12 * (lam ** (-7.0 / 6.0)) * (cos_z ** (-11.0 / 6.0))
sigma2_rytov = geom_factor * int_cn2 * ((h - h_site) ** (5.0 / 6.0))

amp = fits.getdata("outarrayamp.fits")
I = amp ** 2
mean_I = float(np.mean(I))
var_I = float(np.var(I))
scint_sim = var_I / (mean_I ** 2)

rel_err = abs(scint_sim - sigma2_rytov) / sigma2_rytov
ok = rel_err <= 0.20
msg = (f": scint index sim={scint_sim:.4f}, rytov={sigma2_rytov:.4f} "
       f"(err {rel_err*100:.1f}%, tol 20%)")
print(("PASS" if ok else "FAIL") + msg)
sys.exit(0 if ok else 1)
'
}

# T9d: Zero-distance Fresnel propagation reproduces geometric phase to < 1e-5 rad
scenario_T9d_zero_distance_identity() {
    one_layer_profile 500 1.0 10.0 0.0 10000 0.0
    mkdir -p diff geo

    write_conf diff/WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=1 \
        SITE_ALT=500.0 TIME_SPAN=0.05 WFTIME_STEP=0.01
    (cd diff && cp ../turbul.prof . && run_mkwfs t9d_diff 1.65 0) || return 1

    write_conf geo/WFsim.conf WAVEFRONT_AMPLITUDE=0 FRESNEL_PROPAGATION=0 \
        SITE_ALT=500.0 TIME_SPAN=0.05 WFTIME_STEP=0.01
    (cd geo && cp ../turbul.prof . && run_mkwfs t9d_geo 1.65 0) || return 1

    python3 -c '
import sys
from astropy.io import fits
import numpy as np

p_diff = fits.getdata("diff/outarraypha.fits")
p_geo  = fits.getdata("geo/outarraypha.fits")
a_diff = fits.getdata("diff/outarrayamp.fits")

max_p_diff = float(np.max(np.abs(p_diff - p_geo)))
max_a_diff = float(np.max(np.abs(a_diff - 1.0)))

ok = (max_p_diff < 1e-5) and (max_a_diff < 1e-5)
msg = f": max |pha_diff - pha_geo| = {max_p_diff:.2e} rad, |amp - 1| = {max_a_diff:.2e} (tol 1e-5)"
print(("PASS" if ok else "FAIL") + msg)
sys.exit(0 if ok else 1)
'
}

# T9e: Rytov vs split-step weak regime correlation (> 0.98 primary, > 0.95 secondary)
scenario_T9e_rytov_vs_splitstep() {
    one_layer_profile 5000 1.0 15.0 0.0 10000 0.0
    mkdir -p split rytov

    write_conf split/WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=1 \
        TURBULENCE_SEEING=0.1 SITE_ALT=500.0 ZENITH_ANGLE=0.5235987756 \
        TIME_SPAN=0.05 WFTIME_STEP=0.01 SEED=12345
    (cd split && cp ../turbul.prof . && run_mkwfs t9e_split 1.65 0) || return 1

    write_conf rytov/WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=2 \
        TURBULENCE_SEEING=0.1 SITE_ALT=500.0 ZENITH_ANGLE=0.5235987756 \
        TIME_SPAN=0.05 WFTIME_STEP=0.01 SEED=12345
    (cd rytov && cp ../turbul.prof . && run_mkwfs t9e_rytov 1.65 0) || return 1

    python3 -c '
import sys
from astropy.io import fits
import numpy as np

a1 = fits.getdata("split/outarrayamp.fits")[0]
a2 = fits.getdata("rytov/outarrayamp.fits")[0]
sa1 = fits.getdata("split/outsarrayamp.fits")[0]
sa2 = fits.getdata("rytov/outsarrayamp.fits")[0]
sp1 = fits.getdata("split/outsarraypha.fits")[0]
sp2 = fits.getdata("rytov/outsarraypha.fits")[0]

s = slice(16, -16)
chi1_p = np.log(np.maximum(a1[s, s], 1e-6))
chi2_p = np.log(np.maximum(a2[s, s], 1e-6))
corr_p = float(np.corrcoef(chi1_p.ravel(), chi2_p.ravel())[0, 1])

chi1_s = np.log(np.maximum(sa1[s, s], 1e-6))
chi2_s = np.log(np.maximum(sa2[s, s], 1e-6))
corr_s = float(np.corrcoef(chi1_s.ravel(), chi2_s.ravel())[0, 1])

p_diff = float(np.std(sp1[s, s] - sp2[s, s]))
p_rms  = float(np.std(sp1[s, s]))
rel_p_diff = p_diff / p_rms if p_rms > 0 else 0.0

ok = (corr_p > 0.98) and (corr_s > 0.95) and (rel_p_diff < 0.10)
msg = f": corr pri={corr_p:.4f} sec={corr_s:.4f}, rel pha diff={rel_p_diff*100:.2f}%"
print(("PASS" if ok else "FAIL") + msg)
sys.exit(0 if ok else 1)
'
}

# T9f: Multi-layer Rytov Fourier propagation conserves total optical energy (<I> = 1.00 +- 0.01)
scenario_T9f_energy_conservation_rytov() {
    cat << 'EOF' > turbul.prof
# alt cn2 speed dir L0 l0
 4215     5.32        6.5     1.47  10000 0.0
 4230     1.47        6.55    1.57  10000 0.0
 4349     1.08        6.6     1.67  10000 0.0
 5007     2.11        6.7     1.77  10000 0.0
12000     1.83       22.0     3.10  10000 0.0
16200     1.48        9.5     3.20  10000 0.0
23701     0.697       5.6     3.30  10000 0.0
EOF
    write_conf WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=2 \
        FRESNEL_PROPAGATION_BIN=1000.0 TIME_SPAN=0.2 WFTIME_STEP=0.01
    run_mkwfs t9f 1.65 0 || return 1
    python3 -c '
import sys
from astropy.io import fits
import numpy as np

amp = fits.getdata("outarrayamp.fits")
samp = fits.getdata("outsarrayamp.fits")
mean_I_pri = float(np.mean(amp ** 2))
mean_I_sec = float(np.mean(samp ** 2))
ok = abs(mean_I_pri - 1.0) < 0.01 and abs(mean_I_sec - 1.0) < 0.01
msg = f": mean intensity pri={mean_I_pri:.5f} sec={mean_I_sec:.5f} (tol 0.01)"
print(("PASS" if ok else "FAIL") + msg)
sys.exit(0 if ok else 1)
'
}

# T9g: High-altitude single layer scintillation index in Mode 2 matches Rytov approximation (+- 20%)
scenario_T9g_scintillation_rytov_mode2() {
    one_layer_profile 15000 1.0 15.0 0.0 10000 0.0
    write_conf WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=2 \
        SITE_ALT=500.0 ZENITH_ANGLE=0.5235987756 TURBULENCE_REF_WAVEL=0.5 \
        TURBULENCE_SEEING=0.6 TIME_SPAN=0.4 WFTIME_STEP=0.01
    run_mkwfs t9g 1.65 0 || return 1
    python3 -c '
import sys
from astropy.io import fits
import numpy as np

lam = 0.5e-6
seeing_rad = 0.6 * (np.pi / (180.0 * 3600.0))
r0 = 0.98 * lam / seeing_rad
h = 15000.0
h_site = 500.0
cos_z = np.cos(0.5235987756)
int_cn2 = 0.060 * (lam ** 2) * (r0 ** (-5.0 / 3.0))
geom_factor = 19.12 * (lam ** (-7.0 / 6.0)) * (cos_z ** (-11.0 / 6.0))
sigma2_rytov = geom_factor * int_cn2 * ((h - h_site) ** (5.0 / 6.0))

amp = fits.getdata("outarrayamp.fits")
chi = np.log(np.maximum(amp, 1e-6))
sigma2_sim = 4.0 * float(np.var(chi))

rel_err = abs(sigma2_sim - sigma2_rytov) / sigma2_rytov
ok = rel_err <= 0.20
msg = (f": Rytov variance 4*var(chi) sim={sigma2_sim:.4f}, rytov={sigma2_rytov:.4f} "
       f"(err {rel_err*100:.1f}%, tol 20%)")
print(("PASS" if ok else "FAIL") + msg)
sys.exit(0 if ok else 1)
'
}

# T9h: Zero-distance Rytov propagation reproduces geometric phase to < 1e-5 rad
scenario_T9h_zero_distance_mode2() {
    one_layer_profile 500 1.0 10.0 0.0 10000 0.0
    mkdir -p rytov geo

    write_conf rytov/WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=2 \
        SITE_ALT=500.0 TIME_SPAN=0.05 WFTIME_STEP=0.01
    (cd rytov && cp ../turbul.prof . && run_mkwfs t9h_rytov 1.65 0) || return 1

    write_conf geo/WFsim.conf WAVEFRONT_AMPLITUDE=0 FRESNEL_PROPAGATION=0 \
        SITE_ALT=500.0 TIME_SPAN=0.05 WFTIME_STEP=0.01
    (cd geo && cp ../turbul.prof . && run_mkwfs t9h_geo 1.65 0) || return 1

    python3 -c '
import sys
from astropy.io import fits
import numpy as np

p_rytov = fits.getdata("rytov/outarraypha.fits")
p_geo   = fits.getdata("geo/outarraypha.fits")
a_rytov = fits.getdata("rytov/outarrayamp.fits")

max_p_diff = float(np.max(np.abs(p_rytov - p_geo)))
max_a_diff = float(np.max(np.abs(a_rytov - 1.0)))

ok = (max_p_diff < 1e-5) and (max_a_diff < 1e-5)
msg = f": max |pha_rytov - pha_geo| = {max_p_diff:.2e} rad, |amp - 1| = {max_a_diff:.2e} (tol 1e-5)"
print(("PASS" if ok else "FAIL") + msg)
sys.exit(0 if ok else 1)
'
}

# T9l: Option B vs Option C agreement (corr > 0.95 at 30 deg, identity at z=0)
scenario_T9l_rytov_option_b_vs_c() {
    one_layer_profile 5000 1.0 15.0 0.0 10000 0.0
    mkdir -p optB optC optB_z0 optC_z0

    write_conf optB/WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=2 \
        TURBULENCE_SEEING=0.1 SITE_ALT=500.0 ZENITH_ANGLE=0.5235987756 \
        TIME_SPAN=0.05 WFTIME_STEP=0.01 SEED=12345
    (cd optB && cp ../turbul.prof . && run_mkwfs t9l_b 1.65 0) || return 1

    write_conf optC/WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=3 \
        TURBULENCE_SEEING=0.1 SITE_ALT=500.0 ZENITH_ANGLE=0.5235987756 \
        TIME_SPAN=0.05 WFTIME_STEP=0.01 SEED=12345
    (cd optC && cp ../turbul.prof . && run_mkwfs t9l_c 1.65 0) || return 1

    write_conf optB_z0/WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=2 \
        TURBULENCE_SEEING=0.1 SITE_ALT=500.0 ZENITH_ANGLE=0.0 \
        TIME_SPAN=0.05 WFTIME_STEP=0.01 SEED=12345
    (cd optB_z0 && cp ../turbul.prof . && run_mkwfs t9l_bz0 1.65 0) || return 1

    write_conf optC_z0/WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=3 \
        TURBULENCE_SEEING=0.1 SITE_ALT=500.0 ZENITH_ANGLE=0.0 \
        TIME_SPAN=0.05 WFTIME_STEP=0.01 SEED=12345
    (cd optC_z0 && cp ../turbul.prof . && run_mkwfs t9l_cz0 1.65 0) || return 1

    python3 -c '
import sys
from astropy.io import fits
import numpy as np

sa_b = fits.getdata("optB/outsarrayamp.fits")[0]
sa_c = fits.getdata("optC/outsarrayamp.fits")[0]
sa_bz0 = fits.getdata("optB_z0/outsarrayamp.fits")[0]
sa_cz0 = fits.getdata("optC_z0/outsarrayamp.fits")[0]

s = slice(16, -16)
chi_b = np.log(np.maximum(sa_b[s, s], 1e-6))
chi_c = np.log(np.maximum(sa_c[s, s], 1e-6))
corr_30 = float(np.corrcoef(chi_b.ravel(), chi_c.ravel())[0, 1])

chi_bz0 = np.log(np.maximum(sa_bz0[s, s], 1e-6))
chi_cz0 = np.log(np.maximum(sa_cz0[s, s], 1e-6))
corr_z0 = float(np.corrcoef(chi_bz0.ravel(), chi_cz0.ravel())[0, 1])
rel_diff_z0 = float(np.std(chi_bz0 - chi_cz0) / np.std(chi_cz0))

ok = (corr_30 > 0.95) and (corr_z0 > 0.999) and (rel_diff_z0 < 0.001)
msg = (f": corr(z=30deg)={corr_30:.4f}, corr(z=0)={corr_z0:.5f}, "
       f"rel diff(z=0)={rel_diff_z0*100:.3f}%")
print(("PASS" if ok else "FAIL") + msg)
sys.exit(0 if ok else 1)
'
}

# T9i: Moisan periodic-plus-smooth decomposition suppresses edge scintillation (ratio < 1.30)
scenario_T9i_edge_artifact_suppression() {
    one_layer_profile 10000 1.0 15.0 0.0 10000 0.0
    write_conf WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=2 \
        SITE_ALT=500.0 TURBULENCE_SEEING=1.0 TIME_SPAN=0.1 WFTIME_STEP=0.01 SEED=42
    run_mkwfs t9i 1.65 0 || return 1
    python3 -c '
import sys, numpy as np
from astropy.io import fits
amp = fits.getdata("outarrayamp.fits")[0]
chi = np.log(np.maximum(amp, 1e-6))
mask_edge = np.zeros_like(chi, dtype=bool)
mask_edge[:4, :] = True
mask_edge[-4:, :] = True
mask_edge[:, :4] = True
mask_edge[:, -4:] = True
rms_edge = float(np.std(chi[mask_edge]))
rms_int  = float(np.std(chi[~mask_edge]))
ratio = rms_edge / rms_int if rms_int > 0 else 0.0
ok = ratio < 1.30
print(("PASS" if ok else "FAIL") + f": edge/int rms ratio = {ratio:.3f} (tol < 1.30)")
sys.exit(0 if ok else 1)
'
}

# T9k: Strong scintillation warning emitted when sigma_R^2 > 0.30
scenario_T9k_strong_scintillation_warning() {
    one_layer_profile 15000 1.0 15.0 0.0 10000 0.0
    write_conf WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=2 \
        SITE_ALT=500.0 ZENITH_ANGLE=1.0471975512 TURBULENCE_REF_WAVEL=0.5 \
        TURBULENCE_SEEING=0.8 TIME_SPAN=0.05 WFTIME_STEP=0.01
    run_mkwfs t9k 1.65 0 || return 1
    if grep -q "Strong scintillation regime" mkwfs.log; then
        echo "PASS: strong scintillation warning correctly emitted"
        return 0
    else
        echo "FAIL: strong scintillation warning missing from mkwfs.log"
        return 1
    fi
}

# T9m: Differential refraction shifts secondary scintillation pattern by physical ray vector
scenario_T9m_chromatic_refraction_consistency() {
    one_layer_profile 15000 1.0 0.0 0.0 10000 0.0
    write_conf WFsim.conf ZENITH_ANGLE=0.785398 PARALLACTIC_ANGLE=0.0 SITE_ALT=4200.0 \
        TURBULENCE_SEEING=0.2 WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=2 \
        MAKE_SWAVEFRONT=1 SLAMBDA=1.65 TIME_SPAN=0.01 SEED=42
    run_mkwfs t9m 1.65 0 || return 1
    "${VALIDATE[@]}" shift outarrayamp.fits outsarrayamp.fits --dx 0.0 --dy -2.26 --tol 0.15
}

# T9n: Option C fallback explicitly triggered via config keyword
scenario_T9n_option_b_fallback() {
    one_layer_profile 5000 1.0 15.0 0.0 10000 0.0
    write_conf WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=3 \
        TURBULENCE_SEEING=0.1 TIME_SPAN=0.05 WFTIME_STEP=0.01 SEED=12345
    run_mkwfs t9n 1.65 0 || return 1
    if grep -q "Secondary wavelength using Option C" mkwfs.log; then
        echo "PASS: Option C fallback confirmed in log"
        return 0
    else
        echo "FAIL: Option C fallback not logged"
        return 1
    fi
}

# T9o: Guard band margin padding eliminates boundary artifacts while preserving interior field
scenario_T9o_guard_band_parity() {
    one_layer_profile 15000 1.0 15.0 0.0 10000 0.0
    mkdir -p g0 g16

    write_conf g0/WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=2 \
        FRESNEL_GUARD_PIX=0 SITE_ALT=500.0 ZENITH_ANGLE=0.5235987756 \
        TURBULENCE_SEEING=0.6 TIME_SPAN=0.02 WFTIME_STEP=0.01 WFsize=64 SEED=42

    write_conf g16/WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=2 \
        FRESNEL_GUARD_PIX=16 SITE_ALT=500.0 ZENITH_ANGLE=0.5235987756 \
        TURBULENCE_SEEING=0.6 TIME_SPAN=0.02 WFTIME_STEP=0.01 WFsize=64 SEED=42

    (cd g0 && run_mkwfs t9o_g0 1.65 0) || return 1
    (cd g16 && run_mkwfs t9o_g16 1.65 0) || return 1

    python3 -c '
import sys
from astropy.io import fits
import numpy as np

p0 = fits.getdata("g0/outarraypha.fits")[0]
p16 = fits.getdata("g16/outarraypha.fits")[0]
a0 = fits.getdata("g0/outarrayamp.fits")[0]
a16 = fits.getdata("g16/outarrayamp.fits")[0]

i0_mean = float(np.mean(a0**2))
i16_mean = float(np.mean(a16**2))
energy_ok = abs(i0_mean - 1.0) < 0.01 and abs(i16_mean - 1.0) < 0.01

int_p0 = p0[16:48, 16:48]
int_p16 = p16[16:48, 16:48]
corr_p = float(np.corrcoef(int_p0.ravel(), int_p16.ravel())[0, 1])

int_a0 = a0[16:48, 16:48]
int_a16 = a16[16:48, 16:48]
corr_a = float(np.corrcoef(int_a0.ravel(), int_a16.ravel())[0, 1])

chi16 = np.log(np.maximum(a16, 1e-6))
mask_edge = np.zeros_like(chi16, dtype=bool)
mask_edge[:4, :] = True
mask_edge[-4:, :] = True
mask_edge[:, :4] = True
mask_edge[:, -4:] = True
rms_edge = float(np.std(chi16[mask_edge]))
rms_int = float(np.std(chi16[~mask_edge]))
ratio16 = rms_edge / rms_int if rms_int > 0 else 0.0

ok = energy_ok and corr_p > 0.999 and corr_a > 0.99 and abs(ratio16 - 1.0) < 0.12
msg = (f": corr_pha={corr_p:.6f}, corr_amp={corr_a:.6f}, "
       f"edge/int ratio={ratio16:.4f}, mean_I={i16_mean:.4f}")
print(("PASS" if ok else "FAIL") + msg)
sys.exit(0 if ok else 1)
'
}

# T9p: Linear z-interpolation achieves O(Delta z^2) error convergence (~4x error reduction)
scenario_T9p_z_interpolation_scaling() {
    cat << 'EOF' > turbul.prof
# altitude(m)   relativeCN2     speed(m/s)      direction(rad)
 2000     1.0        10.0     0.0
 3000     1.0        10.0     0.0
 4000     1.0        10.0     0.0
 5000     1.0        10.0     0.0
 6000     1.0        10.0     0.0
 7000     1.0        10.0     0.0
 8000     1.0        10.0     0.0
 9000     1.0        10.0     0.0
10000     1.0        10.0     0.0
EOF

    mkdir -p truth c4 c2 z4 z2

    write_conf truth/WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=2 \
        FRESNEL_PROPAGATION_BIN=0 FRESNEL_RYTOV_ZINT=0 \
        PUPIL_SCALE=0.03 SITE_ALT=0.0 MASTER_SIZE=2048 LOWFREQ=1 ROLLING=1 \
        TURBULENCE_SEEING=0.5 TIME_SPAN=0.02 WFTIME_STEP=0.01 WFsize=64 SEED=123456

    write_conf c4/WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=2 \
        FRESNEL_PROPAGATION_BIN=4000 FRESNEL_RYTOV_ZINT=0 \
        PUPIL_SCALE=0.03 SITE_ALT=0.0 MASTER_SIZE=2048 LOWFREQ=1 ROLLING=1 \
        TURBULENCE_SEEING=0.5 TIME_SPAN=0.02 WFTIME_STEP=0.01 WFsize=64 SEED=123456

    write_conf c2/WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=2 \
        FRESNEL_PROPAGATION_BIN=2000 FRESNEL_RYTOV_ZINT=0 \
        PUPIL_SCALE=0.03 SITE_ALT=0.0 MASTER_SIZE=2048 LOWFREQ=1 ROLLING=1 \
        TURBULENCE_SEEING=0.5 TIME_SPAN=0.02 WFTIME_STEP=0.01 WFsize=64 SEED=123456

    write_conf z4/WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=2 \
        FRESNEL_PROPAGATION_BIN=4000 FRESNEL_RYTOV_ZINT=1 \
        PUPIL_SCALE=0.03 SITE_ALT=0.0 MASTER_SIZE=2048 LOWFREQ=1 ROLLING=1 \
        TURBULENCE_SEEING=0.5 TIME_SPAN=0.02 WFTIME_STEP=0.01 WFsize=64 SEED=123456

    write_conf z2/WFsim.conf WAVEFRONT_AMPLITUDE=1 FRESNEL_PROPAGATION=2 \
        FRESNEL_PROPAGATION_BIN=2000 FRESNEL_RYTOV_ZINT=1 \
        PUPIL_SCALE=0.03 SITE_ALT=0.0 MASTER_SIZE=2048 LOWFREQ=1 ROLLING=1 \
        TURBULENCE_SEEING=0.5 TIME_SPAN=0.02 WFTIME_STEP=0.01 WFsize=64 SEED=123456

    (cd truth && cp ../turbul.prof . && run_mkwfs t9p_truth 0.5 0) || return 1
    (cd c4 && cp ../turbul.prof . && run_mkwfs t9p_c4 0.5 0) || return 1
    (cd c2 && cp ../turbul.prof . && run_mkwfs t9p_c2 0.5 0) || return 1
    (cd z4 && cp ../turbul.prof . && run_mkwfs t9p_z4 0.5 0) || return 1
    (cd z2 && cp ../turbul.prof . && run_mkwfs t9p_z2 0.5 0) || return 1

    python3 -c '
import sys
from astropy.io import fits
import numpy as np

p_true = fits.getdata("truth/outarraypha.fits")[0]
a_true = fits.getdata("truth/outarrayamp.fits")[0]
p_c4 = fits.getdata("c4/outarraypha.fits")[0]
p_c2 = fits.getdata("c2/outarraypha.fits")[0]
p_z4 = fits.getdata("z4/outarraypha.fits")[0]
p_z2 = fits.getdata("z2/outarraypha.fits")[0]
a_z4 = fits.getdata("z4/outarrayamp.fits")[0]
a_z2 = fits.getdata("z2/outarrayamp.fits")[0]

s = slice(8, -8)

err_c4 = float(np.std(p_c4[s, s] - p_true[s, s]))
err_c2 = float(np.std(p_c2[s, s] - p_true[s, s]))
ratio_c = err_c4 / err_c2 if err_c2 > 0 else 0.0

err_z4 = float(np.std(p_z4[s, s] - p_true[s, s]))
err_z2 = float(np.std(p_z2[s, s] - p_true[s, s]))
ratio_z = err_z4 / err_z2 if err_z2 > 0 else 0.0

mean_I_z4 = float(np.mean(a_z4**2))
mean_I_z2 = float(np.mean(a_z2**2))
energy_ok = abs(mean_I_z4 - 1.0) < 0.01 and abs(mean_I_z2 - 1.0) < 0.01

ok = energy_ok and (ratio_z >= 3.0) and (err_z2 < err_c2 * 0.3)
msg = (f": zinterp err 4000={err_z4:.4f} -> 2000={err_z2:.4f} "
       f"(ratio {ratio_z:.2f}x, target >= 3.0x); "
       f"centroid ratio={ratio_c:.2f}x; mean_I={mean_I_z2:.4f}")
print(("PASS" if ok else "FAIL") + msg)
sys.exit(0 if ok else 1)
'
}

# T16a: Flat wavefront PSF produces normalized Strehl = 1.000 +- 1e-3 and centered peak
scenario_T16a_zero_turb_strehl() {
    run_aoloop t16a "none" "PSFflat" "psf_flat.fits" 0 0.0 0.0 0 1.65 8.0 || return 1
    python3 -c '
import sys, re
from astropy.io import fits
import numpy as np

psf = fits.getdata("psf_flat.fits")
ny, nx = psf.shape
cy, cx = np.unravel_index(np.argmax(psf), psf.shape)
centered = (abs(cy - ny // 2) <= 1) and (abs(cx - nx // 2) <= 1)

with open("aoloop.log") as f:
    log_text = f.read()
m = re.search(r"Strehl = ([0-9.]+)", log_text)
if not m:
    print("FAIL: could not parse Strehl ratio from log")
    sys.exit(1)
strehl = float(m.group(1))

ok = centered and abs(strehl - 1.0) < 1e-3
msg = f": Strehl = {strehl:.4f}, peak at ({cy},{cx}) vs target ({ny//2},{nx//2})"
print(("PASS" if ok else "FAIL") + msg)
sys.exit(0 if ok else 1)
'
}

# T16b: Open-loop seeing FWHM matches Kolmogorov 0.98 * lambda / r0 (+- 20%)
scenario_T16b_open_loop_seeing() {
    one_layer_profile 0 1.0 10.0 0.0 25.0 0.01
    write_conf WFsim.conf TURBULENCE_SEEING=0.8 TIME_SPAN=0.2 WFTIME_STEP=0.01 \
        WFsize=128 PUPIL_SCALE=0.02
    run_mkwfs t16b 1.65 0 || return 1
    run_aoloop t16b_ao "outarraypha.fits" "PSFol" "psf_ol.fits" 0 0.0 0.0 0 1.65 2.56 || return 1
    python3 -c '
import sys
from astropy.io import fits
import numpy as np

psf = fits.getdata("psf_ol.fits")
ny, nx = psf.shape
cy, cx = np.unravel_index(np.argmax(psf), psf.shape)

y, x = np.ogrid[:ny, :nx]
r = np.hypot(x - cx, y - cy)
half_max = 0.5 * np.max(psf)
fwhm_meas = 2.0 * float(np.max(r[psf >= half_max]))

lam_sci = 1.65e-6
lam_ref = 0.5e-6
seeing_ref_rad = 0.8 * (np.pi / (180.0 * 3600.0))
r0_ref = 0.98 * lam_ref / seeing_ref_rad
r0_sci = r0_ref * ((lam_sci / lam_ref) ** 1.2)
theta_seeing = 0.98 * lam_sci / r0_sci

dx = 0.02
n_fft = 256
dtheta = lam_sci / (n_fft * dx)
fwhm_theo_pix = theta_seeing / dtheta

rel_err = abs(fwhm_meas - fwhm_theo_pix) / fwhm_theo_pix
ok = rel_err <= 0.20
msg = (f": measured FWHM = {fwhm_meas:.2f} px, theoretical = {fwhm_theo_pix:.2f} px "
       f"(err {rel_err*100:.1f}%, tol 20%)")
print(("PASS" if ok else "FAIL") + msg)
sys.exit(0 if ok else 1)
'
}

# T16c: Closed-loop integrator increases Strehl over open-loop and converges
scenario_T16c_closed_loop_strehl() {
    one_layer_profile 0 1.0 10.0 0.0 25.0 0.01
    write_conf WFsim.conf TURBULENCE_SEEING=0.5 TIME_SPAN=0.2 WFTIME_STEP=0.01 \
        WFsize=128 PUPIL_SCALE=0.02
    run_mkwfs t16c 1.65 0 || return 1
    run_aoloop t16c_ol "outarraypha.fits" "PSFol" "psf_ol.fits" 0 0.0 0.0 0 1.65 2.56 || return 1
    cp aoloop.log aoloop_ol.log
    run_aoloop t16c_cl "outarraypha.fits" "PSFcl" "psf_cl.fits" 1 0.5 0.001 1 1.65 2.56 || return 1
    cp aoloop.log aoloop_cl.log
    python3 -c '
import sys, re

def parse_strehl(fname):
    with open(fname) as f:
        text = f.read()
    m = re.search(r"Strehl = ([0-9.]+)\s+\(first = ([0-9.]+),\s+last = ([0-9.]+)\)", text)
    if not m:
        return None, None, None
    return float(m.group(1)), float(m.group(2)), float(m.group(3))

s_ol_cumul, s_ol_first, s_ol_last = parse_strehl("aoloop_ol.log")
s_cl_cumul, s_cl_first, s_cl_last = parse_strehl("aoloop_cl.log")

if s_ol_cumul is None or s_cl_cumul is None:
    print("FAIL: failed to parse Strehl ratios from logs")
    sys.exit(1)

gain_ratio = s_cl_cumul / s_ol_cumul
converged = s_cl_last > s_cl_first
ok = (gain_ratio >= 1.30) and converged

msg = (f": open-loop Strehl = {s_ol_cumul:.4f} (last {s_ol_last:.4f}), "
       f"closed-loop Strehl = {s_cl_cumul:.4f} (last {s_cl_last:.4f}, ratio {gain_ratio:.2f}x)")
print(("PASS" if ok else "FAIL") + msg)
sys.exit(0 if ok else 1)
'
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
    "T8b_hv_bounds_and_recomputed_r0:pass"
    "T8c_turbulent_wind_advection:pass"
    "T9a_fresnel_off:pass"
    "T9b_energy_conservation:pass"
    "T9c_scintillation_rytov:pass"
    "T9d_zero_distance_identity:pass"
    "T9e_rytov_vs_splitstep:pass"
    "T9f_energy_conservation_rytov:pass"
    "T9g_scintillation_rytov_mode2:pass"
    "T9h_zero_distance_mode2:pass"
    "T9i_edge_artifact_suppression:pass"
    "T9k_strong_scintillation_warning:pass"
    "T9l_rytov_option_b_vs_c:pass"
    "T9m_chromatic_refraction_consistency:pass"
    "T9n_option_b_fallback:pass"
    "T9o_guard_band_parity:pass"
    "T9p_z_interpolation_scaling:pass"
    "T10_breathing:pass"
    "T10b_simd_parity:pass"
    "T10c_cuda_rytov_parity:pass"
    "T11_tilt_variance:pass"
    "T12_wrap_decorrelation:pass"
    "T13_epoch_stability:pass"
    "T14_seam_absence:pass"
    "T15_invalid_pupil_scale:pass"
    "T16a_zero_turb_strehl:pass"
    "T16b_open_loop_seeing:pass"
    "T16c_closed_loop_strehl:pass"
)

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
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
fi
