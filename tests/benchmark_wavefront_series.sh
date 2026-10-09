#!/usr/bin/env bash
#
# tests/benchmark_wavefront_series.sh
# End-to-end runtime benchmark for synthesizing wavefront series with mkwfs.
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Locate milk-fpsexec-atmturb-mkwfs
MKWFS_EXEC="${MILK_MKWFS_EXEC:-}"
if [[ -z "$MKWFS_EXEC" ]]; then
    CANDIDATES=(
        "$REPO_ROOT/_build/milk-fpsexec-atmturb-mkwfs"
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
    exit 1
fi

TEST_TMPDIR="$(mktemp -d /tmp/bench_mkwfs_XXXXXX)"
MILK_SHM_DIR="${MILK_SHM_DIR:-/milk/shm}"
trap 'rm -rf "$TEST_TMPDIR"; rm -f "$MILK_SHM_DIR"/twfsbench_*.fps.shm' EXIT

cd "$TEST_TMPDIR"

# Generate 5-layer Hufnagel-Valley profile
cat > turbul.prof << 'EOF'
# altitude(m)   relativeCN2     speed(m/s)      direction(rad)
200.0           5.0             8.0             0.2
1500.0          2.5             12.0            0.8
5000.0          1.8             22.0            1.5
9000.0          1.0             32.0            2.1
14000.0         0.5             18.0            2.9
EOF

echo "########################################################################"
echo "#  milkatmturb End-to-End Wavefront Series Benchmark (mkwfs)           #"
echo "########################################################################"
echo " Executable: $MKWFS_EXEC"
echo ""

run_series_benchmark() {
    local name="$1"
    local wfsize="$2"
    local msize="$3"
    local nbframes="$4"
    local fresnel="$5"

    local tstep="0.005"
    local tspan=$(awk "BEGIN { printf \"%.6f\", $nbframes * $tstep }")

    cat > WFsim.conf << EOF
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
SWF_WRITE2DISK             0
SWF_FILE_PREFIX            ./swf_
SHM_SOUTPUT                0
SHM_SPREFIX                shmswf
SHM_SOUTPUTM               0

WFsize                     $wfsize
PUPIL_SCALE                0.020000

REALTIME                   0
REALTIMEFACTOR             1.0
WFTIME_STEP                $tstep
TIME_SPAN                  $tspan
NB_TSPAN                   1
SIMTDELAY                  0
WAITFORSEM                 0
WAITSEMIMNAME              

SKIP_EXISTING              0
WF_RAW_SIZE                $wfsize
MASTER_SIZE                $msize
WAVEFRONT_AMPLITUDE        $([[ "$fresnel" -gt 0 ]] && echo 1 || echo 0)
FRESNEL_PROPAGATION        $fresnel
FRESNEL_PROPAGATION_BIN    1000.0
EOF

    echo "========================================================================"
    echo " Scenario: $name (${wfsize}x${wfsize} pupil, $nbframes frames, Fresnel=$fresnel)"
    echo "========================================================================"
    printf "%-20s | %10s | %12s | %12s | %8s\n" \
           "Backend" "Wall Time" "Frame Rate" "Throughput" "Speedup"
    echo "---------------------+------------+--------------+--------------+---------"

    local backends=("SCALAR" "AVX2")
    if grep -qE "avx512f" /proc/cpuinfo 2>/dev/null; then
        backends+=("AVX512")
    fi
    if [[ "$fresnel" -eq 0 ]] && command -v nvidia-smi >/dev/null 2>&1 && \
       nvidia-smi >/dev/null 2>&1; then
        backends+=("CUDA")
    fi

    local t_scalar=0.0
    for backend in "${backends[@]}"; do
        local fps_id="twfsbench_${wfsize}_${backend}"
        local start_ns=$(date +%s%N)
        local out
        out=$(ATMTURB_SIMD="$backend" "$MKWFS_EXEC" -n "$fps_id" exec 0.7 0 2>&1)
        local end_ns=$(date +%s%N)
        rm -f "$MILK_SHM_DIR/${fps_id}.fps.shm"

        local dur_s=$(awk "BEGIN { printf \"%.3f\", ($end_ns - $start_ns) / 1000000000 }")
        local fps=$(awk "BEGIN { printf \"%.1f\", $nbframes / $dur_s }")
        local mps=$(awk "BEGIN { printf \"%.2f\", ($nbframes * $wfsize * $wfsize) / $dur_s / 1e6 }")

        local ktime="N/A"
        local kfps="N/A"
        if echo "$out" | grep -q "Rendered .* frames"; then
            ktime=$(echo "$out" | grep "Rendered .* frames" | sed -E 's/.* in ([0-9.]+) s.*/\1/')
            kfps=$(echo "$out" | grep "Rendered .* frames" | sed -E 's/.*\(([0-9.]+) fps.*/\1/')
        fi

        local speedup="1.00x"
        if [[ "$backend" == "SCALAR" ]]; then
            t_scalar="$dur_s"
        else
            speedup=$(awk "BEGIN { printf \"%.2fx\", $t_scalar / $dur_s }")
        fi

        printf "%-20s | %8s s | %10s /s | %8s MP/s | %8s (render: %ss, %s fps)\n" \
               "$backend" "$dur_s" "$fps" "$mps" "$speedup" "$ktime" "$kfps"
    done
    echo ""
}

# compare_modes_benchmark <name> <wfsize> <msize> <nbframes>
compare_modes_benchmark() {
    local title="$1"
    local wfsize="$2"
    local msize="$3"
    local nbframes="$4"

    local tstep="0.005"
    local tspan=$(awk "BEGIN { printf \"%.6f\", $nbframes * $tstep }")

    echo "========================================================================"
    echo " Propagation Mode Comparison: $title (${wfsize}x${wfsize}, $nbframes frames)"
    echo "========================================================================"
    printf "%-26s | %10s | %12s | %12s | %10s\n" \
           "Propagation Engine" "Wall Time" "Frame Rate" "Throughput" "Rel. Speed"
    echo "---------------------------+------------+--------------+--------------+-----------"

    declare -A mode_labels=(
        [0]="Geometric (mode 0)"
        [1]="Split-Step Fresnel (mode 1)"
        [2]="Rytov Option B (mode 2)"
        [3]="Rytov Option C (mode 3)"
    )

    local t_base=1.0
    for m in 1 0 2 3; do
        local amp=$([[ $m -gt 0 ]] && echo 1 || echo 0)
        cat > WFsim.conf << EOF
TURBULENCE_REF_WAVEL       0.500000
TURBULENCE_SEEING          0.60000
TURBULENCE_PROF_FILE       turbul.prof
ZENITH_ANGLE               0.0
SOURCE_XPOS                0.0
SOURCE_YPOS                0.0
WFOUTPUT                   0
SHM_OUTPUT                 0
MAKE_SWAVEFRONT            1
SLAMBDA                    1.650
SWF_WRITE2DISK             0
SHM_SOUTPUT                0
WFsize                     $wfsize
PUPIL_SCALE                0.020000
WFTIME_STEP                $tstep
TIME_SPAN                  $tspan
MASTER_SIZE                $msize
WAVEFRONT_AMPLITUDE        $amp
FRESNEL_PROPAGATION        $m
FRESNEL_PROPAGATION_BIN    1000.0
EOF

        local fps_id="twfsbench_cmp_${wfsize}_${m}"
        local start_ns=$(date +%s%N)
        local out
        out=$("$MKWFS_EXEC" -n "$fps_id" exec 1.65 0 2>&1)
        local end_ns=$(date +%s%N)
        rm -f "$MILK_SHM_DIR/${fps_id}.fps.shm"

        local dur_s=$(awk "BEGIN { printf \"%.3f\", ($end_ns - $start_ns) / 1000000000 }")
        local fps=$(awk "BEGIN { printf \"%.1f\", $nbframes / $dur_s }")
        local mps=$(awk "BEGIN { printf \"%.2f\", ($nbframes * $wfsize * $wfsize) / $dur_s / 1e6 }")

        if [[ $m -eq 1 ]]; then
            t_base="$dur_s"
        fi
        local sp=$(awk "BEGIN { printf \"%.2fx\", $t_base / $dur_s }")

        printf "%-26s | %8s s | %10s /s | %8s MP/s | %10s\n" \
               "${mode_labels[$m]}" "$dur_s" "$fps" "$mps" "$sp"
    done
    echo ""
}

# 1. Fast series: 64x64 pupil, 200 frames
run_series_benchmark "Small Pupil Series" 64 256 200 0

# 2. Production series: 256x256 pupil, 100 frames
run_series_benchmark "Medium Pupil Series" 256 1024 100 0

# 3. ELT scale series: 512x512 pupil, 50 frames
run_series_benchmark "Large Pupil Series" 512 2048 50 0

# 4. Multi-mode comparison on production grid (256x256, 50 frames)
compare_modes_benchmark "Production Grid (256x256)" 256 1024 50

# 5. Multi-mode comparison on ELT scale grid (512x512, 20 frames)
compare_modes_benchmark "ELT Scale Grid (512x512)" 512 2048 20

echo "Wavefront series benchmark completed successfully!"
