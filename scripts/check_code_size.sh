#!/usr/bin/env bash
#
# scripts/check_code_size.sh
# Ratchet-based code size and function scope enforcement for milkatmturb.
#
# Enforces:
#   - File length: <= 1000 lines (soft limit: 600 lines)
#   - Function body length: <= 150 lines (soft limit: 60 lines)
#   - main() body length: <= 80 lines (soft limit: 40 lines)
#
# Exemptions:
#   - Preceded by '/* size-exempt: <reason> */'
#   - Vendored code: src/AtmosphericTurbulence/nrlmsise-00.20131225/
#   - Legacy items tracked in scripts/code_size_baseline.txt
#
# Usage:
#   ./scripts/check_code_size.sh                  # Verify code against baseline
#   ./scripts/check_code_size.sh --verbose        # Show warnings for soft limits
#   ./scripts/check_code_size.sh --update-baseline # Update baseline with current counts
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BASELINE_FILE="$SCRIPT_DIR/code_size_baseline.txt"

FILE_HARD_LIMIT=1000
FILE_SOFT_LIMIT=600
FUNC_HARD_LIMIT=150
FUNC_SOFT_LIMIT=60
MAIN_HARD_LIMIT=80
MAIN_SOFT_LIMIT=40

VERBOSE=0
UPDATE_BASELINE=0

for arg in "$@"; do
    case "$arg" in
        -v|--verbose)
            VERBOSE=1
            ;;
        --update-baseline)
            UPDATE_BASELINE=1
            ;;
        -h|--help)
            echo "Usage: $0 [--verbose] [--update-baseline]"
            exit 0
            ;;
        *)
            echo "Unknown argument: $arg" >&2
            exit 1
            ;;
    esac
done

cd "$ROOT_DIR"

TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

AWK_PARSER="$TMP_DIR/parse_c_code.awk"
cat > "$AWK_PARSER" << 'AWK_EOF'
BEGIN {
    in_fn = 0;
    fn_name = "";
    fn_start = 0;
    brace_depth = 0;
    in_proto = 0;
    proto_accum = "";
    exempt = 0;
}

/size-exempt:/ {
    exempt = 1;
}

/^[ \t]*#/ { next; }

!in_fn && brace_depth == 0 {
    if ($0 ~ /^[a-zA-Z_]/) {
        proto_accum = $0;
        in_proto = 1;
    } else if (in_proto) {
        proto_accum = proto_accum " " $0;
    }

    if (in_proto && $0 ~ /;\s*$/) {
        in_proto = 0;
        proto_accum = "";
        exempt = 0;
    }
}

/^\{/ {
    if (in_proto && proto_accum ~ /\(/ && proto_accum !~ /^(struct|enum|union|typedef)\s+/) {
        in_fn = 1;
        in_proto = 0;
        fn_start = FNR;
        brace_depth = 1;

        sub(/\(.*$/, "", proto_accum);
        match(proto_accum, /[a-zA-Z0-9_]+[ \t]*$/);
        fn_name = substr(proto_accum, RSTART, RLENGTH);
        gsub(/[ \t]/, "", fn_name);
        is_exempt = exempt;
        exempt = 0;
        next;
    } else {
        in_proto = 0;
        proto_accum = "";
        exempt = 0;
    }
}

in_fn {
    t = $0;
    gsub(/"([^"\\]|\\.)*"/, "", t);
    gsub(/'([^'\\]|\\.)*'/, "", t);
    n_open = gsub(/\{/, "{", t);
    n_close = gsub(/\}/, "}", t);
    brace_depth += (n_open - n_close);

    if (brace_depth <= 0 || (FNR > fn_start && /^}/)) {
        fn_len = FNR - fn_start + 1;
        limit = (fn_name == "main") ? main_limit : func_limit;
        if (fn_len > limit && !is_exempt) {
            print "FUNC\t" FILENAME ":" fn_name "\t" fn_len;
        }
        in_fn = 0;
        fn_name = "";
        proto_accum = "";
        brace_depth = 0;
        is_exempt = 0;
    }
}
AWK_EOF

# List tracked C and header files, ignoring third_party / vendored code and symlinks
RAW_FILES=$(git ls-files 'src/*.[ch]' | \
    grep -Ev '(nrlmsise-00.20131225|cjson|third_party)' | sort)
FILES=""
for f in $RAW_FILES; do
    if [ ! -L "$f" ]; then
        FILES="$FILES $f"
    fi
done

if [ "$UPDATE_BASELINE" -eq 1 ]; then
    echo "Updating code size baseline: $BASELINE_FILE"
    {
        echo "# milkatmturb Code Size Ratchet Baseline"
        echo "# Format: 2 lines per entry: <TYPE> <MAX_LINES> followed by <TARGET>"
        echo "# Generated on $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
        echo "# Update via: ./scripts/check_code_size.sh --update-baseline"
        echo ""
        echo "# --- File Lengths (> 1000 lines) ---"
        for f in $FILES; do
            flen=$(wc -l < "$f")
            if [ "$flen" -gt "$FILE_HARD_LIMIT" ]; then
                echo "FILE $flen"
                echo "$f"
            fi
        done

        echo ""
        echo "# --- Function Lengths (> 150 lines, or main > 80 lines) ---"
        TMP_FUNCS="$TMP_DIR/raw_funcs.tsv"
        for f in $FILES; do
            awk -v func_limit="$FUNC_HARD_LIMIT" -v main_limit="$MAIN_HARD_LIMIT" \
                -f "$AWK_PARSER" "$f" || true
        done | sort -k2,2 > "$TMP_FUNCS"

        while IFS=$'\t' read -r type target flen; do
            echo "$type $flen"
            echo "$target"
        done < "$TMP_FUNCS"
    } > "$BASELINE_FILE"
    echo "Baseline updated successfully ($(grep -cE '^(FILE|FUNC) ' "$BASELINE_FILE") items)."
    exit 0
fi

if [ ! -f "$BASELINE_FILE" ]; then
    echo "ERROR: Baseline file not found: $BASELINE_FILE" >&2
    echo "Run '$0 --update-baseline' to create it." >&2
    exit 1
fi

CURRENT_METRICS="$TMP_DIR/current_metrics.tsv"
for f in $FILES; do
    flen=$(wc -l < "$f")
    if [ "$flen" -gt "$FILE_HARD_LIMIT" ]; then
        echo -e "FILE\t$f\t$flen"
    fi
    awk -v func_limit="$FUNC_HARD_LIMIT" -v main_limit="$MAIN_HARD_LIMIT" \
        -f "$AWK_PARSER" "$f" || true
done > "$CURRENT_METRICS"

# Compare current metrics against baseline using awk
VERIFY_SCRIPT="$TMP_DIR/verify_ratchet.awk"
cat > "$VERIFY_SCRIPT" << 'VERIFY_AWK'
BEGIN {
    FS = "\t";
    errors = 0;
    improvements = 0;
}

# First file: baseline
FILENAME == ARGV[1] {
    if ($0 ~ /^[ \t]*#/ || NF == 0) next;
    if ($0 ~ /^(FILE|FUNC) /) {
        split($0, parts, " ");
        type = parts[1];
        limit = parts[2] + 0;
        if ((getline target) > 0) {
            gsub(/^[ \t]+|[ \t]+$/, "", target);
            key = type "\t" target;
            baseline_limit[key] = limit;
            in_baseline[key] = 1;
        }
    }
    next;
}

# Second file: current metrics
FILENAME == ARGV[2] {
    if ($0 ~ /^[ \t]*#/ || NF < 3) next;
    type = $1;
    target = $2;
    curr_len = $3 + 0;
    key = type "\t" target;

    if (!(key in in_baseline)) {
        printf("ERROR [NEW VIOLATION]: %s '%s' exceeds hard limit (%d lines)\n",
               type, target, curr_len);
        errors++;
    } else {
        base_len = baseline_limit[key];
        if (curr_len > base_len) {
            printf("ERROR [REGRESSION]: %s '%s' grew beyond baseline (%d > %d lines)\n",
                   type, target, curr_len, base_len);
            errors++;
        } else if (curr_len < base_len) {
            if (verbose) {
                printf("NOTICE [IMPROVED]: %s '%s' shrunk (%d -> %d lines)\n",
                       type, target, base_len, curr_len);
            }
            improvements++;
        }
    }
}

END {
    if (errors > 0) {
        printf("\nFAILED: %d code size violation(s) or regression(s) found.\n", errors);
        exit 1;
    } else {
        if (improvements > 0) {
            printf("PASS: All limits within baseline. %d item(s) improved!\n", improvements);
            printf("      Run './scripts/check_code_size.sh --update-baseline' ");
            printf("to tighten ratchet.\n");
        } else {
            printf("PASS: All files and functions comply with size constraints and baseline.\n");
        }
        exit 0;
    }
}
VERIFY_AWK

awk -v verbose="$VERBOSE" -f "$VERIFY_SCRIPT" "$BASELINE_FILE" "$CURRENT_METRICS"
