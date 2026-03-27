#!/bin/bash
# =============================================================================
# validate_binaries.sh: Run the full validation pipeline:
#   1. Evolve test binaries on a candidate branch
#   2. Compare against baseline files
#
# Usage:
#   ./validate_binaries.sh <candidate_branch> [baseline_branch] [metallicities]
#                          [--loose] [--rtol VALUE] [--atol VALUE]
#
# Positional arguments:
#   candidate_branch   Branch or tag to validate (required)
#   baseline_branch    Branch or tag to compare against (default: main)
#   metallicities      Space-separated list of Z values, quoted
#                      (default: "2 1 0.45 0.2 0.1 0.01 0.001 0.0001")
#
# Tolerance flags (passed through to compare_runs.py):
#   --loose            Use relaxed floating-point tolerances
#                      (rtol=1e-12, atol=1e-15 unless overridden)
#   --rtol VALUE       Set explicit relative tolerance as per np.allclose
#   --atol VALUE       Set explicit absolute tolerance as per np.allclose
#
#   --rtol and --atol can be combined with --loose (explicit values take
#   precedence over the --loose defaults) or used on their own without --loose.
#
# Examples:
#   ./validate_binaries.sh feature/new-SN                      # compare vs main, exact
#   ./validate_binaries.sh feature/new-SN v2.1.0               # compare vs v2.1.0, exact
#   ./validate_binaries.sh feature/new-SN main "1 0.45"        # subset of metallicities
#   ./validate_binaries.sh feature/new-SN --loose              # relaxed tolerances
#   ./validate_binaries.sh feature/new-SN main --rtol 1e-8     # custom rtol, default atol
#   ./validate_binaries.sh feature/new-SN main "1 0.45" --loose --atol 1e-10
#
# Prerequisites:
#   Run generate_baseline.sh first to create baseline files.
#
# Output:
#   outputs/<branch>/comparison_<Z>Zsun.txt — per-metallicity comparison reports
#   outputs/<branch>/comparison_summary.txt — overall summary
# =============================================================================

set -euo pipefail

# ── Parse arguments ───────────────────────────────────────────────────────

CANDIDATE_BRANCH=${1:?Usage: ./validate_binaries.sh <candidate_branch> [baseline_branch] [metallicities] [--loose] [--rtol VALUE] [--atol VALUE]}
BASELINE_BRANCH=${2:-main}
METALLICITIES=${3:-"2 1 0.45 0.2 0.1 0.01 0.001 0.0001"}
shift $(( $# < 3 ? $# : 3 ))

LOOSE=false
RTOL=""
ATOL=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --loose)
            LOOSE=true
            shift
            ;;
        --rtol)
            RTOL="${2:?--rtol requires a value}"
            shift 2
            ;;
        --atol)
            ATOL="${2:?--atol requires a value}"
            shift 2
            ;;
        *)
            echo "ERROR: Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

DEV_TOOLS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="${DEV_TOOLS_DIR}/script_data"
SRC_DIR="${SCRIPT_DIR}/src"
SAFE_CANDIDATE="${CANDIDATE_BRANCH//\//_}"
SAFE_BASELINE="${BASELINE_BRANCH//\//_}"

BASELINE_DIR="$SCRIPT_DIR/baselines/${SAFE_BASELINE}"
BINARY_OUTPUT_DIR="$SCRIPT_DIR/output/binary_star_tests/${SAFE_CANDIDATE}"
SUMMARY_FILE="$BINARY_OUTPUT_DIR/comparison_summary.txt"

# ── Build compare_runs.py flags ───────────────────────────────────────────
COMPARE_FLAGS=""
if [ "$LOOSE" = "true" ]; then
    COMPARE_FLAGS="$COMPARE_FLAGS --loose"
fi
if [ -n "$RTOL" ]; then
    COMPARE_FLAGS="$COMPARE_FLAGS --rtol $RTOL"
fi
if [ -n "$ATOL" ]; then
    COMPARE_FLAGS="$COMPARE_FLAGS --atol $ATOL"
fi

# Build a human-readable tolerance label for the summary
if [ -n "$RTOL" ] || [ -n "$ATOL" ]; then
    TOL_LABEL="rtol=${RTOL:-default}, atol=${ATOL:-default}"
    if [ "$LOOSE" = "true" ]; then
        TOL_LABEL="$TOL_LABEL (--loose)"
    fi
elif [ "$LOOSE" = "true" ]; then
    TOL_LABEL="--loose (rtol=1e-12, atol=1e-15)"
else
    TOL_LABEL="EXACT (rtol=0, atol=0)"
fi

echo "============================================================"
echo "  POSYDON Binary Validation"
echo "  Candidate:  $CANDIDATE_BRANCH"
echo "  Baseline:   $BASELINE_BRANCH"
echo "  Metallicities: $METALLICITIES"
echo "============================================================"

# ── Verify baseline exists ────────────────────────────────────────────────
if [ ! -d "$BASELINE_DIR" ]; then
    echo "ERROR: Baseline directory not found: $BASELINE_DIR" >&2
    echo "Run generate_baseline.sh first:" >&2
    echo "  ./generate_baseline.sh $BASELINE_BRANCH" >&2
    exit 1
fi

# Check that at least one baseline file exists
BASELINE_COUNT=0
for Z in $METALLICITIES; do
    if [ -f "$BASELINE_DIR/baseline_${Z}Zsun.h5" ]; then
        BASELINE_COUNT=$((BASELINE_COUNT + 1))
    fi
done
if [ $BASELINE_COUNT -eq 0 ]; then
    echo "ERROR: No baseline files found in $BASELINE_DIR for requested metallicities." >&2
    exit 1
fi
echo "  Found $BASELINE_COUNT baseline file(s)."

# ── Step 1: Evolve binaries on candidate branch ──────────────────────────
echo ""
echo "Step 1: Evolving binaries on candidate branch '$CANDIDATE_BRANCH'..."
"$DEV_TOOLS_DIR/run_test_suite.sh" "$CANDIDATE_BRANCH" "" "$METALLICITIES"

# ── Step 2: Compare each metallicity ─────────────────────────────────────
echo ""
echo "Step 2: Comparing results..."

TOTAL=0
PASS=0
FAIL=0
SKIP=0

# Initialize summary
mkdir -p "$BINARY_OUTPUT_DIR"
cat > "$SUMMARY_FILE" << EOF
POSYDON Binary Validation — Comparison Summary
================================================
Candidate branch: $CANDIDATE_BRANCH
Baseline branch:  $BASELINE_BRANCH
Tolerances:       $TOL_LABEL
Date:             $(date -u '+%Y-%m-%d %H:%M:%S UTC')
================================================

EOF

for Z in $METALLICITIES; do
    TOTAL=$((TOTAL + 1))

    BASELINE_FILE="$BASELINE_DIR/baseline_${Z}Zsun.h5"
    CANDIDATE_FILE="$BINARY_OUTPUT_DIR/candidate_${Z}Zsun.h5"
    COMPARISON_FILE="$BINARY_OUTPUT_DIR/comparison_${Z}Zsun.txt"

    echo ""
    echo "--- Z = ${Z} Zsun ---"

    if [ ! -f "$BASELINE_FILE" ]; then
        echo "  SKIP: No baseline file for Z=${Z}"
        echo "Z = ${Z} Zsun: SKIPPED (no baseline)" >> "$SUMMARY_FILE"
        SKIP=$((SKIP + 1))
        continue
    fi

    if [ ! -f "$CANDIDATE_FILE" ]; then
        echo "  FAIL: No candidate file for Z=${Z}"
        echo "Z = ${Z} Zsun: FAIL (no candidate output)" >> "$SUMMARY_FILE"
        FAIL=$((FAIL + 1))
        continue
    fi

    # $COMPARE_FLAGS is intentionally unquoted so it word-splits into
    # separate arguments for compare_runs.py.
    if python "$SRC_DIR/compare_runs.py" "$BASELINE_FILE" "$CANDIDATE_FILE" \
        $COMPARE_FLAGS \
        2>&1 | tee "$COMPARISON_FILE"; then
        echo "  PASS: No differences"
        echo "Z = ${Z} Zsun: PASS" >> "$SUMMARY_FILE"
        PASS=$((PASS + 1))
    else
        echo "  DIFFERENCES DETECTED — see $COMPARISON_FILE"
        echo "Z = ${Z} Zsun: DIFFERENCES DETECTED (see comparison_${Z}Zsun.txt)" >> "$SUMMARY_FILE"
        FAIL=$((FAIL + 1))
    fi
done

# ── Final Summary ─────────────────────────────────────────────────────────
cat >> "$SUMMARY_FILE" << EOF

================================================
TOTAL: $TOTAL | PASS: $PASS | FAIL: $FAIL | SKIP: $SKIP
EOF

echo ""
echo "============================================================"
echo "  Validation Summary"
echo "  TOTAL: $TOTAL | PASS: $PASS | FAIL: $FAIL | SKIP: $SKIP"
echo "  Full summary: $SUMMARY_FILE"
echo "============================================================"

if [ $FAIL -gt 0 ]; then
    exit 1
fi
exit 0
