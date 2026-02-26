#!/bin/bash
# =============================================================================
# validate_binaries.sh — Run the full validation pipeline:
#   1. Evolve test binaries on a candidate branch
#   2. Compare against baseline files
#
# Usage:
#   ./validate_binaries.sh <candidate_branch> [baseline_branch] [metallicities]
#
# Examples:
#   ./validate_binaries.sh feature/new-SN                    # compare vs main baseline
#   ./validate_binaries.sh feature/new-SN v2.1.0             # compare vs v2.1.0 baseline
#   ./validate_binaries.sh feature/new-SN main "1 0.45"      # subset of metallicities
#
# Prerequisites:
#   Run generate_baseline.sh first to create baseline files.
#
# Output:
#   outputs/<branch>/comparison_<Z>Zsun.txt — per-metallicity comparison reports
#   outputs/<branch>/comparison_summary.txt — overall summary
# =============================================================================

set -euo pipefail

CANDIDATE_BRANCH=${1:?Usage: ./validate_binaries.sh <candidate_branch> [baseline_branch] [metallicities]}
BASELINE_BRANCH=${2:-main}
METALLICITIES=${3:-"2 1 0.45 0.2 0.1 0.01 0.001 0.0001"}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAFE_CANDIDATE="${CANDIDATE_BRANCH//\//_}"
SAFE_BASELINE="${BASELINE_BRANCH//\//_}"

BASELINE_DIR="$SCRIPT_DIR/baselines/${SAFE_BASELINE}"
OUTPUT_DIR="$SCRIPT_DIR/outputs/${SAFE_CANDIDATE}"
SUMMARY_FILE="$OUTPUT_DIR/comparison_summary.txt"

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
"$SCRIPT_DIR/evolve_binaries.sh" "$CANDIDATE_BRANCH" "" "$METALLICITIES"

# ── Step 2: Compare each metallicity ─────────────────────────────────────
echo ""
echo "Step 2: Comparing results..."

TOTAL=0
PASS=0
FAIL=0
SKIP=0

# Initialize summary
mkdir -p "$OUTPUT_DIR"
cat > "$SUMMARY_FILE" << EOF
POSYDON Binary Validation — Comparison Summary
================================================
Candidate branch: $CANDIDATE_BRANCH
Baseline branch:  $BASELINE_BRANCH
Date:             $(date -u '+%Y-%m-%d %H:%M:%S UTC')
================================================

EOF

for Z in $METALLICITIES; do
    TOTAL=$((TOTAL + 1))

    BASELINE_FILE="$BASELINE_DIR/baseline_${Z}Zsun.h5"
    CANDIDATE_FILE="$OUTPUT_DIR/candidate_${Z}Zsun.h5"
    COMPARISON_FILE="$OUTPUT_DIR/comparison_${Z}Zsun.txt"

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

    if python "$SCRIPT_DIR/compare_runs.py" "$BASELINE_FILE" "$CANDIDATE_FILE" \
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