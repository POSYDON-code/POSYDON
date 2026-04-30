#!/bin/bash
# =============================================================================
# generate_baseline.sh — Generate baseline HDF5 files from a designated branch.
#
# This runs the binary validation suite against a chosen branch (or commit)
# and saves the results as the baseline for future comparisons.
#
# Usage:
#   ./generate_baseline.sh <branch> [sha] [metallicities]
#   ./generate_baseline.sh --promote <branch> [metallicities]
#
# Examples:
#   ./generate_baseline.sh main                          # evolve + save baseline, all Z
#   ./generate_baseline.sh v2.1.0                        # baseline from a release tag
#   ./generate_baseline.sh main abc123f                  # baseline from a specific commit
#   ./generate_baseline.sh main "" "1 0.45"              # baseline for subset of Z
#   ./generate_baseline.sh --promote main                # promote existing outputs to baseline
#   ./generate_baseline.sh --promote main "1 0.45"       # promote subset of existing outputs
#
# Output:
#   baselines/<branch>/baseline_<Z>Zsun.h5  — one file per metallicity
#   baselines/<branch>/baseline_info.txt    — records branch, commit SHA, date
# =============================================================================

set -euo pipefail

# ── Parse arguments ───────────────────────────────────────────────────────
PROMOTE=false
if [ "${1:-}" = "--promote" ]; then
    PROMOTE=true
    shift
fi

BRANCH=${1:-main}
if [ "$PROMOTE" = true ]; then
    SHA=""
    METALLICITIES=${2:-"2 1 0.45 0.2 0.1 0.01 0.001 0.0001"}
else
    SHA=${2:-}
    METALLICITIES=${3:-"2 1 0.45 0.2 0.1 0.01 0.001 0.0001"}
fi

DEV_TOOLS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="${DEV_TOOLS_DIR}/script_data"
SAFE_BRANCH="${BRANCH//\//_}"
BASELINE_DIR="$SCRIPT_DIR/baselines/${SAFE_BRANCH}"
BINARY_CANDIDATE_DIR="$SCRIPT_DIR/output/binary_star_tests/${SAFE_BRANCH}"

echo "============================================================"
echo "  POSYDON Binary Validation — Generating Baseline"
echo "  Branch:        $BRANCH"
if [ "$PROMOTE" = true ]; then
    echo "  Mode:          --promote (using existing outputs)"
else
    echo "  SHA:           ${SHA:-HEAD}"
fi
echo "  Metallicities: $METALLICITIES"
echo "  Output dir:    $BASELINE_DIR"
echo "============================================================"

# ── Step 1: Evolve binaries (skip if --promote) ──────────────────────────
if [ "$PROMOTE" = true ]; then
    echo ""
    echo "Step 1: SKIPPED (--promote: using existing outputs in $BINARY_CANDIDATE_DIR)"

    if [ ! -d "$BINARY_CANDIDATE_DIR" ]; then
        echo "ERROR: No outputs found at $BINARY_CANDIDATE_DIR" >&2
        echo "Run evolve_binaries.sh first, or drop --promote to evolve from scratch." >&2
        exit 1
    fi
else
    echo ""
    echo "Step 1: Evolving binaries on branch '$BRANCH'..."
    "${DEV_TOOLS_DIR}/run_test_suite.sh" "$BRANCH" "$SHA" "$METALLICITIES"
fi

# ── Step 2: Copy results into the baselines directory ────────────────────
echo ""
echo "Step 2: Copying results to baseline directory..."

mkdir -p "$BASELINE_DIR"

COPIED=0

for Z in $METALLICITIES; do
    SRC="$BINARY_CANDIDATE_DIR/candidate_${Z}Zsun.h5"
    DST="$BASELINE_DIR/baseline_${Z}Zsun.h5"

    if [ -f "$SRC" ]; then
        cp "$SRC" "$DST"
        echo "  Saved: $DST"
        COPIED=$((COPIED + 1))
    else
        echo "  WARNING: Missing output for Z=${Z}: $SRC" >&2
    fi
done

# ── Step 3: Record baseline metadata ─────────────────────────────────────
CLONE_DIR="$SCRIPT_DIR/workdirs/POSYDON_${SAFE_BRANCH}/POSYDON"
ACTUAL_SHA=""
if [ -d "$CLONE_DIR" ]; then
    ACTUAL_SHA=$(cd "$CLONE_DIR" && git rev-parse HEAD 2>/dev/null || echo "unknown")
fi

# Extract PATH_TO_POSYDON_DATA from the first available baseline HDF5 file
POSYDON_DATA_PATH="unknown"
for Z in $METALLICITIES; do
    H5="$BASELINE_DIR/baseline_${Z}Zsun.h5"
    if [ -f "$H5" ]; then
        POSYDON_DATA_PATH=$(python -c "
import pandas as pd
with pd.HDFStore('$H5', mode='r') as s:
    print(s['/metadata']['path_to_posydon_data'].iloc[0])
" 2>&1) || {
            echo "  WARNING: Could not read POSYDON data path from $H5"
            echo "  ($POSYDON_DATA_PATH)"
            POSYDON_DATA_PATH="unknown"
        }
        break
    fi
done

# Check completeness from HDF5 metadata
INCOMPLETE=""
for Z in $METALLICITIES; do
    DST="$BASELINE_DIR/baseline_${Z}Zsun.h5"
    if [ -f "$DST" ]; then
        MISSING=$(python3 -c "
import pandas as pd, sys
try:
    with pd.HDFStore('$DST', mode='r') as s:
        m = s['/metadata']
        n = int(m['n_missing'].iloc[0])
        if n > 0:
            print(f'Z={Z}: {n} missing — {m[\"missing_ids\"].iloc[0]}')
except Exception as e:
    print(f'Z={Z}: could not read metadata ({e})', file=sys.stderr)
" 2>/dev/null)
        if [ -n "$MISSING" ]; then
            echo "  ⚠️  $MISSING"
            INCOMPLETE="${INCOMPLETE}  ${MISSING}\n"
        fi
    fi
done

INFO_FILE="$BASELINE_DIR/baseline_info.txt"
cat > "$INFO_FILE" << EOF
POSYDON Binary Validation Baseline
===================================
Branch:        $BRANCH
Commit SHA:    ${ACTUAL_SHA:-unknown}
Requested SHA: ${SHA:-HEAD}
Mode:          $([ "$PROMOTE" = true ] && echo "promoted from existing outputs" || echo "evolved from scratch")
Generated:     $(date -u '+%Y-%m-%d %H:%M:%S UTC')
Metallicities: $METALLICITIES
Files:         $COPIED
POSYDON data:  $POSYDON_DATA_PATH
EOF

if [ -n "$INCOMPLETE" ]; then
    printf "\nINCOMPLETE BASELINES:\n%b\n" "$INCOMPLETE" >> "$INFO_FILE"
    echo ""
    echo "⚠️  WARNING: Some baselines have missing binaries. See $INFO_FILE"
fi

echo ""
echo "============================================================"
echo "  Baseline generated: $COPIED file(s)"
echo "  Info: $INFO_FILE"
echo "  Directory: $BASELINE_DIR"
echo "============================================================"

if [ $COPIED -eq 0 ]; then
    echo "ERROR: No baseline files were created!" >&2
    exit 1
fi
