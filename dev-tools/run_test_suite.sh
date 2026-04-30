#!/bin/bash

# =============================================================================
# evolve_binaries.sh — Clone a POSYDON branch, install it, and run the
# binary validation suite at all requested metallicities.
#
# Usage:
#   ./evolve_binaries.sh <branch> [sha] [metallicities]
#
# Examples:
#   ./evolve_binaries.sh main                          # all metallicities
#   ./evolve_binaries.sh feature/my-fix abc123f        # specific commit
#   ./evolve_binaries.sh main "" "1 0.45 0.1"          # subset of metallicities
#
# Output structure:
#   outputs/<branch>/candidate_<Z>Zsun.h5   — evolution results per metallicity
#   logs/<branch>/evolve_<Z>Zsun.log        — log per metallicity
#   workdirs/POSYDON_<branch>/              — cloned repo + conda env
#
# The resolved commit SHA and branch name are recorded in each HDF5 file's
# /metadata table for provenance tracking in comparison reports.
# =============================================================================

set -euo pipefail
# Load git if needed
if ! command -v git >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        module load git
    fi
fi

# ── Configuration ──────────────────────────────────────────────────────────
ALL_METALLICITIES="2 1 0.45 0.2 0.1 0.01 0.001 0.0001"

BRANCH=${1:-main}
SHA=${2:-}
METALLICITIES=${3:-$ALL_METALLICITIES}

REPO_URL="https://github.com/POSYDON-code/POSYDON"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/script_data"

# Sanitize branch name for filesystem
SAFE_BRANCH="${BRANCH//\//_}"


# Directories (all relative to SCRIPT_DIR, the dev-tools root)
WORK_DIR="$SCRIPT_DIR/workdirs/POSYDON_${SAFE_BRANCH}"
BINARY_OUTPUT_DIR="$SCRIPT_DIR/output/binary_star_tests/${SAFE_BRANCH}"
POP_OUTPUT_DIR="$SCRIPT_DIR/output/population_tests/${SAFE_BRANCH}"
LOG_DIR="$SCRIPT_DIR/logs/${SAFE_BRANCH}"
CLONE_DIR="$WORK_DIR/POSYDON"

mkdir -p ${LOG_DIR} ${BINARY_OUTPUT_DIR}

# ── Conda Setup ────────────────────────────────────────────────────────────
echo "🔧 Initializing conda"
CONDA_SH=""
for candidate in \
    "$HOME/miniconda3/etc/profile.d/conda.sh" \
    "$HOME/anaconda3/etc/profile.d/conda.sh" \
    "/opt/homebrew/Caskroom/miniconda/base/etc/profile.d/conda.sh"; do
    if [ -f "$candidate" ]; then
        CONDA_SH="$candidate"
        break
    fi
done

if [ -z "$CONDA_SH" ]; then
    if command -v conda >/dev/null 2>&1; then
        CONDA_SH="$(conda info --base)/etc/profile.d/conda.sh"
    else
        echo "ERROR: Could not find conda installation." >&2
        exit 1
    fi
fi
source "$CONDA_SH"

# ── Clone Repository ──────────────────────────────────────────────────────
if [ -d "$WORK_DIR" ]; then
    echo "🗑️ Removing existing work directory: $WORK_DIR"
    rm -rf "$WORK_DIR"
fi
mkdir -p "$WORK_DIR"

echo "🔄 Cloning POSYDON repository (branch: $BRANCH)"
if ! git clone -b "$BRANCH" "$REPO_URL" "$CLONE_DIR" 2>&1 | sed 's/^/    /'; then
    echo "ERROR: Failed to clone branch '$BRANCH'." >&2
    exit 1
fi

# if SHA is provided, checkout that commit
if [[ -n "$SHA" ]]; then
    echo "🔄 Checking out commit: $SHA"
    cd "$CLONE_DIR"
    if ! git checkout "$SHA" 2>&1 | sed 's/^/  /'; then
        echo -e "\033[31mError: Failed to checkout commit '$SHA'. Please check if the commit exists.\033[0m"
        exit 1
    fi
    cd -
fi

# ── Create Conda Environment ─────────────────────────────────────────────
ENV_PREFIX="$WORK_DIR/conda_env"

echo "🐍 Creating conda environment at $ENV_PREFIX"
conda create --prefix="$ENV_PREFIX" python=3.11 -y -q 2>&1 | sed 's/^/    /'
conda activate "$ENV_PREFIX"

echo "📦 Installing POSYDON"
pip install -e "$CLONE_DIR" -q 2>&1 | sed 's/^/    /'

# ── Run Suite for Each Metallicity ────────────────────────────────────────
SUITE_SCRIPT="$SCRIPT_DIR/src/binaries_suite.py"
FAILED=0

# Resolve the actual commit SHA for metadata
ACTUAL_SHA=$(cd "$CLONE_DIR" && git rev-parse HEAD 2>/dev/null || echo "unknown")
echo "  Resolved SHA: $ACTUAL_SHA"

# override environment's PATH_TO_POSYDON variable to point to the
# current branch's clone for these tests
#PATH_TO_POSYDON=$CLONE_DIR

# copy this branch's default .ini file to perform tests
DEFAULT_INI="${CLONE_DIR}/posydon/popsyn/population_params_default.ini"
TEST_INI="${SCRIPT_DIR}/inlists/${SAFE_BRANCH}_test_params.ini"
cp $DEFAULT_INI $TEST_INI
sed -i 's/^\([[:space:]]*\)entropy *= *.*/\1entropy = 0/' $TEST_INI

for Z in $METALLICITIES; do
    OUTPUT_FILE="$BINARY_OUTPUT_DIR/candidate_${Z}Zsun.h5"
    LOG_FILE="$LOG_DIR/evolve_${Z}Zsun.log"

    echo ""
    echo "============================================================"
    echo "  🚀 Evolving binaries for Z = ${Z} Zsun"
    echo "  Output: $OUTPUT_FILE"
    echo "  Log:    $LOG_FILE"
    echo "============================================================"

    #python "$SUITE_SCRIPT" \
    #    --metallicity "$Z" \
    #    --output "$OUTPUT_FILE" \
    #    --branch "$BRANCH" \
    #    --sha "$ACTUAL_SHA" \
    #    2>&1 | tee "$LOG_FILE"
    #EXIT_CODE=${PIPESTATUS[0]}

    #if [ $EXIT_CODE -eq 137 ]; then
    #    echo "ERROR: Process killed (likely OOM) for Z=${Z}. Exit code 137 (SIGKILL)." >&2
    #    echo "  Consider increasing job memory." >&2
    #    FAILED=$((FAILED + 1))
    #elif [ $EXIT_CODE -ne 0 ]; then
    #    echo "WARNING: Suite failed for Z=${Z} (exit code $EXIT_CODE). Check $LOG_FILE" >&2
    #    FAILED=$((FAILED + 1))
    #elif [ ! -f "$OUTPUT_FILE" ]; then
    #    echo "WARNING: Output file not created for Z=${Z}" >&2
    #    FAILED=$((FAILED + 1))
    #else
    #    echo "  Z=${Z} Zsun complete."
    #fi

done

# Run population synthesis tests and capture output (stdout and stderr)
SUITE_SCRIPT="$SCRIPT_DIR/src/popsynth_suite.py"
FAILED=0
LOG_FILE="$LOG_DIR/evolve_populations.log"

MULTIZ_INI="${SCRIPT_DIR}/inlists/${SAFE_BRANCH}_test_multiZ_params.ini"
cp $TEST_INI $MULTIZ_INI
sed -i 's/^\([[:space:]]*\)metallicities *= *\[.*\]/\1metallicities = [2., 1., 1e-4]/' $MULTIZ_INI

echo ""
echo "============================================================"
echo "  🚀 Evolving populations"
echo "  Output: $POP_OUTPUT_DIR"
echo "  Log:    $LOG_FILE"
echo "============================================================"

python "$SUITE_SCRIPT" \
    --output "$POP_OUTPUT_DIR" \
    --ini "$TEST_INI" \
    --multiz "$MULTIZ_INI" \
    2>&1 | tee "$LOG_FILE"
EXIT_CODE=${PIPESTATUS[0]}

if [ $EXIT_CODE -eq 137 ]; then
    echo "ERROR: Process killed (likely OOM) for popsynth_suite.py. Exit code 137 (SIGKILL)." >&2
    echo "  Consider increasing job memory." >&2
    FAILED=$((FAILED + 1))
elif [ $EXIT_CODE -ne 0 ]; then
    echo "WARNING: popsynth_suite.py failed (exit code $EXIT_CODE). Check $LOG_FILE" >&2
    FAILED=$((FAILED + 1))
elif [ ! -f "$OUTPUT_FILE" ]; then
    echo "WARNING: Output file not created for popsynth_suite.py" >&2
    FAILED=$((FAILED + 1))
else
    echo "  popsynth_suite.py completed."
fi

# ── Deactivate Environment ────────────────────────────────────────────────
conda deactivate

echo ""
echo "============================================================"
if [ $FAILED -eq 0 ]; then
    echo "✅ All metallicities completed successfully."
else
    echo "Completed with $FAILED failure(s)."
fi
echo "  Outputs in: $BINARY_OUTPUT_DIR/"
echo "============================================================"
echo ""
echo "============================================================"
if [ $FAILED -eq 0 ]; then
    echo "✅ All population synthesis tests completed successfully."
else
    echo "Completed with $FAILED failure(s)."
fi
echo "  Outputs in: $POP_OUTPUT_DIR/"
echo "============================================================"


exit $FAILED
