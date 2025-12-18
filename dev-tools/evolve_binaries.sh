#!/bin/bash

# Script usage: ./evolve_binaries.sh <branch>
# This script clones the POSYDON repo to the specified branch (defaults to 'main'),
# copies evolve_binaries.py, runs it, and saves output to evolve_binaries.out

# Set default branch to 'main' if not provided
BRANCH=${1:-main}
REPO_URL="https://github.com/POSYDON-code/POSYDON"

if [[ -n "$2" ]]; then
    SHA=$2
    WORK_DIR="POSYDON_${BRANCH}_${SHA}"
else
    WORK_DIR="POSYDON_$BRANCH"
fi

# Remove existing directory if it exists
if [ -d "$WORK_DIR" ]; then
    echo "🗑️  Removing existing directory: $WORK_DIR"
    rm -rf "$WORK_DIR"
fi

echo "📁 Creating working directory: $WORK_DIR"
# Create the working directory
mkdir -p "$WORK_DIR"

FULL_PATH="$(realpath "$WORK_DIR")"
CLONE_DIR="$FULL_PATH/POSYDON"

OUTPUT_DIR="$FULL_PATH/outputs"
LOG_DIR="$FULL_PATH/logs"

SAFE_BRANCH="${BRANCH//\//_}"
OUTPUT_FILE="$OUTPUT_DIR/candidate_${SAFE_BRANCH}.h5"
LOG_FILE="$LOG_DIR/evolve_${SAFE_BRANCH}.log"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

echo "📋 Copying script_data folder"
# copy the script_data folder
cp -r "./script_data" "$WORK_DIR"

cd "$WORK_DIR"

# Initialize conda for bash
echo "🔧 Initializing conda"
# Source conda's shell integration
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
elif [ -f "/opt/homebrew/Caskroom/miniconda/base/etc/profile.d/conda.sh" ]; then
    source "/opt/homebrew/Caskroom/miniconda/base/etc/profile.d/conda.sh"
else
    echo -e "\033[31mError: Could not find conda installation. Please check your conda setup.\033[0m"
    exit 1
fi

# Clone the repository to the specified branch
echo "🔄 Cloning POSYDON repository (branch: $BRANCH)"
if ! git clone -b "$BRANCH" "$REPO_URL" "$CLONE_DIR" 2>&1 | sed 's/^/  /'; then
    echo -e "\033[31mError: Failed to clone branch '$BRANCH'. Please check if the branch exists.\033[0m"
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

# Create conda environment for POSYDON v2
echo "🐍 Creating conda environment"
conda create --prefix="$FULL_PATH/conda_env" python=3.11 -y -q 2>&1 | sed 's/^/  /'

echo "⚡ Activating conda environment"
conda activate "$FULL_PATH/conda_env"

# install POSYDON manually
echo "📦 Installing POSYDON"
pip install -e "$CLONE_DIR" -q 2>&1 | sed 's/^/  /'

echo "🚀 Running evolve_binaries.py"
# # Run the Python script and capture output (stdout and stderr)
python script_data/1Zsun_binaries_suite.py --output "$OUTPUT_FILE" 2>&1 | tee "$LOG_FILE" 

if [ ! -f "$OUTPUT_FILE" ]; then
    echo "ERROR: Results file was not created: $OUTPUT_FILE"
    exit 2
fi

if [ $? -ne 0 ]; then
    echo "ERROR: Python script exited with an error. Check $LOG_FILE for details."
    exit 3
fi

echo -e "✅ Script completed. Output saved to \n$OUTPUT_FILE"
