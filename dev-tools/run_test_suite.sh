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

# copy branch's default .ini file for testing and make one into a multi-metallicity .ini file for pop synth tests
cp $CLONE_DIR/posydon/popsyn/population_params_default.ini $FULL_PATH/script_data/inlists/default_test_params.ini
cp $CLONE_DIR/posydon/popsyn/population_params_default.ini $FULL_PATH/script_data/inlists/multiZ_test_params.ini
sed -i 's/metallicities *= *\[1\.\]/metallicities = [1., 0.1]/' $FULL_PATH/script_data/inlists/multiZ_test_params.ini

# override environment's PATH_TO_POSYDON variable to point to the current branch's clone
# for these tests
PATH_TO_POSYDON=$FULL_PATH

# Run binary evolution tests and capture output (stdout and stderr)
OUT_DIR=$FULL_PATH/script_data/output/binary_star_tests
echo "🚀 Running binaries_suite.py"
python script_data/src/binaries_suite.py > $OUT_DIR/evolve_binaries_$BRANCH.out 2>&1
echo -e "✅ Script completed. Output saved to \n$OUT_DIR/evolve_binaries_$BRANCH.out"

# Run population synthesis tests and capture output (stdout and stderr)
OUT_DIR=$FULL_PATH/script_data/output/population_tests
echo "🚀 Running popsynth_suite.py"
python script_data/src/popsynth_suite.py > $OUT_DIR/evolve_pop_$BRANCH.out 2>&1
echo -e "✅ Script completed. Output saved to \n$OUT_DIR/evolve_binaries_$BRANCH.out"
