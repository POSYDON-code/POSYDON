#!/bin/bash

# A script for validating the outputs of 100 binaries,
# which can be compared to a baseline to monitor changes to the code.

# script usage: ./validate_binaries.sh --branch candidate_branch

BRANCH=$1
SUFFIX=$2

# run candidate binaries and save to file
./evolve_binaries.sh "$BRANCH"

# compare quantitative, qualitative, warnings/errors, structured output
# create outputs/comparison_branchname.txt
SAFE_BRANCH="${BRANCH//\//_}"
CANDIDATE_FILE="outputs/candidate_${SAFE_BRANCH}.h5"
COMPARISON_FILE="outputs/comparison_${SAFE_BRANCH}${SUFFIX:+_$SUFFIX}.txt"
python compare_runs.py baseline.h5 "$CANDIDATE_FILE" > "$COMPARISON_FILE"

if [ $? -ne 0 ]; then
    echo "Error: compare_runs.py failed. Check $COMPARISON_FILE"
    exit 1
fi

echo "Binary evolution comparison saved to $COMPARISON_FILE"
