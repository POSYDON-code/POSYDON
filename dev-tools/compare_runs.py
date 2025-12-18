#!/usr/bin/env python3

"""
Compare evolution outcomes of set of binaries saved to file with script_data/1Zsun_binaries_suite.py
Used for validation of the branch.

Author: Elizabeth Teng
"""


import sys
import h5py
import numpy as np

def compare_datasets(base, cand, path="/"):
    """
    Recursively compare HDF5 datasets.
    Returns a list of strings with differences.
    """
    differences = []
    
    for key in base.keys():
        base_item = base[key]
        if key not in cand:
            differences.append(f"{path}{key} missing in candidate")
            continue
        
        cand_item = cand[key]

        if isinstance(base_item, h5py.Dataset):
            base_data = base_item[()]
            cand_data = cand_item[()]
            
            if np.issubdtype(base_data.dtype, np.number):
                if not np.allclose(base_data, cand_data, rtol=1e-5, atol=1e-8):
                    differences.append(f"{path}{key} numeric mismatch")
            else:
                if not np.array_equal(base_data, cand_data):
                    differences.append(f"{path}{key} non-numeric mismatch")
        
        elif isinstance(base_item, h5py.Group):
            differences.extend(compare_datasets(base_item, cand_item, path=f"{path}{key}/"))
    
    # Check for extra keys in candidate
    for key in cand.keys():
        if key not in base:
            differences.append(f"{path}{key} extra in candidate")
    
    return differences


def main():
    if len(sys.argv) != 3:
        print("Usage: compare_runs.py baseline.h5 candidate.h5")
        sys.exit(1)
    
    baseline_file = sys.argv[1]
    candidate_file = sys.argv[2]

    with h5py.File(baseline_file, 'r') as base, h5py.File(candidate_file, 'r') as cand:
        diffs = compare_datasets(base, cand)

    if diffs:
        print("Differences found:")
        for diff in diffs:
            print(" -", diff)
        sys.exit(1)
    else:
        print("No differences detected between baseline and candidate.")


if __name__ == "__main__":
    main()
