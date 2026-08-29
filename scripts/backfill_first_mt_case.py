#!/usr/bin/env python
"""Backfill the ``first_mt_case`` column into existing PSyGrid HDF5 files.

The column was introduced in the grid-creation step; grids created before then
lack it. This script derives it, for each model, from the existing
``termination_flag_2`` string using the same logic as grid creation, so that
new and legacy grids agree exactly.

The value of ``first_mt_case`` is the Maltsev+25 MT class of the first mass
transfer episode (whichever star is the donor): 'case_A', 'case_B' or
'case_C'. Models without a qualifying MT interaction keep the ``termination
flag_2`` token verbatim (e.g. 'no_RLOF', 'initial_RLOF', 'not_converged'),
matching the stored flag.

Usage
-----
    python backfill_first_mt_case.py grid1.h5 [grid2.h5 ...]

Each grid is modified in place. Grids that already contain the column are
skipped.

"""
import sys

import numpy as np

from posydon.grids.psygrid import PSyGrid
from posydon.utils.common_functions import mt_class_from_cumulative


def backfill_first_mt_case(grid_path):
    """Add the ``first_mt_case`` column to a single grid file."""
    grid = PSyGrid(grid_path, lazy=False)
    if "first_mt_case" in grid.final_values.dtype.names:
        print(f"{grid_path}: column 'first_mt_case' already exists, skipping")
        grid.close()
        return False
    if "termination_flag_2" not in grid.final_values.dtype.names:
        grid.close()
        raise ValueError(f"{grid_path}: no 'termination_flag_2' column found")

    values = np.asarray([
        mt_class_from_cumulative(tf2, retain_flag_if_no_mt=True)
        for tf2 in grid.final_values["termination_flag_2"]
    ])
    grid.add_column("first_mt_case", values)
    grid.close()
    print(f"{grid_path}: added 'first_mt_case' for {len(values)} models")
    return True


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for path in sys.argv[1:]:
        backfill_first_mt_case(path)
