"""Helper function(s) for tests requiring MESA data

used in:
    - posydon/unit_tests/girds/test_psygrid.py
    - posydon/unit_tests/utils/test_compress_mesa_files.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

import os
from posydon.grids.io import SINGLE_OUTPUT_FILE, BINARY_OUTPUT_FILE
from posydon.grids.psygrid import (DEFAULT_BINARY_HISTORY_COLS,
                                   DEFAULT_STAR_HISTORY_COLS,
                                   DEFAULT_SINGLE_HISTORY_COLS,
                                   DEFAULT_PROFILE_COLS)

# helper functions
def add_MESA_run_files(path, idx, binary_run=True, with_histories=True,\
                       with_profiles=True):
    """Create files of one MESA run

    Parameters
    ----------
    path : str
        The path to the directory containing all the runs.
    idx : int
        Run index. (1, 3, 4, 5, 8, 11, 14, and 17 have a special meaning)
    binary_run : bool (default: True)
        If `True` files of a binary run are created, otherwise run of a single
        star.
    with_histories : bool (default: True)
        If `True` creates history files.
    with_profiles : bool (default: True)
        If `True` creates final profile files.

    Returns
    -------
    str
        Path to the added run.

    """
    if not os.path.exists(path):
        raise NotADirectoryError(f"Wrong path: {path}")
    if not isinstance(idx, int):
        raise TypeError("idx must be an integer")
    elif idx < 0:
        raise ValueError("idx cannot be negative")
    if not isinstance(binary_run, bool):
        raise TypeError("binary_run must be a boolean")
    if not isinstance(with_histories, bool):
        raise TypeError("with_histories must be a boolean")
    # create files expected from a MESA run
    idx_str = "{:.4f}".format(1.0+0.1*idx)
    # get run directory
    if binary_run:
        run_name = f"Zbase_0.0142_m1_{idx_str}_m2_0.9000_initial_z_"\
                   +f"1.4200e-02_initial_period_in_days_{idx_str}e-01_"\
                   +f"grid_index_{idx}"
        if path[-2:] == 'v1': # v1 name
            run_name = f"m1_{idx_str}_m2_0.9000_initial_period_in_days_"\
                   +f"{idx_str}e-01_grid_index_{idx}"
    else:
        run_name = f"Zbase_0.0142_initial_mass_{idx_str}_initial_z_"\
                   +f"1.4200e-02_grid_index_{idx}"
        if path[-2:] == 'v1': # v1 name
            run_name = f"initial_mass_{idx_str}_grid_index_{idx}"
    run_path = os.path.join(path, run_name)
    os.mkdir(run_path)
    # get MESA output file
    if binary_run:
        OUTPUT_FILE = BINARY_OUTPUT_FILE
    else:
        OUTPUT_FILE = SINGLE_OUTPUT_FILE
    with open(os.path.join(run_path, OUTPUT_FILE), "w") as out_file:
        out_file.write(f"Test{idx}")
    # get inlist_grid_points file: containing masses and period
    with open(os.path.join(run_path, "inlist_grid_points"), "w") as\
         inlist_file:
        if binary_run:
            inlist_file.write("&binary_controls\n")
            if idx != 20:
                inlist_file.write(f"\tm1 = 1.{idx}d0\n")
                inlist_file.write("\tm2 = 0.9d0\n")
                inlist_file.write(f"\tinitial_period_in_days = 1.{idx}d-1\n")
            if idx == 1:
                # add value to get Z set in case of getting ignored
                inlist_file.write(f"\tZ = 0.0142d0\n")
                # add fake value
                inlist_file.write(f"\ttest = 0.0\n")
            inlist_file.write("/ ! end of binary_controls namelist")
        else:
            inlist_file.write("&controls\n")
            inlist_file.write(f"\tinitial_mass = 1.{idx}d0\n")
            inlist_file.write("\tinitial_z = 0.0142d0\n")
            inlist_file.write(f"\tZbase = 0.0142d0\n")
            inlist_file.write("/ ! end of star_controls namelist")
    if binary_run and with_histories and (idx != 3):
        # get binary_history.data
        with open(os.path.join(run_path, "binary_history.data"), "w")\
             as binary_history_file:
            for j in range(5):
                # initial values are read from header, rest does not matter
                if j == 1:
                    names = "initial_don_mass initial_acc_mass "\
                            +"initial_period_days\n"
                    binary_history_file.write(names)
                elif j == 2:
                    vals = f"             1.{idx}              0.9 "\
                           +f"            1.{idx}E-01\n"
                    binary_history_file.write(vals)
                else:
                    binary_history_file.write(f"Test HEADER{j}\n")
            # table with column headers and two rows of all the required
            # columns
            row1 = "\n"
            row2 = "\n"
            for j,c in enumerate(DEFAULT_BINARY_HISTORY_COLS):
                if (((idx == 4) or (idx == 11)) and
                    (c in ["model_number", "age"])):
                    # skip special columns
                    continue
                l = len(c)
                fmt = " {:" + f"{l}" + "}"
                binary_history_file.write(fmt.format(c))
                fmt = " {:>" + f"{l}" + "}"
                if ((idx == 1) and (c == "star_1_mass")):
                    # get a high mass star for stop_before_carbon_depletion
                    row1 += fmt.format(j+100)
                else:
                    row1 += fmt.format(j)
                row2 += fmt.format(j+1)
            if idx == 11:
                # add special columns at end for header and first row
                for j,c in enumerate(["model_number", "age"]):
                    l = len(c)
                    fmt = " {:" + f"{l}" + "}"
                    binary_history_file.write(fmt.format(c))
                    fmt = " {:>" + f"{l}" + "}"
                    row1 += fmt.format(j)
            if idx != 14:
                # add data rows
                binary_history_file.write(row1)
                binary_history_file.write(row2)
            if ((idx == 8) and (len(DEFAULT_BINARY_HISTORY_COLS) > 1)):
                # add frist two columns of third row
                l = len(DEFAULT_BINARY_HISTORY_COLS[0])
                fmt = " {:>" + f"{l}" + "}"
                binary_history_file.write("\n"+fmt.format(2))
                l = len(DEFAULT_BINARY_HISTORY_COLS[1])
                fmt = " {:>" + f"{l}" + "}"
                binary_history_file.write(fmt.format(3))
    # get LOGS directories
    if binary_run: # (none, one, two, none, one, two, ...)
        limit_log_dir = 3
    else: # (none, one, none, one, ...)
        limit_log_dir = 2
    for k in range(1, idx%limit_log_dir+1):
        if ((idx==17) and (k==1)):
            continue
        logs_path = os.path.join(run_path, f"LOGS{k}")
        os.mkdir(logs_path)
        if with_histories:
            # get history.data
            with open(os.path.join(logs_path, "history.data"), "w") as\
                 history_file:
                for j in range(5):
                    # initial values are read from header, rest does not matter
                    if j == 1:
                        names = "initial_Z initial_Y initial_m\n"
                        history_file.write(names)
                    elif j == 2:
                        vals = f"   0.0142      0.25       1.{idx}\n"
                        history_file.write(vals)
                    else:
                        history_file.write(f"Test HEADER{j}\n")
                # table with column headers and two rows of all the required
                # columns (model_number and star_age need to align with the
                # binary_history.data where they are the first two columns)
                row1 = "\n"
                row2 = "\n"
                if binary_run:
                    cols = ["model_number", "star_age"]\
                           + DEFAULT_STAR_HISTORY_COLS
                else:
                    cols = DEFAULT_SINGLE_HISTORY_COLS
                for j,c in enumerate(cols):
                    if (((idx == 5) or (idx == 11)) and
                        (c in ["model_number", "star_age"])):
                        # skip extra columns
                        continue
                    l = len(c)
                    fmt = " {:" + f"{l}" + "}"
                    history_file.write(fmt.format(c))
                    fmt = " {:>" + f"{l}" + "}"
                    row1 += fmt.format(j)
                    if "center_h" in c:
                        # creates he depletion
                        row2 += fmt.format(0)
                    else:
                        row2 += fmt.format(j+1)
                if (binary_run and (idx == 11)):
                    # add extra columns at end for header and first row
                    for j,c in enumerate(["model_number", "star_age"]):
                        l = len(c)
                        fmt = " {:" + f"{l}" + "}"
                        history_file.write(fmt.format(c))
                        fmt = " {:>" + f"{l}" + "}"
                        row1 += fmt.format(j)
                if idx != 14:
                    # add data rows
                    history_file.write(row1)
                    history_file.write(row2)
                if ((idx == 8) and (len(cols) > 1)):
                    # add frist two columns of third row
                    l = len(cols[0])
                    fmt = " {:>" + f"{l}" + "}"
                    history_file.write("\n"+fmt.format(2))
                    l = len(cols[1])
                    fmt = " {:>" + f"{l}" + "}"
                    history_file.write(fmt.format(3))
        if with_profiles:
            # get final_profile.data
            with open(os.path.join(logs_path, "final_profile.data"), "w")\
                 as final_profile_file:
                for j in range(5):
                    # the first 5 lines are skipped
                    final_profile_file.write(f"Test HEADER{j}\n")
                # table with column headers and two rows of all the required
                # columns
                row1 = "\n"
                row2 = "\n"
                for j,c in enumerate(DEFAULT_PROFILE_COLS):
                    l = len(c)
                    fmt = " {:" + f"{l}" + "}"
                    final_profile_file.write(fmt.format(c))
                    fmt = " {:>" + f"{l}" + "}"
                    row1 += fmt.format(j)
                    row2 += fmt.format(j+1)
                final_profile_file.write(row1)
                final_profile_file.write(row2)
    return run_path


def get_MESA_dir(path, idx, n_runs=20):
    """Create files of one MESA run

    Parameters
    ----------
    path : str
        The path to the directory containing all the new directories and files.
    idx : int
        MESA collection index. Positive ones contain binary runs, negative ones
        single star runs, and 0 will be empty.
    n_runs : int (default: 20)
        Positive number of runs to create as a range(n_runs).

    Returns
    -------
    str
        Path to the added MESA directory.

    """
    if not os.path.exists(path):
        raise NotADirectoryError(f"Wrong path: {path}")
    if not isinstance(idx, int):
        raise TypeError("idx must be an integer")
    if not isinstance(n_runs, int):
        raise TypeError("n_runs must be an integer")
    elif n_runs <= 0:
        raise ValueError("n_runs must be positive")
    # get MESA directory
    MESA_data_path = os.path.join(path, f"MESA_data_index{idx}")
    os.mkdir(MESA_data_path)
    # check for binary or single run
    if idx > 0:
        binary_run = True
    elif idx < 0:
        binary_run = False
    else:
        return MESA_data_path
    # create the individual runs
    for run_idx in range(n_runs):
        add_MESA_run_files(MESA_data_path, run_idx, binary_run=binary_run)
    return MESA_data_path
