"""Helper function(s) for tests requiring a POSYDON Population

used in:
    - posydon/unit_tests/popsyn/test_synthetic_population.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

import os
import h5py
import numpy as np
import pandas as pd

from posydon.popsyn.synthetic_population import Population


# helper functions

def make_ini(tmp_path,content=None):
    """
    Create a minimal dummy .ini file inside the pytest tmp_path using os.path.

    Parameters
    ----------
    tmp_path : pathlib.Path
        pytest temporary directory.
    content : str, optional
        Content to write to the ini file.

    Returns
    -------
    str
        Path (string) to the created dummy ini file.
    """
    dir_path = str(tmp_path)
    ini_path = os.path.join(dir_path, "dummy.ini")

    if content is None:
        content = "[DEFAULT]\nkey=value\n"

    with open(ini_path, "w") as f:
        f.write(content)

    return str(ini_path)

def make_test_pop(
    tmp_path,
    filename="test_population.h5",
    oneline_rows=None,
    history_rows=None,
    metallicity=0.02):
    """
    Create a minimally valid synthetic population HDF5 file and return a
    fully initialized Population object. This centralizes and standardizes
    population generation across unit tests.

    Parameters
    ----------
    tmp_path : Path-like
        Directory in which the HDF5 file will be created.
    filename : str
        Name of the file to write.
    oneline_rows : list[dict], optional
        Rows for the /oneline table.
    history_rows : list[dict], optional
        Rows for the /history table.
    metallicity : float
        Metallicty value for oneline and mass_per_metallicity tables.

    Returns
    -------
    Population
        A fully initialized Population instance.
    """

    # history and oneline tables
    
    if history_rows is None:
        history_rows = [{"binary_index": 0, "event": "start", "time": 0.0},
                        {"binary_index": 0, "event": "end",   "time": 1.0},
                        {"binary_index": 1, "event": "start", "time": 0.0},
                        {"binary_index": 1, "event": "end",   "time": 1.0}]

    if oneline_rows is None:
        oneline_rows = [{"binary_index": 0,
                        "S1_mass_i": 1.0,
                        "S2_mass_i": 1.0,
                        "state_i": "initial",
                        "metallicity": metallicity,
                        "interp_class_HMS_HMS": "initial_MT",
                        "mt_history_HMS_HMS": "Stable"},
                        {"binary_index": 1,
                        "S1_mass_i": 1.0,
                        "S2_mass_i": 1.0,
                        "state_i": "initial",
                        "metallicity": metallicity,
                        "interp_class_HMS_HMS": "no_MT",
                        "mt_history_HMS_HMS": None}]

    # Convert to DataFrames
    oneline_df = pd.DataFrame(oneline_rows).sort_values("binary_index")
    history_df = pd.DataFrame(history_rows).sort_values(["binary_index", "time"])

    # history_lengths = number of rows per binary_index
    history_lengths_df = history_df.groupby("binary_index").size().to_frame("length")

    # ini_parameters
    ini_df = pd.DataFrame({
        "Parameter": ["metallicity", "number_of_binaries"],
        "Value": [metallicity, len(oneline_df)],
    })

    # mass_per_metallicity
    mass_df = pd.DataFrame(
        {"simulated_mass": [1.0], "number_of_systems": [len(oneline_df)]},
        index=[metallicity]
    )

    # Write HDF5 file using pandas/HDFStore (Population expects PyTables layout)
    fpath = os.path.join(tmp_path, filename)
    with pd.HDFStore(fpath, "w") as store:
        store.put("oneline", oneline_df, format="table")
        store.put("history", history_df, format="table")
        store.put("history_lengths", history_lengths_df, format="table")
        store.put("ini_parameters", ini_df, format="table")
        store.put("mass_per_metallicity", mass_df, format="table")

    # Return fully initialized Population object
    return Population(fpath)
