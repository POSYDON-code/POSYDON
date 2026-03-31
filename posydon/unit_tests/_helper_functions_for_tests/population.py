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

from posydon.popsyn.rate_calculation import (
    DEFAULT_SFH_MODEL,
    get_cosmic_time_from_redshift,
    get_redshift_bin_centers,
)
from posydon.popsyn.synthetic_population import Population, Rates, TransientPopulation

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

    # ini_parameters – include all keys that _load_ini_params expects
    ini_params = {
        "metallicities": [metallicity],
        "number_of_binaries": len(oneline_df),
        "binary_fraction_scheme": "const",
        "binary_fraction_const": 0.7,
        "star_formation": "burst",
        "max_simulation_time": 13800000000.0,
        "primary_mass_scheme": "Kroupa2001",
        "primary_mass_min": 0.01,
        "primary_mass_max": 200.0,
        "secondary_mass_scheme": "flat_mass_ratio",
        "secondary_mass_min": 0.0005,
        "secondary_mass_max": 200.0,
        "orbital_scheme": "period",
        "orbital_period_scheme": "Sana+12_period_extended",
        "orbital_period_min": 0.35,
        "orbital_period_max": 6000.0,
        "orbital_separation_scheme": "log_uniform",
        "orbital_separation_min": 5.0,
        "orbital_separation_max": 100000.0,
        "eccentricity_scheme": "zero",
        "posydon_version": "test",
    }
    ini_df = pd.DataFrame({k: [v] for k, v in ini_params.items()})


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

def make_test_transient_pop(
    tmp_path,
    transient_name="test_transient",
    filename="test_population.h5",
    oneline_rows=None,
    history_rows=None,
    transient_rows=None,
    metallicity=0.02):
    """
    Create a minimally valid TransientPopulation HDF5 file and return
    a fully initialized TransientPopulation object.

    Builds on make_test_pop by adding a /transients/{transient_name} table.

    Parameters
    ----------
    tmp_path : Path-like
        Directory in which the HDF5 file will be created.
    transient_name : str
        Name for the transient population key.
    filename : str
        Name of the file to write.
    oneline_rows : list[dict], optional
        Rows for the /oneline table.
    history_rows : list[dict], optional
        Rows for the /history table.
    transient_rows : list[dict], optional
        Rows for the /transients/{transient_name} table.
    metallicity : float
        Metallicity value.

    Returns
    -------
    TransientPopulation
        A fully initialized TransientPopulation instance.
    """
    pop = make_test_pop(
        tmp_path,
        filename=filename,
        oneline_rows=oneline_rows,
        history_rows=history_rows,
        metallicity=metallicity,
    )

    if transient_rows is None:
        transient_rows = [
            {"time": 100.0, "metallicity": metallicity, "channel": "ch_A"},
            {"time": 200.0, "metallicity": metallicity, "channel": "ch_B"},
        ]

    transient_df = pd.DataFrame(transient_rows)

    with pd.HDFStore(pop.filename, "a") as store:
        store.append(
            "transients/" + transient_name,
            transient_df,
            format="table",
            min_itemsize={"channel": 100},
        )

    return TransientPopulation(
        pop.filename, transient_name, verbose=False
    )


def make_test_rates(
    tmp_path,
    transient_name="test_transient",
    SFH_identifier="test_SFH",
    filename="test_population.h5",
    oneline_rows=None,
    history_rows=None,
    transient_rows=None,
    metallicity=0.02,
    MODEL=None):
    """
    Create a minimally valid Rates HDF5 file and return
    a fully initialized Rates object.

    Builds on make_test_transient_pop by adding the rates tables:
    MODEL, weights, z_events, and birth.

    Parameters
    ----------
    tmp_path : Path-like
        Directory in which the HDF5 file will be created.
    transient_name : str
        Name for the transient population key.
    SFH_identifier : str
        Name for the star formation history identifier.
    filename : str
        Name of the file to write.
    oneline_rows : list[dict], optional
        Rows for the /oneline table.
    history_rows : list[dict], optional
        Rows for the /history table.
    transient_rows : list[dict], optional
        Rows for the /transients/{transient_name} table.
    metallicity : float
        Metallicity value.
    MODEL : dict, optional
        The SFH model dict. If None, uses DEFAULT_SFH_MODEL.

    Returns
    -------
    Rates
        A fully initialized Rates instance.
    """
    tpop = make_test_transient_pop(
        tmp_path,
        transient_name=transient_name,
        filename=filename,
        oneline_rows=oneline_rows,
        history_rows=history_rows,
        transient_rows=transient_rows,
        metallicity=metallicity,
    )

    if MODEL is None:
        MODEL = dict(DEFAULT_SFH_MODEL)
    else:
        # Merge user overrides into a copy of the defaults
        merged = dict(DEFAULT_SFH_MODEL)
        merged.update(MODEL)
        MODEL = merged

    n_transients = len(pd.read_hdf(tpop.filename, "transients/" + transient_name))

    # Compute birth bins from MODEL
    z_birth = get_redshift_bin_centers(MODEL["delta_t"])
    t_birth = get_cosmic_time_from_redshift(z_birth)
    nr_of_birth_bins = len(z_birth)

    base_path = "/transients/" + transient_name + "/rates/" + SFH_identifier + "/"

    with pd.HDFStore(tpop.filename, "a") as store:
        # MODEL table — mirror _write_MODEL_data logic for dlogZ lists
        if (MODEL["dlogZ"] is not None) and (not isinstance(MODEL["dlogZ"], float)):
            store.put(base_path + "MODEL", pd.DataFrame(MODEL))
        else:
            store.put(base_path + "MODEL", pd.DataFrame(MODEL, index=[0]))

        # birth table
        store.put(base_path + "birth", pd.DataFrame({"z": z_birth, "t": t_birth}))

        # weights table: (n_transients x nr_of_birth_bins) with small dummy values
        weights = np.full((n_transients, nr_of_birth_bins), 1e-10)
        store.append(
            base_path + "weights",
            pd.DataFrame(data=weights, index=np.arange(n_transients)),
            format="table",
        )

        # z_events table: same shape, dummy redshift values
        z_events = np.full((n_transients, nr_of_birth_bins), 0.1)
        store.append(
            base_path + "z_events",
            pd.DataFrame(data=z_events, index=np.arange(n_transients)),
            format="table",
        )

    return Rates(
        tpop.filename, transient_name, SFH_identifier, verbose=False
    )
