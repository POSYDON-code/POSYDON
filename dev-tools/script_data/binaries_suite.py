#!/usr/bin/env python3
"""
Script to evolve a set of test binaries at a given metallicity.
Used for validation of POSYDON branches.

Usage:
    python binaries_suite.py --output results.h5 --metallicity 1
    python binaries_suite.py --output results.h5 --metallicity 0.45 --verbose

Authors: Max Briel, Elizabeth Teng
"""

import argparse
import os
import sys
import warnings
import pandas as pd

from posydon.binary_evol.binarystar import BinaryStar, SingleStar
from posydon.binary_evol.simulationproperties import SimulationProperties
from posydon.popsyn.io import simprop_kwargs_from_ini
from posydon.utils.common_functions import orbital_separation_from_period

AVAILABLE_METALLICITIES = [2., 1., 0.45, 0.2, 0.1, 0.01, 0.001, 0.0001]

# Display settings
TARGET_ROWS = 12
LINE_LENGTH = 80
COLUMNS_TO_SHOW = ['step_names', 'state', 'event', 
                   'S1_state', 'S1_mass',
                   'S2_state', 'S2_mass', 'orbital_period']

def load_inlist(metallicity, verbose, ini_path=None):
    """Load simulation properties from ini file and configure for given metallicity.

    Args:
        metallicity: float, metallicity in solar units
        verbose: bool
        ini_path: str, path to ini file (auto-detected if None)

    Returns:
        SimulationProperties object with loaded steps
    """
    if ini_path is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        ini_path = os.path.join(script_dir, 'binaries_params.ini')

    if not os.path.exists(ini_path):
        raise FileNotFoundError(f"INI file not found: {ini_path}")

    sim_kwargs = simprop_kwargs_from_ini(ini_path, verbose=verbose)

    metallicity_kwargs = {'metallicity': metallicity, 'verbose': verbose}

    metallicity_steps = ['step_HMS_HMS', 'step_CO_HeMS', 'step_CO_HMS_RLO',
                        'step_CO_HeMS_RLO', 'step_detached', 'step_disrupted',
                        'step_merged', 'step_initially_single']
    for step_name in metallicity_steps:
        if step_name in sim_kwargs:
            sim_kwargs[step_name][1].update(metallicity_kwargs)

    sim_prop = SimulationProperties(**sim_kwargs)
    sim_prop.load_steps(verbose=verbose)
    return sim_prop

def write_binary_to_screen(binary):
    """Writes a binary DataFrame to screen."""
    df = binary.to_df(**{'extra_columns':{'step_names':'str'}})

    # Filter to only existing columns
    available_columns = [col for col in COLUMNS_TO_SHOW if col in df.columns]
    df_filtered = df[available_columns].reset_index(drop=True)

    print("=" * LINE_LENGTH)

    # Print the DataFrame
    df_string = df_filtered.to_string(index=True, float_format='%.3f')
    print(df_string)

    # Add empty lines to reach exactly 10 rows of output
    current_rows = len(df_filtered) + 1 # add one for header

    if current_rows < TARGET_ROWS:
        for _ in range(TARGET_ROWS - current_rows):
            print("")

    print("-" * LINE_LENGTH)


def print_failed_binary(binary,e,  max_error_lines=3):
    """Print information about a binary that failed to evolve."""

    print("=" * LINE_LENGTH)
    print(f"🚨 Binary Evolution Failed!")
    print(f"Exception: {type(e).__name__}")
    print(f"Message: {e}")

    # Get the binary's current state and limit output
    try:
        df = binary.to_df(**{'extra_columns':{'step_names':'str'}})
        if len(df) > 0:
            # Select only the desired columns

            available_columns = [col for col in COLUMNS_TO_SHOW if col in df.columns]
            df_filtered = df[available_columns].reset_index(drop=True)

            # Limit to max_error_lines
            if len(df_filtered) > max_error_lines:
                df_filtered = df_filtered.tail(max_error_lines)
                print(f"\nShowing last {max_error_lines} evolution steps before failure:")
            else:
                print(f"\nEvolution steps before failure ({len(df_filtered)} steps):")

            df_string = df_filtered.to_string(index=True, float_format='%.3f')
            print(df_string)

            current_rows = len(df_filtered) + 1 + 5  # add one for header
            empty_lines_needed = TARGET_ROWS - current_rows
            for _ in range(max(0, empty_lines_needed)):
                print("")
        else:
            print("\nNo evolution steps recorded before failure.")
    except Exception as inner_e:
        print(f"\nCould not retrieve binary state: {inner_e}")

    print("-" * LINE_LENGTH)

def evolve_binary(binary, h5file, binary_id):
    """
    Evolves a single binary, prints its evolution, and saves to HDF5.

    Args:
        binary: BinaryStar object
        h5file: open pd.HDFStore object for writing
        binary_id: unique identifier for this binary
    """

    # Capture warnings during evolution
    captured_warnings = []

    def warning_handler(message, category, filename, lineno, file=None, line=None):
        captured_warnings.append({
            "binary_id": int(binary_id),
            "category": category.__name__,
            "message": str(message),
            "filename": filename,
            "lineno": lineno
        })

    # Set up warning capture
    old_showwarning = warnings.showwarning
    warnings.showwarning = warning_handler
    
    print(f"Binary {binary_id}")
    evolution_df = None 

    try:
        binary.evolve()
        write_binary_to_screen(binary)
        evolution_df = binary.to_df(extra_columns={'step_names':'str'})

    except Exception as e:
        print_failed_binary(binary, e)

        # If evolution fails, create a minimal df with error info
        evolution_df = pd.DataFrame([{
            "binary_id": int(binary_id),
            "exception_type": type(e).__name__,
            "exception_message": str(e)
        }])

    finally:
        warnings.showwarning = old_showwarning
        
        # Ensure we always have a dataframe
        if evolution_df is not None:
            # Decode bytes columns if needed
            for col in evolution_df.select_dtypes([object]):
                if evolution_df[col].apply(lambda x: isinstance(x, bytes)).any():
                    evolution_df[col] = evolution_df[col].apply(
                        lambda x: x.decode('utf-8') if isinstance(x, bytes) else x
                    )

            # Always ensure binary_id exists
            if "binary_id" not in evolution_df.columns:
                evolution_df["binary_id"] = int(binary_id)

            # Defragment
            evolution_df = evolution_df.copy()

            # Determine min_itemsize from the dataframe we're actually saving
            string_cols = evolution_df.select_dtypes([object]).columns
            min_itemsize = {col: 100 for col in string_cols}
            h5file.append("evolution", evolution_df, format="table",
                          data_columns=True, min_itemsize=min_itemsize)

        # Save warnings
        if captured_warnings:
            warn_df = pd.DataFrame(captured_warnings)
            # Ensure consistent string column sizes for warnings table
            warn_string_cols = warn_df.select_dtypes([object]).columns
            warn_min_itemsize = {col: 200 for col in warn_string_cols}
            h5file.append("warnings", warn_df, format="table",
                          min_itemsize=warn_min_itemsize)
            print(f"⚠️  {len(captured_warnings)} warning(s) raised during evolution:")
            for i, w in enumerate(captured_warnings[:3], 1):
                print(f"   {i}. {w['category']}: {w['message'][:80]}")
            if len(captured_warnings) > 3:
                print(f"   ... and {len(captured_warnings) - 3} more warning(s)")
        else:
            print(f"  No warning(s) raised during evolution")

        print(f"✅ Finished binary {binary_id}")
        print("=" * LINE_LENGTH)
        
        
        
def get_test_binaries(metallicity, sim_prop):
    """Return the list of test binaries as (star1_kwargs, star2_kwargs, binary_kwargs, description) tuples.

    All binaries use the specified metallicity.
    """
    Z = metallicity

    binaries = [
        # 0: Failing binary in matching
        ({'mass': 11.948472796094759, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [231.97383621190582, 5.927334890264575, 1.5990566013567014, 6.137994236518587]},
         {'mass': 7.636958434479617, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 190925.99636740884, 'eccentricity': 0.0},
         "Failing binary in matching"),

        # 1: Failing binary in matching
        ({'mass': 30.169861921689556, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [77.96834852144123, 0.05021460132555987, 2.3146518208348152, 1.733054979982291]},
         {'mass': 10.972734402996027, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 20479.71919353725, 'eccentricity': 0.0},
         "Failing binary in matching"),

        # 2: Flipped S1 and S2 (near-equal mass)
        ({'mass': 9.474917413943635, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [133.5713935237759, 4.398754864537542, 2.703102872841114, 1.4633904612711142]},
         {'mass': 9.311073918196263, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 18.605997832086413, 'eccentricity': 0.0},
         "Flipped S1 and S2 (near-equal mass)"),

        # 3: Flipped S1 and S2 (high mass ratio)
        ({'mass': 10.438541, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 1.400713, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 4.190728383757787, 1.1521129697118118, 5.015343794234789]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 9.824025, 'eccentricity': 0.0},
         "Flipped S1 and S2 (high mass ratio)"),

        # 4: Flipped S1 and S2
        ({'mass': 9.845907, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 4.190728383757787, 1.1521129697118118, 5.015343794234789]},
         {'mass': 9.611029, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 3.820571, 'eccentricity': 0.0},
         "Flipped S1 and S2"),

        # 5: Normal binary evolution (high mass)
        ({'mass': 30.845907, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 4.190728383757787, 1.1521129697118118, 5.015343794234789]},
         {'mass': 30.611029, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 30.820571, 'eccentricity': 0.0},
         "Normal binary evolution (high mass)"),

        # 6: Normal binary (wide orbit)
        ({'mass': 9.213534679594247, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [327.5906384501521, 1.7707176050073297, 1.573225822966838, 1.6757313876001914]},
         {'mass': 7.209878522799272, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 63123.74544474666, 'eccentricity': 0.0},
         "Normal binary (wide orbit)"),

        # 7: Normal binary (near-equal mass, close)
        ({'mass': 9.561158487732602, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [317.5423844462847, 2.9095984678057603, 1.754121288652108, 2.3693917842468784]},
         {'mass': 9.382732464319286, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 27.77657038557851, 'eccentricity': 0.0},
         "Normal binary (near-equal mass, close)"),

        # 8: Normal binary
        ({'mass': 7.552858, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [40.91509926587841, 2.6295454150818256, 1.6718337470964977, 6.0408769315244895]},
         {'mass': 6.742063, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 17.957531550841225, 'eccentricity': 0.0},
         "Normal binary"),

        # 9: High BH spin options
        ({'mass': 31.616785, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [10, 4.190728383757787, 1.1521129697118118, 5.015343794234789]},
         {'mass': 26.874267, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 501.99252706449792, 'eccentricity': 0.0},
         "High BH spin options"),

        # 10: Original a>1 spin error
        ({'mass': 18.107506844123645, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [528.2970725443025, 4.190728383757787, 1.1521129697118118, 5.015343794234789]},
         {'mass': 15.641392951875442, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 151.99252706449792, 'eccentricity': 0.0},
         "Original a>1 spin error"),

        # 11: FIXED disrupted crash
        ({'mass': 52.967313, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 36.306444, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 12.877004, 'eccentricity': 0.0},
         "FIXED disrupted crash"),

        # 12: FIXED error with SN type
        ({'mass': 17.782576, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 3.273864, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 4513.150157, 'eccentricity': 0.0},
         "FIXED error with SN type"),

        # 13: FIXED oRLO2 looping
        ({'mass': 170.638207, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [4.921294, 4.31745, 1.777768, 3.509656]},
         {'mass': 37.917852, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 113.352736, 'eccentricity': 0.0},
         "FIXED oRLO2 looping"),

        # 14: Redirect to step_CO_HeMS
        ({'mass': 8.333579, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [17.125568, 4.101834, 0.917541, 3.961291]},
         {'mass': 8.208376, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 66.870417, 'eccentricity': 0.0},
         "Redirect to step_CO_HeMS"),

        # 15: FIXED oRLO2 looping
        ({'mass': 16.921378, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [268.837139, 5.773527, 2.568105, 2.519068]},
         {'mass': 16.286318, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 37.958768, 'eccentricity': 0.0},
         "FIXED oRLO2 looping"),

        # 16: FIXED step_detached failure
        ({'mass': 19.787769, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [24.464803, 0.666314, 1.954698, 5.598975]},
         {'mass': 7.638741, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 3007.865561, 'eccentricity': 0.0},
         "FIXED step_detached failure"),

        # 17: Disrupted binary
        ({'mass': 16.921378, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [268.837139, 5.773527, 2.568105, 2.519068]},
         {'mass': 16.286318, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 3007.865561, 'eccentricity': 0.0},
         "Disrupted binary"),

        # 18: FIXED Detached binary failure (low mass)
        ({'mass': 9, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 0.8, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 4513.150157, 'eccentricity': 0.0},
         "FIXED Detached binary failure (low mass)"),

        # 19: FIXED SN_TYPE = None crash
        ({'mass': 17.782576, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 3.273864, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 4513.150157, 'eccentricity': 0.0},
         "FIXED SN_TYPE = None crash"),

        # 20: FIXED SN_TYPE errors (low mass primary)
        ({'mass': 6.782576, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 3.273864, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 4513.150157, 'eccentricity': 0.0},
         "FIXED SN_TYPE errors (low mass primary)"),

        # 21: FIXED SN_TYPE errors (high mass)
        ({'mass': 40.638207, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [30.921294, 4.31745, 1.777768, 3.509656]},
         {'mass': 37.917852, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 2113.352736, 'eccentricity': 0.0},
         "FIXED SN_TYPE errors (high mass)"),

        # 22: FIXED ECSN errors
        ({'mass': 12.376778, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [80, 4.31745, 1.777768, 3.509656]},
         {'mass': 9.711216, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 79.83702, 'eccentricity': 0.0},
         "FIXED ECSN errors"),

        # 23: Interpolator masses (close)
        ({'mass': 7.592921, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 5.038679, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 5.537807, 'eccentricity': 0.0},
         "Interpolator masses (close)"),

        # 24: Interpolator masses (both kicked)
        ({'mass': 38.741115, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [21.113771, 2.060135, 2.224789, 4.089729]},
         {'mass': 27.776178, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [282.712103, 0.296252, 1.628433, 5.623812]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 93.387072, 'eccentricity': 0.0},
         "Interpolator masses (both kicked)"),

        # 25: FIXED NaN spin (very high mass, wide)
        ({'mass': 70.066924, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 34.183110, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 5.931492e+03, 'eccentricity': 0.0},
         "FIXED NaN spin (very high mass, wide)"),

        # 26: FIXED NaN spin
        ({'mass': 28.837286, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 6.874867, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 35.609894, 'eccentricity': 0.0},
         "FIXED NaN spin"),

        # 27: oRLO2 issue (near-equal mass)
        ({'mass': 29.580210, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 28.814626, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 40.437993, 'eccentricity': 0.0},
         "oRLO2 issue (near-equal mass)"),

        # 28: oRLO2 issue (high mass ratio)
        ({'mass': 67.126795, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 19.622908, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 1484.768582, 'eccentricity': 0.0},
         "oRLO2 issue (high mass ratio)"),

        # 29: oRLO2 issue (very high mass, near-equal)
        ({'mass': 58.947503, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 56.660506, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 2011.300659, 'eccentricity': 0.0},
         "oRLO2 issue (very high mass, near-equal)"),

        # 30: oRLO2 issue (extreme mass, kicked)
        ({'mass': 170.638207, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [47.979957374424956, 5.317304576107798, 2.7259013166068145, 4.700929589520818]},
         {'mass': 37.917852, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 113.352736, 'eccentricity': 0.0},
         "oRLO2 issue (extreme mass, kicked)"),

        # 31: oRLO2 issue (very high mass, close)
        ({'mass': 109.540207, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 84.344530, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 5.651896, 'eccentricity': 0.0},
         "oRLO2 issue (very high mass, close)"),

        # 32: Redirect (extreme mass ratio)
        ({'mass': 13.889634, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 0.490231, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 14513.150157, 'eccentricity': 0.0},
         "Redirect (extreme mass ratio)"),

        # 33: Redirect (low mass secondary)
        ({'mass': 9, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 0.8, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 4513.150157, 'eccentricity': 0.0},
         "Redirect (low mass secondary)"),

        # 34: Max time (very high mass, wide)
        ({'mass': 103.07996766780799, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.2965418610971261, 2.0789170290719117, 3.207488023705968]},
         {'mass': 83.66522615073987, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 1449.1101985875678, 'eccentricity': 0.0},
         "Max time (very high mass, wide)"),

        # 35: Max time
        ({'mass': 8.860934140643465, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [11.818027275431337, 2.812412688633058, 0.4998731824233789, 2.9272630485628643]},
         {'mass': 8.584716012668551, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 20.82030114750744, 'eccentricity': 0.0},
         "Max time"),

        # 36: PR421
        ({'mass': 24.035366, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 23.187355, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 18.865029, 'eccentricity': 0.0},
         "PR421"),

        # 37: CE class
        ({'mass': 33.964274, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 28.98149, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 82.370989, 'eccentricity': 0.0},
         "CE class"),

        # 38: PR574 - stepCE fix
        ({'mass': 29.580210, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 28.814626 * 0.4, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 300.437993, 'eccentricity': 0.0},
         "PR574 - stepCE fix"),

        # 39: e_ZAMS error
        ({'mass': 8.161885721822461, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 3.5907829421526154, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 36.873457164644144, 'eccentricity': 0.0},
         "e_ZAMS error"),

        # 40: e_ZAMS error (eccentric)
        ({'mass': 35.24755025317775, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [19.755993125895806, 0.37149222852233904, 1.6588846085306563, 1.434617029858906]},
         {'mass': 30.000450298072902, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 24060.02101364665, 'eccentricity': 0.8085077857996965},
         "e_ZAMS error (eccentric)"),

        # 41: e_ZAMS error (high mass ratio)
        ({'mass': 11.862930493162692, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 1.4739109294156703, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 4111.083887312003, 'eccentricity': 0.0},
         "e_ZAMS error (high mass ratio)"),

        # 42: e_ZAMS error (extreme mass ratio)
        ({'mass': 8.527361341212108, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 0.7061748406821822, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 2521.1927287891444, 'eccentricity': 0.0},
         "e_ZAMS error (extreme mass ratio)"),

        # 43: e_ZAMS error
        ({'mass': 13.661942533447398, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'mass': 4.466151109802313, 'state': 'H-rich_Core_H_burning', 'metallicity': Z,
          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]},
         {'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
          'orbital_period': 3110.1346707516914, 'eccentricity': 0.0},
         "e_ZAMS error"),
    ]

    return binaries


def evolve_binaries(metallicity, verbose, output_path, ini_path=None):
    """Evolves the test binary suite at the given metallicity and saves results.

    Args:
        metallicity: float, metallicity in solar units
        verbose: bool
        output_path: str, path to save HDF5 output
        ini_path: str, path to ini file (auto-detected if None)
    """
    print(f"{'=' * LINE_LENGTH}")
    print(f"  Evolving test binaries at Z = {metallicity} Zsun")
    print(f"  Output: {output_path}")
    print(f"{'=' * LINE_LENGTH}\n")

    sim_prop = load_inlist(metallicity, verbose, ini_path)
    test_binaries = get_test_binaries(metallicity, sim_prop)

    with pd.HDFStore(output_path, mode="w") as h5file:
        # Save metadata
        meta_df = pd.DataFrame([{
            'metallicity': metallicity,
            'n_binaries': len(test_binaries),
        }])
        h5file.put("metadata", meta_df, format="table")

        for binary_id, (s1_kw, s2_kw, bin_kw, description) in enumerate(test_binaries):
            print(f"\n[{binary_id}/{len(test_binaries)-1}] {description}")

            star_1 = SingleStar(**s1_kw)
            star_2 = SingleStar(**s2_kw)

            # Add separation from period if not explicitly provided
            if 'separation' not in bin_kw and 'orbital_period' in bin_kw:
                bin_kw['separation'] = orbital_separation_from_period(
                    bin_kw['orbital_period'], star_1.mass, star_2.mass
                )

            binary = BinaryStar(star_1, star_2, **bin_kw, properties=sim_prop)
            evolve_binary(binary, h5file, binary_id)

    print(f"\n{'=' * LINE_LENGTH}")
    print(f"  All {len(test_binaries)} binaries complete. Results saved to {output_path}")
    print(f"{'=' * LINE_LENGTH}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Evolve test binaries for POSYDON branch validation.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--verbose', '-v', action='store_true', default=False,
                        help='Enable verbose output')
    parser.add_argument('--output', type=str, required=True,
                        help='Path to save HDF5 output')
    parser.add_argument('--metallicity', '-Z', type=float, default=1.0,
                        help=f'Metallicity in solar units. Available: {AVAILABLE_METALLICITIES}')
    parser.add_argument('--ini', type=str, default=None,
                        help='Path to params ini file (auto-detected if not given)')
    args = parser.parse_args()

    if args.metallicity not in AVAILABLE_METALLICITIES:
        print(f"WARNING: Metallicity {args.metallicity} not in standard set {AVAILABLE_METALLICITIES}.")
        print(f"Proceeding anyway, but POSYDON grids may not exist for this value.")

    evolve_binaries(
        metallicity=args.metallicity,
        verbose=args.verbose,
        output_path=args.output,
        ini_path=args.ini,
    )