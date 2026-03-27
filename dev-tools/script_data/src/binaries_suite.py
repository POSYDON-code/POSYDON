"""
Script to evolve a few binaries.
Used for validation of the branch.

Author: Max Briel
"""

import argparse
import os
import shutil
import signal
import warnings

import numpy as np
import pandas as pd
from binary_test_cases import get_test_binaries
from formatting import AVAILABLE_METALLICITIES, LINE_LENGTH
from utils import print_failed_binary, print_warnings, write_binary_to_screen

from posydon.binary_evol.binarystar import BinaryStar, SingleStar
from posydon.binary_evol.simulationproperties import SimulationProperties
from posydon.config import PATH_TO_POSYDON, PATH_TO_POSYDON_DATA
from posydon.utils.common_functions import orbital_separation_from_period

#base_dir =os.path.dirname(PATH_TO_POSYDON)
#script_dir = os.path.join(PATH_TO_POSYDON, "dev-tools/script_data/")
#path_to_default_params = os.path.join(script_dir, "inlists/default_test_params.ini")

def load_inlist(ini_path, metallicity, verbose):

    if ini_path is None:
        # copy the .ini file from the POSYDON version installed in current environment
        default_ini_fn = os.path.join(PATH_TO_POSYDON, "posydon/popsyn/population_params_default.ini")
        test_ini_fn = os.path.join(PATH_TO_POSYDON, "dev-tools/script_data/inlists/test_params.ini")
        ini_path = shutil.copyfile(default_ini_fn, test_ini_fn)

    print(f"Reading inlist: {ini_path}")
    sim_prop = SimulationProperties.from_ini(ini_path)

    # TODO: try to create/pass RNG with a fixed seed
    RNG = np.random.default_rng(0)
    try:
        sim_prop.load_steps(verbose=verbose, RNG=RNG, metallicity=metallicity)
    except TypeError as e:
        sim_prop.load_steps(verbose=verbose, metallicity=metallicity)

    return sim_prop

def evolve_binary(binary, binary_id):

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

    print(f"Binary {binary_id}")
    evolution_df = None
    error_df = None

    # Set up warning capture
    old_showwarning = warnings.showwarning
    warnings.showwarning = warning_handler

    try:
        binary.evolve()
        # Display the evolution summary for successful evolution
        write_binary_to_screen(binary)
        evolution_df = binary.to_df(extra_columns={'step_names':'str'})

        # Show warnings if any were captured
        print_warnings(captured_warnings)
        print("=" * LINE_LENGTH)

    except Exception as e:

        # turn off binary alarm in case of exception
        signal.alarm(0)

        print_failed_binary(binary, e)
        error_df = pd.DataFrame([{
            "binary_id": int(binary_id),
            "exception_type": type(e).__name__,
            "exception_message": str(e)
        }])

        # Show warnings if any were captured before the exception
        print_warnings(captured_warnings)

        print("=" * LINE_LENGTH)
    finally:
        # Always turn off binary alarm and restore warning handler
        signal.alarm(0)
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

            # Defragment the DataFrame from POSYDON's column-by-column construction
            evolution_df = evolution_df.copy()

        # Save warnings
        if captured_warnings:
            print(f"⚠️  {len(captured_warnings)} warning(s) raised during evolution:")
            for i, w in enumerate(captured_warnings[:3], 1):
                print(f"   {i}. {w['category']}: {w['message'][:80]}")
            if len(captured_warnings) > 3:
                print(f"   ... and {len(captured_warnings) - 3} more warning(s)")
        else:
            print(f"  No warning(s) raised during evolution")

        print(f"✅ Finished binary {binary_id}")
        print("=" * LINE_LENGTH)

    return evolution_df, error_df, captured_warnings

def create_binary(s1_kw, s2_kw, bin_kw, sim_prop):

        star_1 = SingleStar(**s1_kw)
        star_2 = SingleStar(**s2_kw)

        # Add separation from period if not explicitly provided
        if 'separation' not in bin_kw and 'orbital_period' in bin_kw:
            bin_kw['separation'] = orbital_separation_from_period(
                bin_kw['orbital_period'], star_1.mass, star_2.mass
            )

        return BinaryStar(star_1, star_2, **bin_kw, properties=sim_prop)

def evolve_binaries(metallicity, output_path, verbose, ini_path=None):
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

    sim_prop = load_inlist(ini_path, metallicity, verbose)
    test_binaries = get_test_binaries(metallicity)

    # Collect all results in memory, then write once at the end.
    # This avoids repeated HDFStore.append() calls, each of which
    # reconciles schemas, checks string sizing, and flushes to disk.
    all_evolution_dfs = []
    all_error_dfs = []
    all_warning_dfs = []

    for binary_id, (s1_kw, s2_kw, bin_kw, description) in enumerate(test_binaries):
        print(f"\n[{binary_id}/{len(test_binaries)-1}] {description}")

        binary = create_binary(s1_kw, s2_kw, bin_kw, sim_prop)

        evo_df, err_df, warn_list = evolve_binary(binary, binary_id)

        if evo_df is not None:
            all_evolution_dfs.append(evo_df)
        if err_df is not None:
            all_error_dfs.append(err_df)
        if warn_list:
            all_warning_dfs.append(pd.DataFrame(warn_list))

    # ── Completeness check ──────────────────────────────────────────
    expected_ids = set(range(len(test_binaries)))
    evolved_ids = set()
    errored_ids = set()

    for df in all_evolution_dfs:
        if 'binary_id' in df.columns:
            evolved_ids.update(df['binary_id'].unique())
    for df in all_error_dfs:
        if 'binary_id' in df.columns:
            errored_ids.update(df['binary_id'].unique())

    accounted_ids = evolved_ids | errored_ids
    missing_ids = sorted(expected_ids - accounted_ids)

    if missing_ids:
        print(f"\n⚠️  WARNING: {len(missing_ids)} binary(ies) unaccounted for: {missing_ids}")
        print(f"  These produced neither evolution output nor a caught error.")

    # ── Single-pass HDF5 write ──────────────────────────────────────────
    with pd.HDFStore(output_path, mode="w") as h5file:
        # Save metadata
        meta_df = pd.DataFrame([{
            'metallicity': metallicity,
            'n_binaries': len(test_binaries),
            'n_evolved': len(evolved_ids),
            'n_errored': len(errored_ids),
            'n_missing': len(missing_ids),
            'missing_ids': str(missing_ids) if missing_ids else '',
            'path_to_posydon_data': PATH_TO_POSYDON_DATA,
        }])
        h5file.put("metadata", meta_df, format="table")

        if all_evolution_dfs:
            combined_evo = pd.concat(all_evolution_dfs, ignore_index=True)
            string_cols = combined_evo.select_dtypes([object]).columns
            min_itemsize = {col: 500 for col in string_cols}
            h5file.put("evolution", combined_evo, format="table",
                       data_columns=True, min_itemsize=min_itemsize)

        if all_error_dfs:
            combined_err = pd.concat(all_error_dfs, ignore_index=True)
            err_string_cols = combined_err.select_dtypes([object]).columns
            err_min_itemsize = {col: 1000 for col in err_string_cols}
            h5file.put("errors", combined_err, format="table",
                       min_itemsize=err_min_itemsize)

        if all_warning_dfs:
            combined_warn = pd.concat(all_warning_dfs, ignore_index=True)
            warn_string_cols = combined_warn.select_dtypes([object]).columns
            warn_min_itemsize = {col: 1000 for col in warn_string_cols}
            h5file.put("warnings", combined_warn, format="table",
                       min_itemsize=warn_min_itemsize)

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
