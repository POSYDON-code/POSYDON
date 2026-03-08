#!/usr/bin/env python3
"""
Script to evolve a few binaries.
Used for validation of the branch.

Author: Max Briel
"""

import argparse
import os
import signal
import warnings

from binary_test_cases import get_test_binaries
from formatting import line_length
from utils import print_failed_binary, print_warnings, write_binary_to_screen

from posydon.binary_evol.simulationproperties import SimulationProperties
from posydon.config import PATH_TO_POSYDON

path_to_default_params = os.path.join(PATH_TO_POSYDON, "dev-tools/script_data/1Zsun_binaries_params.ini")

def load_inlist(verbose):

    sim_prop = SimulationProperties.from_ini(path_to_default_params)
    sim_prop.load_steps(verbose=verbose, metallicity=1.0)

    return sim_prop

def evolve_binary(binary):

    # Capture warnings during evolution
    captured_warnings = []

    def warning_handler(message, category, filename, lineno, file=None, line=None):
        captured_warnings.append({
            'message': str(message),
            'category': category.__name__,
            'filename': filename,
            'lineno': lineno
        })

    # Set up warning capture
    old_showwarning = warnings.showwarning
    warnings.showwarning = warning_handler

    try:
        binary.evolve()
        # Display the evolution summary for successful evolution
        write_binary_to_screen(binary)

        # Show warnings if any were captured
        print_warnings(captured_warnings)
        print("=" * line_length)

    except Exception as e:

        # turn off binary alarm in case of exception
        signal.alarm(0)

        print_failed_binary(binary, e)

        # Show warnings if any were captured before the exception
        print_warnings(captured_warnings)

        print("=" * line_length)
    finally:
        # Always turn off binary alarm and restore warning handler
        signal.alarm(0)
        warnings.showwarning = old_showwarning


def evolve_binaries(verbose):
    """Evolves a few binaries to validate their output
    """
    sim_prop = load_inlist(verbose)
    test_binaries = get_test_binaries(sim_prop)

    for binary in test_binaries:
        evolve_binary(binary)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Evolve binaries for validation.')
    parser.add_argument('--verbose', '-v', action='store_true', default=False,
                        help='Enable verbose output (default: False)')
    args = parser.parse_args()

    evolve_binaries(verbose=args.verbose)
