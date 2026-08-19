"""Configuration handling for the posydon grid setup."""

import os

import pandas

from posydon.grids.psygrid import PSyGrid
from posydon.utils import configfile
from posydon.utils.posydonwarning import Pwarn


def parse_config_file(inifile_path: str) -> tuple[dict, dict, dict, dict]:
    """Parse an inifile into the four configuration sections.

    Args:
        inifile_path: Path to the inifile to parse.

    Returns:
        The (run_parameters, slurm, mesa_inlists, mesa_extras) dictionaries.
    """
    return configfile.parse_inifile(inifile_path)


def normalize_defaults(run_parameters: dict, mesa_inlists: dict) -> None:
    """Add missing default values to the parsed dictionaries in place.

    Args:
        run_parameters: The run_parameters section of the inifile.
        mesa_inlists: The mesa_inlists section of the inifile.
    """
    if 'scenario' not in mesa_inlists.keys():
        mesa_inlists['scenario'] = None

    if 'keep_profiles' not in run_parameters.keys():
        run_parameters['keep_profiles'] = False

    if 'keep_photos' not in run_parameters.keys():
        run_parameters['keep_photos'] = False


def validate_config(run_parameters: dict, grid_type: str) -> None:
    """Validate the parsed run_parameters against the grid type.

    Args:
        run_parameters: The run_parameters section of the inifile.
        grid_type: Either 'fixed' or 'dynamic'.

    Raises:
        ValueError: If the grid path does not exist or a dynamic grid is
            missing its psycris inifile.
    """
    if ((not os.path.isfile(run_parameters['grid'])) and (not os.path.isdir(run_parameters['grid']))):
        raise ValueError("Supplied grid does not exist, please check your path and try again")

    if ('psycris_inifile' not in run_parameters.keys()) and (grid_type == 'dynamic'):
        raise ValueError("Please add psycris inifile to the [run_parameters] section of the inifile.")


def read_grid_file(filepath: str) -> tuple[pandas.DataFrame, str]:
    """Read a grid from a csv, h5, or directory of MESA runs.

    Args:
        filepath: Path to the grid file or directory of MESA runs.

    Returns:
        The grid DataFrame and the name of the (fixed) grid file.

    Raises:
        ValueError: If the grid format is not recognized.
    """
    if '.csv' in filepath:
        grid_df = pandas.read_csv(filepath)
        fixgrid_file_name = filepath
    elif '.h5' in filepath:
        psy_grid = PSyGrid()
        psy_grid.load(filepath)
        grid_df = psy_grid.get_pandas_initial_final()
        psy_grid.close()
        fixgrid_file_name = filepath
    elif os.path.isdir(filepath):
        PSyGrid().create(filepath, "./fixed_grid_results.h5", slim=True)
        psy_grid = PSyGrid()
        psy_grid.load("./fixed_grid_results.h5")
        grid_df = psy_grid.get_pandas_initial_final()
        psy_grid.close()
        fixgrid_file_name = os.path.join(os.getcwd(), "fixed_grid_results.h5")
    else:
        raise ValueError('Grid format not recognized, please feed in an acceptable format: csv')

    return grid_df, fixgrid_file_name


def resolve_extras(mesa_extras: dict) -> dict:
    """Enforce the mesa, then posydon, then user order for binary extras.

    Args:
        mesa_extras: The mesa_extras section of the inifile.

    Returns:
        The mesa_extras dictionary with lower-precedence binary extras set to
        None.
    """
    extras_files_types = sorted(set([k.split('_')[0] for k in mesa_extras.keys() if ('binary_extras' in k)]))
    print("WE ARE USING THE EXTRA FILE FROM TYPE {0}".format(extras_files_types[-1]))
    for k in mesa_extras.keys():
        if ('binary_extras' in k) and (extras_files_types[-1] not in k):
            Pwarn("Section mesa_extras value {0} is being set to".format(k) + " None", "ReplaceValueWarning")
            mesa_extras[k] = None

    return mesa_extras
