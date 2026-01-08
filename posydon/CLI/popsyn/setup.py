"""Utility functions for POSYDON population synthesis command line interface."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import os

from posydon.CLI.io import (
    create_bash_submit_script,
    create_python_scripts,
    create_slurm_scripts,
)
from posydon.grids.SN_MODELS import get_SN_MODEL_NAME
from posydon.popsyn.io import binarypop_kwargs_from_ini, simprop_kwargs_from_ini
from posydon.utils.common_functions import convert_metallicity_to_string


def check_SN_MODEL_validity(ini_file, verbose_on_fail=True):
    '''Checks if the step_SN model is valid for this script

    Parameters
    ----------
    ini_file : str
        the path to the ini file
    verbose_on_fail : bool (default: True)
        if `True` rerun the model selection with verbose mode if this failed

    Returns
    -------
    bool
        True if the model is valid or use_interp_values=False, False otherwise
    '''

    simprop_kwargs = simprop_kwargs_from_ini(ini_file)
    step_SN_MODEL = simprop_kwargs['step_SN'][1]
    # always allow the use of non-interpolation values
    if step_SN_MODEL['use_interp_values'] == False:
        return True
    # step_SN MODEL check
    SN_MODEL_NAME_SEL = get_SN_MODEL_NAME(step_SN_MODEL)

    if SN_MODEL_NAME_SEL is None:
        if verbose_on_fail:
            get_SN_MODEL_NAME(step_SN_MODEL, verbose=True)
        return False
    else:
        return True

def validate_ini_file(ini_file):
    '''Validates the ini file for population synthesis

    Parameters
    ----------
    ini_file : str
        the path to the ini file

    Raises
    ------
    FileNotFoundError
        if the ini file does not exist
    ValueError
        if the step_SN MODEL is not valid for this script
    '''
    if not os.path.exists(ini_file):
        raise FileNotFoundError(f'File {ini_file} not found')

    if not check_SN_MODEL_validity(ini_file):
        raise ValueError("The step_SN MODEL is not valid for this script. "
                         "Please check the MODEL in the ini file")

################################################
### Main function to setup population synthesis
################################################

def setup_popsyn_function(args):
    '''Function to setup the population synthesis run

    Parameters
    ----------
    args : argparse.Namespace
        the arguments passed to the function
    '''

    # validate the ini file
    validate_ini_file(args.ini_file)

    synpop_params = binarypop_kwargs_from_ini(args.ini_file)
    metallicities = synpop_params['metallicities']
    if synpop_params['number_of_binaries'] / args.job_array < 1:
        raise ValueError("The number of binaries is less than the job array"
                         " length. Please increase the number of binaries or"
                         " decrease the job array length")

    # create the run_metallicity and merge_metallicity scripts
    create_python_scripts(args.ini_file)

    for met in metallicities:
        # create log directory
        str_met = convert_metallicity_to_string(met)
        os.makedirs(f'{str_met}_logs', exist_ok=True)

        # create SLURM array submission script for this metallicity
        create_slurm_scripts(met, args)

    # create bash submission script for all metallicities at once
    create_bash_submit_script('slurm_submit.sh', metallicities)
