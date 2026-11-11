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
from posydon.popsyn.io import (
    DEFAULT_INI_PATH,
    binarypop_kwargs_from_ini,
    parse_inifile,
    resolved_inifile_name,
    simprop_kwargs_from_ini,
    write_resolved_inifile,
)
from posydon.utils.common_functions import convert_metallicity_to_string
from posydon.utils.posydonwarning import Pwarn


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
        True if the model is valid, use_interp_values=False, or [step_SN]
        uses a class that is not a StepSN; False otherwise
    '''
    # imported here, since the step classes are only needed for this check
    from posydon.binary_evol.SN.step_SN import StepSN

    SN_class, SN_kwargs = simprop_kwargs_from_ini(ini_file)['step_SN']
    # the SN model check only applies to POSYDON's SN step; a custom class
    # does not use the SN models of the grids
    if not (isinstance(SN_class, type) and issubclass(SN_class, StepSN)):
        return True
    # check the values the step runs with: a section, which swaps in a
    # subclass of StepSN, only lists some keys and the class fills in the rest
    step_SN_MODEL = {**SN_class.DEFAULT_KWARGS, **SN_kwargs}
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

def _step_sections(ini_file):
    '''Return the names of all step sections in the ini file(s)

    Parameters
    ----------
    ini_file : str or list of str
        the path(s) to the ini file(s)

    Returns
    -------
    set of str
        the section names starting with 'step_', with or without an import
    '''
    return {section for section in parse_inifile(ini_file).sections()
            if section.startswith('step_')}

def check_flow_steps(ini_file):
    '''Checks that the flow and the configured steps match

    Every step the flow chart can send a binary to must have a section in the
    ini file (or in the default ini file it is stacked on). Step sections that
    are neither used by the flow nor part of the default ini file raise a
    warning, since they usually point to a typo in a section name.

    Parameters
    ----------
    ini_file : str
        the path to the ini file

    Raises
    ------
    ValueError
        if the ini file has no [flow] section, or if the flow references
        steps that are not configured in the ini file
    '''
    simprop_kwargs = simprop_kwargs_from_ini(ini_file)
    if 'flow' not in simprop_kwargs:
        raise ValueError("The ini file has no [flow] section.")

    flow_function, flow_kwargs = simprop_kwargs['flow']
    flow_steps = set(flow_function(**flow_kwargs).values())
    configured_steps = {key for key in simprop_kwargs
                        if key.startswith('step_')}

    default_steps = set(parse_inifile(DEFAULT_INI_PATH,
                                      use_defaults=False).sections())

    missing_steps = sorted(flow_steps - configured_steps)
    if missing_steps:
        raise ValueError("The flow uses steps that are not configured in the "
                         f"ini file: {missing_steps}. Add a section for each "
                         "of them.")

    # simprop_kwargs only holds sections with an import key; a misspelled
    # section that only sets parameters is found in the ini file itself
    unused_steps = sorted((configured_steps | _step_sections(ini_file))
                          - flow_steps - default_steps)
    if unused_steps:
        Pwarn(f"The ini file configures steps that the flow never uses: "
              f"{unused_steps}. Check the section names for typos.",
              "StepWarning")

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
        if the step_SN MODEL is not valid for this script, or if the flow
        and the configured steps do not match
    '''
    if not os.path.exists(ini_file):
        raise FileNotFoundError(f'File {ini_file} not found')

    check_flow_steps(ini_file)

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

    # write all parameters of the run (defaults + user ini file) to one file;
    # the jobs read this file, so later changes of the defaults or of the
    # user ini file do not change the run
    resolved_ini = os.path.abspath(resolved_inifile_name(args.ini_file))
    write_resolved_inifile(args.ini_file, resolved_ini)
    print(f"Wrote the resolved ini file with all parameters: {resolved_ini}")

    synpop_params = binarypop_kwargs_from_ini(resolved_ini, use_defaults=False)
    # warn if mass ratio q = M2/M1 could fall below 0.05
    if synpop_params['secondary_mass_scheme'] == 'flat_mass_ratio':
        q_min_possible = (synpop_params['secondary_mass_min']
                          / synpop_params['primary_mass_max'])
        if q_min_possible < 0.05:
            Pwarn(
                f"With secondary_mass_min="
                f"{synpop_params['secondary_mass_min']} and "
                f"primary_mass_max="
                f"{synpop_params['primary_mass_max']}, the mass "
                f"ratio q=M2/M1 can be as low as {q_min_possible:.4f}. "
                f"Some binaries with q<0.05 might fall outside the POSYDON "
                f"default grids.",
                "InappropriateValueWarning"
            )

    metallicities = synpop_params['metallicities']
    if synpop_params['number_of_binaries'] / args.job_array < 1:
        raise ValueError("The number of binaries is less than the job array"
                         " length. Please increase the number of binaries or"
                         " decrease the job array length")

    # create the run_metallicity and merge_metallicity scripts
    create_python_scripts(resolved_ini)

    for met in metallicities:
        # create log directory
        str_met = convert_metallicity_to_string(met)
        os.makedirs(f'{str_met}_logs', exist_ok=True)

        # create SLURM array submission script for this metallicity
        create_slurm_scripts(met, args)

    # create bash submission script for all metallicities at once
    create_bash_submit_script('slurm_submit.sh', metallicities)
