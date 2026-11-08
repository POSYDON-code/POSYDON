"""Module for input/output operations in the CLI of Posydon."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import json
import os
import textwrap

from posydon.CLI.schedulers import (
    ArraySpec,
    JobSpec,
    get_scheduler,
    normalize_scheduler_name,
    parse_option_strings,
)
from posydon.config import PATH_TO_POSYDON, PATH_TO_POSYDON_DATA
from posydon.utils.common_functions import convert_metallicity_to_string
from posydon.utils.posydonwarning import Pwarn

# Output Utility Functions #

# ANSI color codes for terminal output
COLOR_RED = '\033[31m'
COLOR_GREEN = '\033[32m'
COLOR_RESET = '\033[0m'

# Utility scripts for printing colored messages
def print_error(message):
    """Print an error message in red.

    Parameters
    ----------
    message : str
        The error message to print
    """
    print(f"{COLOR_RED}{message}{COLOR_RESET}")

def print_success(message):
    """Print a success message in green.

    Parameters
    ----------
    message : str
        The success message to print
    """
    print(f"{COLOR_GREEN}{message}{COLOR_RESET}")

def print_separator_line():
    """Print a separator line."""
    print("-" * 80)

def print_separator_section():
    """Print a section separator line."""
    print("#" * 80)


def clear_previous_lines(num_lines):
    """Clear the specified number of previous lines in the terminal.

    Parameters
    ----------
    num_lines : int
        Number of lines to clear
    """
    # ANSI escape code to move cursor up one line
    UP = '\033[1A'
    # ANSI escape code to clear line
    CLEAR = '\033[K'

    for _ in range(num_lines):
        # Move up and clear line
        print(UP, end=CLEAR)


# Script Creation Functions #

## Main population setup

def create_merge_script_text(ini_file):
    text = textwrap.dedent(f'''\
        from posydon.popsyn.binarypopulation import BinaryPopulation
        from posydon.popsyn.io import binarypop_kwargs_from_ini
        from posydon.utils.common_functions import convert_metallicity_to_string
        import argparse
        import os

        if __name__ == "__main__":
            parser = argparse.ArgumentParser()
            parser.add_argument("metallicity", type=float)
            args = parser.parse_args()
            ini_kw = binarypop_kwargs_from_ini("{ini_file}")
            ini_kw["metallicity"] = args.metallicity
            str_met = convert_metallicity_to_string(args.metallicity)
            ini_kw["temp_directory"] = str_met + "_Zsun_" + ini_kw["temp_directory"]
            synpop = BinaryPopulation(**ini_kw)
            path_to_batch = ini_kw["temp_directory"]
            tmp_files = [os.path.join(path_to_batch, f) for f in os.listdir(path_to_batch)
                         if os.path.isfile(os.path.join(path_to_batch, f))]
            tmp_files = sorted(tmp_files, key=lambda x: int(x.split(".h5")[0].split(".")[-1]))
            synpop.combine_saved_files(str_met + "_Zsun_population.h5", tmp_files)
            print("done")
            if len(os.listdir(path_to_batch)) == 0:
                os.rmdir(path_to_batch)
    ''')
    return text

def create_run_script_text(ini_file):
    text = textwrap.dedent(f'''\
        from posydon.popsyn.binarypopulation import BinaryPopulation
        from posydon.popsyn.io import binarypop_kwargs_from_ini
        from posydon.utils.common_functions import convert_metallicity_to_string
        from posydon.binary_evol.simulationproperties import SimulationProperties
        import argparse

        if __name__ == "__main__":
            parser = argparse.ArgumentParser()
            parser.add_argument('metallicity', type=float)
            args = parser.parse_args()
            ini_kw = binarypop_kwargs_from_ini('{ini_file}')
            ini_kw['metallicity'] = args.metallicity
            str_met = convert_metallicity_to_string(args.metallicity)
            ini_kw['temp_directory'] = str_met+'_Zsun_' + ini_kw['temp_directory']
            sim_props = SimulationProperties.from_ini('{ini_file}')
            synpop = BinaryPopulation(population_properties=sim_props, **ini_kw)
            synpop.evolve()
    ''')
    return text

def create_run_script(ini_file):
    '''Creates a run script for the population synthesis run.

    This function creates a script that will run the POSYDON population at a specific metallicity.
    The script takes a metallicity as a command-line input.
    It uses that metallicity to run the BinaryPopulation with.

    Parameters
    ----------
    ini_file : str
        Path to the .ini file containing population synthesis parameters.
    '''

    filename ='run_metallicity.py'
    if os.path.exists(filename): # pragma: no cover
        Pwarn('Replace '+filename, "OverwriteWarning")
    with open(filename, mode='w') as file:
        file.write(create_run_script_text(ini_file))

def create_merge_script(ini_file):
    '''Creates a merge script for the population synthesis run.

    The created merge script will combine the "evolution.combined.h5" files
    from the different job arrays in the "temp_directory" into a
    single "{MET}_Zsun_population.h5" file in the run directory.

    Parameters
    ----------
    ini_file : str
        Path to the .ini file containing population synthesis parameters.
    '''

    filename='merge_metallicity.py'
    if os.path.exists(filename): # pragma: no cover
        Pwarn('Replace '+filename, "OverwriteWarning")
    with open(filename, mode='w') as file:
        file.write(create_merge_script_text(ini_file))

## Scheduler independent job descriptions ##

# Settings of a population run, used by `posydon-popsyn check` for rescues
RUN_SETTINGS_FILE = 'popsyn_jobs.json'
RUN_SETTINGS_KEYS = ('job_array', 'walltime', 'merge_walltime', 'mem_per_cpu',
                     'partition', 'account', 'email', 'max_concurrent_jobs',
                     'path_to_posydon', 'path_to_posydon_data')


def array_script_name(str_met, scheduler):
    """File name of the job array script, e.g. "1e+00_Zsun_slurm_array.slurm"."""
    return scheduler.file_name(f"{str_met}_Zsun_{scheduler.name}_array")


def merge_script_name(str_met, scheduler):
    """File name of the merge job script."""
    return scheduler.file_name(f"{str_met}_Zsun_merge_popsyn")


def rescue_script_name(str_met, scheduler):
    """File name of the rescue job script."""
    return scheduler.file_name(f"{str_met}_Zsun_rescue")


def submit_script_name(scheduler):
    """File name of the bash script submitting all jobs of a run."""
    return f"{scheduler.name}_submit.sh"


def resubmit_script_name(scheduler):
    """File name of the bash script submitting the rescue jobs."""
    return f"resubmit_{scheduler.name}.sh"


def popsyn_array_job(metallicity, job_array_length, walltime, mem_per_cpu,
                     partition, account, email, max_concurrent_jobs,
                     path_to_posydon, path_to_posydon_data, task_ids=None):
    '''Job running the population synthesis of one metallicity as an array.

    Parameters
    ----------
    metallicity : float
        The metallicity in solar units
    job_array_length : int
        The length of the job array (number of tasks)
    walltime : str
        Time limit of each task
    mem_per_cpu : str
        Memory per CPU (e.g. "4G")
    partition : str or None
        Queue/partition
    account : str or None
        Account to charge
    email : str or None
        Email address for notifications on failure
    max_concurrent_jobs : int or None
        Maximum number of tasks running at the same time
    path_to_posydon : str
        Path to the POSYDON installation
    path_to_posydon_data : str
        Path to the POSYDON_DATA directory
    task_ids : list of int, optional
        Only run these tasks (rescue job); the array size stays
        `job_array_length`

    Returns
    -------
    JobSpec
    '''
    str_met = convert_metallicity_to_string(metallicity)
    rescue = task_ids is not None
    name = f"{str_met}_popsyn_rescue" if rescue else f"{str_met}_popsyn"
    prefix = "rescue" if rescue else "popsyn"
    return JobSpec(
        name=name, walltime=walltime, mem_per_cpu=mem_per_cpu,
        account=account, queue=partition, email=email, mail_type="FAIL",
        output=f"./{str_met}_logs/{prefix}_job.out",
        setup=[f"export PATH_TO_POSYDON={path_to_posydon}",
               f"export PATH_TO_POSYDON_DATA={path_to_posydon_data}"],
        array=ArraySpec(
            command=f"python ./run_metallicity.py {metallicity}",
            n_tasks=None if rescue else job_array_length,
            task_ids=task_ids, task_count=job_array_length,
            max_concurrent=max_concurrent_jobs,
            log_file=f"./{str_met}_logs/{prefix}_{{job_id}}_{{task_id}}.out"))


def popsyn_merge_job(metallicity, merge_walltime, mem_per_cpu, partition,
                     account, email, path_to_posydon, path_to_posydon_data):
    '''Job merging the array results of one metallicity.

    See `popsyn_array_job` for the parameters.

    Returns
    -------
    JobSpec
    '''
    str_met = convert_metallicity_to_string(metallicity)
    # Override walltime for debug-cpu partition
    if partition == 'debug-cpu':
        merge_walltime = '00:14:00'
    return JobSpec(
        name=f"{str_met}_Zsun_merge", walltime=merge_walltime,
        mem_per_cpu=mem_per_cpu, account=account, queue=partition,
        email=email, mail_type="FAIL",
        output=f"./{str_met}_logs/popsyn_merge.out",
        setup=[f"export PATH_TO_POSYDON={path_to_posydon}",
               f"export PATH_TO_POSYDON_DATA={path_to_posydon_data}"],
        launched_command=f"python ./merge_metallicity.py {metallicity}")


def _valid_option(args, name, types=(str, int, float)):
    '''Get an option from `args` if it is set to a value of `types`.'''
    value = getattr(args, name, None)
    if isinstance(value, bool) or not isinstance(value, types):
        return None
    return value


def scheduler_from_args(args, name=None):
    '''Create the scheduler selected by the command line arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Uses `scheduler`, `scheduler_option` (list of "KEY=VALUE") and
        `exclude`, if present.
    name : str, optional
        Scheduler to use instead of `args.scheduler`.

    Returns
    -------
    Scheduler
    '''
    if name is None:
        name = normalize_scheduler_name(_valid_option(args, 'scheduler'))
    option_strings = getattr(args, 'scheduler_option', None)
    if not isinstance(option_strings, (list, tuple)):
        option_strings = []
    options = parse_option_strings(option_strings)
    exclude = _valid_option(args, 'exclude', (str,))
    if exclude is not None:
        options['exclude'] = exclude
    return get_scheduler(name, **options)


def scheduler_for_run(args, settings=None):
    '''The scheduler of an existing population run.

    The scheduler is `args.scheduler` if given, otherwise the one saved in
    the run's settings, otherwise SLURM (runs set up with older versions).
    Options saved in the settings are used if the scheduler is the same, and
    options given in `args` override them.

    Parameters
    ----------
    args : argparse.Namespace
        Command line arguments (`scheduler`, `scheduler_option`, `exclude`)
    settings : dict, optional
        The run settings (see `read_run_settings`)

    Returns
    -------
    Scheduler
    '''
    settings = settings or {}
    name = normalize_scheduler_name(_valid_option(args, 'scheduler'),
                                    default=None)
    if name is None:
        name = normalize_scheduler_name(settings.get('scheduler'))
    options = {}
    if name == settings.get('scheduler'):
        options.update(settings.get('scheduler_options') or {})
    options.update(scheduler_from_args(args, name).given_options())
    return get_scheduler(name, **options)


## Population Synthesis Script Creation ##

def create_python_scripts(ini_file):
    '''Creates the run and merge scripts for population synthesis.

    Parameters
    ----------
    ini_file : str
        Path to the .ini file containing population synthesis parameters.
    '''
    create_run_script(ini_file)
    create_merge_script(ini_file)
    print("Created run script and merge script")


def create_job_scripts(metallicity, args, scheduler,
                       path_to_posydon=PATH_TO_POSYDON,
                       path_to_posydon_data=None):
    '''Creates the job array and merge scripts for one metallicity.

    Parameters
    ----------
    metallicity : float
        The metallicity in solar units (e.g., 0.02 for Z=0.02)
    args : argparse.Namespace
        the arguments passed to `posydon-popsyn setup`
    scheduler : Scheduler
        The job scheduler backend
    path_to_posydon : str
        Path to the POSYDON installation
    path_to_posydon_data : str, optional
        Path to the POSYDON_DATA directory (default: from the config)

    Returns
    -------
    tuple of str
        The paths of the array and merge scripts
    '''
    if path_to_posydon_data is None:
        path_to_posydon_data = os.path.dirname(PATH_TO_POSYDON_DATA)
    str_met = convert_metallicity_to_string(metallicity)
    array_job = popsyn_array_job(
        metallicity, args.job_array, args.walltime, args.mem_per_cpu,
        args.partition, args.account, args.email, args.max_concurrent_jobs,
        path_to_posydon, path_to_posydon_data)
    merge_job = popsyn_merge_job(
        metallicity, args.merge_walltime, args.mem_per_cpu, args.partition,
        args.account, args.email, path_to_posydon, path_to_posydon_data)
    array_script = scheduler.write(array_job,
                                   array_script_name(str_met, scheduler))
    merge_script = scheduler.write(merge_job,
                                   merge_script_name(str_met, scheduler))
    print(f"{scheduler.name.upper()} script created for metallicity {str_met}")
    return array_script, merge_script


def write_run_settings(args, scheduler, folder='.',
                       path_to_posydon=PATH_TO_POSYDON,
                       path_to_posydon_data=None):
    '''Save the job settings of a population run for later rescues.

    Parameters
    ----------
    args : argparse.Namespace
        the arguments passed to `posydon-popsyn setup`
    scheduler : Scheduler
        The job scheduler backend
    folder : str
        The run folder
    path_to_posydon, path_to_posydon_data : str
        Paths written into the job scripts

    Returns
    -------
    str
        Path to the settings file
    '''
    if path_to_posydon_data is None:
        path_to_posydon_data = os.path.dirname(PATH_TO_POSYDON_DATA)
    settings = {key: _valid_option(args, key) for key in RUN_SETTINGS_KEYS}
    settings['path_to_posydon'] = path_to_posydon
    settings['path_to_posydon_data'] = path_to_posydon_data
    settings['scheduler'] = scheduler.name
    settings['scheduler_options'] = scheduler.given_options()
    path = os.path.join(folder, RUN_SETTINGS_FILE)
    with open(path, 'w') as f:
        json.dump(settings, f, indent=2, sort_keys=True)
    return path


def read_run_settings(run_folder):
    '''Read the job settings saved by `write_run_settings`.

    Returns
    -------
    dict or None
        The settings, or None if the run folder has no settings file (e.g.
        runs set up with older POSYDON versions)
    '''
    path = os.path.join(run_folder, RUN_SETTINGS_FILE)
    if not os.path.isfile(path):
        return None
    with open(path, 'r') as f:
        return json.load(f)


def create_bash_submit_script(filename, metallicities, scheduler=None):
    '''Creates the bash script submitting all jobs of a population run.

    Parameters
    ----------
    filename : str
        The name of the bash submission script to create
    metallicities : list of float
        The list of metallicities in solar units
    scheduler : Scheduler, optional
        The job scheduler backend (default: SLURM)
    '''
    if scheduler is None:
        scheduler = get_scheduler('slurm')
    if os.path.exists(filename): # pragma: no cover
        Pwarn('Replace '+filename, "OverwriteWarning")
    with open(filename, mode='w') as file:
        file.write('#!/bin/bash\n')
        for met in metallicities:
            str_met = convert_metallicity_to_string(met)
            file.write(scheduler.submit_line(
                'array', array_script_name(str_met, scheduler)))
            file.write("echo '"+str_met+" job array submitted as '${array}\n")
            file.write(scheduler.submit_line(
                'merge', merge_script_name(str_met, scheduler),
                after='array'))
            file.write("echo '"+str_met+" merge job submitted as '${merge}\n")

    print_success("Setup complete. You can now submit the jobs using "
                  f"'./{os.path.basename(filename)}'")


## SLURM specific wrappers (kept for backwards compatibility) ##

def create_slurm_array(metallicity,
                       job_array_length,
                       partition,
                       email,
                       walltime,
                       account,
                       mem_per_cpu,
                       max_concurrent_jobs,
                       exclude,
                       path_to_posydon,
                       path_to_posydon_data):
    '''Creates the SLURM array script "{str_met}_Zsun_slurm_array.slurm".

    See `popsyn_array_job` for the parameters; `exclude` is a comma
    separated list of nodes to exclude.
    '''
    scheduler = get_scheduler('slurm', exclude=exclude)
    job = popsyn_array_job(metallicity, job_array_length, walltime,
                           mem_per_cpu, partition, account, email,
                           max_concurrent_jobs, path_to_posydon,
                           path_to_posydon_data)
    str_met = convert_metallicity_to_string(metallicity)
    scheduler.write(job, array_script_name(str_met, scheduler))


def create_slurm_merge(metallicity,
                       partition,
                       email,
                       merge_walltime,
                       account,
                       mem_per_cpu,
                       path_to_posydon,
                       path_to_posydon_data):
    '''Creates the SLURM merge script "{str_met}_Zsun_merge_popsyn.slurm".

    See `popsyn_merge_job` for the parameters.
    '''
    scheduler = get_scheduler('slurm')
    job = popsyn_merge_job(metallicity, merge_walltime, mem_per_cpu,
                           partition, account, email, path_to_posydon,
                           path_to_posydon_data)
    str_met = convert_metallicity_to_string(metallicity)
    scheduler.write(job, merge_script_name(str_met, scheduler))


def create_slurm_rescue(metallicity,
                        missing_indices,
                        job_array_length,
                        partition,
                        email,
                        walltime,
                        account,
                        mem_per_cpu,
                        max_concurrent_jobs,
                        exclude,
                        path_to_posydon,
                        path_to_posydon_data):
    '''Creates the SLURM rescue script "{str_met}_Zsun_rescue.slurm".

    See `popsyn_array_job` for the parameters; only the tasks in
    `missing_indices` are run.
    '''
    scheduler = get_scheduler('slurm', exclude=exclude)
    job = popsyn_array_job(metallicity, job_array_length, walltime,
                           mem_per_cpu, partition, account, email,
                           max_concurrent_jobs, path_to_posydon,
                           path_to_posydon_data, task_ids=missing_indices)
    str_met = convert_metallicity_to_string(metallicity)
    scheduler.write(job, rescue_script_name(str_met, scheduler))


def create_slurm_scripts(metallicity, args):
    '''Creates the SLURM scripts for population synthesis.

    Parameters
    ----------
    metallicity : float
        The metallicity in solar units (e.g., 0.02 for Z=0.02)
    args : argparse.Namespace
        the arguments passed to the function
    '''
    create_job_scripts(metallicity, args, scheduler_from_args(args, 'slurm'))


## Main functions for rescue

def _legacy_slurm_settings(run_folder, str_met):
    """Read the settings of a run without settings file from its SLURM
    array script (runs set up with older POSYDON versions)."""
    slurm_script = os.path.join(run_folder, f'{str_met}_Zsun_slurm_array.slurm')
    with open(slurm_script, 'r') as f:
        lines = f.readlines()

    settings = {key: None for key in RUN_SETTINGS_KEYS}
    settings['scheduler'] = 'slurm'
    settings['scheduler_options'] = {}
    for line in lines:
        if line.startswith('#SBATCH --array='):
            array_range = line.split('=')[1].strip()
            if '%' in array_range:
                tmp_array_range = array_range.split('%')[0]
                settings['max_concurrent_jobs'] = int(array_range.split('%')[1])
            else:
                tmp_array_range = array_range
            if '-' in tmp_array_range:
                start, end = map(int, tmp_array_range.split('-'))
                settings['job_array'] = end - start + 1
        elif line.startswith("#SBATCH --time="):
            settings['walltime'] = line.split('=')[1].strip()
        elif line.startswith("#SBATCH --mem-per-cpu="):
            settings['mem_per_cpu'] = line.split('=')[1].strip()
        elif line.startswith("#SBATCH --partition="):
            settings['partition'] = line.split('=')[1].strip()
        elif line.startswith("#SBATCH --account="):
            settings['account'] = line.split('=')[1].strip()
        elif line.startswith("#SBATCH --mail-user="):
            settings['email'] = line.split('=')[1].strip()
        elif line.startswith("#SBATCH --exclude="):
            settings['scheduler_options']['exclude'] = \
                line.split('=')[1].strip()
        elif line.startswith("export PATH_TO_POSYDON="):
            settings['path_to_posydon'] = line.split('=')[1].strip()
        elif line.startswith("export PATH_TO_POSYDON_DATA="):
            settings['path_to_posydon_data'] = line.split('=')[1].strip()
    return settings


def create_batch_rescue_script(args, batch_status):
    """Generate a script to resubmit failed runs.

    The settings of the original run are read from its settings file
    (`RUN_SETTINGS_FILE`) or, for older runs, from its SLURM array script.
    Arguments given on the command line override them.

    Parameters
    ----------
    args : argparse.Namespace
        The arguments passed to the CLI
    batch_status : dict
        Dictionary with batch status information

    Returns
    -------
    str
        Path to the generated rescue script
    """
    print("Generating rescue script for failed runs...")

    run_folder = args.run_folder
    metallicity = batch_status['metallicity']
    str_met = convert_metallicity_to_string(metallicity)

    # Extract missing indices from batch_status
    missing_indices = batch_status.get('missing_indices', [])
    if not missing_indices:
        raise ValueError("No missing indices found in batch status.\n"
                         "Cannot generate rescue script.\n"
                         "Please run the check function to resubmit the merge jobs.\n\n"
                         "posydon-popsyn check <run_folder>\n")

    settings = read_run_settings(run_folder)
    if settings is None:
        settings = _legacy_slurm_settings(run_folder, str_met)

    # Override with command-line arguments if provided
    for key in ('walltime', 'mem_per_cpu', 'partition', 'account', 'email',
                'max_concurrent_jobs'):
        value = _valid_option(args, key)
        if value is not None:
            settings[key] = value

    scheduler = scheduler_for_run(args, settings)

    job = popsyn_array_job(
        metallicity, settings['job_array'], settings['walltime'],
        settings['mem_per_cpu'], settings['partition'], settings['account'],
        settings['email'], settings['max_concurrent_jobs'],
        settings['path_to_posydon'], settings['path_to_posydon_data'],
        task_ids=sorted(missing_indices))
    rescue_script = os.path.join(run_folder,
                                 rescue_script_name(str_met, scheduler))
    scheduler.write(job, rescue_script)
    print_success(f"Rescue script generated: {rescue_script}")
    return rescue_script


def create_bash_submit_rescue_script(filename, rescue_scripts, scheduler=None):
    """Generate a bash script to submit rescue scripts.

    Parameters
    ----------
    filename : str
        The run folder in which to create the bash submission script
    rescue_scripts : list of str
        The list of rescue script paths
    scheduler : Scheduler, optional
        The job scheduler backend (default: SLURM)

    Returns
    -------
    str
        Path to the generated script
    """
    if scheduler is None:
        scheduler = get_scheduler('slurm')

    # Generate a resubmit script
    print("GENERATING A RESUBMIT SCRIPT ....................\n")
    resubmit_sh_file = os.path.join(filename, resubmit_script_name(scheduler))
    with open(resubmit_sh_file, 'w') as f:
        f.write("#!/bin/bash\n")
        for script in rescue_scripts:
            script_basename = os.path.basename(script)
            met = script_basename.split('_')[0]
            f.write(scheduler.submit_line('resubmit', script_basename))
            f.write("echo 'Rescue script submitted as '${resubmit}\n")
            f.write(scheduler.submit_line(
                'merge', merge_script_name(met, scheduler), after='resubmit'))
            f.write("echo 'Merge job submitted as '${merge}\n")

    print(f"Rescue scripts ready to be resubmitted by 'sh {resubmit_sh_file}'")
    return resubmit_sh_file
