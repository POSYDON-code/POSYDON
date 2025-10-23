"""Module for input/output operations in the CLI of Posydon."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import os
import textwrap

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

def create_slurm_array(metallicity,
                        job_array_length,
                        partition,
                        email,
                        walltime,
                        account,
                        mem_per_cpu,
                        path_to_posydon,
                        path_to_posydon_data):
    '''Creates the slurm array script for population synthesis job arrays.

    Creates a SLURM submission script file named
    "{str_met}_Zsun_slurm_array.slurm" where str_met is the metallicity
    converted to string format.

    Parameters
    ----------
    metallicity : float
        The metallicity in solar units (e.g., 0.02 for Z=0.02)
    job_array_length : int
        The length of the job array (number of jobs)
    partition : str, optional
        SLURM partition to submit the job to
    email : str, optional
        Email address for job notifications
    walltime : str
        Time limit for the job in SLURM format (e.g., "24:00:00")
    account : str, optional
        SLURM account to charge the job to
    mem_per_cpu : str
        Memory per CPU allocation (e.g., "4G")
    path_to_posydon : str
        Path to the POSYDON installation
    path_to_posydon_data : str
        Path to the POSYDON_DATA directory
    '''
    str_met = convert_metallicity_to_string(metallicity)
    job_array_length = job_array_length - 1  # 0 already included

    # Build optional SBATCH directives
    optional_directives = []
    if account is not None:
        optional_directives.append(f"#SBATCH --account={account}")
    if partition is not None:
        optional_directives.append(f"#SBATCH --partition={partition}")
    if email is not None:
        optional_directives.extend([
            "#SBATCH --mail-type=FAIL",
            f"#SBATCH --mail-user={email}"
        ])

    optional_section = "\n".join(optional_directives)
    if optional_section:
        optional_section += "\n"

    text = textwrap.dedent(f'''\
        #!/bin/bash
        #SBATCH --array=0-{job_array_length}
        #SBATCH --job-name={str_met}_popsyn
        #SBATCH --output=./{str_met}_logs/popsyn_%A_%a.out
        #SBATCH --time={walltime}
        #SBATCH --mem-per-cpu={mem_per_cpu}
        {optional_section}
        export PATH_TO_POSYDON={path_to_posydon}
        export PATH_TO_POSYDON_DATA={path_to_posydon_data}
        srun python ./run_metallicity.py {metallicity}
    ''')

    filename = f"{str_met}_Zsun_slurm_array.slurm"
    if os.path.exists(filename): # pragma: no cover
        Pwarn('Replace ' + filename, "OverwriteWarning")
    with open(filename, mode='w') as file:
        file.write(text)

def create_slurm_merge(metallicity,
                       partition,
                       email,
                       merge_walltime,
                       account,
                       mem_per_cpu,
                       path_to_posydon,
                       path_to_posydon_data):
    '''Creates the slurm submit script for merging population synthesis results.

    Creates a SLURM submission script file named
    "{str_met}_Zsun_merge_popsyn.slurm" where str_met is the metallicity
    converted to string format.

    Parameters
    ----------
    metallicity : float
        The metallicity in solar units (e.g., 0.02 for Z=0.02)
    partition : str, optional
        SLURM partition to submit the job to
    email : str, optional
        Email address for job notifications
    merge_walltime : str
        Time limit for the merge job in SLURM format (e.g., "12:00:00")
    account : str, optional
        SLURM account to charge the job to
    mem_per_cpu : str
        Memory per CPU allocation (e.g., "4G")
    path_to_posydon : str
        Path to the POSYDON installation
    path_to_posydon_data : str
        Path to the POSYDON_DATA directory
    '''
    str_met = convert_metallicity_to_string(metallicity)

    # Override walltime for debug-cpu partition
    if partition == 'debug-cpu':
        merge_walltime = '00:14:00'

    # Build optional SBATCH directives
    optional_directives = []
    if account is not None:
        optional_directives.append(f"#SBATCH --account={account}")
    if partition is not None:
        optional_directives.append(f"#SBATCH --partition={partition}")
    if email is not None:
        optional_directives.extend([
            "#SBATCH --mail-type=FAIL",
            f"#SBATCH --mail-user={email}"
        ])

    optional_section = "\n".join(optional_directives)
    if optional_section:
        optional_section += "\n"

    text = textwrap.dedent(f'''\
        #!/bin/bash
        #SBATCH --job-name={str_met}_Zsun_merge
        #SBATCH --output=./{str_met}_logs/popsyn_merge.out
        #SBATCH --mem-per-cpu={mem_per_cpu}
        #SBATCH --time={merge_walltime}
        {optional_section}
        export PATH_TO_POSYDON={path_to_posydon}
        export PATH_TO_POSYDON_DATA={path_to_posydon_data}
        srun python ./merge_metallicity.py {metallicity}
    ''')

    filename = f'{str_met}_Zsun_merge_popsyn.slurm'
    if os.path.exists(filename): # pragma: no cover
        Pwarn('Replace ' + filename, "OverwriteWarning")
    with open(filename, mode='w') as file:
        file.write(text)

## Rescue Script Creation Functions

def create_slurm_rescue(metallicity,
                        missing_indices,
                        job_array_length,
                        partition,
                        email,
                        walltime,
                        account,
                        mem_per_cpu,
                        path_to_posydon,
                        path_to_posydon_data):
    '''Creates the slurm rescue script for resubmitting failed population synthesis jobs.

    Creates a SLURM submission script file named
    "{str_met}_Zsun_rescue.slurm" where str_met is the metallicity
    converted to string format.

    Parameters
    ----------
    metallicity : float
        The metallicity in solar units (e.g., 0.02 for Z=0.02)
    missing_indices : list of int
        List of failed job indices to resubmit
    job_array_length : int
        The original length of the job array (number of jobs)
    partition : str, optional
        SLURM partition to submit the job to
    email : str, optional
        Email address for job notifications
    walltime : str
        Time limit for the job in SLURM format (e.g., "24:00:00")
    account : str, optional
        SLURM account to charge the job to
    mem_per_cpu : str
        Memory per CPU allocation (e.g., "4G")
    path_to_posydon : str
        Path to the POSYDON installation
    path_to_posydon_data : str
        Path to the POSYDON_DATA directory
    '''
    str_met = convert_metallicity_to_string(metallicity)

    # Format job array string for SBATCH
    job_array_str = ','.join(map(str, sorted(missing_indices)))

    # Build optional SBATCH directives
    optional_directives = []
    if account is not None:
        optional_directives.append(f"#SBATCH --account={account}")
    if partition is not None:
        optional_directives.append(f"#SBATCH --partition={partition}")
    if email is not None:
        optional_directives.extend([
            "#SBATCH --mail-type=FAIL",
            f"#SBATCH --mail-user={email}"
        ])

    optional_section = "\n".join(optional_directives)
    if optional_section:
        optional_section += "\n"

    text = textwrap.dedent(f'''\
        #!/bin/bash
        #SBATCH --array={job_array_str}
        #SBATCH --job-name={str_met}_popsyn_rescue
        #SBATCH --output=./{str_met}_logs/rescue_%A_%a.out
        #SBATCH --time={walltime}
        #SBATCH --mem-per-cpu={mem_per_cpu}
        {optional_section}
        export PATH_TO_POSYDON={path_to_posydon}
        export PATH_TO_POSYDON_DATA={path_to_posydon_data}

        export SLURM_ARRAY_TASK_COUNT={job_array_length}
        export SLURM_ARRAY_TASK_MIN=0

        srun python ./run_metallicity.py {metallicity}
    ''')

    filename = f'{str_met}_Zsun_rescue.slurm'
    if os.path.exists(filename): # pragma: no cover
        Pwarn('Replace ' + filename, "OverwriteWarning")
    with open(filename, mode='w') as file:
        file.write(text)

## Main Functions for Population Synthesis Setup ##

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

def create_slurm_scripts(metallicity, args): # pragma: no cover
    '''Creates the slurm scripts for population synthesis.

    Parameters
    ----------
    metallicity : float
        The metallicity in solar units (e.g., 0.02 for Z=0.02)
    args : argparse.Namespace
        the arguments passed to the function
    '''
    create_slurm_array(metallicity, args.job_array, args.partition, args.email,
                       args.walltime, args.account, args.mem_per_cpu,
                       PATH_TO_POSYDON,
                       os.path.dirname(PATH_TO_POSYDON_DATA))

    create_slurm_merge(metallicity, args.partition, args.email,
                       args.merge_walltime, args.account, args.mem_per_cpu,
                       PATH_TO_POSYDON,
                       os.path.dirname(PATH_TO_POSYDON_DATA))

    print(f"SLURM script created for metallicity {convert_metallicity_to_string(metallicity)}")

def create_bash_submit_script(filename, metallicities):
    '''Creates the bash submission script for all SLURM jobs.

    Parameters
    ----------
    filename : str
        The name of the bash submission script to create
    metallicities : list of float
        The list of metallicities in solar units
    '''

    if os.path.exists(filename): # pragma: no cover
        Pwarn('Replace '+filename, "OverwriteWarning")
    with open(filename, mode='w') as file:
        file.write('#!/bin/bash\n')
        for met in metallicities:
            str_met = convert_metallicity_to_string(met)
            file.write("array=$(sbatch --parsable "
                       f"{str_met}_Zsun_slurm_array.slurm)\n")
            file.write("echo '"+str_met+" job array submitted as '${array}\n")
            file.write("merge=$(sbatch --parsable --dependency=afterok:${array} "
                       "--kill-on-invalid-dep=yes "
                       f'{str_met}_Zsun_merge_popsyn.slurm)\n')
            file.write("echo '"+str_met+" merge job submitted as '${merge}\n")

    print_success("Setup complete. You can now submit the jobs using './slurm_submit.sh'")

## Main functions for rescue
def create_batch_rescue_script(args, batch_status):
    """Generate a script to resubmit failed runs.

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
    slurm_script = os.path.join(run_folder, f'{str_met}_Zsun_slurm_array.slurm')

    # Extract missing indices from batch_status
    missing_indices = batch_status.get('missing_indices', [])
    if not missing_indices:
        raise ValueError("No missing indices found in batch status.\n"
                         "Cannot generate rescue script.\n"
                         "Please run the check function to resubmit the merge jobs.\n\n"
                         "posydon-popsyn check <run_folder>\n")

    # Parse the current slurm array script
    with open(slurm_script, 'r') as f:
        lines = f.readlines()

    # Extract parameters from the current script
    job_array_length = None
    account = None
    partition = None
    email = None
    walltime = None
    mem_per_cpu = None
    path_to_posydon = None
    path_to_posydon_data = None

    for line in lines:
        if line.startswith('#SBATCH --array='):
            array_range = line.split('=')[1].strip()
            if '-' in array_range:
                start, end = map(int, array_range.split('-'))
                job_array_length = end - start + 1
        elif line.startswith("#SBATCH --time="):
            walltime = line.split('=')[1].strip()
        elif line.startswith("#SBATCH --mem-per-cpu="):
            mem_per_cpu = line.split('=')[1].strip()
        elif line.startswith("#SBATCH --partition="):
            partition = line.split('=')[1].strip()
        elif line.startswith("#SBATCH --account="):
            account = line.split('=')[1].strip()
        elif line.startswith("#SBATCH --mail-user="):
            email = line.split('=')[1].strip()
        elif line.startswith("export PATH_TO_POSYDON="):
            path_to_posydon = line.split('=')[1].strip()
        elif line.startswith("export PATH_TO_POSYDON_DATA="):
            path_to_posydon_data = line.split('=')[1].strip()

    # Override with command-line arguments if provided
    if args.walltime is not None:
        walltime = args.walltime
    if args.mem_per_cpu is not None:
        mem_per_cpu = args.mem_per_cpu
    if args.partition is not None:
        partition = args.partition
    if args.account is not None:
        account = args.account
    if args.email is not None:
        email = args.email

    # Create the rescue script
    create_slurm_rescue(
        metallicity=metallicity,
        missing_indices=missing_indices,
        job_array_length=job_array_length,
        partition=partition,
        email=email,
        walltime=walltime,
        account=account,
        mem_per_cpu=mem_per_cpu,
        path_to_posydon=path_to_posydon,
        path_to_posydon_data=path_to_posydon_data
    )

    rescue_script = os.path.join(run_folder, f'{str_met}_Zsun_rescue.slurm')
    print_success(f"Rescue script generated: {rescue_script}")
    return rescue_script

def create_bash_submit_rescue_script(filename, rescue_scripts):
    """Generate a bash script to submit rescue scripts.

    Parameters
    ----------
    filename : str
        The name of the bash submission script to create
    rescue_scripts : list of str
        The list of rescue script paths
    """
    MERGE_SCRIPT_PATTERN = "{met}_Zsun_merge_popsyn.slurm"
    RESUBMIT_SCRIPT = "resubmit_slurm.sh"

    # Generate a resubmit_slurm.sh script
    print("GENERATING A RESUBMIT SCRIPT ....................\n")
    resubmit_sh_file = os.path.join(filename, RESUBMIT_SCRIPT)
    with open(resubmit_sh_file, 'w') as f:
        f.write("#!/bin/bash\n")
        for script in rescue_scripts:
            script_basename = os.path.basename(script)
            met = script_basename.split('_')[0]
            merge_script = MERGE_SCRIPT_PATTERN.format(met=met)
            f.write(f'resubmit=$(sbatch --parsable {script_basename})\n')
            f.write("echo 'Rescue script submitted as '${resubmit}\n")
            f.write("merge=$(sbatch --parsable "
                    "--dependency=afterok:${resubmit} "
                    "--kill-on-invalid-dep=yes "
                    f"{merge_script})\n")
            f.write("echo 'Merge job submitted as '${merge}\n")

    print(f"Rescue scripts ready to be resubmitted by 'sh {resubmit_sh_file}'")
    return resubmit_sh_file
