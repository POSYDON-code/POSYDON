"""Module for checking the status of a population synthesis run.
"""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import glob
import os
import subprocess

import numpy as np

from posydon.CLI.io import (
    clear_previous_lines,
    create_bash_submit_rescue_script,
    create_batch_rescue_script,
    print_error,
    print_separator_line,
    print_success,
)
from posydon.popsyn.io import binarypop_kwargs_from_ini
from posydon.popsyn.synthetic_population import Population
from posydon.utils.common_functions import convert_metallicity_to_string
from posydon.utils.posydonwarning import Pwarn

# File naming patterns
MERGE_SCRIPT_PATTERN = "{met}_Zsun_merge_popsyn.slurm"
ARRAY_SCRIPT_PATTERN = "{met}_Zsun_slurm_array.slurm"
RESCUE_SCRIPT_PATTERN = "{met}_Zsun_rescue.slurm"
RESUBMIT_SCRIPT = "resubmit_slurm.sh"

def get_ini_file(args):
    '''Find and select the INI file for the population synthesis run.

    Parameters
    ----------
    args : argparse.Namespace
        Command line arguments that needs to containing the run_folder path.

    Returns
    -------
    str
        Path to the selected INI file.

    Raises
    ------
    FileNotFoundError
        If no INI file is found in the run folder.
    '''

    # Find and select the INI file
    ini_files = glob.glob(os.path.join(args.run_folder, '*.ini'))
    if not ini_files:
        raise FileNotFoundError("No INI file found in the run folder.")

    # Handle multiple INI files
    if len(ini_files) > 1:
        print("Multiple INI files found:\n")
        for idx, file in enumerate(ini_files):
            print(f"{idx}: {file}")

        print("")
        choice = input("Enter the index of the INI file you want to use: ")
        try:
            selected_index = int(choice)
            if not (0 <= selected_index < len(ini_files)):
                print("Invalid index; using the first INI file.")
                selected_ini = ini_files[0]
            else:
                selected_ini = ini_files[selected_index]
        except ValueError:
            print("Invalid input; using the first INI file.")
            selected_ini = ini_files[0]
    else:
        selected_ini = ini_files[0]

    print(f"\nUsing INI file:\n{selected_ini}")
    return selected_ini

def validate_run_folder(run_folder):
    """Validate that the run folder exists and is not empty.

    Parameters
    ----------
    run_folder : str
        Path to the run folder to validate

    """
    if not os.path.exists(run_folder):
        raise FileNotFoundError(
            f"Run folder '{run_folder}' does not exist.\n"
            "Please provide a valid path to a population run folder."
        )

    try:
        folder_contents = os.listdir(run_folder)
    except (OSError, PermissionError) as e: # pragma: no cover
        raise FileNotFoundError(
            f"Cannot access run folder '{run_folder}': {e}\n"
            "Please check folder permissions."
        )

    if not folder_contents:
        raise ValueError(
            f"Run folder '{run_folder}' is empty.\n"
            "This folder does not contain any population run files."
        )

def get_binary_params(ini_file):
    '''Read the binary population parameters from the INI file

    Parameters
    ----------
    ini_file : str
        The path to the INI file

    Returns
    -------
    int
        The number of metallicities
    int
        The number of binaries
    list of floats
        The list of metallicities in solar units
    dict
        The dictionary of population synthesis parameters from the INI file
    '''
    # Read the population synthesis parameters
    synpop_params = binarypop_kwargs_from_ini(ini_file)
    metallicities = synpop_params.get('metallicities', [])
    num_metallicities = len(metallicities)
    number_of_binaries = synpop_params.get('number_of_binaries', 0)
    return num_metallicities, number_of_binaries, metallicities, synpop_params

def get_run_configuration(args):
    """Get and validate the run configuration from the run folder.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments containing run_folder path

    Returns
    -------
    tuple
        (ini_file, num_metallicities, number_of_binaries, metallicities, synpop_params)

    Raises
    ------
    FileNotFoundError
        If folder validation or INI file retrieval fails
    """
    # Validate folder
    validate_run_folder(args.run_folder)

    # Get INI file
    ini_file = get_ini_file(args)

    # Get binary parameters
    num_metallicities, number_of_binaries, metallicities, synpop_params = \
        get_binary_params(ini_file)

    print("REQUESTED PARAMETERS")
    print(f"# metallicities: {num_metallicities}")
    print(f"# binaries:      {number_of_binaries}\n")

    return ini_file, num_metallicities, number_of_binaries, metallicities, synpop_params

# Check the population files and binary counts

def check_population_files(run_folder, metallicities):
    """Check if all merged population files exist.

    Parameters
    ----------
    run_folder : str
        Path to the folder where the population run is located
    metallicities : list of floats
        List of metallicities to check in solar units

    Returns
    -------
    bool
        True if all files exist, False otherwise
    dict
        Dictionary with metallicity as key and existence status as value
    """
    print("POPULATION FILES CHECK ....................\n")
    print("MET\t\tPOP_FILE")

    all_exist = True
    status_dict = {}

    for met in metallicities:
        str_met = convert_metallicity_to_string(met)
        merged_file = os.path.join(run_folder, f'{str_met}_Zsun_population.h5')

        if not os.path.exists(merged_file):
            print_error(f"{str_met}\t\tNO")
            all_exist = False
            status_dict[met] = False
        else:
            print_success(f"{str_met} \t\tOK")
            status_dict[met] = True

    if not all_exist:
        print_error("POPULATION FILES CHECK ....................ERROR")
    else:
        clear_previous_lines(len(metallicities) + 3)
        print_success("POPULATION FILES CHECK ....................OK")

    return all_exist, status_dict

def check_binary_counts(run_folder, metallicities, expected_count):
    """Check if each population file has the expected number of binaries.

    Parameters
    ----------
    run_folder : str
        Path to the folder where the population run is located
    metallicities : list of floats
        List of metallicities to check in solar units
    expected_count : int
        Expected number of binaries per file

    Returns
    -------
    bool
        True if all counts match, False otherwise
    dict
        Dictionary with metallicity as key and binary count as value
    """
    print("BINARY COUNT CHECK     ....................\n")
    print("MET\t\tEXPECTED\tFOUND\tSTATUS")

    all_match = True
    counts_dict = {}

    for met in metallicities:
        str_met = convert_metallicity_to_string(met)
        merged_file = os.path.join(run_folder, f'{str_met}_Zsun_population.h5')

        try:
            pop = Population(merged_file)
            num_binaries = pop.number_of_systems
            counts_dict[met] = num_binaries

            if num_binaries != expected_count:
                print_error(f"{str_met}\t\t{expected_count}\t\t{num_binaries}"
                             "\tMISMATCH")
                all_match = False
            else:
                print_success(f"{str_met}\t\t{expected_count}\t\t{num_binaries}"
                               "\tOK")
        except Exception as e:
            print_error(f"{str_met}\t\t{expected_count}\t\t-"
                         f"\tERROR: {str(e)}")
            all_match = False

    if not all_match:
        print_error("BINARY COUNT CHECK     ....................ERROR")
    else:
        clear_previous_lines(len(metallicities) + 3)
        print_success("BINARY COUNT CHECK     ....................OK")

    return all_match, counts_dict

def check_run_status(run_folder, metallicities, number_of_binaries):
    """Check whether the population files exist and if they have
    the expected number of binary counts.

    Parameters
    ----------
    run_folder : str
        Path to the run folder
    metallicities : list
        List of metallicities to check
    number_of_binaries : int
        Expected number of binaries

    Returns
    -------
    tuple
        (files_exist, counts_match, file_status)
    """
    print("Checking the status of the population run")
    print_separator_line()

    files_exist, file_status = check_population_files(run_folder, metallicities)

    counts_match = False
    if files_exist:
        counts_match, binary_counts = check_binary_counts(
            run_folder, metallicities, number_of_binaries
        )

    return files_exist, counts_match, file_status


# Check individual batch for reasons of the failed runs

def get_expected_batch_count(run_folder, str_met):
    """Parse SLURM script to get expected batch count.

    Parameters
    ----------
    run_folder : str
        Path to the run folder
    str_met : str
        String representation of metallicity

    Returns
    -------
    int or None
        Expected batch count, or None if not found
    """
    slurm_script = os.path.join(run_folder, ARRAY_SCRIPT_PATTERN.format(met=str_met))

    if not os.path.exists(slurm_script):
        return None

    with open(slurm_script, 'r') as f:
        for line in f:
            if line.startswith('#SBATCH --array='):
                array_range = line.split('=')[1].strip()
                if '-' in array_range:
                    start, end = map(int, array_range.split('-'))
                    return end - start + 1
    return None

def find_missing_batch_indices(batch_folder, expected_count):
    """Find which batch indices are missing.

    Parameters
    ----------
    batch_folder : str
        Path to the batch folder
    expected_count : int
        Expected number of batches

    Returns
    -------
    set
        Set of missing batch indices
    """
    batch_files = glob.glob(os.path.join(batch_folder, 'evolution.combined.*'))
    found_indices = set()
    for batch_file in batch_files:
        file_name = os.path.basename(batch_file)
        if 'evolution.combined' in file_name: # pragma: no cover
            idx_str = file_name.split('.')[-2]
            found_indices.add(int(idx_str))

    return set(range(expected_count)) - found_indices

def print_batch_status(str_met, expected_count, actual_count, missing_indices):
    """Print the batch comparison status.

    Parameters
    ----------
    str_met : str
        String representation of metallicity
    expected_count : int or None
        Expected number of batches
    actual_count : int
        Actual number of batches found
    missing_indices : set
        Set of missing batch indices

    Returns
    -------
    str
        Status string: 'unknown_expected_count', 'incomplete', 'complete', or 'extra_files'
    """
    if expected_count is None:
        print_error(f"{str_met}\t\t?\t{actual_count}\t UNKNOWN EXPECTED COUNT")
        return "unknown_expected_count"

    if actual_count < expected_count:
        missing_count = expected_count - actual_count
        print_error(f"{str_met}\t\t{expected_count}\t{actual_count}"
                    f"\t{missing_count} MISSING")
        # Show missing batch indices
        if len(missing_indices) <= 10:
            print(f"  Missing batches: {sorted(missing_indices)}")
        else:
            sample = sorted(list(missing_indices))[:10]
            print(f"  Missing batches include: {sample} and "
                  f"{len(missing_indices)-10} more")
        return "incomplete"

    elif actual_count == expected_count:
        print_success(f"{str_met}\t\t{expected_count}\t{actual_count}"
                       f"\tCOMPLETE")
        return "complete"

    else:  # actual_count > expected_count
        extra_count = actual_count - expected_count
        print_error(f"{str_met}\t\t{expected_count}\t{actual_count}"
                    f"\t{extra_count} EXTRA")
        return "extra_files"

def select_job_id(run_folder, str_met):
    """Find and select a job ID from available log files.

    Parameters
    ----------
    run_folder : str
        Path to the run folder
    str_met : str
        String representation of metallicity

    Returns
    -------
    int or None
        Selected job ID, or None if no logs found
    """
    jobIDs = glob.glob(os.path.join(run_folder, f'{str_met}_logs/popsyn_*.out'))

    if len(jobIDs) == 0:
        print("\n\033[33mNo log files found. Cannot determine failure reasons.\033[0m")
        print("This may indicate that the jobs were never submitted or the log directory is missing.")
        return None

    # Extract unique job IDs
    jobIDs = [
        int(os.path.basename(job).split('_')[1].split('.')[0])
        for job in jobIDs
    ]
    jobIDs = np.unique(jobIDs)

    # If multiple job IDs, let user select
    if len(jobIDs) > 1:
        print("Please select a job ID to use: ")
        for i, job_id in enumerate(jobIDs):
            print(f"{i}: {job_id}")

        selected_job_idx = None
        while selected_job_idx is None: # pragma: no cover
            try:
                idx = int(input("\nEnter the index to the job ID: "))
                if 0 <= idx < len(jobIDs):
                    selected_job_idx = idx
                    return jobIDs[idx]
                else: # pragma: no cover
                    print_error("Invalid selection. Please try again.")
            except ValueError: # pragma: no cover
                print_error("Please enter a valid number.")
    else:
        return jobIDs[0]

def read_batch_log_file(log_file_path, batch_index, str_met, jobID):
    """Read and analyze a batch log file to determine failure reason.

    Parameters
    ----------
    log_file_path : str
        Path to the log file
    batch_index : int
        Index of the batch being analyzed
    str_met : str
        String representation of the metallicity
    jobID : int
        SLURM job ID

    Returns
    -------
    None
        Prints the analysis results directly
    """
    # Check if the log file exists
    if not os.path.exists(log_file_path):
        print_error(f"Batch {batch_index}: Log file not found")
        return

    try:
        with open(log_file_path, 'r') as f:
            lines = f.readlines()

        if not lines:
            print(f"Batch {batch_index}: <empty file>")
            return

        # Check for time limit exceeded
        if len(lines) >= 3:
            last_3_lines = [line.strip() for line in lines[-3:]]
            if "DUE TO TIME LIMIT" in str(last_3_lines):
                print(f"Batch {batch_index}: Wall time exceeded")
            else:
                print(f"{str_met}_logs/popsyn_{jobID}_{batch_index}.out:")
                for i, line in enumerate(last_3_lines):
                    print(f"  {i+1}: {line}")
        else:
            # File has fewer than 3 lines
            print(f"Batch {batch_index}: File contains {len(lines)} line(s):")
            if "DUE TO TIME LIMIT" in lines[0]:
                print("  Wall time exceeded")
            else:
                print(f"{str_met}_logs/popsyn_{jobID}_{batch_index}.out:")
                for i, line in enumerate(lines):
                    print(f"  {i+1}: {line.strip()}")

    except Exception as e: # pragma: no cover
        print_error(f"Batch {batch_index}: Error reading log file - {str(e)}")

def analyze_missing_batch_logs(run_folder, str_met, missing_indices):
    """Analyze log files for missing batches to determine failure reasons.

    Parameters
    ----------
    run_folder : str
        Path to the run folder
    str_met : str
        String representation of metallicity
    missing_indices : set
        Set of missing batch indices
    """
    if len(missing_indices) == 0:
        return

    jobID = select_job_id(run_folder, str_met)
    if jobID is None:
        return

    print(f"\nREADING THE LOGS OF {jobID} ....")
    print("Reasons for missing files:")

    for index in missing_indices:
        log_file_path = os.path.join(run_folder,
                                    f'{str_met}_logs/popsyn_{jobID}_{index}.out')
        read_batch_log_file(log_file_path, index, str_met, jobID)

def check_batch(run_folder, metallicity, batch_folder_name):
    """Check batch files for a specific metallicity when the population file is missing.

    Parameters
    ----------
    run_folder : str
        Path to the folder where the population run is located
    metallicity : float
        Metallicity to check
    batch_folder_name : str
        Name of the folder containing batch files

    Returns
    -------
    dict
        Dictionary with batch information including:
        - status: str ('complete', 'incomplete', 'folder_missing', etc.)
        - expected_count: int or None
        - found_count: int
        - metallicity: float
        - batch_folder: str
        - missing_indices: set or None
    """
    str_met = convert_metallicity_to_string(metallicity)

    # Get expected batch count from SLURM script
    expected_batch_count = get_expected_batch_count(run_folder, str_met)

    # Check if batch folder exists
    batch_folder = os.path.join(run_folder, f'{str_met}_Zsun_{batch_folder_name}')

    if not os.path.exists(batch_folder):
        print_error(f"{str_met}\t\t-\t-\tNO BATCH DIR")
        return {
            "status": "folder_missing",
            "expected_count": expected_batch_count,
            "found_count": 0,
            "metallicity": metallicity,
            "batch_folder": batch_folder,
            "missing_indices": None,
        }

    # Count actual batch files
    batch_files = glob.glob(os.path.join(batch_folder, 'evolution.combined.*'))
    actual_count = len(batch_files)

    # Find missing batch indices if incomplete
    missing_indices = set()
    if expected_batch_count is not None and actual_count < expected_batch_count:
        missing_indices = find_missing_batch_indices(batch_folder, expected_batch_count)

    # Print status and determine overall status
    status = print_batch_status(str_met, expected_batch_count, actual_count, missing_indices)

    # Analyze logs for missing batches
    analyze_missing_batch_logs(run_folder, str_met, missing_indices)

    return {
        "status": status,
        "expected_count": expected_batch_count,
        "found_count": actual_count,
        "metallicity": metallicity,
        "batch_folder": batch_folder,
        "missing_indices": missing_indices if status == "incomplete" else None,
    }


# Get the status of all batches for missing populations

def get_batches_status(run_folder, missing_files, synpop_params):
    """Get the status of batch files for missing populations.

    Parameters
    ----------
    run_folder : str
        Path to the run folder
    missing_files : dict
        Dictionary of metallicities with missing files
    synpop_params : dict
        Population synthesis parameters

    Returns
    -------
    dict
        Batch status information for each metallicity
    """
    batch_status = {}

    if missing_files:
        print("\nChecking batch files for missing populations:")
        for met in missing_files:
            print_separator_line()
            batch_status[met] = check_batch(
                run_folder,
                met,
                synpop_params.get('temp_directory', 'batches')
            )
            print_separator_line()

    return batch_status

def get_user_confirmation(prompt, valid_yes=None, valid_no=None):
    """Get user confirmation with validation.

    Parameters
    ----------
    prompt : str
        The prompt to show the user
    valid_yes : list, optional
        List of strings considered as "yes" (default: ['yes', 'y'])
    valid_no : list, optional
        List of strings considered as "no" (default: ['no', 'n'])

    Returns
    -------
    bool
        True if user confirmed, False otherwise
    """
    if valid_yes is None:
        valid_yes = ['yes', 'y']
    if valid_no is None:
        valid_no = ['no', 'n']

    choice = input(prompt).strip().lower()

    if choice in valid_yes:
        return True
    elif choice in valid_no:
        return False
    else:
        print_error(f"Unrecognized input '{choice}'. Treating as 'no'.")
        return False

def submit_slurm_job(script_path):
    """Submit a SLURM job and handle errors.

    Parameters
    ----------
    script_path : str
        Path to the SLURM script to submit
    description : str, optional
        Description of the job for error messages

    Returns
    -------
    bool
        True if submission succeeded, False otherwise
    """
    try:
        result = subprocess.run(
            ['sbatch', script_path],
            check=True,
            capture_output=True,
            text=True
        )
        job_id = result.stdout.strip()
        print_success(f"  Job submitted: {job_id}")
        return True
    except subprocess.CalledProcessError as e:
        print_error(f"Failed to submit job: {e.stderr.strip()}")
        return False
    except FileNotFoundError:
        print_error("sbatch command not found. Is SLURM installed?")
        return False

def handle_batches_complete(args, missing_files, batch_status):
    """Handle the case where all batch files are complete.

    Parameters
    ----------
    args : argparse.Namespace
        Command line arguments
    missing_files : dict
        Dictionary of metallicities with missing files
    batch_status : dict
        Batch status information

    Returns
    -------
    bool
        True if merge jobs were handled (resubmitted or user declined), False if batches incomplete
    """
    has_incomplete_batches = any(
        status["status"] != "complete"
        for status in batch_status.values()
    )

    if not has_incomplete_batches:
        print_success("\nAll batch files are complete.")
        print("Please resubmit the merge jobs to generate the population files.")

        if get_user_confirmation("Would you like to resubmit the merge jobs? (yes/no): "):
            all_succeeded = True
            for met in missing_files:
                str_met = convert_metallicity_to_string(met)
                script_path = os.path.join(
                    args.run_folder,
                    MERGE_SCRIPT_PATTERN.format(met=str_met)
                )
                if not submit_slurm_job(script_path): #pragma: no cover
                    all_succeeded = False

            if all_succeeded:
                print_success("All merge jobs submitted successfully.")
            else:
                print_error("Some merge jobs failed to submit. Please check the errors above.") #pragma: no cover
        else:
            print("Merge jobs not resubmitted.")

        return True

    return False

###########################################
## Main function for checking a population
###########################################
def check_popsyn_function(args):
    """Function to check the status of a population run.

    Parameters
    ----------
    args : argparse.Namespace
        The arguments passed to the function
    """
    print(f"Checking the status of the population run in {args.run_folder}\n")

    # Get and validate configuration
    ini_file, num_metallicities, number_of_binaries, metallicities, synpop_params = \
        get_run_configuration(args)

    # Check run status
    files_exist, counts_match, file_status = check_run_status(
        args.run_folder, metallicities, number_of_binaries
    )

    # If all checks passed, we're done
    if files_exist and counts_match:
        print_success("\nAll checks passed successfully.")
        return 0

    # Otherwise, investigate further
    print_error("\nOne or more checks failed.")
    print("We will attempt to rescue the failed runs.\n")
    print_separator_line()

    # Check batch status for missing files
    missing_files = {met: status for met, status in file_status.items() if not status}
    batch_status = get_batches_status(args.run_folder, missing_files, synpop_params)

    # Handle complete batches (just need merge)
    if handle_batches_complete(args, missing_files, batch_status):
        return 2

    # If we get here, some batches are incomplete
    print("\nOne or more batch files are incomplete.\n"
          "We need to resubmit some of the batch jobs.\n")

    # Check if batch folders are missing
    if any(status["status"] == "folder_missing" for status in batch_status.values()):
        print_error("One or more batch folders are missing.")
        print("Cannot generate rescue scripts without the batch folders.")
        print("Please ensure that the batch folders exist and try again.")
        return 1

    # Generate rescue scripts for incomplete batches
    print("We can generate rescue scripts for the incomplete batches.\n",
          "These scripts will resubmit only the missing batch jobs.\n"
          "You can change the script parameters by giving them as arguments\n"
          "to the 'posydon-popsyn rescue' command.\n")

    if get_user_confirmation("Do you want to create rescue scripts for the incomplete batches? "
                          "(yes/no): "):
        # Create rescue scripts
        rescue_scripts = []
        for _, status in batch_status.items():
            print_separator_line()
            print("\nGENERATING A RESCUE SCRIPT ....................")
            rescue_scripts.append(create_batch_rescue_script(args, status))
            print_separator_line()

        resubmit_sh_file = create_bash_submit_rescue_script(args.run_folder, rescue_scripts)

        if get_user_confirmation("Do you want to submit the rescue scripts now? (yes/no): "):
            os.system(f'sh {resubmit_sh_file}')
            print("Rescue scripts submitted.")
            return 0
        else:
            print("Please submit the rescue scripts using the following command:")
            print("\n")
            print(f"sh {RESUBMIT_SCRIPT}")
            print("\n")
            return 2

    else:
        print("Rescue scripts not created. Exiting.")

    return 2
