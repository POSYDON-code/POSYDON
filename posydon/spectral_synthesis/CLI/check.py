"""Check function for the spectral synthesis CLI.
"""

__authors__ = [
    "Eirini Kasdagli <kasdaglie@ufl.edu>",
]

import glob
import os

from posydon.spectral_synthesis.CLI.io import (
    create_bash_submit_rescue_script,
    create_slurm_rescue,
    print_error,
    print_separator_line,
    print_success,
)


def get_expected_batch_count(run_folder):
    """Parse SLURM array script to get expected batch count.
    Mirrors get_expected_batch_count in posydon/CLI/popsyn/check.py [5].
    """
    slurm_script = os.path.join(run_folder, 'spectra_slurm_array.slurm')
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


def find_missing_batch_indices(temp_directory, expected_count):
    """Find which batch indices are missing their spectra output.
    Mirrors find_missing_batch_indices in posydon/CLI/popsyn/check.py [5].
    """
    spectra_files = glob.glob(
        os.path.join(temp_directory, 'batch_*_spectra*.h5')
    )
    found_indices = set()
    for f in spectra_files:
        try:
            idx = int(os.path.basename(f).split('_')[1])
            found_indices.add(idx)
        except (ValueError, IndexError):
            pass

    return set(range(expected_count)) - found_indices


def get_user_confirmation(prompt):
    """Get user confirmation.
    Mirrors get_user_confirmation in posydon/CLI/popsyn/check.py [5].
    """
    choice = input(prompt).strip().lower()
    if choice in ['yes', 'y']:
        return True
    elif choice in ['no', 'n']:
        return False
    else:
        print_error(f"Unrecognized input '{choice}'. Treating as 'no'.")
        return False


def check_spectra_function(args):
    """Check the status of a spectral synthesis run.
    Mirrors check_popsyn_function in posydon/CLI/popsyn/check.py [5].

    1. Checks how many batches completed
    2. Reports missing batch IDs
    3. Generates rescue scripts if needed
    """
    run_folder = os.path.abspath(args.run_folder)

    if not os.path.exists(run_folder):
        raise FileNotFoundError(f'Run folder not found: {run_folder}')

    print(f'Checking spectral synthesis run in {run_folder}\n')
    print_separator_line()

    # get expected batch count from SLURM script
    # mirrors: get_expected_batch_count in check.py [5]
    expected_count = args.job_array or get_expected_batch_count(run_folder)
    if expected_count is None:
        print_error(
            'Could not determine expected batch count. '
            'Pass -j to specify manually.'
        )
        return 1

    # find temp_directory from run folder
    temp_directory = os.path.join(run_folder, 'spectra_batches')

    # find completed spectra files
    spectra_files = sorted(glob.glob(
        os.path.join(temp_directory, 'batch_*_spectra*.h5')
    ))
    completed_count = len(spectra_files)

    # find missing batches
    # mirrors: find_missing_batch_indices in check.py [5]
    missing_indices = find_missing_batch_indices(temp_directory, expected_count)
    missing_count = len(missing_indices)

    print(f'Total batches    : {expected_count}')
    print(f'Completed        : {completed_count}')
    print(f'Incomplete/failed: {missing_count}')
    print_separator_line()

    # check if combined output already exists
    output_files = glob.glob(os.path.join(run_folder, '*_spectra.h5'))
    if output_files:
        print_success(f'Combined output found: {output_files[0]}')
        return 0

    if not missing_indices:
        print_success('All batches completed successfully!')
        if get_user_confirmation(
            'Would you like to resubmit the merge job? (yes/no): '
        ):
            os.system('sbatch spectra_merge.slurm')
        return 0

    # mirrors: print missing batch indices in check.py [5]
    if len(missing_indices) <= 10:
        print(f'Failed batch IDs : {sorted(missing_indices)}')
    else:
        sample = sorted(list(missing_indices))[:10]
        print(f'Failed batch IDs include: {sample} and '
              f'{missing_count - 10} more')

    print_separator_line()

    if not get_user_confirmation(
        'Would you like to create rescue scripts for the failed batches? (yes/no): '
    ):
        print('Rescue scripts not created.')
        return 2

    # mirrors: create rescue scripts in check.py [5]
    partition   = args.partition   or 'hpg-default'
    walltime    = args.walltime    or '02:00:00'
    mem_per_cpu = args.mem_per_cpu or '8G'

    create_slurm_rescue(
        missing_indices=missing_indices,
        num_batches=expected_count,
        partition=partition,
        email=args.email,
        walltime=walltime,
        account=args.account,
        mem_per_cpu=mem_per_cpu,
    )

    # mirrors: create_bash_submit_rescue_script in check.py [5]
    resubmit_path = create_bash_submit_rescue_script(
        'resubmit_rescue.sh',
        'spectra_rescue.slurm'
    )

    if get_user_confirmation(
        'Would you like to submit the rescue scripts now? (yes/no): '
    ):
        os.system(f'sh {resubmit_path}')
        print_success('Rescue scripts submitted.')
        return 0
    else:
        print(f'Submit manually with: sh {resubmit_path}')
        return 2
