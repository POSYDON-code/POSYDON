"""Module for input/output operations in the CLI of the spectral synthesis."""

__authors__ = [
    "Eirini Kasdagli <kasdaglie@ufl.edu>",
]

import os
import textwrap

from posydon.utils.posydonwarning import Pwarn

# mirrors: COLOR codes in posydon/CLI/io.py [3]
COLOR_RED   = '\033[31m'
COLOR_GREEN = '\033[32m'
COLOR_RESET = '\033[0m'


# ──────────────────────────────────────────────
# Output Utility Functions
# mirrors: print_error, print_success in posydon/CLI/io.py [3]
# ──────────────────────────────────────────────

def print_error(message):
    """Print an error message in red.
    Mirrors print_error in posydon/CLI/io.py [3].
    """
    print(f"{COLOR_RED}{message}{COLOR_RESET}")


def print_success(message):
    """Print a success message in green.
    Mirrors print_success in posydon/CLI/io.py [3].
    """
    print(f"{COLOR_GREEN}{message}{COLOR_RESET}")


def print_separator_line():
    """Print a separator line.
    Mirrors print_separator_line in posydon/CLI/io.py [3].
    """
    print("-" * 80)


# ──────────────────────────────────────────────
# Script Text Creation
# mirrors: create_run_script_text and create_merge_script_text in io.py [3]
# ──────────────────────────────────────────────

def create_run_script_text(ini_file, temp_directory, num_batches):
    """Create the text for the run script.
    Reads all kwargs directly from the ini file at runtime.
    Mirrors create_run_script_text in posydon/CLI/io.py [3].
    """
    text = textwrap.dedent(f'''\
        import os
        from posydon.spectral_synthesis.run_spectral_synthesis import (
            SpectralSynthesisRunner
        )
        from posydon.spectral_synthesis.io import spectral_kwargs_from_ini

        if __name__ == "__main__":
            # mirrors: ini_kw = binarypop_kwargs_from_ini in run_metallicity.py [3]
            spectral_kwargs = spectral_kwargs_from_ini("{ini_file}")

            # mirrors: JOB_ID and RANK from SLURM in BinaryPopulation [2]
            JOB_ID = int(os.environ["SLURM_ARRAY_TASK_ID"])
            RANK   = int(os.environ.get("SLURM_PROCID", 0))
            size   = int(os.environ.get("SLURM_NTASKS", 1))

            runner = SpectralSynthesisRunner(
                population_file=spectral_kwargs.pop("population_file"),
                temp_directory="{temp_directory}",
                num_batches={num_batches},
                verbose=True,
                JOB_ID=JOB_ID,
                RANK=RANK,
                size=size,
                **spectral_kwargs
            )
            runner.make_spectra()
        ''')
    return text


def create_merge_script_text(ini_file, temp_directory, num_batches):
    """Create the text for the merge script.
    Mirrors create_merge_script_text in posydon/CLI/io.py [3].
    """
    text = textwrap.dedent(f'''\
        from posydon.spectral_synthesis.run_spectral_synthesis import (
            SpectralSynthesisRunner
        )
        from posydon.spectral_synthesis.io import spectral_kwargs_from_ini

        if __name__ == "__main__":
            # mirrors: ini_kw = binarypop_kwargs_from_ini in merge_metallicity.py [3]
            spectral_kwargs = spectral_kwargs_from_ini("{ini_file}")

            runner = SpectralSynthesisRunner(
                population_file=spectral_kwargs.pop("population_file"),
                temp_directory="{temp_directory}",
                num_batches={num_batches},
                verbose=True,
            )
            runner.combine_batches()
        ''')
    return text


# ──────────────────────────────────────────────
# Script File Creation
# mirrors: create_run_script, create_merge_script in posydon/CLI/io.py [3]
# ──────────────────────────────────────────────

def create_run_script(ini_file, temp_directory, num_batches):
    """Create the run script file for the spectral synthesis run.
    Mirrors create_run_script in posydon/CLI/io.py [3].
    """
    filename = 'run_spectral_batch.py'
    if os.path.exists(filename):
        Pwarn(f'Replace {filename}', "OverwriteWarning")
    with open(filename, mode='w') as f:
        f.write(create_run_script_text(ini_file, temp_directory, num_batches))
    print(f'Created run script -> {filename}')


def create_merge_script(ini_file, temp_directory, num_batches):
    """Create the merge script file for the spectral synthesis run.
    Mirrors create_merge_script in posydon/CLI/io.py [3].
    """
    filename = 'merge_spectral_batches.py'
    if os.path.exists(filename):
        Pwarn(f'Replace {filename}', "OverwriteWarning")
    with open(filename, mode='w') as f:
        f.write(create_merge_script_text(ini_file, temp_directory, num_batches))
    print(f'Created merge script -> {filename}')


def create_python_scripts(ini_file, temp_directory, num_batches):
    """Create run and merge scripts for the spectral synthesis run.
    Mirrors create_python_scripts in posydon/CLI/io.py [3].
    """
    create_run_script(ini_file, temp_directory, num_batches)
    create_merge_script(ini_file, temp_directory, num_batches)
    print('Created run script and merge script')


# ──────────────────────────────────────────────
# SLURM Script Creation
# mirrors: create_slurm_array, create_slurm_merge in posydon/CLI/io.py [3]
# ──────────────────────────────────────────────

def create_slurm_array(num_batches, partition, email,
                       walltime, account, mem_per_cpu):
    """Create the SLURM array script for the spectral synthesis run.
    Mirrors create_slurm_array in posydon/CLI/io.py [3].
    """
    # mirrors: job_array_length = job_array_length - 1 in io.py [3]
    job_array_length = num_batches - 1

    # mirrors: optional_directives pattern in posydon/CLI/io.py [3]
    optional_directives = []
    if account is not None:
        optional_directives.append(f'#SBATCH --account={account}')
    if partition is not None:
        optional_directives.append(f'#SBATCH --partition={partition}')
    if email is not None:
        optional_directives.extend([
            '#SBATCH --mail-type=FAIL',
            f'#SBATCH --mail-user={email}'
        ])

    optional_section = '\n'.join(optional_directives)
    if optional_section:
        optional_section += '\n'

    text_pre = textwrap.dedent(f'''\
        #!/bin/bash
        #SBATCH --array=0-{job_array_length}
        #SBATCH --job-name=spectra_batch
        #SBATCH --output=./spectra_logs/batch_%A_%a.out
        #SBATCH --time={walltime}
        #SBATCH --mem-per-cpu={mem_per_cpu}
        ''')

    text_post = textwrap.dedent('''\
        srun python ./run_spectral_batch.py
        ''')

    text = text_pre + optional_section + text_post

    filename = 'spectra_slurm_array.slurm'
    if os.path.exists(filename):
        Pwarn(f'Replace {filename}', "OverwriteWarning")
    with open(filename, mode='w') as f:
        f.write(text)
    print(f'Created SLURM array script -> {filename}')


def create_slurm_merge(partition, email, merge_walltime, account, mem_per_cpu):
    """Create the SLURM merge script for the spectral synthesis run.
    Mirrors create_slurm_merge in posydon/CLI/io.py [3].
    """
    # mirrors: optional_directives pattern in posydon/CLI/io.py [3]
    optional_directives = []
    if account is not None:
        optional_directives.append(f'#SBATCH --account={account}')
    if partition is not None:
        optional_directives.append(f'#SBATCH --partition={partition}')
    if email is not None:
        optional_directives.extend([
            '#SBATCH --mail-type=FAIL',
            f'#SBATCH --mail-user={email}'
        ])

    optional_section = '\n'.join(optional_directives)
    if optional_section:
        optional_section += '\n'

    text_pre = textwrap.dedent(f'''\
        #!/bin/bash
        #SBATCH --job-name=spectra_merge
        #SBATCH --output=./spectra_logs/merge.out
        #SBATCH --mem-per-cpu={mem_per_cpu}
        #SBATCH --time={merge_walltime}
        ''')

    text_post = textwrap.dedent('''\
        srun python ./merge_spectral_batches.py
        ''')

    text = text_pre + optional_section + text_post

    filename = 'spectra_merge.slurm'
    if os.path.exists(filename):
        Pwarn(f'Replace {filename}', "OverwriteWarning")
    with open(filename, mode='w') as f:
        f.write(text)
    print(f'Created SLURM merge script -> {filename}')


def create_slurm_rescue(missing_indices, num_batches, partition,
                        email, walltime, account, mem_per_cpu):
    """Create the SLURM rescue script for failed spectral synthesis jobs.
    Mirrors create_slurm_rescue in posydon/CLI/io.py [3].
    """
    # mirrors: job_array_str in create_slurm_rescue [3]
    job_array_str = ','.join(map(str, sorted(missing_indices)))

    optional_directives = []
    if account is not None:
        optional_directives.append(f'#SBATCH --account={account}')
    if partition is not None:
        optional_directives.append(f'#SBATCH --partition={partition}')
    if email is not None:
        optional_directives.extend([
            '#SBATCH --mail-type=FAIL',
            f'#SBATCH --mail-user={email}'
        ])

    optional_section = '\n'.join(optional_directives)
    if optional_section:
        optional_section += '\n'

    text_pre = textwrap.dedent(f'''\
        #!/bin/bash
        #SBATCH --array={job_array_str}
        #SBATCH --job-name=spectra_rescue
        #SBATCH --output=./spectra_logs/rescue_%A_%a.out
        #SBATCH --time={walltime}
        #SBATCH --mem-per-cpu={mem_per_cpu}
        ''')

    text_post = textwrap.dedent(f'''\
        export SLURM_ARRAY_TASK_COUNT={num_batches}
        export SLURM_ARRAY_TASK_MIN=0

        srun python ./run_spectral_batch.py
        ''')

    text = text_pre + optional_section + text_post

    filename = 'spectra_rescue.slurm'
    if os.path.exists(filename):
        Pwarn(f'Replace {filename}', "OverwriteWarning")
    with open(filename, mode='w') as f:
        f.write(text)
    print(f'Created rescue script -> {filename}')


# ──────────────────────────────────────────────
# Bash Submission Scripts
# mirrors: create_bash_submit_script in posydon/CLI/io.py [3]
# ──────────────────────────────────────────────

def create_bash_submit_script(filename):
    """Create the bash submission script.
    Mirrors create_bash_submit_script in posydon/CLI/io.py [3].
    """
    if os.path.exists(filename):
        Pwarn(f'Replace {filename}', "OverwriteWarning")
    with open(filename, mode='w') as f:
        f.write('#!/bin/bash\n')
        f.write('array=$(sbatch --parsable spectra_slurm_array.slurm)\n')
        f.write("echo 'Batch job array submitted as '${array}\n")
        f.write('merge=$(sbatch --parsable --dependency=afterok:${array} '
                '--kill-on-invalid-dep=yes spectra_merge.slurm)\n')
        f.write("echo 'Merge job submitted as '${merge}\n")

    os.chmod(filename, 0o755)
    print_success(
        f"Setup complete. You can now submit the jobs using './{filename}'"
    )


def create_bash_submit_rescue_script(filename, rescue_script):
    """Create the bash rescue submission script.
    Mirrors create_bash_submit_rescue_script in posydon/CLI/io.py [3].
    """
    with open(filename, mode='w') as f:
        f.write('#!/bin/bash\n')
        f.write(f'resubmit=$(sbatch --parsable {rescue_script})\n')
        f.write("echo 'Rescue script submitted as '${resubmit}\n")
        f.write('merge=$(sbatch --parsable --dependency=afterok:${resubmit} '
                '--kill-on-invalid-dep=yes spectra_merge.slurm)\n')
        f.write("echo 'Merge job submitted as '${merge}\n")

    os.chmod(filename, 0o755)
    print(f'Rescue scripts ready to be resubmitted by sh {filename}')
    return filename
