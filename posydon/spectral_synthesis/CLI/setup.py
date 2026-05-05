"""Setup function for the spectral synthesis CLI.
Mirrors posydon/CLI/popsyn/setup.py [4].
"""

__authors__ = [
    "Eirini Kasdagli <kasdaglie@ufl.edu>",
]

import os
from posydon.spectral_synthesis.CLI.io import (
    create_bash_submit_script,
    create_merge_script,
    create_run_script,
    create_slurm_array,
    create_slurm_merge,
    print_separator_line,
)
from posydon.spectral_synthesis.run_spectral_synthesis import SpectralSynthesisRunner


def setup_spectra_function(args):
    """Setup a spectral synthesis run.
    Mirrors setup_popsyn_function in posydon/CLI/popsyn/setup.py [4].

    1. Creates batch h5 files from the population
    2. Creates the run and merge python scripts
    3. Creates SLURM array and merge scripts
    4. Creates a bash submit script
    """
    population_file = os.path.abspath(args.population_file)
    temp_directory  = os.path.abspath(args.temp_directory)
    num_batches     = args.job_array

    # mirrors: log directory creation in setup_popsyn_function [4]
    os.makedirs('spectra_logs', exist_ok=True)
    os.makedirs(temp_directory, exist_ok=True)

    print_separator_line()

    # Step 1 — create population batches
    # mirrors: BinaryPopulation init in setup_popsyn_function [4]
    print(f'Creating {num_batches} batches from {population_file}...')
    runner = SpectralSynthesisRunner(
        population_file=population_file,
        temp_directory=temp_directory,
        num_batches=num_batches,
        verbose=True,
    )
    runner.create_population_batches()

    print_separator_line()

    # Step 2 — create python run and merge scripts
    # mirrors: create_python_scripts in setup_popsyn_function [4]
    create_run_script(population_file, temp_directory, num_batches)
    create_merge_script(population_file, temp_directory, num_batches)

    print_separator_line()

    # Step 3 — create SLURM scripts
    # mirrors: create_slurm_scripts in setup_popsyn_function [4]
    create_slurm_array(
        num_batches=num_batches,
        partition=args.partition,
        email=args.email,
        walltime=args.walltime,
        account=args.account,
        mem_per_cpu=args.mem_per_cpu,
    )
    create_slurm_merge(
        partition=args.partition,
        email=args.email,
        merge_walltime=args.merge_walltime,
        account=args.account,
        mem_per_cpu=args.mem_per_cpu,
    )

    print_separator_line()

    # Step 4 — create bash submit script
    # mirrors: create_bash_submit_script in setup_popsyn_function [4]
    create_bash_submit_script('slurm_submit.sh')