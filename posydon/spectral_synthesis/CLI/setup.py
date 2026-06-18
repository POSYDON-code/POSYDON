"""Setup function for the spectral synthesis CLI.
Mirrors posydon/CLI/popsyn/setup.py [4].
"""

__authors__ = [
    "Eirini Kasdagli <kasdaglie@ufl.edu>",
]

import os

from posydon.spectral_synthesis.CLI.io import (
    create_bash_submit_script,
    create_python_scripts,
    create_slurm_array,
    create_slurm_merge,
    print_separator_line,
)
from posydon.spectral_synthesis.run_spectral_synthesis import SpectralSynthesisRunner


def setup_spectra_function(args):
    """Setup a spectral synthesis run from an ini file.
    Mirrors setup_popsyn_function in posydon/CLI/popsyn/setup.py [4].
    """
    from posydon.spectral_synthesis.io import (
        spectral_kwargs_from_ini,
        validate_spectral_ini,
    )

    # mirrors: validate_ini_file in setup_popsyn_function [4]
    validate_spectral_ini(args.ini_file)

    # mirrors: binarypop_kwargs_from_ini in setup_popsyn_function [4]
    spectral_kwargs = spectral_kwargs_from_ini(args.ini_file)

    population_file = os.path.abspath(spectral_kwargs['population_file'])
    temp_directory  = os.path.abspath(args.temp_directory)
    num_batches     = args.job_array

    os.makedirs('spectra_logs', exist_ok=True)
    os.makedirs(temp_directory, exist_ok=True)

    print_separator_line()
    print(f'Creating {num_batches} batches from {population_file}...')

    # mirrors: BinaryPopulation init in setup_popsyn_function [4]
    runner = SpectralSynthesisRunner(
        population_file=population_file,
        temp_directory=temp_directory,
        num_batches=num_batches,
        verbose=True,
        **{k: v for k, v in spectral_kwargs.items()
           if k != 'population_file'}
    )
    runner.create_population_batches()

    print_separator_line()

    # mirrors: create_python_scripts in setup_popsyn_function [4]
    create_python_scripts(
        ini_file=args.ini_file,
        temp_directory=temp_directory,
        num_batches=num_batches,
    )

    print_separator_line()

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

    # mirrors: create_bash_submit_script in setup_popsyn_function [4]
    create_bash_submit_script('slurm_submit.sh')
