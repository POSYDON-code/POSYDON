"""Setup orchestrator for the posydon grid setup."""

import argparse
import os
import shutil
import subprocess

from posydon.active_learning.psy_cris.utils import parse_inifile
from posydon.CLI.grids.config import (
    normalize_defaults,
    parse_config_file,
    read_grid_file,
    resolve_extras,
    validate_config,
)
from posydon.CLI.grids.executables import make_executables
from posydon.CLI.grids.inlists import construct_static_inlist
from posydon.CLI.grids.log import setup_logger
from posydon.CLI.grids.scenario import setup_inlists_from_scenario
from posydon.CLI.grids.submission import (
    construct_command_line,
    write_shell_submission_script,
    write_slurm_submission_scripts,
)


def parse_commandline() -> argparse.Namespace:
    """Parse the arguments given on the command-line.

    Returns:
        The parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inifile",
                        help="Name of ini file of params",
                        required=True)
    parser.add_argument("--grid-type",
                        help="Either you are supplying a grid "
                        "of points to run MESA on (fixed) or you are supplying a pre computed MESA "
                        "grid and want to sample new points to run MESA on (dynamic).",
                        required=True)
    parser.add_argument("--run-directory",
                        help="Path where executable will be made and MESA "
                             "simulation output will be placed", default=os.getcwd())
    parser.add_argument("--submission-type",
                        help="Options include creating a shell script or a slurm script",
                        default='shell')
    parser.add_argument("-n", "--nproc",
                        help="number of processors", type=int, default=1)
    parser.add_argument("--verbose", action="store_true", default=False,
                        help="Run in Verbose Mode")

    args = parser.parse_args()

    if args.grid_type not in ['fixed', 'dynamic']:
        raise parser.error("--grid-type must be either fixed or dynamic")

    if args.submission_type not in ['slurm', 'shell']:
        raise parser.error('--submission-type must be either slurm of shell')

    return args


def validate_input(args: argparse.Namespace) -> None:
    """Validate that the environment is ready to run a grid.

    Args:
        args: The parsed command-line arguments.

    Raises:
        ValueError: If MESA_DIR is not defined in the environment.
    """
    try:
        os.environ['MESA_DIR']
    except KeyError:
        raise ValueError("MESA_DIR must be defined in your environment "
                          "before you can run a grid os MESA runs")


def find_run_grid_executable() -> str:
    """Locate the posydon-run-grid executable in the current PATH.

    Returns:
        The decoded path to the posydon-run-grid executable.

    Raises:
        ValueError: If the posydon-run-grid executable cannot be located.
    """
    proc = subprocess.Popen(['which', 'posydon-run-grid'],
                                 stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE
                             )
    (path_to_run_grid_exec, err) = proc.communicate()
    if not path_to_run_grid_exec:
        raise ValueError('Cannot locate posydon-run-grid executable in your path')
    else:
        return path_to_run_grid_exec.decode('utf-8').strip('\n')


def run_setup(args: argparse.Namespace) -> None:
    """Run the full grid setup flow based on the parsed arguments.

    Args:
        args: The parsed command-line arguments.
    """
    setup_logger(args.verbose)
    validate_input(args)
    path_to_run_grid_exec = find_run_grid_executable()

    run_parameters, slurm, mesa_inlists, mesa_extras = parse_config_file(args.inifile)
    normalize_defaults(run_parameters, mesa_inlists)
    validate_config(run_parameters, args.grid_type)

    if mesa_inlists['scenario'] is not None:
        setup_inlists_from_scenario(source=mesa_inlists['scenario'][0],
                                    gitcommit=mesa_inlists['scenario'][1],
                                    system_type=mesa_inlists['scenario'][2],
                                    mesa_inlists=mesa_inlists,
                                    mesa_extras=mesa_extras)

    grid_df, fixgrid_file_name = read_grid_file(run_parameters['grid'])

    mesa_extras = resolve_extras(mesa_extras)

    binary_exe, star1_exe, star2_exe = make_executables(mesa_extras=mesa_extras,
                                                        working_directory=args.run_directory)

    if args.grid_type == "dynamic":
        dynamic_grid_params = parse_inifile(run_parameters["psycris_inifile"])
        mesa_params_to_run_grid_over = dynamic_grid_params["posydon_dynamic_sampling_kwargs"]["mesa_column_names"]
        inlist_star1_formation, inlist_star2_formation, inlist_binary_project, inlist_star1_binary, \
            inlist_star2_binary = construct_static_inlist(mesa_inlists,
                                                          grid_parameters=mesa_params_to_run_grid_over,
                                                          working_directory=args.run_directory)
    else:
        inlist_star1_formation, inlist_star2_formation, inlist_binary_project, inlist_star1_binary, \
            inlist_star2_binary = construct_static_inlist(mesa_inlists,
                                                          grid_parameters=grid_df.columns,
                                                          working_directory=args.run_directory)

    column_lists_folder = os.path.join(args.run_directory, 'column_lists')
    if os.path.exists(column_lists_folder):
        shutil.rmtree(column_lists_folder)
    os.makedirs(column_lists_folder)

    star_history_columns = os.path.join(column_lists_folder, 'history_columns.list')
    binary_history_columns = os.path.join(column_lists_folder, 'binary_history_columns.list')
    profile_columns = os.path.join(column_lists_folder, 'profile_columns.list')

    shutil.copy(mesa_inlists['star_history_columns'], star_history_columns)
    shutil.copy(mesa_inlists['binary_history_columns'], binary_history_columns)
    shutil.copy(mesa_inlists['profile_columns'], profile_columns)

    if slurm['job_array']:
        command_line = construct_command_line(1,
                                              run_parameters['grid'],
                                              binary_exe,
                                              star1_exe,
                                              star2_exe,
                                              inlist_binary_project,
                                              inlist_star1_binary,
                                              inlist_star2_binary,
                                              inlist_star1_formation,
                                              inlist_star2_formation,
                                              star_history_columns,
                                              binary_history_columns,
                                              profile_columns,
                                              args.run_directory,
                                              'fixed',
                                              path_to_run_grid_exec,
                                              keep_profiles=run_parameters['keep_profiles'],
                                              keep_photos=run_parameters['keep_photos'])
        command_line += ' --grid-point-index $SLURM_ARRAY_TASK_ID'
    else:
        command_line = construct_command_line(slurm['number_of_mpi_tasks']*slurm['number_of_nodes'],
                                              fixgrid_file_name,
                                              binary_exe,
                                              star1_exe,
                                              star2_exe,
                                              inlist_binary_project,
                                              inlist_star1_binary,
                                              inlist_star2_binary,
                                              inlist_star1_formation,
                                              inlist_star2_formation,
                                              star_history_columns,
                                              binary_history_columns,
                                              profile_columns,
                                              args.run_directory,
                                              args.grid_type,
                                              path_to_run_grid_exec,
                                              psycris_inifile=run_parameters["psycris_inifile"],
                                              keep_profiles=run_parameters['keep_profiles'],
                                              keep_photos=run_parameters['keep_photos'])
    if args.submission_type == 'slurm':
        command_line += ' --job_end $SLURM_JOB_END_TIME'
    if 'work_dir' in slurm.keys() and not(slurm['work_dir'] == ''):
        command_line += ' --temporary-directory '+slurm['work_dir']

    if args.submission_type == 'shell':
        write_shell_submission_script(command_line, slurm, grid_df, args.run_directory)
    elif args.submission_type == 'slurm':
        write_slurm_submission_scripts(command_line, slurm, grid_df, args.run_directory)


def main() -> None:
    """Entry point for the posydon grid setup."""
    args = parse_commandline()
    run_setup(args)
