"""Submission script generation for the posydon-grid setup CLI."""

import os

import pandas

from posydon.utils.posydonwarning import Pwarn


def _run_shell(command: str) -> None:
    """Run a shell command via os.system.

    Args:
        command: The shell command to execute.
    """
    os.system(command)


def construct_command_line(
    number_of_mpi_processes: int,
    path_to_grid: str,
    binary_exe: str,
    star1_exe: str,
    star2_exe: str,
    inlist_binary_project: str,
    inlist_star1_binary: str,
    inlist_star2_binary: str,
    inlist_star1_formation,
    inlist_star2_formation,
    star_history_columns: str,
    binary_history_columns: str,
    profile_columns: str,
    run_directory: str,
    grid_type: str,
    path_to_run_grid_exec: str,
    psycris_inifile=None,
    keep_profiles: bool = False,
    keep_photos: bool = False,
) -> str:
    """Based on the inifile construct the command line call to posydon-run-grid.

    Args:
        number_of_mpi_processes: Number of MPI processes for dynamic grids.
        path_to_grid: Path to the grid file.
        binary_exe: Path to the compiled binary executable.
        star1_exe: Path to the compiled star1 executable.
        star2_exe: Path to the compiled star2 executable.
        inlist_binary_project: Path to the binary project inlist.
        inlist_star1_binary: Path to the star1 binary inlist.
        inlist_star2_binary: Path to the star2 binary inlist.
        inlist_star1_formation: Path(s) to the star1 formation inlists.
        inlist_star2_formation: Path(s) to the star2 formation inlists.
        star_history_columns: Path to the star history columns list.
        binary_history_columns: Path to the binary history columns list.
        profile_columns: Path to the profile columns list.
        run_directory: Directory to write the grid output to.
        grid_type: Either 'fixed' or 'dynamic'.
        path_to_run_grid_exec: Path to the posydon-run-grid executable.
        psycris_inifile: Path to the psycris inifile for dynamic grids.
        keep_profiles: Whether to keep profile files.
        keep_photos: Whether to keep photo files.

    Returns:
        The constructed command line string.
    """
    if grid_type == "fixed":
        command_line = 'python {15} --mesa-grid {1} --mesa-binary-executable {2} '
    elif grid_type == "dynamic":
        command_line = 'mpirun --bind-to none -np {0} python -m mpi4py {15} --mesa-grid {1} --mesa-binary-executable {2} '
    else:
        raise ValueError("grid_type can either be fixed or dynamic not anything else")
    command_line += '--mesa-star1-executable {3} --mesa-star2-executable {4} --mesa-binary-inlist-project {5} '
    command_line += '--mesa-binary-inlist1 {6} --mesa-binary-inlist2 {7} --mesa-star1-inlist-project {8} '
    command_line += '--mesa-star2-inlist-project {9} --mesa-star-history-columns {10} '
    command_line += '--mesa-binary-history-columns {11} --mesa-profile-columns {12} '
    command_line += '--output-directory {13} --grid-type {14} '
    command_line += '--psycris-inifile {16}'
    if keep_profiles:
        command_line += ' --keep_profiles'
    if keep_photos:
        command_line += ' --keep_photos'
    command_line = command_line.format(number_of_mpi_processes,
                                       path_to_grid,
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
                                       run_directory,
                                       grid_type,
                                       path_to_run_grid_exec,
                                       psycris_inifile)
    return command_line


def write_shell_submission_script(
    command_line: str, slurm: dict, grid_df: pandas.DataFrame, run_directory: str
) -> None:
    """Write the shell submission script grid_command.sh to the current directory.

    Args:
        command_line: The posydon-run-grid command line to embed in the script.
        slurm: The parsed [slurm] inifile section.
        grid_df: The grid points DataFrame.
        run_directory: Unused; kept for interface consistency.
    """
    if os.path.exists('grid_command.sh'):
        Pwarn('Replace grid_command.sh', "OverwriteWarning")
    with open('grid_command.sh', 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write('export OMP_NUM_THREADS={0}\n\n'.format(slurm['number_of_cpus_per_task']))
        f.write('export MESASDK_ROOT={0}\n'.format(os.environ['MESASDK_ROOT']))
        f.write('source $MESASDK_ROOT/bin/mesasdk_init.sh\n')
        f.write('export MESA_DIR={0}\n\n'.format(os.environ['MESA_DIR']))
        if slurm['job_array']:
            f.write('for SLURM_ARRAY_TASK_ID in ')
            for i in range(len(grid_df)):
                f.write('{0} '.format(i))
            f.write('; do ')
        f.write(command_line)
        if slurm['job_array']:
            f.write(' ; done\n\n')
        f.write('compress-mesa .\n')
        if 'newgroup' in slurm.keys():
            f.write('echo "Change group to {0}"\n'.format(slurm['newgroup']))
            f.write('chgrp -fR {0} .\n'.format(slurm['newgroup']))
            f.write('echo "Change group permission to rwX at least"\n')
            f.write('chmod -fR g+rwX .\n')
        f.write('\necho "Done."')
    _run_shell("chmod 755 grid_command.sh")


def write_slurm_submission_scripts(
    command_line: str, slurm: dict, grid_df: pandas.DataFrame, run_directory: str
) -> tuple[str, str]:
    """Write the slurm submission, cleanup and run scripts to the current directory.

    Args:
        command_line: The posydon-run-grid command line to embed in the scripts.
        slurm: The parsed [slurm] inifile section.
        grid_df: The grid points DataFrame.
        run_directory: Unused; kept for interface consistency.

    Returns:
        Tuple of the grid submission script name and the run script name.
    """
    if slurm['job_array']:
        grid_script = 'job_array_grid_submit.slurm'
        if os.path.exists(grid_script):
            Pwarn('Replace '+grid_script, "OverwriteWarning")
        with open(grid_script, 'w') as f:
            f.write('#!/bin/bash\n')

            f.write('#SBATCH --account={0}\n'.format(slurm['account']))
            f.write('#SBATCH --partition={0}\n'.format(slurm['partition']))
            f.write('#SBATCH -N 1\n')
            f.write('#SBATCH --array=0-{0}\n'.format(len(grid_df)-1))
            f.write('#SBATCH --cpus-per-task {0}\n'.format(slurm['number_of_cpus_per_task']))
            f.write('#SBATCH --ntasks-per-node 1\n')
            f.write('#SBATCH --time={0}\n'.format(slurm['walltime']))
            f.write('#SBATCH --job-name="mesa_grid_\${SLURM_ARRAY_TASK_ID}"\n')
            f.write('#SBATCH --output=mesa_grid.%A_%a.out\n')
            f.write('#SBATCH --mail-type=ALL\n')
            f.write('#SBATCH --mail-user={0}\n'.format(slurm['email']))
            f.write('#SBATCH --mem-per-cpu=4G\n\n')

            f.write('export OMP_NUM_THREADS={0}\n\n'.format(slurm['number_of_cpus_per_task']))

            f.write('export MESASDK_ROOT={0}\n'.format(os.environ['MESASDK_ROOT']))
            f.write('source $MESASDK_ROOT/bin/mesasdk_init.sh\n')
            f.write('export MESA_DIR={0}\n\n\n'.format(os.environ['MESA_DIR']))
            f.write(command_line)
    else:
        grid_script = 'mpi_grid_submit.slurm'
        if os.path.exists(grid_script):
            Pwarn('Replace '+grid_script, "OverwriteWarning")
        with open(grid_script, 'w') as f:
            f.write('#!/bin/bash\n')

            f.write('#SBATCH --account={0}\n'.format(slurm['account']))
            f.write('#SBATCH --partition={0}\n'.format(slurm['partition']))
            f.write('#SBATCH -N {0}\n'.format(slurm['number_of_nodes']))
            f.write('#SBATCH --cpus-per-task {0}\n'.format(slurm['number_of_cpus_per_task']))
            f.write('#SBATCH --ntasks-per-node {0}\n'.format(slurm['number_of_mpi_tasks']))
            f.write('#SBATCH --time={0}\n'.format(slurm['walltime']))
            f.write('#SBATCH --output="mesa_grid.out"\n')
            f.write('#SBATCH --mail-type=ALL\n')
            f.write('#SBATCH --mail-user={0}\n'.format(slurm['email']))
            f.write('#SBATCH --mem-per-cpu=4G\n\n')

            f.write('export OMP_NUM_THREADS={0}\n\n'.format(slurm['number_of_cpus_per_task']))

            f.write('export MESASDK_ROOT={0}\n'.format(os.environ['MESASDK_ROOT']))
            f.write('source $MESASDK_ROOT/bin/mesasdk_init.sh\n')
            f.write('export MESA_DIR={0}\n\n\n'.format(os.environ['MESA_DIR']))
            f.write(command_line)
    if os.path.exists('cleanup.slurm'):
        Pwarn('Replace cleanup.slurm', "OverwriteWarning")
    with open('cleanup.slurm', 'w') as f:
        f.write('#!/bin/bash\n')

        f.write('#SBATCH --account={0}\n'.format(slurm['account']))
        f.write('#SBATCH --partition={0}\n'.format(slurm['partition']))
        f.write('#SBATCH -N 1\n')
        f.write('#SBATCH --cpus-per-task 1\n')
        f.write('#SBATCH --ntasks-per-node 1\n')
        f.write('#SBATCH --time={0}\n'.format(slurm['walltime']))
        f.write('#SBATCH --job-name="mesa_grid_cleanup"\n')
        f.write('#SBATCH --output=mesa_cleanup.out\n')
        f.write('#SBATCH --mail-type=ALL\n')
        f.write('#SBATCH --mail-user={0}\n'.format(slurm['email']))
        f.write('#SBATCH --mem-per-cpu=4G\n\n')

        f.write('compress-mesa .\n')
        if 'newgroup' in slurm.keys():
            f.write('echo "Change group to {0}"\n'.format(slurm['newgroup']))
            f.write('chgrp -fR {0} .\n'.format(slurm['newgroup']))
            f.write('echo "Change group permission to rwX at least"\n')
            f.write('chmod -fR g+rwX .\n')
        f.write('\necho "Done."')
    if os.path.exists('run_grid.sh'):
        Pwarn('Replace run_grid.sh', "OverwriteWarning")
    with open('run_grid.sh', 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('ID_GRID=$(sbatch --parsable {0})\n'.format(grid_script))
        f.write('echo "{0}'.format(grid_script)+' submitted as "${ID_GRID}\n')
        f.write('ID_cleanup=$(sbatch --parsable --dependency=afterany:${ID_GRID} '
                '--kill-on-invalid-dep=yes cleanup.slurm)\n')
        f.write('echo "cleanup.slurm submitted as "${ID_cleanup}\n')
    _run_shell("chmod 755 run_grid.sh")
    return grid_script, "run_grid.sh"
