"""Write the ``posydon-run-grid`` command line and the scripts that run it.

Everything a run needs is passed to ``posydon-run-grid`` on one command
line (``command_line``). The scripts around it depend on the backend of
``[submission]`` (or ``--submission-type`` and ``job_array`` of a legacy
``[slurm]`` section)::

    backend      job_array  scripts written              how points run
    shell        True       grid_command.sh              a loop, one by one
    shell        False      grid_command.sh              one call over all
    slurm_array  -          job_array_grid_submit.slurm  one array task each
                            cleanup.slurm, run_grid.sh
    slurm_mpi    -          mpi_grid_submit.slurm        one call over all
                            cleanup.slurm, run_grid.sh

``run_grid.sh`` submits the grid script and then ``cleanup.slurm``, which
waits for it, compresses the output (``compress-mesa``) and, with
``newgroup``, hands the files to that group.

A job array runs a fixed grid, one point per task. A dynamic grid needs
one MPI job (``posydon-run-grid`` samples new points on rank 0 and runs
them on the others), so a dynamic grid with a job array is refused (fixed
S5: before, it silently ran as a fixed grid over the CSV).

All text is written exactly as today, including what the tracker lists: a
missing star inlist is passed as the literal ``None`` and the job name
holds a literal ``\\${SLURM_ARRAY_TASK_ID}`` (C11), only the ``.sh``
scripts are made executable (C12). A shell script without a job array ends
the command with a newline before ``compress-mesa`` (fixed S13: before,
the two were glued into one line).
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import shutil
from pathlib import Path

from posydon.utils.posydonwarning import Pwarn

# the options of posydon-run-grid, in today's order, and the name of the
# value each one takes
RUN_GRID_OPTIONS = (
    ("--mesa-grid", "grid"),
    ("--mesa-binary-executable", "binary_executable"),
    ("--mesa-star1-executable", "star1_executable"),
    ("--mesa-star2-executable", "star2_executable"),
    ("--mesa-binary-inlist-project", "binary_project"),
    ("--mesa-binary-inlist1", "binary_star1"),
    ("--mesa-binary-inlist2", "binary_star2"),
    ("--mesa-star1-inlist-project", "star1_formation"),
    ("--mesa-star2-inlist-project", "star2_formation"),
    ("--mesa-star-history-columns", "history_columns"),
    ("--mesa-binary-history-columns", "binary_history_columns"),
    ("--mesa-profile-columns", "profile_columns"),
    ("--output-directory", "run_directory"),
    ("--grid-type", "mode"),
    ("--psycris-inifile", "psycris_inifile"),
)
GRID_SCRIPTS = {"shell": "grid_command.sh",
                "slurm_array": "job_array_grid_submit.slurm",
                "slurm_mpi": "mpi_grid_submit.slurm"}
CLEANUP_SCRIPT = "cleanup.slurm"
RUN_SCRIPT = "run_grid.sh"


def find_run_grid():
    """The path of the ``posydon-run-grid`` executable.

    Raises
    ------
    ValueError
        If it is not on the ``PATH``.
    """
    path = shutil.which("posydon-run-grid")
    if path is None:
        raise ValueError("Cannot locate posydon-run-grid executable in "
                         "your path")
    return path


def one_point_per_task(submission):
    """Whether each grid point is its own call of ``posydon-run-grid``."""
    return submission.backend == "slurm_array" or \
        (submission.backend == "shell" and submission.job_array)


def run_grid_call(values, run_grid, mode="fixed", mpi_processes=1,
                  keep_profiles=False, keep_photos=False):
    """One call of ``posydon-run-grid``, as today's ``construct_command_line``.

    Parameters
    ----------
    values : dict
        A value for each name in ``RUN_GRID_OPTIONS``; ``None`` is written
        as ``None``.
    run_grid : str
        Path of ``posydon-run-grid``.
    mode : str
        ``fixed`` (one ``python`` process) or ``dynamic`` (``mpirun``).

    Examples
    --------
    >>> values = {name: name.upper() for _, name in RUN_GRID_OPTIONS}
    >>> run_grid_call(values, "run-grid").split(" --")[:2]
    ['python run-grid', 'mesa-grid GRID']
    """
    if mode == "fixed":
        start = "python {}".format(run_grid)
    elif mode == "dynamic":
        start = ("mpirun --bind-to none -np {} python -m mpi4py {}"
                 .format(mpi_processes, run_grid))
    else:
        raise ValueError("grid_type can either be fixed or dynamic not "
                         "anything else")
    values = dict(values, mode=mode)
    parts = [start] + ["{} {}".format(option, values[name])
                       for option, name in RUN_GRID_OPTIONS]
    if keep_profiles:
        parts.append("--keep_profiles")
    if keep_photos:
        parts.append("--keep_photos")
    return " ".join(parts)


def command_line(config, inlists, executables, column_lists, run_directory,
                 run_grid):
    """The full command line that runs the grid (or one point of it).

    Parameters
    ----------
    config : GridConfig
    inlists : BuiltInlists
        From ``build.build_inlists``.
    executables : tuple of str
        Binary, star 1 and star 2, from ``extras.compile_executables``.
    column_lists : dict
        From ``extras.stage_column_lists``.
    run_directory : str or Path
        Written as given.
    run_grid : str
        Path of ``posydon-run-grid`` (``find_run_grid``).
    """
    submission = config.submission
    values = {"grid": config.grid.points,
              "binary_executable": executables[0],
              "star1_executable": executables[1],
              "star2_executable": executables[2],
              "binary_project": inlists.binary_project,
              "binary_star1": inlists.binary_star1,
              "binary_star2": inlists.binary_star2,
              "star1_formation": inlists.star1_formation,
              "star2_formation": inlists.star2_formation,
              "history_columns": column_lists["history_columns"],
              "binary_history_columns":
                  column_lists["binary_history_columns"],
              "profile_columns": column_lists["profile_columns"],
              "run_directory": run_directory,
              "psycris_inifile": config.grid.psycris_inifile}
    keep = {"keep_profiles": config.grid.keep_profiles,
            "keep_photos": config.grid.keep_photos}
    if one_point_per_task(submission):
        if config.grid.mode == "dynamic":
            raise ValueError("a dynamic grid runs as one MPI job; it cannot "
                             "use a job array (S5)")
        # a fixed grid, one point per task; psy-cris is not passed
        values["psycris_inifile"] = None
        command = run_grid_call(values, run_grid, **keep)
        command += " --grid-point-index $SLURM_ARRAY_TASK_ID"
    else:
        command = run_grid_call(
            values, run_grid, mode=config.grid.mode,
            mpi_processes=submission.mpi_tasks * submission.nodes, **keep)
    if submission.backend != "shell":
        command += " --job_end $SLURM_JOB_END_TIME"
    if submission.work_dir not in (None, ""):
        command += " --temporary-directory " + submission.work_dir
    return command


###############################################################################
# the scripts
###############################################################################
def _environment(submission, mesasdk_root, mesa_dir):
    """Threads, MESA SDK and MESA_DIR, as at the top of each grid script."""
    return ("export OMP_NUM_THREADS={}\n\n"
            "export MESASDK_ROOT={}\n"
            "source $MESASDK_ROOT/bin/mesasdk_init.sh\n"
            "export MESA_DIR={}\n\n".format(submission.cpus_per_task,
                                            mesasdk_root, mesa_dir))


def _finish(submission):
    """Compress the output, hand it to ``newgroup``, say done."""
    text = "compress-mesa .\n"
    if submission.newgroup is not None:
        text += ('echo "Change group to {0}"\n'
                 'chgrp -fR {0} .\n'
                 'echo "Change group permission to rwX at least"\n'
                 'chmod -fR g+rwX .\n'.format(submission.newgroup))
    return text + '\necho "Done."'


def _sbatch(submission, lines):
    """``#!/bin/bash`` and the ``#SBATCH`` header around ``lines``."""
    head = ["--account={}".format(submission.account),
            "--partition={}".format(submission.partition)]
    tail = ["--mail-type=ALL",
            "--mail-user={}".format(submission.email),
            "--mem-per-cpu={}".format(submission.mem_per_cpu)]
    return "#!/bin/bash\n" + "".join(
        "#SBATCH {}\n".format(line) for line in head + lines + tail) + "\n"


def shell_script(config, command, n_points, mesasdk_root, mesa_dir):
    """``grid_command.sh``: the whole grid in one shell script."""
    submission = config.submission
    text = "#!/bin/bash\n\n" + _environment(submission, mesasdk_root,
                                             mesa_dir)
    if submission.job_array:
        text += ("for SLURM_ARRAY_TASK_ID in "
                 + "".join("{} ".format(i) for i in range(n_points))
                 + "; do " + command + " ; done\n\n")
    else:
        # the command on its own line (S13: before, compress-mesa was glued
        # to its last argument)
        text += command + "\n\n"
    return text + _finish(submission)


def slurm_grid_script(config, command, n_points, mesasdk_root, mesa_dir):
    """The SLURM script of the grid: a job array or one MPI job."""
    submission = config.submission
    if submission.backend == "slurm_array":
        header = _sbatch(submission, [
            "-N 1",
            "--array=0-{}".format(n_points - 1),
            "--cpus-per-task {}".format(submission.cpus_per_task),
            "--ntasks-per-node 1",
            "--time={}".format(submission.walltime),
            # a literal \$ (C11)
            '--job-name="mesa_grid_\\${SLURM_ARRAY_TASK_ID}"',
            "--output=mesa_grid.%A_%a.out"])
    else:
        header = _sbatch(submission, [
            "-N {}".format(submission.nodes),
            "--cpus-per-task {}".format(submission.cpus_per_task),
            "--ntasks-per-node {}".format(submission.mpi_tasks),
            "--time={}".format(submission.walltime),
            '--output="mesa_grid.out"'])
    return header + _environment(submission, mesasdk_root, mesa_dir) \
        + "\n" + command


def cleanup_script(config):
    """``cleanup.slurm``: runs after the grid, on one CPU."""
    submission = config.submission
    return _sbatch(submission, [
        "-N 1",
        "--cpus-per-task 1",
        "--ntasks-per-node 1",
        "--time={}".format(submission.walltime),
        '--job-name="mesa_grid_cleanup"',
        "--output=mesa_cleanup.out"]) + _finish(submission)


def run_script(grid_script):
    """``run_grid.sh``: submit the grid, then the cleanup after it."""
    return ("#!/bin/bash\n"
            "ID_GRID=$(sbatch --parsable {0})\n"
            'echo "{0} submitted as "${{ID_GRID}}\n'
            "ID_cleanup=$(sbatch --parsable --dependency=afterany:${{ID_GRID}}"
            " --kill-on-invalid-dep=yes {1})\n"
            'echo "{1} submitted as "${{ID_cleanup}}\n'
            .format(grid_script, CLEANUP_SCRIPT))


def _write(path, text, executable=False):
    if path.exists():
        Pwarn("Replace " + path.name, "OverwriteWarning")
    path.write_text(text)
    if executable:
        path.chmod(0o755)
    return path


def write_scripts(config, command, n_points, folder, mesasdk_root,
                  mesa_dir):
    """Write the scripts of the backend into ``folder``.

    Parameters
    ----------
    config : GridConfig
    command : str
        From ``command_line``.
    n_points : int
        Number of grid points (rows of the grid CSV).
    folder : str or Path
        Where the scripts go. Today this is the current directory, not the
        run directory (S9).
    mesasdk_root, mesa_dir : str
        Written into the scripts; today ``$MESASDK_ROOT`` and ``$MESA_DIR``.

    Returns
    -------
    list of Path
        The scripts, in the order they are written.
    """
    folder = Path(folder)
    backend = config.submission.backend
    grid_script = GRID_SCRIPTS[backend]
    if backend == "shell":
        return [_write(folder / grid_script,
                       shell_script(config, command, n_points, mesasdk_root,
                                    mesa_dir), executable=True)]
    return [_write(folder / grid_script,
                   slurm_grid_script(config, command, n_points,
                                     mesasdk_root, mesa_dir)),
            _write(folder / CLEANUP_SCRIPT, cleanup_script(config)),
            _write(folder / RUN_SCRIPT, run_script(grid_script),
                   executable=True)]
