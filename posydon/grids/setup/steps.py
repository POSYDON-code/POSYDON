"""Set up a grid of MESA runs from an ini file: the steps, in order.

This is what ``posydon-grid setup`` and ``explain`` run, and what a script
or notebook can call directly (``setup_grid("grid.ini")``). The work is
done by the other modules of ``posydon.grids.setup``:

    1. read and check the ini file     config.read_config (new format) or
                                       legacy.read_legacy_config
    2. read the grid points            (CSV, PSyGrid .h5 or a MESA grid folder)
    3. get the scenario's inlists      resolve.resolve_scenario; defaults of a
                                       new-format ini: defaults.fill_defaults
    4. stage sources and compile       extras.stage_extras, compile_executables
    5. write the inlists               build.build_inlists
    6. stage the column lists          extras.stage_column_lists
    7. write the submission scripts    submit.command_line, write_scripts
    8. record what it was set up from  record.write_record (setup_record.ini)

``prepare`` does steps 1-3, which check everything; ``explain_grid``
stops there, without step 2, and says where each inlist value comes from;
``setup_grid`` goes on to write the run directory. Only reading the grid
points can write: a MESA grid folder is first converted into
``fixed_grid_results.h5`` (in the current directory, as today).

An ini with a ``[grid]`` section is a new-format ini: it says the mode, run
directory and backend itself. Any other ini is a legacy ini, which needs
the grid type (and takes the run directory and submission type) from the
command line.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import configparser
import os
from dataclasses import dataclass
from pathlib import Path

import pandas

from posydon.grids.setup.build import build_inlists
from posydon.grids.setup.config import ConfigError, read_config
from posydon.grids.setup.defaults import fill_defaults
from posydon.grids.setup.explain import explain_text
from posydon.grids.setup.extras import (
    COLUMN_LIST_FILES,
    column_list_sources,
    compile_executables,
    keep_highest_extras,
    merge_extras,
    place,
    replaced_column_lists,
    stage_column_lists,
    stage_extras,
    whole_file_extras,
)
from posydon.grids.setup.legacy import read_legacy_config
from posydon.grids.setup.record import write_record
from posydon.grids.setup.resolve import resolve_scenario
from posydon.grids.setup.submit import (
    command_line,
    find_run_grid,
    one_point_per_task,
    write_scripts,
)
from posydon.grids.setup.support import support_root
from posydon.utils.posydonwarning import Pwarn

SUBMISSION_TYPES = ("shell", "slurm")


def is_new_format(path):
    """Whether an ini file is in the new format (it has a [grid] section)."""
    parser = configparser.ConfigParser(interpolation=None, strict=False)
    parser.read(path)
    return parser.has_section("grid")


def read_any_config(inifile, grid_type=None, run_directory=None,
                    submission_type=None, overrides=False):
    """Read the ini in either format.

    Parameters
    ----------
    inifile : str or Path
    grid_type, run_directory, submission_type : str or None
        The command-line options. A legacy ini needs ``grid_type``; the
        other two default to the current directory and ``shell``.
    overrides : bool
        For a new-format ini: let the options override ``[grid] mode``,
        ``run_directory`` and ``[submission] backend`` (the deprecated
        ``posydon-setup-grid`` does); otherwise they are refused, since the
        ini says all of it.
    """
    if not is_new_format(inifile):
        if grid_type is None:
            raise ValueError("--grid-type is needed for a legacy ini file")
        return read_legacy_config(
            inifile, grid_type=grid_type,
            submission_type=submission_type or 'shell',
            run_directory=run_directory or os.getcwd(),
            # grid_params/ paths: from POSYDON-MESA-INLISTS at a pinned
            # commit, since POSYDON no longer ships grid_params/
            support=support_root)
    given = [option for option, value in (
        ("--grid-type", grid_type), ("--run-directory", run_directory),
        ("--submission-type", submission_type)) if value is not None]
    if given and not overrides:
        raise ValueError("{} only {} to legacy ini files; a new-format "
                         "ini says this in [grid] mode, run_directory and "
                         "[submission] backend".format(
                             ", ".join(given),
                             "applies" if len(given) == 1 else "apply"))
    config = read_config(inifile)
    if grid_type is not None:
        config.grid.mode = grid_type
    if run_directory is not None:
        config.grid.run_directory = Path(run_directory)
    if submission_type == 'shell':
        config.submission.backend = 'shell'
    elif submission_type == 'slurm' and \
            config.submission.backend == 'shell':
        config.submission.backend = 'slurm_array' \
            if config.submission.job_array else 'slurm_mpi'
    return config


def check_config(config):
    """Stop on settings that would fail later, before anything is written."""
    grid = config.grid
    if not (grid.points.is_file() or grid.points.is_dir()):
        raise ValueError("Supplied grid does not exist, please check your "
                         "path and try again")
    if grid.mode == "dynamic" and grid.psycris_inifile is None:
        raise ValueError("Please add psycris inifile to the [run_parameters] "
                         "section of the inifile.")
    # named as the ini says it: [slurm] job_array in a legacy ini,
    # [submission] backend in a new-format one
    legacy = config.source_format == "legacy"
    where = ("slurm", "job_array") if legacy else ("submission", "backend")
    if grid.mode == "fixed" and not one_point_per_task(config.submission):
        # posydon-run-grid runs a fixed grid one point per call
        # (--grid-point-index); without it nothing would run (tracker S4, R1)
        raise ConfigError(*where, (
            "a fixed grid needs job_array = True" if legacy else
            "a fixed grid runs one point per task: use slurm_array, or "
            "shell with job_array = true"))
    if grid.mode == "dynamic" and one_point_per_task(config.submission):
        # posydon-run-grid samples a dynamic grid on rank 0 of one MPI job;
        # a job array cannot run it (tracker S5: it silently ran as fixed)
        raise ConfigError(*where, (
            "a dynamic grid runs as one MPI job: set job_array = False (new "
            "format: [submission] backend = slurm_mpi)" if legacy else
            "a dynamic grid runs as one MPI job: use slurm_mpi, or shell "
            "with job_array = false"))
    if config.scenario.is_set and config.grid.type is None:
        raise ConfigError("mesa_inlists", "scenario",
                          "unknown system type {!r}".format(
                              config.scenario.system_type))
    for note in config.notes:
        # pointing legacy users to the new format is done by the command
        if not note.startswith("legacy ini format"):
            Pwarn(note, "InappropriateValueWarning")


def read_grid_points(grid):
    """The grid points as a DataFrame, and the grid path for MPI runs."""
    path = str(grid.points)
    if '.csv' in path:
        return pandas.read_csv(path), path
    # PSyGrid pulls in matplotlib; only import it when needed (tracker C14)
    from posydon.grids.psygrid import PSyGrid
    if '.h5' in path:
        psy_grid = PSyGrid()
        psy_grid.load(path)
        grid_df = psy_grid.get_pandas_initial_final()
        psy_grid.close()
        return grid_df, path
    if os.path.isdir(path):
        h5_path = os.path.join(os.getcwd(), "fixed_grid_results.h5")
        PSyGrid().create(path, h5_path, slim=True)
        psy_grid = PSyGrid()
        psy_grid.load(h5_path)
        grid_df = psy_grid.get_pandas_initial_final()
        psy_grid.close()
        return grid_df, h5_path
    raise ValueError('Grid format not recognized, please feed in an '
                     'acceptable format: csv')


def grid_parameters(config, grid_df):
    """The names the grid varies: the CSV columns, or psy-cris's."""
    if config.grid.mode == "dynamic":
        from posydon.active_learning.psy_cris.utils import parse_inifile
        psycris = parse_inifile(str(config.grid.psycris_inifile))
        return psycris["posydon_dynamic_sampling_kwargs"]["mesa_column_names"]
    return grid_df.columns


def chosen_extras(config, scenario):
    """The extras of the ini, the scenario and the inlist layers' whole
    files, one file per place.

    A dropped key is reported only when the ini set it: every key of a
    legacy ini (today's warning), an ``[extras]`` key of a new-format one.
    MESA's stock files filled in by ``fill_defaults``, and scenario files a
    user file replaces, are dropped silently: that is what they are for.

    A whole file of ``[mesa_inlists]`` or ``[posydon_inlists]`` loses to the
    scenario's file at the same place, as the column lists do (before, it
    replaced the release's ``mesa_star1_extras``, so star 1 and star 2
    compiled it while the binary's star took the release's).
    """
    scenario_extras = scenario.extras if scenario is not None else {}
    extras = merge_extras(config.extras, scenario_extras)
    # the places the scenario fills, star 2's copy of star 1's included
    taken = {place(key, path) for key, path in
             merge_extras({}, scenario_extras).items()}
    replaced = []
    for key, path in whole_file_extras(config).items():
        if not key.startswith("user_") and place(key, path) in taken:
            if key.startswith("posydon_") and path not in replaced:
                replaced.append(path)
            continue
        extras[key] = path
    for path in replaced:
        Pwarn("[posydon_inlists] {} is not used: the scenario's file is "
              "used instead".format(path), "ReplaceValueWarning")
    kept = keep_highest_extras(extras)
    for key in extras:
        if key in kept:
            continue
        if config.source_format == "legacy":
            Pwarn("Section mesa_extras value {0} is being set to None"
                  .format(key), "ReplaceValueWarning")
        elif key in config.extras and key not in config.default_extras:
            Pwarn("[extras] {} is not used: a higher layer's file goes "
                  "to the same place".format(key), "ReplaceValueWarning")
    return kept


def _environment(*names):
    for name in names:
        if name not in os.environ:
            raise ValueError("{} must be defined in your environment "
                             "before you can set up a grid of MESA runs"
                             .format(name))


@dataclass
class Prepared:
    """Everything read and checked for a setup; nothing written yet."""

    config: object
    run_directory: str
    grid_df: object
    resolved: object
    extras: dict
    columns: dict
    parameters: object

    @property
    def scenario(self):
        return self.resolved.files if self.resolved is not None else None


def prepare(inifile, grid_type=None, run_directory=None,
            submission_type=None, overrides=False, read_points=True):
    """Steps 1-3: read and check everything.

    Parameters are those of ``read_any_config``, and

    read_points : bool
        Read the grid points (step 2). ``explain_grid`` does not need them:
        it checks that they exist, but does not read them, so it writes
        nothing even for a MESA grid folder (which ``read_grid_points``
        converts into an HDF5 file first). Without them, ``grid_df`` and
        ``parameters`` are None.

    Returns
    -------
    Prepared
    """
    _environment("MESA_DIR")
    config = read_any_config(inifile, grid_type, run_directory,
                             submission_type, overrides)
    # written into the inlists and scripts as given (legacy: --run-directory)
    run_dir = run_directory or str(config.grid.run_directory)
    check_config(config)
    grid_df = parameters = None
    if read_points:
        grid_df, mpi_grid = read_grid_points(config.grid)
        if not one_point_per_task(config.submission):
            config.grid.points = Path(mpi_grid)
        parameters = grid_parameters(config, grid_df)
    resolved = resolve_scenario(config)
    scenario = resolved.files if resolved is not None else None
    fill_defaults(config, resolved, os.environ['MESA_DIR'])
    extras = chosen_extras(config, scenario)
    columns = column_list_sources(config, scenario)
    for name in replaced_column_lists(config, scenario):
        # the scenario's column lists win over the ini's (tracker C7)
        Pwarn("{} in the ini file is not used: the scenario's {} is "
              "used instead".format(name, columns[name]),
              "ReplaceValueWarning")
    missing = [name for name in COLUMN_LIST_FILES if name not in columns]
    if missing:
        raise ConfigError("mesa_inlists", None, "no file for the column "
                          "lists {}".format(", ".join(missing)))
    return Prepared(config=config, run_directory=run_dir, grid_df=grid_df,
                    resolved=resolved, extras=extras, columns=columns,
                    parameters=parameters)


def explain_grid(inifile, **options):
    """Where each inlist value of the grid comes from; nothing is written.

    Options are those of ``read_any_config``.
    """
    ready = prepare(inifile, read_points=False, **options)
    roots = [ready.resolved.repo] if ready.resolved is not None else []
    return explain_text(ready.config, ready.scenario, roots)


def setup_grid(inifile, **options):
    """Set up the run directory of a grid; return it.

    Options are those of ``read_any_config``.
    """
    _environment("MESA_DIR", "MESASDK_ROOT")
    run_grid = find_run_grid()
    ready = prepare(inifile, **options)
    run_directory = ready.run_directory

    # 4-6. write the run directory
    stage_extras(ready.extras, run_directory)
    executables = compile_executables(run_directory)
    inlists = build_inlists(ready.config, run_directory, ready.parameters,
                            scenario=ready.scenario)
    staged_columns = stage_column_lists(ready.columns, run_directory)

    # 7. the scripts that run it, next to the run (tracker S9)
    command = command_line(ready.config, inlists, executables,
                           staged_columns, run_directory, run_grid)
    write_scripts(ready.config, command, len(ready.grid_df), run_directory,
                  os.environ['MESASDK_ROOT'], os.environ['MESA_DIR'])

    # 8. what it was set up from: ini, POSYDON, inlist commit (tracker C13)
    from posydon import __version__ as posydon_version
    write_record(run_directory, ready.config, ready.resolved, inifile,
                 os.environ['MESA_DIR'], posydon_version)
    return Path(run_directory)
