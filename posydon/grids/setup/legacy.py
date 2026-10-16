"""Read today's grid ini files into a ``GridConfig``.

Today's ini files have four sections: ``[slurm]``, ``[mesa_inlists]``,
``[mesa_extras]`` and ``[run_parameters]``. They are read with the same
function ``posydon-setup-grid`` uses today
(``posydon.utils.configfile.parse_inifile``), so every value is understood
exactly as before, and then mapped onto the new settings::

    [run_parameters] grid, psycris_inifile, keep_*   -> config.grid
    [mesa_inlists]   scenario                        -> config.scenario
                     zams_filename*, single_star_grid -> config.grid
                     <target>_mesa_defaults           -> mesa_inlists layer
                     <target>_posydon_defaults        \\
                     <target>_special                  > posydon_inlists layer
                     *_columns                        /
                     <target>_user                    -> user_inlists layer
                     history_*, final_*, binary_*     -> config.output
    [mesa_extras]                                     -> config.extras
    [slurm]                                           -> config.submission

Three things are not in the file and come from the command line instead:
``--grid-type`` (fixed or dynamic), ``--submission-type`` (shell or slurm)
and ``--run-directory``.

The grid type is worked out from the scenario's system type and
``single_star_grid``; ``LEGACY_GRID_TYPES`` lists the combinations. The
non-existent system type ``MS-MS``, used by the main v2 HMS-HMS grids, is
read as HMS-HMS with a note.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import os
import re
from pathlib import Path

from posydon.grids.setup.config import ConfigError, GridConfig, SubmissionSection
from posydon.utils.configfile import parse_inifile

# (system type in the scenario, single_star_grid) -> grid type
LEGACY_GRID_TYPES = {
    ("HMS-HMS", False): "HMS-HMS",
    ("MS-MS", False): "HMS-HMS",
    ("CO-H_star", False): "CO-HMS_RLO",
    ("CO-He_star", False): "CO-HeMS",
    ("HeMS-HMS", False): "HeMS-HMS",
    ("HMS-HMS", True): "single_HMS",
    ("MS-MS", True): "single_HMS",
    ("CO-He_star", True): "single_HeMS",
}

# suffix of an inlist key -> the layer its files belong to
LEGACY_LAYERS = {
    "mesa_defaults": "mesa_inlists",
    "posydon_defaults": "posydon_inlists",
    "special": "posydon_inlists",
    "user": "user_inlists",
}

# "star1_formation_controls_user" -> star1, formation, controls, user
INLIST_KEY = re.compile(
    r"^(?P<who>binary|star[12])_(?P<formation>formation_)?"
    r"(?P<part>controls|job)_(?P<suffix>"
    + "|".join(LEGACY_LAYERS) + r")$")

COLUMN_KEYS = {"star_history_columns": "history_columns",
               "binary_history_columns": "binary_history_columns",
               "profile_columns": "profile_columns"}

OUTPUT_KEYS = ("history_interval", "binary_history", "history_star1",
               "final_profile_star1", "final_model_star1", "history_star2",
               "final_profile_star2", "final_model_star2")

# keys that are read but need no mapping of their own
KNOWN_KEYS = {"posydon_github_root", "mesa_dir"}


def split_scenario_ref(ref):
    """Split a legacy ``'branch-commit'`` into branch and commit.

    The commit is everything after the *last* dash, so branch names that
    contain dashes work.

    Examples
    --------
    >>> split_scenario_ref("main-18710e943edd926a5653c4cdb7d6e18e5bdb35a2")
    ('main', '18710e943edd926a5653c4cdb7d6e18e5bdb35a2')
    >>> split_scenario_ref("CO-HeMS_v1_rerun-0e736a8")
    ('CO-HeMS_v1_rerun', '0e736a8')
    """
    branch, dash, commit = ref.rpartition("-")
    if not dash or not branch or not re.fullmatch(r"[0-9a-f]{7,40}",
                                                  commit):
        raise ConfigError("mesa_inlists", "scenario", "expected "
                          "'branch-commit' with a 7 to 40 digit hex commit, "
                          "got {!r}".format(ref))
    return branch, commit


def _as_list(value):
    return value if isinstance(value, list) else [value]


def _read_with_mesa_dir(path, mesa_dir):
    """Run parse_inifile, which needs $MESA_DIR for ${MESA_DIR} values."""
    if mesa_dir is None and "MESA_DIR" not in os.environ:
        raise ConfigError("mesa_extras", None, "MESA_DIR is not set; pass "
                          "mesa_dir or set the environment variable")
    previous = os.environ.get("MESA_DIR")
    if mesa_dir is not None:
        os.environ["MESA_DIR"] = str(mesa_dir)
    try:
        return parse_inifile(str(path))
    except KeyError as error:
        raise ConfigError(str(error).strip("'"), None, "missing section")
    finally:
        if mesa_dir is not None:
            if previous is None:
                del os.environ["MESA_DIR"]
            else:
                os.environ["MESA_DIR"] = previous


def _read_inlist_keys(mesa_inlists, config, single_star):
    """Sort the [mesa_inlists] file keys into layers and targets.

    Formation keys hold either a list (one file per step) or one file that
    applies to every step. The number of steps is only known once the
    scenario is resolved, so a single file goes to the target
    ``<star>_formation_every_step_<part>``.
    """
    for key, value in mesa_inlists.items():
        match = INLIST_KEY.match(key)
        if match is None or value is None:
            continue
        who, part = match.group("who"), match.group("part")
        files = config.inlists[LEGACY_LAYERS[match.group("suffix")]]
        # single-star grids run star 1's own inlists as formation steps
        is_steps = match.group("formation") or (single_star and
                                                who == "star1")
        if not is_steps:
            for path in _as_list(value):
                files.add("{}_{}".format(who, part), Path(path))
        elif isinstance(value, list):
            for step, path in enumerate(value):
                files.add("{}_formation_step{}_{}".format(who, step, part),
                          Path(path))
        else:
            files.add("{}_formation_every_step_{}".format(who, part),
                      Path(value))


def read_legacy_config(path, grid_type="fixed", submission_type=None,
                       run_directory=None, mesa_dir=None):
    """Read a legacy grid ini file into a ``GridConfig``.

    Parameters
    ----------
    path : str or Path
        The ini file.
    grid_type : str
        ``--grid-type`` of ``posydon-setup-grid``: fixed or dynamic.
    submission_type : str or None
        ``--submission-type``: shell or slurm. ``None`` picks slurm (job
        array or MPI, from ``job_array``).
    run_directory : str, Path or None
        ``--run-directory``.
    mesa_dir : str, Path or None
        Value for ``${MESA_DIR}`` in the file; defaults to ``$MESA_DIR``.

    Returns
    -------
    GridConfig
        With ``source_format = "legacy"``, the raw sections in ``legacy``
        and anything worth telling the user in ``notes``.
    """
    path = Path(path)
    if not path.is_file():
        raise ConfigError("file", None, "{} not found".format(path))
    run_parameters, slurm, mesa_inlists, mesa_extras = \
        _read_with_mesa_dir(path, mesa_dir)
    config = GridConfig(source_format="legacy")
    config.legacy = {"run_parameters": run_parameters, "slurm": slurm,
                     "mesa_inlists": mesa_inlists, "mesa_extras": mesa_extras}
    config.notes.append("legacy ini format; see the grid documentation for "
                        "the new [grid]/[scenario]/[user_inlists] format")
    single_star = bool(mesa_inlists.get("single_star_grid", False))

    # [mesa_inlists] scenario -> scenario and grid type
    if mesa_inlists.get("scenario") is not None:
        source, ref, system_type = mesa_inlists["scenario"]
        branch, commit = split_scenario_ref(ref)
        config.scenario.source = source
        config.scenario.branch = branch
        config.scenario.commit = commit
        config.scenario.system_type = system_type
        grid_type_name = LEGACY_GRID_TYPES.get((system_type, single_star))
        if grid_type_name is None:
            config.notes.append("unknown system type {!r} with "
                                "single_star_grid = {}; grid type not set"
                                .format(system_type, single_star))
        elif system_type == "MS-MS":
            config.notes.append("system type 'MS-MS' does not exist; read "
                                "as {}".format(grid_type_name))
        config.grid.type = grid_type_name
    else:
        config.notes.append("no scenario: the grid type is not stated in "
                            "the file")

    # [mesa_inlists] files, ZAMS models and output options
    _read_inlist_keys(mesa_inlists, config, single_star)
    for key, name in COLUMN_KEYS.items():
        if mesa_inlists.get(key) is not None:
            config.inlists["posydon_inlists"].whole_files[name] = \
                Path(mesa_inlists[key])
    zams_both = mesa_inlists.get("zams_filename")
    for star in ("1", "2"):
        zams = zams_both if zams_both is not None else \
            mesa_inlists.get("zams_filename_" + star)
        if zams is not None:
            setattr(config.grid, "zams_star" + star, Path(zams))
    output = config.output
    output.history_interval = mesa_inlists.get("history_interval", 1)
    output.binary_history = bool(mesa_inlists.get("binary_history", True))
    for star in ("star1", "star2"):
        # e.g. history_star1, final_profile_star1, final_model_star1
        chosen = [option for option in ("history", "final_profile",
                                        "final_model")
                  if mesa_inlists.get("{}_{}".format(option, star))]
        setattr(output, star, tuple(chosen))
    ignored = [key for key in mesa_inlists
               if INLIST_KEY.match(key) is None and key not in COLUMN_KEYS
               and key not in OUTPUT_KEYS and key not in KNOWN_KEYS
               and key not in ("scenario", "single_star_grid",
                               "zams_filename", "zams_filename_1",
                               "zams_filename_2")]
    if ignored:
        config.notes.append("[mesa_inlists] keys not used: {}".format(
            ", ".join(ignored)))

    # [run_parameters] -> grid
    if "grid" not in run_parameters:
        raise ConfigError("run_parameters", "grid", "missing")
    config.grid.points = Path(run_parameters["grid"])
    config.grid.mode = grid_type
    if run_parameters.get("psycris_inifile") is not None:
        config.grid.psycris_inifile = Path(run_parameters["psycris_inifile"])
    config.grid.keep_profiles = bool(run_parameters.get("keep_profiles",
                                                        False))
    config.grid.keep_photos = bool(run_parameters.get("keep_photos", False))
    if run_directory is not None:
        config.grid.run_directory = Path(run_directory)

    # [mesa_extras] -> extras, in the order of the file
    config.extras = {key: Path(value) for key, value in mesa_extras.items()
                     if key != "mesa_dir" and value is not None}

    # [slurm] -> submission
    if submission_type == "shell":
        backend = "shell"
    else:
        backend = "slurm_array" if slurm.get("job_array") else "slurm_mpi"
    defaults = SubmissionSection()
    config.submission = SubmissionSection(
        backend=backend,
        partition=slurm.get("partition"),
        account=slurm.get("account"),
        walltime=slurm.get("walltime"),
        email=slurm.get("email"),
        nodes=slurm.get("number_of_nodes", defaults.nodes),
        mpi_tasks=slurm.get("number_of_mpi_tasks", defaults.mpi_tasks),
        cpus_per_task=slurm.get("number_of_cpus_per_task",
                                defaults.cpus_per_task),
        work_dir=slurm.get("work_dir"),
        newgroup=slurm.get("newgroup"))
    return config
