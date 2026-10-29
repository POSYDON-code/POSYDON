"""The settings of a grid, read from an ini file and checked.

Everything ``posydon-setup-grid`` needs is collected in one ``GridConfig``.
Each part of it comes from one section of the ini file::

    [grid]            -> GridSection        what to run, where
    [scenario]        -> ScenarioSection    which POSYDON inlists to use
    [mesa_inlists]    -> InlistFiles        \\ files per layer
    [user_inlists]    -> InlistFiles        /  and target
    [output]          -> OutputSection      history, profiles, final models
    [submission]      -> SubmissionSection  shell or SLURM settings
    [extras]          -> dict               extra Fortran sources by name

The ``posydon_inlists`` layer between them is filled by the scenario (the
release's files); a new-format ini has no section for it.

``read_config`` reads this new format. Today's ini files (with
``[mesa_inlists]`` holding ``*_mesa_defaults`` keys) are read by
``posydon.grids.setup.legacy`` into the same ``GridConfig``.

Every problem is reported as a ``ConfigError`` that names the section and
key, e.g. ``[grid] type: unknown grid type 'MS-MS'; expected one of ...``.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import configparser
import re
from dataclasses import dataclass, field
from pathlib import Path

GRID_MODES = ("fixed", "dynamic")
# [scenario] keys that used to make the setup check out the inlist
# repository with git; it now reads inlist_repo as it is
GIT_KEYS = ("source", "ref", "commit", "branch")
BACKENDS = ("shell", "slurm_array", "slurm_mpi")

# targets one inlist file can feed, and the shorthands for a whole file
INLIST_TARGETS = ("binary_controls", "binary_job", "star1_controls",
                  "star1_job", "star2_controls", "star2_job")
INLIST_SHORTHANDS = {"binary": ("binary_controls", "binary_job"),
                     "star1": ("star1_controls", "star1_job"),
                     "star2": ("star2_controls", "star2_job")}
# files that are used whole instead of merged: column lists, and Fortran
# sources (each source has one place: [user_inlists]; MESA's own come from
# $MESA_DIR, and the release replaces them)
COLUMN_FILES = ("history_columns", "binary_history_columns",
                "profile_columns")
SOURCE_FILES = ("run_binary_extras", "run_star_extras", "run_star1_extras",
                "run_star2_extras")
WHOLE_FILES = COLUMN_FILES + SOURCE_FILES
# what [extras] of a new-format ini takes: MESA's work files
EXTRAS_KEYS = ("makefile_binary", "makefile_star", "binary_run", "star_run")
OUTPUT_OPTIONS = ("history", "final_profile", "final_model")


class SetupError(ValueError):
    """A mistake the user can fix: in the ini, a file or the environment.

    ``posydon-grid`` prints it as one line, ``posydon-grid: error: ...``,
    and exits with status 2; any other exception is a bug and keeps its
    traceback.
    """


class ConfigError(SetupError):
    """A problem in the ini file, naming the section and key."""

    def __init__(self, section, key, message):
        self.section = section
        self.key = key
        where = "[{}]".format(section) if key is None \
            else "[{}] {}".format(section, key)
        super().__init__("{}: {}".format(where, message))


@dataclass
class GridSection:
    """What to run and where."""

    type: str = None
    mode: str = "fixed"
    points: Path = None
    run_directory: Path = None
    zams_star1: Path = None
    zams_star2: Path = None
    psycris_inifile: Path = None
    keep_profiles: bool = False
    keep_photos: bool = False


@dataclass
class ScenarioSection:
    """Which POSYDON inlists to use; empty when no scenario is set.

    A new-format ini sets ``release`` (+ ``rerun``, ``mesa``) and
    ``inlist_repo``; the setup reads that folder as it is. ``source``,
    ``branch``, ``commit`` and ``system_type`` come from a legacy ini's
    ``scenario = [source, 'branch-commit', system type]``, whose
    repository is still checked out at the commit.
    """

    release: str = None
    rerun: str = None
    mesa: str = None
    inlist_repo: Path = None
    source: str = "posydon"
    branch: str = None
    commit: str = None
    system_type: str = None

    @property
    def is_set(self):
        """Whether the ini names a scenario at all."""
        return any(value is not None for value in
                   (self.release, self.commit, self.branch,
                    self.system_type))


@dataclass
class InlistFiles:
    """The files of one layer: per target a list, plus whole files.

    ``targets`` maps a target (``star1_controls``, ...,
    ``star1_formation_step0_controls``) to its files, applied in order.
    ``whole_files`` maps a name from ``WHOLE_FILES`` to one file.
    """

    targets: dict = field(default_factory=dict)
    whole_files: dict = field(default_factory=dict)

    def add(self, target, path):
        """Append a file to a target."""
        self.targets.setdefault(target, []).append(path)


@dataclass
class OutputSection:
    """What MESA writes: history, final profiles and final models."""

    history_interval: int = 1
    binary_history: bool = True
    star1: tuple = ("history", "final_model")
    star2: tuple = ("history", "final_model")


@dataclass
class SubmissionSection:
    """How the grid is submitted.

    ``job_array`` only matters for the shell backend: run the points one
    after another (True) or as one ``mpirun`` (False). The SLURM backends
    say it in their name.
    """

    backend: str = "shell"
    job_array: bool = True
    partition: str = None
    account: str = None
    walltime: str = None
    email: str = None
    nodes: int = 1
    mpi_tasks: int = 1
    cpus_per_task: int = 1
    mem_per_cpu: str = "4G"
    work_dir: str = None
    newgroup: str = None


@dataclass
class GridConfig:
    """All settings of a grid."""

    grid: GridSection = field(default_factory=GridSection)
    scenario: ScenarioSection = field(default_factory=ScenarioSection)
    inlists: dict = field(default_factory=lambda: {
        layer: InlistFiles() for layer in
        ("mesa_inlists", "posydon_inlists", "user_inlists")})
    output: OutputSection = field(default_factory=OutputSection)
    submission: SubmissionSection = field(default_factory=SubmissionSection)
    extras: dict = field(default_factory=dict)
    # where the settings came from ("new" or "legacy") and what to tell the
    # user about them, e.g. deprecated keys
    source_format: str = "new"
    notes: list = field(default_factory=list)
    # the extras keys defaults.fill_defaults filled in with MESA's stock
    # files, because the ini does not name them
    default_extras: tuple = ()
    # a legacy ini's sections as read by posydon.utils.configfile, kept
    # because their key order decides today's layering
    legacy: dict = field(default_factory=dict)


###############################################################################
# reading values
###############################################################################
def _choice(section, key, value, allowed):
    if value not in allowed:
        raise ConfigError(section, key, "unknown value {!r}; expected one of "
                          "{}".format(value, ", ".join(allowed)))
    return value


def _bool(section, key, value):
    text = str(value).strip().lower()
    if text in ("true", "yes", "1", "on"):
        return True
    if text in ("false", "no", "0", "off"):
        return False
    raise ConfigError(section, key,
                      "expected true or false, got {!r}".format(value))


def _int(section, key, value):
    try:
        return int(value)
    except ValueError:
        raise ConfigError(section, key,
                          "expected an integer, got {!r}".format(value))


def _list(value):
    return [item.strip() for item in value.split(",") if item.strip()]


def _path(base, value):
    """A path from the ini, relative to the ini's folder."""
    path = Path(value).expanduser()
    return path if path.is_absolute() else base / path


###############################################################################
# one parser per section
###############################################################################
def _read_grid(options, base):
    grid = GridSection()
    for key, value in options.items():
        if key == "type":
            # the names come from the registry; imported here because
            # grid_types.py imports ConfigError from this module
            from posydon.grids.setup.grid_types import get_grid_type
            grid.type = get_grid_type(value).name
        elif key == "mode":
            grid.mode = _choice("grid", key, value, GRID_MODES)
        elif key in ("points", "run_directory", "psycris_inifile"):
            setattr(grid, key, _path(base, value))
        elif key == "zams":
            grid.zams_star1 = grid.zams_star2 = _path(base, value)
        elif key in ("zams_star1", "zams_star2"):
            setattr(grid, key, _path(base, value))
        elif key in ("keep_profiles", "keep_photos"):
            setattr(grid, key, _bool("grid", key, value))
        else:
            raise ConfigError("grid", key, "unknown key")
    for key in ("type", "points"):
        if getattr(grid, key) is None:
            raise ConfigError("grid", key, "missing")
    if grid.run_directory is None:
        grid.run_directory = base
    if grid.mode == "dynamic" and grid.psycris_inifile is None:
        raise ConfigError("grid", "psycris_inifile",
                          "needed for mode = dynamic")
    return grid


def _read_scenario(options, base):
    scenario = ScenarioSection()
    for key, value in options.items():
        if key in GIT_KEYS:
            # the release names the version; keeping the clone up to date
            # is done outside the setup
            raise ConfigError("scenario", key, "no longer used: the setup "
                              "reads inlist_repo as it is, and the release "
                              "names the version. Remove this key, and "
                              "update the clone yourself (git pull)")
        elif key in ("release", "rerun", "mesa"):
            setattr(scenario, key, value)
        elif key == "inlist_repo":
            scenario.inlist_repo = _path(base, value)
        else:
            raise ConfigError("scenario", key, "unknown key")
    if scenario.rerun is not None and scenario.release is None:
        raise ConfigError("scenario", "rerun", "a rerun needs a release")
    if scenario.release is None and (scenario.inlist_repo is not None or
                                     scenario.mesa is not None):
        raise ConfigError("scenario", "release", "missing; a scenario names "
                          "a release of the inlist repository")
    return scenario


def _inlist_entries(key, paths):
    """The (target, path) pairs one key feeds, or None for an unknown key."""
    if key in INLIST_TARGETS:
        return [(key, path) for path in paths]
    if key in INLIST_SHORTHANDS:
        return [(target, path) for target in INLIST_SHORTHANDS[key]
                for path in paths]
    if key in ("star1_formation", "star2_formation"):
        # file i is formation step i
        steps = [("{}_step{}".format(key, index), path)
                 for index, path in enumerate(paths)]
    elif re.fullmatch(r"star[12]_formation_step\d+", key):
        # extra files on top of one formation step
        steps = [(key, path) for path in paths]
    else:
        return None
    # a formation step file feeds both its controls and its job
    return [("{}_{}".format(step, part), path)
            for step, path in steps for part in ("controls", "job")]


def _read_inlist_files(section, options, base):
    files = InlistFiles()
    for key, value in options.items():
        paths = [_path(base, item) for item in _list(value)]
        if not paths:
            raise ConfigError(section, key, "no file given")
        if key in WHOLE_FILES:
            if key in SOURCE_FILES and section != "user_inlists":
                # a file here could never win: MESA's stock sources are
                # filled in, and the release replaces them
                raise ConfigError(section, key, "Fortran sources go in "
                                  "[user_inlists]; MESA's own come from "
                                  "$MESA_DIR")
            if len(paths) > 1:
                raise ConfigError(section, key, "expected one file")
            files.whole_files[key] = paths[0]
            continue
        entries = _inlist_entries(key, paths)
        if entries is None:
            raise ConfigError(section, key, "unknown key")
        for target, path in entries:
            files.add(target, path)
    return files


def _read_extras(options, base):
    extras = {}
    for key, value in options.items():
        if key not in EXTRAS_KEYS:
            raise ConfigError("extras", key, "not an [extras] key: [extras] "
                              "takes MESA's work files ({}); Fortran sources "
                              "go in [user_inlists] (run_star_extras, "
                              "run_binary_extras, run_star1_extras, "
                              "run_star2_extras)".format(
                                  ", ".join(EXTRAS_KEYS)))
        extras[key] = _path(base, value)
    return extras


def _read_output(options):
    output = OutputSection()
    for key, value in options.items():
        if key == "history_interval":
            output.history_interval = _int("output", key, value)
        elif key in ("star1", "star2"):
            chosen = _list(value)
            if chosen == ["none"]:
                chosen = []
            for option in chosen:
                _choice("output", key, option, OUTPUT_OPTIONS)
            setattr(output, key, tuple(chosen))
        elif key == "binary":
            chosen = _list(value)
            if chosen not in (["history"], ["none"]):
                raise ConfigError("output", key,
                                  "expected history or none")
            output.binary_history = chosen == ["history"]
        else:
            raise ConfigError("output", key, "unknown key")
    return output


def _read_submission(options):
    submission = SubmissionSection()
    for key, value in options.items():
        if key == "backend":
            submission.backend = _choice("submission", key, value, BACKENDS)
        elif key == "job_array":
            submission.job_array = _bool("submission", key, value)
        elif key in ("nodes", "mpi_tasks", "cpus_per_task"):
            setattr(submission, key, _int("submission", key, value))
        elif key in ("partition", "account", "walltime", "email",
                     "mem_per_cpu", "work_dir", "newgroup"):
            setattr(submission, key, value)
        else:
            raise ConfigError("submission", key, "unknown key")
    return submission


###############################################################################
# the whole file
###############################################################################
SECTIONS = ("grid", "scenario", "mesa_inlists", "user_inlists", "output",
            "submission", "extras")
# the inlist layers an ini can fill; posydon_inlists comes from [scenario]
INI_LAYERS = ("mesa_inlists", "user_inlists")


def read_config(path):
    """Read and check a new-format grid ini file.

    Parameters
    ----------
    path : str or Path
        The ini file. Relative paths inside it are relative to its folder,
        and ``${section:key}`` interpolation works as in configparser.

    Returns
    -------
    GridConfig
        The checked settings.

    Raises
    ------
    ConfigError
        For a missing file, an unknown section or key, a missing required
        key or an invalid value.
    """
    # parse the file
    path = Path(path)
    if not path.is_file():
        raise ConfigError("file", None, "{} not found".format(path))
    parser = configparser.ConfigParser(
        interpolation=configparser.ExtendedInterpolation(),
        inline_comment_prefixes=(";",))  # "key = value  ; comment"
    parser.optionxform = str  # keep the case of keys
    parser.read(path)

    # check the sections before reading any of them
    if parser.has_section("posydon_inlists"):
        raise ConfigError("posydon_inlists", None, "no longer used: "
                          "POSYDON's inlists come from [scenario] (the "
                          "release); put your own files in [user_inlists]")
    unknown = [name for name in parser.sections() if name not in SECTIONS]
    if unknown:
        raise ConfigError(unknown[0], None, "unknown section; expected "
                          "{}".format(", ".join(SECTIONS)))
    if not parser.has_section("grid"):
        raise ConfigError("grid", None, "missing")

    # read each section; a missing optional section gives the defaults
    def options(name):
        return dict(parser.items(name)) if parser.has_section(name) else {}

    # relative paths in the ini are relative to its folder
    base = path.resolve().parent
    config = GridConfig(
        grid=_read_grid(options("grid"), base),
        scenario=_read_scenario(options("scenario"), base),
        output=_read_output(options("output")),
        submission=_read_submission(options("submission")),
        extras=_read_extras(options("extras"), base))
    for layer in INI_LAYERS:
        config.inlists[layer] = _read_inlist_files(layer, options(layer),
                                                   base)
    return config
