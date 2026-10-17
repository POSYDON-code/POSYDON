"""The settings of a grid, read from an ini file and checked.

Everything ``posydon-setup-grid`` needs is collected in one ``GridConfig``.
Each part of it comes from one section of the ini file::

    [grid]            -> GridSection        what to run, where
    [scenario]        -> ScenarioSection    which posydon_inlists to use
    [mesa_inlists]    -> InlistFiles        \\
    [posydon_inlists] -> InlistFiles         > files per layer and target
    [user_inlists]    -> InlistFiles        /
    [output]          -> OutputSection      history, profiles, final models
    [submission]      -> SubmissionSection  shell or SLURM settings
    [extras]          -> dict               extra Fortran sources by name

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

GRID_TYPES = ("HMS-HMS", "CO-HMS_RLO", "CO-HeMS", "single_HMS",
              "single_HeMS", "HeMS-HMS")
GRID_MODES = ("fixed", "dynamic")
SCENARIO_SOURCES = ("posydon", "user", "path")
BACKENDS = ("shell", "slurm_array", "slurm_mpi")

# targets one inlist file can feed, and the shorthands for a whole file
INLIST_TARGETS = ("binary_controls", "binary_job", "star1_controls",
                  "star1_job", "star2_controls", "star2_job")
INLIST_SHORTHANDS = {"binary": ("binary_controls", "binary_job"),
                     "star1": ("star1_controls", "star1_job"),
                     "star2": ("star2_controls", "star2_job")}
# files that are used whole instead of merged
WHOLE_FILES = ("history_columns", "binary_history_columns",
               "profile_columns", "run_binary_extras", "run_star_extras")
OUTPUT_OPTIONS = ("history", "final_profile", "final_model")


class ConfigError(ValueError):
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
    """Which ``posydon_inlists`` to use; empty when no scenario is set."""

    source: str = "posydon"
    release: str = None
    rerun: str = None
    mesa: str = None
    ref: str = None
    commit: str = None
    branch: str = None
    system_type: str = None
    inlist_repo: Path = None

    @property
    def is_set(self):
        """Whether the ini names a scenario at all."""
        return any(value is not None for value in
                   (self.release, self.ref, self.commit, self.branch,
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
            grid.type = _choice("grid", key, value, GRID_TYPES)
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
        if key == "source":
            scenario.source = _choice("scenario", key, value,
                                      SCENARIO_SOURCES)
        elif key == "commit":
            if not re.fullmatch(r"[0-9a-f]{7,40}", value):
                raise ConfigError("scenario", key, "expected a git commit "
                                  "hash (7 to 40 hex digits), got "
                                  "{!r}".format(value))
            scenario.commit = value
        elif key in ("release", "rerun", "mesa", "ref", "branch"):
            setattr(scenario, key, value)
        elif key == "inlist_repo":
            scenario.inlist_repo = _path(base, value)
        else:
            raise ConfigError("scenario", key, "unknown key")
    if scenario.rerun is not None and scenario.release is None:
        raise ConfigError("scenario", "rerun", "a rerun needs a release")
    if scenario.release is not None and scenario.source != "path" and \
            scenario.ref is None and scenario.commit is None:
        # no hidden default: the ini says which version of the inlist
        # repository it was set up with
        raise ConfigError("scenario", "release", "a release needs ref (a "
                          "branch or tag of the inlist repository) or "
                          "commit")
    return scenario


def _read_inlist_files(section, options, base):
    files = InlistFiles()
    for key, value in options.items():
        paths = [_path(base, item) for item in _list(value)]
        if not paths:
            raise ConfigError(section, key, "no file given")
        step = re.fullmatch(r"(star[12])_formation_step(\d+)", key)
        if key in INLIST_TARGETS:
            for path in paths:
                files.add(key, path)
        elif key in INLIST_SHORTHANDS:
            for target in INLIST_SHORTHANDS[key]:
                for path in paths:
                    files.add(target, path)
        elif key in ("star1_formation", "star2_formation"):
            # one file per formation step
            for index, path in enumerate(paths):
                for part in ("controls", "job"):
                    files.add("{}_step{}_{}".format(key, index, part), path)
        elif step is not None:
            # extra files on top of one formation step
            star, index = step.groups()
            for path in paths:
                for part in ("controls", "job"):
                    files.add("{}_formation_step{}_{}".format(
                        star, index, part), path)
        elif key in WHOLE_FILES:
            if len(paths) > 1:
                raise ConfigError(section, key, "expected one file")
            files.whole_files[key] = paths[0]
        else:
            raise ConfigError(section, key, "unknown key")
    return files


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
SECTIONS = ("grid", "scenario", "mesa_inlists", "posydon_inlists",
            "user_inlists", "output", "submission", "extras")


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
    path = Path(path)
    if not path.is_file():
        raise ConfigError("file", None, "{} not found".format(path))
    parser = configparser.ConfigParser(
        interpolation=configparser.ExtendedInterpolation(),
        inline_comment_prefixes=(";",))  # "key = value  ; comment"
    parser.optionxform = str  # keep the case of keys
    parser.read(path)
    unknown = [name for name in parser.sections() if name not in SECTIONS]
    if unknown:
        raise ConfigError(unknown[0], None, "unknown section; expected "
                          "{}".format(", ".join(SECTIONS)))
    if not parser.has_section("grid"):
        raise ConfigError("grid", None, "missing")

    def options(name):
        return dict(parser.items(name)) if parser.has_section(name) else {}

    base = path.resolve().parent
    config = GridConfig(
        grid=_read_grid(options("grid"), base),
        scenario=_read_scenario(options("scenario"), base),
        output=_read_output(options("output")),
        submission=_read_submission(options("submission")),
        extras={key: _path(base, value)
                for key, value in options("extras").items()})
    for layer in config.inlists:
        config.inlists[layer] = _read_inlist_files(layer, options(layer),
                                                   base)
    return config
