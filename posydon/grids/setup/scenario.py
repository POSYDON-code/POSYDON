"""Get the inlist repository at the scenario's commit and list its files.

A scenario points at one commit of POSYDON-MESA-INLISTS. Setting it up has
two parts:

1. ``fetch_inlist_repo``: make sure a local clone exists and has that
   commit checked out. This runs the same git commands as today's
   ``posydon-setup-grid`` (clone if missing, check out the branch, pull,
   check out the commit), but with ``git -C <clone>`` instead of changing
   the working directory, and it stops with an error if the commit cannot
   be checked out instead of carrying on with whatever is there.

2. ``scenario_files``: list the files the grid type takes from that
   commit, as ``posydon_inlists``. For every grid type::

       default_common_inlists/binary/inlist_project -> binary_controls, _job
       default_common_inlists/binary/inlist1        -> star1_controls, _job
       default_common_inlists/binary/inlist2        -> star2_controls, _job
       default_common_inlists/*.list                -> column lists
       default_common_inlists/binary/src/*.f        -> extras

   and on top, from the grid type's own folder (``GridType.folder``)::

       <folder>/binary/inlist_project, inlist1, inlist2
       <folder>/*.list, <folder>/src/*.f
       <folder>/star{1,2}_formation/*step*          -> formation steps

   HeMS-HMS only takes ``binary/inlist1`` and ``star1_formation``. The
   single-star grids add the single-star flag (see ``grid_types``).

The result is data (``ScenarioFiles``); nothing is read or merged here.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from posydon.grids.setup.config import ConfigError, InlistFiles
from posydon.grids.setup.inlist import read_inlist

REPOSITORIES = {
    "posydon": ("https://github.com/POSYDON-code/POSYDON-MESA-INLISTS.git",
                ".posydon_mesa_inlists"),
    "user": ("https://github.com/POSYDON-code/USER-MESA-INLISTS.git",
             ".user_mesa_inlists"),
}
MESA_VERSION = "r11701"
COMMON = "default_common_inlists"
COLUMN_LISTS = {"history_columns": "history_columns.list",
                "binary_history_columns": "binary_history_columns.list",
                "profile_columns": "profile_columns.list"}
# the single-star flag that run_star_extras.f reads
SINGLE_STAR_FLAG = {"x_logical_ctrl(1)": ".true."}


###############################################################################
# 1. the local clone
###############################################################################
def inlist_repo_path(scenario, home=None):
    """Where the clone of the scenario's repository lives.

    ``[scenario] inlist_repo`` if set, otherwise ``~/.posydon_mesa_inlists``
    (or ``~/.user_mesa_inlists`` for ``source = user``), as today.
    """
    if scenario.inlist_repo is not None:
        return Path(scenario.inlist_repo)
    if scenario.source not in REPOSITORIES:
        raise ConfigError("scenario", "source", "no repository is known "
                          "for source {!r}; set inlist_repo".format(
                              scenario.source))
    home = Path.home() if home is None else Path(home)
    return home / REPOSITORIES[scenario.source][1]


def _git(repo, *args):
    """Run git in ``repo``; return the CompletedProcess (never raises)."""
    return subprocess.run(["git", "-C", str(repo), *args],
                          capture_output=True, text=True)


def fetch_inlist_repo(scenario, home=None):
    """Clone or update the inlist repository and check out the commit.

    Parameters
    ----------
    scenario : ScenarioSection
        Needs ``commit``; ``branch`` is checked out and pulled first, as
        today, so that a commit newer than the clone is found.
    home : str, Path or None
        Folder for the default clone location; defaults to ``~``.

    Returns
    -------
    Path
        The clone, with ``scenario.commit`` checked out.

    Raises
    ------
    ConfigError
        If the clone cannot be made or the commit cannot be checked out.
    """
    repo = inlist_repo_path(scenario, home=home)
    if not repo.is_dir():
        url = REPOSITORIES.get(scenario.source, (None,))[0]
        result = subprocess.run(["git", "clone", url, str(repo)],
                                capture_output=True, text=True)
        if result.returncode != 0:
            raise ConfigError("scenario", "source", "could not clone {} "
                              "into {}: {}".format(url, repo,
                                                   result.stderr.strip()))
    if scenario.branch is not None:
        # failures here are not fatal: the branch may only exist remotely,
        # and the pull needs a network
        _git(repo, "checkout", scenario.branch)
        _git(repo, "pull")
    result = _git(repo, "checkout", scenario.commit)
    if result.returncode != 0:
        raise ConfigError("scenario", "commit", "could not check out {} in "
                          "{}: {}".format(scenario.commit, repo,
                                          result.stderr.strip()))
    return repo


###############################################################################
# 2. the files a grid type takes from the repository
###############################################################################
@dataclass
class ScenarioFiles:
    """What the scenario contributes, as the posydon_inlists layer.

    ``inlists`` holds the files per target and the column lists.
    ``flags`` holds values that are not in any file (the single-star flag),
    per target. ``extras`` maps the legacy ``[mesa_extras]`` key names to the
    Fortran sources the scenario sets. ``zams_star1`` is a ZAMS model named
    in the common ``inlist1`` (``None`` if it names none).
    ``formation_stars`` lists the stars the scenario builds from formation
    steps.
    """

    inlists: InlistFiles = field(default_factory=InlistFiles)
    flags: dict = field(default_factory=dict)
    extras: dict = field(default_factory=dict)
    zams_star1: Path = None
    zams_star2: Path = None
    formation_stars: tuple = ()


def _step_order(path):
    """Sort key putting ``inlist_step2`` before ``inlist_step10``."""
    return [int(part) if part.isdigit() else part
            for part in re.split(r"(\d+)", path.name)]


def _formation_steps(folder, star):
    """The formation step files of a star, in order (``*step*``).

    Numbers in the names are compared as numbers (before, a plain sort put
    ``inlist_step10`` before ``inlist_step2``).
    """
    return sorted((folder / "{}_formation".format(star)).glob("*step*"),
                  key=_step_order)


def _zams_from_inlist(root, inlist):
    """The ZAMS model named in an inlist, looked up in ZAMS_models/."""
    if not inlist.is_file():
        return None
    for values in read_inlist(inlist).values():
        name = values.get("zams_filename")
        if name is not None:
            path = root / "ZAMS_models" / Path(name.strip("'\"")).name
            return path if path.is_file() else None
    return None


def _add_star(files, star, inlist):
    for part in ("controls", "job"):
        files.add("{}_{}".format(star, part), inlist)


def _add_every_step(files, star, inlist):
    for part in ("controls", "job"):
        files.add("{}_formation_every_step_{}".format(star, part), inlist)


def _add_formation(result, star, steps, common_inlist1):
    """Formation steps: the common inlist1 under every step, then each."""
    _add_every_step(result.inlists, star, common_inlist1)
    for index, step in enumerate(steps):
        for part in ("controls", "job"):
            result.inlists.add("{}_formation_step{}_{}".format(
                star, index, part), step)
    result.formation_stars += (star,)


def scenario_files(repo, grid_type):
    """List the files a grid type takes from a checked-out repository.

    Parameters
    ----------
    repo : str or Path
        The clone, at the scenario's commit.
    grid_type : GridType
        From ``posydon.grids.setup.grid_types``.

    Returns
    -------
    ScenarioFiles
    """
    root = Path(repo) / MESA_VERSION
    common = root / COMMON
    result = ScenarioFiles()
    files = result.inlists

    # every grid type: default_common_inlists. A single-star grid evolves
    # star 1 alone, through its formation step(s): its inlist1 feeds those
    # steps, and the binary and star 2 take no files.
    common_inlist1 = common / "binary" / "inlist1"
    if not grid_type.single_star:
        files.add("binary_controls", common / "binary" / "inlist_project")
        files.add("binary_job", common / "binary" / "inlist_project")
        _add_star(files, "star1", common_inlist1)
        _add_star(files, "star2", common / "binary" / "inlist2")
    for name, file_name in COLUMN_LISTS.items():
        files.whole_files[name] = common / file_name
    src = common / "binary" / "src"
    result.extras = {"posydon_binary_extras": src / "run_binary_extras.f",
                     "posydon_star_binary_extras": src / "run_star_extras.f",
                     "mesa_star1_extras": src / "run_star_extras.f"}
    result.zams_star1 = _zams_from_inlist(root, common / "binary" /
                                          "inlist1")
    result.zams_star2 = _zams_from_inlist(root, common / "binary" /
                                          "inlist2")

    # the grid type's own folder
    folder = root / grid_type.folder if grid_type.folder else None
    if grid_type.single_star:
        steps = _formation_steps(folder, "star1") if folder else []
        if steps:
            _add_formation(result, "star1", steps, common_inlist1)
            # the flag as an extra, last step
            last = len(steps)
            result.flags["star1_formation_step{}_controls".format(last)] = \
                dict(SINGLE_STAR_FLAG)
            result.flags["star1_formation_step{}_job".format(last)] = {}
        else:
            # one step: the common inlist1 and the flag
            _add_every_step(files, "star1", common_inlist1)
            result.flags["star1_formation_every_step_controls"] = \
                dict(SINGLE_STAR_FLAG)
    elif folder is not None and grid_type.name == "HeMS-HMS":
        if (folder / "binary" / "inlist1").is_file():
            _add_star(files, "star1", folder / "binary" / "inlist1")
        steps = _formation_steps(folder, "star1")
        if steps:
            _add_formation(result, "star1", steps, common_inlist1)
    elif folder is not None:
        if (folder / "binary" / "inlist_project").is_file():
            files.add("binary_controls", folder / "binary" / "inlist_project")
            files.add("binary_job", folder / "binary" / "inlist_project")
        for star, name in (("star1", "inlist1"), ("star2", "inlist2")):
            if (folder / "binary" / name).is_file():
                _add_star(files, star, folder / "binary" / name)
        for name, file_name in COLUMN_LISTS.items():
            if (folder / file_name).is_file():
                files.whole_files[name] = folder / file_name
        for key, file_name in (("user_binary_extras", "run_binary_extras.f"),
                               ("user_star_binary_extras",
                                "run_star_extras.f")):
            if (folder / "src" / file_name).is_file():
                result.extras[key] = folder / "src" / file_name
        for star in ("star1", "star2"):
            steps = _formation_steps(folder, star)
            if steps:
                _add_formation(result, star, steps, common_inlist1)
    return result
