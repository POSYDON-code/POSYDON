"""Turn the ini's ``[scenario]`` into the files of ``posydon_inlists``.

There are two kinds of scenario:

- **release** (new-format ini): ``release`` (+ ``rerun``) of the
  restructured inlist repository, with ``[grid] type``. The setup reads the
  folder ``inlist_repo`` as it is: it never runs git to change it. Which
  version of the repository that folder holds is the user's business
  (``git pull``, a branch, a fork); the release names the physics, and the
  files are read from its manifests (``manifest.release_files``).
- **legacy**: a commit of the inlist repository and a system type, as in
  today's ini files (``scenario = ['posydon', '<branch>-<commit>',
  '<system type>']``). The repository is checked out at the commit and the
  files are found by today's rules (``scenario.scenario_files``).

``resolve_scenario`` returns the files with where they came from, for the
run directory's ``setup_record.ini``: the folder, and, read without
changing anything, its commit, whether the files used differ from that
commit (``dirty``), its ``origin`` and a fingerprint of the files used.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import hashlib
from dataclasses import dataclass
from pathlib import Path

from posydon.grids.setup.config import ConfigError
from posydon.grids.setup.grid_types import get_grid_type
from posydon.grids.setup.manifest import (
    MESA_DEFAULTS,
    find_release,
    release_files,
)
from posydon.grids.setup.scenario import (
    REPOSITORIES,
    ScenarioFiles,
    fetch_inlist_repo,
    git,
    inlist_repo_path,
    scenario_files,
)


@dataclass
class ResolvedScenario:
    """The files of the scenario and where they came from.

    ``commit``, ``dirty`` and ``remote`` are None when ``repo`` is not a
    git repository of its own (or git is not installed). ``fingerprint`` is
    the sha256 of the files used, ``n_files`` how many there are.
    """

    files: ScenarioFiles
    repo: Path
    commit: str = None
    release: object = None
    dirty: bool = None
    remote: str = None
    fingerprint: str = None
    n_files: int = 0


###############################################################################
# where the files came from, read without changing anything
###############################################################################
def files_used(files, release=None, grid_type=None, rerun=None):
    """Every file a scenario reads, sorted.

    The scenario's inlists, column lists, Fortran sources and ZAMS models;
    for a release also its ``release.ini``, the manifests and the MESA
    defaults of its installation (those that exist).
    """
    paths = {Path(path) for targets in files.inlists.targets.values()
             for path in targets}
    paths |= {Path(path) for path in files.inlists.whole_files.values()}
    paths |= {Path(path) for path in files.extras.values()}
    paths |= {Path(path) for path in (files.zams_star1, files.zams_star2)
              if path is not None}
    if release is not None:
        paths.add(release.folder / "release.ini")
        paths.add(release.folder / grid_type / "scenario.ini")
        if rerun is not None:
            paths.add(release.folder / "reruns" / rerun / "scenario.ini")
        defaults = release.folder.parents[1] / "mesa" / release.mesa / \
            "defaults"
        paths |= {defaults / rel for rel in MESA_DEFAULTS.values()
                  if (defaults / rel).is_file()}
    return sorted(paths)


def _relative(repo, path):
    try:
        return Path(path).resolve().relative_to(Path(repo).resolve())
    except ValueError:
        return None


def fingerprint(repo, paths):
    """sha256 over the files' paths (relative to ``repo``) and contents.

    Two setups with the same fingerprint read identical files, from
    whatever folder, branch or machine. A file that does not exist counts
    as ``missing`` (a legacy scenario lists files that may not be there).
    """
    lines = []
    for path in paths:
        rel = _relative(repo, path)
        name = rel.as_posix() if rel is not None else str(path)
        content = hashlib.sha256(Path(path).read_bytes()).hexdigest() \
            if Path(path).is_file() else "missing"
        lines.append("{}\0{}\n".format(name, content))
    return hashlib.sha256("".join(sorted(lines)).encode()).hexdigest()


def git_state(repo, paths, untracked=True):
    """``(commit, dirty, remote)`` of a folder, read-only.

    All None unless ``repo`` is the top of a git repository (a folder
    inside some other repository would report that one's commit). ``dirty``
    says whether any of ``paths`` differs from the commit, untracked files
    included unless ``untracked`` is False.
    """
    try:
        top = git(repo, "rev-parse", "--show-toplevel")
    except OSError:  # no git
        return None, None, None
    if top.returncode != 0 or \
            Path(top.stdout.strip()).resolve() != Path(repo).resolve():
        return None, None, None
    commit = git(repo, "rev-parse", "HEAD").stdout.strip() or None
    rels = [rel.as_posix() for rel in (_relative(repo, path)
                                       for path in paths)
            if rel is not None]
    status = git(repo, "status", "--porcelain", "--untracked-files=" +
                 ("all" if untracked else "no"), "--", *rels)
    dirty = bool(status.stdout.strip()) if status.returncode == 0 else None
    remote = git(repo, "remote", "get-url", "origin")
    return commit, dirty, \
        remote.stdout.strip() if remote.returncode == 0 else None


def _resolved(files, repo, release=None, grid_type=None, rerun=None):
    paths = files_used(files, release, grid_type, rerun)
    commit, dirty, remote = git_state(repo, paths)
    return ResolvedScenario(files=files, repo=repo, commit=commit,
                            release=release, dirty=dirty, remote=remote,
                            fingerprint=fingerprint(repo, paths),
                            n_files=len(paths))


###############################################################################
# the scenario
###############################################################################
def release_repo(scenario, home=None):
    """The folder of the inlist repository, as it is.

    Parameters
    ----------
    scenario : ScenarioSection
        ``inlist_repo``, or the default ``~/.posydon_mesa_inlists``.
    home : str, Path or None
        Folder for the default location; defaults to ``~``.

    Raises
    ------
    ConfigError
        If the folder does not exist; the message says how to clone it.
    """
    repo = inlist_repo_path(scenario, home=home)
    if not repo.is_dir():
        raise ConfigError("scenario", "inlist_repo", "{0} does not exist. "
                          "The setup does not clone or update the inlist "
                          "repository; clone it there yourself:\n"
                          "  git clone --branch v2_sources {1} {0}".format(
                              repo, REPOSITORIES["posydon"][0]))
    return repo


def resolve_scenario(config, home=None):
    """The scenario's files, or ``None`` if the ini sets no scenario.

    Parameters
    ----------
    config : GridConfig
        ``config.grid.type`` names the grid type.
    home : str, Path or None
        Folder for the default inlist repository; defaults to ``~``.

    Returns
    -------
    ResolvedScenario or None
    """
    scenario = config.scenario
    if not scenario.is_set:
        return None
    if config.grid.type is None:
        raise ConfigError("grid", "type", "needed to use a scenario")
    grid_type = get_grid_type(config.grid.type)
    if scenario.release is None:
        if scenario.commit is None:
            raise ConfigError("scenario", "release", "missing; without a "
                              "release, a (legacy) scenario needs a commit")
        repo = fetch_inlist_repo(scenario, home=home)
        return _resolved(scenario_files(repo, grid_type), repo)
    repo = release_repo(scenario, home=home)
    release = find_release(repo, scenario.release, scenario.mesa)
    files = release_files(release, grid_type.name, scenario.rerun)
    return _resolved(files, repo, release, grid_type.name, scenario.rerun)
