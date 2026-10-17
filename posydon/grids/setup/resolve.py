"""Turn the ini's ``[scenario]`` into the files of ``posydon_inlists``.

There are two kinds of scenario:

- **legacy**: a commit of the inlist repository and a system type, as in
  today's ini files (``scenario = ['posydon', '<branch>-<commit>',
  '<system type>']``). The repository is checked out at the commit and the
  files are found by today's rules (``scenario.scenario_files``).
- **release**: ``release`` (+ ``rerun``) of the restructured repository,
  with ``[grid] type``. The repository is checked out at ``ref`` (a branch
  or tag) or ``commit``, and the files are read from the release's
  manifests (``manifest.release_files``). There is no default ref: the ini
  says which version of the inlist repository it uses.

With ``source = path`` the repository at ``inlist_repo`` is used as it is,
without git; useful while editing the inlists.

``resolve_scenario`` returns the files with the commit they came from, so
the run directory can record it.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

from dataclasses import dataclass
from pathlib import Path

from posydon.grids.setup.config import ConfigError
from posydon.grids.setup.grid_types import get_grid_type
from posydon.grids.setup.manifest import find_release, release_files
from posydon.grids.setup.scenario import (
    REPOSITORIES,
    ScenarioFiles,
    clone_if_missing,
    fetch_inlist_repo,
    git,
    inlist_repo_path,
    scenario_files,
)


@dataclass
class ResolvedScenario:
    """The files of the scenario and where they came from."""

    files: ScenarioFiles
    repo: Path
    commit: str = None
    release: object = None


def _head(repo):
    """The checked-out commit, or None outside a git repository."""
    result = git(repo, "rev-parse", "HEAD")
    return result.stdout.strip() if result.returncode == 0 else None


def checkout_release_repo(scenario, home=None):
    """Make the inlist repository available at the scenario's version.

    Parameters
    ----------
    scenario : ScenarioSection
        With ``ref`` and/or ``commit`` (``commit`` wins), or
        ``source = path``.
    home : str, Path or None
        Folder for the default clone location; defaults to ``~``.

    Returns
    -------
    tuple
        The repository folder and the commit checked out (``None`` for a
        ``path`` source that is not a git repository).

    Raises
    ------
    ConfigError
        If the repository cannot be cloned, the ref or commit is unknown,
        or a ``path`` source is not at the given commit.
    """
    repo = inlist_repo_path(scenario, home=home)
    if scenario.source == "path":
        if not repo.is_dir():
            raise ConfigError("scenario", "inlist_repo", "{} does not "
                              "exist".format(repo))
        head = _head(repo)
        if scenario.commit is not None and \
                (head is None or not head.startswith(scenario.commit)):
            raise ConfigError("scenario", "commit", "{} is at {}, not {}"
                              .format(repo, head, scenario.commit))
        return repo, head
    clone_if_missing(REPOSITORIES[scenario.source][0], repo)
    # a failed fetch (no network) is fine if the version is already here
    fetched = git(repo, "fetch", "--tags", "origin").returncode == 0
    if scenario.commit is not None:
        key, names = "commit", [scenario.commit]
    else:
        key, names = "ref", ["origin/" + scenario.ref, scenario.ref]
    target = None
    for name in names:
        result = git(repo, "rev-parse", "--verify", "--quiet",
                      name + "^{commit}")
        if result.returncode == 0:
            target = result.stdout.strip()
            break
    if target is None:
        raise ConfigError("scenario", key, "{!r} is not in {}{}".format(
            names[-1], repo, "" if fetched else " (and fetching failed)"))
    result = git(repo, "checkout", "--quiet", "--detach", target)
    if result.returncode != 0:
        raise ConfigError("scenario", key, "could not check out {} in {}: "
                          "{}".format(target, repo, result.stderr.strip()))
    return repo, target


def resolve_scenario(config, home=None):
    """The scenario's files, or ``None`` if the ini sets no scenario.

    Parameters
    ----------
    config : GridConfig
        ``config.grid.type`` names the grid type.
    home : str, Path or None
        Folder for the default clone location; defaults to ``~``.

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
        return ResolvedScenario(files=scenario_files(repo, grid_type),
                                repo=repo, commit=_head(repo))
    repo, commit = checkout_release_repo(scenario, home=home)
    release = find_release(repo, scenario.release, scenario.mesa)
    files = release_files(release, grid_type.name, scenario.rerun)
    return ResolvedScenario(files=files, repo=repo, commit=commit,
                            release=release)
