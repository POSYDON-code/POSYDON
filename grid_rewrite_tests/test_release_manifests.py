"""The restructured inlist repo against today's scenario rules.

POSYDON-MESA-INLISTS (branch ``v2_sources``) keeps the inlists of every v2
production commit in ``v2_sources/inlists/<hash>_<branch>/`` and the new
layout in ``r11701/releases/v2.0/``, where each of those commits is
``v2.0`` plus a rerun. These tests check that reading the new layout
(``manifest.release_files``) gives what today's rules
(``scenario.scenario_files``) give on the commit itself:

- for every archived commit and grid type, every target merged key by key,
  the column lists and the Fortran sources byte for byte, the ZAMS models
  and the stars built from formation steps;
- for every grid type, ``build_inlists`` writes byte-identical files from
  either.

Needs a clone of the inlist repo: ``POSYDON_INLIST_REPO=<path>``, at the
commit in ``INLIST_REPO_COMMIT`` (CI clones it).
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import os
import subprocess
from pathlib import Path

import pytest
from test_setup_grid_baseline import CASES as GOLDEN_CASES
from test_setup_grid_baseline import inlist_origin  # noqa: F401
from test_setup_grid_baseline import GRID_COLUMNS, build_ini, make_setup_env

from posydon.grids.setup.build import build_inlists
from posydon.grids.setup.grid_types import GRID_TYPES
from posydon.grids.setup.inlist import read_inlist
from posydon.grids.setup.legacy import read_legacy_config
from posydon.grids.setup.manifest import (
    find_release,
    list_reruns,
    release_files,
)
from posydon.grids.setup.scenario import scenario_files

INLIST_REPO = os.environ.get("POSYDON_INLIST_REPO", "")
# v2_sources after "Start releases/development from ..."
INLIST_REPO_COMMIT = "8482d7d2d6f4abd3374437053aa6cf2b11ff70ca"
BASE = "a060b7b_main"
# rerun names that are not the archive folder's branch name
RERUN_NAMES = {"18710e9_main": "main_18710e9",
               "77640b5_seth_lbv_thermo_dedt_hepulse":
                   "seth_lbv_thermo_dedt_hepulse_77640b5",
               "c4f90ce_development": "development_c4f90ce"}

pytestmark = pytest.mark.skipif(
    not INLIST_REPO, reason="set POSYDON_INLIST_REPO to a clone of "
    "POSYDON-MESA-INLISTS at {}".format(INLIST_REPO_COMMIT[:7]))


def archived_commits():
    if not INLIST_REPO:
        return []
    return sorted(path.name for path in
                  (Path(INLIST_REPO) / "v2_sources" / "inlists").iterdir()
                  if path.is_dir())


def rerun_of(commit):
    if commit == BASE:
        return None
    return RERUN_NAMES.get(commit, commit.split("_", 1)[1])


@pytest.fixture(scope="module")
def repo():
    head = subprocess.run(["git", "-C", INLIST_REPO, "rev-parse", "HEAD"],
                          capture_output=True, text=True).stdout.strip()
    assert head == INLIST_REPO_COMMIT, \
        "POSYDON_INLIST_REPO is at {}, expected {}".format(head,
                                                          INLIST_REPO_COMMIT)
    return Path(INLIST_REPO)


@pytest.fixture(scope="module")
def old_roots(repo, tmp_path_factory):
    """Each archived commit as a repo root that scenario_files can read."""
    roots = {}
    base = tmp_path_factory.mktemp("old_layout")
    for commit in archived_commits():
        r11701 = base / commit / "r11701"
        r11701.mkdir(parents=True)
        for folder in (repo / "v2_sources" / "inlists" / commit).iterdir():
            (r11701 / folder.name).symlink_to(folder)
        (r11701 / "ZAMS_models").symlink_to(
            repo / "v2_sources" / "zams_models")
        roots[commit] = base / commit
    return roots


SECTIONS = {"binary": ("binary_controls", "binary_job")}


def merged(files):
    """Every target of a ScenarioFiles, merged key by key."""
    out = {}
    for target in set(files.inlists.targets) | set(files.flags):
        base, part = target.rsplit("_", 1)
        section = SECTIONS.get(base.split("_formation")[0],
                               ("controls", "star_job"))[part == "job"]
        values = {}
        for path in files.inlists.targets.get(target, []):
            values.update(read_inlist(path).get(section, {}))
        values.update(files.flags.get(target, {}))
        out[target] = values
    return out


def contents(paths):
    return {key: Path(path).read_bytes() for key, path in paths.items()}


def pairs():
    for commit in archived_commits():
        for grid in GRID_TYPES:
            if grid == "HeMS-HMS" and not (
                    Path(INLIST_REPO) / "v2_sources" / "inlists" / commit /
                    "HeMS-HMS").is_dir():
                continue  # v2.0 has HeMS-HMS; this commit did not
            yield pytest.param(commit, grid, id="{}-{}".format(commit, grid))


@pytest.mark.parametrize("commit, grid", list(pairs()))
def test_release_reproduces_commit(commit, grid, repo, old_roots):
    old = scenario_files(old_roots[commit], GRID_TYPES[grid])
    new = release_files(find_release(repo, "v2.0"), grid, rerun_of(commit))
    old_merged, new_merged = merged(old), merged(new)
    # a step holding only the flag has an empty &star_job either way
    for target in set(old_merged) ^ set(new_merged):
        assert old_merged.get(target, {}) == new_merged.get(target, {})
    for target in set(old_merged) & set(new_merged):
        assert new_merged[target] == old_merged[target], target
    assert contents(new.inlists.whole_files) == \
        contents(old.inlists.whole_files)
    assert contents(new.extras) == contents(old.extras)
    assert new.zams_star1 == old.zams_star1 is None
    assert new.zams_star2 == old.zams_star2 is None
    assert new.formation_stars == old.formation_stars


def test_every_commit_is_checked(repo):
    """18 archived commits: v2.0 itself and 17 reruns."""
    commits = archived_commits()
    assert len(commits) == 18
    reruns = list_reruns(find_release(repo, "v2.0"))
    assert sorted(rerun_of(c) for c in commits if c != BASE) == reruns


SCENARIO_CASES = [case for case in GOLDEN_CASES
                  if case.name.endswith("__scenario")]


@pytest.mark.parametrize("case", SCENARIO_CASES,
                         ids=[case.name for case in SCENARIO_CASES])
def test_build_inlists_is_byte_identical(case, repo, old_roots, tmp_path,
                                         monkeypatch,
                                         inlist_origin):  # noqa: F811
    """The inlists written from v2.0 equal those written from a060b7b."""
    env = make_setup_env(tmp_path, inlist_origin)
    monkeypatch.setenv("MESA_DIR", str(env["mesa"]))
    ini = tmp_path / "grid.ini"
    ini.write_text(build_ini(case, env))
    config = read_legacy_config(ini, mesa_dir=env["mesa"])
    grid = GRID_TYPES[config.grid.type]
    commit = "35dd7a4_jeff_new_grid" if grid.name == "HeMS-HMS" else BASE
    params = GRID_COLUMNS[case.columns].split(",")
    written = {}
    for side, scenario in (
            ("old", scenario_files(old_roots[commit], grid)),
            ("new", release_files(find_release(repo, "v2.0"), grid.name,
                                  rerun_of(commit)))):
        run = tmp_path / side
        build_inlists(config, run, params, scenario=scenario)
        written[side] = {
            path.relative_to(run).as_posix():
                path.read_text().replace(str(run), "<RUN>")
            for path in sorted(run.rglob("inlist*"))}
    assert len(written["old"]) >= 3  # binary/inlist* at least
    assert written["new"] == written["old"]
