"""The new posydon.grids.setup modules against today's code, piece by piece.

Each test runs a function of today's bin/posydon-setup-grid (loaded as a
module) and its replacement on the same input, and compares the results.
They complement the end-to-end goldens in test_setup_grid_baseline.py.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import importlib.machinery
import importlib.util
import shutil
import subprocess
from pathlib import Path

import pytest
from test_setup_grid_baseline import SETUP_GRID, inlist_origin  # noqa: F401

from posydon.grids.setup import legacy
from posydon.grids.setup.config import GridConfig, ScenarioSection
from posydon.grids.setup.grid_types import GRID_TYPES
from posydon.grids.setup.inlist import read_inlist
from posydon.grids.setup.layers import expand_every_step, formation_step_count
from posydon.grids.setup.scenario import fetch_inlist_repo, scenario_files


def load_setup_grid():
    """bin/posydon-setup-grid as a module (its __main__ block is skipped)."""
    loader = importlib.machinery.SourceFileLoader("legacy_setup_grid",
                                                  str(SETUP_GRID))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def setup_grid():
    return load_setup_grid()


@pytest.fixture
def home(tmp_path, inlist_origin, monkeypatch):  # noqa: F811
    """A HOME with a clone of the fixture inlist repo, as today expects."""
    origin, hashes = inlist_origin
    home = tmp_path / "home"
    home.mkdir()
    subprocess.run(["git", "clone", "-q", str(origin),
                    str(home / ".posydon_mesa_inlists")], check=True,
                   capture_output=True)
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.chdir(tmp_path)
    return home, hashes


# grid type -> (scenario branch, legacy system type)
CASES = {
    "HMS-HMS": ("main", "HMS-HMS"),
    "CO-HMS_RLO": ("main", "CO-H_star"),
    "CO-HeMS": ("main", "CO-He_star"),
    "HeMS-HMS": ("jeff_new_grid", "HeMS-HMS"),
    "single_HMS": ("main", "HMS-HMS"),
}


@pytest.mark.parametrize("name", list(CASES))
def test_scenario_files_match_find_inlist_from_scenario(name, home,
                                                        setup_grid):
    home_dir, hashes = home
    branch, system_type = CASES[name]
    grid_type = GRID_TYPES[name]
    commit = hashes[branch]

    # today: the function fills the module-level mesa_inlists/mesa_extras
    setup_grid.mesa_inlists = {"single_star_grid": grid_type.single_star}
    setup_grid.mesa_extras = {}
    setup_grid.find_inlist_from_scenario(
        source="posydon", gitcommit="{}-{}".format(branch, commit),
        system_type=system_type)
    old = dict(setup_grid.mesa_inlists)
    old_extras = dict(setup_grid.mesa_extras)
    special = old.pop("star1_controls_special", None)

    # new
    scenario = ScenarioSection(branch=branch, commit=commit)
    repo = fetch_inlist_repo(scenario, home=home_dir)
    new = scenario_files(repo, grid_type)

    # inlist files per target: today's posydon_defaults then _user keys,
    # read the way legacy.py reads them, against the new posydon layer
    expected = GridConfig()
    legacy._read_inlist_keys(old, expected, grid_type.single_star)
    targets = {}
    for layer in ("posydon_inlists", "user_inlists"):
        for target, paths in expected.inlists[layer].targets.items():
            targets.setdefault(target, []).extend(paths)
    if grid_type.single_star:
        # today also lists binary and star 2 files for single-star grids,
        # but construct_static_inlist never reads them in that mode
        targets = {target: paths for target, paths in targets.items()
                   if not target.startswith(("binary_", "star2_"))}
    # "one file per step" (today) and "every step" (new) mean the same
    old_targets, new_targets = targets, dict(new.inlists.targets)
    for star in ("star1", "star2"):
        n_steps = formation_step_count(star, [old_targets, new_targets])
        old_targets = expand_every_step(old_targets, star, n_steps)
        new_targets = expand_every_step(new_targets, star, n_steps)
    assert new_targets == old_targets

    # the single-star flag: a file today, values now
    if special is None:
        assert new.flags == {}
    else:
        assert new.flags == {"star1_formation_every_step_controls":
                             read_inlist(special)["controls"]}

    # column lists, extras and ZAMS models
    assert new.inlists.whole_files == {
        name: Path(old[key]) for key, name in legacy.COLUMN_KEYS.items()}
    assert {key: Path(value) for key, value in new.extras.items()} == \
        {key: Path(value) for key, value in old_extras.items()}
    assert new.zams_star1 == old.get("zams_filename_1")
    assert new.zams_star2 == old.get("zams_filename_2")


def test_single_HeMS_crashes_today(home, setup_grid):
    """The reason single_HeMS is compared with the v2 grids instead."""
    _, hashes = home
    setup_grid.mesa_inlists = {"single_star_grid": True}
    setup_grid.mesa_extras = {}
    with pytest.raises(KeyError, match="star1_controls_special"):
        setup_grid.find_inlist_from_scenario(
            source="posydon", gitcommit="main-" + hashes["main"],
            system_type="CO-He_star")


@pytest.mark.skipif(shutil.which("git") is None, reason="needs git")
def test_new_code_leaves_the_clone_clean(home):
    home_dir, hashes = home
    scenario = ScenarioSection(branch="main", commit=hashes["main"])
    repo = fetch_inlist_repo(scenario, home=home_dir)
    scenario_files(repo, GRID_TYPES["single_HMS"])
    status = subprocess.run(["git", "-C", str(repo), "status",
                             "--porcelain"], capture_output=True, text=True)
    assert status.stdout == ""
