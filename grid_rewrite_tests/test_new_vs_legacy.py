"""The new posydon.grids.setup modules against today's code, piece by piece.

Each test runs a function of the pre-rewrite bin/posydon-setup-grid (kept as
legacy_posydon_setup_grid.py and loaded as a module) and its replacement on the same input, and compares the results.
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
from test_setup_grid_baseline import LEGACY_SETUP_GRID  # noqa: F401
from test_setup_grid_baseline import inlist_origin  # noqa: F401

from posydon.grids.setup import legacy
from posydon.grids.setup.config import GridConfig, ScenarioSection
from posydon.grids.setup.grid_types import GRID_TYPES
from posydon.grids.setup.inlist import read_inlist
from posydon.grids.setup.layers import expand_every_step, formation_step_count
from posydon.grids.setup.scenario import fetch_inlist_repo, scenario_files


def load_setup_grid():
    """The pre-rewrite posydon-setup-grid as a module (its __main__ block is skipped)."""
    loader = importlib.machinery.SourceFileLoader("legacy_setup_grid",
                                                  str(LEGACY_SETUP_GRID))
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


###############################################################################
# build.py against construct_static_inlist, for every golden case
###############################################################################
from test_setup_grid_baseline import CASE_IDS
from test_setup_grid_baseline import CASES as GOLDEN_CASES  # noqa: E402
from test_setup_grid_baseline import GRID_COLUMNS, build_ini, make_setup_env

from posydon.grids.setup.build import build_inlists  # noqa: E402
from posydon.grids.setup.grid_types import get_grid_type  # noqa: E402
from posydon.grids.setup.legacy import read_legacy_config  # noqa: E402
from posydon.utils import configfile  # noqa: E402

# keys set by the fixture user inlists (grid_rewrite_tests/setup_grid_data/
# extra/user_inlist*); with a scenario they never took effect today (S6)
USER_KEYS = {"fm", "mesh_delta_coeff", "pgstar_flag"}


def legacy_build(setup_grid, ini, run, params, monkeypatch, home):
    """main() of today's posydon-setup-grid, up to construct_static_inlist."""
    for folder in ("binary", "star1", "star2"):
        (run / folder).mkdir(parents=True, exist_ok=True)
    _, _, mesa_inlists, mesa_extras = configfile.parse_inifile(str(ini))
    mesa_inlists.setdefault("scenario", None)
    setup_grid.mesa_inlists = mesa_inlists
    setup_grid.mesa_extras = mesa_extras
    if mesa_inlists["scenario"] is not None:
        monkeypatch.setenv("HOME", str(home))
        monkeypatch.chdir(run)
        setup_grid.find_inlist_from_scenario(*mesa_inlists["scenario"])
    return setup_grid.construct_static_inlist(mesa_inlists,
                                              grid_parameters=params,
                                              working_directory=str(run))


def written_inlists(run):
    return {path.relative_to(run).as_posix(): path.read_text()
            for path in sorted(run.glob("binary/inlist*")) +
            sorted(run.glob("star[12]/inlist_step*"))}


@pytest.mark.parametrize("case", GOLDEN_CASES, ids=CASE_IDS)
def test_build_matches_construct_static_inlist(case, tmp_path, monkeypatch,
                                               inlist_origin,  # noqa: F811
                                               setup_grid):
    if case.name == "single_HeMS__scenario":
        pytest.skip("crashes today (S1); compared with the v2 grids instead")
    env = make_setup_env(tmp_path, inlist_origin)
    monkeypatch.setenv("MESA_DIR", str(env["mesa"]))
    ini = env["run"] / "grid.ini"
    ini.write_text(build_ini(case, env))
    params = GRID_COLUMNS[case.columns].split(",")
    if case.grid_type == "dynamic":
        params = ["m1", "m2", "initial_period_in_days"]

    old_run, new_run = tmp_path / "old", tmp_path / "new"
    old = legacy_build(setup_grid, ini, old_run, params, monkeypatch,
                       env["home"])
    config = read_legacy_config(ini, grid_type=case.grid_type,
                                mesa_dir=env["mesa"])
    scenario = None
    if config.scenario.is_set:
        repo = fetch_inlist_repo(config.scenario, home=env["home"])
        scenario = scenario_files(repo, get_grid_type(config.grid.type))
    new = build_inlists(config, new_run, params, scenario=scenario)

    if case.name == "HeMS-HMS__scenario":
        # S2: today zams_filename keeps star 1 from forming; now it is formed
        # in steps and binary/inlist1 loads the formed model instead of the
        # ZAMS file. The binary and star 2 do not change.
        assert old[0] is None and new.star1_formation is not None
        old_files, new_files = written_inlists(old_run), {
            name: text.replace(str(new_run), str(old_run))
            for name, text in written_inlists(new_run).items()}
        assert sorted(set(new_files) - set(old_files)) == [
            "star1/inlist_step0", "star1/inlist_step1"]
        for name in ("binary/inlist_project", "binary/inlist2"):
            assert new_files[name] == old_files[name], name
        old1, new1 = (read_inlist_text(old_files["binary/inlist1"]),
                      read_inlist_text(new_files["binary/inlist1"]))
        # the ini's ZAMS model gives way to MESA's default, as in the v2
        # CO-HeMS runs; unused, since star 1 loads its formed model
        assert old1["controls"].pop("zams_filename") != \
            new1["controls"].pop("zams_filename") == "'zams_z2m2_y28.data'"
        assert new1["controls"] == old1["controls"]
        assert {key: value for key, value in new1["star_job"].items()
                if old1["star_job"].get(key) != value} == {
            "load_saved_model": ".true.",
            "saved_model_name": "'initial_star1_step1.mod'"}
        return

    # what posydon-run-grid is told
    def same(path):
        return None if path is None else path.replace(str(new_run),
                                                       str(old_run))
    assert (same(new.star1_formation), same(new.star2_formation),
            same(new.binary_project), same(new.binary_star1),
            same(new.binary_star2)) == old

    # the files, byte for byte
    old_files = written_inlists(old_run)
    new_files = {name: text.replace(str(new_run), str(old_run))
                 for name, text in written_inlists(new_run).items()}
    assert set(new_files) == set(old_files)
    if case.extra_inlists and case.use_scenario:
        # S6: today the scenario's files win over the user's inlists, and
        # when the ini names a *_user file, the grid type's own folder file
        # takes that key's earlier place, so the common inlists override it.
        # Now: common < grid-type folder < user. Changed keys must be the
        # user's, or keys whose new value is the grid-type folder's.
        folder_values = {}
        for paths in scenario.inlists.targets.values():
            for path in paths:
                if "default_common_inlists" not in str(path):
                    for values in read_inlist(path).values():
                        folder_values.update(values)
        for name in old_files:
            old_values, new_values = (read_inlist_text(old_files[name]),
                                      read_inlist_text(new_files[name]))
            changed = {key: value for section in new_values
                       for key, value in new_values[section].items()
                       if old_values[section].get(key) != value}
            unexplained = {key for key, value in changed.items()
                           if key not in USER_KEYS and
                           folder_values.get(key) != value}
            assert not unexplained, (name, unexplained)
        assert any(old_files[name] != new_files[name] for name in old_files)
    elif case.name == "single_HeMS__no_scenario":
        # S3: without a ZAMS model the formation steps no longer carry
        # MESA's default zams_filename; nothing else may change
        for name in old_files:
            old_values, new_values = (read_inlist_text(old_files[name]),
                                      read_inlist_text(new_files[name]))
            for values in old_values.values():
                values.pop("zams_filename", None)
            assert new_values == old_values, name
        assert any(old_files[name] != new_files[name] for name in old_files)
    else:
        for name in old_files:
            assert new_files[name] == old_files[name], name


def read_inlist_text(text):
    from posydon.grids.setup.inlist import parse_inlist
    return parse_inlist(text)


###############################################################################
# extras.py against make_executables and the column-list copy of main()
###############################################################################
from posydon.grids.setup.extras import column_list_sources  # noqa: E402
from posydon.grids.setup.extras import (
    compile_executables,
    keep_highest_extras,
    merge_extras,
    stage_column_lists,
    stage_extras,
)


def legacy_stage(setup_grid, ini, run, monkeypatch, home):
    """main() of today's posydon-setup-grid: extras, make, column lists."""
    run.mkdir(parents=True, exist_ok=True)
    _, _, mesa_inlists, mesa_extras = configfile.parse_inifile(str(ini))
    mesa_inlists.setdefault("scenario", None)
    setup_grid.mesa_inlists = mesa_inlists
    setup_grid.mesa_extras = mesa_extras
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.chdir(run)
    if mesa_inlists["scenario"] is not None:
        setup_grid.find_inlist_from_scenario(*mesa_inlists["scenario"])
    monkeypatch.chdir(run)  # mk is written to and run from the cwd
    kinds = sorted(set(key.split('_')[0] for key in mesa_extras
                       if 'binary_extras' in key))
    for key in mesa_extras:
        if 'binary_extras' in key and kinds[-1] not in key:
            mesa_extras[key] = None
    executables = setup_grid.make_executables(mesa_extras=mesa_extras,
                                              working_directory=str(run))
    folder = run / "column_lists"
    folder.mkdir()
    for key, name in (("star_history_columns", "history_columns.list"),
                      ("binary_history_columns",
                       "binary_history_columns.list"),
                      ("profile_columns", "profile_columns.list")):
        shutil.copy(mesa_inlists[key], folder / name)
    return executables


def staged_files(run):
    names = ["mk"] + [path.relative_to(run).as_posix() for path in
                      sorted(run.glob("*/src/*")) +
                      sorted(run.glob("*/make/*")) +
                      sorted(run.glob("column_lists/*"))]
    names += [path.name for path in run.iterdir()
              if path.is_file() and path.name not in ("mk", "grid.ini")]
    return {name: (run / name).read_bytes() for name in names}


@pytest.mark.parametrize("case", GOLDEN_CASES, ids=CASE_IDS)
def test_extras_match_make_executables(case, tmp_path, monkeypatch,
                                       inlist_origin,  # noqa: F811
                                       setup_grid):
    if case.name == "single_HeMS__scenario":
        pytest.skip("crashes today (S1)")
    env = make_setup_env(tmp_path, inlist_origin)
    monkeypatch.setenv("MESA_DIR", str(env["mesa"]))
    ini = tmp_path / "grid.ini"
    ini.write_text(build_ini(case, env))
    old_run, new_run = tmp_path / "old", tmp_path / "new"
    old_exes = legacy_stage(setup_grid, ini, old_run, monkeypatch,
                            env["home"])

    config = read_legacy_config(ini, mesa_dir=env["mesa"])
    scenario = None
    if config.scenario.is_set:
        repo = fetch_inlist_repo(config.scenario, home=env["home"])
        scenario = scenario_files(repo, get_grid_type(config.grid.type))
    extras = keep_highest_extras(merge_extras(
        config.extras, scenario.extras if scenario else None))
    new_run.mkdir()
    monkeypatch.chdir(tmp_path)  # the new code must not depend on the cwd
    stage_extras(extras, new_run)
    new_exes = compile_executables(new_run)
    stage_column_lists(column_list_sources(config, scenario), new_run)

    assert [path.replace(str(new_run), str(old_run))
            for path in new_exes] == list(old_exes)
    old_files = staged_files(old_run)
    new_files = {name: data.replace(str(new_run).encode(),
                                    str(old_run).encode())
                 for name, data in staged_files(new_run).items()}
    # C15: mk now starts with a #! line
    assert new_files["mk"].startswith(b"#!/bin/sh\n")
    new_files["mk"] = new_files["mk"][len(b"#!/bin/sh\n"):]
    if scenario is not None:
        # C1: star 2 gets the scenario's run_star_extras.f, as star 1 does
        name = "star2/src/run_star_extras.f"
        assert new_files[name] == new_files["star1/src/run_star_extras.f"]
        assert new_files[name] != old_files[name]
        new_files[name] = old_files[name]
    assert new_files == old_files


###############################################################################
# submit.py against the command line and scripts of main()
###############################################################################
from test_setup_grid_baseline import GRID_ROWS, run_setup  # noqa: E402

from posydon.grids.setup.submit import (  # noqa: E402
    command_line,
    write_scripts,
)

SCRIPTS = ["grid_command.sh", "job_array_grid_submit.slurm",
           "mpi_grid_submit.slurm", "cleanup.slurm", "run_grid.sh"]


@pytest.mark.parametrize("case", GOLDEN_CASES, ids=CASE_IDS)
def test_scripts_match_main(case, tmp_path, monkeypatch,
                            inlist_origin):  # noqa: F811
    env = make_setup_env(tmp_path, inlist_origin)
    old_run = env["run"]
    # the pre-rewrite script, run from old_run
    result = run_setup(case, env, script=LEGACY_SETUP_GRID)
    if result.returncode != 0:
        pytest.skip("fails today: " + result.stderr.strip().splitlines()[-1])
    monkeypatch.setenv("MESA_DIR", str(env["mesa"]))

    new_run = tmp_path / "new"
    new_run.mkdir()
    config = read_legacy_config(old_run / "grid.ini",
                                grid_type=case.grid_type,
                                submission_type=case.submission,
                                mesa_dir=env["mesa"])
    scenario = None
    if config.scenario.is_set:
        repo = fetch_inlist_repo(config.scenario, home=env["home"])
        scenario = scenario_files(repo, get_grid_type(config.grid.type))
    params = GRID_COLUMNS[case.columns].split(",")
    if case.grid_type == "dynamic":
        params = ["m1", "m2", "initial_period_in_days"]
    built = build_inlists(config, new_run, params, scenario=scenario)
    stage_extras(keep_highest_extras(merge_extras(
        config.extras, scenario.extras if scenario else None)), new_run)
    monkeypatch.chdir(tmp_path)  # the new code must not depend on the cwd
    executables = compile_executables(new_run)
    columns = stage_column_lists(column_list_sources(config, scenario),
                                 new_run)
    if case.name == "HMS-HMS__scenario__dynamic_job_array":
        # S5: today this writes a fixed-grid command; now it is refused
        assert "--grid-type fixed" in (old_run /
                                       "job_array_grid_submit.slurm"
                                       ).read_text()
        with pytest.raises(ValueError, match="one MPI job"):
            command_line(config, built, executables, columns, str(new_run),
                         str(env["bin"] / "posydon-run-grid"))
        return
    command = command_line(config, built, executables, columns,
                           str(new_run), str(env["bin"] / "posydon-run-grid"))
    scripts = tmp_path / "scripts"
    scripts.mkdir()
    write_scripts(config, command, len(GRID_ROWS[case.columns]), scripts,
                  "/opt/mesasdk-stub", str(env["mesa"]))

    old = {name: (old_run / name) for name in SCRIPTS
           if (old_run / name).exists()}
    new = {path.name: path for path in scripts.iterdir()}
    assert set(new) == set(old)
    for name in old:
        new_text = new[name].read_text().replace(str(new_run), str(old_run))
        if case.name == "HeMS-HMS__scenario" and \
                "--mesa-star1-inlist-project" in new_text:
            # S2: star 1 is formed now, so posydon-run-grid gets its steps
            steps = " {0}/star1/inlist_step0 {0}/star1/inlist_step1".format(
                old_run)
            assert "--mesa-star1-inlist-project " + steps + " " in new_text
            new_text = new_text.replace(steps, "None")
        assert new_text == old[name].read_text(), name
        assert new[name].stat().st_mode & 0o777 == \
            old[name].stat().st_mode & 0o777, name
