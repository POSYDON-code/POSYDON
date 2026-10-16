"""Unit tests of posydon/grids/setup/scenario.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

# import other needed code for the tests, which is not already imported in the
# module you like to test
import shutil
import subprocess
from dataclasses import is_dataclass
from inspect import isroutine
from pathlib import Path

import pytest

# import the module which will be tested
import posydon.grids.setup.scenario as totest
from posydon.grids.setup.config import ConfigError, ScenarioSection
from posydon.grids.setup.grid_types import GRID_TYPES

GIT_ENV = {"GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t.invalid",
           "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@t.invalid"}


def make_repo(root, folders=("CO-H_star", "CO-He_star", "HeMS-HMS"),
              zams=False):
    """A fake POSYDON-MESA-INLISTS tree with the files setup looks for."""
    base = root / "r11701"
    common = base / "default_common_inlists"
    files = ["binary/inlist_project", "binary/inlist1", "binary/inlist2",
             "binary/src/run_binary_extras.f", "binary/src/run_star_extras.f",
             "history_columns.list", "binary_history_columns.list",
             "profile_columns.list"]
    for name in files:
        (common / name).parent.mkdir(parents=True, exist_ok=True)
        (common / name).write_text("&controls\n/\n")
    if zams:
        (common / "binary" / "inlist1").write_text(
            "&controls\n zams_filename = 'models/zams_z1.data'\n/\n")
        (base / "ZAMS_models").mkdir()
        (base / "ZAMS_models" / "zams_z1.data").write_text("")
    layouts = {
        "CO-H_star": ["binary/inlist_project", "binary/inlist1",
                      "binary/inlist2"],
        "CO-He_star": ["binary/inlist_project", "binary/inlist1",
                       "star1_formation/inlist_step1",
                       "star1_formation/inlist_step2",
                       "star1_formation/src/run_star_extras.f"],
        "HeMS-HMS": ["binary/inlist1", "star1_formation/inlist_step1",
                     "star1_formation/inlist_step2"],
    }
    for folder in folders:
        for name in layouts[folder]:
            path = base / folder / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("&controls\n/\n")
    return root


def git(repo, *args):
    env = dict(GIT_ENV, PATH="/usr/bin:/bin")
    return subprocess.run(["git", "-C", str(repo), *args], check=True,
                          capture_output=True, text=True,
                          env=env).stdout.strip()


@pytest.fixture
def origin(tmp_path):
    """A git repository with two commits on main."""
    if shutil.which("git") is None:
        pytest.skip("git is needed")
    repo = make_repo(tmp_path / "origin")
    git(repo, "init", "-q", "-b", "main")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "first")
    first = git(repo, "rev-parse", "HEAD")
    (repo / "r11701" / "new_file").write_text("")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "second")
    return repo, first, git(repo, "rev-parse", "HEAD")


class TestElements:
    def test_dir(self):
        elements = {'COLUMN_LISTS', 'COMMON', 'ConfigError', 'InlistFiles',
                    'MESA_VERSION', 'Path', 'REPOSITORIES',
                    'SINGLE_STAR_FLAG', 'ScenarioFiles', '__authors__',
                    '__builtins__', '__cached__', '__doc__', '__file__',
                    '__loader__', '__name__', '__package__', '__spec__',
                    '_add_every_step', '_add_formation', '_add_star', '_git',
                    '_formation_steps', '_zams_from_inlist', 'dataclass',
                    'fetch_inlist_repo', 'field', 'inlist_repo_path',
                    'read_inlist', 'scenario_files', 'subprocess'}
        assert set(dir(totest)) == elements

    def test_instance_constants(self):
        assert totest.MESA_VERSION == "r11701"
        assert totest.SINGLE_STAR_FLAG == {"x_logical_ctrl(1)": ".true."}
        assert totest.REPOSITORIES["posydon"][1] == ".posydon_mesa_inlists"

    def test_instance_ScenarioFiles(self):
        assert is_dataclass(totest.ScenarioFiles)

    @pytest.mark.parametrize("name", ['inlist_repo_path',
                                      'fetch_inlist_repo', 'scenario_files'])
    def test_instance_functions(self, name):
        assert isroutine(getattr(totest, name))


class TestInlistRepoPath:
    def test_default_locations(self, tmp_path):
        assert totest.inlist_repo_path(ScenarioSection(), home=tmp_path) == \
            tmp_path / ".posydon_mesa_inlists"
        assert totest.inlist_repo_path(ScenarioSection(source="user"),
                                       home=tmp_path) == \
            tmp_path / ".user_mesa_inlists"

    def test_home_default(self):
        assert totest.inlist_repo_path(ScenarioSection()) == \
            Path.home() / ".posydon_mesa_inlists"

    def test_inlist_repo_setting(self, tmp_path):
        scenario = ScenarioSection(inlist_repo=tmp_path / "mine")
        assert totest.inlist_repo_path(scenario) == tmp_path / "mine"

    def test_unknown_source(self):
        with pytest.raises(ConfigError, match="no repository is known"):
            totest.inlist_repo_path(ScenarioSection(source="path"))


class TestFetchInlistRepo:
    def test_existing_clone_is_checked_out(self, origin, tmp_path):
        repo, first, _ = origin
        clone = tmp_path / "clone"
        subprocess.run(["git", "clone", "-q", str(repo), str(clone)],
                       check=True)
        scenario = ScenarioSection(branch="main", commit=first,
                                   inlist_repo=clone)
        assert totest.fetch_inlist_repo(scenario) == clone
        assert git(clone, "rev-parse", "HEAD") == first
        assert not (clone / "r11701" / "new_file").exists()

    def test_pull_brings_new_commits(self, origin, tmp_path):
        repo, first, _ = origin
        clone = tmp_path / "clone"
        subprocess.run(["git", "clone", "-q", str(repo), str(clone)],
                       check=True)
        (repo / "r11701" / "later").write_text("")
        git(repo, "add", "-A")
        git(repo, "commit", "-q", "-m", "third")
        third = git(repo, "rev-parse", "HEAD")
        scenario = ScenarioSection(branch="main", commit=third,
                                   inlist_repo=clone)
        totest.fetch_inlist_repo(scenario)
        assert git(clone, "rev-parse", "HEAD") == third

    def test_no_branch(self, origin, tmp_path):
        repo, first, _ = origin
        clone = tmp_path / "clone"
        subprocess.run(["git", "clone", "-q", str(repo), str(clone)],
                       check=True)
        totest.fetch_inlist_repo(ScenarioSection(commit=first,
                                                 inlist_repo=clone))
        assert git(clone, "rev-parse", "HEAD") == first

    def test_missing_clone_is_cloned(self, origin, tmp_path, monkeypatch):
        repo, first, _ = origin
        monkeypatch.setitem(totest.REPOSITORIES, "posydon",
                            (str(repo), ".posydon_mesa_inlists"))
        scenario = ScenarioSection(branch="main", commit=first)
        clone = totest.fetch_inlist_repo(scenario, home=tmp_path)
        assert clone == tmp_path / ".posydon_mesa_inlists"
        assert git(clone, "rev-parse", "HEAD") == first

    def test_clone_fails(self, tmp_path, monkeypatch):
        monkeypatch.setitem(totest.REPOSITORIES, "posydon",
                            (str(tmp_path / "missing"),
                             ".posydon_mesa_inlists"))
        with pytest.raises(ConfigError, match="could not clone"):
            totest.fetch_inlist_repo(ScenarioSection(commit="1234567"),
                                     home=tmp_path)

    def test_unknown_commit(self, origin, tmp_path):
        repo, _, _ = origin
        clone = tmp_path / "clone"
        subprocess.run(["git", "clone", "-q", str(repo), str(clone)],
                       check=True)
        scenario = ScenarioSection(branch="main", commit="1234567",
                                   inlist_repo=clone)
        with pytest.raises(ConfigError,
                           match=r"\[scenario\] commit: could not check out"):
            totest.fetch_inlist_repo(scenario)


def rel(path):
    """A path from r11701/ on."""
    parts = Path(path).parts
    return "/".join(parts[parts.index("r11701") + 1:])


def short(result):
    return {target: [rel(path) for path in paths]
            for target, paths in result.inlists.targets.items()}


C = "default_common_inlists/binary/"


class TestScenarioFiles:
    def test_HMS_HMS(self, tmp_path):
        result = totest.scenario_files(make_repo(tmp_path),
                                       GRID_TYPES["HMS-HMS"])
        assert short(result) == {
            "binary_controls": [C + "inlist_project"],
            "binary_job": [C + "inlist_project"],
            "star1_controls": [C + "inlist1"], "star1_job": [C + "inlist1"],
            "star2_controls": [C + "inlist2"], "star2_job": [C + "inlist2"]}
        assert {name: rel(path) for name, path in
                result.inlists.whole_files.items()} == {
            "history_columns": "default_common_inlists/history_columns.list",
            "binary_history_columns":
                "default_common_inlists/binary_history_columns.list",
            "profile_columns": "default_common_inlists/profile_columns.list"}
        assert {key: rel(path) for key, path in result.extras.items()} == {
            "posydon_binary_extras": C + "src/run_binary_extras.f",
            "posydon_star_binary_extras": C + "src/run_star_extras.f",
            "mesa_star1_extras": C + "src/run_star_extras.f"}
        assert result.flags == {}
        assert result.formation_stars == ()
        assert result.zams_star1 is None and result.zams_star2 is None

    def test_CO_HMS_RLO(self, tmp_path):
        result = totest.scenario_files(make_repo(tmp_path),
                                       GRID_TYPES["CO-HMS_RLO"])
        assert short(result)["binary_controls"] == [
            C + "inlist_project", "CO-H_star/binary/inlist_project"]
        assert short(result)["star2_job"] == [C + "inlist2",
                                              "CO-H_star/binary/inlist2"]
        assert "user_binary_extras" not in result.extras

    def test_CO_HeMS(self, tmp_path):
        result = totest.scenario_files(make_repo(tmp_path),
                                       GRID_TYPES["CO-HeMS"])
        found = short(result)
        assert found["star1_controls"] == [C + "inlist1",
                                           "CO-He_star/binary/inlist1"]
        assert found["star2_controls"] == [C + "inlist2"]
        assert found["star1_formation_every_step_controls"] == \
            [C + "inlist1"]
        assert found["star1_formation_step0_controls"] == \
            ["CO-He_star/star1_formation/inlist_step1"]
        assert found["star1_formation_step1_job"] == \
            ["CO-He_star/star1_formation/inlist_step2"]
        assert result.formation_stars == ("star1",)

    def test_folder_columns_and_sources(self, tmp_path):
        repo = make_repo(tmp_path)
        folder = repo / "r11701" / "CO-H_star"
        (folder / "src").mkdir()
        (folder / "src" / "run_binary_extras.f").write_text("")
        (folder / "src" / "run_star_extras.f").write_text("")
        (folder / "profile_columns.list").write_text("")
        (folder / "star2_formation").mkdir()
        (folder / "star2_formation" / "inlist_step1").write_text("")
        result = totest.scenario_files(repo, GRID_TYPES["CO-HMS_RLO"])
        assert rel(result.inlists.whole_files["profile_columns"]) == \
            "CO-H_star/profile_columns.list"
        assert rel(result.extras["user_binary_extras"]) == \
            "CO-H_star/src/run_binary_extras.f"
        assert rel(result.extras["user_star_binary_extras"]) == \
            "CO-H_star/src/run_star_extras.f"
        assert result.formation_stars == ("star2",)
        assert short(result)["star2_formation_every_step_job"] == \
            [C + "inlist1"]

    def test_empty_folder(self, tmp_path):
        result = totest.scenario_files(make_repo(tmp_path, folders=()),
                                       GRID_TYPES["CO-HMS_RLO"])
        assert short(result) == short(totest.scenario_files(
            make_repo(tmp_path / "b"), GRID_TYPES["HMS-HMS"]))

    def test_HeMS_HMS(self, tmp_path):
        result = totest.scenario_files(make_repo(tmp_path),
                                       GRID_TYPES["HeMS-HMS"])
        found = short(result)
        assert found["star1_controls"] == [C + "inlist1",
                                           "HeMS-HMS/binary/inlist1"]
        assert found["binary_controls"] == [C + "inlist_project"]
        assert found["star1_formation_step1_controls"] == \
            ["HeMS-HMS/star1_formation/inlist_step2"]
        assert result.formation_stars == ("star1",)

    def test_HeMS_HMS_without_files(self, tmp_path):
        repo = make_repo(tmp_path, folders=())
        (repo / "r11701" / "HeMS-HMS").mkdir()
        result = totest.scenario_files(repo, GRID_TYPES["HeMS-HMS"])
        assert short(result)["star1_controls"] == [C + "inlist1"]
        assert result.formation_stars == ()

    def test_single_HMS(self, tmp_path):
        result = totest.scenario_files(make_repo(tmp_path),
                                       GRID_TYPES["single_HMS"])
        assert short(result) == {
            "star1_formation_every_step_controls": [C + "inlist1"],
            "star1_formation_every_step_job": [C + "inlist1"]}
        assert result.flags == {"star1_formation_every_step_controls":
                                {"x_logical_ctrl(1)": ".true."}}

    def test_single_HeMS(self, tmp_path):
        """As the v2 grids: two He-star steps, then the flag as step 2."""
        result = totest.scenario_files(make_repo(tmp_path),
                                       GRID_TYPES["single_HeMS"])
        assert short(result) == {
            "star1_formation_every_step_controls": [C + "inlist1"],
            "star1_formation_every_step_job": [C + "inlist1"],
            "star1_formation_step0_controls":
                ["CO-He_star/star1_formation/inlist_step1"],
            "star1_formation_step0_job":
                ["CO-He_star/star1_formation/inlist_step1"],
            "star1_formation_step1_controls":
                ["CO-He_star/star1_formation/inlist_step2"],
            "star1_formation_step1_job":
                ["CO-He_star/star1_formation/inlist_step2"]}
        assert result.flags == {
            "star1_formation_step2_controls": {"x_logical_ctrl(1)": ".true."},
            "star1_formation_step2_job": {}}

    def test_zams_named_in_common_inlist1(self, tmp_path):
        result = totest.scenario_files(make_repo(tmp_path, zams=True),
                                       GRID_TYPES["HMS-HMS"])
        assert rel(result.zams_star1) == "ZAMS_models/zams_z1.data"
        assert result.zams_star2 is None

    def test_zams_file_missing(self, tmp_path):
        repo = make_repo(tmp_path, zams=True)
        (repo / "r11701" / "ZAMS_models" / "zams_z1.data").unlink()
        result = totest.scenario_files(repo, GRID_TYPES["HMS-HMS"])
        assert result.zams_star1 is None

    def test_zams_without_inlist(self, tmp_path):
        repo = make_repo(tmp_path)
        (repo / "r11701" / "default_common_inlists" / "binary" /
         "inlist2").unlink()
        result = totest.scenario_files(repo, GRID_TYPES["HMS-HMS"])
        assert result.zams_star2 is None
