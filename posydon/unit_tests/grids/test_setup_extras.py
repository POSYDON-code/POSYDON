"""Unit tests of posydon/grids/setup/extras.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

# import other needed code for the tests, which is not already imported in the
# module you like to test
import doctest
from inspect import isroutine
from pathlib import Path

import pytest

# import the module which will be tested
import posydon.grids.setup.extras as totest
from posydon.grids.setup.config import GridConfig
from posydon.grids.setup.scenario import ScenarioFiles


@pytest.fixture
def sources(tmp_path):
    folder = tmp_path / "sources"
    folder.mkdir()
    for name in ("run_binary_extras.f", "run_star_extras.f",
                 "binary_run.f", "run.f", "makefile", "notes.txt"):
        (folder / name).write_text("! " + name + "\n")
    return folder


class TestElements:
    def test_dir(self):
        elements = {'COLUMN_LIST_FILES', 'Path', 'SOURCE_FOLDERS',
                    '__authors__', '__builtins__', '__cached__', '__doc__',
                    '__file__', '__loader__', '__name__', '__package__',
                    '__spec__', 'column_list_sources', 'compile_executables',
                    'destinations', 'keep_highest_binary_extras',
                    'merge_extras', 'mk_lines', 'os', 'pick_whole_file',
                    'shutil', 'stage_column_lists', 'stage_extras',
                    'subprocess'}
        assert set(dir(totest)) == elements

    @pytest.mark.parametrize("name", [
        'merge_extras', 'keep_highest_binary_extras', 'destinations',
        'mk_lines', 'stage_extras', 'compile_executables',
        'column_list_sources', 'stage_column_lists'])
    def test_instance_functions(self, name):
        assert isroutine(getattr(totest, name))

    def test_doctests(self):
        result = doctest.testmod(totest)
        assert result.failed == 0
        assert result.attempted > 0


class TestMergeExtras:
    def test_scenario_replaces_in_place_and_appends(self):
        extras = totest.merge_extras({"a": "1", "b": "2"},
                                     {"b": "3", "c": "4"})
        assert extras == {"a": Path("1"), "b": Path("3"), "c": Path("4")}
        assert list(extras) == ["a", "b", "c"]

    def test_no_scenario(self):
        assert totest.merge_extras({"a": "1"}) == {"a": Path("1")}


class TestKeepHighestBinaryExtras:
    def test_user_wins(self):
        extras = {"mesa_binary_extras": 1, "mesa_star_binary_extras": 2,
                  "posydon_star_binary_extras": 3, "user_binary_extras": 4,
                  "star_run": 5}
        assert totest.keep_highest_binary_extras(extras) == {
            "user_binary_extras": 4, "star_run": 5}

    def test_no_binary_extras(self):
        extras = {"star_run": 1}
        assert totest.keep_highest_binary_extras(extras) == extras


class TestDestinations:
    @pytest.mark.parametrize("key, expected", [
        ("posydon_binary_extras", ["binary/src/"]),
        ("binary_run", ["binary/src/"]),
        ("star_run", ["star1/src/", "star2/src/"]),
        ("mesa_star1_extras", ["star1/src/"]),
        ("user_star2_extras", ["star2/src/"]),
        ("makefile_binary", ["binary/make/makefile_binary"]),
        ("makefile_star", ["star1/make/makefile_star",
                           "star2/make/makefile_star"]),
        ("mesa_dir", []),
        ("something_else", [""]),
    ])
    def test_destinations(self, key, expected):
        assert totest.destinations(key) == expected


class TestStageAndCompile:
    def test_stage(self, sources, tmp_path):
        run = tmp_path / "run"
        stale = run / "binary" / "src" / "old.f"
        stale.parent.mkdir(parents=True)
        stale.write_text("")
        extras = {"makefile_binary": sources / "makefile",
                  "makefile_star": sources / "makefile",
                  "posydon_binary_extras": sources / "run_binary_extras.f",
                  "binary_run": sources / "binary_run.f",
                  "star_run": sources / "run.f",
                  "mesa_star1_extras": sources / "run_star_extras.f",
                  "mesa_star2_extras": None,
                  "other": sources / "notes.txt",
                  "mesa_dir": sources}
        mk = totest.stage_extras(extras, run)
        assert not stale.exists()
        found = sorted(path.relative_to(run).as_posix()
                       for path in run.rglob("*") if path.is_file())
        assert found == ["binary/make/makefile_binary",
                         "binary/src/binary_run.f",
                         "binary/src/run_binary_extras.f", "mk", "notes.txt",
                         "star1/make/makefile_star", "star1/src/run.f",
                         "star1/src/run_star_extras.f",
                         "star2/make/makefile_star", "star2/src/run.f"]
        assert mk.read_text() == (
            "#!/bin/sh\n"
            "cd {0}/binary/make\nmake -f makefile_binary\n"
            "cd {0}/star1/make\nmake -f makefile_star\n"
            "cd {0}/star2/make\nmake -f makefile_star\n".format(run))
        assert mk.stat().st_mode & 0o111

    def test_mk_runs_on_its_own(self, sources, tmp_path):
        """C15: ./mk works without sh in front."""
        run = tmp_path / "run"
        mk = totest.stage_extras({"other": sources / "notes.txt"}, run)
        mk.write_text(mk.read_text() + "touch ran\n")
        result = totest.subprocess.run(["./mk"], cwd=run)
        assert result.returncode == 0
        assert (run / "ran").exists()

    def test_mk_lines_without_makefiles(self):
        assert totest.mk_lines({"star_run": "x"}, "/run") == []

    def test_compile_runs_mk_in_the_run_directory(self, tmp_path,
                                                  monkeypatch):
        run = tmp_path / "run"
        run.mkdir()
        (run / "mk").write_text("touch compiled\n")
        monkeypatch.chdir(tmp_path)
        assert totest.compile_executables(run) == (
            str(run / "binary" / "binary"), str(run / "star1" / "star"),
            str(run / "star2" / "star"))
        assert (run / "compiled").exists()
        assert not (tmp_path / "compiled").exists()


class TestColumnLists:
    def config(self, sources):
        config = GridConfig()
        whole = config.inlists["posydon_inlists"].whole_files
        for name in totest.COLUMN_LIST_FILES:
            whole[name] = sources / "notes.txt"
        return config

    def test_ini_files(self, sources):
        chosen = totest.column_list_sources(self.config(sources))
        assert chosen == {name: sources / "notes.txt"
                          for name in totest.COLUMN_LIST_FILES}

    def test_scenario_and_user_win(self, sources):
        config = self.config(sources)
        config.inlists["user_inlists"].whole_files["profile_columns"] = \
            sources / "run.f"
        scenario = ScenarioFiles()
        scenario.inlists.whole_files["history_columns"] = sources / "makefile"
        chosen = totest.column_list_sources(config, scenario)
        assert chosen["history_columns"] == sources / "makefile"
        assert chosen["binary_history_columns"] == sources / "notes.txt"
        assert chosen["profile_columns"] == sources / "run.f"

    def test_missing(self):
        assert totest.column_list_sources(GridConfig()) == {}

    def test_stage_new_folder(self, sources, tmp_path):
        staged = totest.stage_column_lists(
            totest.column_list_sources(self.config(sources)), tmp_path)
        assert all(Path(path).is_file() for path in staged.values())

    def test_stage(self, sources, tmp_path):
        run = tmp_path / "run"
        (run / "column_lists").mkdir(parents=True)
        (run / "column_lists" / "stale.list").write_text("")
        staged = totest.stage_column_lists(
            totest.column_list_sources(self.config(sources)), run)
        assert staged == {name: str(run / "column_lists" / file_name)
                          for name, file_name in
                          totest.COLUMN_LIST_FILES.items()}
        assert sorted(path.name for path in
                      (run / "column_lists").iterdir()) == \
            sorted(totest.COLUMN_LIST_FILES.values())
