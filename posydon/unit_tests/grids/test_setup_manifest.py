"""Unit tests of posydon/grids/setup/manifest.py
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
import posydon.grids.setup.manifest as totest
from posydon.grids.setup.config import ConfigError

GRID_MANIFEST = """\
[grid_type]
name = CO-HeMS
star1 = formation   ; built in steps
star2 = point_mass

[inlists]
binary = common/binary/inlist_project
star1 = common/binary/inlist1, CO-HeMS/binary/inlist1
star1_formation_every_step = common/binary/inlist1
star1_formation_step0 = CO-HeMS/star1_formation/inlist_step1

[columns]
history_columns = common/history_columns.list

[extras]
posydon_star_binary_extras = common/binary/src/run_star_extras.f

[zams]
star2 = ZAMS_models/zams.data
"""
RERUN_MANIFEST = """\
[rerun]
name = hot

[overlays]
common/binary/inlist1 = common/binary/inlist1

[replace]
common/binary/src/run_star_extras.f = common/binary/src/run_star_extras.f
"""


def write(path, text=""):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


@pytest.fixture
def repo(tmp_path):
    root = tmp_path / "repo"
    release = root / "r11701" / "releases" / "v2.0"
    write(release / "release.ini", "[release]\nmesa_version = r11701\n"
          "mesa = default\nfrozen = true\n")
    for rel in ("common/binary/inlist_project", "common/binary/inlist1",
                "CO-HeMS/binary/inlist1",
                "CO-HeMS/star1_formation/inlist_step1",
                "common/history_columns.list",
                "common/binary/src/run_star_extras.f"):
        write(release / rel, rel)
    write(release / "CO-HeMS" / "scenario.ini", GRID_MANIFEST)
    rerun = release / "reruns" / "hot"
    write(rerun / "scenario.ini", RERUN_MANIFEST)
    write(rerun / "common" / "binary" / "inlist1", "&controls\n x = 1\n/\n")
    write(rerun / "common" / "binary" / "src" / "run_star_extras.f", "new")
    (release / "reruns" / "empty").mkdir()
    for rel in totest.MESA_DEFAULTS.values():
        write(root / "r11701" / "mesa" / "default" / "defaults" / rel)
    return root


class TestElements:
    def test_dir(self):
        elements = {'ConfigError', 'MESA_DEFAULTS', 'PARTS', 'Path',
                    'Release', 'STARTS', 'ScenarioFiles', '__authors__',
                    '__builtins__', '__cached__', '__doc__', '__file__',
                    '__loader__', '__name__', '__package__', '__spec__',
                    '_paths', '_read_ini', '_rerun_files', '_section',
                    'configparser', 'dataclass', 'find_release',
                    'list_reruns', 'mesa_installation', 'release_files'}
        assert set(dir(totest)) == elements

    @pytest.mark.parametrize("name", [
        'find_release', 'list_reruns', 'release_files',
        'mesa_installation'])
    def test_instance_functions(self, name):
        assert isroutine(getattr(totest, name))

    def test_doctests(self):
        assert doctest.testmod(totest).failed == 0


class TestFindRelease:
    def test_inferred_mesa_version(self, repo):
        release = totest.find_release(repo, "v2.0")
        assert release == totest.Release(
            name="v2.0", folder=repo / "r11701" / "releases" / "v2.0",
            mesa_version="r11701", mesa="default", frozen=True)

    def test_given_mesa_version_and_defaults(self, repo):
        folder = repo / "r11701" / "releases" / "dev"
        write(folder / "release.ini", "[release]\n")
        release = totest.find_release(str(repo), "dev", "r11701")
        assert (release.mesa_version, release.mesa, release.frozen) == \
            ("r11701", "default", False)

    def test_missing(self, repo):
        with pytest.raises(ConfigError, match="no release 'v9'"):
            totest.find_release(repo, "v9")

    def test_missing_release_ini(self, repo):
        (repo / "r11701" / "releases" / "bare").mkdir()
        with pytest.raises(ConfigError, match="release.ini not found"):
            totest.find_release(repo, "bare")

    def test_ambiguous(self, repo):
        write(repo / "r24" / "releases" / "v2.0" / "release.ini",
              "[release]\n")
        with pytest.raises(ConfigError, match=r"several MESA versions "
                           r"\(r11701, r24\); set mesa"):
            totest.find_release(repo, "v2.0")


class TestReruns:
    def test_list(self, repo):
        assert totest.list_reruns(totest.find_release(repo, "v2.0")) == \
            ["hot"]

    def test_none(self, repo):
        folder = repo / "r11701" / "releases" / "dev"
        write(folder / "release.ini", "[release]\n")
        assert totest.list_reruns(totest.find_release(repo, "dev")) == []

    def test_unknown_rerun(self, repo):
        release = totest.find_release(repo, "v2.0")
        with pytest.raises(ConfigError, match=r"no rerun 'cold'.*has: hot"):
            totest.release_files(release, "CO-HeMS", "cold")

    def test_unknown_rerun_without_reruns(self, repo):
        folder = repo / "r11701" / "releases" / "dev"
        write(folder / "release.ini", "[release]\n")
        write(folder / "CO-HeMS" / "scenario.ini", "[grid_type]\n")
        with pytest.raises(ConfigError, match="has: none"):
            totest.release_files(totest.find_release(repo, "dev"),
                                 "CO-HeMS", "cold")


class TestReleaseFiles:
    def test_release(self, repo):
        release = totest.find_release(repo, "v2.0")
        folder = release.folder
        files = totest.release_files(release, "CO-HeMS")
        assert files.inlists.targets == {
            "binary_controls": [folder / "common/binary/inlist_project"],
            "binary_job": [folder / "common/binary/inlist_project"],
            "star1_controls": [folder / "common/binary/inlist1",
                               folder / "CO-HeMS/binary/inlist1"],
            "star1_job": [folder / "common/binary/inlist1",
                          folder / "CO-HeMS/binary/inlist1"],
            "star1_formation_every_step_controls":
                [folder / "common/binary/inlist1"],
            "star1_formation_every_step_job":
                [folder / "common/binary/inlist1"],
            "star1_formation_step0_controls":
                [folder / "CO-HeMS/star1_formation/inlist_step1"],
            "star1_formation_step0_job":
                [folder / "CO-HeMS/star1_formation/inlist_step1"]}
        assert files.inlists.whole_files == {
            "history_columns": folder / "common/history_columns.list"}
        assert files.extras == {"posydon_star_binary_extras":
                                folder / "common/binary/src/run_star_extras.f"}
        assert files.zams_star1 is None
        assert files.zams_star2 == repo / "r11701" / "ZAMS_models/zams.data"
        assert files.formation_stars == ("star1",)
        assert files.flags == {}

    def test_rerun(self, repo):
        release = totest.find_release(repo, "v2.0")
        rerun = release.folder / "reruns" / "hot"
        files = totest.release_files(release, "CO-HeMS", "hot")
        assert files.inlists.targets["star1_controls"] == [
            release.folder / "common/binary/inlist1",
            rerun / "common/binary/inlist1",
            release.folder / "CO-HeMS/binary/inlist1"]
        assert files.inlists.targets["star1_formation_every_step_job"] == [
            release.folder / "common/binary/inlist1",
            rerun / "common/binary/inlist1"]
        assert files.extras["posydon_star_binary_extras"] == \
            rerun / "common/binary/src/run_star_extras.f"

    def test_unknown_grid_type(self, repo):
        with pytest.raises(ConfigError, match="scenario.ini not found"):
            totest.release_files(totest.find_release(repo, "v2.0"),
                                 "HMS-HMS")

    def test_missing_file(self, repo):
        release = totest.find_release(repo, "v2.0")
        (release.folder / "CO-HeMS" / "binary" / "inlist1").unlink()
        with pytest.raises(ConfigError, match="names CO-HeMS/binary/inlist1"
                           ", which does not exist"):
            totest.release_files(release, "CO-HeMS")

    def test_bad_start(self, repo):
        release = totest.find_release(repo, "v2.0")
        manifest = release.folder / "CO-HeMS" / "scenario.ini"
        manifest.write_text(GRID_MANIFEST.replace("point_mass", "bh"))
        with pytest.raises(ConfigError, match=r"star2 = 'bh'"):
            totest.release_files(release, "CO-HeMS")

    def test_paths(self):
        assert totest._paths("a, b,, c ") == ["a", "b", "c"]


class TestMesaInstallation:
    def test_defaults(self, repo):
        release = totest.find_release(repo, "v2.0")
        folder = repo / "r11701" / "mesa" / "default" / "defaults"
        assert totest.mesa_installation(release) == {
            name: folder / rel for name, rel in totest.MESA_DEFAULTS.items()}

    def test_missing(self, repo):
        release = totest.find_release(repo, "v2.0")
        release.mesa = "patched"
        with pytest.raises(ConfigError, match="installation 'patched' is "
                           "missing"):
            totest.mesa_installation(release)


def test_paths_are_paths(repo):
    files = totest.release_files(totest.find_release(repo, "v2.0"),
                                 "CO-HeMS")
    assert all(isinstance(path, Path) for paths in
               files.inlists.targets.values() for path in paths)
