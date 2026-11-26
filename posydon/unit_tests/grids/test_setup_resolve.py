"""Unit tests of posydon/grids/setup/resolve.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

# import other needed code for the tests, which is not already imported in the
# module you like to test
import shutil
import subprocess
from inspect import isroutine

import pytest

# import the module which will be tested
import posydon.grids.setup.resolve as totest
from posydon.grids.setup.config import ConfigError, GridConfig
from posydon.grids.setup.scenario import ScenarioFiles

GIT_ENV = {"GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t",
           "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@t"}
MANIFEST = """\
[grid_type]
star1 = formation
star2 = point_mass

[inlists]
star1 = common/binary/inlist1
"""

pytestmark = pytest.mark.skipif(shutil.which("git") is None,
                                reason="needs git")


def git(repo, *args):
    result = subprocess.run(["git", "-C", str(repo), *args],
                            capture_output=True, text=True,
                            env={**GIT_ENV, "PATH": "/usr/bin:/bin"})
    assert result.returncode == 0, result.stderr
    return result.stdout.strip()


def write(path, text=""):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


@pytest.fixture
def origin(tmp_path):
    """A repository with release v2.0 (commit 'first', tag v2) and a
    second commit on main."""
    repo = tmp_path / "origin"
    release = repo / "r11701" / "releases" / "v2.0"
    write(release / "release.ini", "[release]\n")
    write(release / "common" / "binary" / "inlist1", "&controls\n/\n")
    write(release / "CO-HeMS" / "scenario.ini", MANIFEST)
    git(tmp_path, "init", "-q", "-b", "main", str(repo))
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "first")
    git(repo, "tag", "v2")
    first = git(repo, "rev-parse", "HEAD")
    write(release / "notes.txt", "second")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "second")
    return {"path": repo, "first": first,
            "second": git(repo, "rev-parse", "HEAD")}


@pytest.fixture
def clone(origin, tmp_path):
    path = tmp_path / "clone"
    git(tmp_path, "clone", "-q", str(origin["path"]), str(path))
    return path


def config_for(repo=None, **scenario):
    config = GridConfig()
    config.grid.type = "CO-HeMS"
    config.scenario.inlist_repo = repo
    for key, value in scenario.items():
        setattr(config.scenario, key, value)
    return config


class TestElements:
    def test_dir(self):
        elements = {'ConfigError', 'MESA_DEFAULTS', 'Path', 'REPOSITORIES',
                    'ResolvedScenario', 'ScenarioFiles', '__authors__',
                    '__builtins__', '__cached__', '__doc__', '__file__',
                    '__loader__', '__name__', '__package__', '__spec__',
                    '_relative', '_resolved', 'dataclass',
                    'fetch_inlist_repo', 'files_used', 'find_release',
                    'fingerprint', 'get_grid_type', 'git', 'git_state',
                    'hashlib', 'inlist_repo_path', 'release_files',
                    'release_repo', 'resolve_scenario', 'scenario_files'}
        assert set(dir(totest)) == elements

    @pytest.mark.parametrize("name", ['release_repo', 'resolve_scenario',
                                      'files_used', 'fingerprint',
                                      'git_state'])
    def test_instance_functions(self, name):
        assert isroutine(getattr(totest, name))


class TestResolveScenario:
    def test_no_scenario(self):
        assert totest.resolve_scenario(GridConfig()) is None

    def test_needs_grid_type(self, clone):
        config = config_for(clone, release="v2.0")
        config.grid.type = None
        with pytest.raises(ConfigError, match=r"\[grid\] type: needed"):
            totest.resolve_scenario(config)

    def test_release_as_it_is(self, clone, origin):
        """Whatever the clone has checked out is read; nothing changes."""
        git(clone, "checkout", "-q", origin["first"])
        resolved = totest.resolve_scenario(config_for(clone,
                                                      release="v2.0"))
        assert resolved.commit == origin["first"]
        assert resolved.dirty is False
        assert resolved.remote == str(origin["path"])
        assert resolved.repo == clone
        assert resolved.release.name == "v2.0"
        assert resolved.files.inlists.targets["star1_controls"] == [
            clone / "r11701/releases/v2.0/common/binary/inlist1"]
        assert resolved.files.formation_stars == ("star1",)
        # release.ini, the manifest and inlist1
        assert resolved.n_files == 3
        assert git(clone, "rev-parse", "HEAD") == origin["first"]

    def test_default_folder(self, origin, tmp_path):
        git(tmp_path, "clone", "-q", str(origin["path"]),
            str(tmp_path / ".posydon_mesa_inlists"))
        resolved = totest.resolve_scenario(config_for(release="v2.0"),
                                           home=tmp_path)
        assert resolved.repo == tmp_path / ".posydon_mesa_inlists"

    def test_legacy(self, clone, origin):
        """A legacy scenario still checks out its commit."""
        resolved = totest.resolve_scenario(config_for(
            clone, branch="main", commit=origin["first"],
            system_type="CO-He_star"))
        assert resolved.commit == origin["first"]
        assert git(clone, "rev-parse", "HEAD") == origin["first"]
        assert resolved.release is None
        assert resolved.fingerprint is not None
        assert "star1_controls" in resolved.files.inlists.targets

    def test_no_release_needs_a_commit(self, clone):
        """Before, git was run with commit None and raised a TypeError."""
        config = config_for(clone, branch="main")
        with pytest.raises(ConfigError, match=r"\[scenario\] release: "
                           "missing"):
            totest.resolve_scenario(config)


class TestReleaseRepo:
    def test_missing(self, tmp_path):
        """The setup no longer clones: it says how to."""
        scenario = config_for(tmp_path / "nothing").scenario
        with pytest.raises(ConfigError, match=r"\[scenario\] inlist_repo: "
                           ".*nothing does not exist(.|\n)*git clone "
                           "--branch v2_sources https://github.com/"):
            totest.release_repo(scenario)

    def test_not_a_git_repository(self, tmp_path):
        folder = tmp_path / "copy"
        folder.mkdir()
        assert totest.release_repo(config_for(folder).scenario) == folder


class TestProvenance:
    def test_files_used(self, clone):
        files = ScenarioFiles(zams_star1=clone / "z1")
        files.inlists.add("star1_controls", clone / "a")
        files.inlists.add("star1_job", clone / "a")
        files.inlists.whole_files["history_columns"] = clone / "h"
        files.extras["user_star1_extras"] = clone / "e"
        assert totest.files_used(files) == [clone / name for name in
                                            ("a", "e", "h", "z1")]

    def test_fingerprint(self, tmp_path):
        for folder in ("one", "two"):
            write(tmp_path / folder / "a", "same")
        one, two = (totest.fingerprint(tmp_path / folder,
                                       [tmp_path / folder / "a"])
                    for folder in ("one", "two"))
        # the same files in another folder: the same fingerprint
        assert one == two and len(one) == 64
        write(tmp_path / "two" / "a", "changed")
        assert totest.fingerprint(tmp_path / "two",
                                  [tmp_path / "two" / "a"]) != one
        # a missing file counts, as missing
        assert totest.fingerprint(tmp_path / "two",
                                  [tmp_path / "two" / "gone"]) != one
        # outside the repository, the full path counts
        assert totest.fingerprint(tmp_path / "other",
                                  [tmp_path / "one" / "a"]) != one

    def test_dirty(self, clone):
        used = clone / "r11701/releases/v2.0/common/binary/inlist1"
        assert totest.git_state(clone, [used])[1] is False
        write(clone / "unused.txt", "x")
        assert totest.git_state(clone, [used])[1] is False
        write(used, "&controls\n x = 1\n/\n")
        assert totest.git_state(clone, [used])[1] is True

    def test_not_a_repository(self, tmp_path):
        assert totest.git_state(tmp_path, []) == (None, None, None)

    def test_inside_another_repository(self, clone):
        folder = clone / "r11701"
        assert totest.git_state(folder, []) == (None, None, None)

    def test_no_remote(self, origin):
        commit, dirty, remote = totest.git_state(origin["path"], [])
        assert (commit, dirty, remote) == (origin["second"], False, None)

    def test_without_git(self, clone, monkeypatch):
        def missing(*args, **options):
            raise FileNotFoundError("git")
        monkeypatch.setattr(totest, "git", missing)
        assert totest.git_state(clone, []) == (None, None, None)
