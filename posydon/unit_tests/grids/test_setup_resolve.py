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
        elements = {'ConfigError', 'Path', 'REPOSITORIES',
                    'ResolvedScenario', 'ScenarioFiles', '__authors__',
                    '__builtins__', '__cached__', '__doc__', '__file__',
                    '__loader__', '__name__', '__package__', '__spec__',
                    '_head', 'checkout_release_repo', 'clone_if_missing',
                    'dataclass', 'fetch_inlist_repo', 'find_release',
                    'get_grid_type', 'git', 'inlist_repo_path',
                    'release_files', 'resolve_scenario', 'scenario_files'}
        assert set(dir(totest)) == elements

    @pytest.mark.parametrize("name", ['checkout_release_repo',
                                      'resolve_scenario'])
    def test_instance_functions(self, name):
        assert isroutine(getattr(totest, name))


class TestResolveScenario:
    def test_no_scenario(self):
        assert totest.resolve_scenario(GridConfig()) is None

    def test_needs_grid_type(self, clone):
        config = config_for(clone, release="v2.0", ref="main")
        config.grid.type = None
        with pytest.raises(ConfigError, match=r"\[grid\] type: needed"):
            totest.resolve_scenario(config)

    def test_release_at_branch(self, clone, origin):
        resolved = totest.resolve_scenario(
            config_for(clone, release="v2.0", ref="main"))
        assert resolved.commit == origin["second"]
        assert resolved.repo == clone
        assert resolved.release.name == "v2.0"
        assert resolved.files.inlists.targets["star1_controls"] == [
            clone / "r11701/releases/v2.0/common/binary/inlist1"]
        assert resolved.files.formation_stars == ("star1",)

    def test_release_at_tag_and_commit(self, clone, origin):
        assert totest.resolve_scenario(config_for(
            clone, release="v2.0", ref="v2")).commit == origin["first"]
        assert totest.resolve_scenario(config_for(
            clone, release="v2.0", ref="main",
            commit=origin["first"][:7])).commit == origin["first"]
        assert not (clone / "r11701/releases/v2.0/notes.txt").exists()

    def test_legacy(self, clone, origin):
        resolved = totest.resolve_scenario(config_for(
            clone, branch="main", commit=origin["first"],
            system_type="CO-He_star"))
        assert resolved.commit == origin["first"]
        assert resolved.release is None
        assert "star1_controls" in resolved.files.inlists.targets

    def test_no_release_needs_a_commit(self, clone):
        """Before, git was run with commit None and raised a TypeError."""
        config = config_for(clone, ref="main")
        with pytest.raises(ConfigError, match=r"\[scenario\] release: "
                           "missing"):
            totest.resolve_scenario(config)


class TestCheckout:
    def test_unknown_ref(self, clone):
        with pytest.raises(ConfigError, match=r"\[scenario\] ref: 'nope' "
                           r"is not in .*clone$"):
            totest.checkout_release_repo(
                config_for(clone, ref="nope").scenario)

    def test_unknown_commit_without_network(self, clone, origin):
        shutil.rmtree(origin["path"])
        with pytest.raises(ConfigError, match=r"commit: 'abcdef1' is not "
                           r"in .*\(and fetching failed\)"):
            totest.checkout_release_repo(
                config_for(clone, commit="abcdef1").scenario)

    def test_known_commit_without_network(self, clone, origin):
        shutil.rmtree(origin["path"])
        repo, commit = totest.checkout_release_repo(
            config_for(clone, commit=origin["first"]).scenario)
        assert commit == origin["first"]

    def test_checkout_fails(self, clone, monkeypatch):
        real = totest.git

        def broken(repo, *args, **options):
            if args[0] == "checkout":
                return subprocess.CompletedProcess(args, 1, "", "locked")
            return real(repo, *args, **options)
        monkeypatch.setattr(totest, "git", broken)
        with pytest.raises(ConfigError, match="could not check out .*locked"):
            totest.checkout_release_repo(
                config_for(clone, ref="main").scenario)

    def test_clones_when_missing(self, origin, tmp_path, monkeypatch):
        monkeypatch.setitem(totest.REPOSITORIES, "posydon",
                            (str(origin["path"]), ".inlists"))
        repo, commit = totest.checkout_release_repo(
            config_for(ref="main").scenario, home=tmp_path)
        assert repo == tmp_path / ".inlists"
        assert commit == origin["second"]

    def test_clone_fails(self, tmp_path, monkeypatch):
        monkeypatch.setitem(totest.REPOSITORIES, "posydon",
                            (str(tmp_path / "missing"), ".inlists"))
        with pytest.raises(ConfigError, match="source: could not clone"):
            totest.checkout_release_repo(config_for(ref="main").scenario,
                                         home=tmp_path)


class TestPathSource:
    def test_as_it_is(self, clone, origin):
        repo, commit = totest.checkout_release_repo(
            config_for(clone, source="path").scenario)
        assert (repo, commit) == (clone, origin["second"])

    def test_commit_must_match(self, clone, origin):
        assert totest.checkout_release_repo(config_for(
            clone, source="path",
            commit=origin["second"][:7]).scenario)[1] == origin["second"]
        with pytest.raises(ConfigError, match="commit: .* is at"):
            totest.checkout_release_repo(config_for(
                clone, source="path", commit=origin["first"]).scenario)

    def test_not_a_git_repository(self, tmp_path):
        folder = tmp_path / "plain"
        folder.mkdir()
        assert totest.checkout_release_repo(config_for(
            folder, source="path").scenario) == (folder, None)
        with pytest.raises(ConfigError, match="is at None"):
            totest.checkout_release_repo(config_for(
                folder, source="path", commit="abcdef1").scenario)

    def test_missing(self, tmp_path):
        with pytest.raises(ConfigError, match="inlist_repo: .* does not "
                           "exist"):
            totest.checkout_release_repo(config_for(
                tmp_path / "nope", source="path").scenario)
