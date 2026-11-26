"""Unit tests of posydon/grids/setup/support.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

# import other needed code for the tests, which is not already imported in the
# module you like to test
import doctest
import shutil
import subprocess
from inspect import isroutine

import pytest

# import the module which will be tested
import posydon.grids.setup.support as totest
from posydon.grids.setup.config import ConfigError

GIT_ENV = {"GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t",
           "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@t",
           "PATH": "/usr/bin:/bin"}


def git(repo, *args):
    result = subprocess.run(["git", "-C", str(repo), *args],
                            capture_output=True, text=True, env=GIT_ENV)
    assert result.returncode == 0, result.stderr
    return result.stdout.strip()


@pytest.fixture
def origin(tmp_path):
    """A repository with every support path, one commit."""
    repo = tmp_path / "origin"
    for path in totest.SUPPORT_PATHS:
        (repo / path).mkdir(parents=True)
        (repo / path / "file.txt").write_text(path)
    git(tmp_path, "init", "-q", str(repo))
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "support")
    return repo, git(repo, "rev-parse", "HEAD")


class TestElements:
    def test_dir(self):
        assert {'GRID_PARAMS', 'OLD_FOLDERS', 'SUPPORT_COMMIT',
                'SUPPORT_PATHS', 'map_grid_params_path',
                'support_root'} <= set(dir(totest))

    @pytest.mark.parametrize("name", ['map_grid_params_path',
                                      'support_root'])
    def test_instance_functions(self, name):
        assert isroutine(getattr(totest, name))

    def test_doctests(self):
        result = doctest.testmod(totest)
        assert result.failed == 0 and result.attempted > 0


class TestMap:
    @pytest.mark.parametrize("path, expected", [
        ("/p/grid_params/defaults/r11701/binary/binary_job.defaults",
         "r11701/mesa/default/defaults/binary/binary_job.defaults"),
        ("/p/grid_params/defaults/r11701/history_columns.list",
         "r11701/mesa/default/columns/history_columns.list"),
        ("/p/grid_params/defaults/r11701/src/run.f",
         "r11701/mesa/default/src/run.f"),
        ("/p/grid_params/pipeline_quest.ini",
         "templates/pipeline_quest.ini"),
        ("/p/grid_params/POSYDON-MESA-INLISTS-public/r11701/"
         "default_common_inlists/history_columns.list",
         "r11701/releases/v2.0/common/history_columns.list"),
        ("/p/grid_params/POSYDON-MESA-INLISTS/r11701/CO-He_star/binary/"
         "inlist1", "r11701/releases/v2.0/CO-HeMS/binary/inlist1"),
        ("/p/grid_params/POSYDON-MESA-INLISTS/r11701/running_scripts/a.ini",
         "r11701/running_scripts/a.ini"),
        ("grid_params/defaults/r11701/star/controls.defaults",
         "r11701/mesa/default/defaults/star/controls.defaults"),
        ("/p/elsewhere/controls.defaults", None),
        ("/p/grid_params/defaults/r11701/notes.txt", None),
        ("/p/grid_params/defaults/r24/star/controls.defaults", None),
        ("/p/grid_params/POSYDON-MESA-INLISTS/r11701/single_He_star/x",
         None),
        ("/p/grid_params/POSYDON-MESA-INLISTS/README.md", None),
    ])
    def test_paths(self, path, expected):
        assert totest.map_grid_params_path(path) == expected


@pytest.mark.skipif(shutil.which("git") is None, reason="needs git")
class TestSupportRoot:
    def test_extracts_once(self, origin, tmp_path):
        repo, commit = origin
        clone = tmp_path / "clone"
        git(tmp_path, "clone", "-q", str(repo), str(clone))
        root = totest.support_root(home=tmp_path, clone=clone, commit=commit)
        assert root == tmp_path / ".posydon_mesa_inlists_support" / commit
        assert (root / "r11701/mesa/file.txt").read_text() == "r11701/mesa"
        assert (root / "templates/file.txt").is_file()
        # the cache is used as it is afterwards
        shutil.rmtree(clone)
        assert totest.support_root(home=tmp_path, clone=clone,
                                   commit=commit) == root

    def test_fetches_a_missing_commit(self, origin, tmp_path):
        repo, first = origin
        clone = tmp_path / "clone"
        git(tmp_path, "clone", "-q", str(repo), str(clone))
        (repo / "templates" / "file.txt").write_text("second")
        git(repo, "commit", "-q", "-am", "second")
        second = git(repo, "rev-parse", "HEAD")
        root = totest.support_root(home=tmp_path, clone=clone, commit=second)
        assert (root / "templates/file.txt").read_text() == "second"

    def test_clones_the_default(self, origin, tmp_path, monkeypatch):
        repo, commit = origin
        monkeypatch.setitem(totest.REPOSITORIES, "posydon",
                            (str(repo), ".inlists"))
        root = totest.support_root(home=tmp_path, commit=commit)
        assert (tmp_path / ".inlists").is_dir()
        assert (root / "r11701/ZAMS_models/file.txt").is_file()

    def test_clone_fails(self, tmp_path, monkeypatch):
        monkeypatch.setitem(totest.REPOSITORIES, "posydon",
                            (str(tmp_path / "missing"), ".inlists"))
        with pytest.raises(ConfigError, match="posydon_github_root: "
                           "grid_params/ is no longer part of POSYDON, and "
                           "could not clone"):
            totest.support_root(home=tmp_path)

    def test_unknown_commit(self, origin, tmp_path):
        repo, _ = origin
        with pytest.raises(ConfigError, match="does not have commit 0000000"):
            totest.support_root(home=tmp_path, clone=repo, commit="0000000")

    def test_lost_race(self, origin, tmp_path, monkeypatch):
        """Another setup extracted the cache while this one was busy."""
        repo, commit = origin
        cache = tmp_path / ".posydon_mesa_inlists_support" / commit
        real_mkdtemp = totest.tempfile.mkdtemp

        def mkdtemp(dir):
            staging = real_mkdtemp(dir=dir)
            cache.mkdir()
            (cache / "theirs").write_text("")
            return staging
        monkeypatch.setattr(totest.tempfile, "mkdtemp", mkdtemp)
        root = totest.support_root(home=tmp_path, clone=repo, commit=commit)
        assert root == cache and (cache / "theirs").is_file()
        assert sorted(p.name for p in cache.parent.iterdir()) == [commit]
