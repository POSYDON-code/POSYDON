"""Tests for posydon.CLI.grids.executables."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import os
from pathlib import Path

import pytest

import posydon.CLI.grids.executables as totest
from posydon.utils.posydonwarning import OverwriteWarning


@pytest.fixture(autouse=True)
def _env(tmp_path, monkeypatch):
    """Always provide HOME/MESA_DIR and restore cwd after every test."""
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("MESA_DIR", str(tmp_path / "mesa"))
    original = os.getcwd()
    yield
    os.chdir(original)


class TestCheckFileExist:
    """Tests for check_file_exist()."""

    def test_run_shell(self, monkeypatch):
        """_run_shell forwards the command to os.system."""
        calls = []
        monkeypatch.setattr(totest.os, "system", calls.append)
        totest._run_shell("chmod 755 mk")
        assert calls == ["chmod 755 mk"]

    def test_existing_file_returns_true(self, tmp_path):
        """An existing file returns True."""
        file_path = tmp_path / "present.txt"
        file_path.write_text("present")
        assert totest.check_file_exist(str(file_path)) is True

    def test_missing_file_no_raise(self, tmp_path, capsys):
        """A missing file returns False without raising when raise_error=False."""
        file_path = tmp_path / "absent.txt"
        assert totest.check_file_exist(str(file_path), raise_error=False) is False
        assert "{0} does not exist".format(file_path) in capsys.readouterr().out

    def test_missing_file_raises(self, tmp_path, capsys):
        """A missing file raises ValueError when raise_error=True."""
        file_path = tmp_path / "absent.txt"
        with pytest.raises(ValueError, match="{0} does not exist".format(file_path)):
            totest.check_file_exist(str(file_path))
        assert "{0} does not exist".format(file_path) in capsys.readouterr().out


class TestMakeExecutables:
    """Tests for make_executables()."""

    @staticmethod
    def _extras(tmp_path):
        """Create dummy extra/makefile files and their mesa_extras mapping."""
        extras = tmp_path / "extras"
        extras.mkdir()
        binary_extras = extras / "run_binary_extras.f"
        binary_extras.write_text("! binary extras")
        binary_run = extras / "run_binary.f"
        binary_run.write_text("! binary run")
        star_run = extras / "run.f"
        star_run.write_text("! star run")
        star1_extras = extras / "star1_extras.f"
        star1_extras.write_text("! star1 extras")
        star2_extras = extras / "star2_extras.f"
        star2_extras.write_text("! star2 extras")
        makefile_binary = extras / "makefile_binary"
        makefile_binary.write_text("! binary makefile")
        makefile_star = extras / "makefile_star"
        makefile_star.write_text("! star makefile")
        other = extras / "other.txt"
        other.write_text("other")
        mesa_extras = {
            "mesa_binary_extras": str(binary_extras),
            "mesa_binary_run": str(binary_run),
            "mesa_star_run": str(star_run),
            "star1_extras": str(star1_extras),
            "star2_extras": str(star2_extras),
            "makefile_binary": str(makefile_binary),
            "makefile_star": str(makefile_star),
            "mesa_dir": "/some/mesa",
            "misc_extra": str(other),
        }
        return mesa_extras

    def test_folder_layout_and_mk_script(self, tmp_path, monkeypatch):
        """Creates the folder layout, writes mk, routes extras, returns exe paths."""
        os.chdir(tmp_path)
        mesa_extras = self._extras(tmp_path)

        system_calls = []
        monkeypatch.setattr(totest, "_run_shell", system_calls.append)

        work = tmp_path / "work"
        result = totest.make_executables(mesa_extras, working_directory=str(work))

        for sub in [
            "star1/src",
            "star1/make",
            "star2/src",
            "star2/make",
            "binary/src",
            "binary/make",
        ]:
            assert os.path.isdir(os.path.join(work, sub))

        assert os.path.isfile(
            os.path.join(work, "binary", "src", "run_binary_extras.f")
        )
        assert os.path.isfile(os.path.join(work, "binary", "src", "run_binary.f"))
        assert os.path.isfile(os.path.join(work, "star1", "src", "run.f"))
        assert os.path.isfile(os.path.join(work, "star2", "src", "run.f"))
        assert os.path.isfile(os.path.join(work, "star1", "src", "star1_extras.f"))
        assert os.path.isfile(os.path.join(work, "star2", "src", "star2_extras.f"))
        assert os.path.isfile(os.path.join(work, "binary", "make", "makefile_binary"))
        assert os.path.isfile(os.path.join(work, "star1", "make", "makefile_star"))
        assert os.path.isfile(os.path.join(work, "star2", "make", "makefile_star"))
        assert os.path.isfile(os.path.join(work, "other.txt"))

        mk = Path(tmp_path, "mk").read_text()
        assert "cd {0}\n".format(os.path.join(work, "binary", "make")) in mk
        assert "make -f makefile_binary\n" in mk
        assert "cd {0}\n".format(os.path.join(work, "star1", "make")) in mk
        assert "cd {0}\n".format(os.path.join(work, "star2", "make")) in mk
        assert "make -f makefile_star\n" in mk

        assert system_calls == ["chmod 755 mk", "./mk"]

        assert result == (
            os.path.join(work, "binary", "binary"),
            os.path.join(work, "star1", "star"),
            os.path.join(work, "star2", "star"),
        )

    def test_existing_folders_are_removed(self, tmp_path, monkeypatch):
        """Pre-existing src/make folders are removed and recreated."""
        os.chdir(tmp_path)
        system_calls = []
        monkeypatch.setattr(totest, "_run_shell", system_calls.append)

        work = tmp_path / "work"
        for sub in [
            "star1/src",
            "star1/make",
            "star2/src",
            "star2/make",
            "binary/src",
            "binary/make",
        ]:
            folder = os.path.join(work, sub)
            os.makedirs(folder)
            Path(folder, "stale.txt").write_text("stale")

        totest.make_executables({}, working_directory=str(work))

        for sub in [
            "star1/src",
            "star1/make",
            "star2/src",
            "star2/make",
            "binary/src",
            "binary/make",
        ]:
            assert os.path.isdir(os.path.join(work, sub))
            assert not os.path.exists(os.path.join(work, sub, "stale.txt"))
        assert system_calls == ["chmod 755 mk", "./mk"]

    def test_mk_already_exists_warns(self, tmp_path, monkeypatch):
        """An existing mk file triggers an OverwriteWarning."""
        os.chdir(tmp_path)
        Path("mk").write_text("old")
        system_calls = []
        monkeypatch.setattr(totest, "_run_shell", system_calls.append)

        with pytest.warns(OverwriteWarning, match="Replace mk"):
            totest.make_executables({}, working_directory=str(tmp_path / "work"))
        assert Path("mk").read_text() == ""
        assert system_calls == ["chmod 755 mk", "./mk"]

    def test_none_value_skipped(self, tmp_path, monkeypatch):
        """A None value in mesa_extras is skipped without copying."""
        os.chdir(tmp_path)
        system_calls = []
        monkeypatch.setattr(totest, "_run_shell", system_calls.append)

        work = tmp_path / "work"
        totest.make_executables({"star1_extras": None}, working_directory=str(work))
        assert os.listdir(os.path.join(work, "star1", "src")) == []
        assert system_calls == ["chmod 755 mk", "./mk"]

    def test_empty_mesa_extras(self, tmp_path, monkeypatch):
        """An empty mesa_extras dict writes an empty mk script."""
        os.chdir(tmp_path)
        system_calls = []
        monkeypatch.setattr(totest, "_run_shell", system_calls.append)

        work = tmp_path / "work"
        result = totest.make_executables({}, working_directory=str(work))
        assert Path("mk").read_text() == ""
        assert result == (
            os.path.join(work, "binary", "binary"),
            os.path.join(work, "star1", "star"),
            os.path.join(work, "star2", "star"),
        )
        assert system_calls == ["chmod 755 mk", "./mk"]
