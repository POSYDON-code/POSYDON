"""Unit tests of posydon/CLI/grid/main.py"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import sys

import pytest

# import the module which will be tested
import posydon.CLI.grid.main as totest
from posydon.grids.setup.config import ConfigError, SetupError


@pytest.fixture
def calls(monkeypatch):
    """Record what the commands would do instead of doing it."""
    seen = []
    monkeypatch.setattr(totest, "setup_grid", lambda inifile, **options:
                        seen.append(("setup", inifile, options)))
    monkeypatch.setattr(totest, "explain_grid", lambda inifile, **options:
                        seen.append(("explain", inifile, options))
                        or "EXPLAINED")
    return seen


LEGACY = {"grid_type": None, "run_directory": None, "submission_type": None}


class TestElements:
    def test_dir(self):
        assert {'build_parser', 'main', 'fail', 'init_command', 'run_command',
                'build_legacy_parser', 'equivalent_command',
                'legacy_main'} <= set(dir(totest))


class TestPosydonGrid:
    def test_types(self, capsys):
        totest.main(["types"])
        assert capsys.readouterr().out.startswith("HMS-HMS ")

    def test_init(self, tmp_path, capsys):
        ini = tmp_path / "grid.ini"
        totest.main(["init", "CO-HeMS", str(ini)])
        assert "type = CO-HeMS" in ini.read_text()
        assert "Wrote" in capsys.readouterr().out
        with pytest.raises(SetupError, match="exists; not overwriting"):
            totest.init_command("CO-HeMS", str(ini))

    def test_setup_new_format(self, calls):
        totest.main(["setup", "grid.ini"])
        assert calls == [("setup", "grid.ini",
                          dict(LEGACY, overrides=False,
                               command="posydon-grid setup grid.ini"))]

    def test_setup_legacy_options(self, calls):
        totest.main(["setup", "old.ini", "--grid-type", "fixed",
                     "--run-directory", "/run", "--submission-type",
                     "slurm"])
        assert calls[0][2] == {"grid_type": "fixed", "run_directory": "/run",
                               "submission_type": "slurm",
                               "overrides": False, "command":
                               "posydon-grid setup old.ini --grid-type "
                               "fixed --run-directory /run "
                               "--submission-type slurm"}

    def test_explain(self, calls, capsys):
        totest.main(["explain", "grid.ini"])
        assert capsys.readouterr().out == "EXPLAINED\n"
        assert calls[0][0] == "explain"

    def test_run_forwards_everything(self, monkeypatch):
        seen = []
        monkeypatch.setattr(totest, "find_run_grid", lambda: "/bin/run-grid")
        monkeypatch.setattr(totest.os, "execv",
                            lambda path, argv: seen.append((path, argv)))
        totest.main(["run", "--grid-type", "fixed", "--help"])
        assert seen == [(sys.executable, [sys.executable, "/bin/run-grid",
                                          "--grid-type", "fixed",
                                          "--help"])]

    def test_needs_a_command(self):
        with pytest.raises(SystemExit):
            totest.main([])

    def test_reads_the_command_line(self, calls, monkeypatch):
        monkeypatch.setattr(totest.sys, "argv",
                            ["posydon-grid", "setup", "grid.ini"])
        totest.main()
        assert calls[0][:2] == ("setup", "grid.ini")

    def test_rejects_bad_options(self):
        with pytest.raises(SystemExit):
            totest.main(["setup", "grid.ini", "--grid-type", "adaptive"])


class TestPosydonSetupGrid:
    def test_setup_warns_and_keeps_overrides(self, calls):
        with pytest.warns(match=r"posydon-setup-grid is deprecated.*'posydon"
                          r"-grid setup grid.ini --grid-type fixed'"):
            totest.legacy_main(["--inifile", "grid.ini", "--grid-type",
                                "fixed"])
        assert calls == [("setup", "grid.ini",
                          dict(LEGACY, grid_type="fixed", overrides=True,
                               command="posydon-setup-grid --inifile "
                               "grid.ini --grid-type fixed"))]

    def test_explain(self, calls, capsys):
        with pytest.warns(match="'posydon-grid explain grid.ini'"):
            totest.legacy_main(["--inifile", "grid.ini", "--explain"])
        assert capsys.readouterr().out == "EXPLAINED\n"

    def test_types_and_init(self, tmp_path, capsys):
        with pytest.warns(match="'posydon-grid types'"):
            totest.legacy_main(["--list-grid-types"])
        assert "single_HeMS" in capsys.readouterr().out
        ini = tmp_path / "my grid.ini"
        with pytest.warns(match="'posydon-grid init HeMS-HMS '"):
            totest.legacy_main(["--inifile", str(ini), "--init",
                                "HeMS-HMS"])
        assert ini.is_file()

    def test_needs_an_inifile(self):
        with pytest.raises(SystemExit):
            totest.legacy_main(["--grid-type", "fixed"])

    def test_equivalent_command_quotes(self):
        args = totest.build_legacy_parser().parse_args(
            ["--inifile", "my grid.ini", "--run-directory", "/r",
             "--submission-type", "slurm"])
        assert totest.equivalent_command(args) == (
            "posydon-grid setup 'my grid.ini' --run-directory /r "
            "--submission-type slurm")


def raising(error):
    def command(*args, **kwargs):
        raise error
    return command


class TestUserErrors:
    """A user's mistake is one line and exit status 2; a bug raises."""

    def one_line(self, capsys, call):
        with pytest.raises(SystemExit) as stop:
            call()
        assert stop.value.code == 2
        captured = capsys.readouterr()
        assert captured.out == ""
        assert captured.err.count("\n") == 1
        return captured.err

    def test_init_onto_an_existing_file(self, tmp_path, capsys):
        ini = tmp_path / "grid.ini"
        ini.write_text("mine")
        err = self.one_line(capsys, lambda: totest.main(
            ["init", "CO-HeMS", str(ini)]))
        assert err == "posydon-grid: error: {} exists; not overwriting " \
            "it\n".format(ini)
        assert ini.read_text() == "mine"

    def test_unknown_grid_type(self, tmp_path, capsys):
        err = self.one_line(capsys, lambda: totest.main(
            ["init", "MS-MS", str(tmp_path / "grid.ini")]))
        assert err.startswith("posydon-grid: error: [grid] type: unknown "
                              "grid type 'MS-MS'")

    def test_missing_environment(self, tmp_path, capsys, monkeypatch):
        monkeypatch.setenv("MESA_DIR", str(tmp_path))
        monkeypatch.delenv("MESASDK_ROOT", raising=False)
        err = self.one_line(capsys, lambda: totest.main(
            ["setup", str(tmp_path / "grid.ini")]))
        assert err.startswith("posydon-grid: error: MESASDK_ROOT must be "
                              "defined in your environment")

    def test_missing_ini(self, tmp_path, capsys, monkeypatch):
        monkeypatch.setenv("MESA_DIR", str(tmp_path))
        err = self.one_line(capsys, lambda: totest.main(
            ["explain", str(tmp_path / "nope.ini")]))
        assert err == "posydon-grid: error: {} not found\n".format(
            tmp_path / "nope.ini")

    def test_config_error(self, capsys, monkeypatch):
        monkeypatch.setattr(totest, "setup_grid", raising(ConfigError(
            "submission", "backend", "a fixed grid runs one point per "
            "task")))
        err = self.one_line(capsys, lambda: totest.main(
            ["setup", "grid.ini"]))
        assert err == "posydon-grid: error: [submission] backend: a fixed " \
            "grid runs one point per task\n"

    @pytest.mark.parametrize("args", [["run", "--grid-type", "fixed"],
                                      ["run", "-h"]])
    def test_run_without_posydon_run_grid(self, args, capsys, monkeypatch,
                                          tmp_path):
        monkeypatch.setenv("PATH", str(tmp_path))
        err = self.one_line(capsys, lambda: totest.main(args))
        assert err == "posydon-grid: error: Cannot locate posydon-run-grid " \
            "executable in your path\n"

    def test_deprecated_command(self, capsys, monkeypatch):
        monkeypatch.setattr(totest, "explain_grid", raising(ConfigError(
            "grid", "type", "needed")))
        with pytest.warns(match="deprecated"):
            err = self.one_line(capsys, lambda: totest.legacy_main(
                ["--inifile", "grid.ini", "--explain"]))
        assert err == "posydon-setup-grid: error: [grid] type: needed\n"

    @pytest.mark.parametrize("error", [ValueError("a bug"),
                                       KeyError("a bug"),
                                       RuntimeError("mk failed")])
    def test_other_errors_raise(self, error, monkeypatch):
        monkeypatch.setattr(totest, "setup_grid", raising(error))
        monkeypatch.setattr(totest, "explain_grid", raising(error))
        with pytest.raises(type(error)):
            totest.main(["setup", "grid.ini"])
        with pytest.raises(type(error)), pytest.warns(match="deprecated"):
            totest.legacy_main(["--inifile", "grid.ini", "--explain"])
