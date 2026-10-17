"""Unit tests of posydon/CLI/grid/main.py"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import sys

import pytest

# import the module which will be tested
import posydon.CLI.grid.main as totest


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
        assert {'build_parser', 'main', 'init_command', 'run_command',
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
        with pytest.raises(ValueError, match="exists; not overwriting"):
            totest.main(["init", "CO-HeMS", str(ini)])

    def test_setup_new_format(self, calls):
        totest.main(["setup", "grid.ini"])
        assert calls == [("setup", "grid.ini",
                          dict(LEGACY, overrides=False))]

    def test_setup_legacy_options(self, calls):
        totest.main(["setup", "old.ini", "--grid-type", "fixed",
                     "--run-directory", "/run", "--submission-type",
                     "slurm"])
        assert calls[0][2] == {"grid_type": "fixed", "run_directory": "/run",
                               "submission_type": "slurm",
                               "overrides": False}

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
                          dict(LEGACY, grid_type="fixed", overrides=True))]

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
