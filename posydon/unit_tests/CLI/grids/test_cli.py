"""Tests for the bin/posydon-grid CLI and its deprecated wrappers."""

import importlib.machinery
import importlib.util
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[4]


def _load_script(name, path):
    loader = importlib.machinery.SourceFileLoader(name, str(path))
    spec = importlib.util.spec_from_loader(name, loader)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def posydon_grid():
    """Load bin/posydon-grid once per module."""
    return _load_script("posydon_grid_cli", REPO_ROOT / "bin" / "posydon-grid")


@pytest.fixture(scope="module")
def setup_wrapper():
    """Load bin/posydon-setup-grid once per module."""
    return _load_script("posydon_setup_grid_wrapper", REPO_ROOT / "bin" / "posydon-setup-grid")


@pytest.fixture(scope="module")
def run_wrapper():
    """Load bin/posydon-run-grid once per module."""
    return _load_script("posydon_run_grid_wrapper", REPO_ROOT / "bin" / "posydon-run-grid")


class TestPosydonGridSetupSubcommand:
    """Tests for the setup subcommand of bin/posydon-grid."""

    def test_valid_fixed(self, posydon_grid, tmp_path):
        args = posydon_grid.parse_commandline(
            ["setup", "--inifile", "in.ini", "--grid-type", "fixed"])
        assert args.command == "setup"
        assert args.inifile == "in.ini"
        assert args.grid_type == "fixed"

    def test_valid_dynamic(self, posydon_grid):
        args = posydon_grid.parse_commandline(
            ["setup", "--inifile", "in.ini", "--grid-type", "dynamic"])
        assert args.command == "setup"
        assert args.grid_type == "dynamic"

    def test_defaults(self, posydon_grid, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        args = posydon_grid.parse_commandline(
            ["setup", "--inifile", "in.ini", "--grid-type", "fixed"])
        assert args.run_directory == str(tmp_path)
        assert args.submission_type == "shell"
        assert args.nproc == 1
        assert args.verbose is False

    def test_all_flags(self, posydon_grid, tmp_path):
        args = posydon_grid.parse_commandline(
            ["setup", "--inifile", "in.ini", "--grid-type", "dynamic",
             "--run-directory", "/tmp/run", "--submission-type", "slurm",
             "-n", "4", "-v"])
        assert args.run_directory == "/tmp/run"
        assert args.submission_type == "slurm"
        assert args.nproc == 4
        assert args.verbose is True

    def test_invalid_grid_type(self, posydon_grid):
        with pytest.raises(SystemExit):
            posydon_grid.parse_commandline(
                ["setup", "--inifile", "in.ini", "--grid-type", "bogus"])

    def test_invalid_submission_type(self, posydon_grid):
        with pytest.raises(SystemExit):
            posydon_grid.parse_commandline(
                ["setup", "--inifile", "in.ini", "--grid-type", "fixed",
                 "--submission-type", "bogus"])


class TestPosydonGridRunSubcommand:
    """Tests for the run subcommand of bin/posydon-grid."""

    def test_valid_fixed(self, posydon_grid, tmp_path):
        args = posydon_grid.parse_commandline(
            ["run", "--mesa-grid", "grid.csv", "--grid-type", "fixed",
             "--output-directory", "/tmp/out", "--mesa-binary-executable", "/exe/bin",
             "--mesa-binary-inlist-project", "/inlist/project",
             "--mesa-binary-inlist1", "/inlist/1",
             "--mesa-binary-inlist2", "/inlist/2",
             "--mesa-star-history-columns", "/cols/history",
             "--mesa-binary-history-columns", "/cols/bhistory",
             "--mesa-profile-columns", "/cols/profile"])
        assert args.command == "run"
        assert args.mesa_grid == "grid.csv"
        assert args.grid_type == "fixed"
        assert args.output_directory == "/tmp/out"

    def test_defaults(self, posydon_grid, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        args = posydon_grid.parse_commandline(
            ["run", "--mesa-grid", "grid.csv", "--grid-type", "fixed",
             "--output-directory", "/tmp/out", "--mesa-binary-executable", "/exe/bin",
             "--mesa-binary-inlist-project", "/inlist/project",
             "--mesa-binary-inlist1", "/inlist/1",
             "--mesa-binary-inlist2", "/inlist/2",
             "--mesa-star-history-columns", "/cols/history",
             "--mesa-binary-history-columns", "/cols/bhistory",
             "--mesa-profile-columns", "/cols/profile"])
        assert args.temporary_directory == str(tmp_path)
        assert args.grid_point_index is None
        assert args.Niter == 200
        assert args.job_before_end == 300
        assert args.keep_profiles is False
        assert args.keep_photos is False

    def test_all_flags(self, posydon_grid):
        args = posydon_grid.parse_commandline(
            ["run", "--mesa-grid", "grid.csv", "--grid-type", "dynamic",
             "--output-directory", "/tmp/out", "--temporary-directory", "/tmp/tmp",
             "--grid-point-index", "3", "--psycris-inifile", "/tmp/psy.ini",
             "--Niter", "50", "--mesa-binary-executable", "/exe/bin",
             "--mesa-star1-executable", "/exe/star1",
             "--mesa-star2-executable", "/exe/star2",
             "--mesa-binary-inlist-project", "/inlist/project",
             "--mesa-binary-inlist1", "/inlist/1",
             "--mesa-binary-inlist2", "/inlist/2",
             "--mesa-star1-inlist-project", "/inlist/star1a", "/inlist/star1b",
             "--mesa-star2-inlist-project", "/inlist/star2a",
             "--mesa-star-history-columns", "/cols/history",
             "--mesa-binary-history-columns", "/cols/bhistory",
             "--mesa-profile-columns", "/cols/profile",
             "--verbose", "--keep_profiles", "--keep_photos",
             "--job_end", "1000", "--job_before_end", "10"])
        assert args.temporary_directory == "/tmp/tmp"
        assert args.grid_point_index == 3
        assert args.psycris_inifile == "/tmp/psy.ini"
        assert args.Niter == 50
        assert args.mesa_star1_inlist_project == ["/inlist/star1a", "/inlist/star1b"]
        assert args.mesa_star2_inlist_project == ["/inlist/star2a"]
        assert args.verbose is True
        assert args.keep_profiles is True
        assert args.keep_photos is True
        assert args.job_end == 1000
        assert args.job_before_end == 10

    def test_invalid_grid_type(self, posydon_grid):
        with pytest.raises(SystemExit):
            posydon_grid.parse_commandline(
                ["run", "--mesa-grid", "grid.csv", "--grid-type", "bogus",
                 "--output-directory", "/tmp/out", "--mesa-binary-executable", "/exe/bin",
                 "--mesa-binary-inlist-project", "/inlist/project",
                 "--mesa-binary-inlist1", "/inlist/1",
                 "--mesa-binary-inlist2", "/inlist/2",
                 "--mesa-star-history-columns", "/cols/history",
                 "--mesa-binary-history-columns", "/cols/bhistory",
                 "--mesa-profile-columns", "/cols/profile"])


class TestPosydonGridMain:
    """Tests for the main() dispatch of bin/posydon-grid."""

    def test_no_subcommand_prints_help(self, posydon_grid, capsys, monkeypatch):
        monkeypatch.setattr(sys, "argv", ["posydon-grid"])
        with pytest.raises(SystemExit) as exc:
            posydon_grid.main()
        assert exc.value.code == 1
        out = capsys.readouterr().out
        assert "usage: posydon-grid" in out

    def test_setup_dispatches(self, posydon_grid, monkeypatch):
        called = []
        monkeypatch.setattr(
            posydon_grid.setup, "main",
            lambda: called.append("setup"))
        monkeypatch.setattr(
            sys, "argv",
            ["posydon-grid", "setup", "--inifile", "in.ini", "--grid-type", "fixed"])
        posydon_grid.main()
        assert called == ["setup"]

    def test_run_dispatches(self, posydon_grid, monkeypatch):
        called = []
        monkeypatch.setattr(
            posydon_grid.run, "main",
            lambda: called.append("run"))
        monkeypatch.setattr(
            sys, "argv",
            ["posydon-grid", "run", "--mesa-grid", "grid.csv", "--grid-type", "fixed",
             "--output-directory", "/tmp/out", "--mesa-binary-executable", "/exe/bin",
             "--mesa-binary-inlist-project", "/inlist/project",
             "--mesa-binary-inlist1", "/inlist/1",
             "--mesa-binary-inlist2", "/inlist/2",
             "--mesa-star-history-columns", "/cols/history",
             "--mesa-binary-history-columns", "/cols/bhistory",
             "--mesa-profile-columns", "/cols/profile"])
        posydon_grid.main()
        assert called == ["run"]


class TestSetupWrapperDeprecation:
    """Tests for the deprecated bin/posydon-setup-grid wrapper."""

    def test_emits_deprecation_warning(self, setup_wrapper, monkeypatch):
        warned = []
        monkeypatch.setattr(
            setup_wrapper,
            "Pwarn",
            lambda *a, **k: warned.append((a, k)),
        )
        monkeypatch.setattr(
            setup_wrapper.setup,
            "parse_commandline",
            lambda: ("args",),
        )
        ran = []
        monkeypatch.setattr(
            setup_wrapper.setup,
            "run_setup",
            lambda args: ran.append(args),
        )
        monkeypatch.setattr(sys, "argv", ["posydon-setup-grid"])
        setup_wrapper.main()
        assert len(warned) == 1
        message, category = warned[0][0]
        assert "posydon-setup-grid is deprecated" in message
        assert category == "DeprecationWarning"
        assert ran == [("args",)]


class TestRunWrapperDeprecation:
    """Tests for the deprecated bin/posydon-run-grid wrapper."""

    def test_emits_deprecation_warning(self, run_wrapper, monkeypatch):
        warned = []
        monkeypatch.setattr(
            run_wrapper,
            "Pwarn",
            lambda *a, **k: warned.append((a, k)),
        )
        ran = []
        monkeypatch.setattr(
            run_wrapper.run,
            "main",
            lambda: ran.append("run"),
        )
        monkeypatch.setattr(sys, "argv", ["posydon-run-grid"])
        run_wrapper.main()
        assert len(warned) == 1
        message, category = warned[0][0]
        assert "posydon-run-grid is deprecated" in message
        assert category == "DeprecationWarning"
        assert ran == ["run"]
