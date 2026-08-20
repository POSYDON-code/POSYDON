"""Unit tests for posydon.CLI.grids.setup."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import os
import sys
from types import SimpleNamespace

import pytest

from posydon.CLI.grids import setup

totest = setup


@pytest.fixture(autouse=True)
def _env(tmp_path, monkeypatch):
    """Always provide HOME/MESA_DIR and restore cwd after every test."""
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("MESA_DIR", str(tmp_path / "mesa"))
    original = os.getcwd()
    yield
    os.chdir(original)


class TestParseCommandline:
    """Tests for parse_commandline()."""

    @pytest.fixture(autouse=True)
    def _reset_argv(self, monkeypatch):
        monkeypatch.setattr(sys, "argv", ["posydon-setup-grid"])

    def test_valid_fixed_shell(self, monkeypatch):
        """Valid fixed grid with shell submission parses successfully."""
        monkeypatch.setattr(
            sys,
            "argv",
            ["posydon-setup-grid", "--inifile", "in.ini", "--grid-type", "fixed"],
        )
        args = totest.parse_commandline()
        assert args.inifile == "in.ini"
        assert args.grid_type == "fixed"

    def test_valid_dynamic_slurm(self, monkeypatch):
        """Valid dynamic grid with slurm submission parses successfully."""
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "posydon-setup-grid",
                "--inifile",
                "in.ini",
                "--grid-type",
                "dynamic",
                "--submission-type",
                "slurm",
                "-n",
                "4",
                "--run-directory",
                "/some/where",
                "--verbose",
            ],
        )
        args = totest.parse_commandline()
        assert args.grid_type == "dynamic"
        assert args.submission_type == "slurm"
        assert args.nproc == 4
        assert args.run_directory == "/some/where"
        assert args.verbose is True

    def test_invalid_grid_type_raises(self, monkeypatch):
        """An unknown grid type triggers parser.error (SystemExit)."""
        monkeypatch.setattr(
            sys,
            "argv",
            ["posydon-setup-grid", "--inifile", "in.ini", "--grid-type", "bogus"],
        )
        with pytest.raises(SystemExit):
            totest.parse_commandline()

    def test_invalid_submission_type_raises(self, monkeypatch):
        """An unknown submission type triggers parser.error (SystemExit)."""
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "posydon-setup-grid",
                "--inifile",
                "in.ini",
                "--grid-type",
                "fixed",
                "--submission-type",
                "bogus",
            ],
        )
        with pytest.raises(SystemExit):
            totest.parse_commandline()

    def test_defaults(self, monkeypatch, tmp_path):
        """run_directory defaults to cwd, submission_type to shell, nproc to 1."""
        os.chdir(tmp_path)
        monkeypatch.setattr(
            sys,
            "argv",
            ["posydon-setup-grid", "--inifile", "in.ini", "--grid-type", "fixed"],
        )
        args = totest.parse_commandline()
        assert args.run_directory == os.getcwd()
        assert args.submission_type == "shell"
        assert args.nproc == 1
        assert args.verbose is False


class TestMainFlow:
    """Tests for main() covering the full setup-grid flow."""

    @staticmethod
    def _which_popen(communicate_result):
        """Build a Popen stand-in returning a fixed communicate() result."""

        class FakePopen:
            def __init__(self, cmd, **kwargs):
                self.cmd = cmd

            def communicate(self):
                return communicate_result

        return FakePopen

    @pytest.fixture
    def grid_inputs(self, tmp_path):
        """Create grid.csv, the three column-list files and a dummy extras file."""
        grid = tmp_path / "grid.csv"
        grid.write_text("m1,m2\n10.0,8.0\n11.0,9.0\n")
        star_hist = tmp_path / "star_history_columns.list"
        star_hist.write_text("star columns\n")
        binary_hist = tmp_path / "binary_history_columns.list"
        binary_hist.write_text("binary columns\n")
        profile_hist = tmp_path / "profile_columns.list"
        profile_hist.write_text("profile columns\n")
        extras = tmp_path / "binary_extras.f"
        extras.write_text("! dummy extras\n")
        return grid, star_hist, binary_hist, profile_hist, extras

    @staticmethod
    def _mock_executables_and_inlists(monkeypatch):
        """Patch make_executables and construct_static_inlist for the happy path."""
        monkeypatch.setattr(
            totest,
            "make_executables",
            lambda mesa_extras, working_directory: (
                "/work/binary/binary",
                "/work/star1/star",
                "/work/star2/star",
            ),
        )
        monkeypatch.setattr(
            totest,
            "construct_static_inlist",
            lambda mesa_inlists, grid_parameters, working_directory: (
                "/work/star1/inlist_step0",
                None,
                "/work/binary/inlist_project",
                "/work/binary/inlist1",
                "/work/binary/inlist2",
            ),
        )

    def test_main_mesa_dir_missing(self, monkeypatch):
        """main raises when MESA_DIR is not defined in the environment."""
        monkeypatch.delenv("MESA_DIR", raising=False)
        monkeypatch.setattr(
            sys,
            "argv",
            ["posydon-setup-grid", "--inifile", "in.ini", "--grid-type", "fixed"],
        )
        with pytest.raises(ValueError, match="MESA_DIR must be defined"):
            totest.main()

    def test_main_which_empty(self, monkeypatch):
        """Missing posydon-run-grid executable raises ValueError."""
        monkeypatch.setattr(
            sys,
            "argv",
            ["posydon-setup-grid", "--inifile", "in.ini", "--grid-type", "fixed"],
        )
        monkeypatch.setattr(
            totest.subprocess,
            "Popen",
            self._which_popen((b"", b"")),
        )
        with pytest.raises(
            ValueError, match="Cannot locate posydon-run-grid executable in your path"
        ):
            totest.main()

    def test_main_fixed_shell_happy_path(self, monkeypatch, tmp_path, grid_inputs):
        """Fixed grid shell submission writes a runnable grid_command.sh."""
        grid, star_hist, binary_hist, profile_hist, extras = grid_inputs
        os.chdir(tmp_path)
        monkeypatch.setenv("MESASDK_ROOT", "/opt/mesasdk")
        monkeypatch.setattr(
            sys,
            "argv",
            ["posydon-setup-grid", "--inifile", "in.ini", "--grid-type", "fixed"],
        )
        monkeypatch.setattr(
            totest.subprocess,
            "Popen",
            self._which_popen((b"/usr/bin/posydon-run-grid\n", b"")),
        )
        monkeypatch.setattr(
            totest,
            "parse_config_file",
            lambda inifile: (
                {"grid": str(grid), "psycris_inifile": "/tmp/psycris.ini"},
                {
                    "job_array": False,
                    "number_of_mpi_tasks": 1,
                    "number_of_nodes": 1,
                    "number_of_cpus_per_task": 1,
                    "work_dir": "",
                },
                {
                    "scenario": None,
                    "star_history_columns": str(star_hist),
                    "binary_history_columns": str(binary_hist),
                    "profile_columns": str(profile_hist),
                },
                {"mesa_binary_extras": str(extras)},
            ),
        )
        self._mock_executables_and_inlists(monkeypatch)
        system_calls = []
        monkeypatch.setattr(totest.os, "system", system_calls.append)

        totest.main()

        grid_command = (tmp_path / "grid_command.sh").read_text()
        assert "#!/bin/bash" in grid_command
        assert "export OMP_NUM_THREADS=1" in grid_command
        assert "export MESASDK_ROOT=/opt/mesasdk" in grid_command
        assert "export MESA_DIR={0}".format(os.environ["MESA_DIR"]) in grid_command
        assert "compress-mesa ." in grid_command
        assert 'echo "Done."' in grid_command

        expected_command = (
            "python /usr/bin/posydon-run-grid --mesa-grid {0}/grid.csv "
            "--mesa-binary-executable /work/binary/binary "
            "--mesa-star1-executable /work/star1/star "
            "--mesa-star2-executable /work/star2/star "
            "--mesa-binary-inlist-project /work/binary/inlist_project "
            "--mesa-binary-inlist1 /work/binary/inlist1 "
            "--mesa-binary-inlist2 /work/binary/inlist2 "
            "--mesa-star1-inlist-project /work/star1/inlist_step0 "
            "--mesa-star2-inlist-project None "
            "--mesa-star-history-columns {0}/column_lists/history_columns.list "
            "--mesa-binary-history-columns {0}/column_lists/binary_history_columns.list "
            "--mesa-profile-columns {0}/column_lists/profile_columns.list "
            "--output-directory {0} --grid-type fixed "
            "--psycris-inifile /tmp/psycris.ini"
        ).format(tmp_path)
        assert expected_command in grid_command
        assert system_calls == ["chmod 755 grid_command.sh"]
        assert (
            tmp_path / "column_lists" / "history_columns.list"
        ).read_text() == "star columns\n"

    def test_main_fixed_shell_job_array(self, monkeypatch, tmp_path, grid_inputs):
        """Fixed grid shell submission with a job array loops over grid rows."""
        grid, star_hist, binary_hist, profile_hist, extras = grid_inputs
        os.chdir(tmp_path)
        monkeypatch.setenv("MESASDK_ROOT", "/opt/mesasdk")
        monkeypatch.setattr(
            sys,
            "argv",
            ["posydon-setup-grid", "--inifile", "in.ini", "--grid-type", "fixed"],
        )
        monkeypatch.setattr(
            totest.subprocess,
            "Popen",
            self._which_popen((b"/usr/bin/posydon-run-grid\n", b"")),
        )
        monkeypatch.setattr(
            totest,
            "parse_config_file",
            lambda inifile: (
                {"grid": str(grid)},
                {"job_array": True, "number_of_cpus_per_task": 1, "work_dir": ""},
                {
                    "scenario": None,
                    "star_history_columns": str(star_hist),
                    "binary_history_columns": str(binary_hist),
                    "profile_columns": str(profile_hist),
                },
                {"mesa_binary_extras": str(extras)},
            ),
        )
        self._mock_executables_and_inlists(monkeypatch)
        system_calls = []
        monkeypatch.setattr(totest.os, "system", system_calls.append)

        totest.main()

        grid_command = (tmp_path / "grid_command.sh").read_text()
        assert "--grid-point-index $SLURM_ARRAY_TASK_ID" in grid_command
        assert "for SLURM_ARRAY_TASK_ID in 0 1 ; do" in grid_command
        assert " ; done" in grid_command
        assert system_calls == ["chmod 755 grid_command.sh"]

    def test_main_fixed_slurm_job_array_scripts(
        self, monkeypatch, tmp_path, grid_inputs
    ):
        """Slurm job array submission writes submit, cleanup and run scripts."""
        grid, star_hist, binary_hist, profile_hist, extras = grid_inputs
        os.chdir(tmp_path)
        monkeypatch.setenv("MESASDK_ROOT", "/opt/mesasdk")
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "posydon-setup-grid",
                "--inifile",
                "in.ini",
                "--grid-type",
                "fixed",
                "--submission-type",
                "slurm",
            ],
        )
        monkeypatch.setattr(
            totest.subprocess,
            "Popen",
            self._which_popen((b"/usr/bin/posydon-run-grid\n", b"")),
        )
        monkeypatch.setattr(
            totest,
            "parse_config_file",
            lambda inifile: (
                {"grid": str(grid)},
                {
                    "job_array": True,
                    "account": "myaccount",
                    "partition": "main",
                    "number_of_cpus_per_task": 2,
                    "walltime": "01:00:00",
                    "email": "a@b.c",
                },
                {
                    "scenario": None,
                    "star_history_columns": str(star_hist),
                    "binary_history_columns": str(binary_hist),
                    "profile_columns": str(profile_hist),
                },
                {"mesa_binary_extras": str(extras)},
            ),
        )
        self._mock_executables_and_inlists(monkeypatch)
        system_calls = []
        monkeypatch.setattr(totest.os, "system", system_calls.append)

        totest.main()

        assert (tmp_path / "job_array_grid_submit.slurm").is_file()
        assert (tmp_path / "cleanup.slurm").is_file()
        assert (tmp_path / "run_grid.sh").is_file()
        assert not (tmp_path / "mpi_grid_submit.slurm").exists()

        submit = (tmp_path / "job_array_grid_submit.slurm").read_text()
        assert "#SBATCH --array=0-1" in submit
        assert "#SBATCH --account=myaccount" in submit
        assert "#SBATCH --partition=main" in submit
        assert "#SBATCH --cpus-per-task 2" in submit
        assert '#SBATCH --job-name="mesa_grid_\\${SLURM_ARRAY_TASK_ID}"' in submit
        assert "export OMP_NUM_THREADS=2" in submit
        assert "--grid-point-index $SLURM_ARRAY_TASK_ID" in submit

        run_grid = (tmp_path / "run_grid.sh").read_text()
        assert "ID_GRID=$(sbatch --parsable job_array_grid_submit.slurm)" in run_grid
        assert "compress-mesa ." in (tmp_path / "cleanup.slurm").read_text()
        assert system_calls == ["chmod 755 run_grid.sh"]

    def test_main_fixed_slurm_mpi_scripts(self, monkeypatch, tmp_path, grid_inputs):
        """Slurm MPI submission writes MPI submit, cleanup and run scripts."""
        grid, star_hist, binary_hist, profile_hist, extras = grid_inputs
        os.chdir(tmp_path)
        monkeypatch.setenv("MESASDK_ROOT", "/opt/mesasdk")
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "posydon-setup-grid",
                "--inifile",
                "in.ini",
                "--grid-type",
                "fixed",
                "--submission-type",
                "slurm",
            ],
        )
        monkeypatch.setattr(
            totest.subprocess,
            "Popen",
            self._which_popen((b"/usr/bin/posydon-run-grid\n", b"")),
        )
        monkeypatch.setattr(
            totest,
            "parse_config_file",
            lambda inifile: (
                {"grid": str(grid), "psycris_inifile": "/tmp/psycris.ini"},
                {
                    "job_array": False,
                    "number_of_mpi_tasks": 2,
                    "number_of_nodes": 2,
                    "number_of_cpus_per_task": 4,
                    "account": "acct",
                    "partition": "part",
                    "walltime": "02:00:00",
                    "email": "x@y.z",
                },
                {
                    "scenario": None,
                    "star_history_columns": str(star_hist),
                    "binary_history_columns": str(binary_hist),
                    "profile_columns": str(profile_hist),
                },
                {"mesa_binary_extras": str(extras)},
            ),
        )
        self._mock_executables_and_inlists(monkeypatch)
        system_calls = []
        monkeypatch.setattr(totest.os, "system", system_calls.append)

        totest.main()

        assert (tmp_path / "mpi_grid_submit.slurm").is_file()
        assert (tmp_path / "cleanup.slurm").is_file()
        assert (tmp_path / "run_grid.sh").is_file()
        assert not (tmp_path / "job_array_grid_submit.slurm").exists()

        mpi_submit = (tmp_path / "mpi_grid_submit.slurm").read_text()
        assert "#SBATCH -N 2" in mpi_submit
        assert "#SBATCH --ntasks-per-node 2" in mpi_submit
        assert "#SBATCH --cpus-per-task 4" in mpi_submit
        assert "--job_end $SLURM_JOB_END_TIME" in mpi_submit
        assert "--output-directory" in mpi_submit
        assert system_calls == ["chmod 755 run_grid.sh"]

    def test_main_psycris_missing_dynamic(self, monkeypatch, tmp_path):
        """Dynamic grid without a psycris inifile raises ValueError."""
        grid = tmp_path / "grid.csv"
        grid.write_text("m1,m2\n10.0,8.0\n")
        os.chdir(tmp_path)
        monkeypatch.setattr(
            sys,
            "argv",
            ["posydon-setup-grid", "--inifile", "in.ini", "--grid-type", "dynamic"],
        )
        monkeypatch.setattr(
            totest.subprocess,
            "Popen",
            self._which_popen((b"/usr/bin/posydon-run-grid\n", b"")),
        )
        monkeypatch.setattr(
            totest,
            "parse_config_file",
            lambda inifile: (
                {"grid": str(grid)},
                {
                    "job_array": False,
                    "number_of_mpi_tasks": 1,
                    "number_of_nodes": 1,
                    "number_of_cpus_per_task": 1,
                },
                {"scenario": None},
                {},
            ),
        )
        with pytest.raises(ValueError, match="Please add psycris inifile"):
            totest.main()

    def test_main_scenario_calls_setup_inlists_from_scenario(
        self, monkeypatch, tmp_path, grid_inputs
    ):
        """A scenario in mesa_inlists triggers setup_inlists_from_scenario."""
        grid, star_hist, binary_hist, profile_hist, extras = grid_inputs
        os.chdir(tmp_path)
        monkeypatch.setenv("MESASDK_ROOT", "/opt/mesasdk")
        monkeypatch.setattr(
            sys,
            "argv",
            ["posydon-setup-grid", "--inifile", "in.ini", "--grid-type", "fixed"],
        )
        monkeypatch.setattr(
            totest.subprocess,
            "Popen",
            self._which_popen((b"/usr/bin/posydon-run-grid\n", b"")),
        )
        monkeypatch.setattr(
            totest,
            "parse_config_file",
            lambda inifile: (
                {"grid": str(grid), "psycris_inifile": "/tmp/psycris.ini"},
                {
                    "job_array": False,
                    "number_of_mpi_tasks": 1,
                    "number_of_nodes": 1,
                    "number_of_cpus_per_task": 1,
                    "work_dir": "",
                },
                {
                    "scenario": ["posydon", "master-abc123", "HeMS-HMS"],
                    "star_history_columns": str(star_hist),
                    "binary_history_columns": str(binary_hist),
                    "profile_columns": str(profile_hist),
                },
                {"mesa_binary_extras": str(extras)},
            ),
        )
        self._mock_executables_and_inlists(monkeypatch)
        scenario_calls = []
        monkeypatch.setattr(
            totest,
            "setup_inlists_from_scenario",
            lambda source, gitcommit, system_type, mesa_inlists, mesa_extras: scenario_calls.append(
                (source, gitcommit, system_type, mesa_inlists, mesa_extras)
            ),
        )
        system_calls = []
        monkeypatch.setattr(totest.os, "system", system_calls.append)

        totest.main()

        assert len(scenario_calls) == 1
        source, gitcommit, system_type, mi, me = scenario_calls[0]
        assert source == "posydon"
        assert gitcommit == "master-abc123"
        assert system_type == "HeMS-HMS"
        assert mi["scenario"] == ["posydon", "master-abc123", "HeMS-HMS"]
        assert me["mesa_binary_extras"] == str(extras)

    def test_main_dynamic_uses_psycris_inifile_columns(
        self, monkeypatch, tmp_path, grid_inputs
    ):
        """Dynamic grids source grid params from the psycris inifile columns."""
        grid, star_hist, binary_hist, profile_hist, extras = grid_inputs
        os.chdir(tmp_path)
        monkeypatch.setenv("MESASDK_ROOT", "/opt/mesasdk")
        monkeypatch.setattr(
            sys,
            "argv",
            ["posydon-setup-grid", "--inifile", "in.ini", "--grid-type", "dynamic"],
        )
        monkeypatch.setattr(
            totest.subprocess,
            "Popen",
            self._which_popen((b"/usr/bin/posydon-run-grid\n", b"")),
        )
        monkeypatch.setattr(
            totest,
            "parse_config_file",
            lambda inifile: (
                {"grid": str(grid), "psycris_inifile": "/tmp/psycris.ini"},
                {
                    "job_array": False,
                    "number_of_mpi_tasks": 1,
                    "number_of_nodes": 1,
                    "number_of_cpus_per_task": 1,
                    "work_dir": "",
                },
                {
                    "scenario": None,
                    "star_history_columns": str(star_hist),
                    "binary_history_columns": str(binary_hist),
                    "profile_columns": str(profile_hist),
                },
                {"mesa_binary_extras": str(extras)},
            ),
        )
        calls = {}

        def fake_construct_static_inlist(mesa_inlists, grid_parameters, working_directory):
            calls["grid_parameters"] = list(grid_parameters)
            return (
                "/work/star1/inlist_step0",
                None,
                "/work/binary/inlist_project",
                "/work/binary/inlist1",
                "/work/binary/inlist2",
            )

        monkeypatch.setattr(totest, "construct_static_inlist", fake_construct_static_inlist)
        monkeypatch.setattr(
            totest,
            "make_executables",
            lambda mesa_extras, working_directory: (
                "/work/binary/binary",
                "/work/star1/star",
                "/work/star2/star",
            ),
        )
        monkeypatch.setattr(
            totest,
            "parse_inifile",
            lambda path: {
                "posydon_dynamic_sampling_kwargs": {"mesa_column_names": ["m1", "m2"]}
            },
        )
        system_calls = []
        monkeypatch.setattr(totest.os, "system", system_calls.append)

        totest.main()

        assert calls["grid_parameters"] == ["m1", "m2"]
        assert system_calls == ["chmod 755 grid_command.sh"]

    def test_main_replaces_existing_column_lists_folder(
        self, monkeypatch, tmp_path, grid_inputs
    ):
        """An existing column_lists folder is removed and recreated."""
        grid, star_hist, binary_hist, profile_hist, extras = grid_inputs
        os.chdir(tmp_path)
        column_lists = tmp_path / "column_lists"
        column_lists.mkdir()
        (column_lists / "old.txt").write_text("old")
        monkeypatch.setenv("MESASDK_ROOT", "/opt/mesasdk")
        monkeypatch.setattr(
            sys,
            "argv",
            ["posydon-setup-grid", "--inifile", "in.ini", "--grid-type", "fixed"],
        )
        monkeypatch.setattr(
            totest.subprocess,
            "Popen",
            self._which_popen((b"/usr/bin/posydon-run-grid\n", b"")),
        )
        monkeypatch.setattr(
            totest,
            "parse_config_file",
            lambda inifile: (
                {"grid": str(grid), "psycris_inifile": "/tmp/psycris.ini"},
                {
                    "job_array": False,
                    "number_of_mpi_tasks": 1,
                    "number_of_nodes": 1,
                    "number_of_cpus_per_task": 1,
                    "work_dir": "",
                },
                {
                    "scenario": None,
                    "star_history_columns": str(star_hist),
                    "binary_history_columns": str(binary_hist),
                    "profile_columns": str(profile_hist),
                },
                {"mesa_binary_extras": str(extras)},
            ),
        )
        self._mock_executables_and_inlists(monkeypatch)
        system_calls = []
        monkeypatch.setattr(totest.os, "system", system_calls.append)

        totest.main()

        assert not (column_lists / "old.txt").exists()
        assert (
            tmp_path / "column_lists" / "history_columns.list"
        ).read_text() == "star columns\n"

    def test_main_shell_appends_temporary_directory(
        self, monkeypatch, tmp_path, grid_inputs
    ):
        """A non-empty slurm work_dir appends --temporary-directory."""
        grid, star_hist, binary_hist, profile_hist, extras = grid_inputs
        os.chdir(tmp_path)
        monkeypatch.setenv("MESASDK_ROOT", "/opt/mesasdk")
        monkeypatch.setattr(
            sys,
            "argv",
            ["posydon-setup-grid", "--inifile", "in.ini", "--grid-type", "fixed"],
        )
        monkeypatch.setattr(
            totest.subprocess,
            "Popen",
            self._which_popen((b"/usr/bin/posydon-run-grid\n", b"")),
        )
        monkeypatch.setattr(
            totest,
            "parse_config_file",
            lambda inifile: (
                {"grid": str(grid), "psycris_inifile": "/tmp/psycris.ini"},
                {
                    "job_array": False,
                    "number_of_mpi_tasks": 1,
                    "number_of_nodes": 1,
                    "number_of_cpus_per_task": 1,
                    "work_dir": "/scratch/x",
                },
                {
                    "scenario": None,
                    "star_history_columns": str(star_hist),
                    "binary_history_columns": str(binary_hist),
                    "profile_columns": str(profile_hist),
                },
                {"mesa_binary_extras": str(extras)},
            ),
        )
        self._mock_executables_and_inlists(monkeypatch)
        system_calls = []
        monkeypatch.setattr(totest.os, "system", system_calls.append)

        totest.main()

        grid_command = (tmp_path / "grid_command.sh").read_text()
        assert "--temporary-directory /scratch/x" in grid_command
        assert system_calls == ["chmod 755 grid_command.sh"]

    def test_run_setup_unknown_submission_type_writes_nothing(
        self, monkeypatch, tmp_path, grid_inputs
    ):
        """An unrecognized submission type writes no scripts."""
        grid, star_hist, binary_hist, profile_hist, extras = grid_inputs
        os.chdir(tmp_path)
        args = SimpleNamespace(
            inifile="in.ini",
            grid_type="fixed",
            run_directory=str(tmp_path),
            submission_type="other",
            verbose=False,
        )
        monkeypatch.setattr(
            totest,
            "find_run_grid_executable",
            lambda: "/usr/bin/posydon-run-grid",
        )
        monkeypatch.setattr(
            totest,
            "parse_config_file",
            lambda inifile: (
                {"grid": str(grid), "psycris_inifile": "/tmp/psycris.ini"},
                {
                    "job_array": False,
                    "number_of_mpi_tasks": 1,
                    "number_of_nodes": 1,
                    "number_of_cpus_per_task": 1,
                    "work_dir": "",
                },
                {
                    "scenario": None,
                    "star_history_columns": str(star_hist),
                    "binary_history_columns": str(binary_hist),
                    "profile_columns": str(profile_hist),
                },
                {"mesa_binary_extras": str(extras)},
            ),
        )
        self._mock_executables_and_inlists(monkeypatch)
        system_calls = []
        monkeypatch.setattr(totest.os, "system", system_calls.append)

        totest.run_setup(args)

        assert not (tmp_path / "grid_command.sh").exists()
        assert not (tmp_path / "job_array_grid_submit.slurm").exists()
        assert system_calls == []
