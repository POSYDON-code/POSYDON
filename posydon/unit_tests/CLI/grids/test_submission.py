"""Tests for posydon.CLI.grids.submission."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import os

import pandas
import pytest

import posydon.CLI.grids.submission as totest
from posydon.utils.posydonwarning import OverwriteWarning


@pytest.fixture(autouse=True)
def _env(tmp_path, monkeypatch):
    """Always provide HOME/MESA_DIR/MESASDK_ROOT and restore cwd after tests."""
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("MESA_DIR", str(tmp_path / "mesa"))
    monkeypatch.setenv("MESASDK_ROOT", "/opt/mesasdk")
    original = os.getcwd()
    yield
    os.chdir(original)


class TestConstructCommandLine:
    """Tests for construct_command_line()."""

    def test_run_shell(self, monkeypatch):
        """_run_shell forwards the command to os.system."""
        calls = []
        monkeypatch.setattr(totest.os, "system", calls.append)
        totest._run_shell("chmod 755 grid_command.sh")
        assert calls == ["chmod 755 grid_command.sh"]

    def test_fixed_grid(self):
        """Fixed grid command line uses the 'python {15} ...' template."""
        command = totest.construct_command_line(
            1,
            "/grid/points.csv",
            "/work/binary/binary",
            "/work/star1/star",
            "/work/star2/star",
            "/work/binary/inlist_project",
            "/work/binary/inlist1",
            "/work/binary/inlist2",
            "/work/star1/inlist_step0",
            None,
            "/work/column_lists/history_columns.list",
            "/work/column_lists/binary_history_columns.list",
            "/work/column_lists/profile_columns.list",
            "/work",
            "fixed",
            "/path/to/posydon-run-grid",
            None,
        )
        expected = (
            "python /path/to/posydon-run-grid --mesa-grid /grid/points.csv "
            "--mesa-binary-executable /work/binary/binary "
            "--mesa-star1-executable /work/star1/star "
            "--mesa-star2-executable /work/star2/star "
            "--mesa-binary-inlist-project /work/binary/inlist_project "
            "--mesa-binary-inlist1 /work/binary/inlist1 "
            "--mesa-binary-inlist2 /work/binary/inlist2 "
            "--mesa-star1-inlist-project /work/star1/inlist_step0 "
            "--mesa-star2-inlist-project None "
            "--mesa-star-history-columns /work/column_lists/history_columns.list "
            "--mesa-binary-history-columns /work/column_lists/binary_history_columns.list "
            "--mesa-profile-columns /work/column_lists/profile_columns.list "
            "--output-directory /work --grid-type fixed --psycris-inifile None"
        )
        assert command == expected

    def test_dynamic_grid(self):
        """Dynamic grid command line wraps the run-grid call in mpirun/mpi4py."""
        command = totest.construct_command_line(
            4,
            "/grid/points.h5",
            "/work/binary/binary",
            "/work/star1/star",
            "/work/star2/star",
            "/work/binary/inlist_project",
            "/work/binary/inlist1",
            "/work/binary/inlist2",
            "/work/star1/inlist_step0",
            None,
            "/work/history.list",
            "/work/binary_history.list",
            "/work/profile.list",
            "/work",
            "dynamic",
            "/path/to/posydon-run-grid",
            "/work/psycris.ini",
        )
        assert command.startswith(
            "mpirun --bind-to none -np 4 python -m mpi4py /path/to/posydon-run-grid "
            "--mesa-grid /grid/points.h5 --mesa-binary-executable /work/binary/binary "
        )
        assert " --grid-type dynamic " in command
        assert " --psycris-inifile /work/psycris.ini" in command

    def test_keep_profiles_and_photos(self):
        """keep_profiles/keep_photos append the corresponding flags."""
        command = totest.construct_command_line(
            1,
            "/grid.csv",
            "/b",
            "/s1",
            "/s2",
            "/p",
            "/i1",
            "/i2",
            None,
            None,
            "/h",
            "/bh",
            "/pr",
            "/out",
            "fixed",
            "/path/to/posydon-run-grid",
            None,
            keep_profiles=True,
            keep_photos=True,
        )
        assert command.endswith(" --keep_profiles --keep_photos")

    def test_invalid_grid_type(self):
        """An unrecognized grid_type raises ValueError."""
        with pytest.raises(
            ValueError, match="grid_type can either be fixed or dynamic not anything else"
        ):
            totest.construct_command_line(
                1,
                "/grid.csv",
                "/b",
                "/s1",
                "/s2",
                "/p",
                "/i1",
                "/i2",
                None,
                None,
                "/h",
                "/bh",
                "/pr",
                "/out",
                "bogus",
                "/path/to/posydon-run-grid",
                None,
            )

    def test_keep_profiles_only(self):
        """keep_profiles alone appends only the --keep_profiles flag."""
        command = totest.construct_command_line(
            1,
            "/grid.csv",
            "/b",
            "/s1",
            "/s2",
            "/p",
            "/i1",
            "/i2",
            None,
            None,
            "/h",
            "/bh",
            "/pr",
            "/out",
            "fixed",
            "/path/to/posydon-run-grid",
            None,
            keep_profiles=True,
        )
        assert command.endswith(" --keep_profiles")
        assert " --keep_photos" not in command

    def test_keep_photos_only(self):
        """keep_photos alone appends only the --keep_photos flag."""
        command = totest.construct_command_line(
            1,
            "/grid.csv",
            "/b",
            "/s1",
            "/s2",
            "/p",
            "/i1",
            "/i2",
            None,
            None,
            "/h",
            "/bh",
            "/pr",
            "/out",
            "fixed",
            "/path/to/posydon-run-grid",
            None,
            keep_photos=True,
        )
        assert command.endswith(" --keep_photos")
        assert " --keep_profiles" not in command


class TestWriteShellSubmissionScript:
    """Tests for write_shell_submission_script()."""

    @staticmethod
    def _grid_df():
        """Two-row grid DataFrame used by the submission writers."""
        return pandas.DataFrame({"m1": [10.0, 11.0], "m2": [8.0, 9.0]})

    def test_shell_happy_path(self, tmp_path, monkeypatch):
        """Fixed grid shell submission writes a runnable grid_command.sh."""
        os.chdir(tmp_path)
        slurm = {
            "job_array": False,
            "number_of_mpi_tasks": 1,
            "number_of_nodes": 1,
            "number_of_cpus_per_task": 1,
            "work_dir": "",
        }
        command_line = (
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
        system_calls = []
        monkeypatch.setattr(totest, "_run_shell", system_calls.append)

        totest.write_shell_submission_script(
            command_line, slurm, self._grid_df(), run_directory=str(tmp_path)
        )

        grid_command = (tmp_path / "grid_command.sh").read_text()
        assert "#!/bin/bash" in grid_command
        assert "export OMP_NUM_THREADS=1" in grid_command
        assert "export MESASDK_ROOT=/opt/mesasdk" in grid_command
        assert "export MESA_DIR={0}".format(os.environ["MESA_DIR"]) in grid_command
        assert "compress-mesa ." in grid_command
        assert 'echo "Done."' in grid_command
        assert command_line in grid_command
        assert system_calls == ["chmod 755 grid_command.sh"]

    def test_shell_job_array(self, tmp_path, monkeypatch):
        """Fixed grid shell submission with a job array loops over grid rows."""
        os.chdir(tmp_path)
        slurm = {"job_array": True, "number_of_cpus_per_task": 1, "work_dir": ""}
        command_line = "python /usr/bin/posydon-run-grid --grid-point-index $SLURM_ARRAY_TASK_ID"
        system_calls = []
        monkeypatch.setattr(totest, "_run_shell", system_calls.append)

        totest.write_shell_submission_script(
            command_line, slurm, self._grid_df(), run_directory=str(tmp_path)
        )

        grid_command = (tmp_path / "grid_command.sh").read_text()
        assert "--grid-point-index $SLURM_ARRAY_TASK_ID" in grid_command
        assert "for SLURM_ARRAY_TASK_ID in 0 1 ; do" in grid_command
        assert " ; done" in grid_command
        assert system_calls == ["chmod 755 grid_command.sh"]

    def test_shell_job_array_empty_grid(self, tmp_path, monkeypatch):
        """A job-array loop over an empty grid writes no task indices."""
        os.chdir(tmp_path)
        slurm = {"job_array": True, "number_of_cpus_per_task": 1, "work_dir": ""}
        command_line = "python /usr/bin/posydon-run-grid"
        system_calls = []
        monkeypatch.setattr(totest, "_run_shell", system_calls.append)
        empty_grid = pandas.DataFrame({"m1": [], "m2": []})

        totest.write_shell_submission_script(
            command_line, slurm, empty_grid, run_directory=str(tmp_path)
        )

        grid_command = (tmp_path / "grid_command.sh").read_text()
        assert "for SLURM_ARRAY_TASK_ID in ; do " in grid_command
        assert " ; done" in grid_command
        assert system_calls == ["chmod 755 grid_command.sh"]

    def test_shell_newgroup(self, tmp_path, monkeypatch):
        """A newgroup key adds the group change lines to grid_command.sh."""
        os.chdir(tmp_path)
        slurm = {
            "job_array": False,
            "number_of_cpus_per_task": 1,
            "newgroup": "mygroup",
        }
        command_line = "python /usr/bin/posydon-run-grid"
        system_calls = []
        monkeypatch.setattr(totest, "_run_shell", system_calls.append)

        totest.write_shell_submission_script(
            command_line, slurm, self._grid_df(), run_directory=str(tmp_path)
        )

        grid_command = (tmp_path / "grid_command.sh").read_text()
        assert 'echo "Change group to mygroup"' in grid_command
        assert "chgrp -fR mygroup ." in grid_command
        assert 'echo "Change group permission to rwX at least"' in grid_command
        assert "chmod -fR g+rwX ." in grid_command
        assert system_calls == ["chmod 755 grid_command.sh"]

    def test_shell_grid_command_exists(self, tmp_path, monkeypatch):
        """An existing grid_command.sh triggers an OverwriteWarning."""
        os.chdir(tmp_path)
        (tmp_path / "grid_command.sh").write_text("old")
        slurm = {"job_array": False, "number_of_cpus_per_task": 1}
        command_line = "python /usr/bin/posydon-run-grid"
        system_calls = []
        monkeypatch.setattr(totest, "_run_shell", system_calls.append)

        with pytest.warns(OverwriteWarning, match="Replace grid_command.sh"):
            totest.write_shell_submission_script(
                command_line, slurm, self._grid_df(), run_directory=str(tmp_path)
            )
        assert not (tmp_path / "grid_command.sh").read_text().startswith("old")
        assert system_calls == ["chmod 755 grid_command.sh"]


class TestWriteSlurmSubmissionScripts:
    """Tests for write_slurm_submission_scripts()."""

    @staticmethod
    def _grid_df():
        """Two-row grid DataFrame used by the submission writers."""
        return pandas.DataFrame({"m1": [10.0, 11.0], "m2": [8.0, 9.0]})

    def test_slurm_job_array_scripts(self, tmp_path, monkeypatch):
        """Slurm job array submission writes submit, cleanup and run scripts."""
        os.chdir(tmp_path)
        slurm = {
            "job_array": True,
            "account": "myaccount",
            "partition": "main",
            "number_of_cpus_per_task": 2,
            "walltime": "01:00:00",
            "email": "a@b.c",
        }
        command_line = "python /usr/bin/posydon-run-grid --grid-point-index $SLURM_ARRAY_TASK_ID"
        system_calls = []
        monkeypatch.setattr(totest, "_run_shell", system_calls.append)

        result = totest.write_slurm_submission_scripts(
            command_line, slurm, self._grid_df(), run_directory=str(tmp_path)
        )

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
        assert result == ("job_array_grid_submit.slurm", "run_grid.sh")

    def test_slurm_mpi_scripts(self, tmp_path, monkeypatch):
        """Slurm MPI submission writes MPI submit, cleanup and run scripts."""
        os.chdir(tmp_path)
        slurm = {
            "job_array": False,
            "number_of_mpi_tasks": 2,
            "number_of_nodes": 2,
            "number_of_cpus_per_task": 4,
            "account": "acct",
            "partition": "part",
            "walltime": "02:00:00",
            "email": "x@y.z",
        }
        command_line = (
            "python /usr/bin/posydon-run-grid --job_end $SLURM_JOB_END_TIME "
            "--output-directory /work"
        )
        system_calls = []
        monkeypatch.setattr(totest, "_run_shell", system_calls.append)

        result = totest.write_slurm_submission_scripts(
            command_line, slurm, self._grid_df(), run_directory=str(tmp_path)
        )

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
        assert result == ("mpi_grid_submit.slurm", "run_grid.sh")

    def test_slurm_cleanup_newgroup(self, tmp_path, monkeypatch):
        """A newgroup key adds the group change lines to cleanup.slurm."""
        os.chdir(tmp_path)
        slurm = {
            "job_array": True,
            "account": "myaccount",
            "partition": "main",
            "number_of_cpus_per_task": 1,
            "walltime": "01:00:00",
            "email": "a@b.c",
            "newgroup": "mygroup",
        }
        command_line = "python /usr/bin/posydon-run-grid"
        system_calls = []
        monkeypatch.setattr(totest, "_run_shell", system_calls.append)

        totest.write_slurm_submission_scripts(
            command_line, slurm, self._grid_df(), run_directory=str(tmp_path)
        )

        cleanup = (tmp_path / "cleanup.slurm").read_text()
        assert 'echo "Change group to mygroup"' in cleanup
        assert "chgrp -fR mygroup ." in cleanup
        assert 'echo "Change group permission to rwX at least"' in cleanup
        assert "chmod -fR g+rwX ." in cleanup
        assert system_calls == ["chmod 755 run_grid.sh"]

    def test_slurm_job_array_existing_files(self, tmp_path, monkeypatch):
        """Existing slurm scripts trigger OverwriteWarnings on overwrite."""
        os.chdir(tmp_path)
        for name in ["job_array_grid_submit.slurm", "cleanup.slurm", "run_grid.sh"]:
            (tmp_path / name).write_text("old")
        slurm = {
            "job_array": True,
            "account": "myaccount",
            "partition": "main",
            "number_of_cpus_per_task": 1,
            "walltime": "01:00:00",
            "email": "a@b.c",
        }
        command_line = "python /usr/bin/posydon-run-grid"
        system_calls = []
        monkeypatch.setattr(totest, "_run_shell", system_calls.append)

        with pytest.warns(OverwriteWarning) as record:
            totest.write_slurm_submission_scripts(
                command_line, slurm, self._grid_df(), run_directory=str(tmp_path)
            )
        messages = sorted(str(w.message).strip(chr(39)) for w in record)
        assert messages == [
            "Replace cleanup.slurm",
            "Replace job_array_grid_submit.slurm",
            "Replace run_grid.sh",
        ]
        assert "old" not in (tmp_path / "cleanup.slurm").read_text()
        assert system_calls == ["chmod 755 run_grid.sh"]

    def test_slurm_mpi_existing_files(self, tmp_path, monkeypatch):
        """Existing MPI slurm scripts trigger OverwriteWarnings on overwrite."""
        os.chdir(tmp_path)
        for name in ["mpi_grid_submit.slurm", "cleanup.slurm", "run_grid.sh"]:
            (tmp_path / name).write_text("old")
        slurm = {
            "job_array": False,
            "number_of_mpi_tasks": 2,
            "number_of_nodes": 2,
            "number_of_cpus_per_task": 4,
            "account": "acct",
            "partition": "part",
            "walltime": "02:00:00",
            "email": "x@y.z",
        }
        command_line = "python /usr/bin/posydon-run-grid"
        system_calls = []
        monkeypatch.setattr(totest, "_run_shell", system_calls.append)

        with pytest.warns(OverwriteWarning) as record:
            totest.write_slurm_submission_scripts(
                command_line, slurm, self._grid_df(), run_directory=str(tmp_path)
            )
        messages = sorted(str(w.message).strip(chr(39)) for w in record)
        assert messages == [
            "Replace cleanup.slurm",
            "Replace mpi_grid_submit.slurm",
            "Replace run_grid.sh",
        ]
        assert not (tmp_path / "mpi_grid_submit.slurm").read_text().startswith("old")
        assert system_calls == ["chmod 755 run_grid.sh"]
