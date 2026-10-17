"""Unit tests of posydon/grids/setup/submit.py
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
import posydon.grids.setup.submit as totest
from posydon.grids.setup.build import BuiltInlists
from posydon.grids.setup.config import GridConfig, SubmissionSection

EXECUTABLES = ("/run/binary/binary", "/run/star1/star", "/run/star2/star")
COLUMN_LISTS = {"history_columns": "/run/column_lists/h.list",
                "binary_history_columns": "/run/column_lists/b.list",
                "profile_columns": "/run/column_lists/p.list"}
INLISTS = BuiltInlists("/run/binary/inlist_project", "/run/binary/inlist1",
                       "/run/binary/inlist2", None, None)
OPTIONS = ("--mesa-binary-executable /run/binary/binary "
           "--mesa-star1-executable /run/star1/star "
           "--mesa-star2-executable /run/star2/star "
           "--mesa-binary-inlist-project /run/binary/inlist_project "
           "--mesa-binary-inlist1 /run/binary/inlist1 "
           "--mesa-binary-inlist2 /run/binary/inlist2 "
           "--mesa-star1-inlist-project None "
           "--mesa-star2-inlist-project None "
           "--mesa-star-history-columns /run/column_lists/h.list "
           "--mesa-binary-history-columns /run/column_lists/b.list "
           "--mesa-profile-columns /run/column_lists/p.list "
           "--output-directory /run ")


def make_config(backend="slurm_array", **changes):
    config = GridConfig()
    config.grid.points = Path("/run/grid.csv")
    config.submission = SubmissionSection(
        backend=backend, account="acc", partition="part", walltime="1:00",
        email="a@b.c", cpus_per_task=4, nodes=2, mpi_tasks=8)
    for key, value in changes.items():
        setattr(config.submission, key, value)
    return config


def command(config):
    return totest.command_line(config, INLISTS, EXECUTABLES, COLUMN_LISTS,
                               "/run", "/bin/run-grid")


class TestElements:
    def test_dir(self):
        elements = {'CLEANUP_SCRIPT', 'GRID_SCRIPTS', 'Path', 'Pwarn',
                    'RUN_GRID_OPTIONS', 'RUN_SCRIPT', '__authors__',
                    '__builtins__', '__cached__', '__doc__', '__file__',
                    '__loader__', '__name__', '__package__', '__spec__',
                    '_environment', '_finish', '_sbatch', '_time', '_write',
                    'cleanup_script', 'command_line', 'find_run_grid',
                    'one_point_per_task', 'run_grid_call', 'run_script',
                    'shell_script', 'shutil', 'slurm_grid_script',
                    'write_scripts'}
        assert set(dir(totest)) == elements

    @pytest.mark.parametrize("name", [
        'find_run_grid', 'one_point_per_task', 'run_grid_call',
        'command_line', 'shell_script', 'slurm_grid_script',
        'cleanup_script', 'run_script', 'write_scripts'])
    def test_instance_functions(self, name):
        assert isroutine(getattr(totest, name))

    def test_doctests(self):
        result = doctest.testmod(totest)
        assert result.failed == 0
        assert result.attempted > 0


class TestFindRunGrid:
    def test_found(self, monkeypatch):
        monkeypatch.setattr(totest.shutil, "which", lambda name: "/bin/x")
        assert totest.find_run_grid() == "/bin/x"

    def test_missing(self, monkeypatch):
        monkeypatch.setattr(totest.shutil, "which", lambda name: None)
        with pytest.raises(ValueError, match="Cannot locate"):
            totest.find_run_grid()


class TestOnePointPerTask:
    @pytest.mark.parametrize("backend, job_array, expected", [
        ("slurm_array", False, True), ("slurm_mpi", True, False),
        ("shell", True, True), ("shell", False, False)])
    def test_backends(self, backend, job_array, expected):
        submission = SubmissionSection(backend=backend, job_array=job_array)
        assert totest.one_point_per_task(submission) is expected


class TestCommandLine:
    def test_job_array(self):
        assert command(make_config()) == (
            "python /bin/run-grid --mesa-grid /run/grid.csv " + OPTIONS
            + "--grid-type fixed --psycris-inifile None "
            "--grid-point-index $SLURM_ARRAY_TASK_ID "
            "--job_end $SLURM_JOB_END_TIME")

    @pytest.mark.parametrize("backend", ["slurm_array", "shell"])
    def test_dynamic_job_array_is_refused(self, backend):
        """S5: before, it silently ran as a fixed grid."""
        config = make_config(backend)
        config.grid.mode = "dynamic"
        config.grid.psycris_inifile = Path("/run/psycris.ini")
        with pytest.raises(ValueError, match="dynamic grid runs as one MPI "
                           "job"):
            command(config)

    def test_mpi_dynamic(self):
        config = make_config("slurm_mpi", work_dir="/scratch")
        config.grid.mode = "dynamic"
        config.grid.psycris_inifile = Path("/run/psycris.ini")
        config.grid.keep_profiles = True
        config.grid.keep_photos = True
        assert command(config) == (
            "mpirun --bind-to none -np 16 python -m mpi4py /bin/run-grid "
            "--mesa-grid /run/grid.csv " + OPTIONS
            + "--grid-type dynamic --psycris-inifile /run/psycris.ini "
            "--keep_profiles --keep_photos --job_end $SLURM_JOB_END_TIME "
            "--temporary-directory /scratch")

    def test_shell_without_job_array(self):
        result = command(make_config("shell", job_array=False))
        assert result.startswith("python /bin/run-grid ")
        assert result.endswith("--grid-type fixed --psycris-inifile None")

    def test_unknown_mode(self):
        with pytest.raises(ValueError, match="fixed or dynamic"):
            totest.run_grid_call({}, "/bin/run-grid", mode="adaptive")


class TestScripts:
    def test_shell_job_array(self):
        config = make_config("shell", newgroup="posydon")
        assert totest.shell_script(config, "CMD", 3, "/sdk", "/mesa") == (
            "#!/bin/bash\n\nexport OMP_NUM_THREADS=4\n\n"
            "export MESASDK_ROOT=/sdk\n"
            "source $MESASDK_ROOT/bin/mesasdk_init.sh\n"
            "export MESA_DIR=/mesa\n\n"
            "for SLURM_ARRAY_TASK_ID in 0 1 2 ; do CMD ; done\n\n"
            "compress-mesa .\n"
            'echo "Change group to posydon"\nchgrp -fR posydon .\n'
            'echo "Change group permission to rwX at least"\n'
            'chmod -fR g+rwX .\n\necho "Done."')

    def test_shell_without_job_array(self):
        config = make_config("shell", job_array=False)
        assert totest.shell_script(config, "CMD", 3, "/sdk", "/mesa")\
            .endswith('export MESA_DIR=/mesa\n\nCMD\n\ncompress-mesa .\n\n'
                      'echo "Done."')

    def test_slurm_array(self):
        assert totest.slurm_grid_script(make_config(), "CMD", 5, "/sdk",
                                         "/mesa") == (
            "#!/bin/bash\n#SBATCH --account=acc\n#SBATCH --partition=part\n"
            "#SBATCH -N 1\n#SBATCH --array=0-4\n#SBATCH --cpus-per-task 4\n"
            "#SBATCH --ntasks-per-node 1\n#SBATCH --time=1:00\n"
            '#SBATCH --job-name="mesa_grid"\n'
            "#SBATCH --output=mesa_grid.%A_%a.out\n#SBATCH --mail-type=ALL\n"
            "#SBATCH --mail-user=a@b.c\n#SBATCH --mem-per-cpu=4G\n\n"
            "export OMP_NUM_THREADS=4\n\nexport MESASDK_ROOT=/sdk\n"
            "source $MESASDK_ROOT/bin/mesasdk_init.sh\n"
            "export MESA_DIR=/mesa\n\n\nCMD")

    def test_slurm_mpi(self):
        text = totest.slurm_grid_script(make_config("slurm_mpi"), "CMD", 5,
                                        "/sdk", "/mesa")
        assert ("#SBATCH -N 2\n#SBATCH --cpus-per-task 4\n"
                "#SBATCH --ntasks-per-node 8\n#SBATCH --time=1:00\n"
                '#SBATCH --output="mesa_grid.out"\n') in text
        assert "--array" not in text

    def test_cleanup(self):
        assert totest.cleanup_script(make_config()) == (
            "#!/bin/bash\n#SBATCH --account=acc\n#SBATCH --partition=part\n"
            "#SBATCH -N 1\n#SBATCH --cpus-per-task 1\n"
            "#SBATCH --ntasks-per-node 1\n#SBATCH --time=1:00\n"
            '#SBATCH --job-name="mesa_grid_cleanup"\n'
            "#SBATCH --output=mesa_cleanup.out\n#SBATCH --mail-type=ALL\n"
            "#SBATCH --mail-user=a@b.c\n#SBATCH --mem-per-cpu=4G\n\n"
            'compress-mesa .\n\necho "Done."')

    def test_unset_options_are_left_out(self):
        """Before, they were written as --account=None and so on."""
        config = make_config(account=None, partition=None, walltime=None,
                             email=None)
        text = totest.slurm_grid_script(config, "CMD", 5, "/sdk", "/mesa")
        assert "None" not in text
        assert "--mail-type" not in text
        assert "--mem-per-cpu=4G" in text
        assert "None" not in totest.cleanup_script(config)
        config.submission.mem_per_cpu = None
        assert "--mem-per-cpu" not in totest.cleanup_script(config)

    def test_run_script(self):
        assert totest.run_script("grid.slurm") == (
            "#!/bin/bash\nID_GRID=$(sbatch --parsable grid.slurm)\n"
            'echo "grid.slurm submitted as "${ID_GRID}\n'
            "ID_cleanup=$(sbatch --parsable --dependency=afterany:${ID_GRID}"
            " --kill-on-invalid-dep=yes cleanup.slurm)\n"
            'echo "cleanup.slurm submitted as "${ID_cleanup}\n')


class TestWriteScripts:
    def test_shell(self, tmp_path):
        written = totest.write_scripts(make_config("shell"), "CMD", 2,
                                       tmp_path, "/sdk", "/mesa")
        assert written == [tmp_path / "grid_command.sh"]
        assert written[0].stat().st_mode & 0o777 == 0o755

    def test_slurm_and_overwrite_warning(self, tmp_path):
        (tmp_path / "cleanup.slurm").write_text("old")
        with pytest.warns(match="Replace cleanup.slurm"):
            written = totest.write_scripts(make_config("slurm_mpi"), "CMD",
                                           2, str(tmp_path), "/sdk", "/mesa")
        assert [path.name for path in written] == [
            "mpi_grid_submit.slurm", "cleanup.slurm", "run_grid.sh"]
        assert written[1].read_text().startswith("#!/bin/bash")
        assert written[2].stat().st_mode & 0o777 == 0o755
