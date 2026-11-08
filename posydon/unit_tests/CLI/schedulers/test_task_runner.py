"""Unit tests of posydon/CLI/schedulers/task_runner.py"""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>",
]

import sys

import pytest

# import the module which will be tested
import posydon.CLI.schedulers.task_runner as totest

# a command writing the task environment to stdout
PRINT_ENV = [sys.executable, "-c",
             "import os, sys; print(' '.join(os.environ['POSYDON_ARRAY_' + k]"
             " for k in ['JOB_ID', 'TASK_ID', 'TASK_COUNT', 'TASK_MIN']),"
             " sys.argv[1:])", "{task_id}", "{job_id}"]


@pytest.fixture
def clean_env(monkeypatch):
    for var in totest.RANK_ENV_VARS + ('POSYDON_ARRAY_JOB_ID', 'PBS_JOBID'):
        monkeypatch.delenv(var, raising=False)
    return monkeypatch


class TestHelpers:

    def test_get_rank(self, clean_env):
        assert totest.get_rank({}) is None
        assert totest.get_rank({'PMI_RANK': '3'}) == 3
        assert totest.get_rank({'PALS_RANKID': '1', 'PMI_RANK': '3'}) == 1
        assert totest.get_rank({'PALS_RANKID': 'x', 'PMI_RANK': '3'}) == 3
        assert totest.get_rank() is None
        clean_env.setenv('OMPI_COMM_WORLD_RANK', '5')
        assert totest.get_rank() == 5

    def test_get_job_id(self, clean_env):
        assert totest.get_job_id({}) == "0"
        assert totest.get_job_id({'PBS_JOBID': '1234.crux-pbs-01'}) == "1234"
        assert totest.get_job_id({'PBS_JOBID': '99[2].server'}) == "99"
        assert totest.get_job_id({'POSYDON_ARRAY_JOB_ID': '7',
                                  'PBS_JOBID': '1234.server'}) == "7"
        assert totest.get_job_id() == "0"

    def test_get_my_tasks(self):
        assert totest.get_my_tasks(0, 3, ntasks=8) == [0, 3, 6]
        assert totest.get_my_tasks(2, 3, ntasks=8) == [2, 5]
        assert totest.get_my_tasks(1, 2, task_ids=[1, 4, 9]) == [4]
        assert totest.get_my_tasks(5, 6, ntasks=3) == []

    def test_parse_task_ids(self):
        assert totest.parse_task_ids("1,2,5-7") == [1, 2, 5, 6, 7]
        assert totest.parse_task_ids("3,,4") == [3, 4]


class TestRunTask:

    def test_run_task_with_log(self, tmp_path):
        log = str(tmp_path / "logs" / "task_{job_id}_{task_id:03d}.out")
        code = totest.run_task(PRINT_ENV, 4, "42", 10, log_file=log)
        assert code == 0
        output = (tmp_path / "logs" / "task_42_004.out").read_text()
        assert output.strip() == "42 4 10 0 ['4', '42']"

    def test_run_task_without_log(self, capfd):
        code = totest.run_task(PRINT_ENV, 1, "5", 2)
        assert code == 0
        assert "5 1 2 0 ['1', '5']" in capfd.readouterr().out

    def test_run_task_log_in_cwd(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        code = totest.run_task([sys.executable, "-c", "exit(3)"], 0, "1", 1,
                               log_file="task_{task_id}.out")
        assert code == 3
        assert (tmp_path / "task_0.out").exists()


class TestMain:

    def test_parse_args_errors(self):
        with pytest.raises(SystemExit):
            totest.parse_args(["--nranks", "1", "--ntasks", "2"])
        with pytest.raises(SystemExit):
            totest.parse_args(["--nranks", "1", "--", "echo"])

    def test_parse_args(self):
        args = totest.parse_args(["--nranks", "2", "--task-ids", "1,3",
                                  "--", "echo", "{task_id}"])
        assert args.task_ids == [1, 3]
        assert args.command == ["echo", "{task_id}"]
        args = totest.parse_args(["--nranks", "2", "--ntasks", "3",
                                  "echo"])
        assert args.command == ["echo"]

    def test_main_ranks(self, tmp_path, clean_env):
        log = str(tmp_path / "t_{task_id}.out")
        clean_env.setenv('PBS_JOBID', '777.crux-pbs')
        for rank in range(2):
            clean_env.setenv('PALS_RANKID', str(rank))
            assert totest.main(["--nranks", "2", "--ntasks", "3", "--log",
                                log, "--"] + PRINT_ENV) == 0
        for task in range(3):
            output = (tmp_path / f"t_{task}.out").read_text()
            assert output.startswith(f"777 {task} 3 0")

    def test_main_task_ids(self, tmp_path, clean_env):
        log = str(tmp_path / "t_{task_id}.out")
        assert totest.main(["--nranks", "1", "--rank", "0", "--task-ids",
                            "2,5", "--log", log, "--"] + PRINT_ENV) == 0
        assert (tmp_path / "t_5.out").read_text().startswith("0 5 6 0")
        assert totest.main(["--nranks", "1", "--task-ids", "2",
                            "--task-count", "10", "--log", log, "--"]
                           + PRINT_ENV) == 0
        assert (tmp_path / "t_2.out").read_text().startswith("0 2 10 0")
        # no tasks given: count 0
        assert totest.main(["--nranks", "1", "--task-ids", ",", "--"]
                           + PRINT_ENV) == 0

    def test_main_failures(self, tmp_path, clean_env, capsys):
        # cannot determine the rank
        assert totest.main(["--nranks", "2", "--ntasks", "2", "--",
                            "true"]) == 1
        assert "cannot determine the MPI rank" in capsys.readouterr().err
        # failing task
        assert totest.main(["--nranks", "1", "--ntasks", "1", "--log",
                            str(tmp_path / "x.out"), "--", sys.executable,
                            "-c", "exit(2)"]) == 1
        assert "task 0 failed with exit code 2" in capsys.readouterr().err


def test_shell_braces_are_kept(tmp_path):
    """Only POSYDON placeholders are replaced in commands."""
    log = str(tmp_path / "t.out")
    code = totest.run_task(["echo", "${FILE}", "{task_id:03d}"], 5, "1", 6,
                           log_file=log)
    assert code == 0
    assert open(log).read().strip() == "${FILE} 005"
