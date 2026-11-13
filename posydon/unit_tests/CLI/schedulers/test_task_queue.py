"""Unit tests of posydon/CLI/schedulers/task_queue.py"""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>",
]

import os
import subprocess
import sys

import pytest

import posydon.CLI.schedulers.task_queue as totest


def counter_files(directory):
    return sorted(name for name in os.listdir(directory)
                  if name.startswith(totest.COUNTER_PREFIX))


class TestTaskQueue:

    def test_single_rank_takes_everything(self, tmp_path):
        queue = totest.TaskQueue(str(tmp_path / "q"), [4, 7, 9, 2], 0, 1)
        assert list(queue.tasks()) == [4, 7, 9, 2]
        assert counter_files(tmp_path / "q") == ["next_4"]

    def test_ranks_share_the_tasks(self, tmp_path):
        directory = str(tmp_path / "q")
        ranks = [totest.TaskQueue(directory, list(range(7)), rank, 3)
                 for rank in range(3)]
        tasks = [queue.tasks() for queue in ranks]
        # each rank starts with its own task; rank 0 creates the counter
        assert [next(t) for t in tasks] == [0, 1, 2]
        assert counter_files(directory) == ["next_3"]
        # then they take the next free task, whoever asks first
        assert next(tasks[2]) == 3
        assert next(tasks[0]) == 4
        assert next(tasks[2]) == 5
        assert next(tasks[1]) == 6
        assert [list(t) for t in tasks] == [[], [], []]

    def test_no_queue_with_enough_ranks(self, tmp_path):
        directory = tmp_path / "q"
        assert list(totest.TaskQueue(str(directory), [5, 6], 0, 4).tasks()) \
            == [5]
        assert list(totest.TaskQueue(str(directory), [5, 6], 3, 4).tasks()) \
            == []
        assert not directory.exists()

    def test_counter_is_created_once(self, tmp_path):
        directory = tmp_path / "q"
        queue = totest.TaskQueue(str(directory), list(range(10)), 0, 2)
        queue.create_counter()
        assert queue.claim() == 2
        (directory / "other_file").write_text("")
        (directory / "next_x").write_text("")
        queue.create_counter()
        assert counter_files(directory) == ["next_3", "next_x"]
        assert queue._counter() == 3

    def test_lost_race(self, tmp_path, monkeypatch):
        queue = totest.TaskQueue(str(tmp_path), list(range(5)), 0, 2)
        queue.create_counter()
        rename = os.rename
        calls = []

        def rename_once_too_late(src, dst):
            calls.append(src)
            if len(calls) == 1:
                # another rank claims position 2 just before us
                rename(src, dst)
                raise FileNotFoundError(src)
            rename(src, dst)

        monkeypatch.setattr(totest.os, "rename", rename_once_too_late)
        monkeypatch.setattr(totest.time, "sleep", lambda seconds: None)
        assert queue.claim() == 3
        assert len(calls) == 2

    def test_waits_for_counter(self, tmp_path, monkeypatch):
        directory = str(tmp_path / "q")
        rank1 = totest.TaskQueue(directory, list(range(4)), 1, 2, poll=0)
        sleeps = []

        def sleep(seconds):
            sleeps.append(seconds)
            if len(sleeps) == 2:
                totest.TaskQueue(directory, list(range(4)), 0, 2)\
                    .create_counter()

        monkeypatch.setattr(totest.time, "sleep", sleep)
        assert rank1.claim() == 2
        assert len(sleeps) == 2

    def test_timeout(self, tmp_path):
        queue = totest.TaskQueue(str(tmp_path / "q"), list(range(4)), 1, 2,
                                 timeout=0, poll=0)
        tasks = queue.tasks()
        assert next(tasks) == 1
        with pytest.raises(TimeoutError, match="no task counter"):
            next(tasks)


def test_parallel_ranks_run_each_task_once(tmp_path):
    """Several task runner processes share the tasks without duplicates."""
    done = tmp_path / "done"
    done.mkdir()
    # each task creates its marker file exclusively: running a task twice
    # makes the second run fail
    task = ("import os, random, sys, time; "
            "time.sleep(random.uniform(0, 0.05)); "
            f"os.close(os.open(os.path.join({str(done)!r}, sys.argv[1]), "
            "os.O_CREAT | os.O_EXCL))")
    env = dict(os.environ, POSYDON_QUEUE_ID="test")
    for var in ("PALS_RANKID", "PMI_RANK", "PMIX_RANK",
                "OMPI_COMM_WORLD_RANK", "MV2_COMM_WORLD_RANK",
                "SLURM_PROCID"):
        env.pop(var, None)
    processes = [subprocess.Popen(
        [sys.executable, "-m", "posydon.CLI.schedulers.task_runner",
         "--rank", str(rank), "--nranks", "4", "--ntasks", "40",
         "--queue", str(tmp_path / "queue"), "--",
         sys.executable, "-c", task, "{task_id}"], env=env)
        for rank in range(4)]
    assert [process.wait(timeout=120) for process in processes] == [0] * 4
    assert sorted(int(name) for name in os.listdir(done)) == list(range(40))
    assert os.listdir(tmp_path / "queue" / "test") == ["next_40"]
