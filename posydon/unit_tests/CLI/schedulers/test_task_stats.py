"""Unit tests of posydon/CLI/schedulers/task_stats.py"""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>",
]

import types

import pytest

import posydon.CLI.schedulers.task_stats as totest


def usage(maxrss_kb, utime=1.0, stime=0.5):
    return types.SimpleNamespace(ru_maxrss=maxrss_kb, ru_utime=utime,
                                 ru_stime=stime)


@pytest.fixture
def stats_dir(tmp_path):
    """Two ranks of job 7: tasks 0-2 finished (task 2 failed), task 3 was
    still running; plus one task of another job and a rerun of task 0."""
    writer = totest.StatsWriter(str(tmp_path), "7", 0)
    writer.start(0, timestamp=100.0)
    writer.end(0, 50.0, 0, usage(1024**2), timestamp=150.0)        # 1 GB
    writer.start(2, timestamp=150.0)
    writer.end(2, 30.0, 1, usage(2 * 1024**2), timestamp=180.0)    # 2 GB
    writer = totest.StatsWriter(str(tmp_path), "7", 1)
    writer.start(1, timestamp=100.0)
    writer.end(1, 70.0, 0, timestamp=170.0)                        # no usage
    writer.start(3, timestamp=170.0)                               # killed
    other = totest.StatsWriter(str(tmp_path), "8", 0)
    other.start(0)
    other.end(0, 5.0, 0, usage(10))
    return str(tmp_path)


class TestWriteRead:

    def test_files(self, stats_dir):
        assert sorted(p.name for p in __import__('pathlib').Path(stats_dir)
                      .iterdir()) == ["stats_7_rank0.csv", "stats_7_rank1.csv",
                                      "stats_8_rank0.csv"]

    def test_read_stats(self, stats_dir):
        tasks = totest.read_stats(stats_dir, job_id="7")
        assert set(tasks) == {("7", 0), ("7", 1), ("7", 2), ("7", 3)}
        assert tasks[("7", 0)]["elapsed_s"] == 50.0
        assert tasks[("7", 0)]["max_rss_kb"] == 1024**2
        assert tasks[("7", 1)]["max_rss_kb"] is None
        assert tasks[("7", 3)]["end"] is None
        assert len(totest.read_stats(stats_dir)) == 5

    def test_end_without_start(self, tmp_path):
        writer = totest.StatsWriter(str(tmp_path), "9", 0)
        writer.end(4, 1.0, 0)
        tasks = totest.read_stats(str(tmp_path))
        assert tasks[("9", 4)]["start"] is None
        assert tasks[("9", 4)]["elapsed_s"] == 1.0


class TestSummary:

    def test_summarize(self, stats_dir):
        summary = totest.summarize(totest.read_stats(stats_dir, "7"))
        assert summary["tasks"] == 4
        assert summary["finished"] == 3
        assert summary["failed"] == 1
        assert summary["unfinished"] == 1
        assert summary["elapsed_s"] == {"min": 30.0, "median": 50.0,
                                        "max": 70.0}
        assert summary["max_rss_kb"]["max"] == 2 * 1024**2

    def test_summarize_nothing_finished(self):
        summary = totest.summarize({("1", 0): {"end": None}})
        assert summary["elapsed_s"] is None
        assert summary["max_rss_kb"] is None
        assert "elapsed" not in totest.format_summary(summary)

    def test_suggest_tasks_per_node(self):
        # 2 GB * 1.3 = 2.6 GB -> 98 by memory, 128 by cores
        assert totest.suggest_tasks_per_node(2 * 1024**2) == 98
        assert totest.suggest_tasks_per_node(2 * 1024**2,
                                             cpus_per_task=4) == 32
        assert totest.suggest_tasks_per_node(0) == 128
        assert totest.suggest_tasks_per_node(1024**3) == 1

    def test_format_summary(self, stats_dir):
        text = totest.format_summary(totest.summarize(
            totest.read_stats(stats_dir, "7")))
        assert "tasks: 4  finished: 3  failed: 1  unfinished: 1" in text
        assert "max 0:01:10" in text
        assert "max 2.00" in text
        assert "suggested tasks_per_node: 98" in text


class TestMain:

    def test_main(self, stats_dir, capsys):
        assert totest.main([stats_dir, "--job", "7", "--cpus-per-task",
                            "4"]) == 0
        assert "suggested tasks_per_node: 32" in capsys.readouterr().out

    def test_main_no_stats(self, tmp_path, capsys):
        assert totest.main([str(tmp_path)]) == 1
        assert "No task statistics found" in capsys.readouterr().err
