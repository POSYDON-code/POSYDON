"""Unit tests of posydon/CLI/schedulers/utils.py"""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>",
]

import pytest

import posydon.CLI.schedulers.utils as totest


class TestPlaceholders:

    def test_replace_values_and_callables(self):
        text = "run {task_id} log_{task_id:04d}_{job_id}.out end={job_end}"
        assert totest.replace_placeholders(
            text, task_id=7, job_id="42", job_end="$END") == \
            "run 7 log_0007_42.out end=$END"
        assert totest.replace_placeholders(
            "{task_id:04d}", task_id=lambda spec: f"<{spec}>") == "<:04d>"

    def test_other_braces_untouched(self):
        text = "chgrp ${FILE} {other} {task_id}"
        assert totest.replace_placeholders(text, job_id=1) == text
        assert totest.replace_placeholders(text, task_id=3) == \
            "chgrp ${FILE} {other} 3"

    def test_has_placeholder(self):
        assert totest.has_placeholder("a {job_end} b", "job_end")
        assert totest.has_placeholder("{task_id:03d}", "task_id")
        assert not totest.has_placeholder("{task_id}", "job_end")
        assert not totest.has_placeholder(None, "task_id")


class TestWalltime:

    @pytest.mark.parametrize("walltime, seconds", [
        ("30", 1800), ("30:15", 1815), ("23:00:00", 82800),
        ("2-00:00:00", 172800), ("1-12", 129600), ("1-01:30", 91800),
        ("'2-00:00:00'", 172800), ('"48:00:00"', 172800), (90, 5400)])
    def test_walltime_to_seconds(self, walltime, seconds):
        assert totest.walltime_to_seconds(walltime) == seconds

    def test_walltime_errors(self):
        with pytest.raises(ValueError, match="Cannot interpret walltime"):
            totest.walltime_to_seconds("one day")
        with pytest.raises(TypeError, match="walltime must be a string"):
            totest.walltime_to_seconds(None)
        with pytest.raises(TypeError, match="walltime must be a string"):
            totest.walltime_to_seconds(True)

    @pytest.mark.parametrize("walltime, expected", [
        ("2-00:00:00", "48:00:00"), ("23:00:00", "23:00:00"),
        ("14", "00:14:00"), ("0-00:00:05", "00:00:05")])
    def test_walltime_to_hms(self, walltime, expected):
        assert totest.walltime_to_hms(walltime) == expected


class TestMemory:

    @pytest.mark.parametrize("memory, gb", [
        ("4G", 4.0), ("9g", 9.0), ("4GB", 4.0), ("2048M", 2.0),
        ("1T", 1024.0), ("1048576K", 1.0), ("1024", 1.0), (2048, 2.0),
        ("'4G'", 4.0)])
    def test_memory_to_gb(self, memory, gb):
        assert totest.memory_to_gb(memory) == pytest.approx(gb)

    def test_memory_errors(self):
        with pytest.raises(ValueError, match="Cannot interpret memory"):
            totest.memory_to_gb("lots")
        with pytest.raises(TypeError, match="memory must be a string"):
            totest.memory_to_gb(None)


@pytest.mark.parametrize("value, expected", [
    ("True", True), ("yes", True), ("1", True), ("false", False),
    ("no", False), (True, True), (0, False)])
def test_to_bool(value, expected):
    assert totest.to_bool(value) is expected
