"""Unit tests of posydon/CLI/schedulers/base.py and __init__.py"""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>",
]

import pytest

import posydon.CLI.schedulers as schedulers
from posydon.CLI.schedulers.base import ArraySpec, JobSpec, Scheduler


class DummyScheduler(Scheduler):
    """Minimal scheduler to test the base class."""

    name = 'dummy'
    extension = 'job'
    submit_program = 'submit'
    options = {'count': (int, 1, "a number"),
               'flag': (bool, False, "a flag"),
               'text': (str, None, "a text")}

    def directives(self, job):
        return [f"#DUMMY {job.name}"]

    def body(self, job):
        return job.commands

    def submit_line(self, variable, script, after=None, condition='ok'):
        self._check_condition(condition)
        deps = self._dependency_variables(after)
        return f"{variable}=$(submit {script} {deps})\n"

    def mpi_launcher(self, n_ranks, ranks_per_node=None, cpus_per_rank=1):
        return f"run -n {n_ranks}"


class TestArraySpec:

    def test_range(self):
        array = ArraySpec(command="cmd", n_tasks=3)
        assert array.ids == [0, 1, 2]
        assert array.is_range
        assert array.size == 3

    def test_task_ids(self):
        array = ArraySpec(command="cmd", task_ids=[5, 1], task_count=10)
        assert array.ids == [1, 5]
        assert not array.is_range
        assert array.size == 10
        assert ArraySpec(command="cmd", task_ids=[2, 4]).size == 5

    def test_errors(self):
        with pytest.raises(ValueError, match="Either n_tasks or task_ids"):
            ArraySpec(command="cmd")
        with pytest.raises(ValueError, match="at least one task"):
            ArraySpec(command="cmd", n_tasks=0)
        with pytest.raises(ValueError, match="at least one task"):
            ArraySpec(command="cmd", task_ids=[])


class TestJobSpec:

    def test_all_commands(self):
        job = JobSpec(name="job", setup=["s"], commands=["c"],
                      launched_command="l",
                      array=ArraySpec(command="a", n_tasks=1))
        assert job.all_commands() == ["s", "c", "l", "a"]
        assert JobSpec(name="job").all_commands() == []


class TestSchedulerBase:

    def test_options(self):
        scheduler = DummyScheduler(count="3", FLAG="yes", text=None)
        assert scheduler.settings == {'count': 3, 'flag': True, 'text': None}
        assert scheduler.given_options() == {'count': 3, 'flag': True}
        assert repr(scheduler) == \
            "DummyScheduler({'count': 3, 'flag': True})"

    def test_unknown_option(self):
        with pytest.warns(Warning, match="not supported by the dummy"):
            scheduler = DummyScheduler(unknown=1)
        assert scheduler.given_options() == {}

    def test_render_write_submit(self, tmp_path):
        scheduler = DummyScheduler()
        job = JobSpec(name="job", commands=["echo hi"])
        assert scheduler.render(job) == "#!/bin/bash\n#DUMMY job\n\necho hi\n"
        assert scheduler.file_name("a") == "a.job"
        path = str(tmp_path / "a.job")
        assert scheduler.write(job, path) == path
        with pytest.warns(Warning, match="Replace"):
            scheduler.write(job, path)
        assert open(path).read() == scheduler.render(job)
        assert scheduler.submit_command("a.job") == ["submit", "a.job"]
        assert scheduler.submit_line("X", "a", after="Y") == \
            "X=$(submit a ['Y'])\n"
        with pytest.raises(ValueError, match="Unknown dependency condition"):
            scheduler.submit_line("X", "a", condition="maybe")


class TestRegistry:

    def test_names(self):
        assert 'slurm' in schedulers.scheduler_names()
        assert schedulers.normalize_scheduler_name(' SLURM ') == 'slurm'
        assert schedulers.normalize_scheduler_name(None) == 'slurm'
        assert schedulers.normalize_scheduler_name('x', default=None) is None
        assert schedulers.normalize_scheduler_name(3, default=None) is None

    def test_get_scheduler(self):
        scheduler = schedulers.get_scheduler('slurm', exclude='n1')
        assert isinstance(scheduler, schedulers.SlurmScheduler)
        assert scheduler.settings['exclude'] == 'n1'
        with pytest.raises(ValueError, match="Unknown scheduler 'lsf'"):
            schedulers.get_scheduler('lsf')

    def test_options_from_mapping(self):
        mapping = {'EXCLUDE': 'n1', 'ACCOUNT': 'a', 'launcher': 'srun -l'}
        assert schedulers.scheduler_options_from_mapping('slurm', mapping) \
            == {'exclude': 'n1', 'launcher': 'srun -l'}

    def test_parse_option_strings(self):
        assert schedulers.parse_option_strings(None) == {}
        assert schedulers.parse_option_strings(
            ["Exclude=n1,n2", "select=system=crux"]) == \
            {'exclude': 'n1,n2', 'select': 'system=crux'}
        with pytest.raises(ValueError, match="KEY=VALUE"):
            schedulers.parse_option_strings(["novalue"])
        with pytest.raises(ValueError, match="KEY=VALUE"):
            schedulers.parse_option_strings(["=value"])
