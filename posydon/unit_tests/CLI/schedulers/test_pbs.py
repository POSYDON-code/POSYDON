"""Unit tests of posydon/CLI/schedulers/pbs.py"""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>",
]

import pytest

import posydon.CLI.schedulers.pbs as totest
from posydon.CLI.schedulers.base import ArraySpec, JobSpec

CRUX = dict(select="system=crux", filesystems="home:eagle")


class TestHelpers:

    @pytest.mark.parametrize("name, expected", [
        ("1e+00_popsyn", "p1e_00_popsyn"), ("mesa_grid", "mesa_grid"),
        ("", "p"), ("a b", "a_b")])
    def test_pbs_job_name(self, name, expected):
        assert totest.pbs_job_name(name) == expected

    @pytest.mark.parametrize("mail_type, expected", [
        (None, 'a'), ('ALL', 'abe'), ('FAIL', 'a'), ('BEGIN,END', 'be'),
        (['END', 'BEGIN'], 'be'), ('NONE', 'n'), ('TIME_LIMIT_80', 'a'),
        ('unknown', 'a')])
    def test_mail_type_to_pbs(self, mail_type, expected):
        assert totest.mail_type_to_pbs(mail_type) == expected


class TestTasksPerNode:

    def test_defaults(self):
        scheduler = totest.PBSScheduler()
        assert scheduler.tasks_per_node(JobSpec(name="a")) == 128
        assert scheduler.tasks_per_node(
            JobSpec(name="a", cpus_per_task=4)) == 32
        # limited by memory: 256 GB / 9 GB
        assert scheduler.tasks_per_node(
            JobSpec(name="a", mem_per_cpu="9G")) == 28
        assert scheduler.tasks_per_node(
            JobSpec(name="a", cpus_per_task=4, mem_per_cpu="4G")) == 16
        assert scheduler.tasks_per_node(
            JobSpec(name="a", mem_per_cpu="0G")) == 128
        assert scheduler.tasks_per_node(
            JobSpec(name="a", mem_per_cpu="512G")) == 1

    def test_options(self):
        assert totest.PBSScheduler(tasks_per_node="0").tasks_per_node(
            JobSpec(name="a")) == 1
        scheduler = totest.PBSScheduler(cores_per_node=64,
                                        memory_per_node="128G")
        assert scheduler.tasks_per_node(
            JobSpec(name="a", mem_per_cpu="4G")) == 32

    def test_array_layout(self):
        scheduler = totest.PBSScheduler()
        job = JobSpec(name="a", mem_per_cpu="9G",
                      array=ArraySpec(command="c", n_tasks=50))
        assert scheduler.array_layout(job) == (50, 28, 2)
        job.array.max_concurrent = 10
        assert scheduler.array_layout(job) == (10, 10, 1)
        job = JobSpec(name="a", cpus_per_task=4, mem_per_cpu="4G",
                      array=ArraySpec(command="c", n_tasks=1000))
        assert totest.PBSScheduler(max_nodes=4).array_layout(job) == \
            (64, 16, 4)
        assert totest.PBSScheduler(max_nodes=0).array_layout(job) == \
            (16, 16, 1)


class TestPBSDirectives:

    def test_full_directives(self):
        scheduler = totest.PBSScheduler(**CRUX)
        job = JobSpec(name="1e+00_merge", walltime="2-00:00:00", nodes=2,
                      account="SkAI", queue="workq", email="me@x.edu",
                      mail_type="ALL", output="logs/job_{job_id}.out")
        assert scheduler.directives(job) == [
            "#PBS -N p1e_00_merge", "#PBS -l select=2:system=crux",
            "#PBS -l place=scatter", "#PBS -l walltime=48:00:00",
            "#PBS -l filesystems=home:eagle", "#PBS -A SkAI",
            "#PBS -q workq", "#PBS -V", "#PBS -j oe",
            "#PBS -o logs/job_job.out", "#PBS -m abe",
            "#PBS -M me@x.edu"]
        # no comments on directive lines (not allowed by PBS)
        for line in scheduler.directives(job):
            assert "#" not in line[1:]

    def test_minimal_directives(self):
        scheduler = totest.PBSScheduler(place="", export_env="no")
        assert scheduler.directives(JobSpec(name="a", email="x")) == [
            "#PBS -N a", "#PBS -l select=1", "#PBS -j oe"]

    def test_array_nodes(self):
        scheduler = totest.PBSScheduler(**CRUX)
        job = JobSpec(name="a", mem_per_cpu="9G",
                      array=ArraySpec(command="c", n_tasks=50))
        assert "#PBS -l select=2:system=crux" in scheduler.directives(job)


class TestPBSBody:

    def test_array_body(self):
        scheduler = totest.PBSScheduler()
        job = JobSpec(name="a", walltime="1:00:00", mem_per_cpu="9G",
                      setup=["export A={job_id}"],
                      array=ArraySpec(
                          command="run --index {task_id} --end {job_end} 'a b'",
                          n_tasks=50,
                          log_file="./logs/popsyn_{job_id}_{task_id}.out"))
        assert scheduler.body(job) == [
            'cd "${PBS_O_WORKDIR:-.}"',
            "export POSYDON_JOB_END=$(( $(date +%s) + 3600 ))",
            "export A=${PBS_JOBID%%.*}",
            "export POSYDON_ARRAY_JOB_ID=${PBS_JOBID%%.*}",
            "export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}",
            "mpiexec -n 50 --ppn 28 --depth 1 --cpu-bind depth python -m "
            "posydon.CLI.schedulers.task_runner --nranks 50 --ntasks 50 "
            "--task-count 50 --log './logs/popsyn_{job_id}_{task_id}.out' "
            "-- run --index {task_id} --end $POSYDON_JOB_END 'a b'"]

    def test_rescue_body(self):
        scheduler = totest.PBSScheduler(
            launcher="srun -n {ntasks} -N {nodes} -c {depth}")
        job = JobSpec(name="a", cpus_per_task=4, array=ArraySpec(
            command="run", task_ids=[7, 3], task_count=10, launched=False))
        assert scheduler.body(job)[-1] == (
            "srun -n 2 -N 1 -c 4 python -m posydon.CLI.schedulers.task_runner"
            " --nranks 2 --task-ids 3,7 --task-count 10 -- run")

    def test_single_job_body(self):
        job = JobSpec(name="a", commands=["compress-mesa ."],
                      launched_command="python merge.py")
        assert totest.PBSScheduler().body(job) == [
            'cd "${PBS_O_WORKDIR:-.}"', "compress-mesa .", "python merge.py"]

    def test_job_end_needs_walltime(self):
        job = JobSpec(name="a", commands=["run {job_end}"])
        with pytest.raises(ValueError, match="need a walltime"):
            totest.PBSScheduler().body(job)


class TestPBSSubmission:

    def test_submit_line(self):
        scheduler = totest.PBSScheduler()
        assert scheduler.submit_line("A", "a.pbs") == "A=$(qsub a.pbs)\n"
        assert scheduler.submit_line("B", "b.pbs", after="A") == \
            "B=$(qsub -W depend=afterok:${A} b.pbs)\n"
        assert scheduler.submit_line("C", "c.pbs", after=["A", "B"],
                                     condition="any") == \
            "C=$(qsub -W depend=afterany:${A}:${B} c.pbs)\n"
        assert scheduler.submit_command("a.pbs") == ["qsub", "a.pbs"]
        assert scheduler.file_name("cleanup") == "cleanup.pbs"

    def test_mpi_launcher(self):
        scheduler = totest.PBSScheduler()
        assert scheduler.mpi_launcher(16, 8, 4) == \
            "mpiexec -n 16 --ppn 8 --depth 4 --cpu-bind depth"
        assert scheduler.mpi_launcher(4, cpus_per_rank=None) == \
            "mpiexec -n 4 --ppn 4 --depth 1 --cpu-bind depth"
