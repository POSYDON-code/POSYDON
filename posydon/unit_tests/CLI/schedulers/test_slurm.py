"""Unit tests of posydon/CLI/schedulers/slurm.py"""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>",
]

import pytest

from posydon.CLI.schedulers.base import ArraySpec, JobSpec
from posydon.CLI.schedulers.slurm import SlurmScheduler


class TestSlurmDirectives:

    def test_all_directives(self):
        scheduler = SlurmScheduler(exclude="n1,n2")
        job = JobSpec(name="job", walltime="1-00:00:00", nodes=2,
                      tasks_per_node=8, cpus_per_task=4, mem_per_cpu="4G",
                      account="acc", queue="part", email="me@x.edu",
                      mail_type="ALL", output="out_{job_id}.log",
                      truncate_output=True)
        assert scheduler.directives(job) == [
            "#SBATCH --job-name=job", "#SBATCH --account=acc",
            "#SBATCH --partition=part", "#SBATCH --nodes=2",
            "#SBATCH --ntasks-per-node=8", "#SBATCH --cpus-per-task=4",
            "#SBATCH --mem-per-cpu=4G", "#SBATCH --time=1-00:00:00",
            "#SBATCH --open-mode=truncate", "#SBATCH --output=out_%j.log",
            "#SBATCH --mail-type=ALL", "#SBATCH --mail-user=me@x.edu",
            "#SBATCH --exclude=n1,n2"]

    def test_minimal_directives(self):
        job = JobSpec(name="job", email="me@x.edu")
        assert SlurmScheduler().directives(job) == ["#SBATCH --job-name=job"]

    def test_array_directives(self):
        scheduler = SlurmScheduler()
        job = JobSpec(name="a", output="ignored.out", array=ArraySpec(
            command="c", n_tasks=10, max_concurrent=3,
            log_file="logs/x_{job_id}_{task_id:04d}.out"))
        directives = scheduler.directives(job)
        assert "#SBATCH --array=0-9%3" in directives
        assert "#SBATCH --output=logs/x_%A_%4a.out" in directives
        job = JobSpec(name="a", output="job.out", array=ArraySpec(
            command="c", task_ids=[5, 2], log_file="x_{task_id:x}.out"))
        directives = scheduler.directives(job)
        assert "#SBATCH --array=2,5" in directives
        assert "#SBATCH --output=x_%a.out" in directives
        job = JobSpec(name="a", output="job.out",
                      array=ArraySpec(command="c", n_tasks=2))
        assert "#SBATCH --output=job.out" in scheduler.directives(job)


class TestSlurmBody:

    def test_array_body(self):
        job = JobSpec(name="a", setup=["export A=1"], commands=["echo {job_id}"],
                      array=ArraySpec(command="run {task_id} {job_end}",
                                      n_tasks=4))
        assert SlurmScheduler().body(job) == [
            "export A=1", "echo $SLURM_ARRAY_JOB_ID",
            "srun run $SLURM_ARRAY_TASK_ID $SLURM_JOB_END_TIME"]

    def test_rescue_body(self):
        job = JobSpec(name="a", array=ArraySpec(
            command="run", task_ids=[1, 3], task_count=10, launched=False))
        assert SlurmScheduler().body(job) == [
            "export SLURM_ARRAY_TASK_COUNT=10",
            "export SLURM_ARRAY_TASK_MIN=0", "run"]
        job = JobSpec(name="a", array=ArraySpec(
            command="run", n_tasks=3, task_count=5, launched=False))
        assert SlurmScheduler().body(job)[0] == \
            "export SLURM_ARRAY_TASK_COUNT=5"

    def test_single_job_body(self):
        job = JobSpec(name="a", commands=["echo {job_id}"],
                      launched_command="python merge.py")
        assert SlurmScheduler(launcher="srun -n 1").body(job) == [
            "echo $SLURM_JOB_ID", "srun -n 1 python merge.py"]

    def test_render(self):
        job = JobSpec(name="a", commands=["echo hi"])
        assert SlurmScheduler().render(job) == \
            "#!/bin/bash\n#SBATCH --job-name=a\n\necho hi\n"


class TestSlurmSubmission:

    def test_submit_line(self):
        scheduler = SlurmScheduler()
        assert scheduler.submit_line("A", "a.slurm") == \
            "A=$(sbatch --parsable a.slurm)\n"
        assert scheduler.submit_line("B", "b.slurm", after="A") == (
            "B=$(sbatch --parsable --dependency=afterok:${A} "
            "--kill-on-invalid-dep=yes b.slurm)\n")
        assert scheduler.submit_line("C", "c.slurm", after=["A", "B"],
                                     condition='any') == (
            "C=$(sbatch --parsable --dependency=afterany:${A}:${B} "
            "--kill-on-invalid-dep=yes c.slurm)\n")
        with pytest.raises(ValueError):
            scheduler.submit_line("C", "c.slurm", condition='notok')

    def test_submit_command_and_mpi(self):
        scheduler = SlurmScheduler()
        assert scheduler.submit_command("a.slurm") == ["sbatch", "a.slurm"]
        assert scheduler.mpi_launcher(8) == "mpirun --bind-to none -np 8"
        assert scheduler.file_name("cleanup") == "cleanup.slurm"
