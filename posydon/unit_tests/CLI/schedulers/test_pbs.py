"""Unit tests of posydon/CLI/schedulers/pbs.py"""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>",
]

import pytest

import posydon.CLI.schedulers.pbs as totest
from posydon.CLI.schedulers.base import ArraySpec, JobSpec

CRUX = dict(site="crux")
# a Crux-sized node without the site's other settings
NODE = dict(cores_per_node=128, memory_per_node="256G", chunk_resources=False)
SEARCH = ("POSYDON_MPIEXEC=$(ls /opt/cray/pals/*/bin/mpiexec 2>/dev/null "
          "| sort -V | tail -n 1)")
MPIEXEC = '"${POSYDON_MPIEXEC:-mpiexec}"'


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
        scheduler = totest.PBSScheduler(**NODE)
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
        scheduler = totest.PBSScheduler(**NODE)
        job = JobSpec(name="a", mem_per_cpu="9G",
                      array=ArraySpec(command="c", n_tasks=50))
        assert scheduler.array_layout(job) == (50, 28, 2)
        job.array.max_concurrent = 10
        assert scheduler.array_layout(job) == (10, 10, 1)
        job = JobSpec(name="a", cpus_per_task=4, mem_per_cpu="4G",
                      array=ArraySpec(command="c", n_tasks=1000))
        assert totest.PBSScheduler(max_nodes=4, **NODE).array_layout(job) \
            == (64, 16, 4)
        assert totest.PBSScheduler(max_nodes=0, **NODE).array_layout(job) \
            == (16, 16, 1)


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
        scheduler = totest.PBSScheduler(place="", export_env="no",
                                        chunk_resources=False)
        assert scheduler.directives(JobSpec(name="a", email="x")) == [
            "#PBS -N a", "#PBS -l select=1", "#PBS -j oe"]

    def test_array_nodes(self):
        scheduler = totest.PBSScheduler(**CRUX)
        job = JobSpec(name="a", mem_per_cpu="9G",
                      array=ArraySpec(command="c", n_tasks=50))
        assert "#PBS -l select=2:system=crux" in scheduler.directives(job)


class TestPBSNewOptions:

    def test_rerunnable(self):
        job = JobSpec(name="a")
        assert "#PBS -r y" in totest.PBSScheduler(rerunnable=True)\
            .directives(job)
        assert "#PBS -r n" in totest.PBSScheduler(rerunnable="no")\
            .directives(job)
        assert not any("-r " in line for line in
                       totest.PBSScheduler().directives(job))

    def test_task_stats(self):
        job = JobSpec(name="a", array=ArraySpec(command="run", n_tasks=2))
        assert "--stats stats/run1 --" in totest.PBSScheduler(
            task_stats="stats/run1", **NODE).body(job)[-1]
        assert "--stats" not in totest.PBSScheduler(**NODE).body(job)[-1]


class TestPBSWorkQueue:

    def grid_job(self, n_tasks, cpus=2):
        return JobSpec(name="a", cpus_per_task=cpus, array=ArraySpec(
            command="run {task_id}", n_tasks=n_tasks))

    def test_dynamic_queue(self):
        scheduler = totest.PBSScheduler(tasks_per_node=64, max_nodes=2)
        lines = scheduler.body(self.grid_job(1000))
        assert lines[-3] == \
            "export POSYDON_QUEUE_ID=${PBS_JOBID%%.*}_$(date +%s)"
        assert lines[-2] == "export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}"
        assert lines[-1].startswith(MPIEXEC + " -n 128 --ppn 64 --depth 2 ")
        assert "--nranks 128 --ntasks 1000 --task-count 1000 " \
            "--queue .task_queue -- run {task_id}" in lines[-1]
        scheduler = totest.PBSScheduler(tasks_per_node=64, max_nodes=2,
                                        queue_dir="/eagle/q")
        assert "--queue /eagle/q --" in scheduler.body(self.grid_job(1000))[-1]

    def test_no_queue_needed(self):
        # one rank per task: nothing to hand out
        lines = totest.PBSScheduler(tasks_per_node=64).body(self.grid_job(64))
        assert "--queue" not in lines[-1]
        assert not any("POSYDON_QUEUE_ID" in line for line in lines)

    def test_static_order(self):
        scheduler = totest.PBSScheduler(tasks_per_node=64, max_nodes=2,
                                        task_order="Static")
        assert scheduler.settings['task_order'] == "static"
        assert "--queue" not in scheduler.body(self.grid_job(1000))[-1]
        with pytest.raises(ValueError, match="task_order"):
            totest.PBSScheduler(task_order="random")

    @pytest.mark.parametrize("ppn, cpus, depth", [
        (56, 1, 2), (28, 1, 4), (64, 2, 2), (32, 4, 4), (16, 4, 8),
        (128, 1, 1), (64, 4, 4)])
    def test_depth(self, ppn, cpus, depth):
        scheduler = totest.PBSScheduler(tasks_per_node=ppn, cores_per_node=128)
        assert scheduler.depth(JobSpec(name="a", cpus_per_task=cpus),
                               ppn) == depth
        lines = scheduler.body(self.grid_job(ppn, cpus))
        assert f"--depth {depth} " in lines[-1]
        # threads of a task: its cores, not the spacing
        assert lines[-2] == \
            f"export OMP_NUM_THREADS=${{OMP_NUM_THREADS:-{cpus}}}"


class TestPBSBody:

    def test_array_body(self):
        scheduler = totest.PBSScheduler(**NODE)
        job = JobSpec(name="a", walltime="1:00:00", mem_per_cpu="9G",
                      setup=["export A={job_id}"],
                      array=ArraySpec(
                          command="run --index {task_id} --end {job_end} 'a b'",
                          n_tasks=50,
                          log_file="./logs/popsyn_{job_id}_{task_id}.out"))
        assert scheduler.body(job) == [
            'cd "${PBS_O_WORKDIR:-.}"',
            SEARCH,
            "export POSYDON_JOB_END=$(( $(date +%s) + 3600 ))",
            "export A=${PBS_JOBID%%.*}",
            "export POSYDON_ARRAY_JOB_ID=${PBS_JOBID%%.*}",
            "export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}",
            MPIEXEC + " -n 50 --ppn 28 --depth 4 --cpu-bind depth python -m "
            "posydon.CLI.schedulers.task_runner --nranks 50 --ntasks 50 "
            "--task-count 50 --log './logs/popsyn_{job_id}_{task_id}.out' "
            "-- run --index {task_id} --end $POSYDON_JOB_END 'a b'"]

    def test_rescue_body(self):
        scheduler = totest.PBSScheduler(
            launcher="srun -n {ntasks} -N {nodes} -c {depth}", **NODE)
        job = JobSpec(name="a", cpus_per_task=4, array=ArraySpec(
            command="run", task_ids=[7, 3], task_count=10, launched=False))
        assert scheduler.body(job)[-1] == (
            "srun -n 2 -N 1 -c 64 python -m posydon.CLI.schedulers.task_runner"
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
        fail = ' || {{ echo "Submission of {0} failed" >&2; exit 1; }}\n'
        assert scheduler.submit_line("A", "a.pbs") == \
            "A=$(qsub a.pbs)" + fail.format("a.pbs")
        assert scheduler.submit_line("B", "b.pbs", after="A") == \
            "B=$(qsub -W depend=afterok:${A} b.pbs)" + fail.format("b.pbs")
        assert scheduler.submit_line("C", "c.pbs", after=["A", "B"],
                                     condition="any") == \
            "C=$(qsub -W depend=afterany:${A}:${B} c.pbs)" \
            + fail.format("c.pbs")
        assert scheduler.submit_command("a.pbs") == ["qsub", "a.pbs"]
        assert scheduler.file_name("cleanup") == "cleanup.pbs"

    def test_submission_script_stops_on_failure(self, tmp_path):
        """A failed submission stops the bash script."""
        import subprocess
        scheduler = totest.PBSScheduler()
        fake = tmp_path / "qsub"
        fake.write_text("#!/bin/bash\nexit 1\n")
        fake.chmod(0o755)
        script = tmp_path / "submit.sh"
        script.write_text("#!/bin/bash\n" + scheduler.submit_line("A", "a.pbs")
                          + "echo reached\n")
        result = subprocess.run(["bash", str(script)], capture_output=True,
                                text=True,
                                env={"PATH": f"{tmp_path}:/usr/bin:/bin"})
        assert result.returncode == 1
        assert "reached" not in result.stdout
        assert "Submission of a.pbs failed" in result.stderr

    def test_mpi_launcher(self):
        scheduler = totest.PBSScheduler()
        assert scheduler.mpi_launcher(16, 8, 4) == \
            MPIEXEC + " -n 16 --ppn 8 --depth 4 --cpu-bind depth"
        assert scheduler.mpi_launcher(4, cpus_per_rank=None) == \
            MPIEXEC + " -n 4 --ppn 4 --depth 1 --cpu-bind depth"
        assert totest.PBSScheduler(mpiexec_search="").mpi_launcher(4) == \
            "mpiexec -n 4 --ppn 4 --depth 1 --cpu-bind depth"


class TestFindMpiexec:

    def array_job(self):
        return JobSpec(name="a", array=ArraySpec(command="run", n_tasks=2))

    def test_search_line_only_when_used(self):
        lines = totest.PBSScheduler(**NODE).body(self.array_job())
        assert lines[1] == SEARCH
        assert lines[-1].startswith(MPIEXEC + " -n 2 ")
        # MPI job using the launcher in its command (e.g. MESA grids)
        scheduler = totest.PBSScheduler()
        job = JobSpec(name="a", commands=[
            scheduler.mpi_launcher(8) + " python -m mpi4py run.py"])
        assert scheduler.body(job)[1] == SEARCH
        # jobs without the launcher, or with an explicit mpiexec
        assert SEARCH not in scheduler.body(JobSpec(name="a",
                                                    commands=["echo"]))
        explicit = totest.PBSScheduler(
            launcher="/opt/cray/pals/1.8/bin/mpiexec -n {ntasks} --ppn {ppn}",
            **NODE)
        lines = explicit.body(self.array_job())
        assert SEARCH not in lines
        assert lines[-1].startswith("/opt/cray/pals/1.8/bin/mpiexec -n 2 ")

    def test_no_search(self):
        lines = totest.PBSScheduler(mpiexec_search="", **NODE).body(
            self.array_job())
        assert not any("POSYDON_MPIEXEC" in line for line in lines)
        assert lines[-1].startswith("mpiexec -n 2 ")

    def test_custom_search(self):
        scheduler = totest.PBSScheduler(mpiexec_search="/sw/mpi/*/mpiexec",
                                        **NODE)
        assert scheduler.body(self.array_job())[1] == (
            "POSYDON_MPIEXEC=$(ls /sw/mpi/*/mpiexec 2>/dev/null | sort -V "
            "| tail -n 1)")

    def test_search_in_bash(self, tmp_path):
        """The job script uses the latest system mpiexec, else PATH's."""
        import subprocess

        def fake(path, name):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f"#!/bin/bash\necho {name}\n")
            path.chmod(0o755)

        fake(tmp_path / "conda" / "bin" / "mpiexec", "conda")
        pattern = str(tmp_path / "pals" / "*" / "bin" / "mpiexec")
        scheduler = totest.PBSScheduler(mpiexec_search=pattern,
                                        launcher="{mpiexec}")
        script = "\n".join(scheduler.body(JobSpec(
            name="a", commands=[scheduler.mpi_launcher(1)])))
        env = {"PATH": f"{tmp_path / 'conda' / 'bin'}:/usr/bin:/bin"}

        def run():
            return subprocess.run(["bash", "-c", script], env=env,
                                  capture_output=True, text=True,
                                  cwd=tmp_path).stdout.strip()

        assert run() == "conda"           # nothing found: mpiexec in PATH
        for version in ("1.8", "1.10", "1.9.2"):
            fake(tmp_path / "pals" / version / "bin" / "mpiexec",
                 f"pals-{version}")
        assert run() == "pals-1.10"       # latest version, not the last
                                          # in alphabetical order


class TestSitesAndPortability:

    def test_crux_site(self):
        scheduler = totest.PBSScheduler(site=" 'CRUX' ")
        settings = scheduler.settings
        assert settings['site'] == "crux"
        assert (settings['cores_per_node'], settings['memory_per_node']) == \
            (128, "256G")
        assert (settings['select'], settings['filesystems']) == \
            ("system=crux", "home:eagle")
        assert settings['chunk_resources'] is False
        assert settings['mpiexec_search'] == totest.PALS_MPIEXEC_SEARCH
        job = JobSpec(name="a", mem_per_cpu="4.5G",
                      array=ArraySpec(command="c", n_tasks=56))
        assert "#PBS -l select=1:system=crux" in scheduler.directives(job)
        assert "--ppn 56 --depth 2 " in scheduler.body(job)[-1]

    def test_explicit_options_win(self):
        scheduler = totest.PBSScheduler(site="crux", cores_per_node=64,
                                        select="system=crux:foo=1")
        assert scheduler.settings['cores_per_node'] == 64
        assert scheduler.settings['select'] == "system=crux:foo=1"
        assert scheduler.settings['memory_per_node'] == "256G"

    def test_unknown_site(self):
        with pytest.raises(ValueError, match="Unknown site 'nowhere'"):
            totest.PBSScheduler(site="nowhere")

    def test_node_size_needed(self):
        job = JobSpec(name="a", array=ArraySpec(command="c", n_tasks=4))
        with pytest.raises(ValueError, match="give cores_per_node, a site"):
            totest.PBSScheduler().body(job)
        job.mem_per_cpu = "4G"
        with pytest.raises(ValueError, match="give memory_per_node, a site"):
            totest.PBSScheduler(cores_per_node=64).body(job)
        with pytest.raises(ValueError, match="cores_per_node and memory"):
            totest.PBSScheduler().body(job)
        # explicit tasks per node: no node size needed, no spreading
        lines = totest.PBSScheduler(tasks_per_node=4).body(job)
        assert "-n 4 --ppn 4 --depth 1 " in lines[-1]

    def test_chunk_resources(self):
        scheduler = totest.PBSScheduler(cores_per_node=64,
                                        memory_per_node="512G",
                                        select="model=x")
        # array job: 64 single-core tasks with 4 GB each per node
        job = JobSpec(name="a", mem_per_cpu="4G",
                      array=ArraySpec(command="c", n_tasks=100))
        assert "#PBS -l select=2:ncpus=64:mpiprocs=64:mem=256gb:model=x" in \
            scheduler.directives(job)
        # spread ranks: the chunk asks for the cores they are bound to
        job = JobSpec(name="a", array=ArraySpec(command="c", n_tasks=16))
        assert "#PBS -l select=1:ncpus=64:mpiprocs=16:model=x" in \
            scheduler.directives(job)
        # single job, e.g. a merge or an MPI job
        job = JobSpec(name="a", nodes=2, tasks_per_node=8, cpus_per_task=2,
                      mem_per_cpu="1.5G")
        assert "#PBS -l select=2:ncpus=16:mpiprocs=8:mem=24gb:model=x" in \
            scheduler.directives(job)
        assert "#PBS -l select=1:ncpus=1:mpiprocs=1" in \
            totest.PBSScheduler().directives(JobSpec(name="a"))
        assert "#PBS -l select=1" in totest.PBSScheduler(
            chunk_resources=False).directives(JobSpec(name="a"))

    @pytest.mark.parametrize("name, start", [
        ("pals", '"${POSYDON_MPIEXEC:-mpiexec}" -n 8 --ppn 4 --depth 2 '
                 '--cpu-bind depth'),
        ("MPICH", "mpiexec -n 8 -ppn 4 -bind-to core:2"),
        ("openmpi", "mpiexec -np 8 --map-by ppr:4:node:PE=2 --bind-to core")])
    def test_launchers(self, name, start):
        scheduler = totest.PBSScheduler(launcher=name, tasks_per_node=4,
                                        cores_per_node=8)
        job = JobSpec(name="a", array=ArraySpec(command="c", n_tasks=8))
        lines = scheduler.body(job)
        assert lines[-1].startswith(start + " python -m ")
        # only the PALS launcher looks for Cray's mpiexec by default
        assert any("POSYDON_MPIEXEC=" in line for line in lines) == \
            (name == "pals")

    def test_launcher_search_and_errors(self):
        scheduler = totest.PBSScheduler(launcher="mpich",
                                        mpiexec_search="/opt/mpich/*/mpiexec",
                                        tasks_per_node=2)
        lines = scheduler.body(JobSpec(name="a", array=ArraySpec(
            command="c", n_tasks=2)))
        assert lines[1].startswith("POSYDON_MPIEXEC=$(ls /opt/mpich/*/mpiexec")
        with pytest.raises(ValueError, match="launcher must be one of"):
            totest.PBSScheduler(launcher="srun")
