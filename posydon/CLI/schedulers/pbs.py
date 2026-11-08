"""PBS Pro backend of POSYDON's job scheduler interface.

PBS sites like ALCF Crux schedule whole nodes (jobs never share a node) and
allow only a few jobs or array subjobs per project. Translating each array
task into its own PBS job would occupy a full node for every single-core
task. Array jobs are therefore "packed": a single PBS job requests enough
nodes to run `tasks_per_node` tasks per node and starts them with the MPI
launcher (`mpiexec`). Each launched process runs
`python -m posydon.CLI.schedulers.task_runner`, which works out the array
indices it is responsible for and runs the task command for each of them,
with one log file per task, like a SLURM array task.
"""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>",
]

import math
import re
import shlex

from posydon.CLI.schedulers.base import Scheduler
from posydon.CLI.schedulers.utils import (
    has_placeholder,
    memory_to_gb,
    replace_placeholders,
    walltime_to_hms,
    walltime_to_seconds,
)

# Default launcher: Cray PALS mpiexec as used on ALCF systems.
DEFAULT_LAUNCHER = ("mpiexec -n {ntasks} --ppn {ppn} --depth {depth} "
                    "--cpu-bind depth")

# Shell expression giving the numeric part of the PBS job ID
# (e.g. "12345.crux-pbs-01" -> "12345")
PBS_JOB_NUMBER = '${PBS_JOBID%%.*}'

# Environment variable holding the end of the job (seconds since the epoch)
JOB_END_VARIABLE = 'POSYDON_JOB_END'


def pbs_job_name(name):
    """Sanitize a job name for PBS (no special characters, starts with a
    letter, at most 236 characters)."""
    clean = re.sub(r'[^A-Za-z0-9_.-]', '_', str(name))
    if not clean or not clean[0].isalpha():
        clean = 'p' + clean
    return clean[:236]


def mail_type_to_pbs(mail_type):
    """Convert SLURM mail types (e.g. "ALL", "FAIL", "BEGIN,END") to the
    mail events of PBS' `-m` option ('abe', 'a', 'be', 'n', ...)."""
    if mail_type is None:
        return 'a'
    if isinstance(mail_type, str):
        mail_type = mail_type.split(',')
    events = ''
    for item in mail_type:
        item = str(item).strip().upper()
        if item == 'NONE':
            return 'n'
        for key, value in (('ALL', 'abe'), ('BEGIN', 'b'), ('END', 'e'),
                           ('FAIL', 'a'), ('REQUEUE', 'a'),
                           ('TIME_LIMIT', 'a'), ('INVALID_DEPEND', 'a')):
            if item.startswith(key):
                events += value
    events = ''.join(event for event in 'abe' if event in events)
    return events or 'a'


class PBSScheduler(Scheduler):
    """Job scripts for PBS Pro (`qsub`), e.g. on ALCF Crux.

    Array jobs are packed onto nodes and started with the MPI launcher and
    POSYDON's task runner (see the module documentation).

    Options
    -------
    select : str
        Additional resources of each chunk in "-l select=N:...", e.g.
        "system=crux".
    filesystems : str
        Value of "-l filesystems=" (required at ALCF, e.g. "home:eagle").
    place : str
        Value of "-l place=" (default: "scatter").
    tasks_per_node : int
        Array tasks per node (default: as many as fit on a node given the
        job's cores and memory per task).
    max_nodes : int
        Maximum number of nodes of an array job; tasks beyond the available
        slots run one after another.
    cores_per_node, memory_per_node : int, str
        Size of a node, used for the default `tasks_per_node` (default: 128
        cores and 256G, an ALCF Crux node).
    launcher : str
        Launcher template with the fields {ntasks}, {ppn}, {depth} and
        {nodes} (default: Cray PALS mpiexec).
    export_env : bool
        Export the environment of qsub to the job ("-V"), as SLURM does by
        default (default: True).
    """

    name = 'pbs'
    extension = 'pbs'
    submit_program = 'qsub'
    options = {
        'select': (str, None, "additional resources of each chunk"),
        'filesystems': (str, None, "filesystems needed by the job"),
        'place': (str, 'scatter', "placement of the chunks"),
        'tasks_per_node': (int, None, "array tasks per node"),
        'max_nodes': (int, None, "maximum number of nodes of array jobs"),
        'cores_per_node': (int, 128, "cores of a node"),
        'memory_per_node': (str, '256G', "memory of a node"),
        'launcher': (str, DEFAULT_LAUNCHER, "launcher template"),
        'export_env': (bool, True, "export the environment of qsub (-V)"),
    }

    def tasks_per_node(self, job):
        """Number of array tasks run on each node for `job`."""
        if self.settings['tasks_per_node'] is not None:
            return max(int(self.settings['tasks_per_node']), 1)
        cpus = max(int(job.cpus_per_task or 1), 1)
        n_tasks = self.settings['cores_per_node'] // cpus
        if job.mem_per_cpu is not None:
            mem_per_task = memory_to_gb(job.mem_per_cpu) * cpus
            if mem_per_task > 0:
                n_tasks = min(n_tasks, int(memory_to_gb(
                    self.settings['memory_per_node']) // mem_per_task))
        return max(n_tasks, 1)

    def array_layout(self, job):
        """Number of launched ranks, ranks per node and nodes of an array
        job."""
        n_ranks = len(job.array.ids)
        if job.array.max_concurrent is not None:
            n_ranks = min(n_ranks, max(int(job.array.max_concurrent), 1))
        ppn = min(self.tasks_per_node(job), n_ranks)
        if self.settings['max_nodes'] is not None:
            n_ranks = min(n_ranks, ppn * max(self.settings['max_nodes'], 1))
        nodes = max(math.ceil(n_ranks / ppn), 1)
        return n_ranks, ppn, nodes

    def directives(self, job):
        """#PBS lines for `job`."""
        if job.array is not None:
            nodes = self.array_layout(job)[2]
        else:
            nodes = job.nodes or 1
        select = f"{int(nodes)}"
        if self.settings['select']:
            select += ':' + self.settings['select'].lstrip(':')
        d = [f"-N {pbs_job_name(job.name)}", f"-l select={select}"]
        if self.settings['place']:
            d.append(f"-l place={self.settings['place']}")
        if job.walltime is not None:
            d.append(f"-l walltime={walltime_to_hms(job.walltime)}")
        if self.settings['filesystems']:
            d.append(f"-l filesystems={self.settings['filesystems']}")
        if job.account is not None:
            d.append(f"-A {job.account}")
        if job.queue is not None:
            d.append(f"-q {job.queue}")
        if self.settings['export_env']:
            d.append("-V")
        d.append("-j oe")
        if job.output is not None:
            # PBS has no job ID pattern for output file names
            d.append("-o " + replace_placeholders(job.output, job_id='job'))
        if job.email is not None and job.mail_type is not None:
            d.append(f"-m {mail_type_to_pbs(job.mail_type)}")
            d.append(f"-M {job.email}")
        return ["#PBS " + directive for directive in d]

    def _fill(self, text):
        return replace_placeholders(text, job_id=PBS_JOB_NUMBER,
                                    job_end=f"${JOB_END_VARIABLE}")

    def body(self, job):
        """Shell lines of the PBS script."""
        # PBS starts jobs in the home directory
        lines = ['cd "${PBS_O_WORKDIR:-.}"']
        if any(has_placeholder(text, 'job_end')
               for text in job.all_commands()):
            if job.walltime is None:
                raise ValueError("Jobs using {job_end} need a walltime.")
            lines.append(f"export {JOB_END_VARIABLE}=$(( $(date +%s) + "
                         f"{walltime_to_seconds(job.walltime)} ))")
        lines += [self._fill(line) for line in job.setup]
        lines += [self._fill(line) for line in job.commands]
        if job.launched_command:
            lines.append(self._fill(job.launched_command))
        if job.array is not None:
            lines += self.array_lines(job)
        return lines

    def array_lines(self, job):
        """Lines launching the tasks of an array job."""
        array = job.array
        n_ranks, ppn, nodes = self.array_layout(job)
        depth = max(int(job.cpus_per_task or 1), 1)
        runner = ["python", "-m", "posydon.CLI.schedulers.task_runner",
                  "--nranks", str(n_ranks)]
        if array.is_range:
            runner += ["--ntasks", str(len(array.ids))]
        else:
            runner += ["--task-ids", ",".join(map(str, array.ids))]
        runner += ["--task-count", str(array.size)]
        if array.log_file is not None:
            runner += ["--log", array.log_file]
        runner_text = " ".join(shlex.quote(part) for part in runner + ["--"])
        # {task_id} and {job_id} are replaced by the task runner, $VARIABLES
        # by the shell
        command = replace_placeholders(array.command,
                                       job_end=f"${JOB_END_VARIABLE}")
        command_text = " ".join(
            part if re.fullmatch(r'[\w./:={}$%+,@-]+', part)
            else shlex.quote(part) for part in shlex.split(command))
        launcher = self.settings['launcher'].format(
            ntasks=n_ranks, ppn=ppn, depth=depth, nodes=nodes)
        return [f"export POSYDON_ARRAY_JOB_ID={PBS_JOB_NUMBER}",
                f"export OMP_NUM_THREADS=${{OMP_NUM_THREADS:-{depth}}}",
                f"{launcher} {runner_text} {command_text}"]

    def submit_line(self, variable, script, after=None, condition='ok'):
        """Bash line submitting `script` with qsub."""
        self._check_condition(condition)
        after = self._dependency_variables(after)
        dependency = ""
        if after:
            ids = ":".join("${" + var + "}" for var in after)
            dependency = f"-W depend=after{condition}:{ids} "
        return f"{variable}=$(qsub {dependency}{script})\n"

    def mpi_launcher(self, n_ranks, ranks_per_node=None, cpus_per_rank=1):
        """MPI launch inside PBS jobs, using the launcher template."""
        ppn = ranks_per_node or n_ranks
        return self.settings['launcher'].format(
            ntasks=n_ranks, ppn=ppn, depth=max(int(cpus_per_rank or 1), 1),
            nodes=max(math.ceil(n_ranks / ppn), 1))
