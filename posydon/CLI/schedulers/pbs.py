"""PBS Pro backend of POSYDON's job scheduler interface.

PBS sites like ALCF Crux schedule whole nodes (jobs never share a node) and
allow only a few jobs or array subjobs per project. Translating each array
task into its own PBS job would occupy a full node for every single-core
task. Array jobs are therefore "packed": a single PBS job requests enough
nodes to run `tasks_per_node` tasks per node and starts them with the MPI
launcher (`mpiexec`). Each launched process runs
`python -m posydon.CLI.schedulers.task_runner`, which runs the task command
for its share of the array indices, with one log file per task, like a
SLURM array task. If there are more tasks than ranks, a rank that becomes
free takes the next task nobody has taken yet (a work queue, see
`posydon.CLI.schedulers.task_queue`), so no rank sits idle while tasks are
left.

The ranks are spread evenly over the cores of a node (`--depth`), e.g. 56
single-core tasks get 2 cores each, using both sockets of a node.
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
DEFAULT_LAUNCHER = ("{mpiexec} -n {ntasks} --ppn {ppn} --depth {depth} "
                    "--cpu-bind depth")

# Where job scripts look for the system's mpiexec (Cray PALS on ALCF
# systems). Python environments often bring their own MPI (e.g. Open MPI
# with mpi4py from conda), whose mpiexec would otherwise come first in PATH.
DEFAULT_MPIEXEC_SEARCH = '/opt/cray/pals/*/bin/mpiexec'

# Shell variable holding the mpiexec found in the job script
MPIEXEC_VARIABLE = 'POSYDON_MPIEXEC'

# Shell expression giving the numeric part of the PBS job ID
# (e.g. "12345.crux-pbs-01" -> "12345")
PBS_JOB_NUMBER = '${PBS_JOBID%%.*}'

# Environment variable holding the end of the job (seconds since the epoch)
JOB_END_VARIABLE = 'POSYDON_JOB_END'

# Ways of assigning array tasks to the ranks of a packed job
TASK_ORDERS = ('dynamic', 'static')


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
        Launcher template with the fields {mpiexec}, {ntasks}, {ppn},
        {depth} and {nodes} (default: Cray PALS mpiexec). {mpiexec} is the
        mpiexec found with `mpiexec_search`. {depth} is the number of cores
        per rank: the cores of a node divided evenly between its ranks, but
        at least the cores per task.
    mpiexec_search : str
        File pattern the job script uses to find the system's mpiexec when
        it starts; the latest version found is used, otherwise the mpiexec
        in PATH (default: "/opt/cray/pals/*/bin/mpiexec"; empty: always the
        mpiexec in PATH).
    export_env : bool
        Export the environment of qsub to the job ("-V"), as SLURM does by
        default (default: True).
    rerunnable : bool
        Write "-r y" (True) or "-r n" (False); not set by default. A
        rerunnable job is restarted from the beginning after e.g. being
        preempted.
    task_stats : str
        Directory in which the task runner records the start, end, CPU time
        and peak memory of each array task (see
        `posydon.CLI.schedulers.task_stats`); off by default.
    task_order : str
        "dynamic" (default): if there are more tasks than ranks, a free rank
        takes the next task nobody has taken yet. "static": rank r runs the
        tasks r, r+nranks, ... assigned at the start.
    queue_dir : str
        Directory of the dynamic task queues, on a shared file system
        (default: ".task_queue" in the submission directory).
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
        'mpiexec_search': (str, DEFAULT_MPIEXEC_SEARCH,
                           "file pattern of the system's mpiexec"),
        'export_env': (bool, True, "export the environment of qsub (-V)"),
        'rerunnable': (bool, None, "rerun the job after e.g. preemption"),
        'task_stats': (str, None, "directory for per-task statistics"),
        'task_order': (str, 'dynamic', "'dynamic' or 'static' assignment "
                       "of array tasks to ranks"),
        'queue_dir': (str, '.task_queue', "directory of the task queues"),
    }

    def __init__(self, **options):
        super().__init__(**options)
        self.settings['task_order'] = str(self.settings['task_order']).lower()
        if self.settings['task_order'] not in TASK_ORDERS:
            raise ValueError(f"task_order must be one of {TASK_ORDERS}, not "
                             f"'{self.settings['task_order']}'.")

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
        if self.settings['rerunnable'] is not None:
            d.append("-r y" if self.settings['rerunnable'] else "-r n")
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
        if any(f"${{{MPIEXEC_VARIABLE}" in line for line in lines):
            lines.insert(1, self.mpiexec_search_line())
        return lines

    def mpiexec(self):
        """Value of the {mpiexec} launcher field."""
        if self.settings['mpiexec_search']:
            return f'"${{{MPIEXEC_VARIABLE}:-mpiexec}}"'
        return "mpiexec"

    def mpiexec_search_line(self):
        """Shell line finding the latest system mpiexec."""
        return (f"{MPIEXEC_VARIABLE}=$(ls {self.settings['mpiexec_search']} "
                "2>/dev/null | sort -V | tail -n 1)")

    def depth(self, job, ppn):
        """Cores per rank: the cores of a node divided evenly between `ppn`
        ranks, but at least the cores per task."""
        cpus = max(int(job.cpus_per_task or 1), 1)
        return max(cpus, self.settings['cores_per_node'] // max(ppn, 1))

    def array_lines(self, job):
        """Lines launching the tasks of an array job."""
        array = job.array
        n_ranks, ppn, nodes = self.array_layout(job)
        cpus = max(int(job.cpus_per_task or 1), 1)
        depth = self.depth(job, ppn)
        dynamic = (self.settings['task_order'] == 'dynamic'
                   and n_ranks < len(array.ids))
        runner = ["python", "-m", "posydon.CLI.schedulers.task_runner",
                  "--nranks", str(n_ranks)]
        if array.is_range:
            runner += ["--ntasks", str(len(array.ids))]
        else:
            runner += ["--task-ids", ",".join(map(str, array.ids))]
        runner += ["--task-count", str(array.size)]
        if array.log_file is not None:
            runner += ["--log", array.log_file]
        if self.settings['task_stats']:
            runner += ["--stats", self.settings['task_stats']]
        if dynamic:
            runner += ["--queue", self.settings['queue_dir']]
        runner_text = " ".join(shlex.quote(part) for part in runner + ["--"])
        # {task_id} and {job_id} are replaced by the task runner, $VARIABLES
        # by the shell
        command = replace_placeholders(array.command,
                                       job_end=f"${JOB_END_VARIABLE}")
        command_text = " ".join(
            part if re.fullmatch(r'[\w./:={}$%+,@-]+', part)
            else shlex.quote(part) for part in shlex.split(command))
        launcher = self.settings['launcher'].format(
            mpiexec=self.mpiexec(), ntasks=n_ranks, ppn=ppn, depth=depth,
            nodes=nodes)
        lines = [f"export POSYDON_ARRAY_JOB_ID={PBS_JOB_NUMBER}"]
        if dynamic:
            # unique for each run of the job (a rerun job keeps its ID)
            lines.append("export POSYDON_QUEUE_ID="
                         f"{PBS_JOB_NUMBER}_$(date +%s)")
        lines += [f"export OMP_NUM_THREADS=${{OMP_NUM_THREADS:-{cpus}}}",
                  f"{launcher} {runner_text} {command_text}"]
        return lines

    def submit_line(self, variable, script, after=None, condition='ok'):
        """Bash line submitting `script` with qsub."""
        self._check_condition(condition)
        after = self._dependency_variables(after)
        dependency = ""
        if after:
            ids = ":".join("${" + var + "}" for var in after)
            dependency = f"-W depend=after{condition}:{ids} "
        return (f"{variable}=$(qsub {dependency}{script})"
                + self._on_failure(script) + "\n")

    def mpi_launcher(self, n_ranks, ranks_per_node=None, cpus_per_rank=1):
        """MPI launch inside PBS jobs, using the launcher template."""
        ppn = ranks_per_node or n_ranks
        return self.settings['launcher'].format(
            mpiexec=self.mpiexec(), ntasks=n_ranks, ppn=ppn,
            depth=max(int(cpus_per_rank or 1), 1),
            nodes=max(math.ceil(n_ranks / ppn), 1))
