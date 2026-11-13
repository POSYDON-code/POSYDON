"""SLURM backend of POSYDON's job scheduler interface."""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>",
]

import re

from posydon.CLI.schedulers.base import Scheduler
from posydon.CLI.schedulers.utils import replace_placeholders


def _slurm_log_pattern(index_letter):
    """Placeholder replacement for SLURM file name patterns (%a, %4a)."""
    def replace(spec):
        match = re.fullmatch(r':0?(\d+)d', spec)
        if match:
            return f"%{match.group(1)}{index_letter}"
        return f"%{index_letter}"
    return replace


class SlurmScheduler(Scheduler):
    """Job scripts for SLURM (`sbatch`).

    Array jobs use SLURM job arrays; each array task is started with the
    task launcher (default: `srun`) unless `ArraySpec.launched` is False.

    Options
    -------
    exclude : str
        Comma separated list of nodes to exclude.
    launcher : str
        Command prefix starting a task (default: "srun").
    """

    name = 'slurm'
    extension = 'slurm'
    submit_program = 'sbatch'
    options = {
        'exclude': (str, None, "comma separated list of nodes to exclude"),
        'launcher': (str, 'srun', "command prefix starting a task"),
    }

    def directives(self, job):
        """#SBATCH lines for `job`."""
        d = []
        d.append(f"--job-name={job.name}")
        if job.account is not None:
            d.append(f"--account={job.account}")
        if job.queue is not None:
            d.append(f"--partition={job.queue}")
        if job.nodes is not None:
            d.append(f"--nodes={job.nodes}")
        if job.tasks_per_node is not None:
            d.append(f"--ntasks-per-node={job.tasks_per_node}")
        if job.cpus_per_task is not None:
            d.append(f"--cpus-per-task={job.cpus_per_task}")
        if job.mem_per_cpu is not None:
            d.append(f"--mem-per-cpu={job.mem_per_cpu}")
        if job.walltime is not None:
            d.append(f"--time={job.walltime}")
        output = job.output
        if job.array is not None:
            array = job.array
            if array.is_range:
                indices = f"0-{len(array.ids) - 1}"
            else:
                indices = ",".join(map(str, array.ids))
            if array.max_concurrent is not None:
                indices += f"%{array.max_concurrent}"
            d.append(f"--array={indices}")
            if array.log_file is not None:
                output = replace_placeholders(
                    array.log_file, job_id=_slurm_log_pattern('A'),
                    task_id=_slurm_log_pattern('a'))
        elif output is not None:
            output = replace_placeholders(output,
                                          job_id=_slurm_log_pattern('j'))
        if job.truncate_output:
            d.append("--open-mode=truncate")
        if output is not None:
            d.append(f"--output={output}")
        if job.email is not None and job.mail_type is not None:
            d.append(f"--mail-type={job.mail_type}")
            d.append(f"--mail-user={job.email}")
        if self.settings['exclude'] is not None:
            d.append(f"--exclude={self.settings['exclude']}")
        return ["#SBATCH " + directive for directive in d]

    def _fill(self, text, array):
        job_id = '$SLURM_ARRAY_JOB_ID' if array else '$SLURM_JOB_ID'
        return replace_placeholders(text, task_id='$SLURM_ARRAY_TASK_ID',
                                    job_id=job_id,
                                    job_end='$SLURM_JOB_END_TIME')

    def body(self, job):
        """Shell lines of the SLURM script."""
        is_array = job.array is not None
        lines = [self._fill(line, is_array) for line in job.setup]
        if is_array and (not job.array.is_range
                         or job.array.size != len(job.array.ids)):
            # tasks of a larger array: let programs know the full array
            lines.append(f"export SLURM_ARRAY_TASK_COUNT={job.array.size}")
            lines.append("export SLURM_ARRAY_TASK_MIN=0")
        lines += [self._fill(line, is_array) for line in job.commands]
        launcher = self.settings['launcher']
        if job.launched_command:
            lines.append(f"{launcher} "
                         + self._fill(job.launched_command, is_array))
        if is_array:
            command = self._fill(job.array.command, True)
            if job.array.launched:
                command = f"{launcher} {command}"
            lines.append(command)
        return lines

    def submit_line(self, variable, script, after=None, condition='ok'):
        """Bash line submitting `script` with sbatch."""
        self._check_condition(condition)
        after = self._dependency_variables(after)
        dependency = ""
        if after:
            ids = ":".join("${" + var + "}" for var in after)
            dependency = (f"--dependency=after{condition}:{ids} "
                          "--kill-on-invalid-dep=yes ")
        return (f"{variable}=$(sbatch --parsable {dependency}{script})"
                + self._on_failure(script) + "\n")

    def mpi_launcher(self, n_ranks, ranks_per_node=None, cpus_per_rank=1):
        """MPI launch used inside SLURM jobs."""
        return f"mpirun --bind-to none -np {n_ranks}"
