"""Scheduler independent description of batch jobs.

POSYDON's command line tools (`posydon-popsyn`, `posydon-setup-grid` and
`posydon-setup-pipeline`) describe each job they need with a `JobSpec`. A
`Scheduler` backend (e.g. SLURM) turns this description into a job script
and provides the bash lines to submit it, including dependencies between
jobs. Adding support for another scheduler therefore only requires a new
`Scheduler` subclass registered in `posydon.CLI.schedulers`.

Commands and log file names may contain the placeholders "{task_id}"
(index of an array task), "{job_id}" (job ID) and "{job_end}" (the time the
job will be stopped, in seconds since the epoch). Format specs are allowed,
e.g. "log_{task_id:04d}.out".
"""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>",
]

import os
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import List, Optional, Sequence

from posydon.CLI.schedulers.utils import to_bool
from posydon.utils.posydonwarning import Pwarn


@dataclass
class ArraySpec:
    """An array of independent tasks within a job.

    Parameters
    ----------
    command : str
        Command run for each task; "{task_id}" is replaced by the task index.
    n_tasks : int, optional
        Number of tasks with the indices 0..n_tasks-1.
    task_ids : sequence of int, optional
        Explicit task indices (e.g. to rerun failed tasks); overrides
        `n_tasks`.
    task_count : int, optional
        Size of the full array the tasks belong to (default: number of
        tasks). Programs see it as the array size, so that rerun tasks split
        the work in the same way as the original array.
    max_concurrent : int, optional
        Maximum number of tasks running at the same time.
    log_file : str, optional
        Log file per task, using "{job_id}" and "{task_id}".
    launched : bool
        Start each task with the scheduler's task launcher (e.g. `srun`).
    """

    command: str
    n_tasks: Optional[int] = None
    task_ids: Optional[Sequence[int]] = None
    task_count: Optional[int] = None
    max_concurrent: Optional[int] = None
    log_file: Optional[str] = None
    launched: bool = True

    def __post_init__(self):
        if self.task_ids is not None:
            self.task_ids = sorted(int(i) for i in self.task_ids)
        elif self.n_tasks is None:
            raise ValueError("Either n_tasks or task_ids is needed.")
        if len(self.ids) < 1:
            raise ValueError("A job array needs at least one task.")

    @property
    def ids(self):
        """The list of task indices."""
        if self.task_ids is not None:
            return list(self.task_ids)
        return list(range(int(self.n_tasks)))

    @property
    def is_range(self):
        """True if the tasks are 0..N-1."""
        return self.task_ids is None

    @property
    def size(self):
        """Size of the full array (see `task_count`)."""
        if self.task_count is not None:
            return int(self.task_count)
        if self.task_ids is not None:
            return max(self.task_ids) + 1
        return int(self.n_tasks)


@dataclass
class JobSpec:
    """Scheduler independent description of a batch job.

    Parameters
    ----------
    name : str
        Job name.
    walltime : str, optional
        Time limit, e.g. "23:00:00" or "2-00:00:00".
    nodes : int, optional
        Number of nodes; for array jobs the nodes per task.
    tasks_per_node : int, optional
        Number of (MPI) tasks per node.
    cpus_per_task : int, optional
        Number of cores per task (or array task).
    mem_per_cpu : str, optional
        Memory per core, e.g. "4G".
    account : str, optional
        Account/project to charge.
    queue : str, optional
        Queue (SLURM: partition).
    email : str, optional
        Address for notifications.
    mail_type : str, optional
        When to send mails, in SLURM's vocabulary ("FAIL", "ALL", "BEGIN",
        "END", "NONE" or comma separated combinations); requires `email`.
    output : str, optional
        Log file of the job (for array jobs see `ArraySpec.log_file`).
    truncate_output : bool
        Overwrite an existing log file instead of appending to it.
    setup : list of str
        Shell lines run first, e.g. exports.
    commands : list of str
        Shell lines run as they are.
    launched_command : str, optional
        Main command, started with the scheduler's task launcher.
    array : ArraySpec, optional
        Array of tasks run by this job.
    """

    name: str
    walltime: Optional[str] = None
    nodes: Optional[int] = None
    tasks_per_node: Optional[int] = None
    cpus_per_task: Optional[int] = None
    mem_per_cpu: Optional[str] = None
    account: Optional[str] = None
    queue: Optional[str] = None
    email: Optional[str] = None
    mail_type: Optional[str] = None
    output: Optional[str] = None
    truncate_output: bool = False
    setup: List[str] = field(default_factory=list)
    commands: List[str] = field(default_factory=list)
    launched_command: Optional[str] = None
    array: Optional[ArraySpec] = None

    def all_commands(self):
        """All command texts of the job (for placeholder checks)."""
        texts = list(self.setup) + list(self.commands)
        if self.launched_command:
            texts.append(self.launched_command)
        if self.array is not None:
            texts.append(self.array.command)
        return texts


class Scheduler(ABC):
    """Base class of job scheduler backends.

    Subclasses define `name`, `extension`, `options` and the abstract
    methods. Options are settings of the computing site rather than of a
    single job (e.g. nodes to exclude); they are given when creating the
    scheduler and validated against `options`, a dictionary mapping option
    names to (type, default, help).
    """

    name = None
    extension = None
    submit_program = None
    options = {}

    def __init__(self, **options):
        self.settings = {key: value[1] for key, value in self.options.items()}
        for key, value in options.items():
            key = key.lower()
            if key not in self.options:
                Pwarn(f"Option '{key}' is not supported by the {self.name} "
                      "scheduler and will be ignored.",
                      "InappropriateValueWarning")
                continue
            if value is None:
                continue
            option_type = self.options[key][0]
            if option_type is bool:
                value = to_bool(value)
            elif not isinstance(value, option_type):
                value = option_type(value)
            self.settings[key] = value

    def given_options(self):
        """The options that differ from their defaults."""
        return {key: value for key, value in self.settings.items()
                if value != self.options[key][1]}

    def __repr__(self):
        return f"{type(self).__name__}({self.given_options()})"

    def file_name(self, stem):
        """File name of a job script, e.g. "cleanup" -> "cleanup.slurm"."""
        return f"{stem}.{self.extension}"

    def render(self, job):
        """The full job script for `job` (ending with a new line)."""
        lines = ["#!/bin/bash"] + self.directives(job) + [""] + self.body(job)
        return "\n".join(lines) + "\n"

    def write(self, job, path):
        """Write the job script for `job` to `path`.

        Returns
        -------
        str
            `path`
        """
        if os.path.exists(path):
            Pwarn(f'Replace {path}', "OverwriteWarning")
        with open(path, 'w') as f:
            f.write(self.render(job))
        return path

    @abstractmethod
    def directives(self, job):
        """Scheduler directives (header lines) for `job`."""

    @abstractmethod
    def body(self, job):
        """Shell lines of the job script after the directives."""

    @abstractmethod
    def submit_line(self, variable, script, after=None, condition='ok'):
        """Bash line submitting `script` and storing its job ID.

        Parameters
        ----------
        variable : str
            Name of the bash variable to store the job ID in.
        script : str
            Path to the job script.
        after : str or list of str, optional
            Names of bash variables holding job IDs this job depends on.
        condition : str
            'ok' (start only if all dependencies succeeded) or 'any' (start
            after the dependencies ended in any way).
        """

    def submit_command(self, script):
        """Command (list of arguments) submitting `script` directly."""
        return [self.submit_program, script]

    @abstractmethod
    def mpi_launcher(self, n_ranks, ranks_per_node=None, cpus_per_rank=1):
        """Command prefix starting an MPI program with `n_ranks` ranks."""

    @staticmethod
    def _dependency_variables(after):
        if after is None:
            return []
        if isinstance(after, str):
            return [after]
        return list(after)

    @staticmethod
    def _check_condition(condition):
        if condition not in ('ok', 'any'):
            raise ValueError(f"Unknown dependency condition '{condition}', "
                             "use 'ok' or 'any'.")
