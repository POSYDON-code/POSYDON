"""Run array tasks inside a packed PBS job.

This module is started once per MPI rank by the job scripts that POSYDON
writes for PBS (see `posydon.CLI.schedulers.pbs`). It emulates a SLURM job array task:
each rank picks the array indices `rank, rank+nranks, rank+2*nranks, ...`
and, for each of them, runs the given command with the environment
variables

- POSYDON_ARRAY_JOB_ID: numeric part of the PBS job ID
- POSYDON_ARRAY_TASK_ID: the array index
- POSYDON_ARRAY_TASK_COUNT: number of tasks of the full array
- POSYDON_ARRAY_TASK_MIN: 0

and "{task_id}"/"{job_id}" in the command and log file name replaced.

Example::

    mpiexec -n 4 --ppn 4 python -m posydon.CLI.schedulers.task_runner \\
        --nranks 4 --ntasks 10 --log logs/task_{task_id}.out \\
        -- python run.py --index {task_id}
"""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>",
]

import argparse
import os
import subprocess
import sys
import time

from posydon.CLI.schedulers.task_stats import StatsWriter
from posydon.CLI.schedulers.utils import replace_placeholders

# environment variables in which MPI launchers store the rank of a process
RANK_ENV_VARS = ('PALS_RANKID', 'PMI_RANK', 'PMIX_RANK',
                 'OMPI_COMM_WORLD_RANK', 'MV2_COMM_WORLD_RANK',
                 'SLURM_PROCID')


def get_rank(environ=None):
    """Get the rank of this process from the launcher's environment.

    Parameters
    ----------
    environ : dict or None
        Environment to search (default: os.environ).

    Returns
    -------
    int or None
        The rank, or None if no rank variable is set.
    """
    if environ is None:
        environ = os.environ
    for var in RANK_ENV_VARS:
        value = environ.get(var)
        if value is not None and str(value).strip().isdigit():
            return int(value)
    return None


def get_job_id(environ=None):
    """Get the (numeric) job ID.

    Parameters
    ----------
    environ : dict or None
        Environment to search (default: os.environ).

    Returns
    -------
    str
        POSYDON_ARRAY_JOB_ID if set, otherwise the part of PBS_JOBID
        before the first '.', otherwise "0".
    """
    if environ is None:
        environ = os.environ
    job_id = environ.get('POSYDON_ARRAY_JOB_ID')
    if job_id:
        return str(job_id)
    pbs_job_id = environ.get('PBS_JOBID')
    if pbs_job_id:
        return pbs_job_id.split('.')[0].split('[')[0]
    return "0"


def get_my_tasks(rank, nranks, ntasks=None, task_ids=None):
    """Array indices to be processed by `rank`.

    Parameters
    ----------
    rank : int
        Rank of this process.
    nranks : int
        Total number of ranks.
    ntasks : int or None
        Number of tasks (indices 0..ntasks-1), used if `task_ids` is None.
    task_ids : list of int or None
        Explicit list of indices.

    Returns
    -------
    list of int
    """
    if task_ids is None:
        task_ids = list(range(int(ntasks)))
    return list(task_ids)[rank::nranks]


def parse_task_ids(text):
    """Parse "1,2,5-7" into [1, 2, 5, 6, 7]."""
    ids = []
    for part in str(text).split(','):
        part = part.strip()
        if not part:
            continue
        if '-' in part:
            start, end = part.split('-', 1)
            ids.extend(range(int(start), int(end) + 1))
        else:
            ids.append(int(part))
    return ids


def _execute(cmd, env, stdout=None):
    """Run `cmd` and return its exit code and resource usage.

    The resource usage (from os.wait4) covers the command and the child
    processes it waited for, e.g. MESA started by posydon-run-grid.
    """
    if stdout is None:
        process = subprocess.Popen(cmd, env=env)
    else:
        process = subprocess.Popen(cmd, env=env, stdout=stdout,
                                   stderr=subprocess.STDOUT)
    _, status, usage = os.wait4(process.pid, 0)
    process.returncode = os.waitstatus_to_exitcode(status)
    return process.returncode, usage


def run_task(command, task_id, job_id, task_count, log_file=None,
             return_usage=False):
    """Run the command of a single array task.

    Parameters
    ----------
    command : list of str
        The command; "{task_id}" and "{job_id}" are replaced.
    task_id : int
        The array index.
    job_id : str
        The job ID.
    task_count : int
        Total number of tasks of the array.
    log_file : str or None
        Log file pattern; None writes to the inherited stdout/stderr.
    return_usage : bool
        If True, also return the resource usage of the task.

    Returns
    -------
    int or tuple (int, resource.struct_rusage)
        The return code of the command (and its resource usage).
    """
    env = os.environ.copy()
    env['POSYDON_ARRAY_JOB_ID'] = str(job_id)
    env['POSYDON_ARRAY_TASK_ID'] = str(task_id)
    env['POSYDON_ARRAY_TASK_COUNT'] = str(task_count)
    env['POSYDON_ARRAY_TASK_MIN'] = '0'
    fields = {'task_id': int(task_id), 'job_id': job_id}
    cmd = [replace_placeholders(part, **fields) for part in command]
    if log_file is None:
        code, usage = _execute(cmd, env)
    else:
        log_path = replace_placeholders(log_file, **fields)
        log_dir = os.path.dirname(log_path)
        if log_dir:
            os.makedirs(log_dir, exist_ok=True)
        with open(log_path, 'w') as log:
            code, usage = _execute(cmd, env, stdout=log)
    if return_usage:
        return code, usage
    return code


def parse_args(argv=None):
    """Parse the command line of the task runner."""
    parser = argparse.ArgumentParser(
        prog='python -m posydon.CLI.schedulers.task_runner',
        description='Run POSYDON array tasks inside a packed PBS job.')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--ntasks', type=int,
                       help='number of array tasks (indices 0..ntasks-1)')
    group.add_argument('--task-ids', type=parse_task_ids,
                       help='comma-separated list of array indices')
    parser.add_argument('--nranks', type=int, required=True,
                        help='number of launched ranks')
    parser.add_argument('--rank', type=int, default=None,
                        help='rank of this process (default: from the '
                             'MPI launcher environment)')
    parser.add_argument('--task-count', type=int, default=None,
                        help='value of POSYDON_ARRAY_TASK_COUNT')
    parser.add_argument('--log', default=None,
                        help='log file pattern using {task_id}/{job_id}')
    parser.add_argument('--stats', default=None, metavar='DIR',
                        help='record start, end, CPU time and peak memory '
                             'of each task in DIR (see '
                             'posydon.CLI.schedulers.task_stats)')
    parser.add_argument('command', nargs=argparse.REMAINDER,
                        help='command to run (after "--")')
    args = parser.parse_args(argv)
    if args.command and args.command[0] == '--':
        args.command = args.command[1:]
    if not args.command:
        parser.error('no command given')
    return args


def main(argv=None):
    """Entry point; returns 0 on success and 1 if any task failed."""
    args = parse_args(argv)
    rank = args.rank if args.rank is not None else get_rank()
    if rank is None:
        if args.nranks > 1:
            print("task_runner: cannot determine the MPI rank; "
                  "use --rank or launch with mpiexec.", file=sys.stderr)
            return 1
        rank = 0
    if args.task_ids is not None:
        task_count = args.task_count or (max(args.task_ids) + 1
                                         if args.task_ids else 0)
    else:
        task_count = args.task_count or args.ntasks
    job_id = get_job_id()
    stats = None
    if args.stats is not None:
        stats = StatsWriter(args.stats, job_id, rank)
    failed = 0
    for task_id in get_my_tasks(rank, args.nranks, ntasks=args.ntasks,
                                task_ids=args.task_ids):
        if stats is not None:
            stats.start(task_id)
        start = time.time()
        code, usage = run_task(args.command, task_id, job_id, task_count,
                               log_file=args.log, return_usage=True)
        if stats is not None:
            stats.end(task_id, time.time() - start, code, usage)
        if code != 0:
            print(f"task_runner: task {task_id} failed with exit code "
                  f"{code}", file=sys.stderr)
            failed += 1
    return 1 if failed else 0


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main())
