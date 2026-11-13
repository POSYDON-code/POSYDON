"""Per-task statistics of packed array jobs.

The task runner (`posydon.CLI.schedulers.task_runner`) can record, for each
array task, when it started and ended, its exit code, CPU time and peak
memory. Each rank writes its own CSV file (`stats_<jobid>_rank<rank>.csv`)
in a statistics directory, so that many ranks never write to the same file.
A "start" line is written when a task begins and an "end" line when it
finishes, so tasks killed by the walltime or by preemption show up as
unfinished.

Summarize a statistics directory with::

    python -m posydon.CLI.schedulers.task_stats stats/
"""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>",
]

import argparse
import csv
import glob
import math
import os
import socket
import statistics
import sys
import time

from posydon.CLI.schedulers.utils import memory_to_gb

STATS_FIELDS = ('event', 'job_id', 'task_id', 'rank', 'host', 'time',
                'elapsed_s', 'exit_code', 'user_cpu_s', 'sys_cpu_s',
                'max_rss_kb')


class StatsWriter:
    """Append task statistics of one rank to its CSV file.

    Parameters
    ----------
    directory : str
        Statistics directory (created if needed).
    job_id : str
        The job ID.
    rank : int
        Rank of the task runner.
    """

    def __init__(self, directory, job_id, rank):
        os.makedirs(directory, exist_ok=True)
        self.job_id = str(job_id)
        self.rank = int(rank)
        self.host = socket.gethostname()
        self.path = os.path.join(directory,
                                 f"stats_{self.job_id}_rank{self.rank}.csv")

    def _write(self, row):
        new_file = not os.path.exists(self.path)
        with open(self.path, 'a', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=STATS_FIELDS)
            if new_file:
                writer.writeheader()
            writer.writerow(row)

    def _row(self, event, task_id, timestamp):
        row = {field: '' for field in STATS_FIELDS}
        row.update(event=event, job_id=self.job_id, task_id=int(task_id),
                   rank=self.rank, host=self.host,
                   time=f"{timestamp:.3f}")
        return row

    def start(self, task_id, timestamp=None):
        """Record the start of a task."""
        if timestamp is None:
            timestamp = time.time()
        self._write(self._row('start', task_id, timestamp))

    def end(self, task_id, elapsed, exit_code, usage=None, timestamp=None):
        """Record the end of a task.

        Parameters
        ----------
        task_id : int
            The array index.
        elapsed : float
            Elapsed time of the task in seconds.
        exit_code : int
            Exit code of the task.
        usage : resource.struct_rusage or None
            Resource usage of the task (from os.wait4).
        timestamp : float or None
            Time of the end (default: now).
        """
        if timestamp is None:
            timestamp = time.time()
        row = self._row('end', task_id, timestamp)
        row.update(elapsed_s=f"{elapsed:.3f}", exit_code=int(exit_code))
        if usage is not None:
            row.update(user_cpu_s=f"{usage.ru_utime:.3f}",
                       sys_cpu_s=f"{usage.ru_stime:.3f}",
                       max_rss_kb=int(usage.ru_maxrss))
        self._write(row)


def _float(value):
    return float(value) if value not in (None, '') else None


def read_stats(directory, job_id=None):
    """Read the statistics files in `directory`.

    Parameters
    ----------
    directory : str
        Statistics directory.
    job_id : str or None
        Only read the files of this job.

    Returns
    -------
    dict
        Maps (job_id, task_id) to a dict with the keys 'host', 'start',
        'end', 'elapsed_s', 'exit_code', 'user_cpu_s', 'sys_cpu_s' and
        'max_rss_kb' (None if unknown, e.g. for unfinished tasks). If a task
        was run several times (e.g. a rerun job), its last run is used.
    """
    pattern = f"stats_{job_id}_rank*.csv" if job_id else "stats_*_rank*.csv"
    tasks = {}
    for path in sorted(glob.glob(os.path.join(directory, pattern))):
        with open(path, newline='') as f:
            for row in csv.DictReader(f):
                key = (row['job_id'], int(row['task_id']))
                if row['event'] == 'start':
                    tasks[key] = {'host': row['host'],
                                  'start': _float(row['time']), 'end': None,
                                  'elapsed_s': None, 'exit_code': None,
                                  'user_cpu_s': None, 'sys_cpu_s': None,
                                  'max_rss_kb': None}
                else:  # 'end'
                    task = tasks.setdefault(key, {'host': row['host'],
                                                  'start': None})
                    task.update(end=_float(row['time']),
                                elapsed_s=_float(row['elapsed_s']),
                                exit_code=int(row['exit_code']),
                                user_cpu_s=_float(row['user_cpu_s']),
                                sys_cpu_s=_float(row['sys_cpu_s']),
                                max_rss_kb=_float(row['max_rss_kb']))
    return tasks


def summarize(tasks):
    """Summarize task statistics (see `read_stats`).

    Returns
    -------
    dict
        Counts ('tasks', 'finished', 'failed', 'unfinished') and, for the
        finished tasks, 'elapsed_s' and 'max_rss_kb', each a dict with
        'min', 'median', 'max' (or None if no task finished).
    """
    finished = [t for t in tasks.values() if t.get('end') is not None]
    summary = {'tasks': len(tasks), 'finished': len(finished),
               'failed': sum(1 for t in finished if t['exit_code'] != 0),
               'unfinished': len(tasks) - len(finished)}
    for key in ('elapsed_s', 'max_rss_kb'):
        values = [t[key] for t in finished if t.get(key) is not None]
        summary[key] = ({'min': min(values),
                         'median': statistics.median(values),
                         'max': max(values)} if values else None)
    return summary


def suggest_tasks_per_node(max_rss_kb, memory_per_node='256G',
                           cores_per_node=128, cpus_per_task=1, safety=1.3):
    """Tasks per node fitting the largest measured peak memory.

    Parameters
    ----------
    max_rss_kb : float
        Largest peak memory of a task in KB.
    memory_per_node : str
        Memory of a node, e.g. "256G".
    cores_per_node : int
        Cores of a node.
    cpus_per_task : int
        Cores used by each task.
    safety : float
        Factor applied to the peak memory.

    Returns
    -------
    int
        Number of tasks per node (at least 1).
    """
    by_cores = cores_per_node // max(int(cpus_per_task), 1)
    per_task_gb = max_rss_kb / 1024.0**2 * safety
    if per_task_gb <= 0:
        return max(by_cores, 1)
    by_memory = math.floor(memory_to_gb(memory_per_node) / per_task_gb)
    return max(min(by_cores, by_memory), 1)


def _hms(seconds):
    seconds = int(round(seconds))
    return f"{seconds // 3600}:{seconds % 3600 // 60:02d}:{seconds % 60:02d}"


def format_summary(summary, memory_per_node='256G', cores_per_node=128,
                   cpus_per_task=1, safety=1.3):
    """Human readable summary (see `summarize`)."""
    lines = [f"tasks: {summary['tasks']}  finished: {summary['finished']}  "
             f"failed: {summary['failed']}  "
             f"unfinished: {summary['unfinished']}"]
    elapsed = summary['elapsed_s']
    if elapsed is not None:
        lines.append(f"elapsed (h:mm:ss): min {_hms(elapsed['min'])}  "
                     f"median {_hms(elapsed['median'])}  "
                     f"max {_hms(elapsed['max'])}")
    rss = summary['max_rss_kb']
    if rss is not None:
        to_gb = 1.0 / 1024.0**2
        lines.append(f"peak memory (GB): min {rss['min'] * to_gb:.2f}  "
                     f"median {rss['median'] * to_gb:.2f}  "
                     f"max {rss['max'] * to_gb:.2f}")
        lines.append("suggested tasks_per_node: " + str(suggest_tasks_per_node(
            rss['max'], memory_per_node, cores_per_node, cpus_per_task,
            safety)) + f" (node: {cores_per_node} cores, {memory_per_node}; "
            f"{cpus_per_task} core(s) per task; peak x {safety})")
    return "\n".join(lines)


def main(argv=None):
    """Print the summary of a statistics directory."""
    parser = argparse.ArgumentParser(
        prog='python -m posydon.CLI.schedulers.task_stats',
        description='Summarize the per-task statistics of packed array jobs.')
    parser.add_argument('directory', help='statistics directory')
    parser.add_argument('--job', default=None, help='only this job ID')
    parser.add_argument('--memory-per-node', default='256G')
    parser.add_argument('--cores-per-node', type=int, default=128)
    parser.add_argument('--cpus-per-task', type=int, default=1)
    parser.add_argument('--safety', type=float, default=1.3,
                        help='factor applied to the peak memory')
    args = parser.parse_args(argv)
    tasks = read_stats(args.directory, job_id=args.job)
    if not tasks:
        print(f"No task statistics found in {args.directory}",
              file=sys.stderr)
        return 1
    print(format_summary(summarize(tasks), args.memory_per_node,
                         args.cores_per_node, args.cpus_per_task,
                         args.safety))
    return 0


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main())
