"""Dynamic assignment of array tasks to the ranks of a packed job.

In a packed PBS job (see `posydon.CLI.schedulers.pbs`), many array tasks
share a fixed number of ranks. Assigning the tasks round-robin when the job
starts leaves ranks idle when task run times differ a lot (e.g. MESA runs
taking minutes or days). With a `TaskQueue`, a rank that becomes free takes
the next task that no other rank has taken yet, so all ranks stay busy until
the tasks run out or the job ends.

The ranks coordinate through a directory on a shared file system: it holds a
single counter file `next_<n>`, where n is the position (in the list of
tasks) of the next task to hand out. A rank claims position n by renaming
the file to `next_<n+1>`. A rename succeeds for exactly one process, also
across the nodes of a shared file system like Lustre, so every task is run
by one rank only. No locks are used, so a rank that dies cannot block the
others.

The directory must be unique for each run of a job (e.g. the job ID and the
start time), so that a rerun job does not continue the counter of an
earlier run.
"""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>",
]

import os
import random
import time

COUNTER_PREFIX = 'next_'


class TaskQueue:
    """Hand out the tasks of a packed array job to its ranks.

    Rank r first runs the task at position r. Afterwards, ranks claim the
    remaining positions (from `nranks` on) one at a time from the shared
    counter. Rank 0 creates the counter when it starts.

    Parameters
    ----------
    directory : str
        Directory of the counter, on a file system shared by all ranks.
    task_ids : list of int
        The array indices, in the order in which they are handed out.
    rank : int
        Rank of this process.
    nranks : int
        Number of ranks.
    timeout : float
        Seconds to wait for the counter to appear before giving up.
    poll : float
        Seconds between checks for the counter.
    """

    def __init__(self, directory, task_ids, rank, nranks, timeout=600.0,
                 poll=0.5):
        self.directory = directory
        self.task_ids = list(task_ids)
        self.rank = int(rank)
        self.nranks = int(nranks)
        self.timeout = timeout
        self.poll = poll

    def _path(self, position):
        return os.path.join(self.directory, f"{COUNTER_PREFIX}{position}")

    def _counter(self):
        """Current position of the counter, or None if there is none."""
        try:
            names = os.listdir(self.directory)
        except FileNotFoundError:
            return None
        positions = [int(name[len(COUNTER_PREFIX):]) for name in names
                     if name.startswith(COUNTER_PREFIX)
                     and name[len(COUNTER_PREFIX):].isdigit()]
        return max(positions) if positions else None

    def create_counter(self):
        """Create the counter at position `nranks` (done by rank 0)."""
        os.makedirs(self.directory, exist_ok=True)
        if self._counter() is None:
            with open(self._path(self.nranks), 'x'):
                pass

    def claim(self):
        """Claim the next task.

        Returns
        -------
        int or None
            The array index of the claimed task, or None if all tasks have
            been handed out.

        Raises
        ------
        TimeoutError
            If the counter does not appear within `timeout` seconds.
        """
        deadline = time.monotonic() + self.timeout
        while True:
            position = self._counter()
            if position is None:
                if time.monotonic() > deadline:
                    raise TimeoutError(
                        f"no task counter in {self.directory} after "
                        f"{self.timeout:g} s")
                time.sleep(self.poll)
                continue
            if position >= len(self.task_ids):
                return None
            try:
                os.rename(self._path(position), self._path(position + 1))
            except FileNotFoundError:
                # another rank claimed this position first
                time.sleep(random.uniform(0.0, 0.05))
                continue
            return self.task_ids[position]

    def tasks(self):
        """Iterate over the tasks of this rank."""
        if self.rank < len(self.task_ids):
            if self.rank == 0 and len(self.task_ids) > self.nranks:
                self.create_counter()
            yield self.task_ids[self.rank]
        if len(self.task_ids) <= self.nranks:
            return
        while True:
            task_id = self.claim()
            if task_id is None:
                return
            yield task_id
