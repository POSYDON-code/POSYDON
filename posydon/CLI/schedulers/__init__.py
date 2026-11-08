"""Job scheduler backends of POSYDON's command line tools.

Tools describe their jobs with `JobSpec` and use a `Scheduler` to write job
scripts and submission lines. Get a scheduler with `get_scheduler`::

    from posydon.CLI.schedulers import JobSpec, ArraySpec, get_scheduler

    scheduler = get_scheduler('slurm', exclude='node01')
    job = JobSpec(name='example', walltime='1:00:00',
                  array=ArraySpec(command='python run.py {task_id}',
                                  n_tasks=10))
    scheduler.write(job, scheduler.file_name('example'))

To support another scheduler, subclass `Scheduler` and add it to
`SCHEDULERS`.
"""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>",
]

from posydon.CLI.schedulers.base import ArraySpec, JobSpec, Scheduler
from posydon.CLI.schedulers.slurm import SlurmScheduler

# registry of the available scheduler backends
SCHEDULERS = {
    SlurmScheduler.name: SlurmScheduler,
}
DEFAULT_SCHEDULER = 'slurm'


def scheduler_names():
    """Names of the available schedulers."""
    return tuple(SCHEDULERS)


def normalize_scheduler_name(name, default=DEFAULT_SCHEDULER):
    """Return a valid scheduler name or `default`.

    Parameters
    ----------
    name : any
        Requested scheduler; anything that is not the name of an available
        scheduler (e.g. None) results in `default`.
    default : str or None
        Value returned for unknown names.
    """
    if isinstance(name, str) and name.strip().lower() in SCHEDULERS:
        return name.strip().lower()
    return default


def get_scheduler(name=DEFAULT_SCHEDULER, **options):
    """Create the scheduler backend `name` with the given options.

    Parameters
    ----------
    name : str
        Name of the scheduler, e.g. 'slurm'.
    **options
        Options of the scheduler (see its `options` attribute).

    Returns
    -------
    Scheduler
    """
    key = normalize_scheduler_name(name, default=None)
    if key is None:
        raise ValueError(f"Unknown scheduler '{name}'. Available: "
                         f"{', '.join(scheduler_names())}")
    return SCHEDULERS[key](**options)


def scheduler_options_from_mapping(name, mapping):
    """Pick the options of scheduler `name` from a configuration mapping.

    Keys are matched case-insensitively, so this works for ini sections like
    `[slurm]` (lower case keys) and the pipeline's `[account]` section
    (upper case keys). Other keys are ignored.

    Parameters
    ----------
    name : str
        Name of the scheduler.
    mapping : dict
        Configuration values.

    Returns
    -------
    dict
        The scheduler options found in `mapping`.
    """
    option_names = SCHEDULERS[normalize_scheduler_name(name)].options
    return {key.lower(): value for key, value in mapping.items()
            if key.lower() in option_names}


def parse_option_strings(option_strings):
    """Parse ["key=value", ...] (e.g. from the command line) into a dict.

    The value is everything after the first "=", so "select=system=crux"
    gives {'select': 'system=crux'}.
    """
    options = {}
    for text in option_strings or []:
        key, sep, value = str(text).partition('=')
        if not sep or not key.strip():
            raise ValueError(f"Scheduler options need the form KEY=VALUE, "
                             f"got '{text}'")
        options[key.strip().lower()] = value.strip()
    return options


__all__ = ['ArraySpec', 'JobSpec', 'Scheduler', 'SlurmScheduler',
           'SCHEDULERS', 'DEFAULT_SCHEDULER', 'get_scheduler',
           'normalize_scheduler_name', 'parse_option_strings',
           'scheduler_names', 'scheduler_options_from_mapping']
