"""Helpers shared by the job scheduler backends."""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>",
]

import re

# Placeholders that tools may use in commands and log file names. Each
# scheduler backend replaces them with its own mechanism, e.g. SLURM uses
# $SLURM_ARRAY_TASK_ID in commands and %a in log file names.
PLACEHOLDER_NAMES = ('task_id', 'job_id', 'job_end')
PLACEHOLDER_RE = re.compile(r'\{(' + '|'.join(PLACEHOLDER_NAMES)
                            + r')(:[^{}]*)?\}')


def replace_placeholders(text, **replacements):
    """Replace "{task_id}", "{job_id}" and "{job_end}" placeholders.

    Only the given placeholders are replaced; everything else, including
    shell constructs like "${VAR}", is left untouched.

    Parameters
    ----------
    text : str
        The text containing placeholders, e.g. "log_{task_id:04d}.out".
    **replacements
        For each placeholder name either a value, which is formatted with the
        placeholder's format spec (e.g. "{task_id:04d}" -> "0007"), or a
        callable taking the format spec (e.g. ":04d" or "") and returning
        the replacement string.

    Returns
    -------
    str
        The text with the placeholders replaced.
    """
    def substitute(match):
        name, spec = match.group(1), match.group(2) or ''
        if name not in replacements:
            return match.group(0)
        value = replacements[name]
        if callable(value):
            return value(spec)
        return format(value, spec[1:])
    return PLACEHOLDER_RE.sub(substitute, text)


def has_placeholder(text, name):
    """Whether `text` contains the placeholder `name` (with any spec)."""
    return any(match.group(1) == name
               for match in PLACEHOLDER_RE.finditer(text or ''))


def walltime_to_seconds(walltime):
    """Convert a walltime to seconds.

    Accepted are the formats of SLURM's `--time` option: "MM", "MM:SS",
    "HH:MM:SS", "D-HH", "D-HH:MM" and "D-HH:MM:SS", where hours may exceed
    24 (PBS style "HH:MM:SS"). Surrounding quotes are ignored.

    Parameters
    ----------
    walltime : str or int
        The walltime. A number is interpreted as minutes (as in SLURM).

    Returns
    -------
    int
        The walltime in seconds.
    """
    if isinstance(walltime, (int, float)) and not isinstance(walltime, bool):
        return int(walltime * 60)
    if not isinstance(walltime, str):
        raise TypeError(f"walltime must be a string, got {type(walltime)}")
    text = walltime.strip().strip('"').strip("'")
    match = re.fullmatch(r'(?:(\d+)-)?(\d+)(?::(\d+))?(?::(\d+))?', text)
    if match is None:
        raise ValueError(f"Cannot interpret walltime '{walltime}'")
    days, first, second, third = match.groups()
    first = int(first)
    if days is not None:
        # D-HH, D-HH:MM, D-HH:MM:SS
        return ((int(days) * 24 + first) * 60 + int(second or 0)) * 60 \
            + int(third or 0)
    if second is None:      # MM
        return first * 60
    if third is None:       # MM:SS
        return first * 60 + int(second)
    return (first * 60 + int(second)) * 60 + int(third)  # HH:MM:SS


def walltime_to_hms(walltime):
    """Convert a walltime to "HH:MM:SS" (hours may exceed 24)."""
    total = walltime_to_seconds(walltime)
    hours, rest = divmod(total, 3600)
    minutes, seconds = divmod(rest, 60)
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"


def memory_to_gb(memory):
    """Convert a memory specification like "4G", "4000M" or "1T" to GB.

    Numbers without unit are interpreted as MB (SLURM's default).
    """
    if isinstance(memory, (int, float)) and not isinstance(memory, bool):
        return float(memory) / 1024.0
    if not isinstance(memory, str):
        raise TypeError(f"memory must be a string, got {type(memory)}")
    match = re.fullmatch(r'\s*([0-9]*\.?[0-9]+)\s*([KMGT]?)B?\s*',
                         memory.strip('"').strip("'"), flags=re.IGNORECASE)
    if match is None:
        raise ValueError(f"Cannot interpret memory '{memory}'")
    value, unit = float(match.group(1)), match.group(2).upper()
    factor = {'K': 1.0 / 1024**2, '': 1.0 / 1024, 'M': 1.0 / 1024,
              'G': 1.0, 'T': 1024.0}[unit]
    return value * factor


def to_bool(value):
    """Interpret a configuration value as boolean."""
    if isinstance(value, str):
        return value.strip().lower() in ('1', 'true', 'yes', 'y', 'on')
    return bool(value)
