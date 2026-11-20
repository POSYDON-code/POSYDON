"""Read MESA inlists into dictionaries and write them back.

An inlist is plain text made of one or more *sections* (Fortran calls them
namelists). Each section starts with ``&name``, lists ``key = value`` lines
and ends with ``/``. Everything after a ``!`` is a comment::

    &controls
        initial_mass = 10d0      ! solar masses
        x_ctrl(1) = 0.5d0
    / ! end of controls section

    &star_job
        pgstar_flag = .false.
    /

``parse_inlist`` turns this into::

    {"controls": {"initial_mass": "10d0", "x_ctrl(1)": "0.5d0"},
     "star_job": {"pgstar_flag": ".false."}}

Values stay strings exactly as written (``10d0`` is not turned into a float),
so writing them back gives MESA the same text. MESA's ``*.defaults`` files
have no ``&name`` lines; their keys go under the ``default_section`` given by
the caller.

Compared with ``posydon.utils.gridutils.clean_inlist_file``, which
``posydon-setup-grid`` uses today:

- a ``!`` inside quotes is kept (``clean_inlist_file`` cuts the line there),
- section names are returned without the leading ``&``,
- values to skip are chosen by the caller: pass ``LEGACY_SKIPPED_VALUES`` to
  drop ``''`` and ``'.'`` like ``clean_inlist_file`` does.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

from pathlib import Path

# values that posydon.utils.gridutils.clean_inlist_file silently drops
LEGACY_SKIPPED_VALUES = ("''", "'.'")


def strip_comment(line):
    """Remove a ``!`` comment from one line, ignoring ``!`` inside quotes.

    Parameters
    ----------
    line : str
        One line of an inlist.

    Returns
    -------
    str
        The line without its comment, stripped of surrounding whitespace.

    Examples
    --------
    >>> strip_comment("initial_mass = 10d0  ! in Msun")
    'initial_mass = 10d0'
    >>> strip_comment("log_directory = 'LOGS!1'")
    "log_directory = 'LOGS!1'"
    """
    quote = None
    for i, char in enumerate(line):
        if quote is not None:
            if char == quote:
                quote = None
        elif char in ("'", '"'):
            quote = char
        elif char == "!":
            return line[:i].strip()
    return line.strip()


def parse_inlist(text, default_section=None, skip_values=()):
    """Parse inlist text into ``{section: {key: value}}``.

    Parameters
    ----------
    text : str
        The content of an inlist or of a MESA ``*.defaults`` file.
    default_section : str or None
        Section for ``key = value`` lines that appear outside any
        ``&name`` ... ``/`` block, as in MESA's ``*.defaults`` files.
    skip_values : iterable of str
        Values to leave out, e.g. ``LEGACY_SKIPPED_VALUES``.

    Returns
    -------
    dict
        Section name (without ``&``) to an ordered dict of key to value.
        Keys keep the order of their first appearance; a key that appears
        twice keeps its last value.
    """
    skip_values = set(skip_values)
    sections = {}
    current = None
    for raw_line in text.splitlines():
        line = strip_comment(raw_line)
        if not line:
            continue
        if line.startswith("&"):
            # a new section starts: "&controls"
            current = line[1:].split()[0]
            sections.setdefault(current, {})
        elif line.startswith("/"):
            # the section ends; following keys belong to no section
            current = None
        elif "=" in line:
            key, value = (part.strip() for part in line.split("=", 1))
            if value in skip_values:
                continue
            section = current if current is not None else default_section
            sections.setdefault(section, {})[key] = value
    return sections


def read_inlist(path, default_section=None, skip_values=()):
    """Read an inlist file; see ``parse_inlist`` for the parameters."""
    text = Path(path).read_text()
    return parse_inlist(text, default_section=default_section,
                        skip_values=skip_values)


def format_value(value):
    """Turn a Python value into the text MESA expects in an inlist.

    Parameters
    ----------
    value : bool, int, float or str
        The value to write. Strings are written as they are, so a string
        value must already carry its quotes, e.g. ``"'LOGS1'"``.

    Returns
    -------
    str
        ``.true.``/``.false.`` for booleans, a plain integer for integers and
        for floats without a fractional part, Fortran double precision
        (``1.2345000000d+01``) for other floats.

    Examples
    --------
    >>> format_value(True)
    '.true.'
    >>> format_value(3.0)
    '3'
    >>> format_value(12.345)
    '1.2345000000d+01'
    """
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if value.is_integer():
            return str(int(value))
        return "{:.10e}".format(value).replace("e", "d")
    return str(value)


def render_inlist_section(name, values, indent="    "):
    """Write one inlist section as text.

    Parameters
    ----------
    name : str
        Section name without ``&``, e.g. ``"controls"``.
    values : dict
        Key to value; values go through ``format_value``.
    indent : str
        Text written before each ``key = value`` line.

    Returns
    -------
    str
        ``&name``, one line per key, and a closing ``/``.

    Examples
    --------
    >>> print(render_inlist_section("star_job", {"pgstar_flag": False}),
    ...       end="")
    &star_job
        pgstar_flag = .false.
    /
    """
    lines = ["&{}".format(name)]
    lines += ["{}{} = {}".format(indent, key, format_value(value))
              for key, value in values.items()]
    lines.append("/")
    return "\n".join(lines) + "\n"


def render_inlist(sections, indent="    "):
    """Write a whole inlist: its sections separated by a blank line.

    Parameters
    ----------
    sections : dict
        Section name to its values, as returned by ``parse_inlist``.
    indent : str
        Passed on to ``render_inlist_section``.

    Returns
    -------
    str
        The inlist text. ``parse_inlist`` reads it back to ``sections``.
    """
    return "\n".join(render_inlist_section(name, values, indent=indent)
                     for name, values in sections.items())
