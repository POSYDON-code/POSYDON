"""Stage the Fortran sources and column lists, and compile.

These files are not merged like inlists: each one is copied whole.

Fortran sources and makefiles come from ``[mesa_extras]`` (legacy ini) or
``[extras]``, plus what the scenario sets (``ScenarioFiles.extras``). Each
key says where its file goes, by today's rules (``destinations``)::

    key contains            copied to
    binary_extras           binary/src/
    binary_run              binary/src/
    star_run                star1/src/ and star2/src/
    star1_extras            star1/src/
    star2_extras            star2/src/
    makefile_binary         binary/make/<key>   + "make -f <key>" in mk
    makefile_star           star1/make/<key>, star2/make/<key>  + mk lines
    mesa_dir                (not a file)
    anything else           the run directory itself

Several ``*binary_extras*`` keys can name the same file name; only the keys
with the highest prefix are kept (``mesa`` < ``posydon`` < ``user``, by
alphabetical order as today; tracker C4).

The column lists (``history_columns``, ``binary_history_columns``,
``profile_columns``) go to ``column_lists/``; the highest layer that names
one wins (``layers.pick_whole_file``).

``mk`` and the copies are written into the run directory. Today ``mk`` is
written to the current working directory instead (tracker S9).
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import os
import shutil
import subprocess
from pathlib import Path

from posydon.grids.setup.layers import pick_whole_file

COLUMN_LIST_FILES = {"history_columns": "history_columns.list",
                     "binary_history_columns": "binary_history_columns.list",
                     "profile_columns": "profile_columns.list"}
# the folders make_executables recreates, in today's order
SOURCE_FOLDERS = ("star1/src", "star1/make", "star2/src", "star2/make",
                  "binary/src", "binary/make")


def merge_extras(config_extras, scenario_extras=None):
    """The extras of the ini with the scenario's applied, as today.

    The scenario's keys replace ini keys of the same name in place and are
    appended otherwise (today it writes into the same dictionary).
    """
    extras = {key: Path(path) for key, path in config_extras.items()}
    for key, path in (scenario_extras or {}).items():
        extras[key] = Path(path)
    return extras


def keep_highest_binary_extras(extras):
    """Drop ``*binary_extras*`` keys whose prefix is not the highest.

    The prefix is the part before the first ``_`` (``mesa``, ``posydon``,
    ``user``); the alphabetically last one is kept, as today.

    Examples
    --------
    >>> keep_highest_binary_extras({"mesa_binary_extras": "a",
    ...                             "posydon_binary_extras": "b",
    ...                             "binary_run": "c"})
    {'posydon_binary_extras': 'b', 'binary_run': 'c'}
    """
    prefixes = sorted({key.split("_")[0] for key in extras
                       if "binary_extras" in key})
    if not prefixes:
        return dict(extras)
    return {key: path for key, path in extras.items()
            if "binary_extras" not in key or prefixes[-1] in key}


def destinations(key):
    """Where a ``[mesa_extras]`` key's file goes, relative to the run.

    Returns
    -------
    list of str
        Folders (ending in ``/``) or file paths; empty for ``mesa_dir``.

    Examples
    --------
    >>> destinations("star_run")
    ['star1/src/', 'star2/src/']
    >>> destinations("makefile_star")
    ['star1/make/makefile_star', 'star2/make/makefile_star']
    """
    if "binary_extras" in key or "binary_run" in key:
        return ["binary/src/"]
    if "star_run" in key:
        return ["star1/src/", "star2/src/"]
    if "star1_extras" in key:
        return ["star1/src/"]
    if "star2_extras" in key:
        return ["star2/src/"]
    if "makefile_binary" in key:
        return ["binary/make/" + key]
    if "makefile_star" in key:
        return ["star1/make/" + key, "star2/make/" + key]
    if key == "mesa_dir":
        return []
    return [""]


def mk_lines(extras, run_directory):
    """The commands of ``mk``: one ``cd`` + ``make`` per makefile."""
    run = str(run_directory)
    lines = []
    for key in extras:
        for destination in destinations(key):
            if "/make/" in destination:
                folder = os.path.join(run, os.path.dirname(destination))
                lines += ["cd {}".format(folder),
                          "make -f {}".format(key)]
    return lines


def stage_extras(extras, run_directory):
    """Copy the sources and makefiles into the run directory, write ``mk``.

    The ``src`` and ``make`` folders are emptied first, as today.

    Returns
    -------
    Path
        The ``mk`` script (executable).
    """
    run = Path(run_directory)
    for folder in SOURCE_FOLDERS:
        path = run / folder
        if path.exists():
            shutil.rmtree(path)
        path.mkdir(parents=True)
    for key, source in extras.items():
        if source is None:
            continue
        for destination in destinations(key):
            target = run / destination
            if destination.endswith("/") or destination == "":
                target = target / Path(source).name
            shutil.copy(source, target)
    mk = run / "mk"
    mk.write_text("".join(line + "\n"
                          for line in mk_lines(extras, run_directory)))
    mk.chmod(0o755)
    return mk


def compile_executables(run_directory):
    """Run ``mk`` in the run directory; return the three executables.

    The make output goes to the terminal, as today. ``mk`` has no ``#!``
    line (today's layout), so it is run through ``sh``, as ``os.system``
    did.
    """
    run = Path(run_directory)
    subprocess.run(["sh", "./mk"], cwd=run, check=False)
    return (str(run / "binary" / "binary"), str(run / "star1" / "star"),
            str(run / "star2" / "star"))


def column_list_sources(config, scenario=None):
    """The file to use for each column list, by layer.

    Within ``posydon_inlists`` the scenario's file replaces the ini's, as
    today; a ``user_inlists`` file replaces both.

    Returns
    -------
    dict
        Column-list name to the chosen path, for the names that have one.
    """
    chosen = {}
    for name in COLUMN_LIST_FILES:
        candidates = {}
        for layer in ("mesa_inlists", "posydon_inlists", "user_inlists"):
            path = config.inlists[layer].whole_files.get(name)
            if layer == "posydon_inlists" and scenario is not None:
                path = scenario.inlists.whole_files.get(name, path)
            candidates[layer] = path
        picked = pick_whole_file(candidates)
        if picked is not None:
            chosen[name] = Path(picked[1])
    return chosen


def stage_column_lists(sources, run_directory):
    """Copy the column lists into ``column_lists/`` (emptied first).

    Returns
    -------
    dict
        Column-list name to its path in the run directory, as given to
        ``posydon-run-grid``.
    """
    folder = os.path.join(str(run_directory), "column_lists")
    if os.path.exists(folder):
        shutil.rmtree(folder)
    os.makedirs(folder)
    staged = {}
    for name, file_name in COLUMN_LIST_FILES.items():
        staged[name] = os.path.join(folder, file_name)
        shutil.copy(sources[name], staged[name])
    return staged
