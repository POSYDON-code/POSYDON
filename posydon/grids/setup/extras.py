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

Several keys can put a file with the same name in the same place (e.g.
``mesa_binary_extras`` and ``posydon_binary_extras`` both give
``binary/src/run_binary_extras.f``); then only the highest layer's is kept,
``mesa`` < ``posydon`` < ``user`` (``keep_highest_extras``; fixed C4:
before, the alphabetically last prefix of all ``*binary_extras*`` keys
won, and every other ``*binary_extras*`` key was dropped, whatever file it
gave). Places are compared by the name MESA compiles (``file_name``), and
files are copied under that name.

``run_binary_extras`` and ``run_star_extras`` in an inlist layer are whole
files too; ``whole_file_extras`` turns them into extras keys of that layer.

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

# extras keys by layer, lowest first: mesa_*, posydon_*, user_*
LAYER_PREFIXES = ("mesa", "posydon", "user")
COLUMN_LIST_FILES = {"history_columns": "history_columns.list",
                     "binary_history_columns": "binary_history_columns.list",
                     "profile_columns": "profile_columns.list"}
# the file MESA compiles for each kind of source (first match wins): a
# source given under another name is copied under this one
SOURCE_NAMES = (("star_binary_extras", "run_star_extras.f"),
                ("binary_extras", "run_binary_extras.f"),
                ("star1_extras", "run_star_extras.f"),
                ("star2_extras", "run_star_extras.f"),
                ("binary_run", "binary_run.f"),
                ("star_run", "run.f"))
# whole files of the inlist layers that are Fortran sources, and the extras
# keys (after the layer prefix) they stand for
WHOLE_FILE_EXTRAS = {
    "run_binary_extras": ("binary_extras",),
    "run_star_extras": ("star_binary_extras", "star1_extras",
                        "star2_extras"),
}
INLIST_LAYERS = ("mesa_inlists", "posydon_inlists", "user_inlists")
# the folders make_executables recreates, in today's order
SOURCE_FOLDERS = ("star1/src", "star1/make", "star2/src", "star2/make",
                  "binary/src", "binary/make")


def merge_extras(config_extras, scenario_extras=None):
    """The extras of the ini with the scenario's applied.

    The scenario's keys replace ini keys of the same name in place and are
    appended otherwise (today it writes into the same dictionary). If the
    scenario gives star 1 its ``run_star_extras.f`` (``mesa_star1_extras``)
    but names none for star 2, star 2 gets the same file (fixed C1: before,
    star 2 kept MESA's stock one from the ini).
    """
    scenario_extras = scenario_extras or {}
    extras = {key: Path(path) for key, path in config_extras.items()}
    for key, path in scenario_extras.items():
        extras[key] = Path(path)
    if "mesa_star1_extras" in scenario_extras and \
            "mesa_star2_extras" not in scenario_extras:
        extras["mesa_star2_extras"] = Path(
            scenario_extras["mesa_star1_extras"])
    return extras


def file_name(key, path):
    """The name a key's file gets in the run directory.

    Sources are copied under the name MESA compiles (``SOURCE_NAMES``), so
    ``user_star1_extras = my_extras.f`` becomes ``star1/src/
    run_star_extras.f`` (before, it was copied as ``my_extras.f``, which
    MESA never compiles, next to the release's ``run_star_extras.f``).
    Other files keep their own name.

    Examples
    --------
    >>> file_name("user_star1_extras", "my_extras.f")
    'run_star_extras.f'
    >>> file_name("posydon_star_binary_extras", "x/a.f")
    'run_star_extras.f'
    >>> file_name("something_else", "x/notes.txt")
    'notes.txt'
    """
    for part, name in SOURCE_NAMES:
        if part in key:
            return name
    return Path(path).name


def whole_file_extras(config):
    """The Fortran sources the inlist layers name as whole files.

    ``run_binary_extras`` and ``run_star_extras`` in ``[mesa_inlists]``,
    ``[posydon_inlists]`` or ``[user_inlists]`` become extras keys of that
    layer: ``run_star_extras`` in ``[user_inlists]`` is
    ``user_star_binary_extras``, ``user_star1_extras`` and
    ``user_star2_extras``, so it is used wherever a star runs. Before, these
    keys were accepted and then never used.

    Returns
    -------
    dict
        Extras key to path, in layer order.
    """
    extras = {}
    for layer, prefix in zip(INLIST_LAYERS, LAYER_PREFIXES):
        whole = config.inlists[layer].whole_files
        for name, parts in WHOLE_FILE_EXTRAS.items():
            if whole.get(name) is not None:
                for part in parts:
                    extras["{}_{}".format(prefix, part)] = Path(whole[name])
    return extras


def keep_highest_extras(extras):
    """For each file the extras put in place, keep only the highest layer's.

    Keys are grouped by where their file goes and its name; within a group,
    the key whose prefix (the part before the first ``_``) is highest in
    ``LAYER_PREFIXES`` is kept. Keys without such a prefix (``binary_run``,
    ``star_run``, makefiles) are kept as they are.

    Examples
    --------
    >>> keep_highest_extras({"mesa_binary_extras": "x/run_binary_extras.f",
    ...                      "mesa_star_binary_extras": "x/run_star_extras.f",
    ...                      "user_binary_extras": "y/run_binary_extras.f",
    ...                      "binary_run": "x/binary_run.f"})
    ... # doctest: +NORMALIZE_WHITESPACE
    {'mesa_star_binary_extras': 'x/run_star_extras.f',
     'user_binary_extras': 'y/run_binary_extras.f',
     'binary_run': 'x/binary_run.f'}
    """
    best = {}
    for key, path in extras.items():
        prefix = key.split("_")[0]
        if prefix not in LAYER_PREFIXES or path is None:
            continue
        place = (tuple(destinations(key)), file_name(key, path))
        rank = LAYER_PREFIXES.index(prefix)
        if place not in best or rank >= best[place][0]:
            best[place] = (rank, key)
    kept = {key for _, key in best.values()}
    return {key: path for key, path in extras.items()
            if key in kept or key.split("_")[0] not in LAYER_PREFIXES
            or path is None}


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
                target = target / file_name(key, source)
            shutil.copy(source, target)
    mk = run / "mk"
    # a #! line, so ./mk runs on its own (C15: before, only through sh)
    mk.write_text("#!/bin/sh\n" + "".join(
        line + "\n" for line in mk_lines(extras, run_directory)))
    mk.chmod(0o755)
    return mk


def compile_executables(run_directory):
    """Run ``mk`` in the run directory; return the three executables.

    The make output goes to the terminal, as today. ``mk`` is run through
    ``sh``, as ``os.system`` did, which works whether or not the run
    directory's file system allows executing it.
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


def replaced_column_lists(config, scenario):
    """Column lists the ini names but the scenario replaces with another file.

    Within ``posydon_inlists`` the scenario's list wins (tracker C7); this
    says which of the ini's lists that ignores, so the user can be told.
    Lists whose files have the same content are not reported.

    Returns
    -------
    list of str
        Column-list names, in ``COLUMN_LIST_FILES`` order.
    """
    if scenario is None:
        return []
    ini = config.inlists["posydon_inlists"].whole_files
    used = scenario.inlists.whole_files
    names = []
    for name in COLUMN_LIST_FILES:
        if name in ini and name in used:
            mine, theirs = Path(ini[name]), Path(used[name])
            if not (mine.is_file() and
                    mine.read_bytes() == theirs.read_bytes()):
                names.append(name)
    return names


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
