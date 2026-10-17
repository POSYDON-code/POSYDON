"""Record what a run directory was set up from (tracker C13).

``posydon-grid setup`` writes ``setup_record.ini`` into the run directory::

    [setup]        the ini, its format, the grid type and mode, when, on
                   which host and with which command
    [posydon]      the POSYDON that ran: version, location, its commit and
                   whether its code differed from it; Python
    [mesa]         $MESA_DIR and $MESASDK_ROOT
    [scenario]     the release (or legacy scenario), and the inlist
                   repository as it was read: folder, origin, commit,
                   whether the files used differed from it (inlist_dirty),
                   and a fingerprint of those files
    [inlists]      per written inlist and namelist, the files stacked into
                   it, lowest layer first: mesa, release (posydon for a
                   legacy ini), user
    [whole_files]  per column list and Fortran source in the run
                   directory, the file it was copied from

A grid can then be traced to its physics even after the ini, POSYDON or
the inlist repository changes; two grids with the same fingerprint read
identical files. ``None`` means git could not tell (not a git repository).
Paths inside the inlist repository are relative to it.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import configparser
import datetime
import io
import platform
import socket
import sys
from pathlib import Path

from posydon.grids.setup.resolve import git_state

RECORD_FILE = "setup_record.ini"
# how the layers are named in [inlists]; with a release, the
# posydon_inlists layer holds only its files
LAYER_NAMES = {"mesa_inlists": "mesa", "posydon_inlists": "posydon",
               "user_inlists": "user", "generated": "generated"}


def posydon_state():
    """The POSYDON that runs: version, location, commit, dirty, Python.

    ``dirty`` looks at tracked files of ``posydon/`` and ``bin/`` only: the
    code that ran, not notes or notebooks next to it.
    """
    import posydon
    root = Path(posydon.__file__).resolve().parents[1]
    commit, dirty, _ = git_state(root, [root / "posydon", root / "bin"],
                                 untracked=False)
    return {"version": str(posydon.__version__), "location": str(root),
            "commit": str(commit), "dirty": str(dirty),
            "python": "{} ({})".format(platform.python_version(),
                                       sys.executable)}


def _shown(path, repo):
    """A path relative to the inlist repository, if it is inside it."""
    if repo is not None:
        try:
            return Path(path).resolve().relative_to(
                Path(repo).resolve()).as_posix()
        except (ValueError, TypeError):
            pass
    return str(path)


def record_text(config, resolved, inifile, mesa_dir, posydon_version,
                mesasdk_root=None, command=None, stacks=None,
                whole_files=None, posydon=None, now=None, host=None):
    """The text of ``setup_record.ini``.

    Parameters
    ----------
    config : GridConfig
    resolved : ResolvedScenario or None
        From ``resolve.resolve_scenario``.
    inifile, mesa_dir, mesasdk_root : str or Path
    posydon_version : str
    command : str or None
        The command line that ran the setup; None when it was called from
        Python.
    stacks : dict or None
        ``BuiltInlists.stacks``: written namelist to its ``(layer,
        source)`` files.
    whole_files : dict or None
        File in the run directory to the file it was copied from.
    posydon : dict or None
        ``posydon_state()``; read now if None.
    now, host : str or None
        When and where; read now if None.

    Returns
    -------
    str
    """
    repo = resolved.repo if resolved is not None else None
    parser = configparser.ConfigParser(interpolation=None)
    parser.optionxform = str
    parser["setup"] = {
        "inifile": str(Path(inifile).resolve()),
        "ini_format": config.source_format,
        "grid_type": str(config.grid.type),
        "mode": config.grid.mode,
        "written": now or datetime.datetime.now().isoformat(
            timespec="seconds"),
        "host": host or socket.gethostname(),
        "command": command or "None (called from Python)",
    }
    parser["posydon"] = dict(posydon or posydon_state(),
                             version=str(posydon_version))
    parser["mesa"] = {"mesa_dir": str(mesa_dir),
                      "mesasdk_root": str(mesasdk_root)}
    if resolved is not None:
        scenario = config.scenario
        values = {}
        if resolved.release is None:
            # legacy: the scenario the repository was checked out at
            values["source"] = scenario.source
        for key in ("release", "rerun", "branch", "commit", "system_type"):
            if getattr(scenario, key) is not None:
                values[key] = getattr(scenario, key)
        if resolved.release is not None:
            values["mesa_version"] = resolved.release.mesa_version
            values["mesa"] = resolved.release.mesa
        # what was read: the folder as it was, and its files
        values["inlist_repo"] = str(resolved.repo)
        values["inlist_remote"] = str(resolved.remote)
        values["inlist_commit"] = str(resolved.commit)
        values["inlist_dirty"] = str(resolved.dirty)
        values["inlist_fingerprint"] = str(resolved.fingerprint)
        values["inlist_files"] = str(resolved.n_files)
        parser["scenario"] = values
    if stacks:
        names = dict(LAYER_NAMES)
        if resolved is not None and resolved.release is not None:
            names["posydon_inlists"] = "release"
        parser["inlists"] = {
            written: "".join("\n{:9s}{}".format(
                names.get(layer, layer), _shown(source, repo))
                for layer, source in sources)
            for written, sources in stacks.items()}
    if whole_files:
        parser["whole_files"] = {
            target: _shown(source, repo)
            for target, source in sorted(whole_files.items())}
    out = io.StringIO()
    out.write("; What this run directory was set up from; written by "
              "posydon-grid setup.\n\n")
    parser.write(out)
    return out.getvalue().rstrip("\n") + "\n"


def write_record(run_directory, *args, **options):
    """Write ``setup_record.ini`` into the run directory; return its path.

    The other arguments are those of ``record_text``.
    """
    path = Path(run_directory) / RECORD_FILE
    path.write_text(record_text(*args, **options))
    return path
