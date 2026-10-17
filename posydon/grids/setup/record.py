"""Record what a run directory was set up from (tracker C13).

``posydon-grid setup`` writes ``setup_record.ini`` into the run directory:
the ini it read, the POSYDON version, ``$MESA_DIR``, the grid type and, for
a scenario, the inlist repository as it was read: its folder, ``origin``,
commit, whether the files used differed from that commit
(``inlist_dirty``) and a fingerprint of those files. A grid can then be
traced to its physics even after the ini or the inlist repository
changes; two grids with the same fingerprint read identical files.
``None`` means git could not tell (not a git repository).
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import configparser
import io
from pathlib import Path

RECORD_FILE = "setup_record.ini"


def record_text(config, resolved, inifile, mesa_dir, posydon_version):
    """The text of ``setup_record.ini``.

    Parameters
    ----------
    config : GridConfig
    resolved : ResolvedScenario or None
        From ``resolve.resolve_scenario``.
    inifile, mesa_dir : str or Path
    posydon_version : str

    Returns
    -------
    str
    """
    parser = configparser.ConfigParser(interpolation=None)
    parser.optionxform = str
    parser["setup"] = {
        "inifile": str(Path(inifile).resolve()),
        "ini_format": config.source_format,
        "posydon_version": str(posydon_version),
        "mesa_dir": str(mesa_dir),
        "grid_type": str(config.grid.type),
        "mode": config.grid.mode,
    }
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
    out = io.StringIO()
    out.write("; What this run directory was set up from; written by "
              "posydon-grid setup.\n\n")
    parser.write(out)
    return out.getvalue().rstrip("\n") + "\n"


def write_record(run_directory, *args):
    """Write ``setup_record.ini`` into the run directory; return its path.

    The arguments after ``run_directory`` are those of ``record_text``.
    """
    path = Path(run_directory) / RECORD_FILE
    path.write_text(record_text(*args))
    return path
