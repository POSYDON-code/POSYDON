"""Record what a run directory was set up from (tracker C13).

``posydon-setup-grid`` writes ``setup_record.ini`` into the run directory:
the ini it read, the POSYDON version, ``$MESA_DIR``, the grid type and, for
a scenario, the inlist repository and the exact commit it was checked out
at. A grid can then be traced to its physics even after the ini or a
branch of the inlist repository changes.
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
        values = {"source": scenario.source}
        for key in ("release", "rerun", "ref", "branch", "commit",
                    "system_type"):
            if getattr(scenario, key) is not None:
                values[key] = getattr(scenario, key)
        values["inlist_repo"] = str(resolved.repo)
        values["inlist_commit"] = str(resolved.commit)
        if resolved.release is not None:
            values["mesa_version"] = resolved.release.mesa_version
            values["mesa"] = resolved.release.mesa
        parser["scenario"] = values
    out = io.StringIO()
    out.write("; What this run directory was set up from; written by "
              "posydon-setup-grid.\n\n")
    parser.write(out)
    return out.getvalue().rstrip("\n") + "\n"


def write_record(run_directory, *args):
    """Write ``setup_record.ini`` into the run directory; return its path.

    The arguments after ``run_directory`` are those of ``record_text``.
    """
    path = Path(run_directory) / RECORD_FILE
    path.write_text(record_text(*args))
    return path
