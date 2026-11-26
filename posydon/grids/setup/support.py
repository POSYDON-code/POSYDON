"""Files that used to ship in POSYDON's ``grid_params/``, from the inlist repo.

POSYDON no longer ships ``grid_params/`` (MESA's defaults, the ini
templates and the POSYDON-MESA-INLISTS submodule). Its files live in
POSYDON-MESA-INLISTS, and POSYDON reads them at one pinned commit,
``SUPPORT_COMMIT``, named here and changed only by editing this file.

``support_root`` extracts the files needed at that commit into a cache
folder once (``git archive`` from the local clone, fetching if the clone
does not have the commit yet) and returns it; the cache is never changed
afterwards. Legacy ini files name the old paths through
``${posydon_github_root}/grid_params/...``; ``map_grid_params_path`` says
where each one lives in the cache::

    grid_params/defaults/r11701/star|binary/*  r11701/mesa/default/defaults/
    grid_params/defaults/r11701/*.list         r11701/mesa/default/columns/
    grid_params/defaults/r11701/src/*          r11701/mesa/default/src/
    grid_params/*.ini                          templates/
    grid_params/POSYDON-MESA-INLISTS[-public]/r11701/
        ZAMS_models/, running_scripts/         the same
        default_common_inlists/                r11701/releases/v2.0/common/
        CO-H_star/, CO-He_star/                .../v2.0/CO-HMS_RLO/, CO-HeMS/
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from posydon.grids.setup.config import ConfigError
from posydon.grids.setup.scenario import (
    REPOSITORIES,
    clone_if_missing,
    git,
)

# POSYDON-MESA-INLISTS, branch v2_sources: the restructured layout with
# r11701/mesa/default (byte-identical to grid_params/defaults/r11701)
SUPPORT_COMMIT = "8482d7d2d6f4abd3374437053aa6cf2b11ff70ca"
# what is extracted from it
SUPPORT_PATHS = ("r11701/mesa", "r11701/ZAMS_models",
                 "r11701/releases/v2.0", "r11701/running_scripts",
                 "templates")
OLD_FOLDERS = {"default_common_inlists": "common",
               "CO-H_star": "CO-HMS_RLO",
               "CO-He_star": "CO-HeMS",
               "HeMS-HMS": "HeMS-HMS"}
GRID_PARAMS = re.compile(r"^(?P<root>.*?)/?grid_params/(?P<rest>.+)$")


def map_grid_params_path(path):
    """Where a ``grid_params/...`` path lives in the support cache.

    Parameters
    ----------
    path : str or Path
        A path containing ``grid_params/``.

    Returns
    -------
    str or None
        The path relative to the cache, or None if ``path`` is not under
        ``grid_params/`` or has no counterpart.

    Examples
    --------
    >>> map_grid_params_path("/x/grid_params/defaults/r11701/star/"
    ...                      "controls.defaults")
    'r11701/mesa/default/defaults/star/controls.defaults'
    >>> map_grid_params_path("/x/grid_params/POSYDON-MESA-INLISTS/r11701/"
    ...                      "ZAMS_models/zams_z1.42m2_y0.2703.data")
    'r11701/ZAMS_models/zams_z1.42m2_y0.2703.data'
    """
    match = GRID_PARAMS.match(str(path))
    if match is None:
        return None
    rest = match.group("rest")
    parts = rest.split("/")
    if parts[:2] == ["defaults", "r11701"] and len(parts) > 2:
        name = "/".join(parts[2:])
        if parts[2] in ("star", "binary"):
            return "r11701/mesa/default/defaults/" + name
        if parts[2] == "src":
            return "r11701/mesa/default/src/" + "/".join(parts[3:])
        if name.endswith(".list"):
            return "r11701/mesa/default/columns/" + name
        return None
    if len(parts) == 1 and rest.endswith(".ini"):
        return "templates/" + rest
    if parts[0] in ("POSYDON-MESA-INLISTS", "POSYDON-MESA-INLISTS-public") \
            and len(parts) > 3 and parts[1] == "r11701":
        folder, inner = parts[2], "/".join(parts[3:])
        if folder in ("ZAMS_models", "running_scripts"):
            return "r11701/{}/{}".format(folder, inner)
        if folder in OLD_FOLDERS:
            return "r11701/releases/v2.0/{}/{}".format(OLD_FOLDERS[folder],
                                                     inner)
    return None


def support_root(home=None, clone=None, commit=SUPPORT_COMMIT):
    """The cache folder with the support files at ``commit``.

    Parameters
    ----------
    home : str, Path or None
        Where the default clone and the cache live; defaults to ``~``.
    clone : str, Path or None
        A clone of POSYDON-MESA-INLISTS; defaults to
        ``<home>/.posydon_mesa_inlists`` (cloned if missing).
    commit : str
        The commit to extract.

    Returns
    -------
    Path
        ``<home>/.posydon_mesa_inlists_support/<commit>``.

    Raises
    ------
    ConfigError
        If the clone cannot be made or does not have the commit.
    """
    home = Path.home() if home is None else Path(home)
    cache = home / ".posydon_mesa_inlists_support" / commit
    if cache.is_dir():
        return cache
    url, name = REPOSITORIES["posydon"]
    clone = home / name if clone is None else Path(clone)
    clone_if_missing(url, clone, "mesa_inlists", "posydon_github_root",
                     "grid_params/ is no longer part of POSYDON, and ")
    if git(clone, "cat-file", "-e", commit + "^{commit}").returncode:
        git(clone, "fetch", "origin")
    archive = git(clone, "archive", "--format=tar", commit, *SUPPORT_PATHS,
                  text=False)
    if archive.returncode != 0:
        raise ConfigError("mesa_inlists", "posydon_github_root",
                          "grid_params/ is no longer part of POSYDON, and "
                          "{} does not have commit {} of "
                          "POSYDON-MESA-INLISTS: {}".format(
                              clone, commit,
                              archive.stderr.decode().strip()))
    cache.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(dir=cache.parent))
    subprocess.run(["tar", "-x", "-C", str(staging)], input=archive.stdout,
                   check=True)
    try:
        staging.rename(cache)
    except OSError:  # another setup extracted it first
        shutil.rmtree(staging)
    return cache
