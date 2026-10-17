"""Fill in what a new-format ini leaves to the release and MESA.

A legacy ini spells out its MESA defaults and Fortran sources. A
new-format ini may leave them out; they are then taken as every v2
production ini did:

- ``mesa_inlists``: MESA's ``.defaults`` of the release's MESA
  installation (``<mesa_version>/mesa/<installation>/defaults/``) for the
  binary, star 1 and star 2, and for every formation step of a star the
  grid type builds in steps. A single-star grid runs star 1 as its
  formation steps, so star 1's defaults go there.
- ``extras``: MESA's own work files in ``$MESA_DIR`` (makefiles,
  ``run.f``, ``binary_run.f`` and the stock extras), under the legacy
  ``[mesa_extras]`` names; the scenario's sources then replace the stock
  extras, as before.

Anything the ini sets is kept as it is.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

from pathlib import Path

from posydon.grids.setup.config import ConfigError
from posydon.grids.setup.grid_types import get_grid_type
from posydon.grids.setup.manifest import mesa_installation

# legacy [mesa_extras] name -> file in $MESA_DIR, in the order of the v2 inis
MESA_EXTRAS = {
    "makefile_binary": "binary/work/make/makefile",
    "makefile_star": "star/work/make/makefile",
    "mesa_binary_extras": "binary/work/src/run_binary_extras.f",
    "mesa_star_binary_extras": "binary/work/src/run_star_extras.f",
    "mesa_star1_extras": "star/work/src/run_star_extras.f",
    "mesa_star2_extras": "star/work/src/run_star_extras.f",
    "binary_run": "binary/work/src/binary_run.f",
    "star_run": "star/work/src/run.f",
}


def mesa_defaults_targets(grid_type, defaults):
    """Which ``.defaults`` file feeds which target, for a grid type.

    Parameters
    ----------
    grid_type : GridType
    defaults : dict
        From ``manifest.mesa_installation``.

    Returns
    -------
    dict
        Target to the file, in the order the v2 inis list them.
    """
    targets = {}
    stars = {"star1": grid_type.star1, "star2": grid_type.star2}
    for star, start in stars.items():
        if start == "formation" and not grid_type.single_star:
            for part in ("job", "controls"):
                targets["{}_formation_every_step_{}".format(star, part)] = \
                    defaults["star_" + part]
    targets["binary_controls"] = defaults["binary_controls"]
    targets["binary_job"] = defaults["binary_job"]
    for star in stars:
        steps = grid_type.single_star and star == "star1"
        for part in ("job", "controls"):
            target = "{}_formation_every_step_{}".format(star, part) \
                if steps else "{}_{}".format(star, part)
            targets[target] = defaults["star_" + part]
    return targets


def fill_defaults(config, resolved, mesa_dir):
    """Fill the ``mesa_inlists`` layer and the extras of a new-format ini.

    Parameters
    ----------
    config : GridConfig
        Changed in place. Legacy configs are left alone.
    resolved : ResolvedScenario or None
        From ``resolve.resolve_scenario``; its release names the MESA
        installation.
    mesa_dir : str or Path
        ``$MESA_DIR``, for the default extras.

    Raises
    ------
    ConfigError
        If MESA defaults are needed but no release names an installation.
    """
    if config.source_format != "new":
        return
    mesa_layer = config.inlists["mesa_inlists"]
    if not mesa_layer.targets:
        release = resolved.release if resolved is not None else None
        if release is None:
            raise ConfigError("mesa_inlists", None, "needed without a "
                              "[scenario] release, which would name the "
                              "MESA defaults")
        targets = mesa_defaults_targets(get_grid_type(config.grid.type),
                                        mesa_installation(release))
        for target, path in targets.items():
            mesa_layer.add(target, path)
    if not config.extras:
        config.extras = {key: Path(mesa_dir) / rel
                         for key, rel in MESA_EXTRAS.items()}
