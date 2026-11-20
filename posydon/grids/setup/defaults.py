"""Fill in what a new-format ini leaves to the release and MESA.

A legacy ini spells out its MESA defaults and Fortran sources. A
new-format ini may leave them out; they are then taken as every v2
production ini did:

- ``mesa_inlists``: MESA's ``.defaults`` of the release's MESA
  installation (``<mesa_version>/mesa/<installation>/defaults/``) for the
  binary, star 1 and star 2, and for every formation step of a star the
  grid type builds in steps, for every target ``[mesa_inlists]`` does not
  name files for. A single-star grid runs star 1 as its formation steps,
  so star 1's defaults go there.
- ``extras``: MESA's own work files in ``$MESA_DIR`` (makefiles,
  ``run.f``, ``binary_run.f`` and the stock extras), under the legacy
  ``[mesa_extras]`` names, for every name the ini's ``[extras]`` does not
  set; the scenario's sources then replace the stock extras, as before.

Anything the ini sets is kept as it is.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

from pathlib import Path

from posydon.grids.setup.config import ConfigError
from posydon.grids.setup.grid_types import get_grid_type
from posydon.grids.setup.layers import formation_step_count
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


def _per_step(targets, config, scenario):
    """Turn every-step targets into one target per formation step.

    So that ``[mesa_inlists] star1_formation_step0`` replaces MESA's
    defaults for that step, as any other target's file does, instead of
    being stacked on top of them.
    """
    layers = [files.targets for files in config.inlists.values()]
    layers += [scenario.inlists.targets, scenario.flags]
    out = {}
    for target, path in targets.items():
        star, every, part = target.partition("_formation_every_step_")
        if not every:
            out[target] = path
            continue
        for step in range(max(formation_step_count(star, layers), 1)):
            out["{}_formation_step{}_{}".format(star, step, part)] = path
    return out


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
    release = resolved.release if resolved is not None else None
    if release is None:
        # without a release, the ini's [mesa_inlists] is all there is
        if not mesa_layer.targets:
            raise ConfigError("mesa_inlists", None, "needed without a "
                              "[scenario] release, which would name the "
                              "MESA defaults")
    else:
        # target by target: an ini that names the defaults of one target
        # still gets MESA's for the others (before, naming any left out
        # all of them)
        grid_type = get_grid_type(config.grid.type)
        targets = mesa_defaults_targets(grid_type,
                                        mesa_installation(release))
        targets = _per_step(targets, config, resolved.files)
        given = set(mesa_layer.targets)
        if grid_type.single_star:
            # a single star runs as its steps: [mesa_inlists] star1 is
            # every step's defaults
            for part in ("controls", "job"):
                if "star1_{}".format(part) in given:
                    given |= {target for target in targets if
                              target.startswith("star1_formation_step") and
                              target.endswith("_" + part)}
        for target, path in targets.items():
            if target not in given:
                mesa_layer.add(target, path)
    # key by key: an ini that names one source still gets MESA's makefiles
    # and the others (before, naming any [extras] key left out all of them,
    # so nothing was compiled)
    filled = tuple(key for key in MESA_EXTRAS if key not in config.extras)
    if filled:
        given = config.extras
        config.extras = {key: given.get(key, Path(mesa_dir) / rel)
                         for key, rel in MESA_EXTRAS.items()}
        config.extras.update({key: path for key, path in given.items()
                              if key not in MESA_EXTRAS})
        config.default_extras = filled
