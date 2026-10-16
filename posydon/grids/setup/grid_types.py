"""The grid types POSYDON can set up, and what makes each one different.

Every grid type is one ``GridType`` entry in ``GRID_TYPES``. It records

- ``folder``: the grid type's own folder in POSYDON-MESA-INLISTS
  (``r11701/<folder>/``), whose files are layered on top of
  ``default_common_inlists``; ``None`` means common inlists only,
- how each star starts: ``"zams"`` (a ZAMS model file), ``"formation"``
  (built from formation steps), ``"point_mass"`` (a compact object, not
  evolved by MESA) or ``None`` (no second star),
- ``single_star``: only star 1 is evolved, without a binary,
- ``single_star_flag``: whether the single-star flag
  ``x_logical_ctrl(1) = .true.`` is set for star 1 (``run_star_extras.f``
  reads it to evolve the star on its own). A star that starts from ZAMS
  gets the flag on its one inlist; a star built from formation steps gets
  it as an extra, last step after the formation (as in the v2
  single_HeMS grids).

The table::

    grid type    folder      star 1     star 2      single star
    HMS-HMS      -           zams       zams        no
    CO-HMS_RLO   CO-H_star   zams       point_mass  no
    CO-HeMS      CO-He_star  formation  point_mass  no
    HeMS-HMS     HeMS-HMS    formation  zams        no
    single_HMS   -           zams       -           yes
    single_HeMS  CO-He_star  formation  -           yes

Adding a grid type is adding a row here; the rest of setup reads the
registry instead of testing grid-type names.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

from dataclasses import dataclass

from posydon.grids.setup.config import ConfigError

STAR_STARTS = ("zams", "formation", "point_mass", None)


@dataclass(frozen=True)
class GridType:
    """What sets one grid type apart."""

    name: str
    description: str
    folder: str = None
    star1: str = "zams"
    star2: str = "zams"
    single_star: bool = False
    single_star_flag: bool = False

    def __post_init__(self):
        for star in ("star1", "star2"):
            if getattr(self, star) not in STAR_STARTS:
                raise ValueError("{}: unknown start {!r} for {}".format(
                    self.name, getattr(self, star), star))
        if self.single_star != (self.star2 is None):
            raise ValueError("{}: a single-star grid has no star 2, and "
                             "a binary grid needs one".format(self.name))


GRID_TYPES = {grid_type.name: grid_type for grid_type in (
    GridType("HMS-HMS", "two hydrogen-rich main-sequence stars"),
    GridType("CO-HMS_RLO", "a compact object and a hydrogen-rich "
             "main-sequence star", folder="CO-H_star", star2="point_mass"),
    GridType("CO-HeMS", "a compact object and a helium main-sequence star",
             folder="CO-He_star", star1="formation", star2="point_mass"),
    GridType("HeMS-HMS", "a helium and a hydrogen-rich main-sequence star",
             folder="HeMS-HMS", star1="formation"),
    GridType("single_HMS", "single hydrogen-rich main-sequence stars",
             star2=None, single_star=True, single_star_flag=True),
    GridType("single_HeMS", "single helium main-sequence stars",
             folder="CO-He_star", star1="formation", star2=None,
             single_star=True, single_star_flag=True),
)}


def get_grid_type(name):
    """Look up a grid type by name.

    Raises
    ------
    ConfigError
        For an unknown name, listing the known ones.

    Examples
    --------
    >>> get_grid_type("CO-HeMS").folder
    'CO-He_star'
    """
    try:
        return GRID_TYPES[name]
    except KeyError:
        raise ConfigError("grid", "type", "unknown grid type {!r}; expected "
                          "one of {}".format(name, ", ".join(GRID_TYPES)))


def describe_grid_types():
    """One line per grid type, for ``posydon-setup-grid --list-grid-types``.

    Examples
    --------
    >>> print(describe_grid_types().splitlines()[0])
    HMS-HMS       two hydrogen-rich main-sequence stars
    """
    return "\n".join("{:13s} {}".format(name, grid_type.description)
                     for name, grid_type in GRID_TYPES.items())
