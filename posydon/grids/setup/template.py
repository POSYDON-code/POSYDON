"""A starter new-format ini for a grid type (``--init``).

The template names every section with comments and fills what the grid
type decides (its type, what MESA writes for each star). What only the user
knows is left as ``<...>`` placeholders: the SLURM account and partition,
and a ZAMS model for stars that start from ZAMS.

The scenario is written out in full, with no hidden default: the release,
a commented-out rerun and the inlist repository's ref.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

from posydon.grids.setup.grid_types import get_grid_type

TEMPLATE = """\
; posydon-grid settings for a {name} grid:
; {description}.
; Written by posydon-grid init {name}. Paths are relative to this file.
; Check with: posydon-grid explain <this file>
; Set up with: posydon-grid setup <this file>

[grid]
type = {name}
mode = fixed                ; fixed | dynamic (dynamic needs psycris_inifile)
points = grid.csv           ; the grid points, one row per MESA run
run_directory = .           ; where executables, inlists and output go
{zams}
[scenario]                  ; which POSYDON inlists to use
release = v2.0
; rerun = matthias_lbv_thermo_dedt_noAgeLimit  ; a physics variant of v2.0
ref = v2_sources            ; branch or tag of POSYDON-MESA-INLISTS;
                            ; pin a production grid with commit = <hash>
; source = path             ; use inlist_repo as it is, without git
; inlist_repo = ~/.posydon_mesa_inlists

[user_inlists]              ; your own MESA inlists, on top of the release
; binary = my_inlist_project
; star1_controls = my_inlist1

[output]
history_interval = 1
star1 = {star1}
star2 = {star2}
binary = history            ; history | none

[submission]
backend = slurm_array       ; shell | slurm_array | slurm_mpi
partition = <partition>
account = <account>
walltime = 2-00:00:00
email = <email>
cpus_per_task = 4
"""


def ini_template(name):
    """The starter ini text for a grid type.

    Examples
    --------
    >>> "star2 = none" in ini_template("CO-HMS_RLO")
    True
    """
    grid_type = get_grid_type(name)
    zams_stars = [star for star in ("star1", "star2")
                  if getattr(grid_type, star) == "zams"]
    if zams_stars:
        key = "zams" if len(zams_stars) == 2 else "zams_" + zams_stars[0]
        zams = ("{} = <ZAMS model, e.g. POSYDON-MESA-INLISTS/r11701/"
                "ZAMS_models/zams_z1.42m2_y0.2703.data>\n".format(key))
    else:
        zams = ""
    outputs = {}
    for star in ("star1", "star2"):
        evolved = getattr(grid_type, star) in ("zams", "formation") and \
            not (grid_type.single_star and star == "star2")
        outputs[star] = "history, final_model" if evolved else "none"
    return TEMPLATE.format(name=grid_type.name,
                           description=grid_type.description, zams=zams,
                           star1=outputs["star1"], star2=outputs["star2"])
