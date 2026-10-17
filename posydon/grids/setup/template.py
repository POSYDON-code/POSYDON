"""A starter new-format ini for a grid type (``posydon-grid init``).

The template lists every option the new format takes for the grid type,
section by section, each with a short note. Options that have a default or
are optional are written commented out (``; key = value``): remove the
``;`` to use one. Options that do not apply to the grid type are left out,
e.g. star 2's inlists in a single-star grid or formation steps for a star
that starts from ZAMS.

What the grid type decides is filled in (its type, what MESA writes for
each star). What only the user knows is left as ``<...>`` placeholders: the
SLURM account and partition, and a ZAMS model for stars that start from
ZAMS. The scenario names the release; the inlist repository's clone is
read as it is.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

from posydon.grids.setup.config import EXTRAS_KEYS
from posydon.grids.setup.grid_types import get_grid_type

# column of the inline comments, and the longest line
COMMENT_COLUMN = 28
WIDTH = 79

ZAMS_NOTE = "e.g. <inlist repo>/r11701/ZAMS_models/zams_z1.42m2_y0.2703.data"


def _with_note(text, note):
    """``text  ; note``, or the note on a line of its own above ``text``
    when the two do not fit in ``WIDTH`` columns."""
    if not note:
        return text + "\n"
    inline = text.ljust(COMMENT_COLUMN - 1) + " ; " + note
    if len(inline) <= WIDTH:
        return inline + "\n"
    return "; " + note + "\n" + text + "\n"


def _line(key, value, note=None, active=True):
    """``key = value  ; note``, commented out unless ``active``."""
    return _with_note("{}{} = {}".format("" if active else "; ", key, value),
                      note)


def _header(name, *notes):
    """``[name]`` and its notes, one comment line each."""
    return "[{}]\n".format(name) + "".join("; " + note + "\n"
                                           for note in notes)


def _grid(grid_type):
    text = _header("grid", "what to run, and where")
    text += _line("type", grid_type.name, "see posydon-grid types")
    text += _line("mode", "fixed", "fixed | dynamic (psy-cris samples)")
    text += _line("points", "grid.csv", "a CSV, a PSyGrid .h5 or a folder "
                  "of runs")
    text += _line("run_directory", ".", "where executables, inlists and "
                  "output go")
    text += _line("psycris_inifile", "<psy-cris ini>", "needed for mode = "
                  "dynamic", active=False)
    zams_stars = [star for star in ("star1", "star2")
                  if getattr(grid_type, star) == "zams"]
    if len(zams_stars) == 2:
        text += _line("zams", "<ZAMS model>", "both stars; " + ZAMS_NOTE)
        text += _line("zams_star1", "<ZAMS model>", "instead of zams, one "
                      "per star", active=False)
        text += _line("zams_star2", "<ZAMS model>", active=False)
    elif zams_stars:
        text += _line("zams_" + zams_stars[0], "<ZAMS model>", ZAMS_NOTE)
    text += _line("keep_profiles", "false", "keep each run's profiles",
                  active=False)
    text += _line("keep_photos", "false", "keep each run's photos",
                  active=False)
    return text


def _scenario():
    text = _header("scenario", "which POSYDON inlists to use, from a "
                   "clone of POSYDON-MESA-INLISTS.",
                   "The setup reads the clone as it is; keeping it up to "
                   "date (git pull) is",
                   "up to you. setup_record.ini notes its commit and the "
                   "files used.")
    text += _line("release", "v2.0", "r11701/releases/<release>/")
    text += _line("rerun", "matthias_lbv_thermo_dedt_noAgeLimit",
                  "a physics variant: reruns/<rerun>/", active=False)
    text += _line("inlist_repo", "~/.posydon_mesa_inlists", "the clone, "
                  "on the v2_sources branch", active=False)
    text += _line("mesa", "r11701", "if the release exists for several "
                  "MESA versions", active=False)
    return text


def _inlist_keys(grid_type):
    """``(key, note)`` of every inlist-layer key that applies."""
    single = grid_type.single_star
    keys = []
    if not single:
        keys += [("binary", "&binary_controls and &binary_job"),
                 ("binary_controls", "only &binary_controls"),
                 ("binary_job", "only &binary_job")]
    keys += [("star1", "&controls and &star_job of star 1" +
              (", every step" if single else "")),
             ("star1_controls", "only &controls of star 1"),
             ("star1_job", "only &star_job of star 1")]
    if grid_type.star2 in ("zams", "formation"):
        keys += [("star2", "&controls and &star_job of star 2"),
                 ("star2_controls", "only &controls of star 2"),
                 ("star2_job", "only &star_job of star 2")]
    for star in ("star1", "star2"):
        if getattr(grid_type, star) == "formation":
            keys += [("{}_formation".format(star),
                      "one file per formation step, in order"),
                     ("{}_formation_step0".format(star),
                      "on top of formation step 0 (one key per step)")]
    keys += [("history_columns", "whole files: replace the release's"),
             ("binary_history_columns", None),
             ("profile_columns", None),
             ("run_star_extras", "Fortran, any file name; " +
              ("star 1" if single else "every star")),
             ]
    if not single:
        keys += [("run_binary_extras", None),
                 ("run_star1_extras", "only star 1, over run_star_extras")]
    if grid_type.star2 in ("zams", "formation"):
        keys += [("run_star2_extras", "only star 2, over run_star_extras")]
    return keys


def _inlists(grid_type):
    keys = _inlist_keys(grid_type)
    text = _header("mesa_inlists", "MESA's defaults; by default those of "
                   "the release's MESA installation.",
                   "Same keys as [user_inlists] but no Fortran sources; a "
                   "file given for", "a namelist replaces MESA's defaults "
                   "there.")
    text += _line(keys[0][0], "<file>", active=False)
    text += "\n" + _header("user_inlists", "Your own MESA inlists, on top "
                           "of the release, key by key.",
                           "Several files for one key: <file>, <file> (the "
                           "last wins).")
    for key, note in keys:
        text += _line(key, "<file>", note, active=False)
    return text


def _output(grid_type):
    text = _header("output", "what MESA writes")
    text += _line("history_interval", "1", "MESA steps between history "
                  "lines")
    stars = ("star1",) if grid_type.single_star else ("star1", "star2")
    for star in stars:
        evolved = getattr(grid_type, star) in ("zams", "formation")
        text += _line(star, "history, final_model" if evolved else "none",
                      "history, final_profile, final_model or none"
                      if star == "star1" else None)
    if not grid_type.single_star:
        text += _line("binary", "history", "history | none")
    return text


def _submission():
    text = _header("submission", "how the grid is run")
    text += _line("backend", "slurm_array", "shell | slurm_array | "
                  "slurm_mpi")
    text += _line("job_array", "true", "shell: one by one, or false: one "
                  "mpirun", active=False)
    text += "; SLURM settings, left out of the scripts if not set:\n"
    text += _line("partition", "<partition>")
    text += _line("account", "<account>")
    text += _line("walltime", "2-00:00:00")
    text += _line("email", "<email>")
    text += _line("cpus_per_task", "4", "OpenMP threads per MESA run")
    text += _line("mem_per_cpu", "4G", active=False)
    text += _line("nodes", "1", "slurm_mpi: nodes, and MPI tasks per node",
                  active=False)
    text += _line("mpi_tasks", "1", active=False)
    text += _line("work_dir", "<scratch folder>", "run each point in a "
                  "temporary folder there", active=False)
    text += _line("newgroup", "<group>", "hand the output to this Unix "
                  "group", active=False)
    return text


def _extras():
    text = _header("extras", "MESA's own work files; each one not set "
                   "comes from $MESA_DIR,", "e.g. makefile_binary: "
                   "binary/work/make/makefile. Fortran sources go in",
                   "[user_inlists].")
    for key in EXTRAS_KEYS:
        text += _line(key, "<file>", active=False)
    return text


def ini_template(name):
    """The starter ini text for a grid type.

    Examples
    --------
    >>> "star2 = none" in ini_template("CO-HMS_RLO")
    True
    >>> "; star1_formation_step0 = <file>" in ini_template("CO-HeMS")
    True
    """
    grid_type = get_grid_type(name)
    head = ("; posydon-grid settings for a {0} grid:\n"
            "; {1}.\n"
            "; Written by posydon-grid init {0}.\n"
            "; Paths are relative to this file.\n"
            "; Every option for this grid type is listed; those starting "
            "with ';' are\n"
            "; left at their default: remove the ';' to set one.\n"
            "; Check with: posydon-grid explain <this file>\n"
            "; Set up with: posydon-grid setup <this file>\n").format(
                grid_type.name, grid_type.description)
    return "\n".join([head, _grid(grid_type), _scenario(),
                      _inlists(grid_type), _output(grid_type),
                      _submission(), _extras()])
