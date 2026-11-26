"""The example ini files of the docs read, and set up a grid.

``docs/_source/components-overview/mesa_grids/examples/`` holds the
new-format ini files the documentation links to. Each one is read, and set
up against the synthetic inlist repository of ``test_setup_matrix``, with
only what points at the user's machine replaced: the clone, the ZAMS model
and the rerun's name.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

from pathlib import Path

import pytest

import posydon.grids.setup.steps as totest
from posydon.grids.setup.config import read_config
from posydon.grids.setup.grid_types import GRID_TYPES
from posydon.unit_tests.grids.test_setup_matrix import (  # noqa: F401
    env,
    inlist,
    write,
)

EXAMPLES = Path(__file__).resolve().parents[3] / "docs" / "_source" / \
    "components-overview" / "mesa_grids" / "examples"
NAMES = sorted(path.stem for path in EXAMPLES.glob("*.ini"))
CLONE = "~/POSYDON-MESA-INLISTS"
ZAMS = CLONE + "/r11701/ZAMS_models/zams_z1.42m2_y0.2703.data"
# the scripts each backend writes
SCRIPTS = {"slurm_array": "job_array_grid_submit.slurm",
           "slurm_mpi": "mpi_grid_submit.slurm", "shell": "grid_command.sh"}


def test_there_are_examples():
    assert EXAMPLES.is_dir()
    assert (EXAMPLES / "README.md").is_file()
    # one per grid type, at least
    assert set(GRID_TYPES) <= set(NAMES)


def test_readme_lists_them_all():
    readme = (EXAMPLES / "README.md").read_text()
    for name in NAMES:
        assert "[{0}.ini]({0}.ini)".format(name) in readme, name


@pytest.mark.parametrize("name", NAMES)
def test_reads(name):
    config = read_config(EXAMPLES / (name + ".ini"))
    grid_type = GRID_TYPES[config.grid.type]
    assert name.split("_")[0] in (grid_type.name, grid_type.name.split(
        "_")[0])
    assert config.scenario.release == "v2.0"
    assert config.scenario.inlist_repo == Path(CLONE).expanduser()
    # a ZAMS model for exactly the stars that start from one
    zams = {star for star in ("star1", "star2")
            if getattr(config.grid, "zams_" + star) is not None}
    assert zams == {star for star in ("star1", "star2")
                    if getattr(grid_type, star) == "zams"}
    for star in zams:
        assert str(getattr(config.grid, "zams_" + star)) == \
            str(Path(ZAMS).expanduser())


@pytest.mark.parametrize("name", NAMES)
def test_sets_up(env, name):  # noqa: F811
    text = (EXAMPLES / (name + ".ini")).read_text()
    text = text.replace(ZAMS, str(env["zams"]))
    text = text.replace("inlist_repo = " + CLONE,
                        "inlist_repo = {}".format(env["repo"]))
    # the synthetic repository's rerun
    text = text.replace("rerun = matthias_lbv_thermo_dedt_noAgeLimit",
                        "rerun = hot")
    assert CLONE not in text
    # stand-ins for the user's own files
    for user_file in ("my_inlist_project", "my_inlist1",
                      "my_inlist2_controls"):
        write(env["tmp"] / user_file, inlist({"u": "1"}, {"u": "1"},
                                             binary="project" in user_file))
    write(env["tmp"] / "my_run_star_extras.f", "! mine\n")
    write(env["tmp"] / "my_history_columns.list", "mine\n")
    ini = write(env["tmp"] / (name + ".ini"), text)
    run = totest.setup_grid(ini)
    config = read_config(ini)
    assert run == env["tmp"] / "run"
    assert (run / SCRIPTS[config.submission.backend]).is_file()
    assert (run / "setup_record.ini").is_file()
