"""A new-format ini sets up the same grid as its legacy ini.

For every grid type, ``posydon-setup-grid`` runs twice on the real
POSYDON-MESA-INLISTS:

- on a legacy ini: ``scenario = ['posydon', '<branch>-<commit>',
  '<system type>']`` at the production commit, MESA defaults named as
  ``${posydon_github_root}/grid_params/defaults/...``, which POSYDON no
  longer ships: the setup reads them from POSYDON-MESA-INLISTS at the
  pinned commit (``support.py``);
- on the equivalent new-format ini: ``[scenario] release = v2.0`` (+ the
  rerun of that commit) at the restructured repo, MESA defaults and
  extras filled in by ``defaults.fill_defaults``.

Every file written into the two run directories must be the same.
HeMS-HMS is the exception: its legacy ini never forms star 1 (tracker S2),
and its production commit (jeff_new_grid) is not on v2_sources, so only
the new ini is run and must form star 1.

Needs ``POSYDON_INLIST_REPO``, as test_release_manifests.py.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import re
import shutil
import subprocess
import sys

import pytest
from test_release_manifests import (
    INLIST_REPO,
    INLIST_REPO_COMMIT,
    REQUIRE_INLIST_REPO,
)
from test_setup_grid_baseline import CASES as GOLDEN_CASES
from test_setup_grid_baseline import inlist_origin  # noqa: F401
from test_setup_grid_baseline import (
    DATA,
    GRID_COLUMNS,
    GRID_ROWS,
    POSYDON_GRID,
    ROOT,
    build_ini,
    make_setup_env,
)

pytestmark = pytest.mark.skipif(
    not INLIST_REPO and not REQUIRE_INLIST_REPO,
    reason="set POSYDON_INLIST_REPO (see test_release_manifests.py)")

# the production commit the legacy inis name (v2.0 itself, no rerun)
MAIN = ("main", "a060b7bf1a4d94d693f77be9a1b0b3a522c1eadf")
ZAMS = "/posydon-tests/ZAMS_models/zams_z1.42m2_y0.2703.data"
CASES = [case for case in GOLDEN_CASES if case.name.endswith("__scenario")]


def clone(target):
    subprocess.run(["git", "clone", "-q", INLIST_REPO, str(target)],
                   check=True)
    return target


def new_ini(case, run, repo):
    """The new-format ini equivalent to the golden case's legacy ini."""
    # HeMS-HMS as Jeff set it up: v2.0 + the rerun of his branch
    rerun = "jeff_new_grid" if case.name.startswith("HeMS-HMS") else None
    out = case.output
    star_output = {
        star: ", ".join(option for option, key in (
            ("history", "history_" + star),
            ("final_profile", "final_profile_" + star),
            ("final_model", "final_model_" + star)) if out[key]) or "none"
        for star in ("star1", "star2")}
    lines = ["[grid]", "type = {}".format(case.name.split("__")[0]),
             "points = {}".format(run / "grid.csv"),
             "run_directory = {}".format(run)]
    if case.zams == "both":
        lines.append("zams = {}".format(ZAMS))
    elif case.zams == "star2":
        lines.append("zams_star2 = {}".format(ZAMS))
    lines += ["", "[scenario]", "release = v2.0",
              "commit = {}".format(INLIST_REPO_COMMIT),
              "inlist_repo = {}".format(repo)]
    if rerun:
        lines.append("rerun = {}".format(rerun))
    lines += ["", "[output]",
              "history_interval = {}".format(out["history_interval"]),
              "star1 = {}".format(star_output["star1"]),
              "star2 = {}".format(star_output["star2"]),
              "binary = {}".format("history" if out["binary_history"]
                                   else "none"),
              "", "[submission]", "backend = slurm_array",
              "partition = test-cpu", "account = test",
              "walltime = 2-00:00:00", "email = tester@posydon.invalid",
              "cpus_per_task = 4", "nodes = 1", "mpi_tasks = 1"]
    return "\n".join(lines) + "\n"


def run_dir_files(run):
    # setup_record.ini names the ini and scenario, which differ by design
    skip = {"grid.csv", "grid.ini", "binary/binary", "star1/star",
            "star2/star", "setup_record.ini"}
    return {path.relative_to(run).as_posix():
            path.read_bytes().replace(str(run).encode(), b"<RUN>")
            for path in sorted(run.rglob("*")) if path.is_file()
            and path.relative_to(run).as_posix() not in skip}


@pytest.mark.parametrize("case", CASES, ids=[c.name for c in CASES])
def test_new_ini_equals_legacy_ini(case, tmp_path,
                                   inlist_origin):  # noqa: F811
    assert INLIST_REPO, "POSYDON_INLIST_REPO is not set"
    grid_type = case.name.split("__")[0]
    env = make_setup_env(tmp_path, inlist_origin)
    shutil.rmtree(env["clone"])
    clone(env["clone"])  # the legacy scenario reads ~/.posydon_mesa_inlists
    branch, commit = MAIN
    runs = {}
    kinds = ("new",) if grid_type == "HeMS-HMS" else ("legacy", "new")
    for kind in kinds:
        run = tmp_path / kind
        run.mkdir()
        (run / "grid.csv").write_text("\n".join(
            [GRID_COLUMNS[case.columns]] + GRID_ROWS[case.columns]) + "\n")
        if kind == "legacy":
            text = re.sub(r"(?m)^scenario = .*$",
                          "scenario = ['posydon', '{}-{}', '{}']".format(
                              branch, commit, case.system_type),
                          build_ini(case, {**env, "run": run}))
            # as the v2 ini files: MESA defaults under grid_params/
            assert not (ROOT / "grid_params").exists()
            text = text.replace(str(DATA / "mesa_defaults"),
                                str(ROOT / "grid_params/defaults/r11701"))
            options = ["--grid-type", case.grid_type, "--run-directory",
                       str(run), "--submission-type", case.submission]
        else:
            text = new_ini(case, run, clone(tmp_path / "inlists_new"))
            options = []  # the ini says it all
        (run / "grid.ini").write_text(text)
        result = subprocess.run(
            [sys.executable, str(POSYDON_GRID), "setup", "grid.ini"]
            + options, cwd=run, env=env["env"], capture_output=True, text=True,
            timeout=300)
        assert result.returncode == 0, (kind, result.stderr[-3000:])
        runs[kind] = run_dir_files(run)
    if grid_type == "HeMS-HMS":
        assert any(name.startswith("star1/inlist_step")
                   for name in runs["new"])
        return
    # inlists, column lists, sources, mk and the submit scripts
    assert len(runs["legacy"]) >= 20
    assert sorted(runs["new"]) == sorted(runs["legacy"])
    for name in runs["legacy"]:
        assert runs["new"][name] == runs["legacy"][name], name
