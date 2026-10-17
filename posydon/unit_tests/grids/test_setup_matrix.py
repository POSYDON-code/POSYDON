"""The configuration matrix of posydon-grid setup, checked end to end.

Every grid type, backend, ZAMS choice, layer and output option of the new
ini format (``grid_rewrite_tests/CONFIG_MATRIX.md``) is set up against a
synthetic inlist repository laid out like POSYDON-MESA-INLISTS v2.0. Each
file of that repository sets a marker key, so the tests can say which files
reached which inlist. The expectations come from the matrix and the ini
documentation, not from the code.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import configparser
from pathlib import Path

import pytest

import posydon.grids.setup.steps as totest
from posydon.grids.setup.config import ConfigError
from posydon.grids.setup.defaults import MESA_EXTRAS
from posydon.grids.setup.inlist import parse_inlist

GRID_TYPES = ("HMS-HMS", "CO-HMS_RLO", "CO-HeMS", "HeMS-HMS", "single_HMS",
              "single_HeMS")
SINGLE = ("single_HMS", "single_HeMS")
# number of formation steps of star 1, by grid type (matrix table)
STAR1_STEPS = {"HMS-HMS": 0, "CO-HMS_RLO": 0, "CO-HeMS": 2, "HeMS-HMS": 2,
               "single_HMS": 1, "single_HeMS": 3}
# how each star starts
STARTS = {"HMS-HMS": ("zams", "zams"), "CO-HMS_RLO": ("zams", "point_mass"),
          "CO-HeMS": ("formation", "point_mass"),
          "HeMS-HMS": ("formation", "zams"), "single_HMS": ("zams", "none"),
          "single_HeMS": ("formation", "none")}


###############################################################################
# a synthetic POSYDON-MESA-INLISTS, v2.0 layout
###############################################################################
def inlist(controls=None, job=None, binary=False):
    names = ("binary_controls", "binary_job") if binary else \
        ("controls", "star_job")
    text = ""
    for name, values in zip(names, (controls or {}, job or {})):
        text += "&{}\n".format(name) + "".join(
            "  {} = {}\n".format(key, value)
            for key, value in values.items()) + "/\n"
    return text


def manifest(name, star1, star2, inlists):
    return ("[grid_type]\nname = {}\nstar1 = {}\nstar2 = {}\n\n[inlists]\n"
            .format(name, star1, star2)
            + "".join("{} = {}\n".format(key, value)
                      for key, value in inlists.items())
            + "\n[columns]\n"
            "history_columns = common/history_columns.list\n"
            "binary_history_columns = common/binary_history_columns.list\n"
            "profile_columns = common/profile_columns.list\n"
            "\n[extras]\n"
            "posydon_binary_extras = common/binary/src/run_binary_extras.f\n"
            "posydon_star_binary_extras = common/binary/src/"
            "run_star_extras.f\n"
            "mesa_star1_extras = common/binary/src/run_star_extras.f\n"
            "\n[zams]\n")


MANIFESTS = {
    "HMS-HMS": manifest("HMS-HMS", "zams", "zams", {
        "binary": "common/binary/inlist_project",
        "star1": "common/binary/inlist1",
        "star2": "common/binary/inlist2"}),
    "CO-HMS_RLO": manifest("CO-HMS_RLO", "zams", "point_mass", {
        "binary": "common/binary/inlist_project, "
                  "CO-HMS_RLO/binary/inlist_project",
        "star1": "common/binary/inlist1, CO-HMS_RLO/binary/inlist1",
        "star2": "common/binary/inlist2, CO-HMS_RLO/binary/inlist2"}),
    "CO-HeMS": manifest("CO-HeMS", "formation", "point_mass", {
        "binary": "common/binary/inlist_project, "
                  "CO-HeMS/binary/inlist_project",
        "star1": "common/binary/inlist1, CO-HeMS/binary/inlist1",
        "star2": "common/binary/inlist2",
        "star1_formation_every_step": "common/binary/inlist1",
        "star1_formation_step0": "CO-HeMS/star1_formation/inlist_step1",
        "star1_formation_step1": "CO-HeMS/star1_formation/inlist_step2"}),
    "HeMS-HMS": manifest("HeMS-HMS", "formation", "zams", {
        "binary": "common/binary/inlist_project",
        "star1": "common/binary/inlist1, HeMS-HMS/binary/inlist1",
        "star2": "common/binary/inlist2",
        "star1_formation_every_step": "common/binary/inlist1",
        "star1_formation_step0": "HeMS-HMS/star1_formation/inlist_step1",
        "star1_formation_step1": "HeMS-HMS/star1_formation/inlist_step2"}),
    "single_HMS": manifest("single_HMS", "zams", "none", {
        "star1_formation_every_step": "common/binary/inlist1, "
                                      "single_HMS/single_star_flag"}),
    "single_HeMS": manifest("single_HeMS", "formation", "none", {
        "star1_formation_every_step": "common/binary/inlist1",
        "star1_formation_step0": "CO-HeMS/star1_formation/inlist_step1",
        "star1_formation_step1": "CO-HeMS/star1_formation/inlist_step2",
        "star1_formation_step2": "single_HeMS/single_star_flag"}),
}

# MESA's defaults: every key a later layer may override, plus one each
# layer leaves alone
MESA_FILES = {
    "binary/binary_controls.defaults": "fm = 1d-2\nm_mesa_bc = 1\n"
                                       "history_interval = 1\n",
    "binary/binary_job.defaults": "evolve_both_stars = .true.\n"
                                  "m_mesa_bj = 1\n",
    "star/controls.defaults": "mesh_delta_coeff = 1d0\nm_mesa_c = 1\n"
                              "initial_mass = 1d0\nzams_filename = ''\n"
                              "x_ctrl(1:num_x_ctrls) = 0d0\n",
    "star/star_job.defaults": "pgstar_flag = .false.\nm_mesa_j = 1\n",
}


def write(path, text=""):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def make_inlist_repo(root):
    """POSYDON-MESA-INLISTS with one release, v2.0, and a rerun."""
    r11701 = root / "r11701"
    for rel, text in MESA_FILES.items():
        write(r11701 / "mesa" / "default" / "defaults" / rel,
              text + "lv = 'mesa_default'\n")
    write(r11701 / "ZAMS_models" / "zams.data", "zams")
    release = r11701 / "releases" / "v2.0"
    write(release / "release.ini",
          "[release]\nname = v2.0\nmesa_version = r11701\nmesa = default\n"
          "frozen = true\n")
    common = release / "common"
    write(common / "binary" / "inlist_project", inlist(
        {"fm": "5d-3", "m_common_project": "1", "lv": "'release'"},
        {"m_common_project_job": "1", "lv": "'release'"}, binary=True))
    write(common / "binary" / "inlist1", inlist(
        {"mesh_delta_coeff": "0.8d0", "m_common_star1": "1",
         "lv": "'release'"},
        {"m_common_star1_job": "1", "lv": "'release'"}))
    write(common / "binary" / "inlist2", inlist(
        {"mesh_delta_coeff": "0.9d0", "m_common_star2": "1",
         "lv": "'release'"},
        {"m_common_star2_job": "1", "lv": "'release'"}))
    for name in ("history_columns", "binary_history_columns",
                 "profile_columns"):
        write(common / (name + ".list"), "release " + name + "\n")
    for name in ("run_binary_extras.f", "run_star_extras.f"):
        write(common / "binary" / "src" / name, "! release " + name + "\n")
    for folder in ("CO-HMS_RLO", "CO-HeMS"):
        write(release / folder / "binary" / "inlist_project", inlist(
            {"m_{}_project".format(folder): "1"}, {}, binary=True))
        write(release / folder / "binary" / "inlist1", inlist(
            {"m_{}_star1".format(folder): "1"}, {}))
    write(release / "CO-HMS_RLO" / "binary" / "inlist2", inlist(
        {"m_CO-HMS_RLO_star2": "1"}, {}))
    write(release / "HeMS-HMS" / "binary" / "inlist1", inlist(
        {"m_HeMS-HMS_star1": "1"}, {}))
    for folder in ("CO-HeMS", "HeMS-HMS"):
        for step in (1, 2):
            write(release / folder / "star1_formation" /
                  "inlist_step{}".format(step), inlist(
                      {"m_{}_step{}".format(folder, step): "1"},
                      {"m_{}_step{}_job".format(folder, step): "1"}))
    for name in SINGLE:
        write(release / name / "single_star_flag",
              inlist({"x_logical_ctrl(1)": ".true."}, {}))
    for name, text in MANIFESTS.items():
        write(release / name / "scenario.ini", text)
    rerun = release / "reruns" / "hot"
    write(rerun / "scenario.ini",
          "[overlays]\ncommon/binary/inlist1 = common/binary/inlist1\n")
    write(rerun / "common" / "binary" / "inlist1",
          inlist({"mesh_delta_coeff": "0.7d0", "m_rerun": "1",
                  "lv": "'rerun'"}, {"lv": "'rerun'"}))
    return root


@pytest.fixture
def env(tmp_path, monkeypatch):
    """The inlist repository, a stub MESA and the environment."""
    repo = make_inlist_repo(tmp_path / "inlists")
    mesa = tmp_path / "mesa"
    for key, rel in MESA_EXTRAS.items():
        # makefiles that build nothing, so mk succeeds
        write(mesa / rel, "all:\n\t@true\n" if "makefile" in rel
              else "! MESA " + rel + "\n")
    monkeypatch.setenv("MESA_DIR", str(mesa))
    monkeypatch.setenv("MESASDK_ROOT", "/sdk")
    monkeypatch.setattr(totest, "find_run_grid", lambda: "/bin/run-grid")
    import posydon.active_learning.psy_cris.utils as utils
    monkeypatch.setattr(utils, "parse_inifile", lambda path: {
        "posydon_dynamic_sampling_kwargs": {
            "mesa_column_names": ["initial_mass", "m2"]}})
    write(tmp_path / "grid.csv", "initial_mass,m2\n1,2\n3,4\n5,6\n")
    write(tmp_path / "psycris.ini", "")
    return {"tmp": tmp_path, "repo": repo, "mesa": mesa,
            "zams": repo / "r11701" / "ZAMS_models" / "zams.data"}


###############################################################################
# writing an ini and reading the run directory
###############################################################################
def ini_text(env, grid_type, grid=None, scenario=None, sections=None,
             submission=None):
    """A new-format ini: the release from env's repository by default."""
    grid_lines = {"type": grid_type, "points": env["tmp"] / "grid.csv",
                  "run_directory": env["tmp"] / "run"}
    grid_lines.update(grid or {})
    scenario_lines = {"source": "path", "release": "v2.0",
                      "inlist_repo": env["repo"]}
    if scenario is not None:
        scenario_lines.update(scenario)
    submission_lines = {"backend": "shell", "job_array": "true"}
    submission_lines.update(submission or {})
    blocks = [("grid", grid_lines), ("scenario", scenario_lines),
              ("submission", submission_lines)]
    blocks += list((sections or {}).items())
    return "".join("[{}]\n".format(name) + "".join(
        "{} = {}\n".format(key, value) for key, value in values.items()
        if value is not None) + "\n" for name, values in blocks)


def setup(env, grid_type, **options):
    ini = write(env["tmp"] / "grid.ini",
                ini_text(env, grid_type, **options))
    run = totest.setup_grid(ini)
    return Run(run)


class Run:
    """A written run directory, read back."""

    def __init__(self, path):
        self.path = Path(path)

    def inlist(self, rel):
        return parse_inlist((self.path / rel).read_text())

    def steps(self, star):
        return sorted((self.path / star).glob("inlist_step*"),
                      key=lambda path: int(path.name[len("inlist_step"):]))

    def text(self, rel):
        return (self.path / rel).read_text()

    def files(self):
        return {path.relative_to(self.path).as_posix()
                for path in self.path.rglob("*") if path.is_file()}

    def star_runs_in(self, grid_type, star="star1"):
        """The inlist that evolves the star to the end of the run."""
        if grid_type in SINGLE:
            return self.steps("star1")[-1].relative_to(self.path).as_posix()
        return "binary/inlist{}".format(star[-1])


def user_file(env, name, controls=None, job=None, binary=False):
    return write(env["tmp"] / "user" / name,
                 inlist(controls, job, binary=binary))


###############################################################################
# the grid types: structure, release files, MESA defaults
###############################################################################
@pytest.mark.parametrize("grid_type", GRID_TYPES)
class TestGridTypes:
    def test_files_written(self, env, grid_type):
        run = setup(env, grid_type)
        assert len(run.steps("star1")) == STAR1_STEPS[grid_type]
        assert run.steps("star2") == []
        for rel in ("binary/inlist_project", "binary/inlist1",
                    "binary/inlist2", "mk", "grid_command.sh",
                    "setup_record.ini", "column_lists/history_columns.list",
                    "column_lists/binary_history_columns.list",
                    "column_lists/profile_columns.list",
                    "binary/src/run_binary_extras.f",
                    "binary/src/run_star_extras.f",
                    "star1/src/run_star_extras.f",
                    "star2/src/run_star_extras.f",
                    "binary/make/makefile_binary",
                    "star1/make/makefile_star", "star2/make/makefile_star",
                    "binary/src/binary_run.f", "star1/src/run.f"):
            assert rel in run.files(), rel

    def test_command_line_names_the_steps(self, env, grid_type):
        run = setup(env, grid_type)
        command = run.text("grid_command.sh")
        steps = " ".join(str(path) for path in run.steps("star1"))
        if steps:
            assert "--mesa-star1-inlist-project  {} ".format(steps) in \
                command
        else:
            assert "--mesa-star1-inlist-project None " in command
        assert "--mesa-star2-inlist-project None " in command
        for _ in range(3):
            assert "for SLURM_ARRAY_TASK_ID in 0 1 2 ; do" in command

    def test_release_files_reach_their_inlists(self, env, grid_type):
        run = setup(env, grid_type)
        star1 = run.inlist(run.star_runs_in(grid_type))
        assert star1["controls"]["m_common_star1"] == "1"
        if grid_type in SINGLE:
            # the single-star flag, in the step that runs the star
            assert star1["controls"]["x_logical_ctrl(1)"] == ".true."
            return
        project = run.inlist("binary/inlist_project")
        assert project["binary_controls"]["m_common_project"] == "1"
        assert project["binary_controls"]["fm"] == "5d-3"
        assert run.inlist("binary/inlist2")["controls"][
            "m_common_star2"] == "1"
        if grid_type in ("CO-HMS_RLO", "CO-HeMS"):
            assert project["binary_controls"][
                "m_{}_project".format(grid_type)] == "1"
            assert star1["controls"]["m_{}_star1".format(grid_type)] == "1"
        if grid_type == "CO-HMS_RLO":
            assert run.inlist("binary/inlist2")["controls"][
                "m_CO-HMS_RLO_star2"] == "1"

    def test_mesa_defaults_under_everything(self, env, grid_type):
        run = setup(env, grid_type)
        written = [run.star_runs_in(grid_type)]
        written += [path.relative_to(run.path).as_posix()
                    for path in run.steps("star1")]
        if grid_type not in SINGLE:
            written.append("binary/inlist2")
            project = run.inlist("binary/inlist_project")
            assert project["binary_controls"]["m_mesa_bc"] == "1"
            assert project["binary_job"]["m_mesa_bj"] == "1"
        for rel in written:
            values = run.inlist(rel)
            assert values["controls"]["m_mesa_c"] == "1", rel
            assert values["star_job"]["m_mesa_j"] == "1", rel
            # MESA's x_ctrl(1:num_x_ctrls) cannot be read back by MESA
            assert "x_ctrl(1:num_x_ctrls)" not in values["controls"], rel

    def test_formation_chain(self, env, grid_type):
        run = setup(env, grid_type)
        steps = run.steps("star1")
        for index, path in enumerate(steps):
            job = parse_inlist(path.read_text())["star_job"]
            assert job["save_model_when_terminate"] == ".true."
            assert job["save_model_filename"] == \
                "'initial_star1_step{}.mod'".format(index)
            if index > 0:
                assert job["load_saved_model"] == ".true."
                assert job["saved_model_name"] == \
                    "'initial_star1_step{}.mod'".format(index - 1)
            if grid_type == "CO-HeMS" and index < 2:
                assert parse_inlist(path.read_text())["controls"][
                    "m_CO-HeMS_step{}".format(index + 1)] == "1"
        if STARTS[grid_type][0] == "formation" and grid_type not in SINGLE:
            job = run.inlist("binary/inlist1")["star_job"]
            assert job["load_saved_model"] == ".true."
            assert job["saved_model_name"] == \
                "'initial_star1_step{}.mod'".format(len(steps) - 1)

    def test_record(self, env, grid_type):
        run = setup(env, grid_type)
        record = configparser.ConfigParser(interpolation=None)
        record.read(run.path / "setup_record.ini")
        assert record["setup"]["grid_type"] == grid_type
        assert record["setup"]["mode"] == "fixed"
        assert record["scenario"]["release"] == "v2.0"

    def test_column_lists_and_extras_from_the_release(self, env, grid_type):
        run = setup(env, grid_type)
        assert run.text("column_lists/history_columns.list") == \
            "release history_columns\n"
        assert run.text("binary/src/run_binary_extras.f") == \
            "! release run_binary_extras.f\n"
        for rel in ("binary/src/run_star_extras.f",
                    "star1/src/run_star_extras.f",
                    "star2/src/run_star_extras.f"):
            assert run.text(rel) == "! release run_star_extras.f\n", rel
        assert run.text("mk").count("make -f") == 3


###############################################################################
# ZAMS models
###############################################################################
@pytest.mark.parametrize("grid_type", GRID_TYPES)
@pytest.mark.parametrize("zams_key", ["zams", "zams_star1", "zams_star2"])
def test_zams(env, grid_type, zams_key):
    run = setup(env, grid_type, grid={zams_key: env["zams"]})
    zams = "'{}'".format(env["zams"])
    given = {"zams": ("star1", "star2"), "zams_star1": ("star1",),
             "zams_star2": ("star2",)}[zams_key]
    for star, start in zip(("star1", "star2"), STARTS[grid_type]):
        if start == "none":
            continue
        rel = run.star_runs_in(grid_type, star)
        value = run.inlist(rel)["controls"].get("zams_filename")
        if start == "zams" and star in given:
            assert value == zams, (star, rel)
        elif start == "formation":
            # a formed star never starts from a ZAMS model
            assert value != zams, (star, rel)
            for path in run.steps(star):
                assert parse_inlist(path.read_text())["controls"].get(
                    "zams_filename") != zams


###############################################################################
# user layers: where each key goes
###############################################################################
def reach(run, grid_type, key):
    """The written inlists (and sections) whose values include ``key``."""
    found = set()
    rels = ["binary/inlist_project", "binary/inlist1", "binary/inlist2"]
    rels += [path.relative_to(run.path).as_posix()
             for path in run.steps("star1")]
    for rel in rels:
        for section, values in run.inlist(rel).items():
            if key in values:
                found.add((rel, section))
    return found


def step_sections(run, sections=("controls", "star_job")):
    return {(path.relative_to(run.path).as_posix(), section)
            for path in run.steps("star1") for section in sections}


@pytest.mark.parametrize("grid_type", GRID_TYPES)
class TestUserInlists:
    def test_binary(self, env, grid_type):
        mine = user_file(env, "project", {"fm": "7d-3", "u_bc": "1"},
                         {"u_bj": "1"}, binary=True)
        run = setup(env, grid_type,
                     sections={"user_inlists": {"binary": mine}})
        if grid_type in SINGLE:
            return  # the binary inlists are not used
        assert reach(run, grid_type, "u_bc") == {
            ("binary/inlist_project", "binary_controls")}
        assert reach(run, grid_type, "u_bj") == {
            ("binary/inlist_project", "binary_job")}
        assert run.inlist("binary/inlist_project")["binary_controls"][
            "fm"] == "7d-3"

    def test_star1(self, env, grid_type):
        mine = user_file(env, "star1", {"mesh_delta_coeff": "0.5d0",
                                        "u_c": "1"}, {"u_j": "1"})
        run = setup(env, grid_type,
                     sections={"user_inlists": {"star1": mine}})
        if grid_type in SINGLE:
            # a single star runs as its steps: the file applies to each
            assert reach(run, grid_type, "u_c") == step_sections(
                run, ("controls",))
            assert reach(run, grid_type, "u_j") == step_sections(
                run, ("star_job",))
        else:
            assert reach(run, grid_type, "u_c") == {
                ("binary/inlist1", "controls")}
            assert reach(run, grid_type, "u_j") == {
                ("binary/inlist1", "star_job")}
        assert run.inlist(run.star_runs_in(grid_type))["controls"][
            "mesh_delta_coeff"] == "0.5d0"

    def test_star1_controls_only(self, env, grid_type):
        mine = user_file(env, "c1", {"u_c": "1"}, {"u_j": "1"})
        run = setup(env, grid_type,
                     sections={"user_inlists": {"star1_controls": mine}})
        assert reach(run, grid_type, "u_j") == set()
        assert {section for _, section in
                reach(run, grid_type, "u_c")} == {"controls"}

    def test_star2(self, env, grid_type):
        mine = user_file(env, "star2", {"u_c2": "1"}, {"u_j2": "1"})
        run = setup(env, grid_type,
                     sections={"user_inlists": {"star2": mine}})
        if grid_type in SINGLE:
            return
        assert reach(run, grid_type, "u_c2") == {
            ("binary/inlist2", "controls")}
        assert reach(run, grid_type, "u_j2") == {
            ("binary/inlist2", "star_job")}

    def test_formation_step(self, env, grid_type):
        if STAR1_STEPS[grid_type] == 0:
            pytest.skip("star 1 is not built in steps")
        mine = user_file(env, "step0", {"u_s0": "1"}, {"u_s0j": "1"})
        run = setup(env, grid_type, sections={
            "user_inlists": {"star1_formation_step0": mine}})
        assert reach(run, grid_type, "u_s0") == {
            ("star1/inlist_step0", "controls")}
        assert reach(run, grid_type, "u_s0j") == {
            ("star1/inlist_step0", "star_job")}

    def test_posydon_inlists_between_release_and_user(self, env, grid_type):
        """[posydon_inlists] adds to the release; [user_inlists] wins."""
        extra = user_file(env, "extra", {"mesh_delta_coeff": "0.6d0",
                                         "p_c": "1"})
        mine = user_file(env, "mine", {"mesh_delta_coeff": "0.5d0"})
        run = setup(env, grid_type, sections={
            "posydon_inlists": {"star1_controls": extra},
            "user_inlists": {"star1_controls": mine}})
        controls = run.inlist(run.star_runs_in(grid_type))["controls"]
        assert controls["p_c"] == "1"
        assert controls["mesh_delta_coeff"] == "0.5d0"

    def test_mesa_inlists_one_target(self, env, grid_type):
        """Naming one target's MESA defaults keeps the others'."""
        mine = write(env["tmp"] / "user" / "bj.defaults", "m_mine = 1\n")
        run = setup(env, grid_type,
                     sections={"mesa_inlists": {"binary_job": mine}})
        assert run.inlist(run.star_runs_in(grid_type))["controls"][
            "m_mesa_c"] == "1"
        if grid_type not in SINGLE:
            job = run.inlist("binary/inlist_project")["binary_job"]
            assert job["m_mine"] == "1"
            assert "m_mesa_bj" not in job

    def test_column_list(self, env, grid_type):
        mine = write(env["tmp"] / "user" / "h.list", "mine\n")
        run = setup(env, grid_type,
                     sections={"user_inlists": {"history_columns": mine}})
        assert run.text("column_lists/history_columns.list") == "mine\n"

    def test_run_star_extras(self, env, grid_type):
        mine = write(env["tmp"] / "user" / "my_extras.f", "! mine\n")
        run = setup(env, grid_type,
                     sections={"user_inlists": {"run_star_extras": mine}})
        for rel in ("binary/src/run_star_extras.f",
                    "star1/src/run_star_extras.f",
                    "star2/src/run_star_extras.f"):
            assert run.text(rel) == "! mine\n", rel
        assert "binary/src/my_extras.f" not in run.files()

    def test_run_binary_extras(self, env, grid_type):
        mine = write(env["tmp"] / "user" / "b.f", "! mine\n")
        run = setup(env, grid_type,
                     sections={"user_inlists": {"run_binary_extras": mine}})
        assert run.text("binary/src/run_binary_extras.f") == "! mine\n"

    def test_rerun_overlay(self, env, grid_type):
        run = setup(env, grid_type, scenario={"rerun": "hot"})
        controls = run.inlist(run.star_runs_in(grid_type))["controls"]
        assert controls["m_rerun"] == "1"
        assert controls["mesh_delta_coeff"] == "0.7d0"


###############################################################################
# output options and wiring
###############################################################################
# single-star grids run star 1 in its steps, but the [output] options are
# written only into binary/inlist1, which they never run; the old
# posydon-setup-grid does the same (golden single_*__scenario.txt)
SINGLE_OUTPUT_XFAIL = pytest.mark.xfail(
    strict=True, reason="[output] does not reach a single star's steps "
    "(legacy behaviour, kept)")


@pytest.mark.parametrize("grid_type", [
    pytest.param(name, marks=SINGLE_OUTPUT_XFAIL) if name in SINGLE else name
    for name in GRID_TYPES])
@pytest.mark.parametrize("star1", ["none", "history",
                                   "history, final_profile, final_model"])
def test_output_options(env, grid_type, star1):
    run = setup(env, grid_type, sections={"output": {
        "star1": star1, "binary": "none", "history_interval": 5}})
    rel = run.star_runs_in(grid_type)
    values = run.inlist(rel)
    assert values["controls"]["do_history_file"] == \
        (".true." if "history" in star1 else ".false."), rel
    assert values["star_job"]["write_profile_when_terminate"] == \
        (".true." if "final_profile" in star1 else ".false."), rel
    if "final_model" in star1:
        assert values["star_job"]["save_model_when_terminate"] == ".true."
        assert values["star_job"]["save_model_filename"] == \
            "'final_star1.mod'", rel
    assert values["controls"]["history_interval"] == "5", rel
    if grid_type not in SINGLE:
        assert run.inlist("binary/inlist_project")["binary_controls"][
            "history_interval"] == "-1"


@pytest.mark.parametrize("grid_type", GRID_TYPES)
def test_wiring(env, grid_type):
    """initial_mass is a grid column and a star control."""
    run = setup(env, grid_type)
    rel = run.star_runs_in(grid_type)
    controls = run.inlist(rel)["controls"]
    if grid_type in SINGLE:
        # posydon-run-grid wraps each step in its own inlist that reads the
        # step and then inlist_grid_points; the step itself is not wired
        assert "read_extra_controls_inlist1" not in controls
        return
    assert controls["read_extra_controls_inlist1"] == ".true.", rel
    assert controls["extra_controls_inlist1_name"] == \
        "'inlist_grid_star1_binary_controls'"
    job = run.inlist("binary/inlist_project")["binary_job"]
    assert job["inlist_names(1)"] == "'{}'".format(
        run.path / "binary" / "inlist1")


###############################################################################
# modes and backends
###############################################################################
@pytest.mark.parametrize("grid_type", ["HMS-HMS", "CO-HeMS", "single_HMS"])
@pytest.mark.parametrize("mode, submission, scripts", [
    ("fixed", {"backend": "shell", "job_array": "true"},
     {"grid_command.sh"}),
    ("fixed", {"backend": "slurm_array"},
     {"job_array_grid_submit.slurm", "cleanup.slurm", "run_grid.sh"}),
    ("dynamic", {"backend": "slurm_mpi", "nodes": 2, "mpi_tasks": 4},
     {"mpi_grid_submit.slurm", "cleanup.slurm", "run_grid.sh"}),
    ("dynamic", {"backend": "shell", "job_array": "false"},
     {"grid_command.sh"}),
])
def test_backends(env, grid_type, mode, submission, scripts):
    grid = {"mode": mode}
    if mode == "dynamic":
        grid["psycris_inifile"] = env["tmp"] / "psycris.ini"
    run = setup(env, grid_type, grid=grid, submission=submission)
    written = {name for name in run.files() if name.endswith(
        (".sh", ".slurm"))}
    assert written == scripts
    for name in scripts:
        assert "None\n" not in run.text(name), name
        assert "=None" not in run.text(name), name
    main = sorted(scripts - {"cleanup.slurm", "run_grid.sh"})[0]
    text = run.text(main)
    if mode == "dynamic":
        assert "mpirun --bind-to none -np 8 " in text \
            if submission["backend"] == "slurm_mpi" else \
            "mpirun --bind-to none -np 1 " in text
        assert "--grid-type dynamic" in text
        assert "--psycris-inifile {}".format(env["tmp"] / "psycris.ini") \
            in text
        assert "--grid-point-index" not in text
    else:
        assert "--grid-type fixed" in text
        assert "--grid-point-index $SLURM_ARRAY_TASK_ID" in text
    if submission["backend"] == "slurm_array":
        assert "#SBATCH --array=0-2\n" in text
    if submission["backend"] != "shell":
        assert "--job_end $SLURM_JOB_END_TIME" in text


def test_slurm_settings(env):
    run = setup(env, "HMS-HMS", grid={"keep_profiles": "true",
                                      "keep_photos": "true"},
                submission={"backend": "slurm_array", "account": "acc",
                            "partition": "part", "walltime": "1:00:00",
                            "email": "a@b.c", "cpus_per_task": 3,
                            "mem_per_cpu": "8G", "work_dir": "/scratch",
                            "newgroup": "grp"})
    text = run.text("job_array_grid_submit.slurm")
    for line in ("--account=acc", "--partition=part", "--time=1:00:00",
                 "--mail-user=a@b.c", "--cpus-per-task 3",
                 "--mem-per-cpu=8G"):
        assert "#SBATCH {}\n".format(line) in text, line
    assert "export OMP_NUM_THREADS=3" in text
    assert "--temporary-directory /scratch" in text
    assert "--keep_profiles --keep_photos" in text
    assert "chgrp -fR grp ." in run.text("cleanup.slurm")


###############################################################################
# refusals
###############################################################################
@pytest.mark.parametrize("grid, submission, scenario, sections, message", [
    ({}, {"backend": "slurm_mpi"}, None, None,
     r"\[submission\] backend: a fixed grid runs one point per task"),
    ({}, {"backend": "shell", "job_array": "false"}, None, None,
     r"\[submission\] backend: a fixed grid"),
    ({"mode": "dynamic", "psycris_inifile": "psycris.ini"},
     {"backend": "slurm_array"}, None, None,
     r"\[submission\] backend: a dynamic grid runs as one MPI job"),
    ({"mode": "dynamic"}, {"backend": "slurm_mpi"}, None, None,
     r"\[grid\] psycris_inifile: needed"),
    ({}, None, {"release": None, "ref": "main"}, None,
     r"\[scenario\] release: missing"),
    ({}, None, {"release": None, "rerun": "hot"}, None,
     r"\[scenario\] rerun: a rerun needs a release"),
    ({}, None, {"rerun": "cold"}, None, r"\[scenario\] rerun: no rerun"),
    ({}, None, {"release": "v9"}, None, r"\[scenario\] release: no release"),
    ({}, None, None, {"user_inlists": {"star3": "x"}},
     r"\[user_inlists\] star3: unknown key"),
    ({}, None, None, {"output": {"star1": "everything"}},
     r"\[output\] star1: unknown value"),
    ({"type": "MS-MS"}, None, None, None, r"\[grid\] type: unknown grid"),
])
def test_refused(env, grid, submission, scenario, sections, message):
    ini = write(env["tmp"] / "grid.ini", ini_text(
        env, "HMS-HMS", grid=grid, scenario=scenario, sections=sections,
        submission=submission))
    with pytest.raises(ConfigError, match=message):
        totest.setup_grid(ini)
    assert not (env["tmp"] / "run").exists()


def test_no_scenario_needs_mesa_defaults_and_columns(env):
    ini = write(env["tmp"] / "grid.ini", ini_text(
        env, "HMS-HMS", scenario={"release": None, "source": None,
                                  "inlist_repo": None}))
    with pytest.raises(ConfigError, match=r"\[mesa_inlists\]: needed"):
        totest.setup_grid(ini)


###############################################################################
# explain
###############################################################################
@pytest.mark.parametrize("grid_type", GRID_TYPES)
def test_explain_names_the_user_layer(env, grid_type):
    mine = user_file(env, "star1", {"mesh_delta_coeff": "0.5d0"})
    ini = write(env["tmp"] / "grid.ini", ini_text(
        env, grid_type, sections={"user_inlists": {"star1": mine}}))
    text = totest.explain_grid(ini)
    assert "mesh_delta_coeff = 0.5d0  [user_inlists] " in text
    assert not (env["tmp"] / "run").exists()


###############################################################################
# every layer filled at once
###############################################################################
LAYER_SETS = {
    "mesa": ("mesa_inlists",),
    "mesa+posydon": ("mesa_inlists", "posydon_inlists"),
    "all": ("mesa_inlists", "posydon_inlists", "user_inlists"),
}
SHORT = {"mesa_inlists": "mesa", "posydon_inlists": "posydon",
         "user_inlists": "user"}


def layer_files(env, layer, grid_type, style):
    """One file per target of a layer; each sets lv and a marker."""
    short = SHORT[layer]
    folder = env["tmp"] / layer

    def make(name, binary=False):
        values = {"lv": "'{}_{}'".format(short, name),
                  "{}_m_{}".format(short, name): "1"}
        if layer == "mesa_inlists":
            # MESA's .defaults format: no namelist lines
            return write(folder / name, "".join(
                "{} = {}\n".format(key, value)
                for key, value in values.items()))
        return write(folder / name, inlist(values, values, binary=binary))

    keys = {}
    if grid_type not in SINGLE:
        keys["binary"] = make("binary", binary=True)
        keys["star2"] = make("star2")
    keys["star1"] = make("star1")
    steps = [make("step{}".format(step))
             for step in range(STAR1_STEPS[grid_type])]
    if steps and style == "list":
        keys["star1_formation"] = ", ".join(str(path) for path in steps)
    else:
        for step, path in enumerate(steps):
            keys["star1_formation_step{}".format(step)] = path
    for name in ("history_columns", "binary_history_columns",
                 "profile_columns"):
        keys[name] = write(folder / (name + ".list"), short + "\n")
    keys["run_star_extras"] = write(folder / "rse.f", "! " + short + "\n")
    keys["run_binary_extras"] = write(folder / "rbe.f", "! " + short + "\n")
    return keys


def written_targets(run, grid_type):
    """Written inlist -> the target name its most specific files use."""
    out = {path.relative_to(run.path).as_posix(): "step{}".format(index)
           for index, path in enumerate(run.steps("star1"))}
    if grid_type not in SINGLE:
        out.update({"binary/inlist_project": "binary",
                    "binary/inlist1": "star1", "binary/inlist2": "star2"})
    return out


@pytest.mark.parametrize("grid_type", GRID_TYPES)
@pytest.mark.parametrize("layers", list(LAYER_SETS))
@pytest.mark.parametrize("rerun", [False, True])
@pytest.mark.parametrize("style", ["steps", "list"])
def test_every_layer_filled(env, grid_type, layers, rerun, style):
    """mesa < release < rerun < posydon < user, for every target at once."""
    given = LAYER_SETS[layers]
    sections = {layer: layer_files(env, layer, grid_type, style)
                for layer in given}
    run = setup(env, grid_type, sections=sections,
                scenario={"rerun": "hot"} if rerun else None)
    for rel, target in written_targets(run, grid_type).items():
        values = run.inlist(rel)
        star1_runs_here = rel.startswith("star1/") and grid_type in SINGLE \
            or rel == "binary/inlist1"
        # every section of the file: the highest layer wins
        if "user_inlists" in given:
            top = "'user_{}'".format(target)
        elif "posydon_inlists" in given and target.startswith("step"):
            # the scenario's files come after the ini's [posydon_inlists]
            # (as for column lists, C7); only a step's own file is more
            # specific than the release's every-step file
            top = "'posydon_{}'".format(target)
        elif rerun and rel != "binary/inlist_project" \
                and rel != "binary/inlist2":
            top = "'rerun'"
        else:
            top = "'release'"
        for section, keys in values.items():
            assert keys["lv"] == top, (rel, section)
            for layer in given:
                short = SHORT[layer]
                # the layer's own file for this target...
                assert keys.get("{}_m_{}".format(short, target)) == "1", \
                    (rel, section, layer)
                # ...and star 1's on every step of a single star only
                has_star1 = keys.get("{}_m_star1".format(short)) == "1"
                assert has_star1 == star1_runs_here, (rel, section, layer)
            # an ini file replaces MESA's defaults for its target
            assert "m_mesa_c" not in keys and "m_mesa_j" not in keys and \
                "m_mesa_bc" not in keys and "m_mesa_bj" not in keys, rel
    top_whole = SHORT[given[-1]] if "user_inlists" in given else "release"
    for name in ("history_columns", "binary_history_columns",
                 "profile_columns"):
        expected = "release {}\n".format(name) if top_whole == "release" \
            else top_whole + "\n"
        assert run.text("column_lists/{}.list".format(name)) == expected
    # Fortran sources follow the same rule as the column lists
    for rel in ("binary/src/run_star_extras.f",
                "star1/src/run_star_extras.f",
                "star2/src/run_star_extras.f"):
        expected = "! release run_star_extras.f\n" \
            if top_whole == "release" else "! user\n"
        assert run.text(rel) == expected, rel
    assert run.text("binary/src/run_binary_extras.f") == (
        "! release run_binary_extras.f\n" if top_whole == "release"
        else "! user\n")
