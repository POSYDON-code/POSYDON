"""Unit tests of posydon/grids/setup/build.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

# import other needed code for the tests, which is not already imported in the
# module you like to test
from dataclasses import is_dataclass
from inspect import isroutine
from pathlib import Path

import pytest

# import the module which will be tested
import posydon.grids.setup.build as totest
from posydon.grids.setup.config import GridConfig, OutputSection
from posydon.grids.setup.inlist import parse_inlist
from posydon.grids.setup.scenario import ScenarioFiles


@pytest.fixture
def files(tmp_path):
    """Small stand-ins for MESA defaults and POSYDON inlists."""
    texts = {
        "controls.defaults": "mesh_delta_coeff = 1d0\n"
                             "x_ctrl(1:num_x_ctrls) = 0d0\n"
                             "initial_mass = 1d0\nzams_filename = ''\n",
        "star_job.defaults": "pgstar_flag = .false.\n"
                             "save_model_when_terminate = .false.\n",
        "binary_controls.defaults": "fm = 1d-2\nm1 = 1d0\n",
        "binary_job.defaults": "inlist_names(1) = 'inlist1'\n"
                               "evolve_both_stars = .true.\n",
        "inlist_project": "&binary_controls\n fm = 5d-3\n/\n"
                          "&binary_job\n evolve_both_stars = .true.\n/\n",
        "inlist1": "&controls\n mesh_delta_coeff = 0.8d0\n"
                   " read_extra_controls_inlist1 = .true.\n/\n"
                   "&star_job\n pgstar_flag = .false.\n/\n",
        "step1": "&controls\n relax_initial_to_xaccrete = .true.\n/\n"
                 "&star_job\n/\n",
        "step2": "&controls\n relax_initial_Z = .true.\n/\n&star_job\n/\n",
    }
    for name, text in texts.items():
        (tmp_path / name).write_text(text)
    return tmp_path


def binary_config(files, single=False):
    """A two-star grid with MESA defaults and POSYDON inlists."""
    config = GridConfig()
    mesa = config.inlists["mesa_inlists"]
    posydon = config.inlists["posydon_inlists"]
    for target, name in (("binary_controls", "binary_controls.defaults"),
                         ("binary_job", "binary_job.defaults"),
                         ("star1_controls", "controls.defaults"),
                         ("star1_job", "star_job.defaults"),
                         ("star2_controls", "controls.defaults"),
                         ("star2_job", "star_job.defaults")):
        mesa.add(target, files / name)
    for target in ("binary_controls", "binary_job"):
        posydon.add(target, files / "inlist_project")
    for target in ("star1_controls", "star1_job", "star2_controls",
                   "star2_job"):
        posydon.add(target, files / "inlist1")
    config.legacy = {"mesa_inlists": {"single_star_grid": single,
                                      "history_interval": 1}}
    return config


class TestElements:
    def test_dir(self):
        elements = {'BuiltInlists', 'LAYER_ORDER', 'LayerStack', 'Path',
                    'STARS', '__authors__', '__builtins__', '__cached__',
                    '__doc__', '__file__', '__loader__', '__name__',
                    '__package__', '__spec__', '_cleanup_for', '_lines',
                    '_merge_dicts', 'build_inlists', 'dataclass',
                    'drop_wiring_keys', 'expand_every_step',
                    'fix_num_x_ctrls', 'formation_step_count',
                    'get_grid_type', 'is_single_star', 'layer_targets',
                    'legacy_zams', 'os', 'render_binary_project',
                    'render_binary_star', 'render_formation_step',
                    'stack_target'}
        assert set(dir(totest)) == elements

    def test_instance_BuiltInlists(self):
        assert is_dataclass(totest.BuiltInlists)

    @pytest.mark.parametrize("name", [
        'is_single_star', 'layer_targets', 'stack_target', 'legacy_zams',
        'render_binary_project', 'render_binary_star',
        'render_formation_step', 'build_inlists'])
    def test_instance_functions(self, name):
        assert isroutine(getattr(totest, name))


class TestIsSingleStar:
    def test_from_grid_type(self):
        config = GridConfig()
        config.grid.type = "single_HeMS"
        assert totest.is_single_star(config) is True
        config.grid.type = "CO-HeMS"
        assert totest.is_single_star(config) is False

    def test_from_legacy_ini(self):
        config = GridConfig()
        config.legacy = {"mesa_inlists": {"single_star_grid": True}}
        assert totest.is_single_star(config) is True

    def test_default(self):
        assert totest.is_single_star(GridConfig()) is False


class TestLayerTargets:
    def test_scenario_files_and_flags_are_added(self, files):
        config = GridConfig()
        config.inlists["posydon_inlists"].add("star1_controls", "ini")
        config.inlists["mesa_inlists"].add(
            "star1_formation_every_step_controls", "defaults")
        scenario = ScenarioFiles()
        scenario.inlists.add("star1_controls", "scenario")
        scenario.inlists.add("star1_formation_step1_controls", "step2")
        scenario.flags = {"star1_formation_every_step_controls": {"a": "1"}}
        layers, flags = totest.layer_targets(config, scenario)
        assert layers["posydon_inlists"]["star1_controls"] == \
            ["ini", "scenario"]
        assert layers["mesa_inlists"] == {
            "star1_formation_step0_controls": ["defaults"],
            "star1_formation_step1_controls": ["defaults"]}
        assert flags == {"star1_formation_step0_controls": {"a": "1"},
                         "star1_formation_step1_controls": {"a": "1"}}

    def test_without_scenario(self):
        config = GridConfig()
        config.inlists["user_inlists"].add("binary_job", "mine")
        layers, flags = totest.layer_targets(config)
        assert layers["user_inlists"] == {"binary_job": ["mine"]}
        assert flags == {}


class TestStackTarget:
    def test_cleanup_and_flags(self, files):
        layers = {"mesa_inlists": {"star1_controls":
                                   [files / "controls.defaults"]},
                  "posydon_inlists": {"star1_controls":
                                      [files / "inlist1"]}}
        flags = {"star1_controls": {"x_logical_ctrl(1)": ".true."}}
        merged = totest.stack_target("star1_controls", layers,
                                     flags).merged()
        assert merged == {"mesh_delta_coeff": "0.8d0",
                          "x_ctrl(1:1)": "0d0", "initial_mass": "1d0",
                          "x_logical_ctrl(1)": ".true."}

    def test_binary_keeps_num_x_ctrls(self, tmp_path):
        path = tmp_path / "p"
        path.write_text("&binary_controls\n x(1:num_x_ctrls) = 0\n"
                        " read_extra_binary_controls_inlist1 = .true.\n/\n")
        merged = totest.stack_target(
            "binary_controls", {"user_inlists": {"binary_controls": [path]}},
            {}).merged()
        assert merged == {"x(1:num_x_ctrls)": "0"}


class TestLegacyZams:
    def test_ini_keys(self):
        config = GridConfig()
        config.legacy = {"mesa_inlists": {"zams_filename_2": "z2"}}
        assert totest.legacy_zams(config) == {"star1": None, "star2": "z2"}

    def test_formed_star_has_no_zams_whatever_zams_filename_says(self):
        """S2: a HeMS-HMS grid with zams_filename still forms star 1."""
        config = GridConfig()
        config.legacy = {"mesa_inlists": {"zams_filename": "z"}}
        scenario = ScenarioFiles(formation_stars=("star1",))
        assert totest.legacy_zams(config, scenario) == {"star1": None,
                                                        "star2": "z"}

    def test_scenario(self):
        config = GridConfig()
        config.legacy = {"mesa_inlists": {"zams_filename_1": "ini1",
                                          "zams_filename_2": "ini2"}}
        scenario = ScenarioFiles(zams_star1=Path("/r/z1"),
                                 zams_star2=Path("/r/z2"),
                                 formation_stars=("star2",))
        assert totest.legacy_zams(config, scenario) == {"star1": "/r/z1",
                                                        "star2": None}

    def test_formed_single_star_has_no_zams(self):
        """S3: the scenario forms star 1 of single_HeMS."""
        config = GridConfig()
        config.grid.type = "single_HeMS"
        config.legacy = {"mesa_inlists": {"zams_filename": "z"}}
        scenario = ScenarioFiles(formation_stars=("star1",))
        assert totest.legacy_zams(config, scenario) == {"star1": None,
                                                        "star2": None}

    def test_single_star_clears_star2(self):
        config = GridConfig()
        config.grid.type = "single_HeMS"
        config.legacy = {"mesa_inlists": {"zams_filename_2": "z2"}}
        assert totest.legacy_zams(config, ScenarioFiles())["star2"] is None

    def test_new_format_single_star_keeps_star2(self):
        config = GridConfig()
        config.grid.type = "single_HMS"
        config.grid.zams_star1 = config.grid.zams_star2 = Path("/z.data")
        assert totest.legacy_zams(config, ScenarioFiles()) == {
            "star1": "/z.data", "star2": "/z.data"}

    def test_new_format(self):
        config = GridConfig()
        config.grid.zams_star1 = Path("/z/one.data")
        assert totest.legacy_zams(config) == {"star1": "/z/one.data",
                                              "star2": None}


class TestRender:
    def test_binary_project(self):
        assert totest.render_binary_project({"a": "1"}, {"b": "2"}) == (
            "&binary_controls\n\n\ta = 1\n"
            "\n/ ! end of binary_controls namelist\n\n"
            "&binary_job\n\n\tb = 2\n"
            "\n/ ! end of binary_job namelist\n                ")

    def test_binary_star(self):
        assert totest.render_binary_star("star2", {"a": "1"}, {}) == (
            "&controls\n\n\ta = 1\n\n\n"
            "\n/ ! end of star2_controls namelist\n\n"
            "&star_job\n\n"
            "\n/ ! end of star2_job namelist\n\n")

    def test_formation_step(self):
        assert totest.render_formation_step("star1", {}, {"b": "2"}) == (
            "&controls\n\n\n\n"
            "\n        / ! end of star1_controls namelist\n\n        "
            "&star_job\n\n\tb = 2\n"
            "\n        / ! end of star1_job namelist\n\n        ")


def read(path):
    return parse_inlist(Path(path).read_text())


class TestBuildInlists:
    def test_binary_grid(self, files, tmp_path):
        config = binary_config(files)
        config.legacy["mesa_inlists"]["zams_filename"] = "/z/zams.data"
        run = tmp_path / "run"
        built = totest.build_inlists(config, run, ["m1", "mesh_delta_coeff",
                                                   "pgstar_flag"])
        assert built == totest.BuiltInlists(
            binary_project=str(run / "binary" / "inlist_project"),
            binary_star1=str(run / "binary" / "inlist1"),
            binary_star2=str(run / "binary" / "inlist2"))
        project = read(built.binary_project)
        assert project["binary_controls"] == {"fm": "5d-3", "m1": "1d0",
                                              "history_interval": "1"}
        assert project["binary_job"] == {
            "evolve_both_stars": ".true.",
            "inlist_names(1)": "'{}'".format(run / "binary" / "inlist1"),
            "inlist_names(2)": "'{}'".format(run / "binary" / "inlist2")}
        star1 = read(built.binary_star1)
        assert star1["controls"]["read_extra_controls_inlist1"] == ".true."
        assert star1["controls"]["extra_controls_inlist1_name"] == \
            "'inlist_grid_star1_binary_controls'"
        assert star1["controls"]["zams_filename"] == "'/z/zams.data'"
        assert star1["controls"]["do_history_file"] == ".true."
        assert star1["star_job"]["extra_star_job_inlist1_name"] == \
            "'inlist_grid_star1_binary_job'"
        assert star1["star_job"]["save_model_filename"] == "'final_star1.mod'"
        assert star1["star_job"]["write_profile_when_terminate"] == \
            ".false."
        assert not (run / "star1").exists()

    def test_no_grid_columns_in_star_inlists(self, files, tmp_path):
        built = totest.build_inlists(binary_config(files), tmp_path, ["m1"])
        assert "read_extra_controls_inlist1" not in \
            read(built.binary_star1)["controls"]

    def test_output_options(self, files, tmp_path):
        config = binary_config(files)
        config.output = OutputSection(binary_history=False,
                                      star1=("final_profile",), star2=())
        built = totest.build_inlists(config, tmp_path, [])
        assert read(built.binary_project)["binary_controls"][
            "history_interval"] == "-1"
        star1 = read(built.binary_star1)
        assert star1["star_job"]["filename_for_profile_when_terminate"] == \
            "'final_profile_star1.data'"
        assert star1["star_job"]["save_model_when_terminate"] == ".false."
        assert star1["controls"]["do_history_file"] == ".false."

    def test_formation_steps(self, files, tmp_path):
        config = binary_config(files)
        mesa = config.inlists["mesa_inlists"]
        mesa.add("star1_formation_every_step_controls",
                 files / "controls.defaults")
        mesa.add("star1_formation_every_step_job",
                 files / "star_job.defaults")
        user = config.inlists["user_inlists"]
        for index, name in enumerate(("step1", "step2")):
            for part in ("controls", "job"):
                user.add("star1_formation_step{}_{}".format(index, part),
                         files / name)
        built = totest.build_inlists(config, tmp_path, [])
        steps = [tmp_path / "star1" / "inlist_step0",
                 tmp_path / "star1" / "inlist_step1"]
        assert built.star1_formation == " {} {}".format(*steps)
        assert built.star2_formation is None
        step0, step1 = read(steps[0]), read(steps[1])
        assert step0["controls"]["relax_initial_to_xaccrete"] == ".true."
        assert step0["star_job"]["save_model_filename"] == \
            "'initial_star1_step0.mod'"
        assert "load_saved_model" not in step0["star_job"]
        assert step1["star_job"]["saved_model_name"] == \
            "'initial_star1_step0.mod'"
        assert step1["star_job"]["save_model_filename"] == \
            "'initial_star1_step1.mod'"
        job = read(built.binary_star1)["star_job"]
        assert job["load_saved_model"] == ".true."
        assert job["saved_model_name"] == "'initial_star1_step1.mod'"

    def test_zams_star_is_not_formed(self, files, tmp_path):
        """Today: a ZAMS model wins over formation steps (S2)."""
        config = binary_config(files)
        config.legacy["mesa_inlists"]["zams_filename"] = "/z.data"
        for target in ("star1_formation_every_step_controls",
                       "star2_formation_every_step_controls"):
            config.inlists["user_inlists"].add(target, files / "step1")
        built = totest.build_inlists(config, tmp_path, [])
        assert built.star1_formation is None
        assert built.star2_formation is None

    def test_star2_formation(self, files, tmp_path):
        config = binary_config(files)
        config.inlists["user_inlists"].add(
            "star2_formation_every_step_controls", files / "step1")
        built = totest.build_inlists(config, tmp_path, [])
        assert built.star2_formation == " {}".format(
            tmp_path / "star2" / "inlist_step0")
        assert read(built.binary_star2)["star_job"]["saved_model_name"] == \
            "'initial_star2_step0.mod'"

    def test_single_star(self, files, tmp_path):
        config = binary_config(files, single=True)
        config.legacy["mesa_inlists"]["zams_filename"] = "/z.data"
        for part, name in (("controls", "controls.defaults"),
                           ("job", "star_job.defaults")):
            config.inlists["mesa_inlists"].add(
                "star1_formation_every_step_{}".format(part), files / name)
        scenario = ScenarioFiles()
        scenario.flags = {"star1_formation_every_step_controls":
                          {"x_logical_ctrl(1)": ".true."}}
        built = totest.build_inlists(config, tmp_path, ["initial_mass"],
                                     scenario=scenario)
        step0 = read(tmp_path / "star1" / "inlist_step0")
        assert step0["controls"]["x_logical_ctrl(1)"] == ".true."
        assert step0["controls"]["zams_filename"] == "'/z.data'"

        # without a ZAMS model, not even MESA's default is kept (S3)
        del config.legacy["mesa_inlists"]["zams_filename"]
        totest.build_inlists(config, tmp_path, ["initial_mass"],
                             scenario=scenario)
        step0 = read(tmp_path / "star1" / "inlist_step0")
        assert "zams_filename" not in step0["controls"]
        assert read(built.binary_project) == {
            "binary_controls": {"history_interval": "1"}, "binary_job": {}}
        assert "read_extra_controls_inlist1" not in \
            read(built.binary_star1)["controls"]
