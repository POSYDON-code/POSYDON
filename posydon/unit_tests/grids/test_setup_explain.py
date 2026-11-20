"""Unit tests of posydon/grids/setup/explain.py and template.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

# import other needed code for the tests, which is not already imported in the
# module you like to test
import configparser
import doctest
import re
from dataclasses import fields
from inspect import isroutine

import pytest

# import the module which will be tested
import posydon.grids.setup.explain as totest
import posydon.grids.setup.template as template
from posydon.grids.setup.config import (
    COLUMN_FILES,
    EXTRAS_KEYS,
    INLIST_SHORTHANDS,
    INLIST_TARGETS,
    GridConfig,
    GridSection,
    OutputSection,
    ScenarioSection,
    SubmissionSection,
    read_config,
)
from posydon.grids.setup.grid_types import GRID_TYPES
from posydon.grids.setup.scenario import ScenarioFiles


def write(path, text):
    path.write_text(text)
    return path


@pytest.fixture
def config(tmp_path):
    defaults = write(tmp_path / "controls.defaults",
                     "mesh_delta_coeff = 1d0\nalpha_semiconvection = 0d0\n")
    common = write(tmp_path / "inlist1",
                   "&controls\n  mesh_delta_coeff = 0.9d0\n/\n")
    mine = write(tmp_path / "mine",
                 "&controls\n  mesh_delta_coeff = 0.8d0\n/\n")
    config = GridConfig()
    config.grid.type = "HMS-HMS"
    config.inlists["mesa_inlists"].add("star1_controls", defaults)
    config.inlists["user_inlists"].add("star1_controls", mine)
    scenario = ScenarioFiles()
    scenario.inlists.add("star1_controls", common)
    return config, scenario, {"defaults": defaults, "common": common,
                              "mine": mine}


class TestElements:
    def test_dir(self):
        assert {'BINARY_TARGETS', 'Path', '_short', 'explain_text',
                'explained_targets',
                'is_single_star', 'layer_targets', 'stack_target'} <= \
            set(dir(totest))

    @pytest.mark.parametrize("name", ['explained_targets', 'explain_text'])
    def test_instance_functions(self, name):
        assert isroutine(getattr(totest, name))


class TestExplain:
    def test_history(self, config):
        config, scenario, files = config
        explained = totest.explained_targets(config, scenario)
        assert list(explained) == ["star1_controls"]
        (key, history), = explained["star1_controls"]
        assert key == "mesh_delta_coeff"
        assert [(e["layer"], e["value"]) for e in history] == [
            ("mesa_inlists", "1d0"), ("posydon_inlists", "0.9d0"),
            ("user_inlists", "0.8d0")]

    def test_text(self, config):
        config, scenario, files = config
        text = totest.explain_text(config, scenario)
        assert text.splitlines()[:4] == [
            "star1_controls",
            "  mesh_delta_coeff = 0.8d0  [user_inlists] {}".format(
                files["mine"]),
            "      over 0.9d0  [posydon_inlists] {}".format(files["common"]),
            "      over 1d0  [mesa_inlists] {}".format(files["defaults"])]

    def test_short_paths_and_same_values(self, config, tmp_path):
        config, scenario, files = config
        files["mine"].write_text("&controls\n  mesh_delta_coeff = 0.9d0\n/\n")
        text = totest.explain_text(config, scenario,
                                   roots=[tmp_path / "elsewhere", tmp_path])
        assert text.splitlines()[1] == ("  mesh_delta_coeff = 0.9d0  "
                                        "[user_inlists] mine  (same value "
                                        "as below)")
        assert text.splitlines()[2] == ("      over 0.9d0  [posydon_inlists]"
                                        " inlist1")

    def test_defaults_only_and_order(self, config):
        config, scenario, files = config
        for target in ("star2_controls", "binary_controls",
                       "star1_formation_step0_controls",
                       "star1_formation_every_step_controls"):
            config.inlists["mesa_inlists"].add(target, files["defaults"])
        explained = totest.explained_targets(config)
        assert list(explained) == ["binary_controls", "star1_controls",
                                   "star2_controls",
                                   "star1_formation_step0_controls"]
        assert "(MESA defaults only)" in totest.explain_text(config)

    def test_single_star_skips_the_binary(self, config):
        config, scenario, files = config
        config.grid.type = "single_HMS"
        config.inlists["user_inlists"].add("star1_formation_step0_controls",
                                           files["mine"])
        assert list(totest.explained_targets(config, scenario)) == [
            "star1_formation_step0_controls"]


class TestTemplate:
    def test_doctests(self):
        assert doctest.testmod(template).failed == 0

    @pytest.mark.parametrize("name", list(GRID_TYPES))
    def test_parses(self, name, tmp_path):
        ini = write(tmp_path / "grid.ini", template.ini_template(name))
        config = read_config(ini)
        assert config.grid.type == name
        assert config.scenario.release == "v2.0"
        assert config.scenario.inlist_repo is None
        assert config.submission.account == "<account>"

    @pytest.mark.parametrize("name", list(GRID_TYPES))
    def test_every_option_is_valid(self, name, tmp_path):
        """Uncommenting every option of the template still reads."""
        text = template.ini_template(name)
        options = re.findall(r"^; ([a-z0-9_]+) = ", text, flags=re.M)
        assert len(options) > 20
        text = re.sub(r"^; ([a-z0-9_]+ = )", r"\1", text, flags=re.M)
        text = re.sub(r"<[^>]*>", "x", text)
        config = read_config(write(tmp_path / "grid.ini", text))
        assert config.scenario.inlist_repo is not None
        assert config.scenario.rerun is not None
        assert config.submission.newgroup == "x"
        assert set(config.extras) == set(EXTRAS_KEYS)

    @pytest.mark.parametrize("name", list(GRID_TYPES))
    def test_lists_every_option(self, name):
        """Each option the new format has for the grid type is in it."""
        grid_type = GRID_TYPES[name]
        text = template.ini_template(name)
        parser = configparser.ConfigParser(interpolation=None,
                                           inline_comment_prefixes=(";",))
        parser.optionxform = str
        parser.read_string(re.sub(r"^; ([a-z0-9_]+ = )", r"\1", text,
                                  flags=re.M))
        listed = {section: set(parser[section])
                  for section in parser.sections()}
        zams = {star for star in ("star1", "star2")
                if getattr(grid_type, star) == "zams"}
        grid = {field.name for field in fields(GridSection)} - {
            "zams_star1", "zams_star2"}
        grid |= {"zams_" + star for star in zams}
        if len(zams) == 2:
            grid.add("zams")
        assert listed["grid"] == grid
        # source, branch, commit and system_type come from legacy inis
        assert listed["scenario"] == {
            field.name for field in fields(ScenarioSection)} - {
            "source", "branch", "commit", "system_type"}
        stars = {"star1"}
        if grid_type.star2 in ("zams", "formation"):
            stars.add("star2")
        # Fortran sources: every star; the binary and star 1 of a binary
        # grid; star 2 when it is evolved
        sources = {"run_star_extras"}
        if not grid_type.single_star:
            sources |= {"run_binary_extras", "run_star1_extras"}
        if "star2" in stars:
            sources.add("run_star2_extras")
        inlists = set(COLUMN_FILES) | sources | {
            key for key in set(INLIST_TARGETS) | set(INLIST_SHORTHANDS)
            if key.split("_")[0] in stars or
            (key.startswith("binary") and not grid_type.single_star)}
        for star in stars:
            if getattr(grid_type, star) == "formation":
                inlists |= {star + "_formation", star + "_formation_step0"}
        assert listed["user_inlists"] == inlists
        output = {"history_interval", "star1"}
        if not grid_type.single_star:
            output |= {"star2", "binary"}
        assert listed["output"] <= {"history_interval", "star1", "star2",
                                    "binary"}
        assert len(fields(OutputSection)) == 4
        assert listed["output"] == output
        assert listed["submission"] == {
            field.name for field in fields(SubmissionSection)}
        assert listed["extras"] == set(EXTRAS_KEYS)

    def test_zams_keys(self):
        assert "\nzams = <" in template.ini_template("HMS-HMS")
        assert "\nzams_star2 = <" in template.ini_template("HeMS-HMS")
        assert "\nzams_star1 = <" in template.ini_template("single_HMS")
        assert "zams" not in template.ini_template(
            "CO-HeMS").split("[scenario]")[0]

    def test_outputs(self):
        text = template.ini_template("single_HeMS")
        assert "star1 = history, final_model" in text
        # a single star has no star 2 and no binary
        assert "\nstar2 =" not in text and "\nbinary =" not in text
        assert "\nstar2 = none" in template.ini_template("CO-HeMS")

    def test_width(self):
        for name in GRID_TYPES:
            assert max(len(line) for line in template.ini_template(
                name).splitlines()) <= template.WIDTH
