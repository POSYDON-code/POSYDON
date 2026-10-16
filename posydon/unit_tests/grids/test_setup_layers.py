"""Unit tests of posydon/grids/setup/layers.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

# import other needed code for the tests, which is not already imported in the
# module you like to test
import doctest
from inspect import isclass, isroutine
from pathlib import Path

import pytest

# import the module which will be tested
import posydon.grids.setup.layers as totest
from posydon.utils.gridutils import clean_inlist_file

MESA_DEFAULTS = Path(__file__).resolve().parents[3] / "grid_params" \
    / "defaults" / "r11701"

POSYDON_INLIST1 = """
&controls
    mesh_delta_coeff = 0.8d0
    read_extra_controls_inlist1 = .true.
    extra_controls_inlist1_name = 'inlist_extra'
    x_ctrl(1:num_x_ctrls) = 0d0
    zams_filename = ''
/
&star_job
    pgstar_flag = .false.
/
"""


POSYDON_INLIST_PROJECT = """
&binary_controls
    fm = 0.002d0
    read_extra_binary_controls_inlist1 = .true.
/
&binary_job
    inlist_names(1) = 'inlist1'
    evolve_both_stars = .true.
/
"""


@pytest.fixture
def posydon_inlist1(tmp_path):
    path = tmp_path / "inlist1"
    path.write_text(POSYDON_INLIST1)
    return path


@pytest.fixture
def posydon_inlist_project(tmp_path):
    path = tmp_path / "inlist_project"
    path.write_text(POSYDON_INLIST_PROJECT)
    return path


class TestElements:
    def test_dir(self):
        elements = {'LAYER_ORDER', 'LEGACY_SKIPPED_VALUES', 'LayerStack',
                    '__authors__', '__builtins__', '__cached__', '__doc__',
                    '__file__', '__loader__', '__name__', '__package__',
                    '__spec__', 'drop_wiring_keys', 'fix_num_x_ctrls',
                    'pick_whole_file', 'read_inlist', 'section_for'}
        assert set(dir(totest)) == elements

    def test_instance_LAYER_ORDER(self):
        assert totest.LAYER_ORDER == ("mesa_inlists", "posydon_inlists",
                                      "user_inlists", "generated")

    def test_instance_LayerStack(self):
        assert isclass(totest.LayerStack)

    @pytest.mark.parametrize("name", ['section_for', 'pick_whole_file',
                                      'drop_wiring_keys', 'fix_num_x_ctrls'])
    def test_instance_functions(self, name):
        assert isroutine(getattr(totest, name))

    def test_doctests(self):
        result = doctest.testmod(totest)
        assert result.failed == 0
        assert result.attempted > 0


class TestSectionFor:
    @pytest.mark.parametrize("target, section", [
        ("binary_controls", "binary_controls"),
        ("binary_job", "binary_job"),
        ("star1_controls", "controls"),
        ("star2_job", "star_job"),
        ("star1_formation_step0_controls", "controls"),
        ("star2_formation_step3_job", "star_job"),
    ])
    def test_known_targets(self, target, section):
        assert totest.section_for(target) == section

    @pytest.mark.parametrize("target", ["binary", "star1", "controls",
                                        "binary_pgstar"])
    def test_unknown_target(self, target):
        with pytest.raises(ValueError, match="unknown target"):
            totest.section_for(target)


class TestLayerStack:
    def test_init(self):
        stack = totest.LayerStack("star2_job")
        assert stack.target == "star2_job"
        assert stack.section == "star_job"
        assert stack.layer_order == totest.LAYER_ORDER
        assert stack.merged() == {}

    def test_init_unknown_target(self):
        with pytest.raises(ValueError, match="unknown target"):
            totest.LayerStack("star1")

    def test_later_layer_wins_first_position_kept(self):
        stack = totest.LayerStack("star1_controls")
        stack.add("mesa_inlists", {"a": "1", "b": "1", "c": "1"})
        stack.add("posydon_inlists", {"b": "2", "d": "2"})
        stack.add("user_inlists", {"a": "3"})
        merged = stack.merged()
        assert merged == {"a": "3", "b": "2", "c": "1", "d": "2"}
        assert list(merged) == ["a", "b", "c", "d"]

    def test_layer_order_beats_add_order(self):
        stack = totest.LayerStack("binary_controls")
        stack.add("user_inlists", {"fm": "0.002d0"})
        stack.add("mesa_inlists", {"fm": "1d0", "fr": "0.01d0"})
        assert stack.merged() == {"fm": "0.002d0", "fr": "0.01d0"}

    def test_add_order_within_a_layer(self):
        stack = totest.LayerStack("binary_controls")
        stack.add("user_inlists", {"fm": "1"}, source="first")
        stack.add("user_inlists", {"fm": "2"}, source="second")
        assert stack.merged() == {"fm": "2"}

    def test_add_copies_values(self):
        values = {"a": "1"}
        stack = totest.LayerStack("binary_job")
        stack.add("generated", values)
        values["a"] = "2"
        assert stack.merged() == {"a": "1"}

    def test_unknown_layer(self):
        stack = totest.LayerStack("binary_job")
        with pytest.raises(ValueError, match="unknown layer 'special'"):
            stack.add("special", {})

    def test_custom_layer_order(self):
        stack = totest.LayerStack("binary_job", layer_order=("b", "a"))
        stack.add("a", {"x": "a"})
        stack.add("b", {"x": "b"})
        assert stack.merged() == {"x": "a"}

    def test_explain(self):
        stack = totest.LayerStack("star1_controls")
        stack.add("mesa_inlists", {"a": "1"}, source="defaults")
        stack.add("posydon_inlists", {"b": "2"}, source="common")
        stack.add("user_inlists", {"a": "3"}, source="mine")
        assert stack.explain("a") == [
            {"layer": "mesa_inlists", "source": "defaults", "value": "1"},
            {"layer": "user_inlists", "source": "mine", "value": "3"}]
        assert stack.explain("missing") == []

    def test_add_file_reads_the_target_section(self, posydon_inlist1):
        stack = totest.LayerStack("star1_job")
        added = stack.add_file("posydon_inlists", posydon_inlist1)
        assert added == {"pgstar_flag": ".false."}
        assert stack.explain("pgstar_flag")[0]["source"] == \
            str(posydon_inlist1)

    def test_add_file_skips_legacy_values(self, posydon_inlist1):
        stack = totest.LayerStack("star1_controls")
        stack.add_file("posydon_inlists", posydon_inlist1)
        assert "zams_filename" not in stack.merged()

    def test_add_file_with_transform(self, posydon_inlist1):
        stack = totest.LayerStack("star1_controls")
        stack.add_file("posydon_inlists", posydon_inlist1,
                       transform=lambda values: fix(drop(values)))
        assert stack.merged() == {"mesh_delta_coeff": "0.8d0",
                                  "x_ctrl(1:1)": "0d0"}

    def test_add_file_without_the_section(self, posydon_inlist1):
        stack = totest.LayerStack("binary_controls")
        assert stack.add_file("posydon_inlists", posydon_inlist1) == {}
        assert stack.merged() == {}

    def test_add_file_defaults_file(self, tmp_path):
        path = tmp_path / "binary_controls.defaults"
        path.write_text("fm = 1d0\nfr = 0.01d0\n")
        stack = totest.LayerStack("binary_controls")
        stack.add_file("mesa_inlists", path)
        assert stack.merged() == {"fm": "1d0", "fr": "0.01d0"}


def drop(values):
    return totest.drop_wiring_keys(values)


def fix(values):
    return totest.fix_num_x_ctrls(values)


class TestPickWholeFile:
    def test_highest_layer_wins(self):
        assert totest.pick_whole_file({"mesa_inlists": "m",
                                       "user_inlists": "u",
                                       "posydon_inlists": "p"}) == \
            ("user_inlists", "u")

    def test_none_is_skipped(self):
        assert totest.pick_whole_file({"mesa_inlists": "m",
                                       "posydon_inlists": None}) == \
            ("mesa_inlists", "m")

    def test_no_file(self):
        assert totest.pick_whole_file({}) is None
        assert totest.pick_whole_file({"user_inlists": None}) is None

    def test_unknown_layer(self):
        with pytest.raises(ValueError, match="unknown layer"):
            totest.pick_whole_file({"special": "x"})

    def test_custom_layer_order(self):
        assert totest.pick_whole_file({"a": "x", "b": "y"},
                                      layer_order=("b", "a")) == ("a", "x")


class TestLegacyKeyRules:
    def test_drop_wiring_keys(self):
        values = {"read_extra_star_job_inlist1": ".true.",
                  "inlist_names(1)": "'inlist1'",
                  "mesh_delta_coeff": "1d0"}
        assert totest.drop_wiring_keys(values) == {"mesh_delta_coeff": "1d0"}

    def test_fix_num_x_ctrls(self):
        values = {"x_logical_ctrl(1:num_x_ctrls)": ".false.",
                  "x_ctrl(2)": "1"}
        assert totest.fix_num_x_ctrls(values) == {
            "x_logical_ctrl(1:1)": ".false.", "x_ctrl(2)": "1"}


def legacy_merge(files, section, star_controls):
    """The merge loop of construct_static_inlist in bin/posydon-setup-grid."""
    final = {}
    for path in files:
        values = clean_inlist_file(str(path), section=section)[section]
        for key, value in values.items():
            if ('read_extra' in key) or ('inlist' in key):
                continue
            if star_controls and 'num_x_ctrls' in key:
                final[key.replace('num_x_ctrls', '1')] = value
            else:
                final[key] = value
    return final


@pytest.mark.skipif(not MESA_DEFAULTS.is_dir(),
                    reason="grid_params/defaults not present")
class TestSameAsTodaysSetup:
    """Stacking reproduces today's merge for MESA defaults + an inlist."""

    @pytest.mark.parametrize("target, defaults, section", [
        ("star1_controls", "star/controls.defaults", "&controls"),
        ("star1_job", "star/star_job.defaults", "&star_job"),
        ("binary_controls", "binary/binary_controls.defaults",
         "&binary_controls"),
        ("binary_job", "binary/binary_job.defaults", "&binary_job"),
    ])
    def test_merge(self, target, defaults, section, posydon_inlist1,
                   posydon_inlist_project):
        inlist = posydon_inlist_project if target.startswith("binary") \
            else posydon_inlist1
        star_controls = target.endswith("_controls") and \
            target.startswith("star")

        def clean(values):
            values = totest.drop_wiring_keys(values)
            return totest.fix_num_x_ctrls(values) if star_controls \
                else values

        stack = totest.LayerStack(target)
        stack.add_file("mesa_inlists", MESA_DEFAULTS / defaults,
                       transform=clean)
        stack.add_file("posydon_inlists", inlist, transform=clean)
        legacy = legacy_merge([MESA_DEFAULTS / defaults, inlist],
                              section, star_controls)
        assert list(stack.merged().items()) == list(legacy.items())
