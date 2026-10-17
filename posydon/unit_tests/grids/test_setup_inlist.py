"""Unit tests of posydon/grids/setup/inlist.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

# import other needed code for the tests, which is not already imported in the
# module you like to test
import doctest
from inspect import isroutine
from pathlib import Path

import pytest

# import the module which will be tested
import posydon.grids.setup.inlist as totest
from posydon.utils.gridutils import clean_inlist_file

# MESA r11701's defaults, as vendored for the grid-setup rewrite tests from
# POSYDON-MESA-INLISTS (r11701/mesa/default/defaults); the comparison with
# today's reader below is skipped once that temporary folder is gone
MESA_DEFAULTS = Path(__file__).resolve().parents[3] / "grid_rewrite_tests" \
    / "setup_grid_data" / "mesa_defaults"

INLIST = """
! a comment line
&controls
    initial_mass = 10d0      ! solar masses
    x_ctrl(1) = 0.5d0
    log_directory = 'LOGS!1'  ! a ! inside quotes
    empty_name = ''
/ ! end of controls section

&star_job
    pgstar_flag = .false.
    pgstar_flag = .true.
/
"""


class TestElements:
    def test_dir(self):
        elements = {'LEGACY_SKIPPED_VALUES', 'Path', '__authors__',
                    '__builtins__', '__cached__', '__doc__', '__file__',
                    '__loader__', '__name__', '__package__', '__spec__',
                    'format_value', 'parse_inlist', 'read_inlist',
                    'render_inlist_section', 'render_inlist', 'strip_comment'}
        assert set(dir(totest)) == elements

    def test_instance_LEGACY_SKIPPED_VALUES(self):
        assert totest.LEGACY_SKIPPED_VALUES == ("''", "'.'")

    @pytest.mark.parametrize("name", ['strip_comment', 'parse_inlist',
                                      'read_inlist', 'format_value',
                                      'render_inlist_section',
                                      'render_inlist'])
    def test_instance_functions(self, name):
        assert isroutine(getattr(totest, name))

    def test_doctests(self):
        result = doctest.testmod(totest)
        assert result.failed == 0
        assert result.attempted > 0


class TestStripComment:
    @pytest.mark.parametrize("line, expected", [
        ("initial_mass = 10d0  ! in Msun", "initial_mass = 10d0"),
        ("! only a comment", ""),
        ("   ", ""),
        ("log_directory = 'LOGS!1'", "log_directory = 'LOGS!1'"),
        ('fmt = "a!b" ! c', 'fmt = "a!b"'),
        ("mixed = '\"' ! comment", "mixed = '\"'"),
    ])
    def test_strip_comment(self, line, expected):
        assert totest.strip_comment(line) == expected


class TestParseInlist:
    def test_sections_keys_and_values(self):
        parsed = totest.parse_inlist(INLIST)
        assert list(parsed) == ["controls", "star_job"]
        assert parsed["controls"] == {"initial_mass": "10d0",
                                      "x_ctrl(1)": "0.5d0",
                                      "log_directory": "'LOGS!1'",
                                      "empty_name": "''"}

    def test_last_value_wins_first_position_kept(self):
        parsed = totest.parse_inlist(
            "&s\n a = 1\n b = 2\n a = 3\n/\n")
        assert parsed["s"] == {"a": "3", "b": "2"}
        assert list(parsed["s"]) == ["a", "b"]

    def test_skip_values(self):
        parsed = totest.parse_inlist(
            INLIST, skip_values=totest.LEGACY_SKIPPED_VALUES)
        assert "empty_name" not in parsed["controls"]

    def test_default_section(self):
        parsed = totest.parse_inlist("a = 1\nb = 2\n",
                                        default_section="controls")
        assert parsed == {"controls": {"a": "1", "b": "2"}}

    def test_keys_after_closing_slash_go_to_default_section(self):
        parsed = totest.parse_inlist("&s\n a = 1\n/\nb = 2\n")
        assert parsed == {"s": {"a": "1"}, None: {"b": "2"}}

    def test_section_line_with_trailing_text(self):
        parsed = totest.parse_inlist("&binary_job  \n a = 1\n/\n")
        assert parsed == {"binary_job": {"a": "1"}}

    def test_lines_without_equal_sign_are_ignored(self):
        parsed = totest.parse_inlist("&s\n just text\n a = 1\n/\n")
        assert parsed == {"s": {"a": "1"}}

    def test_value_with_equal_sign(self):
        parsed = totest.parse_inlist("&s\n fmt = 'a=b'\n/\n")
        assert parsed == {"s": {"fmt": "'a=b'"}}


class TestReadInlist:
    def test_read_file(self, tmp_path):
        path = tmp_path / "inlist"
        path.write_text(INLIST)
        assert totest.read_inlist(path) == totest.parse_inlist(INLIST)

    def test_options_are_passed_on(self, tmp_path):
        path = tmp_path / "controls.defaults"
        path.write_text("a = 1\nb = ''\n")
        assert totest.read_inlist(
            path, default_section="controls",
            skip_values=totest.LEGACY_SKIPPED_VALUES) == {"controls":
                                                          {"a": "1"}}

    @pytest.mark.skipif(not MESA_DEFAULTS.is_dir(),
                        reason="vendored MESA defaults not present")
    @pytest.mark.parametrize("name, section", [
        ("star/controls.defaults", "&controls"),
        ("star/star_job.defaults", "&star_job"),
        ("binary/binary_controls.defaults", "&binary_controls"),
        ("binary/binary_job.defaults", "&binary_job"),
    ])
    def test_same_as_clean_inlist_file_on_mesa_defaults(self, name, section):
        """The legacy reader and this one agree on real MESA defaults."""
        path = MESA_DEFAULTS / name
        legacy = clean_inlist_file(str(path), section=section)
        new = totest.read_inlist(
            path, default_section=section,
            skip_values=totest.LEGACY_SKIPPED_VALUES)
        assert list(new[section].items()) == list(legacy[section].items())


class TestFormatValue:
    @pytest.mark.parametrize("value, expected", [
        (True, ".true."),
        (False, ".false."),
        (3, "3"),
        (3.0, "3"),
        (12.345, "1.2345000000d+01"),
        (1e-5, "1.0000000000d-05"),
        ("'LOGS1'", "'LOGS1'"),
        ("1d-3", "1d-3"),
    ])
    def test_format_value(self, value, expected):
        assert totest.format_value(value) == expected


class TestRender:
    def test_render_inlist_section(self):
        text = totest.render_inlist_section("controls", {"a": 1.5, "b": "'x'"},
                                      indent="\t")
        assert text == "&controls\n\ta = 1.5000000000d+00\n\tb = 'x'\n/\n"

    def test_round_trip(self):
        parsed = totest.parse_inlist(INLIST)
        text = totest.render_inlist(parsed)
        assert totest.parse_inlist(text) == parsed

    def test_empty_section(self):
        assert totest.render_inlist_section("binary_job", {}) == \
            "&binary_job\n/\n"
