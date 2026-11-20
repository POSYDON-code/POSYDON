"""Unit tests of posydon/grids/setup/grid_types.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

# import other needed code for the tests, which is not already imported in the
# module you like to test
import doctest
from dataclasses import FrozenInstanceError, is_dataclass
from inspect import isroutine

import pytest

# import the module which will be tested
import posydon.grids.setup.grid_types as totest
from posydon.grids.setup import config, legacy


class TestElements:
    def test_dir(self):
        elements = {'ConfigError', 'GRID_TYPES', 'GridType', 'STAR_STARTS',
                    '__authors__', '__builtins__', '__cached__', '__doc__',
                    '__file__', '__loader__', '__name__', '__package__',
                    '__spec__', 'dataclass', 'describe_grid_types',
                    'get_grid_type'}
        assert set(dir(totest)) == elements

    def test_instance_GridType(self):
        assert is_dataclass(totest.GridType)

    @pytest.mark.parametrize("name", ['get_grid_type',
                                      'describe_grid_types'])
    def test_instance_functions(self, name):
        assert isroutine(getattr(totest, name))

    def test_doctests(self):
        result = doctest.testmod(totest)
        assert result.failed == 0
        assert result.attempted > 0


class TestRegistry:
    @pytest.mark.parametrize("name", list(totest.GRID_TYPES))
    def test_config_accepts_every_name(self, name, tmp_path):
        """config.py validates [grid] type against this registry."""
        grid = config._read_grid({"type": name, "points": "g.csv"}, tmp_path)
        assert grid.type == name

    def test_legacy_mapping_points_into_the_registry(self):
        assert set(legacy.LEGACY_GRID_TYPES.values()) <= \
            set(totest.GRID_TYPES)

    def test_legacy_folders_match(self):
        """The legacy system type is the folder, except HMS-HMS/MS-MS."""
        for (system_type, single), name in \
                legacy.LEGACY_GRID_TYPES.items():
            grid_type = totest.GRID_TYPES[name]
            assert grid_type.single_star == single
            expected = None if system_type in ("HMS-HMS", "MS-MS") \
                else system_type
            assert grid_type.folder == expected

    @pytest.mark.parametrize("name, folder, star1, star2, single", [
        ("HMS-HMS", None, "zams", "zams", False),
        ("CO-HMS_RLO", "CO-H_star", "zams", "point_mass", False),
        ("CO-HeMS", "CO-He_star", "formation", "point_mass", False),
        ("HeMS-HMS", "HeMS-HMS", "formation", "zams", False),
        ("single_HMS", None, "zams", None, True),
        ("single_HeMS", "CO-He_star", "formation", None, True),
    ])
    def test_table(self, name, folder, star1, star2, single):
        grid_type = totest.GRID_TYPES[name]
        assert (grid_type.folder, grid_type.star1, grid_type.star2,
                grid_type.single_star) == (folder, star1, star2, single)
        assert grid_type.single_star_flag == single
        assert grid_type.name == name
        assert grid_type.description


class TestGridType:
    def test_frozen(self):
        with pytest.raises(FrozenInstanceError):
            totest.GRID_TYPES["HMS-HMS"].folder = "x"

    def test_unknown_start(self):
        with pytest.raises(ValueError, match="unknown start 'merger'"):
            totest.GridType("X", "x", star2="merger")

    @pytest.mark.parametrize("star2, single", [(None, False),
                                               ("zams", True)])
    def test_single_star_needs_no_star2(self, star2, single):
        with pytest.raises(ValueError, match="single-star grid"):
            totest.GridType("X", "x", star2=star2, single_star=single)


class TestLookup:
    def test_get_grid_type(self):
        assert totest.get_grid_type("single_HMS").single_star is True

    def test_unknown(self):
        with pytest.raises(config.ConfigError,
                           match=r"\[grid\] type: unknown grid type 'MS-MS'"):
            totest.get_grid_type("MS-MS")

    def test_describe(self):
        lines = totest.describe_grid_types().splitlines()
        assert len(lines) == 6
        assert lines[2].startswith("CO-HeMS       a compact object")
