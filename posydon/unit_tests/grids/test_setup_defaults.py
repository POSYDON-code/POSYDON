"""Unit tests of posydon/grids/setup/defaults.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

# import other needed code for the tests, which is not already imported in the
# module you like to test
from inspect import isroutine
from pathlib import Path

import pytest

# import the module which will be tested
import posydon.grids.setup.defaults as totest
from posydon.grids.setup.config import ConfigError, GridConfig
from posydon.grids.setup.grid_types import GRID_TYPES
from posydon.grids.setup.manifest import MESA_DEFAULTS, Release
from posydon.grids.setup.resolve import ResolvedScenario
from posydon.grids.setup.scenario import ScenarioFiles

DEFAULTS = {name: Path("/d") / name for name in MESA_DEFAULTS}


@pytest.fixture
def release(tmp_path):
    for rel in MESA_DEFAULTS.values():
        path = tmp_path / "r11701" / "mesa" / "default" / "defaults" / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("")
    return Release(name="v2.0",
                   folder=tmp_path / "r11701" / "releases" / "v2.0",
                   mesa_version="r11701", mesa="default", frozen=True)


class TestElements:
    def test_dir(self):
        elements = {'ConfigError', 'MESA_EXTRAS', 'Path', '__authors__',
                    '__builtins__', '__cached__', '__doc__', '__file__',
                    '__loader__', '__name__', '__package__', '__spec__',
                    'fill_defaults', 'get_grid_type', 'mesa_defaults_targets',
                    'mesa_installation'}
        assert set(dir(totest)) == elements

    @pytest.mark.parametrize("name", ['fill_defaults',
                                      'mesa_defaults_targets'])
    def test_instance_functions(self, name):
        assert isroutine(getattr(totest, name))


class TestMesaDefaultsTargets:
    def test_binary_from_zams(self):
        targets = totest.mesa_defaults_targets(GRID_TYPES["HMS-HMS"],
                                               DEFAULTS)
        assert list(targets) == ["binary_controls", "binary_job",
                                 "star1_job", "star1_controls",
                                 "star2_job", "star2_controls"]
        assert targets["star2_controls"] == DEFAULTS["star_controls"]

    def test_formed_star(self):
        targets = totest.mesa_defaults_targets(GRID_TYPES["CO-HeMS"],
                                               DEFAULTS)
        assert list(targets)[:2] == ["star1_formation_every_step_job",
                                     "star1_formation_every_step_controls"]
        assert "star1_controls" in targets

    @pytest.mark.parametrize("name", ["single_HMS", "single_HeMS"])
    def test_single_star(self, name):
        targets = totest.mesa_defaults_targets(GRID_TYPES[name], DEFAULTS)
        assert list(targets) == [
            "binary_controls", "binary_job",
            "star1_formation_every_step_job",
            "star1_formation_every_step_controls",
            "star2_job", "star2_controls"]


class TestFillDefaults:
    def config(self, grid_type="CO-HeMS"):
        config = GridConfig()
        config.grid.type = grid_type
        return config

    def test_legacy_is_left_alone(self):
        config = self.config()
        config.source_format = "legacy"
        totest.fill_defaults(config, None, "/mesa")
        assert config.extras == {}
        assert config.inlists["mesa_inlists"].targets == {}

    def test_fills_both(self, release):
        config = self.config()
        resolved = ResolvedScenario(files=ScenarioFiles(), repo=Path("/r"),
                                    release=release)
        totest.fill_defaults(config, resolved, "/mesa")
        folder = release.folder.parents[1] / "mesa" / "default" / "defaults"
        assert config.inlists["mesa_inlists"].targets["binary_job"] == \
            [folder / "binary" / "binary_job.defaults"]
        assert config.extras["star_run"] == Path("/mesa/star/work/src/run.f")
        assert list(config.extras) == list(totest.MESA_EXTRAS)
        assert config.default_extras

    def test_keeps_what_the_ini_sets(self):
        config = self.config()
        config.inlists["mesa_inlists"].add("binary_job", Path("/mine"))
        config.extras = {"star_run": Path("/my/run.f")}
        totest.fill_defaults(config, None, "/mesa")
        assert config.inlists["mesa_inlists"].targets == {
            "binary_job": [Path("/mine")]}
        assert config.extras == {"star_run": Path("/my/run.f")}
        assert not config.default_extras

    @pytest.mark.parametrize("resolved", [
        None, ResolvedScenario(files=ScenarioFiles(), repo=Path("/r"))])
    def test_needs_a_release(self, resolved):
        with pytest.raises(ConfigError, match=r"\[mesa_inlists\]: needed "
                           "without a \\[scenario\\] release"):
            totest.fill_defaults(self.config(), resolved, "/mesa")
