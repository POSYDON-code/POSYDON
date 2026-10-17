"""Unit tests of posydon/grids/setup/record.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

# import other needed code for the tests, which is not already imported in the
# module you like to test
import configparser
from inspect import isroutine
from pathlib import Path

import pytest

# import the module which will be tested
import posydon.grids.setup.record as totest
from posydon.grids.setup.config import GridConfig
from posydon.grids.setup.manifest import Release
from posydon.grids.setup.resolve import ResolvedScenario
from posydon.grids.setup.scenario import ScenarioFiles


def parse(text):
    parser = configparser.ConfigParser(interpolation=None)
    parser.optionxform = str
    parser.read_string(text)
    return {name: dict(parser[name]) for name in parser.sections()}


class TestElements:
    def test_dir(self):
        assert {'RECORD_FILE', 'record_text', 'write_record'} <= \
            set(dir(totest))

    @pytest.mark.parametrize("name", ['record_text', 'write_record'])
    def test_instance_functions(self, name):
        assert isroutine(getattr(totest, name))


class TestRecord:
    def test_no_scenario(self, tmp_path):
        config = GridConfig()
        config.grid.type = "HMS-HMS"
        text = totest.record_text(config, None, tmp_path / "grid.ini",
                                  "/mesa", "2.2")
        assert text.startswith("; What this run directory")
        assert parse(text) == {"setup": {
            "inifile": str((tmp_path / "grid.ini").resolve()),
            "ini_format": "new", "posydon_version": "2.2",
            "mesa_dir": "/mesa", "grid_type": "HMS-HMS", "mode": "fixed"}}

    def test_release(self, tmp_path):
        config = GridConfig()
        config.grid.type = "CO-HeMS"
        config.scenario.release = "v2.0"
        config.scenario.rerun = "hot"
        release = Release(name="v2.0", folder=Path("/r"),
                          mesa_version="r11701", mesa="default", frozen=True)
        resolved = ResolvedScenario(
            files=ScenarioFiles(), repo=Path("/repo"), commit="abc123",
            release=release, dirty=False, remote="https://x/inlists",
            fingerprint="f" * 64, n_files=12)
        path = totest.write_record(tmp_path, config, resolved, "grid.ini",
                                   "/mesa", "2.2")
        assert path == tmp_path / "setup_record.ini"
        assert parse(path.read_text())["scenario"] == {
            "release": "v2.0", "rerun": "hot", "mesa_version": "r11701",
            "mesa": "default", "inlist_repo": "/repo",
            "inlist_remote": "https://x/inlists", "inlist_commit": "abc123",
            "inlist_dirty": "False", "inlist_fingerprint": "f" * 64,
            "inlist_files": "12"}

    def test_legacy(self):
        config = GridConfig(source_format="legacy")
        config.grid.type = "HMS-HMS"
        config.scenario.branch = "main"
        config.scenario.commit = "18710e9"
        config.scenario.system_type = "MS-MS"
        resolved = ResolvedScenario(files=ScenarioFiles(),
                                    repo=Path("/repo"), commit="18710e9abc")
        record = parse(totest.record_text(config, resolved, "g.ini", "/m",
                                          "2.2"))
        assert record["setup"]["ini_format"] == "legacy"
        assert record["scenario"] == {
            "source": "posydon", "branch": "main", "commit": "18710e9",
            "system_type": "MS-MS", "inlist_repo": "/repo",
            "inlist_remote": "None", "inlist_commit": "18710e9abc",
            "inlist_dirty": "None", "inlist_fingerprint": "None",
            "inlist_files": "0"}
