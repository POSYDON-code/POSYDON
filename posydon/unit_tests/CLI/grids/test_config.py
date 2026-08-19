"""Unit tests for posydon.CLI.grids.config."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import os
import re

import pandas
import pytest

from posydon.CLI.grids import config

totest = config


@pytest.fixture(autouse=True)
def _env(tmp_path, monkeypatch):
    """Always provide HOME/MESA_DIR and restore cwd after every test."""
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("MESA_DIR", str(tmp_path / "mesa"))
    original = os.getcwd()
    yield
    os.chdir(original)


class TestParseConfigFile:
    """Tests for parse_config_file()."""

    def test_returns_parse_inifile_tuple(self, monkeypatch):
        """parse_config_file returns the 4-tuple from configfile.parse_inifile."""
        expected = ({"a": 1}, {"b": 2}, {"c": 3}, {"d": 4})
        monkeypatch.setattr(totest.configfile, "parse_inifile", lambda path: expected)
        assert totest.parse_config_file("some.ini") == expected


class TestNormalizeDefaults:
    """Tests for normalize_defaults()."""

    def test_adds_defaults_when_absent(self):
        """Missing keys get their default values added in place."""
        run_parameters = {}
        mesa_inlists = {}
        totest.normalize_defaults(run_parameters, mesa_inlists)
        assert mesa_inlists["scenario"] is None
        assert run_parameters["keep_profiles"] is False
        assert run_parameters["keep_photos"] is False

    def test_leaves_existing_values_alone(self):
        """Existing keys are not overwritten by normalize_defaults."""
        run_parameters = {"keep_profiles": True, "keep_photos": True}
        mesa_inlists = {"scenario": "HeMS-HMS"}
        totest.normalize_defaults(run_parameters, mesa_inlists)
        assert mesa_inlists == {"scenario": "HeMS-HMS"}
        assert run_parameters == {"keep_profiles": True, "keep_photos": True}


class TestValidateConfig:
    """Tests for validate_config()."""

    def test_missing_grid_raises(self, tmp_path):
        """A non-existent grid path raises the exact ValueError."""
        with pytest.raises(
            ValueError,
            match="Supplied grid does not exist, please check your path and try again",
        ):
            totest.validate_config({"grid": str(tmp_path / "missing.csv")}, "fixed")

    def test_directory_grid_is_valid(self, tmp_path):
        """A directory grid path passes validation."""
        grid_dir = tmp_path / "grid_dir"
        grid_dir.mkdir()
        totest.validate_config({"grid": str(grid_dir)}, "fixed")

    def test_dynamic_without_psycris_raises(self, tmp_path):
        """A dynamic grid without a psycris inifile raises the exact ValueError."""
        grid = tmp_path / "grid.csv"
        grid.write_text("m1\n10.0\n")
        with pytest.raises(
            ValueError,
            match=re.escape(
                "Please add psycris inifile to the [run_parameters] section of the inifile."
            ),
        ):
            totest.validate_config({"grid": str(grid)}, "dynamic")

    def test_fixed_without_psycris_is_valid(self, tmp_path):
        """A fixed grid does not require a psycris inifile."""
        grid = tmp_path / "grid.csv"
        grid.write_text("m1\n10.0\n")
        totest.validate_config({"grid": str(grid)}, "fixed")

    def test_dynamic_with_psycris_is_valid(self, tmp_path):
        """A dynamic grid with a psycris inifile passes validation."""
        grid = tmp_path / "grid.csv"
        grid.write_text("m1\n10.0\n")
        totest.validate_config(
            {"grid": str(grid), "psycris_inifile": "x.ini"}, "dynamic"
        )


class TestReadGridFile:
    """Tests for read_grid_file()."""

    def test_csv(self, tmp_path):
        """A csv grid is read with pandas and its path is returned as-is."""
        grid = tmp_path / "grid.csv"
        grid.write_text("m1,m2\n10.0,8.0\n11.0,9.0\n")
        df, name = totest.read_grid_file(str(grid))
        assert list(df.columns) == ["m1", "m2"]
        assert len(df) == 2
        assert name == str(grid)

    def test_h5(self, tmp_path, monkeypatch):
        """An h5 grid is loaded via PSyGrid and its path is returned as-is."""
        grid_df = pandas.DataFrame({"m1": [1.0, 2.0]})

        class FakeGrid:
            def load(self, path):
                pass

            def get_pandas_initial_final(self):
                return grid_df

            def close(self):
                pass

        monkeypatch.setattr(totest, "PSyGrid", FakeGrid)
        path = str(tmp_path / "grid.h5")
        df, name = totest.read_grid_file(path)
        assert df.equals(grid_df)
        assert name == path

    def test_directory(self, tmp_path, monkeypatch):
        """A directory grid is slim-created into ./fixed_grid_results.h5."""
        grid_dir = tmp_path / "grid_dir"
        grid_dir.mkdir()
        grid_df = pandas.DataFrame({"m1": [1.0, 2.0]})

        class FakeGrid:
            def create(self, *args, **kwargs):
                return self

            def load(self, path):
                pass

            def get_pandas_initial_final(self):
                return grid_df

            def close(self):
                pass

        monkeypatch.setattr(totest, "PSyGrid", FakeGrid)
        os.chdir(tmp_path)
        df, name = totest.read_grid_file(str(grid_dir))
        assert df.equals(grid_df)
        assert name == os.path.join(os.getcwd(), "fixed_grid_results.h5")

    def test_unknown_format_raises(self, tmp_path):
        """An unrecognized grid format raises the exact ValueError."""
        bogus = tmp_path / "grid.txt"
        bogus.write_text("data")
        with pytest.raises(
            ValueError,
            match="Grid format not recognized, please feed in an acceptable format: csv",
        ):
            totest.read_grid_file(str(bogus))


class TestResolveExtras:
    """Tests for resolve_extras()."""

    def test_lower_precedence_types_set_to_none(self, monkeypatch, capsys):
        """Only the highest-precedence binary_extras type survives."""
        warnings_calls = []
        monkeypatch.setattr(
            totest,
            "Pwarn",
            lambda message, category: warnings_calls.append((message, category)),
        )
        mesa_extras = {
            "mesa_binary_extras": "mesa.f",
            "posydon_binary_extras": "posydon.f",
            "user_binary_extras": "user.f",
            "mesa_star_run": "star.f",
        }
        result = totest.resolve_extras(mesa_extras)
        assert result["mesa_binary_extras"] is None
        assert result["posydon_binary_extras"] is None
        assert result["user_binary_extras"] == "user.f"
        assert result["mesa_star_run"] == "star.f"
        assert warnings_calls == [
            (
                "Section mesa_extras value mesa_binary_extras is being set to None",
                "ReplaceValueWarning",
            ),
            (
                "Section mesa_extras value posydon_binary_extras is being set to None",
                "ReplaceValueWarning",
            ),
        ]
        assert capsys.readouterr().out == "WE ARE USING THE EXTRA FILE FROM TYPE user\n"

    def test_only_highest_precedence_no_change(self, monkeypatch, capsys):
        """A single binary_extras type is left untouched."""
        warnings_calls = []
        monkeypatch.setattr(
            totest,
            "Pwarn",
            lambda message, category: warnings_calls.append((message, category)),
        )
        mesa_extras = {"user_binary_extras": "user.f"}
        result = totest.resolve_extras(mesa_extras)
        assert result == {"user_binary_extras": "user.f"}
        assert warnings_calls == []
        assert capsys.readouterr().out == "WE ARE USING THE EXTRA FILE FROM TYPE user\n"
