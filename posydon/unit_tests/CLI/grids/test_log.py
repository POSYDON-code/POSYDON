"""Unit tests for posydon.CLI.grids.log."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import logging
import os

import pytest

from posydon.CLI.grids import log

totest = log


@pytest.fixture(autouse=True)
def _env(tmp_path, monkeypatch):
    """Always provide HOME/MESA_DIR and restore cwd after every test."""
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("MESA_DIR", str(tmp_path / "mesa"))
    original = os.getcwd()
    yield
    os.chdir(original)


class TestColoredFormatter:
    """Tests for the ColoredFormatter class."""

    @staticmethod
    def _make_record(level, levelname=None, msg="message"):
        record = logging.LogRecord(
            name="posydon.CLI.grids",
            level=level,
            pathname=__file__,
            lineno=1,
            msg=msg,
            args=(),
            exc_info=None,
        )
        if levelname is not None:
            record.levelname = levelname
        return record

    def test_info_level_gets_cyan(self):
        """An INFO record is formatted with the cyan color code."""
        formatter = totest.ColoredFormatter()
        formatted = formatter.format(self._make_record(logging.INFO, "INFO"))
        assert formatted == "[\x1b[96mINFO\x1b[0m] message"

    def test_debug_level_gets_gray(self):
        """A DEBUG record is formatted with the gray color code."""
        formatter = totest.ColoredFormatter()
        formatted = formatter.format(self._make_record(logging.DEBUG, "DEBUG"))
        assert formatted == "[\x1b[90mDEBUG\x1b[0m] message"

    def test_unknown_levelname_falls_back_to_reset(self):
        """An unknown level name falls back to the RESET color."""
        formatter = totest.ColoredFormatter()
        formatted = formatter.format(self._make_record(logging.DEBUG, "FOO"))
        assert formatted == "[\x1b[0mFOO\x1b[0m] message"


class TestSetupLogger:
    """Tests for setup_logger()."""

    def test_verbose_true_sets_debug_level(self):
        """setup_logger(verbose=True) sets the root logger level to DEBUG."""
        totest.setup_logger(verbose=True)
        assert logging.getLogger().level == logging.DEBUG

    def test_verbose_false_sets_info_level(self):
        """setup_logger(verbose=False) sets the root logger level to INFO."""
        totest.setup_logger(verbose=False)
        assert logging.getLogger().level == logging.INFO

    def test_logger_singleton(self):
        """The module logger is logging.getLogger('posydon.CLI.grids')."""
        assert totest.logger is logging.getLogger("posydon.CLI.grids")

    def test_logger_level_toggles_with_setup_logger(self):
        """The singleton logger's effective level follows setup_logger()."""
        totest.setup_logger(verbose=True)
        assert totest.logger.getEffectiveLevel() == logging.DEBUG
        totest.setup_logger(verbose=False)
        assert totest.logger.getEffectiveLevel() == logging.INFO
