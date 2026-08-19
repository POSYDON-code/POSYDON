"""Logging utilities for the posydon.CLI.grids package."""

import logging

# ANSI color codes
GREEN = '\033[92m'
GRAY = '\033[90m'
CYAN = '\033[96m'
YELLOW = '\033[93m'
MAGENTA = '\033[95m'
RED = '\033[91m'
RESET = '\033[0m'
BOLD = '\033[1m'

# Setup logger
logger = logging.getLogger('posydon.CLI.grids')


class ColoredFormatter(logging.Formatter):
    """Custom formatter that adds colors to log level names."""

    LEVEL_COLORS = {
        'DEBUG': GRAY,
        'INFO': CYAN,
        'WARNING': YELLOW,
        'ERROR': RED,
        'CRITICAL': RED
    }

    def format(self, record: logging.LogRecord) -> str:
        """Format the record with a colored level name.

        Args:
            record: The log record to format.

        Returns:
            The formatted string with the colored level name.
        """
        color = self.LEVEL_COLORS.get(record.levelname, RESET)
        formatted = f'[{color}{record.levelname}{RESET}] {record.getMessage()}'
        return formatted


def setup_logger(verbose: bool = False) -> None:
    """Setup logging configuration based on verbosity level.

    Args:
        verbose: If True set DEBUG level, otherwise INFO level (default).
    """
    level = logging.DEBUG if verbose else logging.INFO
    handler = logging.StreamHandler()
    handler.setFormatter(ColoredFormatter())

    logging.basicConfig(
        level=level,
        handlers=[handler],
        force=True
    )
