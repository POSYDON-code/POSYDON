
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
logger = logging.getLogger(__name__)

class ColoredFormatter(logging.Formatter):
    """Custom formatter that adds colors to log level names."""

    LEVEL_COLORS = {
        'DEBUG': GRAY,
        'INFO': CYAN,
        'WARNING': YELLOW,
        'ERROR': RED,
        'CRITICAL': RED
    }

    def format(self, record):
        # Get the color for this log level
        color = self.LEVEL_COLORS.get(record.levelname, RESET)

        # Create the formatted message with colored level name
        formatted = f'[{color}{record.levelname}{RESET}] {record.getMessage()}'
        return formatted

def setup_logger(verbose=False):
    """Setup logging configuration based on verbosity level.

    Parameters
    ----------
    verbose : bool, optional
        If False, set INFO level (default)
        If True, set DEBUG level for detailed output
    """
    level = logging.DEBUG if verbose else logging.INFO
    handler = logging.StreamHandler()
    handler.setFormatter(ColoredFormatter())

    logging.basicConfig(
        level=level,
        handlers=[handler],
        force=True
    )
