import argparse
import glob
import logging
import os
import shutil
import subprocess

import pandas as pd

from posydon.active_learning.psy_cris.utils import parse_inifile
from posydon.grids.psygrid import PSyGrid
from posydon.utils import configfile
from posydon.utils import gridutils as utils

# Define column types and their filenames
column_types = {'star_history_columns'  :'history_columns.list',
                'binary_history_columns':'binary_history_columns.list',
                'profile_columns'       :'profile_columns.list'}
column_filenames = ['history_columns.list', 'binary_history_columns.list', 'profile_columns.list']

# Define extras keys
extras_keys = ['makefile_binary', 'makefile_star', 'binary_run',
               'star_run', 'run_binary_extras', 'run_star_binary_extras', 'run_star1_extras', 'run_star2_extras',]

# define inlist keys
inlist_keys = ['binary_controls', 'binary_job',
               'star1_controls', 'star1_job',
               'star2_controls', 'star2_job']

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

def setup_logger(verbose=None):
    """Setup logging configuration based on verbosity level.

    Parameters
    ----------
    verbose : str or None, optional
        Verbose logging option:
        None = No verbose logging (INFO level to console only)
        "console" = DEBUG level output to console
        filename = DEBUG level output to specified file
        Default is None.
    """
    if verbose is None:
        level = logging.INFO
        # Only setup console handler for INFO level
        handler = logging.StreamHandler()
        handler.setFormatter(ColoredFormatter())

        # Configure the root logger
        logging.basicConfig(
            level=level,
            handlers=[handler],
            force=True
        )
    elif verbose == "console":
        level = logging.DEBUG
        print(f"{CYAN}DEBUG verbosity enabled: Detailed output will be shown on console.{RESET}")

        # Create a custom handler with colored formatter
        handler = logging.StreamHandler()
        handler.setFormatter(ColoredFormatter())

        # Configure the root logger
        logging.basicConfig(
            level=level,
            handlers=[handler],
            force=True
        )
    else:
        # verbose is a filename
        level = logging.DEBUG
        print(f"{CYAN}DEBUG verbosity enabled: Detailed output will be logged to {verbose}.{RESET}")

        # Create both file and console handlers
        file_handler = logging.FileHandler(verbose)
        console_handler = logging.StreamHandler()

        # Use colored formatter for console, plain formatter for file
        console_handler.setFormatter(ColoredFormatter())
        file_handler.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s',
                                   datefmt='%Y-%m-%d %H:%M'))

        # Configure the root logger
        logging.basicConfig(
            level=level,
            handlers=[file_handler, console_handler],
            force=True
        )


def check_file_exist(file_path, raise_error=True):
    """Check if a file exists at the given path

    Parameters
    ----------
    file_path : str
        Path to the file to check
    raise_error : bool, optional
        If True, raise ValueError when file doesn't exist.
        If False, return boolean. Default is True.

    Returns
    -------
    bool
        True if file exists, False otherwise (only when raise_error=False)
    """
    exists = os.path.exists(file_path)

    if not exists:
        if raise_error:
            print(f"File {file_path} does not exist")
            raise ValueError(f"File {file_path} does not exist")
        else:
            return False

    return True

def setup_inlist_repository(inlist_repository, MESA_version):
    """Setup the inlist repository by creating it if it does not exist

    Parameters
    ----------
    inlist_repository : str
        Path to the inlist repository where we will store inlists
    base : str
        Path to the base to use for the run
    """
    POSYDON_inlist_URL = 'https://github.com/POSYDON-code/POSYDON-MESA-INLISTS.git'
    logger.info("Loading repository for inlists.")

    if inlist_repository is None:
        # check if inlist_repository is provided
        logger.info("No inlist_repository provided in inifile.")
        logger.info('Using your home folder as the inlist repository')
        inlist_repository = os.path.expanduser('~')

    # check if the inlist repository path exists
    if not os.path.exists(inlist_repository):
        logger.debug(f"Creating inlist repository at {inlist_repository}")
        os.makedirs(inlist_repository)

    # check if it contains anything and if not, clone the repo
    if os.listdir(inlist_repository):
        logger.debug(f"Files found in inlist repository: {os.listdir(inlist_repository)}")
        logger.info("Assuming the POSYDON inlist repository is already cloned into this folder!")
    else:
        out = subprocess.run(['git', 'clone', POSYDON_inlist_URL, inlist_repository],
                               capture_output=True,
                               text=True,
                               check=True,)
        if out.stderr:
            logger.error(out.stderr)
        else:
            logger.debug("Cloned the POSYDON inlist repository successfully.")

    # update the repository
    logger.debug("Updating the inlist repository to the latest version.")
    # TODO: Re-enable git pull after testing
    logger.debug("Currently disabled for testing purposes")
    # out = subprocess.run(['git', 'pull'],
    #                cwd=inlist_repository,
    #                capture_output=True,
    #                text=True,
    #                check=True,)

    # check if the base is available as a folder in the repository
    version_root_path = os.path.join(inlist_repository, MESA_version)
    if not os.path.exists(version_root_path):
        logger.error(version_root_path)
        raise ValueError("The provided MESA version does not exist in the inlist repository, please check your provided MESA version and try again.")

    return version_root_path

def setup_MESA_defaults(path_to_version):
    """Setup the MESA default base inlists, extras and columns

    Parameters
    ----------
    path_to_version : str
        Path to the MESA version in the inlist repository
        (root directory of the version)

    Returns
    -------
    MESA_default_inlists : dict
        Dictionary of MESA default inlists paths
    MESA_default_extras : dict
        Dictionary of MESA default extras paths
    MESA_default_columns : dict
        Dictionary of MESA default column files paths
    """
    MESA_DIR = os.environ['MESA_DIR']
    logger.debug(f"Setting up MESA defaults from MESA_DIR: {MESA_DIR}")
    #----------------------------------
    #            Inlists
    #----------------------------------
    # Common inlists
    # TODO: These are currently stored in the POSYDON inlist repository:
    # The default MESA inlists with r11701 has a bug and we needed to fix it.
    # We can add this to our changed MESA version, when we release it?
    # Then this can be reverted to the MESA default inlists!
    MESA_default_inlists = {}

    MESA_defaults_inlists_path = os.path.join(path_to_version,
                                              'MESA_defaults',
                                              'inlists')
    MESA_default_inlists['binary_controls'] = [os.path.join(MESA_defaults_inlists_path,
                                                            'binary',
                                                            'binary_controls.defaults')]
    MESA_default_inlists['binary_job'] = [os.path.join(MESA_defaults_inlists_path,
                                                       'binary',
                                                      'binary_job.defaults')]
    MESA_default_inlists['star1_controls'] = [os.path.join(MESA_defaults_inlists_path,
                                                          'star',
                                                          'controls.defaults')]
    MESA_default_inlists['star1_job'] = [os.path.join(MESA_defaults_inlists_path,
                                                     'star',
                                                     'star_job.defaults')]
    MESA_default_inlists['star2_controls'] = [os.path.join(MESA_defaults_inlists_path,
                                                          'star',
                                                          'controls.defaults')]
    MESA_default_inlists['star2_job'] = [os.path.join(MESA_defaults_inlists_path,
                                                     'star',
                                                     'star_job.defaults')]

    #----------------------------------
    #              EXTRAS
    #----------------------------------
    MESA_default_extras = {}

    # Helper to build MESA work directory paths
    def mesa_path(module, *parts):
        return os.path.join(MESA_DIR, module, 'work', *parts)

    # Makefiles
    MESA_default_extras['makefile_binary'] = mesa_path('binary',
                                                       'make',
                                                       'makefile')
    MESA_default_extras['makefile_star'] = mesa_path('star',
                                                     'make',
                                                     'makefile')

    # Run files
    MESA_default_extras['star_run'] = mesa_path('star',
                                                'src',
                                                'run.f')
    MESA_default_extras['binary_run'] = mesa_path('binary',
                                                  'src',
                                                  'binary_run.f')

    # Extras files for binary evolution
    MESA_default_extras['run_binary_extras'] = mesa_path('binary',
                                                     'src',
                                                     'run_binary_extras.f')
    MESA_default_extras['run_star_binary_extras'] = mesa_path('binary',
                                                         'src',
                                                         'run_star_extras.f')

    # star1_extras and star2_extras are needed for pre-MS formation steps (if any).
    # During binary evolution, star_binary_extras is used for both stars.
    # #Both stars use the same single-star module extras file
    star_extras_path = mesa_path('star', 'src', 'run_star_extras.f')
    MESA_default_extras['run_star1_extras'] = star_extras_path
    MESA_default_extras['run_star2_extras'] = star_extras_path

    # Verify all extras files exist
    for _, path in MESA_default_extras.items():
        check_file_exist(path)

    #----------------------------------
    #             Columns
    #----------------------------------

    # Column files from MESA defaults
    MESA_default_columns = {
        'star_history_columns': os.path.join(MESA_DIR, 'star',
                                             'defaults', 'history_columns.list'),
        'binary_history_columns': os.path.join(MESA_DIR, 'binary',
                                               'defaults', 'binary_history_columns.list'),
        'profile_columns': os.path.join(MESA_DIR, 'star',
                                        'defaults', 'profile_columns.list')
    }

    # Verify all column files exist
    for _, path in MESA_default_columns.items():
        check_file_exist(path)

    return MESA_default_inlists, MESA_default_extras, MESA_default_columns

def setup_POSYDON(path_to_version, base, system_type):
    """Setup the POSYDON inlists, extras and columns

    Parameters
    ----------
    path_to_version : str
        Path to the POSYDON version in the inlist repository
        (root directory of the version)
    base : str
        Path to the POSYDON version to use for the run.
        If "MESA", returns "empty" dictionaries.
    system_type : str
        Type of binary system

    Returns
    -------
    POSYDON_inlists : dict
        Dictionary of POSYDON inlists paths
    POSYDON_extras : dict
        Dictionary of POSYDON extras paths
    POSYDON_columns : dict
        Dictionary of POSYDON column files paths
    """
    # If user wants to use MESA base only, return empty dictionaries
    if base == "MESA" or base[0] == "MESA":
        POSYDON_columns = {name: None for name in column_types}
        POSYDON_inlists = {}
        POSYDON_extras = {key: None for key in extras_keys}
        return POSYDON_inlists, POSYDON_extras, POSYDON_columns

    # Setup POSYDON base path
    if len(base) == 2:
        POSYDON_path = os.path.join(path_to_version, base[0], base[1])
    else:
        POSYDON_path = os.path.join(path_to_version, base[0], base[1], base[2])
    check_file_exist(POSYDON_path)

    logger.debug(f"Setting up POSYDON configuration: {POSYDON_path}")

    #----------------------------------
    #            Inlists
    #----------------------------------
    POSYDON_inlists = {}
    # Common inlists
    # TODOL these are not all needed for single stars.
    common_inlists_path = os.path.join(POSYDON_path, 'common_inlists')
    POSYDON_inlists['binary_controls'] = [os.path.join(common_inlists_path, 'inlist_project')]
    POSYDON_inlists['binary_job'] = [os.path.join(common_inlists_path, 'inlist_project')]
    # setup star1 inlists for binaries
    POSYDON_inlists['star1_controls'] = [os.path.join(common_inlists_path, 'inlist1')]
    POSYDON_inlists['star1_job'] = [os.path.join(common_inlists_path, 'inlist1')]
    # setup star2 inlists for binaries
    POSYDON_inlists['star2_controls'] = [os.path.join(common_inlists_path, 'inlist2')]
    POSYDON_inlists['star2_job'] = [os.path.join(common_inlists_path, 'inlist2')]


    if system_type == 'single_HMS':
        # setup the paths to extra single star inlists
        single_star_inlist_path = os.path.join(POSYDON_path,
                                               'single_HMS',
                                               'single_star_inlist')
        POSYDON_inlists['star1_controls'].append(single_star_inlist_path)
    elif system_type == 'single_HeMS':
        # setup the paths to extra single star He inlists
        # The HeMS single star inlists have two steps
        # step 1: create HeMS star
        # step 2: evolve HeMS star
        helium_star_inlist_step1 = os.path.join(POSYDON_path,
                                                 'single_HeMS',
                                                 'inlist_step1')
        helium_star_inlist_step2 = os.path.join(POSYDON_path,
                                                 'single_HeMS',
                                                 'inlist_step2')
        # We need to also include the HMS single star inlist to set up
        # the single star evolution
        single_star_inlist_path = os.path.join(POSYDON_path,
                                               'single_HMS',
                                               'single_star_inlist')

        single_helium_inlists = [helium_star_inlist_step1,
                                 helium_star_inlist_step2,
                                 single_star_inlist_path]

        # the helium star setup steps contain control & job in the same inlist files
        POSYDON_inlists['star1_controls'].extend(single_helium_inlists)
        # We don't need to add the single star inlist again for the hob,
        # since it only contains controls section in the file
        POSYDON_inlists['star1_job'].extend(single_helium_inlists)

    elif system_type in ['HMS-HeMS', 'HeMS-HeMS']:
        pass
    elif system_type in ['CO-HMS', 'CO-HeMS']:
        pass
    elif system_type in ['HMS-HMS']:
        # the common inlists are sufficient
        pass
    else:
        raise ValueError(f"System type {system_type} not recognized.")

    #----------------------------------
    #            Extras
    #----------------------------------

    POSYDON_extras = {}
    POSYDON_extras['run_binary_extras'] = os.path.join(POSYDON_path,
                                                   'extras_files',
                                                   'run_binary_extras.f')
    POSYDON_extras['run_star_binary_extras'] = os.path.join(POSYDON_path,
                                                       'extras_files',
                                                       'run_star_extras.f')
    POSYDON_extras['run_star1_extras'] = os.path.join(POSYDON_path,
                                                 'extras_files',
                                                 'run_star_extras.f')

    #----------------------------------
    #             Columns
    #----------------------------------

    # Setup POSYDON columns
    POSYDON_columns = {}
    for name, filename in column_types.items():
        file = os.path.join(POSYDON_path, 'column_files', filename)
        # only add if the file exists
        if check_file_exist(file, raise_error=False):
            POSYDON_columns[name] = file
        else:
            POSYDON_columns[name] = None

    return POSYDON_inlists, POSYDON_extras, POSYDON_columns

def setup_user(user_mesa_inlists, user_mesa_extras):
    """Separates out user inlists, extras and columns

    Parameters
    ----------
    user_mesa_inlists : dict
        Dictionary of user inlists paths from the inifile
    user_mesa_extras : dict
        Dictionary of user extras paths from the inifile


    Returns
    -------
    user_mesa_inlists : dict
        Dictionary of user inlists paths
    user_mesa_extras : dict
        Dictionary of user extras paths
    user_columns : dict
        Dictionary of user column files paths
    """

    #----------------------------------
    #            Inlists
    #----------------------------------
    # separate out inlists from user inlists
    user_inlists = {}
    for key in inlist_keys:
        if key not in user_mesa_inlists.keys():
            user_inlists[key] = []
        else:
            user_inlists[key] = [user_mesa_inlists[key]]
            check_file_exist(user_mesa_inlists[key])

    #----------------------------------
    #            Extras
    #----------------------------------
    # separate out extras from user inlists

    user_extras = {}
    for key in extras_keys:
        if key in user_mesa_extras.keys():
            user_extras[key] = user_mesa_extras[key]
            check_file_exist(user_extras[key])
        else:
            user_extras[key] = None

    #----------------------------------
    #             Columns
    #----------------------------------
    # separate out columns from user inlists

    user_columns = {}

    # separate out columns from user inlists
    for name in column_types:
        if name in user_mesa_inlists.keys():
            user_columns[name] = user_mesa_inlists[name]
            check_file_exist(user_columns[name])
        else:
            user_columns[name] = None

    return user_inlists, user_extras, user_columns


def resolve_configuration(keys, MESA_defaults, POSYDON_config, user_config, title="Configuration Priority"):
    """Resolve final configuration to use based on priority:
    user_config > POSYDON_config > MESA_defaults

    Parameters
    ----------
    keys : list
        List of keys to iterate over for resolution
    MESA_defaults : dict
        Dictionary of MESA default paths
    POSYDON_config : dict
        Dictionary of POSYDON configuration paths
    user_config : dict
        Dictionary of user configuration paths
    title : str, optional
        Title for logging output

    Returns
    -------
    final_config : dict
        Dictionary of final configuration paths to use
    """
    final_config = {}

    for key in keys:
        if user_config.get(key) is not None:
            final_config[key] = user_config[key]
        elif POSYDON_config.get(key) is not None:
            final_config[key] = POSYDON_config[key]
        else:
            final_config[key] = MESA_defaults.get(key)

    # Log at DEBUG level with detailed table
    if logger.isEnabledFor(logging.DEBUG):
        print_priority_table(keys, MESA_defaults, POSYDON_config, user_config, final_config, title)

    return final_config


def resolve_columns(MESA_default_columns, POSYDON_columns, user_columns):
    """Resolve final columns to use based on priority:
    user_columns > POSYDON_columns > MESA_default_columns

    Parameters
    ----------
    MESA_default_columns : dict
        Dictionary of MESA default column files paths
    POSYDON_columns : dict
        Dictionary of POSYDON column files paths
    user_columns : dict
        Dictionary of user column files paths

    Returns
    -------
    final_columns : dict
        Dictionary of final column files paths to use
    """
    final_columns = resolve_configuration(column_types, MESA_default_columns,
                                          POSYDON_columns, user_columns,
                                          title="Column Files Priority")

    # Log at INFO level which layer is used
    logger.info(f"{BOLD}Column Files:{RESET}")
    for name, filename in column_types.items():
        if user_columns.get(name) is not None:
            logger.info(f"  {name}: {MAGENTA}user{RESET}")
        elif POSYDON_columns.get(name) is not None:
            logger.info(f"  {name}: {YELLOW}POSYDON{RESET}")
        else:
            logger.info(f"  {name}: {CYAN}MESA{RESET}")

    return final_columns


def resolve_extras(MESA_default_extras, POSYDON_extras, user_extras):
    """Resolve final extras to use based on priority:
    user_extras > POSYDON_extras > MESA_default_extras

    Parameters
    ----------
    MESA_default_extras : dict
        Dictionary of MESA default extras paths
    POSYDON_extras : dict
        Dictionary of POSYDON extras paths
    user_extras : dict
        Dictionary of user extras paths

    Returns
    -------
    final_extras : dict
        Dictionary of final extras paths to use
    """
    final_extras = resolve_configuration(extras_keys, MESA_default_extras,
                                         POSYDON_extras, user_extras,
                                         title="EXTRAS Files Priority")

    # Log at INFO level which layer is used
    logger.info(f"{BOLD}EXTRAS Files:{RESET}")
    for key in extras_keys:
        if user_extras.get(key) is not None:
            logger.info(f"  {key}: {MAGENTA}user{RESET}")
        elif POSYDON_extras.get(key) is not None:
            logger.info(f"  {key}: {YELLOW}POSYDON{RESET}")
        else:
            logger.info(f"  {key}: {CYAN}MESA{RESET}")

    return final_extras


def print_priority_table(keys, MESA_defaults, POSYDON_config, user_config, final_config, title="Configuration Priority"):
    """Log a visual table showing which configuration layer is used for each key.

    Parameters
    ----------
    keys : list
        List of keys to display
    MESA_defaults : dict
        Dictionary of MESA default paths
    POSYDON_config : dict
        Dictionary of POSYDON configuration paths
    user_config : dict
        Dictionary of user configuration paths
    final_config : dict
        Dictionary of final resolved configuration paths
    title : str, optional
        Title for the table
    """

    # Find the longest key name for formatting
    max_key_len = max(len(str(key)) for key in keys)
    col_width = 12

    # Print title
    logger.debug(f"{BOLD}{title}{RESET}")
    logger.debug("=" * (max_key_len + col_width * 3 + 4))

    # Print header with colored column names
    header = f"{'Key':<{max_key_len}}  {CYAN}{'MESA':^{col_width}}{RESET}{YELLOW}{'POSYDON':^{col_width}}{RESET}{MAGENTA}{'user':^{col_width}}{RESET}"
    logger.debug(f"{BOLD}{header}{RESET}")
    logger.debug("-" * (max_key_len + col_width * 3 + 4))

    # Print each row
    for key in keys:
        # Check which configs have this key
        has_mesa = MESA_defaults.get(key) is not None
        has_posydon = POSYDON_config.get(key) is not None
        has_user = user_config.get(key) is not None

        # Determine which one is used (priority: user > POSYDON > MESA)
        used = 'user' if has_user else ('POSYDON' if has_posydon else ('MESA' if has_mesa else None))

        # Format each column
        mesa_mark = 'x' if has_mesa else ' '
        posydon_mark = 'x' if has_posydon else ' '
        user_mark = 'x' if has_user else ' '

        # Apply colors
        if used == 'MESA':
            mesa_str = f"{GREEN}{mesa_mark}{RESET}"
            posydon_str = f"{GRAY}{posydon_mark}{RESET}"
            user_str = f"{GRAY}{user_mark}{RESET}"
        elif used == 'POSYDON':
            mesa_str = f"{GRAY}{mesa_mark}{RESET}"
            posydon_str = f"{GREEN}{posydon_mark}{RESET}"
            user_str = f"{GRAY}{user_mark}{RESET}"
        elif used == 'user':
            mesa_str = f"{GRAY}{mesa_mark}{RESET}"
            posydon_str = f"{GRAY}{posydon_mark}{RESET}"
            user_str = f"{GREEN}{user_mark}{RESET}"
        else:
            mesa_str = f"{GRAY}{mesa_mark}{RESET}"
            posydon_str = f"{GRAY}{posydon_mark}{RESET}"
            user_str = f"{GRAY}{user_mark}{RESET}"

        # Log row
        logger.debug(f"{key:<{max_key_len}}  {mesa_str:^{col_width+9}}{posydon_str:^{col_width+9}}{user_str:^{col_width+9}}")

    logger.debug("=" * (max_key_len + col_width * 3 + 4))
    logger.debug(f"{GREEN}Green{RESET} = used, {GRAY}Gray{RESET} = available but not used.")

def print_inlist_stacking_table(keys, MESA_defaults, POSYDON_config, user_config, title="Inlist Stacking"):
    """Print a visual table showing how inlists are stacked for each key.

    Unlike extras and columns which use replacement, inlists stack together where
    each layer contributes files. The final inlist for each key contains entries
    from user config (highest priority), then POSYDON config, then MESA defaults (lowest priority).

    Parameters
    ----------
    keys : list
        List of inlist keys to display
    MESA_defaults : dict
        Dictionary of MESA default inlist paths (each value is a list)
    POSYDON_config : dict
        Dictionary of POSYDON configuration inlist paths (each value is a list)
    user_config : dict
        Dictionary of user configuration inlist paths (each value is a list)
    title : str, optional
        Title for the table
    """
    # Find the longest key name for formatting
    max_key_len = max(len(str(key)) for key in keys)

    # Print title
    print(f"{BOLD}{title}{RESET}")
    print("=" * 80)
    print(f"{BOLD}{'Key':<{max_key_len}}  Layer Stack{RESET}")
    print("-" * 80)

    # Print each row
    for key in keys:
        mesa_list = MESA_defaults.get(key, []) or []
        posydon_list = POSYDON_config.get(key, []) or []
        user_list = user_config.get(key, []) or []

        # Count total files
        total = len(mesa_list) + len(posydon_list) + len(user_list)

        if total == 0:
            print(f"{key:<{max_key_len}}  {GRAY}(no inlists){RESET}")
            continue

        # Print key with count
        print(f"{BOLD}{key:<{max_key_len}}{RESET}  ({total} file{'s' if total != 1 else ''})")

        # Print each layer with indentation (highest priority first)
        layer_num = 1

        for inlist_path in user_list:
            filename = os.path.basename(inlist_path)
            print(f"{'':>{max_key_len}}  {MAGENTA}[{layer_num}]{RESET} {GRAY}user:{RESET}    {filename}")
            layer_num += 1

        for inlist_path in posydon_list:
            filename = os.path.basename(inlist_path)
            print(f"{'':>{max_key_len}}  {YELLOW}[{layer_num}]{RESET} {GRAY}POSYDON:{RESET} {filename}")
            layer_num += 1

        for inlist_path in mesa_list:
            filename = os.path.basename(inlist_path)
            print(f"{'':>{max_key_len}}  {CYAN}[{layer_num}]{RESET} {GRAY}MESA:{RESET}    {filename}")
            layer_num += 1

        print()  # Blank line between keys

    print("=" * 80)
    print(f"Note: Files at the top override parameters from files below")
    print(f"{MAGENTA}user{RESET} (highest priority) → {YELLOW}POSYDON{RESET} (config) → {CYAN}MESA{RESET} (base)")


def print_inlist_parameter_override_table(key, layer_params, final_params, show_details=False):
    """Log a table showing which layer each parameter comes from (supports all layers).

    Parameters
    ----------
    key : str
        The inlist key (e.g., 'binary_controls', 'star1_controls')
    layer_params : dict
        Dictionary mapping layer names to parameter dictionaries for this key
        Format: {'MESA': {params}, 'POSYDON': {params}, 'user': {params}, 'grid': {params}, 'output': {params}}
    final_params : dict
        Final merged parameters
    show_details : bool, optional
        If True, show detailed parameter-by-parameter table. Default is False.
    """
    if not show_details:
        return

    # Get parameters for each layer for this specific key
    mesa_params = layer_params.get('MESA', {}).get(key, {})
    posydon_params = layer_params.get('POSYDON', {}).get(key, {})
    user_params = layer_params.get('user', {}).get(key, {})
    grid_params = layer_params.get('grid', {}).get(key, {})
    output_params = layer_params.get('output', {}).get(key, {})

    # Get all unique parameter names
    all_params = sorted(set(
        list(mesa_params.keys()) +
        list(posydon_params.keys()) +
        list(user_params.keys()) +
        list(grid_params.keys()) +
        list(output_params.keys())
    ))

    if not all_params:
        return

    # Find the longest parameter name for formatting
    max_param_len = max(len(str(param)) for param in all_params)
    max_param_len = max(max_param_len, 15)  # Minimum width
    col_width = 10

    # Color mapping for each layer
    layer_colors = {
        'MESA': CYAN,
        'POSYDON': YELLOW,
        'user': MAGENTA,
        'grid': '\033[94m',      # Blue
        'output': '\033[92m'      # Green
    }

    # Determine which layers have parameters for this key
    active_layers = []
    for layer_name in ['MESA', 'POSYDON', 'user', 'grid', 'output']:
        layer_data = layer_params.get(layer_name, {}).get(key, {})
        if layer_data:
            active_layers.append(layer_name)

    # Count overridden parameters
    overridden_count = 0
    for param in all_params:
        count = sum(1 for layer_name in active_layers
                   if param in layer_params.get(layer_name, {}).get(key, {}))
        if count > 1:
            overridden_count += 1

    # Log summary header
    logger.debug(f"  {BOLD}Detailed Parameters:{RESET} {len(all_params)} total, {overridden_count} overridden")
    total_width = max_param_len + col_width * len(active_layers) + len(active_layers) * 2
    logger.debug("  " + "=" * total_width)

    # Log header with colored column names
    header_parts = [f"{'Parameter':<{max_param_len}}"]
    for layer_name in active_layers:
        color = layer_colors.get(layer_name, RESET)
        header_parts.append(f"{color}{layer_name:^{col_width}}{RESET}")
    header = "  ".join(header_parts)
    logger.debug(f"  {BOLD}{header}{RESET}")
    logger.debug("  " + "-" * total_width)

    # Log each parameter row
    for param in all_params:
        # Determine which layer provides the final value (check in priority order: highest to lowest)
        # Priority: output > grid > user > POSYDON > MESA
        used_layer = None
        for layer_name in ['output', 'grid', 'user', 'POSYDON', 'MESA']:
            if param in layer_params.get(layer_name, {}).get(key, {}):
                used_layer = layer_name
                break  # Found the highest priority layer with this parameter

        row_parts = [f"{param:<{max_param_len}}"]
        for layer_name in active_layers:
            has_param = param in layer_params.get(layer_name, {}).get(key, {})
            mark = 'x' if has_param else ' '
            color = layer_colors.get(layer_name, RESET)

            # Green if this layer provides the final value, gray if available but not used
            if has_param and layer_name == used_layer:
                row_parts.append(f"{GREEN}{mark:^{col_width}}{RESET}")
            elif has_param:
                row_parts.append(f"{GRAY}{mark:^{col_width}}{RESET}")
            else:
                row_parts.append(f"{mark:^{col_width}}")

        logger.debug("  " + "  ".join(row_parts))

    logger.debug("  " + "=" * total_width)
    logger.debug(f"  {GREEN}Green{RESET} = used, {GRAY}Gray{RESET} = available but not used")




def print_inlist_summary_table_v2(all_keys, layer_counts):
    """Log a summary table showing parameter counts per section at each layer.

    This version supports multiple layers including grid and output configurations.
    Logs at INFO level.

    Parameters
    ----------
    all_keys : list
        List of inlist keys (sections)
    layer_counts : dict
        Dictionary mapping layer names to count dictionaries
        Format: {'MESA': {key: count}, 'POSYDON': {key: count}, ...}
    """
    # Color mapping for each layer
    layer_colors = {
        'MESA': CYAN,
        'POSYDON': YELLOW,
        'user': MAGENTA,
        'grid': '\033[94m',      # Blue
        'output': '\033[92m'      # Green (reuse)
    }

    # Find the longest key name for formatting
    max_key_len = max(len(str(key)) for key in all_keys)
    max_key_len = max(max_key_len, 15)  # Minimum width
    col_width = 10

    # Get layers that have any parameters
    active_layers = [layer for layer in ['MESA', 'POSYDON', 'user', 'grid', 'output']
                     if layer in layer_counts and any(layer_counts[layer].values())]

    num_layers = len(active_layers)
    total_width = max_key_len + col_width * num_layers + (num_layers + 4)

    logger.info(f"{BOLD}Parameter Count Summary{RESET}")
    logger.info("=" * total_width)

    # Print header with colored column names
    header_parts = [f"{'Section':<{max_key_len}}"]
    for layer in active_layers:
        color = layer_colors.get(layer, RESET)
        header_parts.append(f"{color}{layer:^{col_width}}{RESET}")
    header = "  ".join(header_parts)
    logger.info(f"{BOLD}{header}{RESET}")
    logger.info("-" * total_width)

    # Print each section row
    for key in all_keys:
        row_parts = [f"{key:<{max_key_len}}"]
        for layer in active_layers:
            count = layer_counts[layer].get(key, 0)
            color = layer_colors.get(layer, RESET)
            row_parts.append(f"{color}{count:^{col_width}}{RESET}")
        logger.info("  ".join(row_parts))

    logger.info("=" * total_width)
    logger.info(f"Layer priority (lowest → highest): {' → '.join(active_layers)}")

def _get_section_from_key(key):
    """Determine the MESA inlist section based on the key name.

    Parameters
    ----------
    key : str
        The inlist key name

    Returns
    -------
    str or None
        The section name (e.g., '&binary_controls') or None
    """
    if 'binary_controls' in key:
        return '&binary_controls'
    elif 'binary_job' in key:
        return '&binary_job'
    elif ('star1_job' in key) or ('star2_job' in key):
        return '&star_job'
    elif ('star1_controls' in key) or ('star2_controls' in key):
        return '&controls'
    else:
        return None

def _process_inlist_layer(inlist_paths, section):
    """Process a single layer of inlist files and return merged parameters.
    Files are processed in order, with later files overriding parameters from
    earlier files.

    Parameters
    ----------
    inlist_paths : list or None
        List of file paths for this layer
    section : str or None
        The MESA section to extract

    Returns
    -------
    dict
        Merged parameters from all files in this layer
    """
    layer_params = {}
    if inlist_paths:
        for file_path in inlist_paths:
            inlist_dict = utils.clean_inlist_file(file_path, section=section)[section]
            layer_params.update(inlist_dict)
    return layer_params

def _clean_inlist_parameters(params_dict):
    """Clean inlist parameters by removing unwanted keys and replacing placeholders.

    Parameters
    ----------
    params_dict : dict
        Dictionary of parameters to clean
    key : str
        The inlist key (used for logging)

    Returns
    -------
    dict
        Cleaned parameters dictionary
    """
    # Remove read_extra and inlist references (any parameter containing these substrings)
    # This catches: read_extra, inlist, read_extra_controls_inlist1, inlist_names, etc.
    cleaned = {k: v for k, v in params_dict.items()
               if not any(substring in k for substring in ['read_extra', 'inlist'])}

    # Replace num_x_ctrls with actual index (e.g., num_x_ctrls -> 1)
    # This is a special default that the default value in the .defaults
    # file in MESA does not work because it is a placeholder
    keys_to_replace = {k: k.replace('num_x_ctrls', '1')
                      for k in cleaned.keys()
                      if 'num_x_ctrls' in k}

    if keys_to_replace:
        for old_key, new_key in keys_to_replace.items():
            cleaned[new_key] = cleaned.pop(old_key)

    return cleaned

def _build_grid_parameter_layer(grid_parameters, final_inlists):
    """Build a layer of inlist parameters for grid-specific configurations.

    Parameters
    ----------
    grid_parameters : list or set
        Collection of grid parameter names
    final_inlists : dict
        Current state of final inlists (to check which sections are affected)

    Returns
    -------
    dict
        Dictionary mapping section names to parameter dictionaries for grid config
    """
    grid_layer = {}

    # Configuration for each section: (read_extra_param, extra_name_param, inlist_filename)
    section_config = {
        'star1_controls': ('read_extra_controls_inlist1',
                          'extra_controls_inlist1_name',
                          'inlist_grid_star1_binary_controls'),
        'star2_controls': ('read_extra_controls_inlist1',
                          'extra_controls_inlist1_name',
                          'inlist_grid_star2_binary_controls'),
        'star1_job': ('read_extra_star_job_inlist1',
                     'extra_star_job_inlist1_name',
                     'inlist_grid_star1_job'),
        'star2_job': ('read_extra_star_job_inlist1',
                     'extra_star_job_inlist1_name',
                     'inlist_grid_star2_job'),
    }

    # Check which sections have grid parameters
    for section, config in section_config.items():
        matching_params = [param for param in grid_parameters
                          if param in final_inlists.get(section, {})]

        if matching_params:
            read_param, name_param, filename = config
            grid_layer[section] = {
                read_param: '.true.',
                name_param: f"'{filename}'"
            }
            # Log at DEBUG level which sections are affected by grid parameters
            logger.debug(f"  Grid parameters affecting {section}: {', '.join(matching_params)}")
        else:
            grid_layer[section] = {}

    # Log at DEBUG level summary
    affected_sections = [sec for sec, params in grid_layer.items() if params]
    if affected_sections:
        logger.debug(f"Grid parameters affect sections: {', '.join(affected_sections)}")

    return grid_layer


def _build_output_controls_layer(output_settings):
    """Build a layer of inlist parameters for output control configurations.

    Parameters
    ----------
    output_settings : dict
        Dictionary of output settings from the configuration file

    Returns
    -------
    dict
        Dictionary mapping section names to parameter dictionaries for output controls
    """
    # Convert boolean to MESA Fortran string
    def to_fortran_bool(value):
        return ".true." if value else ".false."

    output_layer = {
        'binary_controls': {},
        'binary_job': {},
        'star1_controls': {},
        'star1_job': {},
        'star2_controls': {},
        'star2_job': {}
    }

    # Configuration: (config_key, section, enabled_param, filename_param, filename_value)
    output_config = [
        ('final_profile_star1', 'star1_job', 'write_profile_when_terminate',
         'filename_for_profile_when_terminate', "'final_profile_star1.data'"),
        ('final_profile_star2', 'star2_job', 'write_profile_when_terminate',
         'filename_for_profile_when_terminate', "'final_profile_star2.data'"),
        ('final_model_star1', 'star1_job', 'save_model_when_terminate',
         'save_model_filename', "'final_star1.mod'"),
        ('final_model_star2', 'star2_job', 'save_model_when_terminate',
         'save_model_filename', "'final_star2.mod'"),
        ('history_star1', 'star1_controls', 'do_history_file', None, None),
        ('history_star2', 'star2_controls', 'do_history_file', None, None),
    ]

    # Process each output configuration
    for config_key, section, enabled_param, filename_param, filename_value in output_config:
        if config_key in output_settings:
            is_enabled = output_settings[config_key]
            output_layer[section][enabled_param] = to_fortran_bool(is_enabled)

            # Set filename parameter if provided and feature is enabled
            if filename_param and is_enabled:
                output_layer[section][filename_param] = filename_value

    # Handle history_interval (applies to all sections)
    if 'history_interval' in output_settings:
        interval = output_settings['history_interval']
        output_layer['binary_controls']['history_interval'] = interval
        output_layer['star1_controls']['history_interval'] = interval
        output_layer['star2_controls']['history_interval'] = interval

    # Disable binary history if requested
    if 'binary_history' in output_settings and not output_settings['binary_history']:
        output_layer['binary_controls']['history_interval'] = "-1"

    # Log at DEBUG level which sections have output control parameters
    affected_sections = [sec for sec, params in output_layer.items() if params]
    if affected_sections:
        logger.debug(f"Output control parameters affect sections: {', '.join(affected_sections)}")
        for section in affected_sections:
            params = output_layer[section]
            if params:
                param_list = ', '.join(params.keys())
                logger.debug(f"  {section}: {param_list}")

    # Handle ZAMS filenames if provided
    if 'zams_filename_1' in output_settings and output_settings['zams_filename_1'] is not None:
       output_layer['star1_controls']['zams_filename'] = f"'{output_settings['zams_filename_1']}'"

    if 'zams_filename_2' in output_settings and output_settings['zams_filename_2'] is not None:
       output_layer['star2_controls']['zams_filename'] = f"'{output_settings['zams_filename_2']}'"

    return output_layer


def resolve_inlists(MESA_default_inlists, POSYDON_inlists,
                    user_inlists, system_type, run_directory, grid_parameters=None,
                    output_settings=None):
    """Resolve final inlists to use based on priority:
    output_settings > grid_parameters > user_inlists > POSYDON_inlists > MESA_default_inlists

    The inlists are stacked lists, so the final inlist for each key
    contains all layers in priority order.

    Examples
    -------
    MESA_default_inlist contains
    - `alpha_semiconvection = 0.1d0`
    POSYDON_inlist contains
    - `alpha_semiconvection = 0.2d0`
    user_inlist contains
    - `alpha_semiconvection = 0.15d0`
    The final inlist will contain
    - `alpha_semiconvection = 0.15d0`

    Parameters
    ----------
    MESA_default_inlists : dict
        Dictionary of MESA default inlists paths
    POSYDON_inlists : dict
        Dictionary of POSYDON inlists paths
    user_inlists : dict
        Dictionary of user inlists paths
    system_type : str
        Type of binary system
    run_directory : str
        Path to the run directory where inlist files will be created
    grid_parameters : list or set, optional
        Collection of grid parameter names. If provided, adds grid configuration layer.
    output_settings : dict, optional
        Dictionary of output settings. If provided, adds output control layer.

    Returns
    -------
    final_inlists : dict
        Dictionary where each key maps to parameter dictionaries
    """
    # Get all unique keys from all dictionaries
    all_keys = sorted(set(list(MESA_default_inlists.keys()) +
                   list(POSYDON_inlists.keys()) +
                   list(user_inlists.keys())))

    # Stack inlists: MESA (base) → POSYDON → user → grid → output (highest priority)
    final_inlists = {}

    # Track parameter counts for summary
    layer_counts = {
        'MESA': {},
        'POSYDON': {},
        'user': {},
        'grid': {},
        'output': {},
        'inlist_names': {}
    }

    # Track layer parameters for detailed printing if requested
    layer_params = {
        'MESA': {},
        'POSYDON': {},
        'user': {},
        'grid': {},
        'output': {},
        'inlist_names': {}
    }

    # First pass: process file-based layers (MESA, POSYDON, user)
    for key in all_keys:
        # Determine the section based on the key name
        section = _get_section_from_key(key)
        if 'single' in system_type and ('binary' in key or 'star2' in key):
            # Skip binary or star2 sections for single star systems
            final_inlists[key] = {}
            layer_counts['MESA'][key] = 0
            layer_counts['POSYDON'][key] = 0
            layer_counts['user'][key] = 0
            layer_params['MESA'][key] = {}
            layer_params['POSYDON'][key] = {}
            layer_params['user'][key] = {}
            continue

        # Process each file-based layer
        mesa_layer_params = _process_inlist_layer(MESA_default_inlists.get(key), section)
        posydon_layer_params = _process_inlist_layer(POSYDON_inlists.get(key), section)
        user_layer_params = _process_inlist_layer(user_inlists.get(key), section)

        # Merge file-based layers (order matters: MESA first, then POSYDON, then user)
        final_inlists[key] = {}
        final_inlists[key].update(mesa_layer_params)
        final_inlists[key].update(posydon_layer_params)
        final_inlists[key].update(user_layer_params)

        # Store counts and parameters for summary
        layer_counts['MESA'][key] = len(mesa_layer_params)
        layer_counts['POSYDON'][key] = len(posydon_layer_params)
        layer_counts['user'][key] = len(user_layer_params)

        layer_params['MESA'][key] = mesa_layer_params
        layer_params['POSYDON'][key] = posydon_layer_params
        layer_params['user'][key] = user_layer_params

    # Clean the final inlist parameters.
    # Needs to happen before adding grid/output layers!
    # because grid_parameters add read_extra parameters that need to be preserved!
    for key in all_keys:
        final_inlists[key] = _clean_inlist_parameters(final_inlists[key])

    # Build grid configuration layer if provided
    if grid_parameters:
        grid_layer_dict = _build_grid_parameter_layer(grid_parameters, final_inlists)
        for key in all_keys:
            grid_params = grid_layer_dict.get(key, {})
            final_inlists[key].update(grid_params)
            layer_counts['grid'][key] = len(grid_params)
            layer_params['grid'][key] = grid_params
    else:
        for key in all_keys:
            layer_counts['grid'][key] = 0
            layer_params['grid'][key] = {}

    # Build output controls layer if provided
    if output_settings:
        output_layer_dict = _build_output_controls_layer(output_settings)
        for key in all_keys:
            output_params = output_layer_dict.get(key, {})
            final_inlists[key].update(output_params)
            layer_counts['output'][key] = len(output_params)
            layer_params['output'][key] = output_params
    else:
        for key in all_keys:
            layer_counts['output'][key] = 0
            layer_params['output'][key] = {}

    # Add inlist_names layer for binary systems
    # This must happen after all other layers to use the constructed run_directory paths
    if 'single' not in system_type.lower():
        inlist_star1_binary = os.path.join(run_directory, 'binary', 'inlist1')
        inlist_star2_binary = os.path.join(run_directory, 'binary', 'inlist2')

        inlist_names_params = {
            'inlist_names(1)': f"'{inlist_star1_binary}'",
            'inlist_names(2)': f"'{inlist_star2_binary}'"
        }
        final_inlists['binary_job'].update(inlist_names_params)

        # Track this layer
        for key in all_keys:
            if key == 'binary_job':
                layer_counts['inlist_names'][key] = len(inlist_names_params)
                layer_params['inlist_names'][key] = inlist_names_params
            else:
                layer_counts['inlist_names'][key] = 0
                layer_params['inlist_names'][key] = {}
    else:
        # Single star systems don't need inlist_names
        for key in all_keys:
            layer_counts['inlist_names'][key] = 0
            layer_params['inlist_names'][key] = {}

    # Log at INFO level: Parameter count summary
    print_inlist_summary_table_v2(all_keys, layer_counts)

    # Log at DEBUG level: Detailed parameter tables
    if logger.isEnabledFor(logging.DEBUG):
        for key in all_keys:
            # Only show sections that have parameters in any layer
            if any(layer_params[layer][key] for layer in layer_params):
                logger.debug(f"{BOLD}═══ {key} ═══{RESET}")
                print_inlist_parameter_override_table(
                    key,
                    layer_params,
                    final_inlists[key],
                    show_details=True
                )

    # Return the final parameter dictionaries
    return final_inlists


def read_grid_file(filepath):
    """Read grid file and return grid parameters as a dictionary.

    Parameters
    ----------
    filepath : str
        Path to the grid file

    Returns
    -------
    number of grid points : int
        Number of grid points in the grid
    grid parameter names : set
        Set of grid parameter names (column names)
    fixgrid_file_name : str
        Path to the fixed grid file used for the grid
    """
    # TODO: I"m not sure what the last option is for.
    # Is it processing runs for a grid stored in a directory?
    if '.csv' in filepath:
        grid_df = pd.read_csv(filepath)
        fixgrid_file_name = filepath
    elif '.h5' in filepath:
        psy_grid = PSyGrid()
        psy_grid.load(filepath)
        grid_df = psy_grid.get_pandas_initial_final()
        psy_grid.close()
        fixgrid_file_name = filepath
    elif os.path.isdir(filepath):
        PSyGrid().create(filepath, "./fixed_grid_results.h5", slim=True)
        psy_grid = PSyGrid()
        psy_grid.load("./fixed_grid_results.h5")
        grid_df = psy_grid.get_pandas_initial_final()
        psy_grid.close()
        fixgrid_file_name = os.path.join(os.getcwd(), "fixed_grid_results.h5")
    else:
        raise ValueError('Grid format not recognized, please feed in an acceptable format: csv')

    grid_parameters = set(grid_df.columns.tolist())
    grid_parameters = [params.lower() for params in grid_parameters]

    return (len(grid_df), grid_parameters, fixgrid_file_name)

def get_additional_user_settings(mesa_inlists):
    """Extract additional settings from configuration for output controls.

    Parameters
    ----------
    mesa_inlists : dict
        Dictionary of user inlists and configuration

    Returns
    -------
    dict
        Dictionary of additional settings
    """
    output_settings_list = ['history_interval', 'binary_history',
                            'final_profile_star1', 'final_profile_star2',
                            'final_model_star1', 'final_model_star2',
                            'history_star1', 'history_star2',
                            'zams_filename_1', 'zams_filename_2']
    output_settings = {}

    for setting in output_settings_list:
        if setting in mesa_inlists:
            output_settings[setting] = mesa_inlists[setting]
        else:
            output_settings[setting] = None

    # Handle single zams_filename case
    if 'zams_filename' in mesa_inlists:
        output_settings['zams_filename_1'] = mesa_inlists['zams_filename']
        output_settings['zams_filename_2'] = mesa_inlists['zams_filename']


    return output_settings



def _write_inlist_section(f, section_name, parameters):
    """Write a single MESA inlist section to file.

    Parameters
    ----------
    f : file object
        File handle opened in binary write mode
    section_name : str
        Name of the section (e.g., 'controls', 'binary_controls')
    parameters : dict
        Dictionary of parameter name -> value pairs
    """
    f.write(f'&{section_name}\n\n'.encode('utf-8'))
    for key, value in parameters.items():
        f.write(f'\t{key} = {value}\n'.encode('utf-8'))
    f.write(f'\n/ ! end of {section_name} namelist\n'.encode('utf-8'))


def _write_binary_inlist(filepath, binary_controls, binary_job):
    """Write the binary inlist_project file.

    Parameters
    ----------
    filepath : str
        Path to write the inlist file
    binary_controls : dict
        Parameters for binary_controls section
    binary_job : dict
        Parameters for binary_job section
    """
    with open(filepath, 'wb') as f:
        _write_inlist_section(f, 'binary_controls', binary_controls)
        f.write(b'\n')
        _write_inlist_section(f, 'binary_job', binary_job)
    logger.debug(f'Wrote inlist: {filepath}')


def _write_star_inlist(filepath, star_controls, star_job):
    """Write a star inlist file (for star1 or star2).

    Parameters
    ----------
    filepath : str
        Path to write the inlist file
    star_controls : dict
        Parameters for controls section
    star_job : dict
        Parameters for star_job section
    """
    with open(filepath, 'wb') as f:
        _write_inlist_section(f, 'controls', star_controls)
        f.write(b'\n')
        _write_inlist_section(f, 'star_job', star_job)
    logger.debug(f'Wrote inlist: {filepath}')


def _create_build_script(path):
    """Create the 'mk' build script for compiling MESA executables.

    Parameters
    ----------
    path : str
        Path to the grid run folder
    """
    mk_filepath = os.path.join(path, 'mk')
    if os.path.exists(mk_filepath):
        logger.warning(f"'mk' file already exists at {mk_filepath}. It will be overwritten.")

    with open(mk_filepath, 'w') as f:
        f.write(f'cd {os.path.join(path, "binary/make")}\n')
        f.write(f'make -f makefile_binary\n')
        f.write(f'cd {os.path.join(path, "star1/make")}\n')
        f.write(f'make -f makefile_star\n')
        f.write(f'cd {os.path.join(path, "star2/make")}\n')
        f.write(f'make -f makefile_star\n')

    logger.debug(f'Created build script: {mk_filepath}')
    subprocess.run(['chmod', '755', mk_filepath])
    out = subprocess.run(['./mk'], shell=True, cwd=path, capture_output=True, text=True)
    if out.returncode != 0:
        logger.error(f"Building the MESA executables has failed!")
        for line in out.stderr.strip().split('\n'):
            logger.error(line)

def _copy_columns(path, final_columns):
    """Copy column list files to the grid run folder.

    Parameters
    ----------
    path : str
        Path to the grid run folder
    final_columns : dict
        Dictionary mapping column types to file paths

    Returns
    -------
    dict
        Dictionary mapping column types to destination paths in the grid run folder
    """

    logger.debug(f'{BOLD}COLUMN LISTS USED:{RESET}')
    out_paths = {}
    for key, value in final_columns.items():
        logger.debug(f"{key}:  {value}")
        dest = os.path.join(path, 'column_lists', column_types[key])
        shutil.copy(value, dest)
        out_paths[key] = dest
    return out_paths

def _copy_extras(path, final_extras):
    """Copy extras files (makefiles, run files) to the grid run folder.

    Parameters
    ----------
    path : str
        Path to the grid run folder
    final_extras : dict
        Dictionary mapping extras keys to file paths
    """
    logger.debug(f'{BOLD}EXTRAS USED:{RESET}')

    # Define destination mapping for each extras key
    extras_destinations = {
        'makefile_binary': [('binary/make', 'makefile_binary')],
        'makefile_star': [('star1/make', 'makefile_star'),
                          ('star2/make', 'makefile_star')],
        'binary_run': [('binary/src', 'binary_run.f')],
        'star_run': [('star1/src', 'run.f'),
                     ('star2/src', 'run.f')],
        'run_binary_extras': [('binary/src', 'run_binary_extras.f')],
        'run_star_binary_extras': [('binary/src', 'run_star_extras.f')],
        'run_star1_extras': [('star1/src', 'run_star_extras.f')],
        'run_star2_extras': [('star2/src', 'run_star_extras.f')],
    }

    for key, value in final_extras.items():
        logger.debug(f"{key}: {value}")
        if key == 'mesa_dir':
            continue

        if key in extras_destinations:
            for subdir, filename in extras_destinations[key]:
                dest = os.path.join(path, subdir, filename)
                shutil.copy(value, dest)
        else:
            logger.warning(f"Unrecognized extras key '{key}'. Copying to root.")
            shutil.copy(value, path)


def setup_grid_run_folder(path, final_columns, final_extras,
                          final_inlists):
    """Set up the grid run folder by:

    1. Creating necessary subdirectories
    2. Copying columns and extras into the right folders
    3. Writing the inlist files
    4. Building MESA executables

    Parameters
    ----------
    path : str
        Path to the grid run folder
    final_columns : dict
        Dictionary mapping column types to file paths
    final_extras : dict
        Dictionary mapping extras keys to file paths
    final_inlists : dict
        Dictionary mapping inlist keys to parameter dictionaries
    """
    # Create directory structure
    subdirs = ['binary', 'binary/make', 'binary/src',
               'star1', 'star1/make', 'star1/src',
               'star2', 'star2/make', 'star2/src',
               'column_lists']

    for subdir in subdirs:
        dir_path = os.path.join(path, subdir)
        os.makedirs(dir_path, exist_ok=True)

    logger.debug(f"Setting up grid run folder at: {path}")

    # Copy columns and extras
    column_paths = _copy_columns(path, final_columns)
    _copy_extras(path, final_extras)

    # Create and run build script
    _create_build_script(path)

    # Write inlist files
    logger.debug(f'{BOLD}Writing MESA inlist files:{RESET}')
    inlist_binary_project = os.path.join(path, 'binary', 'inlist_project')
    inlist_star1_binary = os.path.join(path, 'binary', 'inlist1')
    inlist_star2_binary = os.path.join(path, 'binary', 'inlist2')

    # inlist_names(1) and inlist_names(2) are now set in resolve_inlists
    # for binary systems, so no need to add them here

    # Write all three inlist files
    _write_binary_inlist(inlist_binary_project,
                         final_inlists['binary_controls'],
                         final_inlists['binary_job'])

    _write_star_inlist(inlist_star1_binary,
                       final_inlists['star1_controls'],
                       final_inlists['star1_job'])

    _write_star_inlist(inlist_star2_binary,
                       final_inlists['star2_controls'],
                       final_inlists['star2_job'])

    logger.info('MESA inlist files written successfully.')

    # Essentials paths created in this functions
    output_paths = {
        'binary_executable': os.path.join(path, 'binary', 'binary'),
        'star1_executable': os.path.join(path, 'star1', 'star'),
        'star2_executable': os.path.join(path, 'star2', 'star'),
        'inlist_binary_project': inlist_binary_project,
        'inlist_star1_binary': inlist_star1_binary,
        'inlist_star2_binary': inlist_star2_binary,
    }
    output_paths.update(column_paths)

    return output_paths



def _write_environment_setup(f, slurm):
    """Write environment setup commands to a script file.

    Parameters
    ----------
    f : file object
        Open file to write to
    slurm : dict
        SLURM configuration dictionary
    """
    f.write(f'export OMP_NUM_THREADS={slurm["number_of_cpus_per_task"]}\n\n')
    f.write(f'export MESASDK_ROOT={os.environ["MESASDK_ROOT"]}\n')
    f.write('source $MESASDK_ROOT/bin/mesasdk_init.sh\n')
    f.write(f'export MESA_DIR={os.environ["MESA_DIR"]}\n\n')


def _write_sbatch_header(f, slurm, job_type='grid', array_size=None):
    """Write SBATCH header directives to a script file.

    Parameters
    ----------
    f : file object
        Open file to write to
    slurm : dict
        SLURM configuration dictionary
    job_type : str
        Type of job ('grid', 'cleanup')
    array_size : int or None
        If provided, creates job array with this size
    """
    f.write('#!/bin/bash\n')
    f.write(f'#SBATCH --account={slurm["account"]}\n')
    f.write(f'#SBATCH --partition={slurm["partition"]}\n')

    if job_type == 'cleanup':
        f.write('#SBATCH -N 1\n')
        f.write('#SBATCH --cpus-per-task 1\n')
        f.write('#SBATCH --ntasks-per-node 1\n')
        f.write(f'#SBATCH --time={slurm["walltime"]}\n')
        f.write('#SBATCH --job-name="mesa_grid_cleanup"\n')
        f.write('#SBATCH --output=mesa_cleanup.out\n')
    elif array_size is not None:  # Job array mode
        f.write('#SBATCH -N 1\n')
        f.write(f'#SBATCH --array=0-{array_size-1}\n')
        f.write(f'#SBATCH --cpus-per-task {slurm["number_of_cpus_per_task"]}\n')
        f.write('#SBATCH --ntasks-per-node 1\n')
        f.write(f'#SBATCH --time={slurm["walltime"]}\n')
        f.write('#SBATCH --job-name="mesa_grid_${SLURM_ARRAY_TASK_ID}"\n')
        f.write('#SBATCH --output=mesa_grid.%A_%a.out\n')
    else:  # MPI mode
        f.write(f'#SBATCH -N {slurm["number_of_nodes"]}\n')
        f.write(f'#SBATCH --cpus-per-task {slurm["number_of_cpus_per_task"]}\n')
        f.write(f'#SBATCH --ntasks-per-node {slurm["number_of_mpi_tasks"]}\n')
        f.write(f'#SBATCH --time={slurm["walltime"]}\n')
        f.write('#SBATCH --output="mesa_grid.out"\n')

    f.write('#SBATCH --mail-type=ALL\n')
    f.write(f'#SBATCH --mail-user={slurm["email"]}\n')
    f.write('#SBATCH --mem-per-cpu=4G\n\n')


def _write_cleanup_commands(f, slurm):
    """Write cleanup commands to a script file.

    Parameters
    ----------
    f : file object
        Open file to write to
    slurm : dict
        SLURM configuration dictionary
    """
    f.write('compress-mesa .\n')
    if 'newgroup' in slurm:
        f.write(f'echo "Change group to {slurm["newgroup"]}"\n')
        f.write(f'chgrp -fR {slurm["newgroup"]} .\n')
        f.write('echo "Change group permission to rwX at least"\n')
        f.write('chmod -fR g+rwX .\n')
    f.write('\necho "Done."')


def generate_submission_scripts(submission_type, command_line, slurm, nr_systems):
    """Generate submission scripts (shell or SLURM) for running the grid.

    Parameters
    ----------
    submission_type : str
        Type of submission script to generate ('shell' or 'slurm')
    command_line : str
        The base command line to execute
    slurm : dict
        Dictionary containing SLURM configuration parameters
    nr_systems : list
        List of systems in the grid (used for job array length)
    """
    if submission_type == 'shell':
        script_name = 'grid_command.sh'
        if os.path.exists(script_name):
            logger.warning(f'Replacing existing script: {script_name}')

        with open(script_name, 'w') as f:
            f.write('#!/bin/bash\n\n')
            _write_environment_setup(f, slurm)

            # Job array loop if needed
            if slurm['job_array']:
                indices = ' '.join(str(i) for i in range(nr_systems))
                f.write(f'for SLURM_ARRAY_TASK_ID in {indices}; do ')

            f.write(command_line)

            if slurm['job_array']:
                f.write(' ; done\n\n')
            else:
                f.write('\n\n')

            _write_cleanup_commands(f, slurm)

        os.system(f"chmod 755 {script_name}")
        logger.debug(f"Created {os.path.abspath(script_name)}")

    elif submission_type == 'slurm':
        # Generate main grid submission script
        grid_script = 'job_array_grid_submit.slurm' if slurm['job_array'] else 'mpi_grid_submit.slurm'
        if os.path.exists(grid_script):
            logger.warning(f'Replace {grid_script}')

        array_size = nr_systems if slurm['job_array'] else None
        with open(grid_script, 'w') as f:
            _write_sbatch_header(f, slurm, job_type='grid', array_size=array_size)
            _write_environment_setup(f, slurm)
            f.write(command_line)

        # Generate cleanup script
        cleanup_script = 'cleanup.slurm'
        if os.path.exists(cleanup_script):
            logger.warning(f'Replace {cleanup_script}')

        with open(cleanup_script, 'w') as f:
            _write_sbatch_header(f, slurm, job_type='cleanup')
            _write_cleanup_commands(f, slurm)

        # Generate wrapper run script
        run_script = 'run_grid.sh'
        if os.path.exists(run_script):
            logger.warning(f'Replace {run_script}')
        with open(run_script, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write(f'ID_GRID=$(sbatch --parsable {grid_script})\n')
            f.write(f'echo "{grid_script} submitted as "${{ID_GRID}}\n')
            f.write('ID_cleanup=$(sbatch --parsable --dependency=afterany:${ID_GRID} '
                    '--kill-on-invalid-dep=yes cleanup.slurm)\n')
            f.write('echo "cleanup.slurm submitted as "${ID_cleanup}\n')

        os.system(f"chmod 755 {run_script}")
        logger.info(f"Created {grid_script}, {cleanup_script}, and {run_script}")


def construct_command_line(number_of_mpi_processes, path_to_grid,
                           binary_exe, star1_exe, star2_exe,
                           inlist_binary_project, inlist_star1_binary, inlist_star2_binary,
                           inlist_star1_formation, inlist_star2_formation,
                           star_history_columns, binary_history_columns, profile_columns,
                           run_directory, grid_type, path_to_run_grid_exec,
                           psycris_inifile=None, keep_profiles=False,
                           keep_photos=False):
    """Based on the inifile construct the command line call to posydon-run-grid
    """
    if grid_type == "fixed":
        command_line = 'python {15} --mesa-grid {1} --mesa-binary-executable {2} '
    elif grid_type == "dynamic":
        command_line = 'mpirun --bind-to none -np {0} python -m mpi4py {15} --mesa-grid {1} --mesa-binary-executable {2} '
    else:
        raise ValueError("grid_type can either be fixed or dynamic not anything else")
    command_line += '--mesa-star1-executable {3} --mesa-star2-executable {4} --mesa-binary-inlist-project {5} '
    command_line += '--mesa-binary-inlist1 {6} --mesa-binary-inlist2 {7} --mesa-star1-inlist-project {8} '
    command_line += '--mesa-star2-inlist-project {9} --mesa-star-history-columns {10} '
    command_line += '--mesa-binary-history-columns {11} --mesa-profile-columns {12} '
    command_line += '--output-directory {13} --grid-type {14} '
    command_line += '--psycris-inifile {16}'
    if keep_profiles:
        command_line += ' --keep_profiles'
    if keep_photos:
        command_line += ' --keep_photos'
    command_line = command_line.format(number_of_mpi_processes,
                                       path_to_grid,
                                       binary_exe,
                                       star1_exe,
                                       star2_exe,
                                       inlist_binary_project,
                                       inlist_star1_binary,
                                       inlist_star2_binary,
                                       inlist_star1_formation,
                                       inlist_star2_formation,
                                       star_history_columns,
                                       binary_history_columns,
                                       profile_columns,
                                       run_directory,
                                       grid_type,
                                       path_to_run_grid_exec,
                                       psycris_inifile)
    return command_line


###############################################################################
# VALIDATION AND ENVIRONMENT HELPERS
###############################################################################

def validate_environment():
    """Validate that required environment variables are set.

    Raises
    ------
    ValueError
        If MESA_DIR is not set in the environment
    """
    if 'MESA_DIR' not in os.environ:
        raise ValueError(
            "MESA_DIR must be defined in your environment "
            "before you can run a grid of MESA runs"
        )


def find_run_grid_executable():
    """Find the posydon-run-grid executable in the system PATH.

    Returns
    -------
    str
        Path to the posydon-run-grid executable

    Raises
    ------
    ValueError
        If the executable cannot be found
    """
    proc = subprocess.Popen(
        ['which', 'posydon-run-grid'],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    (path_to_exec, err) = proc.communicate()

    if not path_to_exec:
        raise ValueError('Cannot locate posydon-run-grid executable in your path')

    return path_to_exec.decode('utf-8').strip('\n')


def validate_inifile(inifile_path):
    """Validate that the inifile exists.

    Parameters
    ----------
    inifile_path : str
        Path to the inifile

    Raises
    ------
    FileNotFoundError
        If the inifile does not exist
    """
    if not os.path.isfile(inifile_path):
        raise FileNotFoundError(
            "The provided inifile does not exist, please check the path and try again"
        )


def validate_and_setup_run_parameters(run_parameters, grid_type):
    """Validate run parameters and set defaults.

    Parameters
    ----------
    run_parameters : dict
        Dictionary of run parameters from the inifile
    grid_type : str
        Type of grid ('fixed' or 'dynamic')

    Returns
    -------
    dict
        Updated run parameters with defaults applied

    Raises
    ------
    ValueError
        If required parameters are missing or invalid
    """
    # Set defaults
    if 'keep_profiles' not in run_parameters:
        run_parameters['keep_profiles'] = False

    if 'keep_photos' not in run_parameters:
        run_parameters['keep_photos'] = False

    # Validate grid path
    grid_path = run_parameters.get('grid')
    if grid_path is None or (not os.path.isfile(grid_path) and not os.path.isdir(grid_path)):
        raise ValueError(
            "Supplied grid does not exist, please check your path and try again"
        )

    # Validate dynamic grid requirements
    if grid_type == 'dynamic' and 'psycris_inifile' not in run_parameters:
        raise ValueError(
            "Please add psycris inifile to the [run_parameters] section of the inifile."
        )

    return run_parameters


def validate_mesa_inlists(user_mesa_inlists):
    """Validate user MESA inlist configuration.

    Parameters
    ----------
    user_mesa_inlists : dict
        Dictionary of user MESA inlist settings

    Raises
    ------
    ValueError
        If required settings are missing
    """
    if ('base' not in user_mesa_inlists
        or user_mesa_inlists['base'] is None
        or user_mesa_inlists['base'] == ""):
        raise ValueError(
            "Please provide a base for the MESA run in the configuration file"
        )


###############################################################################
# CONFIGURATION BUILDING
###############################################################################

def build_configuration_stack(user_mesa_inlists, user_mesa_extras, run_parameters, run_directory):
    """Build the complete configuration stack for the MESA grid.

    This function orchestrates the layering of configurations:
    MESA defaults → POSYDON defaults → User settings → Grid parameters

    Parameters
    ----------
    user_mesa_inlists : dict
        User MESA inlist settings from the inifile
    user_mesa_extras : dict
        User MESA extras settings from the inifile
    run_parameters : dict
        Run parameters from the inifile
    run_directory : str
        Path to the run directory where files will be created

    Returns
    -------
    tuple
        (final_columns, final_extras, final_inlists, nr_systems, grid_parameters, fixgrid_file_name)
    """
    # Setup the inlist repository
    MESA_version_root_path = setup_inlist_repository(
        user_mesa_inlists.get('inlist_repository', None),
        user_mesa_inlists['mesa_version']
    )

    # Setup MESA defaults (always needed)
    MESA_default_inlists, \
    MESA_default_extras, \
    MESA_default_columns = setup_MESA_defaults(MESA_version_root_path)

    # Setup POSYDON configuration (handles MESA base internally)
    POSYDON_inlists, \
    POSYDON_extras, \
    POSYDON_columns = setup_POSYDON(
        MESA_version_root_path,
        user_mesa_inlists['base'],
        user_mesa_inlists['system_type']
    )

    # Extract user inlists, extras and columns
    user_inlists, \
    user_extras, \
    user_columns = setup_user(user_mesa_inlists, user_mesa_extras)

    # Build final columns dictionary
    final_columns = resolve_columns(
        MESA_default_columns,
        POSYDON_columns,
        user_columns
    )

    # Build final extras dictionary
    final_extras = resolve_extras(
        MESA_default_extras,
        POSYDON_extras,
        user_extras
    )

    # Read grid to get grid parameters
    nr_systems, grid_parameters, fixgrid_file_name = read_grid_file(run_parameters['grid'])

    # Extract output settings from user configuration
    user_output_settings = get_additional_user_settings(mesa_inlists=user_mesa_inlists)

    # Stack all inlist layers together
    final_inlists = resolve_inlists(
        MESA_default_inlists,
        POSYDON_inlists,
        user_inlists,
        system_type=user_mesa_inlists['system_type'],
        run_directory=run_directory,
        grid_parameters=grid_parameters,
        output_settings=user_output_settings
    )

    return final_columns, final_extras, final_inlists, nr_systems, grid_parameters, fixgrid_file_name


###############################################################################
# COMMAND LINE BUILDING
###############################################################################

def build_command_line_for_grid(slurm, output_paths, fixgrid_file_name,
                                 run_directory, run_parameters,
                                 path_to_run_grid_exec, submission_type):
    """Build the command line for running the grid.

    Parameters
    ----------
    slurm : dict
        SLURM configuration dictionary
    output_paths : dict
        Dictionary of output paths from setup_grid_run_folder
    fixgrid_file_name : str
        Path to the grid file
    run_directory : str
        Directory for output
    run_parameters : dict
        Run parameters from inifile
    path_to_run_grid_exec : str
        Path to the posydon-run-grid executable
    submission_type : str
        Type of submission ('shell' or 'slurm')

    Returns
    -------
    str
        The complete command line string
    """
    if slurm['job_array']:
        command_line = construct_command_line(
            1,
            fixgrid_file_name,
            output_paths['binary_executable'],
            output_paths['star1_executable'],
            output_paths['star2_executable'],
            output_paths['inlist_binary_project'],
            output_paths['inlist_star1_binary'],
            output_paths['inlist_star2_binary'],
            None,  # inlist_star1_formation
            None,  # inlist_star2_formation
            output_paths['star_history_columns'],
            output_paths['binary_history_columns'],
            output_paths['profile_columns'],
            run_directory,
            'fixed',
            path_to_run_grid_exec,
            keep_profiles=run_parameters['keep_profiles'],
            keep_photos=run_parameters['keep_photos']
        )
        command_line += ' --grid-point-index $SLURM_ARRAY_TASK_ID'
    else:
        # MPI mode (for future dynamic grids)
        command_line = construct_command_line(
            slurm.get('number_of_mpi_tasks', 1) * slurm.get('number_of_nodes', 1),
            fixgrid_file_name,
            output_paths['binary_executable'],
            output_paths['star1_executable'],
            output_paths['star2_executable'],
            output_paths['inlist_binary_project'],
            output_paths['inlist_star1_binary'],
            output_paths['inlist_star2_binary'],
            None,  # inlist_star1_formation
            None,  # inlist_star2_formation
            output_paths['star_history_columns'],
            output_paths['binary_history_columns'],
            output_paths['profile_columns'],
            run_directory,
            'fixed',
            path_to_run_grid_exec,
            keep_profiles=run_parameters['keep_profiles'],
            keep_photos=run_parameters['keep_photos']
        )

    # Add SLURM-specific options
    if submission_type == 'slurm':
        command_line += ' --job_end $SLURM_JOB_END_TIME'

    # Add work directory if specified
    if slurm.get('work_dir') and slurm['work_dir'] != '':
        command_line += f' --temporary-directory {slurm["work_dir"]}'

    return command_line


###############################################################################
# MAIN SETUP ENTRY POINT
###############################################################################

def run_setup(args):
    """Main entry point for the grid setup process.

    This function orchestrates the entire setup workflow:
    1. Validate environment and inputs
    2. Parse configuration file
    3. Build configuration stack (MESA → POSYDON → User)
    4. Setup run directory with all necessary files
    5. Generate submission scripts

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments containing:
        - inifile: Path to the configuration file
        - grid_type: 'fixed' or 'dynamic'
        - run_directory: Output directory
        - submission_type: 'shell' or 'slurm'
        - nproc: Number of processors
        - verbose: Verbosity setting
    """
    from posydon.utils import configfile

    # Setup logging based on verbosity level
    setup_logger(args.verbose)

    # Validate environment
    validate_environment()

    # Validate inifile exists
    validate_inifile(args.inifile)

    # Find the run grid executable
    path_to_run_grid_exec = find_run_grid_executable()

    # Parse the configuration file
    run_parameters, slurm, user_mesa_inlists, user_mesa_extras = configfile.parse_inifile(args.inifile)

    # Validate and setup run parameters
    run_parameters = validate_and_setup_run_parameters(run_parameters, args.grid_type)

    # Validate MESA inlist configuration
    validate_mesa_inlists(user_mesa_inlists)

    logger.debug(f'Base provided in inifile:\n{user_mesa_inlists["base"]}')

    # Build the complete configuration stack
    final_columns, final_extras, final_inlists, \
    nr_systems, grid_parameters, fixgrid_file_name = build_configuration_stack(
        user_mesa_inlists,
        user_mesa_extras,
        run_parameters,
        args.run_directory
    )

    # Setup the run directory with all necessary files
    output_paths = setup_grid_run_folder(
        args.run_directory,
        final_columns,
        final_extras,
        final_inlists
    )

    # Build the command line
    command_line = build_command_line_for_grid(
        slurm,
        output_paths,
        fixgrid_file_name,
        args.run_directory,
        run_parameters,
        path_to_run_grid_exec,
        args.submission_type
    )

    # Generate submission scripts
    generate_submission_scripts(args.submission_type, command_line, slurm, nr_systems)

    logger.info("Setup complete! You can now submit your grid to the cluster.")

    if args.submission_type == 'slurm':
        logger.info("To submit, run the following command:")
        logger.info("  sbatch submit_slurm.sh")
