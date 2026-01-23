import argparse
import glob
import os
import shutil
import subprocess

import pandas as pd

from posydon.active_learning.psy_cris.utils import parse_inifile
from posydon.grids.psygrid import PSyGrid
from posydon.utils import configfile
from posydon.utils import gridutils as utils
from posydon.utils.posydonwarning import Pwarn

# Define column types and their filenames
column_types = {'star_history_columns'  :'history_columns.list',
                'binary_history_columns':'binary_history_columns.list',
                'profile_columns'       :'profile_columns.list'}
column_filenames = ['history_columns.list', 'binary_history_columns.list', 'profile_columns.list']

# Define extras keys
extras_keys = ['makefile_binary', 'makefile_star', 'binary_run',
               'star_run', 'binary_extras', 'star_binary_extras', 'star1_extras', 'star2_extras',]

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
RESET = '\033[0m'
BOLD = '\033[1m'

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
    print("We are setting up your inlist repository now")

    # check if the inlist repository path exists
    if not os.path.exists(inlist_repository):
        print(f"Creating inlist repository at {inlist_repository}")
        os.makedirs(inlist_repository)

    # check if it contains anything and if not, clone the repo
    if os.listdir(inlist_repository):
        print(os.listdir(inlist_repository))
        print("Files found in inlist repository, assuming POSYDON inlist repository is already cloned!")
    else:
        out = subprocess.run(['git', 'clone', POSYDON_inlist_URL, inlist_repository],
                               capture_output=True,
                               text=True,
                               check=True,)
        print(out.stdout)

    # update the repository
    print("Updating inlist repository")
    # TODO: Re-enable git pull after testing
    print("Currently disabled for testing purposes")
    # out = subprocess.run(['git', 'pull'],
    #                cwd=inlist_repository,
    #                capture_output=True,
    #                text=True,
    #                check=True,)
    # print(out.stdout)

    # check if the base is available as a folder in the repository
    version_root_path = os.path.join(inlist_repository, MESA_version)
    if not os.path.exists(version_root_path):
        print(version_root_path)
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
    MESA_default_extras['binary_extras'] = mesa_path('binary',
                                                     'src',
                                                     'run_binary_extras.f')
    MESA_default_extras['star_binary_extras'] = mesa_path('binary',
                                                         'src',
                                                         'run_star_extras.f')

    # star1_extras and star2_extras are needed for pre-MS formation steps (if any).
    # During binary evolution, star_binary_extras is used for both stars.
    # #Both stars use the same single-star module extras file
    star_extras_path = mesa_path('star', 'src', 'run_star_extras.f')
    MESA_default_extras['star1_extras'] = star_extras_path
    MESA_default_extras['star2_extras'] = star_extras_path

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
    if base == "MESA":
        POSYDON_columns = {name: None for name in column_types}
        POSYDON_inlists = {}
        POSYDON_extras = {key: None for key in extras_keys}
        return POSYDON_inlists, POSYDON_extras, POSYDON_columns

    # Setup POSYDON base path
    POSYDON_path = os.path.join(path_to_version, base)
    check_file_exist(POSYDON_path)

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
        # We need to also include the HMS single star inlist to set up the single star evolution
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
    POSYDON_extras['binary_extras'] = os.path.join(POSYDON_path,
                                                   'extras_files',
                                                   'run_binary_extras.f')
    POSYDON_extras['star_binary_extras'] = os.path.join(POSYDON_path,
                                                       'extras_files',
                                                       'run_star_extras.f')
    POSYDON_extras['star1_extras'] = os.path.join(POSYDON_path,
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


def resolve_configuration(keys, MESA_defaults, POSYDON_config, user_config, verbose=False):
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
    verbose : bool, optional
        If True, print a visual priority table. Default is False.

    Returns
    -------
    final_config : dict
        Dictionary of final configuration paths to use
    """
    if verbose:
        print_priority_table(keys, MESA_defaults, POSYDON_config, user_config)

    final_config = {}

    for key in keys:
        if user_config.get(key) is not None:
            final_config[key] = user_config[key]
        elif POSYDON_config.get(key) is not None:
            final_config[key] = POSYDON_config[key]
        else:
            final_config[key] = MESA_defaults.get(key)

    return final_config


def resolve_columns(MESA_default_columns, POSYDON_columns, user_columns, verbose=False):
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
    verbose : bool, optional
        If True, print a visual priority table. Default is False.

    Returns
    -------
    final_columns : dict
        Dictionary of final column files paths to use
    """
    if verbose:
        print_priority_table(column_types, MESA_default_columns,
                           POSYDON_columns, user_columns,
                           title="Column Files Priority")
    return resolve_configuration(column_types, MESA_default_columns,
                                POSYDON_columns, user_columns, verbose=False)


def resolve_extras(MESA_default_extras, POSYDON_extras, user_extras, verbose=False):
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
    verbose : bool, optional
        If True, print a visual priority table. Default is False.

    Returns
    -------
    final_extras : dict
        Dictionary of final extras paths to use
    """
    if verbose:
        print_priority_table(extras_keys, MESA_default_extras,
                           POSYDON_extras, user_extras,
                           title="EXTRAS Files Priority")
    return resolve_configuration(extras_keys, MESA_default_extras,
                                POSYDON_extras, user_extras, verbose=False)


def print_priority_table(keys, MESA_defaults, POSYDON_config, user_config, title="Configuration Priority"):
    """Print a visual table showing which configuration layer is used for each key.

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
    title : str, optional
        Title for the table
    """

    # Find the longest key name for formatting
    max_key_len = max(len(str(key)) for key in keys)
    col_width = 12

    # Print title
    print(f"\n{BOLD}{title}{RESET}")
    print("=" * (max_key_len + col_width * 3 + 4))

    # Print header with colored column names
    header = f"{'Key':<{max_key_len}}  {CYAN}{'MESA':^{col_width}}{RESET}{YELLOW}{'POSYDON':^{col_width}}{RESET}{MAGENTA}{'user':^{col_width}}{RESET}"
    print(f"{BOLD}{header}{RESET}")
    print("-" * (max_key_len + col_width * 3 + 4))

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

        # Print row
        print(f"{key:<{max_key_len}}  {mesa_str:^{col_width+9}}{posydon_str:^{col_width+9}}{user_str:^{col_width+9}}")

    print("=" * (max_key_len + col_width * 3 + 4))
    print(f"{GREEN}Green{RESET} = used, {GRAY}Gray{RESET} = available but not used\n")

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
    print(f"\n{BOLD}{title}{RESET}")
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
    print(f"{MAGENTA}user{RESET} (highest priority) → {YELLOW}POSYDON{RESET} (config) → {CYAN}MESA{RESET} (base)\n")

def print_inlist_parameter_override_table(key, mesa_params, posydon_params, user_params, final_params, show_details=False):
    """Print a table showing which layer each parameter comes from, similar to extras/columns tables.

    Parameters
    ----------
    key : str
        The inlist key (e.g., 'binary_controls', 'star1_controls')
    mesa_params : dict
        Parameters from MESA defaults
    posydon_params : dict
        Parameters from POSYDON config
    user_params : dict
        Parameters from user config
    final_params : dict
        Final merged parameters
    show_details : bool, optional
        If True, show detailed parameter-by-parameter table. Default is False.
    """
    # Get all unique parameter names
    all_params = sorted(set(list(mesa_params.keys()) +
                           list(posydon_params.keys()) +
                           list(user_params.keys())))

    if not all_params:
        return

    # Only show detailed table if requested
    if not show_details:
        return

    # Find the longest parameter name for formatting
    max_param_len = max(len(str(param)) for param in all_params)
    max_param_len = max(max_param_len, 15)  # Minimum width
    col_width = 12

    # Print summary header
    overridden_count = sum(1 for p in all_params if (p in user_params and (p in mesa_params or p in posydon_params)) or
                          (p in posydon_params and p in mesa_params))
    print(f"\n  {BOLD}Detailed Parameters:{RESET} {len(all_params)} total, {overridden_count} overridden")
    print("  " + "=" * (max_param_len + col_width * 3 + 4))

    # Print header with colored column names
    header = f"  {'Parameter':<{max_param_len}}  {CYAN}{'MESA':^{col_width}}{RESET}{YELLOW}{'POSYDON':^{col_width}}{RESET}{MAGENTA}{'user':^{col_width}}{RESET}"
    print(f"{BOLD}{header}{RESET}")
    print("  " + "-" * (max_param_len + col_width * 3 + 4))

    # Print each parameter row
    for param in all_params:
        # Check which configs have this parameter
        has_mesa = param in mesa_params
        has_posydon = param in posydon_params
        has_user = param in user_params

        # Determine which one is used (priority: user > POSYDON > MESA)
        used = 'user' if has_user else ('POSYDON' if has_posydon else ('MESA' if has_mesa else None))

        # Format each column
        mesa_mark = 'x' if has_mesa else ' '
        posydon_mark = 'x' if has_posydon else ' '
        user_mark = 'x' if has_user else ' '

        # Apply colors - green for used, gray for available but not used
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

        # Print row
        print(f"  {param:<{max_param_len}}  {mesa_str:^{col_width+9}}{posydon_str:^{col_width+9}}{user_str:^{col_width+9}}")

    print("  " + "=" * (max_param_len + col_width * 3 + 4))
    print(f"  {GREEN}Green{RESET} = used, {GRAY}Gray{RESET} = available but not used")


def print_inlist_parameter_override_table_v2(key, layer_params, final_params, show_details=False):
    """Print a table showing which layer each parameter comes from (supports all layers).

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

    # Print summary header
    print(f"\n  {BOLD}Detailed Parameters:{RESET} {len(all_params)} total, {overridden_count} overridden")
    total_width = max_param_len + col_width * len(active_layers) + len(active_layers) * 2
    print("  " + "=" * total_width)

    # Print header with colored column names
    header_parts = [f"{'Parameter':<{max_param_len}}"]
    for layer_name in active_layers:
        color = layer_colors.get(layer_name, RESET)
        header_parts.append(f"{color}{layer_name:^{col_width}}{RESET}")
    header = "  ".join(header_parts)
    print(f"  {BOLD}{header}{RESET}")
    print("  " + "-" * total_width)

    # Print each parameter row
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

        print("  " + "  ".join(row_parts))

    print("  " + "=" * total_width)
    print(f"  {GREEN}Green{RESET} = used, {GRAY}Gray{RESET} = available but not used")




def print_inlist_summary_table_v2(all_keys, layer_counts):
    """Print a summary table showing parameter counts per section at each layer.

    This version supports multiple layers including grid and output configurations.

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
    total_width = max_key_len + col_width * num_layers + (num_layers + 1)

    print(f"\n{BOLD}Parameter Count Summary{RESET}")
    print("=" * total_width)

    # Print header with colored column names
    header_parts = [f"{'Section':<{max_key_len}}"]
    for layer in active_layers:
        color = layer_colors.get(layer, RESET)
        header_parts.append(f"{color}{layer:^{col_width}}{RESET}")
    header = "  ".join(header_parts)
    print(f"{BOLD}{header}{RESET}")
    print("-" * total_width)

    # Print each section row
    for key in all_keys:
        row_parts = [f"{key:<{max_key_len}}"]
        for layer in active_layers:
            count = layer_counts[layer].get(key, 0)
            color = layer_colors.get(layer, RESET)
            row_parts.append(f"{color}{count:^{col_width}}{RESET}")
        print("  ".join(row_parts))

    print("=" * total_width)
    print(f"\nLayer priority (lowest → highest): {' → '.join(active_layers)}\n")

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
        else:
            grid_layer[section] = {}

    #
    print('Adding grid parameters to sections:',
          ', '.join([sec for sec, params in grid_layer.items() if params]))

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


    print('Adding output control parameters to sections:',
              ', '.join([sec for sec, params in output_layer.items() if params]))

    # Handle ZAMS filenames if provided
    if 'zams_filename_1' in output_settings and output_settings['zams_filename_1'] is not None:
       output_layer['star1_controls']['zams_filename'] = f"'{output_settings['zams_filename_1']}'"

    if 'zams_filename_2' in output_settings and output_settings['zams_filename_2'] is not None:
       output_layer['star2_controls']['zams_filename'] = f"'{output_settings['zams_filename_2']}'"

    return output_layer


def resolve_inlists(MESA_default_inlists, POSYDON_inlists, user_inlists, system_type,
                     grid_parameters=None, output_settings=None, verbose=False, show_details=False):
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
    grid_parameters : list or set, optional
        Collection of grid parameter names. If provided, adds grid configuration layer.
    output_settings : dict, optional
        Dictionary of output settings. If provided, adds output control layer.
    verbose : bool, optional
        If True, print visual stacking and parameter count summary. Default is False.
    show_details : bool, optional
        If True, print detailed parameter-by-parameter tables. Default is False.

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
        'output': {}
    }

    # Track layer parameters for detailed printing if requested
    layer_params = {
        'MESA': {},
        'POSYDON': {},
        'user': {},
        'grid': {},
        'output': {}
    }

    # First pass: process file-based layers (MESA, POSYDON, user)
    for key in all_keys:
        # Determine the section based on the key name
        section = _get_section_from_key(key)

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


    if show_details:
        for key in all_keys:
            # Only show sections that have parameters in any layer
            if any(layer_params[layer][key] for layer in layer_params):
                print(f"\n{BOLD}═══ {key} ═══{RESET}")
                print_inlist_parameter_override_table_v2(
                    key,
                    layer_params,
                    final_inlists[key],
                    show_details=True
                )


    if verbose:
        print_inlist_summary_table_v2(all_keys, layer_counts)

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
    print(f'Wrote inlist: {filepath}')


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
    print(f'Wrote inlist: {filepath}')


def _create_build_script(path):
    """Create the 'mk' build script for compiling MESA executables.

    Parameters
    ----------
    path : str
        Path to the grid run folder
    """
    mk_filepath = os.path.join(path, 'mk')
    if os.path.exists(mk_filepath):
        print(f"Warning: 'mk' file already exists. It will be overwritten.")

    with open(mk_filepath, 'w') as f:
        f.write(f'cd {os.path.join(path, "binary/make")}\n')
        f.write(f'make -f makefile_binary\n')
        f.write(f'cd {os.path.join(path, "star1/make")}\n')
        f.write(f'make -f makefile_star\n')
        f.write(f'cd {os.path.join(path, "star2/make")}\n')
        f.write(f'make -f makefile_star\n')

    print(f'Created build script: {mk_filepath}')
    subprocess.run(['chmod', '755', mk_filepath])
    subprocess.run(['./mk'], shell=True, cwd=path)

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
    print('Using the following configuration layers for columns:')
    out_paths = {}
    for key, value in final_columns.items():
        print(f"  - {key}: {value}")
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
    print('Using the following configuration layers for extras:')

    # Define destination mapping for each extras key
    extras_destinations = {
        'makefile_binary': [('binary/make', 'makefile_binary')],
        'makefile_star': [('star1/make', 'makefile_star'),
                          ('star2/make', 'makefile_star')],
        'binary_run': [('binary/src', 'binary_run.f')],
        'star_run': [('star1/src', 'run.f'),
                     ('star2/src', 'run.f')],
        'binary_extras': [('binary/src', 'run_binary_extras.f')],
        'star_binary_extras': [('binary/src', 'run_star_extras.f')],
        'star1_extras': [('star1/src', 'run_star_extras.f')],
        'star2_extras': [('star2/src', 'run_star_extras.f')],
    }

    for key, value in final_extras.items():
        print(f"  - {key}: {value}")

        if key == 'mesa_dir':
            continue

        if key in extras_destinations:
            for subdir, filename in extras_destinations[key]:
                dest = os.path.join(path, subdir, filename)
                shutil.copy(value, dest)
        else:
            print(f"Warning: Unrecognized extras key '{key}'. Copying to root.")
            shutil.copy(value, path)


def setup_grid_run_folder(path, final_columns, final_extras, final_inlists, verbose=False):
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
    verbose : bool, optional
        If True, print additional information
    """
    # Create directory structure
    subdirs = ['binary', 'binary/make', 'binary/src',
               'star1', 'star1/make', 'star1/src',
               'star2', 'star2/make', 'star2/src',
               'column_lists']

    for subdir in subdirs:
        dir_path = os.path.join(path, subdir)
        os.makedirs(dir_path, exist_ok=True)

    print(f"\nSetting up grid run folder at: {path}")

    # Copy columns and extras
    column_paths = _copy_columns(path, final_columns)
    _copy_extras(path, final_extras)

    # Create and run build script
    _create_build_script(path)

    # Write inlist files
    print('\nWriting MESA inlist files:')
    inlist_binary_project = os.path.join(path, 'binary', 'inlist_project')
    inlist_star1_binary = os.path.join(path, 'binary', 'inlist1')
    inlist_star2_binary = os.path.join(path, 'binary', 'inlist2')

    # Add inlist names to binary_job section
    final_inlists['binary_job']['inlist_names(1)'] = f"'{inlist_star1_binary}'"
    final_inlists['binary_job']['inlist_names(2)'] = f"'{inlist_star2_binary}'"

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
            Pwarn(f'Replace {script_name}', "OverwriteWarning")

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
        print(f"Created {script_name}")

    elif submission_type == 'slurm':
        # Generate main grid submission script
        grid_script = 'job_array_grid_submit.slurm' if slurm['job_array'] else 'mpi_grid_submit.slurm'
        if os.path.exists(grid_script):
            Pwarn(f'Replace {grid_script}', "OverwriteWarning")

        array_size = nr_systems if slurm['job_array'] else None
        with open(grid_script, 'w') as f:
            _write_sbatch_header(f, slurm, job_type='grid', array_size=array_size)
            _write_environment_setup(f, slurm)
            f.write(command_line)

        # Generate cleanup script
        cleanup_script = 'cleanup.slurm'
        if os.path.exists(cleanup_script):
            Pwarn(f'Replace {cleanup_script}', "OverwriteWarning")

        with open(cleanup_script, 'w') as f:
            _write_sbatch_header(f, slurm, job_type='cleanup')
            _write_cleanup_commands(f, slurm)

        # Generate wrapper run script
        run_script = 'run_grid.sh'
        if os.path.exists(run_script):
            Pwarn(f'Replace {run_script}', "OverwriteWarning")

        with open(run_script, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write(f'ID_GRID=$(sbatch --parsable {grid_script})\n')
            f.write(f'echo "{grid_script} submitted as "${{ID_GRID}}\n')
            f.write('ID_cleanup=$(sbatch --parsable --dependency=afterany:${ID_GRID} '
                    '--kill-on-invalid-dep=yes cleanup.slurm)\n')
            f.write('echo "cleanup.slurm submitted as "${ID_cleanup}\n')

        os.system(f"chmod 755 {run_script}")
        print(f"Created {grid_script}, {cleanup_script}, and {run_script}")


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
