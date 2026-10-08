''''
This module provides the setup process for POSYDON MESA grids.

This function orchestrates the entire setup workflow:
1. Validate environment and inputs
2. Parse configuration file
3. Build configuration stack (MESA → POSYDON → User)
4. Setup run directory with all necessary files
5. Generate submission scripts
'''
import os
import shutil
import subprocess

import pandas as pd

from posydon.CLI.grids.inlist_manipulation import (
    Inlist,
    InlistManager,
    InlistSection,
    MESAInlists,
)
from posydon.CLI.log import RESET, logger, setup_logger
from posydon.grids.psygrid import PSyGrid
from posydon.utils import configfile

COLUMNS_FILES = {'star_history_columns':'history_columns.list',
                 'binary_history_columns':'binary_history_columns.list',
                 'profile_columns':'profile_columns.list'}

EXTRAS_FILES = ['makefile_binary', 'makefile_star', 'binary_run',
               'star_run', 'run_binary_extras', 'run_star_binary_extras', 'run_star1_extras', 'run_star2_extras',]

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

def validate_input(args):

    logger.debug("Validating input arguments and environment variables.")

    if 'MESA_DIR' not in os.environ:
        raise ValueError(
            "MESA_DIR must be defined in your environment "
            "before you can run a grid of MESA runs"
        )

    inifile_path = args.inifile
    if not os.path.isfile(inifile_path):
        raise FileNotFoundError(
            "The provided inifile does not exist, please check the path and try again"
        )
    logger.debug("Done")
    logger.debug('')


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

def read_configuration_file(inifile_path, grid_type):
    """Read and parse the configuration file for the grid setup.

    Parameters
    ----------
    inifile_path : str
        Path to the configuration file

    Returns
    -------
    tuple
        A tuple containing:
        - run_parameters: Dictionary of parameters for running the grid
        - slurm: Boolean indicating if SLURM is used for submission
        - user_inlist_path: Path to the user's inlist directory
        - user_inlist_extras: Additional inlist parameters from the user
    """
    logger.debug("Reading configuration file at:")
    logger.debug(f"{inifile_path}")

    config_data = configfile.parse_inifile(inifile_path)
    # unpack the configuration data into the expected variables
    run_parameters = config_data[0]
    slurm = config_data[1]
    user_inlists = config_data[2]
    user_extras = config_data[3]

    # validate the configuration data
    if 'keep_profiles' not in run_parameters:
        run_parameters['keep_profiles'] = False

    if 'keep_photo' not in run_parameters:
        run_parameters['keep_photo'] = False

    # Check if the grid exists
    grid_path = run_parameters.get('grid')
    if grid_path is None or (not os.path.isfile(grid_path) and not os.path.isdir(grid_path)):
        logger.error(grid_path)
        raise ValueError(
            "Supplied grid does not exist, please check your path and try again"
        )
    # Validate dynamic grid requirements
    if grid_type == 'dynamic' and 'psycris_inifile' not in run_parameters:
        logger.error(run_parameters)
        raise ValueError(
            "Please add psycris inifile to the [run_parameters] section of the inifile."
        )

    # user_inlists checks
    if ('base' not in user_inlists
        or user_inlists['base'] is None
        or user_inlists['base'] == ''):
            logger.error(user_inlists)
            raise ValueError(
                "Please provide a base ofr the MESA inlists in the configuration file under the [user_inlists] section."
            )
    if 'inlist_repository' not in user_inlists:
        user_inlists['inlist_repository'] = None
    if 'MESA_version' not in user_inlists:
        user_inlists['MESA_version'] = 'r11701'
    if 'repo_URL' not in user_inlists:
        user_inlists['repo_URL'] = 'https://github.com/POSYDON-code/POSYDON-MESA-INLISTS.git'



    logger.debug("Configuration file read successfully.")
    logger.debug('')

    return run_parameters, slurm, user_inlists, user_extras

def resolve_files(mesa_files, posydon_files, user_files, file_keys):
    """Generic file resolver: priority is USER > POSYDON > MESA

    Parameters
    ----------
    mesa_files : dict
    posydon_files : dict
    user_files : dict
    file_keys : list
        Keys to resolve

    Returns
    -------
    dict
        Final resolved files
    """
    final_files = {}
    for key in file_keys:
        if user_files.get(key) is not None:
            final_files[key] = user_files[key]
        elif posydon_files.get(key) is not None:
            final_files[key] = posydon_files[key]
        else:
            final_files[key] = mesa_files.get(key)
    return final_files


def setup_posydon_inlist_repository(inlist_repository,
                                    MESA_version,
                                    repo_URL,
                                    ):
    logger.info("Setting up POSYDON MESA inlist repository.")
    logger.info('This can take a few minutes depending on your internet connection.')
    if inlist_repository is None:
        # check if inlist_repository is provided
        inlist_repository = os.path.expanduser('~')

    # check if the inlist repository path exists
    if not os.path.exists(inlist_repository):
        os.makedirs(inlist_repository)

    if os.listdir(inlist_repository):
        out = subprocess.run(['git', 'remote', '-v'],
                       cwd=inlist_repository,
                       check=True,
                       capture_output=True,
                       text=True,)
        if repo_URL not in out.stdout:
            raise ValueError(f"The provided inlist repository path is not empty and does not contain the correct POSYDON inlist repository, please check the path and try again.")
        logger.debug('Existing POSYDON inlist repository found.')
        logger.debug('Files in folder:')
        logger.debug(os.listdir(inlist_repository))
        logger.debug("Files in folder. Assuming the POSYDON inlist repository is already cloned into this folder!")
    else:
        out = subprocess.run(['git', 'clone', repo_URL, inlist_repository],
                               capture_output=True,
                               text=True,
                               check=True,)
        if out.stderr:
            logger.error(out.stderr)
        else:
            logger.info("Cloned the POSYDON inlist repository successfully.")

    try:
        out = subprocess.run(['git', 'pull'],
                   cwd=inlist_repository,
                   capture_output=True,
                   text=True,
                   check=True,)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error pulling the latest changes: {e.stderr}")

    # check if the base is available as a folder in the repository
    version_root_path = os.path.join(inlist_repository, MESA_version)
    if not os.path.exists(version_root_path):
        logger.error(version_root_path)
        raise ValueError("The provided MESA version does not exist in the inlist repository, please check your provided MESA version and try again.")

    logger.info("POSYDON MESA inlist repository setup complete.")
    logger.debug("POSYDON MESA inlist repository setup complete:")
    logger.debug(version_root_path)
    logger.debug('')
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
    path_to_MESA_inlists = os.path.join(path_to_version, 'MESA_defaults', 'inlists')
    MESA_inlists = MESAInlists(path_to_MESA_inlists)

    #----------------------------------
    #              EXTRAS
    #----------------------------------

    MESA_extras = {}

    # Helper to build MESA work directory paths
    def mesa_path(module, *parts):
        return os.path.join(MESA_DIR, module, 'work', *parts)

    # Makefiles
    MESA_extras['makefile_binary'] = mesa_path('binary',
                                                       'make',
                                                       'makefile')
    MESA_extras['makefile_star'] = mesa_path('star',
                                                     'make',
                                                     'makefile')

    # Run files
    MESA_extras['star_run'] = mesa_path('star',
                                                'src',
                                                'run.f')
    MESA_extras['binary_run'] = mesa_path('binary',
                                                  'src',
                                                  'binary_run.f')

    # Extras files for binary evolution
    MESA_extras['run_binary_extras'] = mesa_path('binary',
                                                     'src',
                                                     'run_binary_extras.f')
    MESA_extras['run_star_binary_extras'] = mesa_path('binary',
                                                         'src',
                                                         'run_star_extras.f')

    # star1_extras and star2_extras are needed for pre-MS formation steps (if any).
    # During binary evolution, star_binary_extras is used for both stars.
    # #Both stars use the same single-star module extras file
    star_extras_path = mesa_path('star', 'src', 'run_star_extras.f')
    MESA_extras['run_star1_extras'] = star_extras_path
    MESA_extras['run_star2_extras'] = star_extras_path

    # Verify all extras files exist
    for _, path in MESA_extras.items():
        check_file_exist(path)

    #----------------------------------
    #             Columns
    #----------------------------------

    # Column files from MESA defaults
    MESA_columns = {
        'star_history_columns': os.path.join(MESA_DIR, 'star',
                                             'defaults', 'history_columns.list'),
        'binary_history_columns': os.path.join(MESA_DIR, 'binary',
                                               'defaults', 'binary_history_columns.list'),
        'profile_columns': os.path.join(MESA_DIR, 'star',
                                        'defaults', 'profile_columns.list')
    }

    # Verify all column files exist
    for _, path in MESA_columns.items():
        check_file_exist(path)

    logger.info("MESA defaults setup complete.")
    logger.debug('MESA inlists:')
    logger.debug(MESA_inlists)
    logger.debug('MESA extras:')
    logger.debug(MESA_extras)
    logger.debug('MESA columns:')
    logger.debug(MESA_columns)
    logger.debug('')
    return MESA_inlists, MESA_extras, MESA_columns


def setup_POSYDON(path_to_version, base, system_type, mesa_inlists, run_directory, ZAMS_filenames):
    """Setup the POSYDON configuration inlists, extras and columns based on the provided base and system type.

    Parameters
    ----------
    path_to_version : str
        Path to the POSYDON inlist repository
    base : list or str
        Base to use for the POSYDON inlists
    system_type : str
        System type to use for the POSYDON inlists
    mesa_inlists : MESAInlists
        The MESA default inlists to use as a base for the POSYDON inlists
    run_directory : str
        Directory where the run will be executed
    ZAMS_filenames : list or str
        Filename for the ZAMS model to be used
    Returns
    -------
    inlists : dict
        Dictionary of POSYDON inlists paths stacked on top of MESA defaults
    POSYDON_extras : dict
        Dictionary of POSYDON extras paths
    POSYDON_columns : dict
        Dictionary of POSYDON column files paths
        """
    if base == "MESA" or base[0] == "MESA":
        POSYDON_columns = {name: None for name in COLUMNS_FILES.keys()}
        POSYDON_inlists = {}
        POSYDON_extras = {key: None for key in EXTRAS_FILES}
        return POSYDON_inlists, POSYDON_extras, POSYDON_columns

    if len(base) == 3:
        POSYDON_inlist_folder = os.path.join(path_to_version, base[0], base[1], base[2])
    else:
        raise ValueError("Base should be a list of 3 elements corresponding to the path in the inlist repository where the POSYDON inlists are located. For example: ['POSYDON', 'DR2', 'dedt_hepulse']")
    check_file_exist(POSYDON_inlist_folder)

    logger.debug(f"Setting up POSYDON configuration: {POSYDON_inlist_folder}")

    logger.info("Setting up POSYDON configuration.")
    logger.info('Requested POSYDON configuration:')
    logger.info(f"Base: {base}")
    logger.info(f"System type: {system_type}")

    if len(ZAMS_filenames) == 2:
        ZAMS_filename_1 = ZAMS_filenames[0]
        ZAMS_filename_2 = ZAMS_filenames[1]
    elif len(ZAMS_filenames) == 1:
        ZAMS_filename_1 = ZAMS_filenames[0]
        ZAMS_filename_2 = ZAMS_filenames[0]
    else:
        raise ValueError("ZAMS_filenames not understood")

    #----------------------------------
    #            Inlists
    #----------------------------------
    inlists = InlistManager()

    # Load common base inlists once
    base_inlist1 = Inlist.from_file(f'{POSYDON_inlist_folder}/base_inlists/inlist1')
    base_inlist2 = Inlist.from_file(f'{POSYDON_inlist_folder}/base_inlists/inlist2')
    base_project_inlist = Inlist.from_file(f'{POSYDON_inlist_folder}/base_inlists/inlist_project')

    # Load HeMS setup steps (used by multiple system types)
    hems_step_0 = Inlist.from_file(f'{POSYDON_inlist_folder}/HeMS_setup_inlists/inlist_step1')
    hems_step_1 = Inlist.from_file(f'{POSYDON_inlist_folder}/HeMS_setup_inlists/inlist_step2')

    # Load single star inlist (used by multiple system types)
    single_star_inlist = Inlist.from_file(f'{POSYDON_inlist_folder}/single/single_star_inlist')

    mesa_base_project_inlist = mesa_inlists.base_binary_inlist.merge(base_project_inlist)
    mesa_base_star_inlist1 = mesa_inlists.base_star_inlist.merge(base_inlist1)
    mesa_base_star_inlist2 = mesa_inlists.base_star_inlist.merge(base_inlist2)

    # remove the zams_filename from the base inlists since we will add it in the POSYDON inlists based on the configuration, and we want to avoid confusion about which zams model is being used.
    mesa_base_star_inlist1.controls.parameters.pop('zams_filename')
    mesa_base_star_inlist2.controls.parameters.pop('zams_filename')
    hems_step_0.controls.parameters.pop('zams_filename', None) # not all steps have zams_filename, so use pop with default value to avoid errors
    hems_step_1.controls.parameters.pop('zams_filename', None)

    if system_type == 'single_HMS':
        combined_inlist = mesa_base_star_inlist1.merge(single_star_inlist)
        # Add save_model_when_terminate to the inlist for single star runs, since we want to save the final model for the single star case.
        parameters = {'save_model_when_terminate': True,
                      'save_model_filename': "'initial_star1_step0.mod'"}
        combined_inlist.add_section(InlistSection(name='star_job',
                                                    parameters=parameters))
        combined_inlist.add_section(InlistSection(name='controls',
                                                  parameters={'zams_filename': f"'{ZAMS_filename_1}'"}))
        inlists.append_star1_inlist(combined_inlist)

        # TODO: Can we remove the lines below??
        # add specific to "save" the initial model loading in the binary star1 job
        parameters = {'create_pre_main_sequence_model': False,
                      'load_saved_model': True,
                      'saved_model_name': "'initial_star1_step0.mod'"}
        tmp_inlist = Inlist(name='binary_star1_inlist')
        tmp_inlist.add_section(InlistSection(name='star_job', parameters=parameters))
        inlists.append_binary_star1_inlist(tmp_inlist)


    elif system_type == 'single_HeMS':
        # link each step
        parameters = {'save_model_when_terminate': True,
                      'save_model_filename': f"'initial_star1_step0.mod'"}
        hems_step_0.add_section(InlistSection(name='star_job', parameters=parameters))

        parameters = {'load_saved_model': True,
                      'saved_model_name': "'initial_star1_step0.mod'",
                      'save_model_when_terminate': True,
                      'save_model_filename': f"'initial_star1_step1.mod'"}
        hems_step_1.add_section(InlistSection(name='star_job', parameters=parameters))

        parameters = {'load_saved_model': True,
                      'saved_model_name': "'initial_star1_step1.mod'",
                      'save_model_when_terminate': True,
                      'save_model_filename': f"'initial_star1_step2.mod'"}
        single_star_inlist.add_section(InlistSection(name='star_job', parameters=parameters))

        inlists.append_star1_inlist(mesa_base_star_inlist1.merge(hems_step_0))
        inlists.append_star1_inlist(mesa_base_star_inlist1.merge(hems_step_1))
        inlists.append_star1_inlist(mesa_base_star_inlist1.merge(single_star_inlist))

        # TODO: Can we remove the lines below??
        # add specific to "save" the initial model loading in the binary star1 job
        parameters = {'create_pre_main_sequence_model': False,
                      'load_saved_model': True,
                      'saved_model_name': "'initial_star1_step2.mod'"}
        tmp_inlist = Inlist(name='binary_star1_inlist')
        tmp_inlist.add_section(InlistSection(name='star_job', parameters=parameters))
        inlists.append_binary_star1_inlist(tmp_inlist)


    elif system_type == 'HMS-HMS':

        inlists.append_binary_inlist(mesa_base_project_inlist)

        tmp_inlist = Inlist(name='binary_star1_inlist')
        tmp_inlist.add_section(InlistSection(name='controls', parameters={'zams_filename': f"'{ZAMS_filename_1}'"}))
        inlists.append_binary_star1_inlist(mesa_base_star_inlist1.merge(tmp_inlist))

        tmp_inlist = Inlist(name='binary_star2_inlist')
        tmp_inlist.add_section(InlistSection(name='controls', parameters={'zams_filename': f"'{ZAMS_filename_2}'"}))
        inlists.append_binary_star2_inlist(mesa_base_star_inlist2.merge(tmp_inlist))


    elif system_type == 'CO-HMS':
        co_hms_inlist1 = Inlist.from_file(f'{POSYDON_inlist_folder}/CO-HMS/inlist1')
        co_hms_inlist2 = Inlist.from_file(f'{POSYDON_inlist_folder}/CO-HMS/inlist2')
        co_hms_project = Inlist.from_file(f'{POSYDON_inlist_folder}/CO-HMS/inlist_project')

        inlists.append_binary_inlist(mesa_base_project_inlist.merge(co_hms_project))

        co_hms_inlist1.add_section(InlistSection(name='controls', parameters={'zams_filename': f"'{ZAMS_filename_1}'"}))
        # TODO: we don't really need the secondary....
        co_hms_inlist2.add_section(InlistSection(name='controls', parameters={'zams_filename': f"'{ZAMS_filename_2}'"}))

        inlists.append_binary_star1_inlist(mesa_base_star_inlist1.merge(co_hms_inlist1))
        inlists.append_binary_star2_inlist(mesa_base_star_inlist2.merge(co_hms_inlist2))

    elif system_type == 'CO-HeMS':
        co_inlist1 = Inlist.from_file(f'{POSYDON_inlist_folder}/CO-HeMS/binary/inlist1')
        co_inlist_project = Inlist.from_file(f'{POSYDON_inlist_folder}/CO-HeMS/binary/inlist_project')

        mesa_HMS_base = mesa_base_star_inlist1.merge(base_inlist1)

        inlists.append_binary_inlist(mesa_base_project_inlist.merge(co_inlist_project))

        parameters = {'save_model_when_terminate': True,
                      'save_model_filename': f"'initial_star1_step0.mod'"}

        hems_step_0.add_section(InlistSection(name='star_job', parameters=parameters))
        inlists.append_star1_inlist(mesa_HMS_base.merge(hems_step_0))

        parameters = {'load_saved_model': True,
                      'saved_model_name': "'initial_star1_step0.mod'",
                      'save_model_when_terminate': True,
                      'save_model_filename': f"'initial_star1_step1.mod'"}

        hems_step_1.add_section(InlistSection(name='star_job', parameters=parameters))
        inlists.append_star1_inlist(mesa_HMS_base.merge(hems_step_1))

        parameters = {'load_saved_model': True,
                      'saved_model_name': "'initial_star1_step1.mod'"}

        co_inlist1.add_section(InlistSection(name='star_job', parameters=parameters))
        inlists.append_binary_star1_inlist(mesa_HMS_base.merge(base_inlist1).merge(co_inlist1))

        inlists.append_binary_star2_inlist(mesa_base_star_inlist2.merge(base_inlist2))

    elif system_type == 'HeMS-HMS':

        mesa_HMS_base = mesa_base_star_inlist1.merge(base_inlist1)
        inlists.append_binary_inlist(mesa_base_project_inlist)

        parameters = {'save_model_when_terminate': True,
                      'save_model_filename': f"'initial_star1_step0.mod'"}
        hems_step_0.add_section(InlistSection(name='star_job', parameters=parameters))
        parameters = {'load_saved_model': True,
                      'saved_model_name': "'initial_star1_step0.mod'",
                      'save_model_when_terminate': True,
                      'save_model_filename': f"'initial_star1_step1.mod'"}
        hems_step_1.add_section(InlistSection(name='star_job', parameters=parameters))

        inlists.append_star1_inlist(mesa_HMS_base.merge(hems_step_0))
        inlists.append_star1_inlist(mesa_HMS_base.merge(hems_step_1))
        inlists.append_binary_star1_inlist(mesa_HMS_base)

        tmp_inlist = Inlist(name='binary_star2_inlist')
        tmp_inlist.add_section(InlistSection(name='controls', parameters={'zams_filename': f"'{ZAMS_filename_2}'"}))

        inlists.append_binary_star2_inlist(mesa_base_star_inlist2.merge(base_inlist2).merge(tmp_inlist))


    else:
        raise ValueError(f"System type {system_type} not recognized. Please check your configuration file and try again.")

    if not 'single' in system_type:
        inlist_star1_binary = os.path.join(run_directory, 'binary', 'inlist1')
        inlist_star2_binary = os.path.join(run_directory, 'binary', 'inlist2')

        inlist_names_params = {
            'inlist_names(1)': f"'{inlist_star1_binary}'",
            'inlist_names(2)': f"'{inlist_star2_binary}'"
        }
        tmp_inlist = Inlist(name='binary_inlist')
        tmp_inlist.add_section(InlistSection(name='binary_job', parameters=inlist_names_params))
        inlists.binary_inlists[0] = inlists.binary_inlists[0].merge(tmp_inlist)

    #----------------------------------
    #            Extras
    #----------------------------------

    POSYDON_extras = {}
    POSYDON_extras['run_binary_extras'] = os.path.join(POSYDON_inlist_folder,
                                                   'extras_files',
                                                   'run_binary_extras.f')
    POSYDON_extras['run_star_binary_extras'] = os.path.join(POSYDON_inlist_folder,
                                                       'extras_files',
                                                       'run_star_extras.f')
    POSYDON_extras['run_star1_extras'] = os.path.join(POSYDON_inlist_folder,
                                                 'extras_files',
                                                 'run_star_extras.f')

    #----------------------------------
    #             Columns
    #----------------------------------

    # Setup POSYDON columns
    POSYDON_columns = {}
    for name, filename in COLUMNS_FILES.items():
        file = os.path.join(POSYDON_inlist_folder, 'column_files', filename)
        # only add if the file exists
        if check_file_exist(file, raise_error=False):
            POSYDON_columns[name] = file
        else:
            POSYDON_columns[name] = None

    logger.info("POSYDON configuration setup complete.")
    logger.debug('POSYDON inlists:')
    logger.debug(inlists)
    logger.debug('POSYDON extras:')
    logger.debug(POSYDON_extras)
    logger.debug('POSYDON columns:')
    logger.debug(POSYDON_columns)
    logger.debug('')

    return inlists, POSYDON_extras, POSYDON_columns

def get_additional_user_settings(user_inlists):
    """Extract additional settings from configuration for output controls.

    Parameters
    ----------
    user_inlists : dict
        Dictionary of user inlists and configuration

    Returns
    -------
    dict
        Dictionary of additional settings
    """

    binary_inlist = Inlist(name='binary_inlists')
    binary_star1_inlist = Inlist(name='binary_star1_inlists')
    binary_star2_inlist = Inlist(name='binary_star2_inlists')


    if 'history_interval' in user_inlists:
        logger.debug('history_interaval found in user inlists, setting history_interval for binary, binary_star1 and binary_star2 inlists to user value')
        interval = user_inlists['history_interval']
        binary_inlist.add_section(InlistSection(name='binary_controls',
                                                parameters={'history_interval': interval}))
        binary_star1_inlist.add_section(InlistSection(name='controls',
                                                parameters={'history_interval': interval}))
        binary_star2_inlist.add_section(InlistSection(name='controls',
                                                parameters={'history_interval': interval}))

    if 'binary_history' in user_inlists and not user_inlists['binary_history']:
        logger.debug('User has set binary_history to False, setting history_interval to -1 to turn off binary history output')
        binary_inlist.add_section(InlistSection(name='binary_controls',
                                                parameters={'history_interval': -1}))


    for star, binary_star_inlist in zip(['star1', 'star2'],
                            (binary_star1_inlist, binary_star2_inlist)):

        final_profile_key = f'final_profile_{star}'
        final_model_key = f'final_model_{star}'
        history_key = f'history_{star}'

        if final_profile_key in user_inlists and user_inlists[final_profile_key]:
            parameters = {'write_profile_when_terminate': True,
                          'filename_for_profile_when_terminate': "'final_profile_star1.data'" if star == 'star1' else "'final_profile_star2.data'"}

            binary_star_inlist.add_section(InlistSection(name=f'star_job',
                                                         parameters=parameters))

        else:
            parameters = {'write_profile_when_terminate': False}
            binary_star_inlist.add_section(InlistSection(name=f'star_job',
                                                         parameters=parameters))

        if final_model_key in user_inlists:
            parameters = {'save_model_when_terminate': True,
                          'save_model_filename': "'final_star1.mod'" if star == 'star1' else "'final_star2.mod'"}

            binary_star_inlist.add_section(InlistSection(name=f'star_job',
                                                         parameters=parameters))

        if history_key in user_inlists:
            binary_star_inlist.add_section(InlistSection(name=f'controls',
                                                    parameters={'do_history_file': user_inlists[history_key]}))

    return (binary_inlist, binary_star1_inlist, binary_star2_inlist)


def setup_user(user_inlists, user_extras, POSYDON_inlists):
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
    logger.debug("Resolving user inlists and merging with POSYDON inlists where applicable.")

    for key in POSYDON_inlists.keys():
        num_steps = len(POSYDON_inlists[key])
        if ('star1_inlists' in key or 'star2_inlists' in key) and not 'binary' in key:
            key = key[:-1] # drop the 's' at the end of star1_inlists and star2_inlists to look for user inlist keys
            if num_steps == 0:
                continue
            elif num_steps == 1:
                # single step case: look for key without or without step
                lookup_key = f"{key}_step0" if not key.endswith("_step0") else key
                if lookup_key in user_inlists and user_inlists[lookup_key] is not None:
                    logger.debug(f"Found user inlist for {key} with key {lookup_key}")
                    # Add the key to the POSYDON_inlists
                    star1_user_inlist = Inlist.from_file(user_inlists[lookup_key])

                elif key in user_inlists and user_inlists[key] is not None:
                    logger.debug(f"Found user inlist for {key} with key {key}")
                    # Add the key to the POSYDON_inlists
                    star1_user_inlist = Inlist.from_file(user_inlists[key])
                else:
                    continue
                POSYDON_inlists[key][0] = POSYDON_inlists[key][0].merge(star1_user_inlist)

            else:
                # multiple steps case: look for key_stepX versions
                for step in range(num_steps):
                    lookup_key = f"{key}_step{step}" if not key.endswith(f"_step{step}") else key
                    if lookup_key in user_inlists and user_inlists[lookup_key] is not None:
                        logger.debug(f"Found user inlist for {key} step {step} with key {lookup_key}")
                        star1_user_inlist = Inlist.from_file(user_inlists[lookup_key])
                        POSYDON_inlists[key][step] = POSYDON_inlists[key][step].merge(star1_user_inlist)
                    else:
                        continue

        else: # binary_inlists and binary_star1/2_inlists
            if key in user_inlists and user_inlists[key] is not None:
                user_inlist = Inlist.from_file(user_inlists[key])
                POSYDON_inlists[key] = POSYDON_inlists[key].merge(user_inlist)
            else:
                continue

    logger.info("User inlists resolved and merged where applicable.")
    logger.debug('Final inlists after merging with user inlists:')
    logger.debug(POSYDON_inlists)
    logger.debug('')

    # Additional user inlist parameters that are set in user_inlists,
    # but do not correspond to full inlist files.

    logger.info("Checking for additional user inlist parameters in configuration file.")

    additional_user_settings = get_additional_user_settings(user_inlists)
    if additional_user_settings:
        binary_inlist, binary_star1_inlist, binary_star2_inlist = additional_user_settings
        if len(POSYDON_inlists['binary_inlists']) == 0:
            POSYDON_inlists['binary_inlists'].append(binary_inlist)
        else:
            POSYDON_inlists['binary_inlists'][0] = POSYDON_inlists['binary_inlists'][0].merge(binary_inlist)
        if len(POSYDON_inlists['binary_star1_inlists']) == 0:
            POSYDON_inlists['binary_star1_inlists'].append(binary_star1_inlist)
        else:
            POSYDON_inlists['binary_star1_inlists'][0] = POSYDON_inlists['binary_star1_inlists'][0].merge(binary_star1_inlist)
        if len(POSYDON_inlists['binary_star2_inlists']) == 0:
            POSYDON_inlists['binary_star2_inlists'].append(binary_star2_inlist)
        else:
            POSYDON_inlists['binary_star2_inlists'][0] = POSYDON_inlists['binary_star2_inlists'][0].merge(binary_star2_inlist)


    #----------------------------------
    #            Extras
    #----------------------------------
    # separate out extras from user inlists
    logger.info("Resolving user extras")

    user_extras = {}
    for key in EXTRAS_FILES:
        if key in user_extras.keys():
            user_extras[key] = user_extras[key]
            check_file_exist(user_extras[key])
        else:
            user_extras[key] = None

    logger.info("User extras resolved.")
    logger.debug('User extras:')
    logger.debug(user_extras)
    logger.debug('')

    #----------------------------------
    #             Columns
    #----------------------------------
    # separate out columns from user inlists
    logger.info("Resolving user columns")

    user_columns = {}

    # separate out columns from user inlists
    for name in COLUMNS_FILES.keys():
        if name in user_inlists.keys():
            user_columns[name] = user_inlists[name]
            check_file_exist(user_columns[name])
        else:
            user_columns[name] = None

    logger.info("User columns resolved.")
    logger.debug('User columns:')
    logger.debug(user_columns)
    logger.debug('')

    return POSYDON_inlists, user_extras, user_columns

def _copy_column_files(path, final_columns):
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

    logger.debug(f'COLUMN LISTS USED:')
    out_paths = {}
    for key, value in final_columns.items():
        logger.debug(f"{key}:  {value}")
        dest = os.path.join(path, 'column_lists', COLUMNS_FILES[key])
        shutil.copy(value, dest)
        out_paths[key] = dest
    return out_paths

def _copy_extras_files(path, final_extras):
    """Copy extras files (makefiles, run files) to the grid run folder.

    Parameters
    ----------
    path : str
        Path to the grid run folder
    final_extras : dict
        Dictionary mapping extras keys to file paths
    """
    logger.debug(f'EXTRAS USED:')

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


def create_run_directory(run_directory, inlists, extras, columns):
        # Create directory structure

    logger.debug(f"Creating run directory at:")
    logger.debug(f'{run_directory}')

    os.makedirs(run_directory, exist_ok=True)

    subdirs = ['binary', 'binary/make', 'binary/src',
               'star1', 'star1/make', 'star1/src',
               'star2', 'star2/make', 'star2/src',
               'column_lists']

    for subdir in subdirs:
        dir_path = os.path.join(run_directory, subdir)
        os.makedirs(dir_path, exist_ok=True)

    # columns
    _copy_column_files(run_directory, columns)
    _copy_extras_files(run_directory, extras)

    inlists.write_inlists(run_directory)

def add_grid_parameters_to_inlists(inlists, grid_parameters):

    nr, parameters, fixgrid_filename = read_grid_file(grid_parameters['grid'])
    print(parameters)
    # TODO: the old code added the extra read to the binary_star1/star2 inlists.
    # I think it should be added to star1/star2 inlists too for the single star case.
    for key in inlists.keys():
        print(key)
        if not 'binary_star1' in key and not 'binary_star2' in key:
            continue
        print('here', key)
        inlist_list = inlists[key]
        for i in range(len(inlist_list)):
            inlist = inlist_list[i]
            for section in inlist.sections.keys():
                for param in parameters:
                    if param in inlist.sections[section].parameters:
                        logger.debug(f'\t\t{param} is in {section} section of {inlist.name} inlist')
                        tmp_key = 'star1' if 'star1' in key else 'star2'
                        binary_key = 'binary_' if 'binary' in key else ''
                        out_params = {f'read_extra_{section}_inlist1': True,
                                      f'extra_{section}_inlist1_name': f"'inlist_grid_{tmp_key}_{binary_key}{section}'"}
                        # add an read_extras to the inlist with the grid parameters
                        tmp_inlist = Inlist(name=f'{key}')
                        tmp_inlist.add_section(InlistSection(name=section,
                                                              parameters=out_params))

                        inlist_list[i] = inlist_list[i].merge(tmp_inlist)

    return inlists


def run_setup(args):
    """Run the setup process for POSYDON MESA grids.


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

    setup_logger(args.verbose)

    validate_input(args)

    run_parameters, slurm, user_inlists_params, user_extras = read_configuration_file(args.inifile, args.grid_type)

    # Setup the run directory
    run_directory = args.run_directory
    if not os.path.exists(run_directory):
        os.makedirs(run_directory)

    # Setup the POSYDON MESA inlist repository
    posydon_inlist_repo_path = setup_posydon_inlist_repository(
        inlist_repository = user_inlists_params['inlist_repository'],
        MESA_version      = user_inlists_params['MESA_version'],
        repo_URL          = user_inlists_params['repo_URL'],
    )

    MESA_inlists, MESA_extras, MESA_columns = setup_MESA_defaults(posydon_inlist_repo_path)

    print(user_inlists_params['zams_filename'])
    if 'zams_filename' in user_inlists_params:
        ZAMS_filenames =  (user_inlists_params['zams_filename'],
                           user_inlists_params['zams_filename'])
    elif 'zams_filename_1' in user_inlists_params and 'zams_filename_2' in user_inlists_params:
        ZAMS_filenames = (user_inlists_params['zams_filename_1'],
                           user_inlists_params['zams_filename_2'])
    else:
        ZAMS_filenames = (None, None)


    # Stacks the MESA inlists
    (POSYDON_inlists,
     POSYDON_extras,
     POSYDON_columns) = setup_POSYDON(posydon_inlist_repo_path,
                                      user_inlists_params['base'],
                                      user_inlists_params['system_type'],
                                      MESA_inlists,
                                      run_directory,
                                      ZAMS_filenames)

    user_inlists, user_extras, user_columns = setup_user(user_inlists_params, user_extras, POSYDON_inlists)

    # add grid parameters from run_parameters to the inlist stack.

    # write/copy extras
    final_extras = resolve_files(MESA_extras, POSYDON_extras, user_extras, EXTRAS_FILES)

    # write/copy columns
    final_columns = resolve_files(MESA_columns, POSYDON_columns, user_columns, COLUMNS_FILES.keys())

    # Add grid_parameters to the POSYDON_inlists if needed, by creating a new inlist with the grid parameters and merging it on top of the existing stack.
    # TODO: remove if statement? This is to mimick the old code, which doesn't add
    # the grid parameter to the single star inlist, but adds it to binary stars.
    # I'm not sure if this is intentional or an oversight, but it seems more consistent to add the grid parameters to the single star inlist as well, since it can be used for both single and binary star runs.
    if not 'single' in user_inlists_params['system_type']:
        POSYDON_inlists = add_grid_parameters_to_inlists(POSYDON_inlists, run_parameters)

    # create the run directory with the final inlists, extras and columns
    create_run_directory(run_directory, POSYDON_inlists, final_extras, final_columns)

    # write/copy inlists
    # Build the configuration stack + write the run directory
        # 1. Get the POSYDON MESA inlist repository
        # 2. Get the MESA defaults inlists
        # 3. Check the config for POSYDON base
        # 4. Add the user inlist base to the stack
    # write/copy the submission scripts


    path_to_run_grid_exe = find_run_grid_executable()
