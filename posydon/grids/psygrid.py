"""Module defining the PSyGrid, the main grid object of POSYDON.

I. CREATING A GRID
------------------
a) One-liner for creating a grid without loading it in memory:

    mygrid = PSyGrid(verbose=True).create("mysimulations/", "mygrids/grid1.h5")

b) Create a grid and keep it for processing:

    mygrid = PSyGrid()
    mygrid.create("MESA_simulations/ZAMS_grid/", "stored_grids/ZAMS_grid.h5")
    print(mygrid)

c) Store only metadata and initial/final values of runs:

    mygrid = PSyGrid().create(in_dir_path, out_file_path, slim=True)

d) Various options (see definition of `create` and `GRIDPROPERTIES` for more):

    mygrid = PSyGrid().create(in_dir_path, out_file_path,
                              overwrite=True,
                              warn="suppress",
                              max_number_of_runs=1000,
                              description="My 1000-run grid.")

    NOTE: you can use a family of HDF5 files in case of file size restrictions.
    For example, using `mygrid.h5.%d` for the path when creating/loading grids,
    implied that split files will be `mygrid.h5.0`, `mygrid.h5.1`, and so on.
    When overwriting the file, all the family is overwritten (effectively, the
    files are emptied to avoid remainder split files from previous operations.)

e) The history/profile data can be downsampled to significantly reduce the size
   of the grid. There are three downsampling procedures: one combining all data
   from histories (binary / star 1 / star 2), and two for each star's profile.
   The parameters of the downsampling of histories and profiles can be set via
   the `GRIDPROPERTIES` dictionary (see below). For example, for history down-
   sampling, you can set a tolerance level (between 0 and 1) which corresponds
   to the maximum interpolation error tolerated with respect to the range of
   the parameters:

   history_DS_error = 0.01.

   Columns to be ignored during the downsampling process are set using the
   `history_DS_exclude` parameter. The corresponding parameters for profile
    downsampling are `profile_DS_error` and `profile_DS_exclude`.

f) The input path can be a list of paths to indicate multiple folders with MESA
   runs. Each of the paths can have wildcard characters to select subfolders or
   even runs. E.g.:
   - ["/grid_part1", "/grid_part2"]
   - "/grid_part*"
   - ["/grid/run1", "/grid/run2*"]

   Note that all provided paths are searched recursively, allowing for hierar-
   chical organizations of the MESA runs. For example, the parent folder could
   look like this:

   part_A/
     run1/
     run2/
   part_B/
     run3/
     run4/
   extra_run1/
   extra_run2/

g) Only one run can be included for a specific set of initial parameters, being
   identified by the folder name (excluding the `_grid_index*` part). E.g.,

   reruns/m1_28_m2_11_initial_z_1.42e-02_initial_period_in_days_2_grid_index_0

   will be preferred over

   main/m1_28_m2_11_initial_z_1.42e-02_initial_period_in_days_2_grid_index_102

   if the PSyGrid is created with the appropriate order of input paths:

       mygrid = PSyGrid(["mygrid/main", "mygrid/reruns"])

   The user is notified through warnings about such substitutions.

h) To include specific columns from MESA data files, set the grid properties of
   the form `*_saved_columns`. For example,

       grid.create("mygrid/", "mygrid.h5",
                   star1_history_saved_columns=["star_mass"])

   By default, the provided columns will be added to the default list. In order
   to select the exact set of columns, provide them in a tuple instead of list.


II. READING DATA FROM GRID
--------------------------
mygrid = PSyGrid("thegrid.h5")  # opens a stored grid

# Indexing...
t1 = mygrid[0]["binary_history"]["age"]
t3 = mygrid[0].binary_history["age"]

tenth_run = mygrid[10]
if tenth_run["binary_history"] is None:
    print("No binary history in 10th run.")


Use `.close()` method to close the HDF5 file, or delete the grid variable. E.g.

    mygrid.close()

    or

    del mygrid


III. OTHER FUNCTIONS
--------------------
a) Printing / logging summary:
    print(mygrid)

    with open("grid_summary.txt", "w") as f:
        f.write(str(mygrid))

b) In `for` loop (`iter` and `next` are supported too):

    for run in mygrid:
        print(run)

c) Getting number of runs:

    n_runs = len(psygrid)

d) Checking if run index is defined in grid:

    if (13 in grid):
        print("My good luck lasted for", grid[13].binary_history.age[-1], "yr")

e) Getting list or set of runs (i.e., PSyGridView objects):

    run_lst = list(psygrid)
    run_set = set(psygrid)

f) Printing parameters used during the creation of the grid

    print(psygrid.config)


WARNING: reducing reading times when accessing a PSyGrid
--------------------------------------------------------
When accessing a table in a specific run, e.g.,

    mygrid[0].history1

    or

    for run in mygrid:
        do_something_with(run.history1)

the table is loaded from the HDF5 file temporarily. This means that:

    myrun.history1.star_mass + myrun.history1.star_mass

loads the table two times! To reduce reading times, store in local variables
all the tables that are going to be needed more than once:

    myhistory1 = myrun.history1
    myhistory1.star_mass + myhistory1.star_mass

"""


__authors__ = [
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
    "Scott Coughlin <scottcoughlin2014@u.northwestern.edu>",
    "Devina Misra <devina.misra@unige.ch>",
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]


import os
import glob
import json
import ast
import warnings
import csv
import h5py
import numpy as np
import pandas as pd
import tqdm

from posydon.utils.common_functions import check_state_of_star_history_array
from posydon.grids.io import (GridReader, read_initial_values,
                              initial_values_from_dirname, EEP_FILE_EXTENSIONS)
from posydon.grids.termination_flags import (get_flags_from_MESA_run,
                                             check_state_from_history,
                                             get_flag_from_MESA_output,
                                             infer_interpolation_class,
                                             get_detected_initial_RLO,
                                             get_nearest_known_initial_RLO)
from posydon.utils.configfile import ConfigFile
from posydon.utils.common_functions import (orbital_separation_from_period,
                                            initialize_empty_array,
                                            infer_star_state)
from posydon.utils.gridutils import (read_MESA_data_file, read_EEP_data_file,
                                     add_field, join_lists, fix_He_core)
from posydon.visualization.plot2D import plot2D
from posydon.visualization.plot1D import plot1D
from posydon.grids.downsampling import TrackDownsampler
from posydon.grids.scrubbing import scrub, keep_after_RLO


HDF5_MEMBER_SIZE = 2**31 - 1            # maximum HDF5 file size when splitting
H5_UNICODE_DTYPE = h5py.string_dtype()  # dtype for unicode strings in the HDF5
H5_REC_STR_DTYPE = "S70"    # string dtype when writing record arrays to HDF5

# valid keys in a run in the HDF5 file
VALID_KEYS = ["binary_history", "history1", "history2", "final_profile1",
              "final_profile2", "initial_values", "final_values"]

# Acceptable values for `warn` argument of PSyGrid.create()
WARN_VALUES = ["end", "normal", "suppress"]

# Number of termination flags per run
N_FLAGS = 4
N_FLAGS_SINGLE = 2
TERMINATION_FLAG_COLUMNS = ["termination_flag_{}".format(i+1)
                            for i in range(N_FLAGS)]
TERMINATION_FLAG_COLUMNS_SINGLE = ["termination_flag_1", "termination_flag_3"]

# Default columns to be included from history and profile tables
DEFAULT_BINARY_HISTORY_COLS = [
    "model_number", "age",
    "star_1_mass", "star_2_mass",
    "period_days", "binary_separation",
    "lg_system_mdot_1", "lg_system_mdot_2",
    "lg_wind_mdot_1", "lg_wind_mdot_2",
    "lg_mstar_dot_1", "lg_mstar_dot_2",
    "lg_mtransfer_rate", "xfer_fraction",
    "rl_relative_overflow_1", "rl_relative_overflow_2",
    "trap_radius", "acc_radius",
    "t_sync_rad_1", "t_sync_conv_1", "t_sync_rad_2", "t_sync_conv_2"
]

DEFAULT_STAR_HISTORY_COLS = [
    "he_core_mass", "c_core_mass", "o_core_mass",
    "he_core_radius", "c_core_radius", "o_core_radius",
    "center_h1", "center_he4", "center_c12", "center_n14", "center_o16",
    "surface_h1", "surface_he4", "surface_c12", "surface_n14", "surface_o16",
    "c12_c12", "center_gamma", "avg_c_in_c_core",
    "surf_avg_omega", "surf_avg_omega_div_omega_crit",
    "log_LH", "log_LHe", "log_LZ", "log_Lnuc", "log_Teff", "log_L", "log_R",
    "log_center_T", "log_center_Rho",
    "total_moment_of_inertia", "spin_parameter", "log_total_angular_momentum",
    "conv_env_top_mass", "conv_env_bot_mass",
    "conv_env_top_radius", "conv_env_bot_radius",
    "conv_env_turnover_time_g", "conv_env_turnover_time_l_b",
    "conv_env_turnover_time_l_t", "envelope_binding_energy",
    "mass_conv_reg_fortides", "thickness_conv_reg_fortides",
    "radius_conv_reg_fortides",
    "lambda_CE_1cent", "lambda_CE_10cent", "lambda_CE_30cent",
    "co_core_mass", "co_core_radius", "lambda_CE_pure_He_star_10cent",
    "log_L_div_Ledd"
]

DEFAULT_SINGLE_HISTORY_COLS = (["model_number", "star_age", "star_mass"]
                               + DEFAULT_STAR_HISTORY_COLS)

DEFAULT_EEP_HISTORY_COLS = ["star_age", "star_mass"]+DEFAULT_STAR_HISTORY_COLS

DEFAULT_PROFILE_COLS = [
    "radius",
    "mass",
    "logRho",
    "omega",
    "energy",
    "x_mass_fraction_H",
    "y_mass_fraction_He",
    "z_mass_fraction_metals",
    "neutral_fraction_H",
    "neutral_fraction_He",
    "avg_charge_He",
]

# Default columns to exclude from computing history/profile downsampling
DEFAULT_HISTORY_DS_EXCLUDE = [
    "model_number",     # not physical
    "age",              # this is the independent variable in binary grids
    "star_age"          # this is the independent variable in single grids
]

DEFAULT_PROFILE_DS_EXCLUDE = [
    "mass",             # this is the independent variable in binary grids
    "star_mass",        # this is the independent variable in single grids
]

# to gain extra compression exclude the follwoing columns
# we support two versions:
# - SMALL: DS_err = 0.01 and exclude COLS+['surface_he4', 'surface_h1']
# - LITE: DS_err = 0.1 and exclude COLS
EXTRA_COLS_DS_EXCLUDE = [
    'acc_radius', 'age', 'avg_charge_He', 'c12_c12', 'c_core_radius',
    'center_c12', 'center_gamma', 'center_n14', 'center_o16', 'co_core_radius',
    'conv_env_bot_mass', 'conv_env_bot_radius', 'conv_env_top_mass',
    'conv_env_top_radius', 'conv_env_turnover_time_g',
    'conv_env_turnover_time_l_b', 'conv_env_turnover_time_l_t', 'energy',
    'envelope_binding_energy', 'he_core_radius', 'lambda_CE_10cent',
    'lambda_CE_1cent', 'lambda_CE_30cent', 'lambda_CE_pure_He_star_10cent',
    'lg_mstar_dot_1', 'lg_mstar_dot_2', 'lg_wind_mdot_1', 'lg_wind_mdot_2',
    'log_LH', 'log_LHe', 'log_LZ', 'log_Lnuc', 'mass',
    'mass_conv_reg_fortides', 'model_number', 'neutral_fraction_H',
    'neutral_fraction_He', 'o_core_radius', 'radius_conv_reg_fortides',
    'rl_relative_overflow_1', 'rl_relative_overflow_2', 'spin_parameter',
    'star_age', 'star_mass', 'surf_avg_omega', 'surf_avg_omega_div_omega_crit',
    'surface_c12', 'surface_h1', 'surface_he4', 'surface_n14', 'surface_o16',
    't_sync_conv_1', 't_sync_conv_2', 't_sync_rad_1', 't_sync_rad_2',
    'thickness_conv_reg_fortides', 'trap_radius', 'x_mass_fraction_H',
    'xfer_fraction', 'y_mass_fraction_He', 'z_mass_fraction_metals',
    'log_L_div_Ledd'
]


# Properties defining the PSyGrid creation. They're saved as metadata.
GRIDPROPERTIES = {
    # file loading parameters
    "description": "",                      # description text
    "max_number_of_runs": None,
    "format": "hdf5",
    "compression": "gzip9",
    # history downsampling parameters
    "history_DS_error": None,
    "history_DS_exclude": DEFAULT_HISTORY_DS_EXCLUDE,
    # profile downsampling parameters
    "profile_DS_error": None,
    "profile_DS_interval": None,
    "profile_DS_exclude": DEFAULT_PROFILE_DS_EXCLUDE,
    # columns to pass to the PSyGrid
    "star1_history_saved_columns": "minimum",
    "star2_history_saved_columns": "minimum",
    "binary_history_saved_columns": "minimum",
    "star1_profile_saved_columns": "minimum",
    "star2_profile_saved_columns": "minimum",
    # Initial/Final value arrays
    "initial_value_columns": None,
    "final_value_columns": None,
    # grid-specific arguments
    "start_at_RLO": False,
    "binary": True,
    "eep": None,    # path to EEP files
    "initial_RLO_fix": False,
    "He_core_fix": True,
    "accept_missing_profile": True,
}


class PSyGrid:
    """Class handling a grid of MESA runs encoded in HDF5 format."""

    def __init__(self, filepath=None, verbose=False):
        """Initialize the PSyGrid, and if `filename` is provided, open it.

        Parameters
        ----------
        verbose : bool
            If `True`, the objects reports by printing to standard output.
        path : str or None
            If not None, it is the path of the HDF5 file to be loaded.

        """
        self._reset()    # initialize all except for `filepath` and `verbose`
        self.filepath = filepath
        self.verbose = verbose
        if self.filepath is not None:
            self.load(self.filepath)

    def _reset(self):
        """(Re)set attributes to defaults, except `filename` and `verbose`."""
        self.close()
        self.hdf5 = None                    # h5py File handle if file is open
        self.MESA_dirs = []
        self.initial_values = None
        self.final_values = None
        self.config = ConfigFile()
        self._make_compression_args()
        self.n_runs = 0
        self.eeps = None

    def _make_compression_args(self):
        if "compression" not in self.config:
            self.compression_args = {}
            return

        compression = self.config["compression"]
        compression = "" if compression is None else compression.lower()
        if compression == "":
            compression_args = {}
        elif compression == "lzf":
            compression_args = {"compression": "lzf"}
        elif compression.startswith("gzip"):
            if compression == "gzip":
                compression_args = {"compression": "gzip"}
            else:
                if len(compression) > 5 or not compression[4].isdigit():
                    raise ValueError("GZIP compression level must be 0-9.")
                compression_level = int(compression[4])
                compression_args = {"compression": "gzip",
                                    "compression_opts": compression_level}
        else:
            raise ValueError("Unknown compression `{}`.".format(compression))
        self.compression_args = compression_args

    def _discover_eeps(self, path):
        self.eeps = {}
        for extension in EEP_FILE_EXTENSIONS:
            searchfor = os.path.join(path, "*" + extension)
            for filename in glob.glob(searchfor):
                identifier = os.path.basename(filename.split(extension)[-2])
                assert identifier not in self.eeps
                self.eeps[identifier] = filename
        if len(self.eeps) == 0:
            self.eeps = None

    def _say(self, something):
        if self.verbose:
            print(something, flush=True)

    def generate_config(self, **grid_kwargs):
        """Store the grid's configuration in a ConfigFile object."""
        self.config = ConfigFile()

        for key in grid_kwargs:
            if key not in GRIDPROPERTIES:
                raise ValueError("`{}` is not a valid parameter name.".
                                 format(key))

        new_dict = {}
        for varname in GRIDPROPERTIES:
            default_value = GRIDPROPERTIES[varname]
            new_dict[varname] = grid_kwargs.get(varname, default_value)

        self.config.update(new_dict)
        self._make_compression_args()

    def create(self, MESA_grid_path, psygrid_path=None, overwrite=False,
               slim=False, warn="end", fmt="posydon", **grid_kwargs):
        """Create a new PSyGrid object from a MESA grid and save it on disk.

        Parameters
        ----------
        MESA_grid_path : str
            The path of the top directory where the input grid is located.
        psygrid_path : str
            The path of the HDF5 file to be created (or overwritten). If not
            provided, it is assumed that it was defined during initialization.
        overwrite : bool
            Whether to overwrite the HDF5 file if it already exists.
        slim : bool
            If `True`, only the metadata and initial/final values are stored.
        warn : str
            How warnings are handled. If "normal", then warnings are shown when
            they occur. If "end", they are collected and shown at the end. If
            "suppress", then warnings are not shown at all.
        **grid_kwargs : dict
            Configuration of the new grid, passed as ddditional arguments.

        """
        self._reset()

        if warn not in WARN_VALUES:
            raise ValueError("`warn` must be in: {}".format(WARN_VALUES))

        if psygrid_path is not None:
            self.filepath = psygrid_path

        if self.filepath is None:
            ValueError("The path of the HDF5 file was not defined.")

        # Create the output directory if it doesn't exist
        dirpath = os.path.dirname(self.filepath)
        if dirpath != "":
            os.makedirs(dirpath, exist_ok=True)

        # Create the configuration object from the keyword arguments
        self.generate_config(**grid_kwargs)

        # Create the PSyGrid HDF5 file. Handle warnings according to `warn`.

        # if already exists, am I allowed to overwrite?
        if not overwrite and os.path.exists(self.filepath):
            raise Exception("File {} already exists.".format(self.filepath))

        # open the HDF5 file and transcribe the grid
        driver_args = {} if "%d" not in self.filepath else {
            "driver": "family", "memb_size": HDF5_MEMBER_SIZE}

        with h5py.File(self.filepath, "w", **driver_args) as hdf5:
            if warn == "normal":
                self._create_psygrid(MESA_grid_path,
                                     hdf5=hdf5, slim=slim, fmt=fmt)
            else:
                collected_warnings = []
                with warnings.catch_warnings(record=True) as caught_warnings:
                    self._create_psygrid(MESA_grid_path,
                                         hdf5=hdf5, slim=slim, fmt=fmt)
                    collected_warnings = caught_warnings

                if warn == "end":
                    for warning_message in collected_warnings:
                        warnings.showwarning(warning_message.message,
                                             warning_message.category,
                                             warning_message.filename,
                                             warning_message.lineno,
                                             line="")
                else:
                    # for consistency
                    assert warn == "suppress"

        self.load()

    def _create_psygrid(self, MESA_path, hdf5, slim=False, fmt="posydon"):
        """Generate a PSyGrid instance from a path to MESA data.

        Parameters
        ----------
        MESA_path : str
            The path of the directory with the MESA runs.
        overwrite : bool
            Whether existing HDF5 file can be overwritten.
        slim : bool
            Whether to include runs' history/profile data in the HDF5 file.

        """
        outpath = self.filepath
        assert outpath is not None  # for consistency - expected to be defined

        binary_grid = self.config["binary"]
        initial_RLO_fix = self.config["initial_RLO_fix"]
        start_at_RLO = self.config["start_at_RLO"]
        eep = self.config["eep"]

        if eep is not None:
            if binary_grid:
                warnings.warn("Selected EEPs, switching to single-star grid.")
                self.config["binary"] = True
                binary_grid = True
            self._discover_eeps(eep)

        # prepare for loss-less compression
        compression_args = self.compression_args

        self._say("Reading input directory structure...")
        grid = GridReader(MESA_path, fmt=fmt, binary=binary_grid,
                          verbose=self.verbose)

        # Determine expected number of runs (may not be as many in the end)
        N_runs = len(grid.runs)
        if not self.config['max_number_of_runs'] is None:
            N_runs = min(N_runs, self.config['max_number_of_runs'])

        # Decide on the columns to be included
        def decide_columns(key_in_config, defaults):
            """Get columns to be saved or `None` if setting is 'all')."""
            inconfig = self.config[key_in_config]
            if inconfig == "all":
                return None
            if inconfig == "minimum":
                return defaults
            if np.iterable(inconfig):
                if isinstance(inconfig, tuple):
                    return list(inconfig)
                return join_lists(defaults, inconfig)
            raise Exception("{} setting not recognized.".format(key_in_config))

        if binary_grid:
            H1_DEFAULTS = DEFAULT_STAR_HISTORY_COLS
        else:
            if self.eeps is None:
                H1_DEFAULTS = DEFAULT_SINGLE_HISTORY_COLS
            else:
                H1_DEFAULTS = DEFAULT_EEP_HISTORY_COLS

        BH_columns = decide_columns("binary_history_saved_columns",
                                    DEFAULT_BINARY_HISTORY_COLS)
        H1_columns = decide_columns("star1_history_saved_columns",
                                    H1_DEFAULTS)
        H2_columns = decide_columns("star2_history_saved_columns",
                                    DEFAULT_STAR_HISTORY_COLS)
        P1_columns = decide_columns("star1_profile_saved_columns",
                                    DEFAULT_PROFILE_COLS)
        P2_columns = decide_columns("star2_profile_saved_columns",
                                    DEFAULT_PROFILE_COLS)

        self._say("Infering column names in history tables...")
        all_history_columns = grid.infer_history_columns(
            BH_columns, H1_columns, H2_columns)
        if len(all_history_columns) == 0:
            raise Exception("PSyGrid object cannot be created: "
                            "No history data in all runs.")

        self._say("Preparing initial/final_values arrays.")
        initial_value_columns = all_history_columns.copy()
        termination_flag_columns = (TERMINATION_FLAG_COLUMNS
                                    if binary_grid else
                                    TERMINATION_FLAG_COLUMNS_SINGLE)
        dtype_initial_values = (
            [(col, 'f8') for col in initial_value_columns]
            + [(col, 'f8') for col in ["X", "Y", "Z"]]
        )
        dtype_final_values = (
            [(col, 'f8') for col in all_history_columns]
            + [(col, H5_UNICODE_DTYPE) for col in termination_flag_columns]
            + ([("interpolation_class", H5_UNICODE_DTYPE)]
               if binary_grid else [])
        )
        self.initial_values = initialize_empty_array(
            np.empty(N_runs, dtype=dtype_initial_values))
        self.final_values = initialize_empty_array(
            np.empty(N_runs, dtype=dtype_final_values))

        # in case lists were used, get a proper `dtype` object
        dtype_initial_values = self.initial_values.dtype
        dtype_final_values = self.final_values.dtype

        self._say('Loading MESA data...')

        #this int array will store the run_index of the i-th run and will be
        # -1 if the run is not included.
        run_included_at = np.full(N_runs, -1, dtype=int)
        run_index = 0
        for i in tqdm.tqdm(range(N_runs)):
            # Select the ith run
            run = grid.runs[i]
            ignore_data = False    # if failed run, do not save any data
            self._say('Processing {}'.format(run.path))

            # Restrict number of runs if limit is set inside config
            if self.config['max_number_of_runs'] is not None:
                if run_index == self.config['max_number_of_runs']:
                    self._say("Maximum number of runs reached.")
                    break

            binary_history = read_MESA_data_file(run.binary_history_path,
                                                 BH_columns)
            # if no binary history, ignore this run
            if binary_grid and binary_history is None:
                ignore_data = True
                ignore_reason = "ignored_no_BH"
                warnings.warn("Ignored MESA run because of missing binary "
                              "history in: {}\n".format(run.path))
                if not initial_RLO_fix:
                    continue

            if ignore_data:
                history1 = None
            elif self.eeps is None:
                history1 = read_MESA_data_file(run.history1_path, H1_columns)
            else:
                try:
                    eep_path = self.eeps[os.path.basename(run.path)]
                except KeyError:
                    warnings.warn("No matching EEP file for `" + run.path
                                  + "`. Ignoring run.")
                    continue
                history1 = read_EEP_data_file(eep_path, H1_columns)
            if self.config["He_core_fix"]:
                history1 = fix_He_core(history1)

            if not binary_grid and history1 is None:
                warnings.warn("Ignored MESA run because of missing "
                              "history in: {}\n".format(run.path))
                ignore_data = True
                ignore_reason = "ignore_no_H1"
                continue

            if ignore_data:
                history2 = None
            else:
                history2 = read_MESA_data_file(run.history2_path, H2_columns)
            if self.config["He_core_fix"]:
                history2 = fix_He_core(history2)

            # scrub histories (unless EEPs are selected or run is ignored)
            if not ignore_data and self.eeps is None:
                # read the model numbers and ages from the histories
                colname = "model_number"
                if history1 is not None:
                    history1_mod = np.int_(read_MESA_data_file(
                        run.history1_path, [colname])[colname])
                    if len(history1_mod) == len(history1) + 1:
                        history1_mod = history1_mod[:-1]
                    history1_age = read_MESA_data_file(
                        run.history1_path, ["star_age"])["star_age"]
                    if len(history1_age) == len(history1) + 1:
                        history1_age = history1_age[:-1]
                else:
                    history1_mod = None
                    history1_age = None

                if history2 is not None:
                    history2_mod = np.int_(read_MESA_data_file(
                        run.history2_path, [colname])[colname])
                    if len(history2_mod) == len(history2) + 1:
                        history2_mod = history2_mod[:-1]
                    history2_age = read_MESA_data_file(
                        run.history2_path, ["star_age"])["star_age"]
                    if len(history2_age) == len(history2) + 1:
                        history2_age = history2_age[:-1]
                else:
                    history2_mod = None
                    history2_age = None

                if binary_history is not None:
                    binary_history_mod = np.int_(read_MESA_data_file(
                        run.binary_history_path, [colname])[colname])
                    if len(binary_history_mod) == len(binary_history) + 1:
                        binary_history_mod = binary_history_mod[:-1]
                    binary_history_age = read_MESA_data_file(
                        run.binary_history_path, ["age"])["age"]
                    if len(binary_history_age) == len(binary_history) + 1:
                        binary_history_age = binary_history_age[:-1]
                else:
                    binary_history_mod = None
                    binary_history_age = None

                # scrub the histories
                if binary_grid:
                    history1, history2, binary_history = scrub(
                        tables=[history1, history2, binary_history],
                        models=[history1_mod, history2_mod,
                                binary_history_mod],
                        ages=[history1_age, history2_age, binary_history_age]
                    )
                else:
                    history1, history2, binary_history = scrub(
                        tables=[history1, None, None],
                        models=[history1_mod, None, None],
                        ages=[history1_age, None, None]
                    )

                # if scubbing wiped all the binary history, discard run
                if binary_grid and len(binary_history) == 0:
                    ignore_data = True
                    ignore_reason = "ignored_scrubbed"
                    warnings.warn("Ignored MESA run because of scrubbed binary"
                                  " history in: {}\n".format(run.path))
                    if not initial_RLO_fix:
                        continue
                if not binary_grid and len(history1) == 0:
                    ignore_data = True
                    warnings.warn("Ignored MESA run because of scrubbed"
                                  " history in: {}\n".format(run.path))
                    continue

                # check whether start at RLO is requested, and chop the history
                if start_at_RLO:
                    kept = keep_after_RLO(binary_history, history1, history2)
                    if kept is None:
                        ignore_data = True
                        ignore_reason = "ignored_no_RLO"
                        warnings.warn("Ignored MESA run because of no RLO"
                                      " in: {}\n".format(run.path))
                        if not initial_RLO_fix:
                            continue
                        binary_history, history1, history2 = None, None, None
                    else:
                        binary_history, history1, history2 = kept

            if ignore_data:
                binary_history, history1, history2 = None, None, None
                final_profile1, final_profile2 = None, None
            else:
                final_profile1 = read_MESA_data_file(
                    run.final_profile1_path, P1_columns)
                final_profile2 = read_MESA_data_file(
                    run.final_profile2_path, P2_columns)
                if not binary_grid and final_profile1 is None:
                    if self.config["accept_missing_profile"]:
                        warnings.warn("Including MESA run despite the missing "
                                      "profile in {}\n".format(run.path))
                        ignore_data = False
                    else:
                        warnings.warn("Ignored MESA run because of missing "
                                      "profile in: {}\n".format(run.path))
                        ignore_data = True
                        ignore_reason = "ignore_no_FP"
                        continue

            if binary_history is not None:
                if binary_history.shape == ():  # if history is only one line
                    initial_BH = binary_history.copy()
                    final_BH = binary_history.copy()
                else:
                    # disentagle arrays with `copy` because the LHSs may change
                    initial_BH = binary_history[0].copy()
                    final_BH = binary_history[-1].copy()
            else:
                initial_BH = None
                final_BH = None

            if history1 is not None:
                if history1.shape == ():  # if history is only one line
                    initial_H1 = history1.copy()
                    final_H1 = history1.copy()
                else:
                    # disentagle arrays with `copy` because the LHSs may change
                    initial_H1 = history1[0].copy()
                    final_H1 = history1[-1].copy()
            else:
                initial_H1 = None
                final_H1 = None

            if history2 is not None:
                if history2.shape == ():  # if history is only one line
                    initial_H2 = history2.copy()
                    final_H2 = history2.copy()
                else:
                    # disentagle arrays with `copy` because the LHSs may change
                    initial_H2 = history2[0].copy()
                    final_H2 = history2[-1].copy()
            else:
                initial_H2 = None
                final_H2 = None

            # get some initial values from the `binary_history.data` header
            # if of course, no RLO fix is applied
            if binary_grid and not (start_at_RLO or ignore_data):
                bh_header = np.genfromtxt(run.binary_history_path,
                                          skip_header=1,
                                          max_rows=1, names=True)

                can_compute_new_separation = True

                if "star_1_mass" in dtype_initial_values.names:
                    init_mass_1 = bh_header["initial_don_mass"]
                    initial_BH["star_1_mass"] = init_mass_1
                else:
                    can_compute_new_separation = False

                if "star_2_mass" in dtype_initial_values.names:
                    init_mass_2 = bh_header["initial_acc_mass"]
                    initial_BH["star_2_mass"] = init_mass_2
                else:
                    can_compute_new_separation = False

                if "period_days" in dtype_initial_values.names:
                    init_period = bh_header["initial_period_days"]
                    initial_BH["period_days"] = init_period
                else:
                    can_compute_new_separation = False

                if ("binary_separation" in dtype_initial_values.names
                        and can_compute_new_separation):
                    init_separation = orbital_separation_from_period(
                        init_period, init_mass_1, init_mass_2)
                    initial_BH["binary_separation"] = init_separation
            elif not binary_grid and not (start_at_RLO or ignore_data):
                # use header to get initial mass in single-star grids
                h1_header = np.genfromtxt(run.history1_path,
                                          skip_header=1,
                                          max_rows=1, names=True)
                if "S1_star_mass" in dtype_initial_values.names:
                    init_mass_1 = h1_header["initial_m"]
                    initial_H1["star_mass"] = init_mass_1

            # get some initial values from the `LOGS1/history.data` header
            addX = "X" in dtype_initial_values.names
            addY = "Y" in dtype_initial_values.names
            addZ = "Z" in dtype_initial_values.names
            if (addX or addY or addZ) and not ignore_data:
                # read abundances from history1 if present, else history2
                if history1 is not None:
                    read_from = run.history1_path
                elif history2 is not None:
                    read_from = run.history2_path
                else:
                    read_from = None

                # if not star history file, NaN values
                if read_from is None:
                    init_X, init_Y, init_Z = np.nan, np.nan, np.nan
                else:
                    star_header = np.genfromtxt(
                        read_from, skip_header=1, max_rows=1, names=True)
                    init_Z = star_header["initial_Z"]
                    init_Y = star_header["initial_Y"]
                    init_X = 1.0 - init_Z - init_Y

                if binary_grid:
                    where_to_add = initial_BH
                else:
                    where_to_add = initial_H1

                if addZ:
                    where_to_add = add_field(where_to_add, [("Z", '<f8')])
                    where_to_add["Z"] = init_Z
                if addY:
                    where_to_add = add_field(where_to_add, [("Y", '<f8')])
                    where_to_add["Y"] = init_Y
                if addX:
                    where_to_add = add_field(where_to_add, [("X", '<f8')])
                    where_to_add["X"] = init_X

            # set initial/final values
            if initial_BH is not None:
                for col in initial_BH.dtype.names:
                    self.initial_values[i][col] = initial_BH[col]
                for col in final_BH.dtype.names:
                    self.final_values[i][col] = final_BH[col]
            if initial_H1 is not None:
                for col in initial_H1.dtype.names:
                    newcol = "S1_" + col
                    self.initial_values[i][newcol] = initial_H1[col]
                for col in final_H1.dtype.names:
                    newcol = "S1_" + col
                    self.final_values[i][newcol] = final_H1[col]
            if initial_H2 is not None:
                for col in initial_H2.dtype.names:
                    newcol = "S2_" + col
                    self.initial_values[i][newcol] = initial_H2[col]
                for col in final_H2.dtype.names:
                    newcol = "S2_" + col
                    self.final_values[i][newcol] = final_H2[col]

            if not ignore_data:
                if addX:
                    self.initial_values["X"] = where_to_add["X"]
                if addY:
                    self.initial_values["Y"] = where_to_add["Y"]
                if addZ:
                    self.initial_values["Z"] = where_to_add["Z"]

            if binary_grid:
                if ignore_data:
                    termination_flags = [ignore_reason] * N_FLAGS
                else:
                    termination_flags = get_flags_from_MESA_run(
                        run.out_txt_path, binary_history=binary_history,
                        history1=history1, history2=history2,
                        start_at_RLO=start_at_RLO)
            else:
                if ignore_data:
                    termination_flags = [ignore_reason] * N_FLAGS_SINGLE
                else:
                    termination_flags = [
                        get_flag_from_MESA_output(run.out_txt_path),
                        check_state_from_history(history1,
                                                 history1["star_mass"])
                    ]

            if ignore_data:
                for colname in self.final_values.dtype.names:
                    if (colname.startswith("termination_flag_")
                            or colname.startswith("interpolation_class")):
                        self.final_values[i][colname] = ignore_reason
                    else:
                        self.final_values[i][colname] = np.nan
                for colname in self.initial_values.dtype.names:
                    self.initial_values[i][colname] = np.nan

                # now update the initial values from the grid point data
                grid_point = read_initial_values(run.path)
                for colname, value in grid_point.items():
                    if colname in self.initial_values.dtype.names:
                        self.initial_values[i][colname] = value
            else:
                for flag, col in zip(termination_flags,
                                     termination_flag_columns):
                    if flag is None:
                        flag = ""
                    self.final_values[i][col] = flag

                if binary_grid:
                    self.final_values[i]["interpolation_class"] = \
                        infer_interpolation_class(*termination_flags[:2])

            if ignore_data:
                # if not fix requested and failed run, do not include it
                continue

            # if data to be included, downsample and store
            if not slim:
                if ignore_data:
                    hdf5.create_group("/grid/run{}/".format(run_index))
                else:
                    # downsample and save
                    histories_DS = downsample_history(
                        binary_history, history1, history2, params=self.config)
                    binary_history_DS, history1_DS, history2_DS = histories_DS

                    final_profile1_DS = downsample_profile(
                        final_profile1, params=self.config)
                    final_profile2_DS = downsample_profile(
                        final_profile2, params=self.config)

                    # store the data
                    arrays_to_save = [binary_history_DS,
                                      history1_DS,
                                      history2_DS,
                                      final_profile1_DS,
                                      final_profile2_DS]
                    keys_to_arrays = ['binary_history',
                                      'history1', 'history2',
                                      'final_profile1',
                                      'final_profile2']

                    for array, name in zip(arrays_to_save, keys_to_arrays):
                        if array is not None:
                            hdf5.create_dataset(
                                "/grid/run{}/{}".format(run_index, name),
                                data=array, **compression_args)

            # consider the run (and the input directory) included
            self.MESA_dirs.append(run.path)
            run_included_at[i] = run_index
            run_index += 1
            #check that new MESA path is added at run_index
            lenMESA_dirs = len(self.MESA_dirs)
            if lenMESA_dirs!=run_index:
                warnings.warn("Non synchronous indexing: " +
                          "run_index={} != ".format(run_index) +
                          "length(MESA_dirs)={}".format(lenMESA_dirs))

        #general fix for termination_flag in case of initial RLO in binaries
        if binary_grid and initial_RLO_fix:
            #create list of already detected initial RLO
            detected_initial_RLO = get_detected_initial_RLO(self)
            colnames = ["termination_flag_1", "termination_flag_2",
                        "interpolation_class"]
            valtoset = ["forced_initial_RLO", "forced_initial_RLO",
                        "initial_MT"]
            for i in range(N_runs):
                flag1 = self.final_values[i]["termination_flag_1"]
                if flag1 != "Terminate because of overflowing initial model":
                    #use grid point data (if existing) to detect initial RLO
                    grid_point = read_initial_values(grid.runs[i].path)
                    if "star_1_mass" in grid_point:
                        mass1 = grid_point["star_1_mass"]
                    else:
                        mass1 = self.initial_values[i]["star_1_mass"]
                        warnings.warn("No star_1_mass in "+grid.runs[i].path)
                    if "star_2_mass" in grid_point:
                        mass2 = grid_point["star_2_mass"]
                    else:
                        mass2 = self.initial_values[i]["star_2_mass"]
                        warnings.warn("No star_2_mass in "+grid.runs[i].path)
                    if "period_days" in grid_point:
                        period = grid_point["period_days"]
                    else:
                        period = self.initial_values[i]["period_days"]
                        warnings.warn("No period_days in "+grid.runs[i].path)
                    nearest = get_nearest_known_initial_RLO(mass1, mass2,
                                                        detected_initial_RLO)
                    if period<nearest["period_days"]:
                        #set values
                        for colname, value in zip(colnames, valtoset):
                            self.final_values[i][colname] = value
                        #copy values from nearest known system
                        for colname in ["termination_flag_3",
                                        "termination_flag_4"]:
                            if colname in nearest:
                                self.final_values[i][colname]=nearest[colname]
                        #reset the initial values to the grid point data
                        for colname, value in grid_point.items():
                            if colname in self.initial_values.dtype.names:
                                self.initial_values[i][colname] = value
                        #add initial RLO system if not added before
                        if run_included_at[i]==-1:
                            if not slim:
                                hdf5.create_group(
                                    "/grid/run{}/".format(run_index))
                            #include the run (and the input directory)
                            self.MESA_dirs.append(grid.runs[i].path)
                            run_included_at[i] = run_index
                            run_index += 1
                            #check that new MESA path is added at run_index
                            lenMESA_dirs = len(self.MESA_dirs)
                            if lenMESA_dirs!=run_index:
                                warnings.warn("Non synchronous indexing: " +
                                  "run_index={} != ".format(run_index) +
                                  "length(MESA_dirs)={}".format(lenMESA_dirs))


        self._say("Storing initial/final values and metadata to HDF5...")
        #create new array of initial and finial values with included runs
        # only and sort it by run_index
        new_initial_values = initialize_empty_array(
            np.empty(run_index, dtype=dtype_initial_values))
        new_final_values = initialize_empty_array(
            np.empty(run_index, dtype=dtype_final_values))
        for i in range(N_runs):
            #only use included runs with a valid index
            if run_included_at[i]>=0:
                #check for index range
                if run_included_at[i]>=run_index:
                    warnings.warn("run {} has a run_index out of ".format(i) +
                        "range: {}>={}".format(run_included_at[i], run_index))
                    continue
                for colname in self.initial_values.dtype.names:
                    value = self.initial_values[i][colname]
                    new_initial_values[run_included_at[i]][colname] = value
                for colname in self.final_values.dtype.names:
                    value = self.final_values[i][colname]
                    new_final_values[run_included_at[i]][colname] = value
        #replace old initial/final value array
        self.initial_values = np.copy(new_initial_values)
        self.final_values = np.copy(new_final_values)

        # Store the full table of initial_values
        hdf5.create_dataset("/grid/initial_values", data=self.initial_values,
                            **compression_args)

        # Store the full table of final_values (including termination flags)
        termination_flag_columns = np.array(termination_flag_columns)

        hdf5.create_dataset("/grid/final_values", data=self.final_values,
                            **compression_args)

        hdf5.attrs["config"] = json.dumps(str(dict(self.config)))
        hdf5.create_dataset(
            "relative_file_paths",
            data=np.asarray(self.MESA_dirs, dtype=H5_UNICODE_DTYPE),
            **compression_args)

    def add_column(self, colname, array, where="final_values", overwrite=True):
        """Add a new numerical column in the final values array."""
        arr = np.asarray(array)

        if where != "final_values":
            raise ValueError("Only adding columns to `final_values` allowed.")

        if len(self) != len(arr):
            raise ValueError("`array` has {} elements but the grid has {} runs"
                             .format(len(arr), len(self)))
        if colname in self.final_values.dtype.names:
            if overwrite:
                self.final_values[colname] = arr
            else:
                raise Exception("Column `{}` already exists in final values.".
                                format(colname))
        else:
            self.final_values = add_field(self.final_values,
                                          [(colname, arr.dtype.descr[0][1])])
            self.final_values[colname] = arr

        self.update_final_values()

    def update_final_values(self):
        """Update the final values in the HDF5 file."""
        self._reload_hdf5_file(writeable=True)
        new_dtype = []
        for dtype in self.final_values.dtype.descr:
            if (dtype[0].startswith("termination_flag")
                    or dtype[0] == "interpolation_class"
                    or "SN_type" in dtype[0] or "_state" in dtype[0]):
                dtype = (dtype[0], H5_REC_STR_DTYPE.replace("U", "S"))
            new_dtype.append(dtype)
        final_values = self.final_values.astype(new_dtype)
        del self.hdf5["/grid/final_values"]
        self.hdf5.create_dataset("/grid/final_values", data=final_values,
                                 **self.compression_args)
        # self.hdf5["/grid/final_values"] = final_values
        self._reload_hdf5_file(writeable=False)

    def _reload_hdf5_file(self, writeable=False):
        driver_args = {} if "%d" not in self.filepath else {
            "driver": "family", "memb_size": HDF5_MEMBER_SIZE}
        self.close()
        mode = "a" if writeable else "r"
        self.hdf5 = h5py.File(self.filepath, mode, **driver_args)

    def load(self, filepath=None):
        """Load up a previously created PSyGrid object from an HDF5 file.

        Parameters
        ----------
        filepath : str
            Location of the HDF5 file to be loaded. If not provided, assume
            it was defined during the initialization (argument: `filepath`).

        """
        self._say("Loading HDF5 grid...")
        # if not filepath defined, take it from the attribute
        if filepath is not None:
            self.filepath = filepath

        if self.filepath is None:
            raise ValueError("The path of the HDF5 file was not defined.")

        driver_args = {} if "%d" not in self.filepath else {
            "driver": "family", "memb_size": HDF5_MEMBER_SIZE}

        self.hdf5 = h5py.File(self.filepath, "r", **driver_args)
        hdf5 = self.hdf5
        # load initial/final_values
        self._say("\tLoading initial/final values...")
        self.initial_values = hdf5['/grid/initial_values'][()]
        self.final_values = hdf5['/grid/final_values'][()]

        # change ASCII to UNICODE in termination flags in `final_values`
        new_dtype = []
        for dtype in self.final_values.dtype.descr:
            if (dtype[0].startswith("termination_flag")
                    or dtype[0] == "interpolation_class"
                    or "SN_type" in dtype[0] or "_state" in dtype[0]):
                dtype = (dtype[0], H5_REC_STR_DTYPE.replace("S", "U"))
            new_dtype.append(dtype)
        self.final_values = self.final_values.astype(new_dtype)

        # load MESA dirs
        self._say("\tAcquiring paths to MESA directories...")
        self.MESA_dirs = list(hdf5["/relative_file_paths"])

        # load config & assert that run indices are 0-indexed & sequential
        # load configuration
        self._say("\tGetting configuration metadata...")
        configuration_string = hdf5.attrs["config"]
        config = ast.literal_eval(ast.literal_eval(configuration_string))
        self.config = ConfigFile()
        self.config.update(config)

        # enumerate the runs included, and ensure that:
        #     (i)   are as many as the MESA_dirs, or 0 (slim psyrid)
        #    (ii)   sequential and starting from 0
        self._say("\tEnumerating runs and checking integrity of grid...")
        n_expected = len(self.MESA_dirs)
        included = np.zeros(n_expected, dtype=bool)
        for key in hdf5["/grid"].keys():
            if key.startswith("run") and key != "run":
                index = int(key[3:])
                if index < 0:
                    raise KeyError("Negative index {} does not make sense".
                                   format(index))
                if index >= n_expected:
                    raise KeyError("More runs than MESA dirs? Gaps?")
                if included[index]:
                    raise KeyError("Duplicate key {}?".format(key))
                included[index] = True

        if not np.any(included):    # if slim grid... no problem!
            self.n_runs = 0
        elif np.all(included):
            self.n_runs = n_expected
        else:
            raise KeyError("Some runs are missing from the HDF5 grid.")

        self._say("\tDone.")

    def close(self):
        """Close the HDF5 file if open."""
        if hasattr(self, "hdf5") and self.hdf5 is not None:
            self.hdf5.close()
            self.hdf5 = None

    def __str__(self):
        """Return the status of the PSyGrid."""
        ret = "PSyGrid instance:\n"

        if self.filepath is None:
            ret += "\tNo HDF5 file path. PSyGrid instance is likely empty.\n"
            return ret

        ret += "\tLoaded from: {}\n\n".format(self.filepath)

        # Check if instance has loaded up the .npz data
        if self.n_runs == 0:
            ret += "\tNo runs found in the grid (empty or 'slim' version).\n"
        else:
            ret += "\t{} runs found. They include:\n".format(self.n_runs)
            # Cycle through runs to count the number of data files

            def one_if_ok(data):
                """1 if data is not None and has at least one row, else 0."""
                return 0 if data is None or len(data) == 0 else 1

            N_binary_history = 0
            N_history1 = 0
            N_history2 = 0
            N_final_profile1 = 0
            N_final_profile2 = 0
            for run in self:
                N_binary_history += one_if_ok(run.binary_history)
                N_history1 += one_if_ok(run.history1)
                N_history2 += one_if_ok(run.history2)
                N_final_profile1 += one_if_ok(run.final_profile1)
                N_final_profile2 += one_if_ok(run.final_profile2)
            ret += "\t\t{} binary_history files.\n".format(N_binary_history)
            ret += "\t\t{} history1 files.\n".format(N_history1)
            ret += "\t\t{} history2 files.\n".format(N_history2)
            ret += "\t\t{} final_profile1 files.\n".format(N_final_profile1)
            ret += "\t\t{} final_profile2 files.\n".format(N_final_profile2)

            if run.binary_history is not None:
                ret += "\nColumns in binary_history: {}\n".format(
                    run.binary_history.dtype.names)
            if run.history1 is not None:
                ret += "\nColumns in history1: {}\n".format(
                    run.history1.dtype.names)
            if run.history2 is not None:
                ret += "\nColumns in history2: {}\n".format(
                    run.history2.dtype.names)
            if run.final_profile1 is not None:
                ret += "\nColumns in final_profile1: {}\n".format(
                    run.final_profile1.dtype.names)
            if run.final_profile2 is not None:
                ret += "\nColumns in final_profile2: {}\n".format(
                    run.final_profile1.dtype.names)
        ret += "\n"

        # Print out initial values array parameters
        if self.initial_values is None:
            ret += "Initial values array is not loaded.\n\n"
        else:
            ret += "Columns in initial values:\n"
            ret += str(self.initial_values.dtype.names) + "\n\n"

        # Print out final values array parameters
        if self.final_values is None:
            ret += "Final values array is not loaded.\n\n"
        else:
            ret += "Columns in final values:\n"
            ret += str(self.final_values.dtype.names) + "\n\n"

        # Print out configuration file parameters
        ret += "Configuration:\n"
        config_string = str(self.config)
        if config_string == "":
            ret += "\t(empty)\n\n"
        else:
            # separate config items with lines that start with a tab
            ret += "\t" + "\n\t".join(config_string.split("\n")) + "\n\n"

        # Print out examples of MESA dirs (if grid is not empty)
        n_dirs = len(self.MESA_dirs)
        if n_dirs != 0 or self.n_runs != 0:
            ret += "Relative paths to MESA run directories:\n"
            if n_dirs == 0:
                ret += "\tNo MESA directory found.\n"
            elif n_dirs == 1:
                ret += "\t{}\n".format(self.MESA_dirs[0])
            elif n_dirs == 2:
                ret += "\t{}\n".format(self.MESA_dirs[0])
                ret += "\t{}\n".format(self.MESA_dirs[1])
            else:
                ret += "\t{}\n".format(self.MESA_dirs[0])
                ret += "\t...({} other directories)\n".format(n_dirs - 2)
                ret += "\t{}\n".format(self.MESA_dirs[-1])

        return ret

    def __getitem__(self, idx):
        """Return a PSyRunView instance for the run with index `idx`."""
        if idx not in self:
            raise KeyError("Index {} out of bounds.".format(idx))
        return PSyRunView(self, idx)

    def get_pandas_initial_final(self):
        """Convert the initial/final values into a single Pandas dataframe."""
        df = pd.DataFrame()
        for key in self.initial_values.dtype.names:
            new_col_name = "initial_" + key
            df[new_col_name] = self.initial_values[key]

        for key in self.final_values.dtype.names:
            new_col_name = "final_" + key
            df[new_col_name] = self.final_values[key]
        return df

    def __len__(self):
        """Return the number of runs in the grid."""
        return self.n_runs

    def __contains__(self, index):
        """Return True if run with index `index` is in the grid."""
        return 0 <= index < self.n_runs

    def __iter__(self):
        """Allow iteration of runs (`for` loops, lists, sets)."""
        return PSyGridIterator(self)

    def __del__(self):
        """Destructor of the object, closing the HDF5 file if opened."""
        # IMPORTANT: do not delete this - it is necessary in interactive mode
        self.close()

    def rerun(self, path_to_file='./', runs_to_rerun=None,
              termination_flags=None, new_mesa_flag=None):
        """Create a CSV file with the PSyGrid initial values to rerun.

        This methods allows you to create a CSV file with the psygrid initial
        values you want to rerun.

        Parameters
        ----------
        path_to_file : str
            The path to the directory where the new `grid.csv` file will be
            saved. If the directory does not exist it will be created.
        runs_to_rerun : integer array
            Array containing the indecies of the psygrid runs you want to rerun
            e.g., runs_to_rerun = np.array([2,3])
        termination_flags : str
            The runs with this termination flag will be rerun.
            e.g. termination_flags='max_number_retries'
        new_mesa_flag : dict
            Dictionary of flags with their value to add as extra columns to the
            `grid.csv`. The user can specify any arbitrary amount of flags.
            e.g. new_mesa_flag = {'varcontrol_target': 0.01}

        """
        # check that the `path_to_file` exists
        if not os.path.exists(path_to_file):
            os.makedirs(path_to_file)

        if runs_to_rerun is not None and termination_flags is None:
            n_runs_to_rerun = len(runs_to_rerun)

            # find the key of initial values to save
            initial_values = {}
            for key in self.initial_values.dtype.names:
                initial_values[key] = []

            # find the value of initial values to save
            for i in runs_to_rerun.tolist():
                for key in initial_values:
                    initial_values[key].append(self.initial_values[key][i])

            # replace star_1_mass, star_2_mass, period_days, Z
            NDIG = 10 # rounding matches initial point rounding
            if 'star_1_mass' in self.initial_values.dtype.names:
                initial_values['m1'] = np.around(initial_values['star_1_mass'],
                                                 NDIG)
            if 'star_2_mass' in self.initial_values.dtype.names:
                initial_values['m2'] = np.around(initial_values['star_2_mass'],
                                                 NDIG)
            if 'period_days' in self.initial_values.dtype.names:
                initial_values['initial_period_in_days'] = np.around(
                                            initial_values['period_days'],NDIG)
            MESA_dir_name = self.MESA_dirs[0].decode("utf-8")
            if  'initial_z' in MESA_dir_name:
                initial_values['initial_z'] = np.around(initial_values['Z'],
                                                        NDIG)
            if  'Zbase' in MESA_dir_name:
                initial_values['Zbase'] = np.around(initial_values['Z'], NDIG)
            if  'new_Z' in MESA_dir_name:
                initial_values['new_Z'] = np.around(initial_values['Z'], NDIG)
            for key in self.initial_values.dtype.names:
                if key not in ['m1', 'm2', 'initial_period_in_days',
                               'Zbase', 'new_Z', 'initial_z']:
                    del initial_values[key]

            # add new_mesa_flag
            if new_mesa_flag is not None:
                for key in new_mesa_flag.keys():
                    initial_values[key] = [new_mesa_flag[key]]*n_runs_to_rerun

            # create the CSV file
            with open(path_to_file+'grid.csv', 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(initial_values.keys())
                for i in range(n_runs_to_rerun):
                    writer.writerow(
                        [initial_values[key][i]for key in initial_values])
        elif termination_flags is not None and runs_to_rerun is None:
            if isinstance(termination_flags, str):
                rerun_flags = [termination_flags]
            else:
                try:
                    rerun_flags = list(termination_flags)
                except TypeError as err:
                    msg = "`termination_flags` should be a str of a list."
                    raise ValueError(msg) from err

            # find the key of initial values to save
            initial_values = {}
            for key in self.initial_values.dtype.names:
                initial_values[key] = []

            # find the value of initial values to save
            n_runs_to_rerun = 0
            n = self.initial_values.shape[0]
            for i in range(n):
                # this should be accessed with key when possible
                if self.final_values[i]['termination_flag_1'] in rerun_flags:
                    n_runs_to_rerun += 1
                    for key in initial_values:
                        initial_values[key].append(self.initial_values[key][i])

            # replace star_1_mass, star_2_mass, period_days, Z
            NDIG = 10 # rounding matches initial point rounding
            if 'star_1_mass' in self.initial_values.dtype.names:
                initial_values['m1'] = np.around(initial_values['star_1_mass'],
                                                 NDIG)
            if 'star_2_mass' in self.initial_values.dtype.names:
                initial_values['m2'] = np.around(initial_values['star_2_mass'],
                                                 NDIG)
            if 'period_days' in self.initial_values.dtype.names:
                initial_values['initial_period_in_days'] = np.around(
                                            initial_values['period_days'],NDIG)
            MESA_dir_name = self.MESA_dirs[0].decode("utf-8")
            if  'initial_z' in MESA_dir_name:
                initial_values['initial_z'] = np.around(initial_values['Z'],
                                                        NDIG)
            if  'Zbase' in MESA_dir_name:
                initial_values['Zbase'] = np.around(initial_values['Z'], NDIG)
            if  'new_Z' in MESA_dir_name:
                initial_values['new_Z'] = np.around(initial_values['Z'], NDIG)
            for key in self.initial_values.dtype.names:
                if key not in ['m1', 'm2', 'initial_period_in_days', 'Zbase',
                               'new_Z', 'initial_z']:
                    del initial_values[key]

            # add new_mesa_flag
            if new_mesa_flag is not None:
                for key in new_mesa_flag.keys():
                    initial_values[key] = [new_mesa_flag[key]]*n_runs_to_rerun

            # create the CSV file
            with open(path_to_file+'grid.csv', 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(initial_values.keys())
                for i in range(n_runs_to_rerun):
                    writer.writerow(
                        [initial_values[key][i]for key in initial_values])
        else:
            raise ValueError("Choose either the runs manually, or "
                             "indicate the termination flag(s).")

    def plot2D(self, x_var_str, y_var_str, z_var_str=None,
               termination_flag='termination_flag_1',
               grid_3D=None, slice_3D_var_str=None, slice_3D_var_range=None,
               grid_4D=None, slice_4D_var_str=None, slice_4D_var_range=None,
               extra_grid=None, slice_at_RLO=False,
               MARKERS_COLORS_LEGENDS=None,
               verbose=False, **kwargs):
        """Plot a 2D slice of x_var_str vs y_var_str of one or more runs.

        Parameters
        ----------
        x_var_str : str
            String of the initial value to plot on the x axis. Allowed strings
            are `psygrid.initial_values.dtype.names`.
        y_var_str : str
            String of the initial value to plot on the y axis. Allowed strings
            are `psygrid.initial_values.dtype.names`.
        z_var_str : str
            String of the initial value to plot on the z axis (displayed as
            a color). Allowed strings are
            `psygrid.final_values.dtype.names`, `psygrid.history1.dtype.names`,
            `psygrid.binary_history.dtype.names`.
        termination_flag : str
            Termination flag to display, allowed values are:
            "termination_flag_1", "termination_flag_2", "termination_flag_3",
            "termination_flag_4", "all".
        grid_3D : bool
            If `True`, the psygrid object is a 3D grid and needs to be sliced.
        slice_3D_var_str : str
            Variable along which the 3D space will be sliced. Allowed values
            are `psygrid.initial_values.dtype.names`.
        slice_3D_var_range : tuple
            Range between which you want to slice the variable slice_3D_var_str
            e.g., `(2.5,3.)`.
        grid_4D : bool
            If `True`, the psygrid object is a 4D grid and needs to be sliced.
        slice_4D_var_str : str
            Variable along which the 4D space will be sliced. Allowed values
            are `psygrid.initial_values.dtype.names`.
        slice_4D_var_range : toople
            Range between which you want to slice the variable slice_4D_var_str
            e.g., `(2.5,3.)`.
        extra_grid : object or array of objects
            If subset of the grid was rerun a or an extention was added, one
            can overlay the new psygrid by passing it here.
        slice_at_RLO : bool
            If `True`, the object plots the tracks until onset of Roche Lobe
            overflow.
        MARKERS_COLORS_LEGENDS : dict
            Each termination flag is associated with a marker shape, size,
            color and label (cf. `MARKERS_COLORS_LEGENDS` in
            `plot_defaults.py`).
        verbose : bool
            If `True`, the object reports by printing to standard output.
        **kwargs : dict
            Dictionary containing extra visualisation options (cf.
            `PLOT_PROPERTIES` in `plot_defaults.py`.

        """
        plot = plot2D(psygrid=self,
                      x_var_str=x_var_str,
                      y_var_str=y_var_str,
                      z_var_str=z_var_str,
                      termination_flag=termination_flag,
                      grid_3D=grid_3D,
                      slice_3D_var_str=slice_3D_var_str,
                      slice_3D_var_range=slice_3D_var_range,
                      grid_4D=grid_4D,
                      slice_4D_var_str=slice_4D_var_str,
                      slice_4D_var_range=slice_4D_var_range,
                      extra_grid=extra_grid,
                      slice_at_RLO=slice_at_RLO,
                      MARKERS_COLORS_LEGENDS=MARKERS_COLORS_LEGENDS,
                      verbose=verbose,
                      **kwargs)
        plot()

    def plot(self, idx, x_var_str, y_var_str, z_var_str=None,
             history='binary_history', verbose=False, **kwargs):
        """Plot a 1D track of x_var_str vs y_var_str.

        Parameters
        ----------
        idx : int or list of int
            Index or indices of runs to plot.
        x_var_str : str
            String of values to plot on the x axis. Allowed strings are the
            one in `psygrid.history.dtype.names` where "history" needs to be
            chosen accordingly.
        y_var_str : str or list of str
            String or list of stringvalues to plot on the y axis. Allowed
            strings are the one in `psygrid.history.dtype.names` where
            "history" needs to be chosen accordingly.
        z_var_str : str
            String of values to plot on the z axis (displayed with a color).
            Allowed strings are the one in `psygrid.history.dtype.names` where
            "history" needs to be chosen accordingly.
        history : str
            The x, y, z variables are read from either: "binary_history",
            "history1", "history2".
        verbose : bool
            If `True`, the object reports by printing to standard output.
        **kwargs : dict
            Dictionary containing extra visualisation options
            (cf. `PLOT_PROPERTIES` in `plot_defaults.py`).

        """
        if isinstance(idx, list):
            runs = []
            for i in idx:
                runs.append(self[i])
        elif isinstance(idx, int):
            runs = [self[idx]]
        else:
            raise ValueError('Invalid idx = {}!'.format(idx))

        plot = plot1D(run=runs,
                      x_var_str=x_var_str,
                      y_var_str=y_var_str,
                      z_var_str=z_var_str,
                      history=history,
                      HR=False,
                      verbose=verbose,
                      **kwargs)
        plot()

    def HR(self, idx, history='history1', states=False, verbose=False,
           **kwargs):
        """Plot the HR diagram of one or more runs.

        Parameters
        ----------
        idx : int or list of int
            Index or indices of runs to plot.
        history : str
            Which history is going to be used. The options are:
            "binary_history", "history1", or "history2".
        states : bool
            If true the HR diagram shows the stellar state with a color map.
        verbose : bool
            If `True`, the object reports by printing to standard output.
        **kwargs : dict
            Dictionary containing extra visualisation options
            (cf. `PLOT_PROPERTIES` in `plot_defaults.py`).

        """
        if isinstance(idx, list):
            runs = []
            for i in idx:
                runs.append(self[i])
        elif isinstance(idx, int):
            runs = [self[idx]]
        else:
            raise ValueError('Invalid idx = {}!'.format(idx))

        if states:
            from posydon.binary_evol.singlestar import SingleStar
            star_states = []
            for run in runs:
                star = SingleStar.from_run(run, history=True, profile=False)
                star_states.append(check_state_of_star_history_array(
                    star, N=len(star.mass_history)))
        else:
            star_states = None

        plot = plot1D(run=runs,
                      x_var_str=None,
                      y_var_str=None,
                      history=history,
                      star_states=star_states,
                      HR=True,
                      verbose=verbose,
                      **kwargs)
        plot()

    def __eq__(self, other, verbose=True):
        """Check the equality (in terms of the data) of two PSyGrid objects."""
        def say(msg):
            if verbose:
                print("COMPARISON:", msg)

        if not isinstance(other, PSyGrid):
            say("Only PSyGrid instances should be compared.")
            return False

        if len(self) != len(other):
            say("The grids do not contain the same number of runs ({} != {})".
                format(len(self), len(other)))

        for prop in PROPERTIES_TO_BE_CONSISTENT:
            val1, val2 = self.config[prop], other.config[prop]
            if val1 != val2:
                say("Property `{}` is not the same ({} != {})".
                    format(prop, val1, val2))
                return False

        # just report differences for parameters that are not data-relevant
        for prop in ["compression", "description"]:
            val1, val2 = self.config[prop], other.config[prop]
            if val1 != val2:
                say("Property `{}` is not the same ({} != {})".
                    format(prop, val1, val2))

        for tablename in ["initial_values", "final_values"]:
            arr1, arr2 = getattr(self, tablename), getattr(other, tablename)
            columns1, columns2 = arr1.dtype.names, arr2.dtype.names

            if np.all(arr1 == arr2):
                continue
            if len(columns1) != len(columns2):
                say("Number of columns in `{}` do not match ({} != {})".
                    format(tablename, len(columns1), len(columns2)))
                return False

            if columns1 != columns2:
                say("Columns in `{}` do not match:\n{}\n!=\n{}\n".format
                    (tablename, "\n".join(columns1), "\n".join(columns2)))
                return False

            for colname in columns1:
                data1, data2 = arr1[colname], arr2[colname]
                if len(data1) != len(data2):
                    say("Column `{}` in `{}` not of same length.".
                        format(colname, tablename))
                if np.any(data1 != data2):
                    for val1, val2 in zip(data1, data2):
                        if ((val1 == val2) or (val1 is None and val2 is None)
                                or (np.isnan(val1) and np.isnan(val2))):
                            continue
                        say("Column `{}` in `{}` is not the same ({} != {}).".
                            format(colname, tablename, val1, val2))
                        return False

        for i, (run1, run2) in enumerate(zip(self, other)):
            for tablename in ["binary_history", "history1", "history2",
                              "final_profile1", "final_profile2"]:
                arr1, arr2 = getattr(run1, tablename), getattr(run2, tablename)
                if arr1 is None and arr2 is None:
                    continue        # tables missing for both runs
                if arr1 is None or arr2 is None:
                    which = "1st" if arr1 is None else "2nd"
                    say("Table `{}` for run `{}` missing in {} grid.".
                        format(tablename, i, which))
                    return False
                if np.all(arr1 != arr2):
                    say("Table `{}` for run `{}` is not the same.".
                        format(tablename, i))

        return True


class PSyGridIterator:
    """Class for iterating runs (i.e., PSyRunView instances) in a PSyGrid."""

    def __init__(self, grid):
        """Initialize by pointing to the PSyGrid object."""
        self._grid = grid
        self._index = 0

    def __next__(self):
        """Return current run and move to the next."""
        if self._index >= len(self._grid):
            raise StopIteration
        item = self._grid[self._index]
        self._index += 1
        return item


class PSyRunView:
    """Access runs in HDF5 grid without opening all data.

    Example
    -------
    mygrid = PSyGrid("thegrid.h5")  # opens a stored grid
    myrun = mygrid[0]               # get a "view" of the first run in the grid
    print(mygrid[1].history1.star_age)

    t1 = myrun["binary_history"]["age"]
    t2 = myrun["binary_history"].age
    t3 = myrun.binary_history["age"]
    t4 = myrun.binary_history.age

    print(mygrid[1].history1.star_age)
    """

    def __init__(self, psygrid, run_index):
        """Initialize by linking to a PSyGrid object and setting run index."""
        self.psygrid = psygrid
        self.index = run_index

    def _hdf5_key(self):
        return "/grid/run{}".format(self.index)

    def __getitem__(self, key):
        """Get a table for a specific run using its name as in ["name"]."""
        if key not in VALID_KEYS:
            raise KeyError("Key {} not in list of valid keys.".format(key))

        if key == "initial_values":
            return self.psygrid.initial_values[self.index]
        elif key == "final_values":
            return self.psygrid.final_values[self.index]

        hdf5 = self.psygrid.hdf5
        if hdf5 is None:
            raise Exception("The HDF5 file is not open.")
        fullkey = self._hdf5_key() + "/" + key
        try:
            return hdf5[fullkey][()]
        except KeyError:
            if key in VALID_KEYS:  # key is valid, so the data is just missing
                return None
            raise KeyError("There is no key '{}' in '{}'".
                           format(fullkey, self.psygrid.filepath))

    def __getattr__(self, key):
        """Enable the ability of using `.table1` instead of ["table1"]."""
        return self[key]

    def __str__(self):
        """Return a summary of the PSyRunView in a string."""
        return "View of the run {} in the file '{}' at key '{}'".format(
            self.index, self.psygrid.filepath, self._hdf5_key())


def downsample_history(bh, h1, h2, params):
    """Downsample history data of a run.

    Parameters
    ----------
    bh : record array
        The binary history table.
    h1 : record array
        The history of star 1.
    h2 : record array
        The history of star 2.
    params : ConfigFile
        A ConfigFile instance holding the downsampling parameters.

    Returns
    -------
    tuple of record arrays
        The three record arrays corresponding to the downsample binary,
        star 1, and star2 histories.

    """
    max_err, exclude = params["history_DS_error"], params["history_DS_exclude"]

    if max_err is None or (bh is None and h1 is None and h2 is None):
        return bh, h1, h2

    # select the data based on which the downsampling will be computed
    number_of_rows = None
    dependent = []
    for history_table in [bh, h1, h2]:
        if history_table is not None:
            if number_of_rows is None:
                number_of_rows = len(history_table)
            else:
                if len(history_table) != number_of_rows:
                    raise Exception("Unequal numbers of rows in histories.")

            for column_name in history_table.dtype.names:
                if column_name not in exclude:
                    dependent.append(history_table[column_name])
    dependent = np.array(dependent).T

    if params["binary"]:
        independent = bh["age"]
    else:
        independent = h1["star_age"]

    # find which rows to keep for having optimum downsampling
    TD = TrackDownsampler(independent, dependent)
    TD.find_downsample(max_err)

    # ... and return only these rows
    def correct_shape(array):
        return array[0] if len(array.shape) == 2 else array

    return (correct_shape(bh[TD.keep]) if bh is not None else None,
            correct_shape(h1[TD.keep]) if h1 is not None else None,
            correct_shape(h2[TD.keep]) if h2 is not None else None)


def downsample_profile(profile, params):
    """Downsample a profile table.

    Parameters
    ----------
    profile : record array
        The profile data with column names defined.
    params : ConfigFile
        A ConfigFile instance holding the downsampling parameters.

    Returns
    -------
    record array
        The downsampled profile.

    """
    max_err, exclude = params["profile_DS_error"], params["profile_DS_exclude"]
    max_interval = params["profile_DS_interval"]

    if max_err is None or profile is None:
        return profile

    # select the data based on which the downsampling will be computed
    dependent = []
    for column_name in profile.dtype.names:
        if column_name not in exclude:
            dependent.append(profile[column_name][::-1])  # invert order
    dependent = np.array(dependent).T
    independent = profile["mass"][::-1]

    # find which rows to keep for having optimum downsampling
    TD = TrackDownsampler(independent, dependent)
    TD.find_downsample(max_err, max_interval)
    # ... and return only these rows
    return profile[TD.keep[::-1]]    # invert order again


PROPERTIES_ALLOWED = {"format": "hdf5"}

PROPERTIES_TO_BE_SET = ["compression", "description"]

PROPERTIES_TO_BE_NONE = {
    "max_number_of_runs": None,
    "star1_history_saved_columns": None,
    "star2_history_saved_columns": None,
    "binary_history_saved_columns": None,
    "star1_profile_saved_columns": None,
    "star2_profile_saved_columns": None,
    "initial_value_columns": None,
    "final_value_columns": None,
}

PROPERTIES_TO_BE_CONSISTENT = ["binary", "eep", "start_at_RLO",
                               "initial_RLO_fix",
                               "He_core_fix", "history_DS_error",
                               "history_DS_exclude", "profile_DS_error",
                               "profile_DS_exclude", "profile_DS_interval"]

ALL_PROPERTIES = (list(PROPERTIES_ALLOWED.keys()) + PROPERTIES_TO_BE_CONSISTENT
                  + list(PROPERTIES_TO_BE_NONE.keys()) + PROPERTIES_TO_BE_SET)


def join_grids(input_paths, output_path,
               compression="gzip9", description="joined", verbose=True):
    """Join two or more PSyGrid HDF5 files into one."""
    def say(something):
        if verbose:
            print(something)

    # open the grids, and figure out the configuration of the new joined grid
    grids = []
    newconfig = None
    initial_dtype = None
    final_dtype = None
    for input_path in input_paths:
        say("Checking {}".format(input_path))
        say("    opening and getting configuration...")
        grid = PSyGrid(input_path)
        curconfig = grid.config

        # check we did not forget about a property, and their consistency
        assert len(ALL_PROPERTIES) == len(curconfig.keys())
        for key in curconfig.keys():
            assert key in ALL_PROPERTIES
        if newconfig is None:   # this is the first grid
            say("    copying configuration...")
            newconfig = curconfig.deepcopy()
            initial_dtype = grid.initial_values.dtype
            final_dtype = grid.final_values.dtype
        else:
            say("    checking consistency of configuration and dtypes...")
            for prop in PROPERTIES_TO_BE_CONSISTENT:
                newvalue, curvalue = newconfig[prop], curconfig[prop]
                if newvalue != curvalue:
                    raise Exception("Inconsistent value for `{}`: {} != {}".
                                    format(prop, curvalue, newvalue))
            if initial_dtype != grid.initial_values.dtype:
                raise Exception("Initial values dtype's do not match.")
            if final_dtype != grid.final_values.dtype:
                raise Exception("Final values dtype's do not match.")
        grids.append(grid)

    say("Finalizing grid properties...")
    # ensure supported property values
    for prop, allowed_value in PROPERTIES_ALLOWED.items():
        if newconfig[prop] != allowed_value:
            raise ValueError("Only grid with `{}`={} can be joined.".
                             format(prop, allowed_value))
    # set fixed property values
    newconfig.update(PROPERTIES_TO_BE_NONE)
    # new values set by used
    assert("compression" in PROPERTIES_TO_BE_SET
           and "description" in PROPERTIES_TO_BE_SET)
    newconfig["compression"] = compression
    newconfig["description"] = description

    # for each set of initial parameters, indicate the grid and the run to use
    # (this handles reruns, substituting older runs).
    say("Inferring which runs to include, from which grid...")
    initial_params = {}
    n_substitutions = 0
    for grid_index, grid in enumerate(grids):
        for run_index, dir_path in enumerate(grid.MESA_dirs):
            # get the parameters part from the dir name
            params_from_path = initial_values_from_dirname(dir_path)
            if params_from_path in initial_params:
                n_substitutions += 1
            initial_params[params_from_path] = (grid_index, run_index)
    say("    {} substituions detected.".format(n_substitutions))
    say("    {} runs to be joined.".format(len(initial_params)))

    say("Opening new file...")
    # open new HDF5 file and start copying runs
    driver_args = {} if "%d" not in output_path else {
        "driver": "family", "memb_size": HDF5_MEMBER_SIZE}

    # prepare for loss-less compression
    compression = newconfig["compression"]
    compression = "" if compression is None else compression.lower()
    if compression == "":
        compression_args = {}
    elif compression == "lzf":
        compression_args = {"compression": "lzf"}
    elif compression.startswith("gzip"):
        if compression == "gzip":
            compression_args = {"compression": "gzip"}
        else:
            if len(compression) > 5 or not compression[4].isdigit():
                raise ValueError("GZIP compression level must be 0-9.")
            compression_level = int(compression[4])
            compression_args = {"compression": "gzip",
                                "compression_opts": compression_level}
    else:
        raise ValueError("Unknown compression `{}`.".format(compression))

    with h5py.File(output_path, "w", **driver_args) as new_grid:
        say("Copying runs...")
        new_mesa_dirs = []
        new_initial_values = []
        new_final_values = []
        new_index = 0
        for grid_index, run_index in initial_params.values():
            grid = grids[grid_index]
            new_mesa_dirs.append(grid.MESA_dirs[run_index])
            new_initial_values.append(grid.initial_values[run_index])
            new_final_values.append(grid.final_values[run_index])

            for subarray in ['binary_history', 'history1', 'history2',
                             'final_profile1', 'final_profile2']:
                # ensure group is created even if no data to be written
                new_grp = "/grid/run{}/".format(new_index)
                if new_grp not in new_grid:
                    new_grid.create_group(new_grp)

                # now copy the data
                new_key = "/grid/run{}/{}".format(new_index, subarray)
                old_key = "/grid/run{}/{}".format(run_index, subarray)
                try:
                    old_data = grid.hdf5[old_key]
                    new_grid.create_dataset(new_key, data=old_data,
                                            **compression_args)
                except KeyError:
                    pass

            new_index += 1

        say("Writing initial/final values and metadata...")

        new_mesa_dirs = np.array(new_mesa_dirs, dtype=H5_UNICODE_DTYPE)
        new_initial_values = np.array(new_initial_values, dtype=initial_dtype)
        new_final_dtype = []
        for dtype in final_dtype.descr:
            if (dtype[0].startswith("termination_flag")
                    or dtype[0] == "interpolation_class"
                    or "SN_type" in dtype[0] or "_state" in dtype[0]):
                dtype = (dtype[0], H5_REC_STR_DTYPE.replace("U", "S"))
            new_final_dtype.append(dtype)
        new_final_values = np.array(new_final_values, dtype=new_final_dtype)

        new_grid.attrs["config"] = json.dumps(str(dict(newconfig)))
        new_grid.create_dataset("relative_file_paths", data=new_mesa_dirs,
                                **compression_args)
        new_grid.create_dataset("/grid/initial_values",
                                data=new_initial_values,
                                **compression_args)
        new_grid.create_dataset("/grid/final_values",
                                data=new_final_values,
                                **compression_args)
    say("Grids have been successfully joined.")
