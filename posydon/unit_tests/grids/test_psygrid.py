"""Unit tests of posydon/grids/psygrid.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.grids.psygrid as totest
# aliases
np = totest.np
os = totest.os
h5py = totest.h5py

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns
from inspect import isclass, isroutine
from posydon.utils.posydonwarning import ReplaceValueWarning

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['ALL_PROPERTIES', 'Catch_POSYDON_Warnings', 'ConfigFile',\
                    'DEFAULT_BINARY_HISTORY_COLS', 'DEFAULT_EEP_HISTORY_COLS',\
                    'DEFAULT_HISTORY_DS_EXCLUDE', 'DEFAULT_PROFILE_COLS',\
                    'DEFAULT_PROFILE_DS_EXCLUDE',\
                    'DEFAULT_SINGLE_HISTORY_COLS',\
                    'DEFAULT_STAR_HISTORY_COLS', 'EEP_FILE_EXTENSIONS',\
                    'EXTRA_COLS_DS_EXCLUDE',\
                    'EXTRA_STAR_COLS_AT_HE_DEPLETION', 'GRIDPROPERTIES',\
                    'GridReader', 'H5_REC_STR_DTYPE', 'H5_UNICODE_DTYPE',\
                    'HDF5_MEMBER_SIZE', 'IgnoreReason', 'N_FLAGS',\
                    'N_FLAGS_SINGLE', 'PROPERTIES_ALLOWED',\
                    'PROPERTIES_TO_BE_CONSISTENT', 'PROPERTIES_TO_BE_NONE',\
                    'PROPERTIES_TO_BE_SET', 'PSyGrid', 'PSyGridIterator',\
                    'PSyRunView', 'Pwarn', 'TERMINATION_FLAG_COLUMNS',\
                    'TERMINATION_FLAG_COLUMNS_SINGLE',\
                    'THRESHOLD_CENTRAL_ABUNDANCE',\
                    'THRESHOLD_CENTRAL_ABUNDANCE_LOOSE_C', 'TrackDownsampler',\
                    'VALID_KEYS', 'WARN_VALUES', '__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__',\
                    'add_field', 'ast', 'check_state_from_history',\
                    'check_state_of_star_history_array', 'downsample_history',\
                    'downsample_profile', 'fix_He_core',\
                    'get_detected_initial_RLO', 'get_flag_from_MESA_output',\
                    'get_flags_from_MESA_run', 'get_i_He_depl',\
                    'get_nearest_known_initial_RLO', 'glob', 'h5py',\
                    'infer_interpolation_class', 'infer_star_state',\
                    'initial_values_from_dirname', 'initialize_empty_array',\
                    'join_grids', 'join_lists', 'json', 'keep_after_RLO',\
                    'keep_till_central_abundance_He_C', 'np',\
                    'orbital_separation_from_period', 'os', 'pd', 'plot1D',\
                    'plot2D', 'read_EEP_data_file', 'read_MESA_data_file',\
                    'read_initial_values', 'scrub', 'tqdm']
        assert dir(totest) == elements, "There might be added or removed "\
                                        + "objects without an update on the "\
                                        + "unit test."

    def test_instance_HDF5_MEMBER_SIZE(self):
        assert isinstance(totest.HDF5_MEMBER_SIZE, (int, float)),\
               "HDF5_MEMBER_SIZE is of type: "\
               + str(type(totest.HDF5_MEMBER_SIZE))

    def test_instance_H5_UNICODE_DTYPE(self):
        assert isinstance(totest.H5_UNICODE_DTYPE, np.dtype),\
               "H5_UNICODE_DTYPE is of type: "\
               + str(type(totest.H5_UNICODE_DTYPE))

    def test_instance_H5_REC_STR_DTYPE(self):
        assert isinstance(totest.H5_REC_STR_DTYPE, str),\
               "H5_REC_STR_DTYPE is of type: "\
               + str(type(totest.H5_REC_STR_DTYPE))

    def test_instance_VALID_KEYS(self):
        assert isinstance(totest.VALID_KEYS, list),\
               "VALID_KEYS is of type: " + str(type(totest.VALID_KEYS))

    def test_instance_WARN_VALUES(self):
        assert isinstance(totest.WARN_VALUES, list),\
               "WARN_VALUES is of type: " + str(type(totest.WARN_VALUES))

    def test_instance_N_FLAGS(self):
        assert isinstance(totest.N_FLAGS, int),\
               "N_FLAGS is of type: " + str(type(totest.N_FLAGS))

    def test_instance_N_FLAGS_SINGLE(self):
        assert isinstance(totest.N_FLAGS_SINGLE, int),\
               "N_FLAGS_SINGLE is of type: " + str(type(totest.N_FLAGS_SINGLE))

    def test_instance_TERMINATION_FLAG_COLUMNS(self):
        assert isinstance(totest.TERMINATION_FLAG_COLUMNS, list),\
               "TERMINATION_FLAG_COLUMNS is of type: "\
               + str(type(totest.TERMINATION_FLAG_COLUMNS))

    def test_instance_TERMINATION_FLAG_COLUMNS_SINGLE(self):
        assert isinstance(totest.TERMINATION_FLAG_COLUMNS_SINGLE, list),\
               "TERMINATION_FLAG_COLUMNS_SINGLE is of type: "\
               + str(type(totest.TERMINATION_FLAG_COLUMNS_SINGLE))

    def test_instance_DEFAULT_BINARY_HISTORY_COLS(self):
        assert isinstance(totest.DEFAULT_BINARY_HISTORY_COLS, list),\
               "DEFAULT_BINARY_HISTORY_COLS is of type: "\
               + str(type(totest.DEFAULT_BINARY_HISTORY_COLS))

    def test_instance_DEFAULT_STAR_HISTORY_COLS(self):
        assert isinstance(totest.DEFAULT_STAR_HISTORY_COLS, list),\
               "DEFAULT_STAR_HISTORY_COLS is of type: "\
               + str(type(totest.DEFAULT_STAR_HISTORY_COLS))

    def test_instance_DEFAULT_SINGLE_HISTORY_COLS(self):
        assert isinstance(totest.DEFAULT_SINGLE_HISTORY_COLS, list),\
               "DEFAULT_SINGLE_HISTORY_COLS is of type: "\
               + str(type(totest.DEFAULT_SINGLE_HISTORY_COLS))

    def test_instance_DEFAULT_EEP_HISTORY_COLS(self):
        assert isinstance(totest.DEFAULT_EEP_HISTORY_COLS, list),\
               "DEFAULT_EEP_HISTORY_COLS is of type: "\
               + str(type(totest.DEFAULT_EEP_HISTORY_COLS))

    def test_instance_DEFAULT_PROFILE_COLS(self):
        assert isinstance(totest.DEFAULT_PROFILE_COLS, list),\
               "DEFAULT_PROFILE_COLS is of type: "\
               + str(type(totest.DEFAULT_PROFILE_COLS))

    def test_instance_EXTRA_STAR_COLS_AT_HE_DEPLETION(self):
        assert isinstance(totest.EXTRA_STAR_COLS_AT_HE_DEPLETION, list),\
               "EXTRA_STAR_COLS_AT_HE_DEPLETION is of type: "\
               + str(type(totest.EXTRA_STAR_COLS_AT_HE_DEPLETION))

    def test_instance_DEFAULT_HISTORY_DS_EXCLUDE(self):
        assert isinstance(totest.DEFAULT_HISTORY_DS_EXCLUDE, list),\
               "DEFAULT_HISTORY_DS_EXCLUDE is of type: "\
               + str(type(totest.DEFAULT_HISTORY_DS_EXCLUDE))

    def test_instance_DEFAULT_PROFILE_DS_EXCLUDE(self):
        assert isinstance(totest.DEFAULT_PROFILE_DS_EXCLUDE, list),\
               "DEFAULT_PROFILE_DS_EXCLUDE is of type: "\
               + str(type(totest.DEFAULT_PROFILE_DS_EXCLUDE))

    def test_instance_EXTRA_COLS_DS_EXCLUDE(self):
        assert isinstance(totest.EXTRA_COLS_DS_EXCLUDE, list),\
               "EXTRA_COLS_DS_EXCLUDE is of type: "\
               + str(type(totest.EXTRA_COLS_DS_EXCLUDE))

    def test_instance_GRIDPROPERTIES(self):
        assert isinstance(totest.GRIDPROPERTIES, dict),\
               "GRIDPROPERTIES is of type: " + str(type(totest.GRIDPROPERTIES))

    def test_instance_PSyGrid(self):
        assert isclass(totest.PSyGrid)

    def test_instance_PSyGridIterator(self):
        assert isclass(totest.PSyGridIterator)

    def test_instance_PSyRunView(self):
        assert isclass(totest.PSyRunView)

    def test_instance_downsample_history(self):
        assert isroutine(totest.downsample_history)

    def test_instance_downsample_profile(self):
        assert isroutine(totest.downsample_profile)

    def test_instance_PROPERTIES_ALLOWED(self):
        assert isinstance(totest.PROPERTIES_ALLOWED, dict),\
               "PROPERTIES_ALLOWED is of type: "\
               + str(type(totest.PROPERTIES_ALLOWED))

    def test_instance_PROPERTIES_TO_BE_SET(self):
        assert isinstance(totest.PROPERTIES_TO_BE_SET, list),\
               "PROPERTIES_TO_BE_SET is of type: "\
               + str(type(totest.PROPERTIES_TO_BE_SET))

    def test_instance_PROPERTIES_TO_BE_NONE(self):
        assert isinstance(totest.PROPERTIES_TO_BE_NONE, dict),\
               "PROPERTIES_TO_BE_NONE is of type: "\
               + str(type(totest.PROPERTIES_TO_BE_NONE))

    def test_instance_PROPERTIES_TO_BE_CONSISTENT(self):
        assert isinstance(totest.PROPERTIES_TO_BE_CONSISTENT, list),\
               "PROPERTIES_TO_BE_CONSISTENT is of type: "\
               + str(type(totest.PROPERTIES_TO_BE_CONSISTENT))

    def test_instance_ALL_PROPERTIES(self):
        assert isinstance(totest.ALL_PROPERTIES, list),\
               "ALL_PROPERTIES is of type: " + str(type(totest.ALL_PROPERTIES))

    def test_instance_join_grids(self):
        assert isroutine(totest.join_grids)


class TestValues:
    # check that the values fit
    def test_value_HDF5_MEMBER_SIZE(self):
        assert totest.HDF5_MEMBER_SIZE == 2147483647

    def test_value_H5_UNICODE_DTYPE(self):
        assert totest.H5_UNICODE_DTYPE == h5py.string_dtype()

    def test_value_H5_REC_STR_DTYPE(self):
        assert "S" in totest.H5_REC_STR_DTYPE

    def test_value_VALID_KEYS(self):
        for v in ["binary_history", "history1", "history2", "final_profile1",\
                  "final_profile2", "initial_values", "final_values"]:
            assert v in totest.VALID_KEYS

    def test_value_WARN_VALUES(self):
        for v in ["end", "normal", "suppress"]:
            assert v in totest.WARN_VALUES

    def test_value_N_FLAGS(self):
        assert totest.N_FLAGS == 4

    def test_value_N_FLAGS_SINGLE(self):
        assert totest.N_FLAGS_SINGLE == 2

    def test_value_TERMINATION_FLAG_COLUMNS(self):
        for i in range(1,5):
            assert f"termination_flag_{i}" in totest.TERMINATION_FLAG_COLUMNS

    def test_value_TERMINATION_FLAG_COLUMNS_SINGLE(self):
        for i in [1, 3]:
            assert f"termination_flag_{i}" in\
                   totest.TERMINATION_FLAG_COLUMNS_SINGLE

    def test_value_DEFAULT_BINARY_HISTORY_COLS(self):
        for v in ["model_number", "age", "star_1_mass", "star_2_mass",\
                  "period_days", "binary_separation", "lg_system_mdot_1",\
                  "lg_system_mdot_2", "lg_wind_mdot_1", "lg_wind_mdot_2",\
                  "lg_mstar_dot_1", "lg_mstar_dot_2", "lg_mtransfer_rate",\
                  "xfer_fraction", "rl_relative_overflow_1",\
                  "rl_relative_overflow_2", "trap_radius", "acc_radius",\
                  "t_sync_rad_1", "t_sync_conv_1", "t_sync_rad_2",\
                  "t_sync_conv_2"]:
            assert v in totest.DEFAULT_BINARY_HISTORY_COLS

    def test_value_DEFAULT_STAR_HISTORY_COLS(self):
        for v in ["he_core_mass", "c_core_mass", "o_core_mass",\
                  "he_core_radius", "c_core_radius", "o_core_radius",\
                  "center_h1", "center_he4", "center_c12", "center_n14",\
                  "center_o16", "surface_h1", "surface_he4", "surface_c12",\
                  "surface_n14", "surface_o16", "c12_c12", "center_gamma",\
                  "avg_c_in_c_core", "surf_avg_omega",\
                  "surf_avg_omega_div_omega_crit", "log_LH", "log_LHe",\
                  "log_LZ", "log_Lnuc", "log_Teff", "log_L", "log_R",\
                  "log_center_T", "log_center_Rho", "total_moment_of_inertia",\
                  "spin_parameter", "log_total_angular_momentum",\
                  "conv_env_top_mass", "conv_env_bot_mass",\
                  "conv_env_top_radius", "conv_env_bot_radius",\
                  "conv_env_turnover_time_g", "conv_env_turnover_time_l_b",\
                  "conv_env_turnover_time_l_t", "envelope_binding_energy",\
                  "mass_conv_reg_fortides", "thickness_conv_reg_fortides",\
                  "radius_conv_reg_fortides", "lambda_CE_1cent",\
                  "lambda_CE_10cent", "lambda_CE_30cent", "co_core_mass",\
                  "co_core_radius", "lambda_CE_pure_He_star_10cent",\
                  "log_L_div_Ledd"]:
            assert v in totest.DEFAULT_STAR_HISTORY_COLS

    def test_value_DEFAULT_SINGLE_HISTORY_COLS(self):
        for v in ["model_number", "star_age", "star_mass"]:
            assert v in totest.DEFAULT_SINGLE_HISTORY_COLS
        for v in totest.DEFAULT_STAR_HISTORY_COLS:
            assert v in totest.DEFAULT_SINGLE_HISTORY_COLS

    def test_value_DEFAULT_EEP_HISTORY_COLS(self):
        for v in ["star_age", "star_mass"]:
            assert v in totest.DEFAULT_EEP_HISTORY_COLS
        for v in totest.DEFAULT_STAR_HISTORY_COLS:
            assert v in totest.DEFAULT_EEP_HISTORY_COLS

    def test_value_DEFAULT_PROFILE_COLS(self):
        for v in ["radius", "mass", "logRho", "omega", "energy",\
                  "x_mass_fraction_H", "y_mass_fraction_He",\
                  "z_mass_fraction_metals", "neutral_fraction_H",\
                  "neutral_fraction_He", "avg_charge_He"]:
            assert v in totest.DEFAULT_PROFILE_COLS

    def test_value_EXTRA_STAR_COLS_AT_HE_DEPLETION(self):
        for v in ["avg_c_in_c_core", "co_core_mass"]:
            assert v in totest.EXTRA_STAR_COLS_AT_HE_DEPLETION

    def test_value_DEFAULT_HISTORY_DS_EXCLUDE(self):
        for v in ["model_number", "age", "star_age"]:
            assert v in totest.DEFAULT_HISTORY_DS_EXCLUDE

    def test_value_DEFAULT_PROFILE_DS_EXCLUDE(self):
        for v in ["mass", "star_mass"]:
            assert v in totest.DEFAULT_PROFILE_DS_EXCLUDE

    def test_value_EXTRA_COLS_DS_EXCLUDE(self):
        for v in ["acc_radius", "age", "avg_charge_He", "c12_c12",\
                  "c_core_radius", "center_c12", "center_gamma", "center_n14",\
                  "center_o16", "co_core_radius", "conv_env_bot_mass",\
                  "conv_env_bot_radius", "conv_env_top_mass",\
                  "conv_env_top_radius", "conv_env_turnover_time_g",\
                  "conv_env_turnover_time_l_b", "conv_env_turnover_time_l_t",\
                  "energy", "envelope_binding_energy", "he_core_radius",\
                  "lambda_CE_10cent", "lambda_CE_1cent", "lambda_CE_30cent",\
                  "lambda_CE_pure_He_star_10cent", "lg_mstar_dot_1",\
                  "lg_mstar_dot_2", "lg_wind_mdot_1", "lg_wind_mdot_2",\
                  "log_LH", "log_LHe", "log_LZ", "log_Lnuc", "mass",\
                  "mass_conv_reg_fortides", "model_number",\
                  "neutral_fraction_H", "neutral_fraction_He",\
                  "o_core_radius", "radius_conv_reg_fortides",\
                  "rl_relative_overflow_1", "rl_relative_overflow_2",\
                  "spin_parameter", "star_age", "star_mass", "surf_avg_omega",\
                  "surf_avg_omega_div_omega_crit", "surface_c12",\
                  "surface_h1", "surface_he4", "surface_n14", "surface_o16",\
                  "t_sync_conv_1", "t_sync_conv_2", "t_sync_rad_1",\
                  "t_sync_rad_2", "thickness_conv_reg_fortides",\
                  "trap_radius", "x_mass_fraction_H", "xfer_fraction",\
                  "y_mass_fraction_He", "z_mass_fraction_metals",\
                  "log_L_div_Ledd"]:
            assert v in totest.EXTRA_COLS_DS_EXCLUDE

    def test_value_GRIDPROPERTIES(self):
        for v in ["description", "format", "compression",\
                  "star1_history_saved_columns",\
                  "star2_history_saved_columns",\
                  "binary_history_saved_columns",\
                  "star1_profile_saved_columns",\
                  "star2_profile_saved_columns", "eep"]:
            assert v in totest.GRIDPROPERTIES
            assert isinstance(totest.GRIDPROPERTIES[v], (str, type(None)))
        for v in ["max_number_of_runs", "history_DS_error",\
                  "profile_DS_error", "profile_DS_interval"]:
            assert v in totest.GRIDPROPERTIES
            assert isinstance(totest.GRIDPROPERTIES[v], (int, float,\
                                                         type(None)))
        for v in ["history_DS_exclude", "profile_DS_exclude",\
                  "initial_value_columns", "final_value_columns"]:
            assert v in totest.GRIDPROPERTIES
            assert isinstance(totest.GRIDPROPERTIES[v], (list, type(None)))
        for v in ["start_at_RLO", "stop_before_carbon_depletion", "binary",\
                  "initial_RLO_fix", "He_core_fix", "accept_missing_profile"]:
            assert v in totest.GRIDPROPERTIES
            assert isinstance(totest.GRIDPROPERTIES[v], (bool, type(None)))

    def test_value_PROPERTIES_ALLOWED(self):
        for v in ["format"]:
            assert v in totest.PROPERTIES_ALLOWED
            assert isinstance(totest.PROPERTIES_ALLOWED[v], str)

    def test_value_PROPERTIES_TO_BE_SET(self):
        for v in ["compression", "description"]:
            assert v in totest.PROPERTIES_TO_BE_SET

    def test_value_PROPERTIES_TO_BE_NONE(self):
        for v in ["max_number_of_runs", "star1_history_saved_columns",\
                  "star2_history_saved_columns",\
                  "binary_history_saved_columns",\
                  "star1_profile_saved_columns",\
                  "star2_profile_saved_columns",\
                  "initial_value_columns", "final_value_columns"]:
            assert v in totest.PROPERTIES_TO_BE_NONE
            assert totest.PROPERTIES_TO_BE_NONE[v] is None

    def test_value_PROPERTIES_TO_BE_CONSISTENT(self):
        for v in ["binary", "eep", "start_at_RLO",\
                  "stop_before_carbon_depletion", "initial_RLO_fix",\
                  "He_core_fix", "accept_missing_profile", "history_DS_error",\
                  "history_DS_exclude", "profile_DS_error",\
                  "profile_DS_exclude", "profile_DS_interval"]:
            assert v in totest.PROPERTIES_TO_BE_CONSISTENT

    def test_value_ALL_PROPERTIES(self):
        for v in totest.PROPERTIES_ALLOWED.keys():
            assert v in totest.ALL_PROPERTIES
        for v in totest.PROPERTIES_TO_BE_CONSISTENT:
            assert v in totest.ALL_PROPERTIES
        for v in totest.PROPERTIES_TO_BE_NONE.keys():
            assert v in totest.ALL_PROPERTIES
        for v in totest.PROPERTIES_TO_BE_SET:
            assert v in totest.ALL_PROPERTIES


class TestFunctions:
    @fixture
    def ConfigFile(self):
        # a ConfigFile object for testing
        configfile = totest.ConfigFile()
        configfile.entries["binary"] = False
        configfile.entries["history_DS_error"] = None
        configfile.entries["history_DS_exclude"] = None
        configfile.entries["profile_DS_error"] = None
        configfile.entries["profile_DS_exclude"] = None
        configfile.entries["profile_DS_interval"] = None
        return configfile

    @fixture
    def star_history(self):
        # a temporary star history for testing
        return np.array([(1.0, 0.2), (1.0e+2, 0.9), (1.0e+3, 0.2)],\
                        dtype=[('star_age', '<f8'), ('center_he4', '<f8')])

    @fixture
    def binary_history(self):
        # a temporary binary history for testing
        return np.array([(1.0, 1.0), (1.1, 1.0e+2), (1.2, 1.0e+3)],\
                        dtype=[('period_days', '<f8'), ('age', '<f8')])

    @fixture
    def profile(self):
        # a temporary profile for testing
        return np.array([(2.0, 1.0e+3), (1.1, 1.0e+2), (0.1, 1.0)],\
                        dtype=[('mass', '<f8'), ('radius', '<f8')])

    @fixture
    def no_path(self, tmp_path):
        # a path which does not exist for testing
        return os.path.join(tmp_path, "does_not_exist.test")

    @fixture
    def h5_out_path(self, tmp_path):
        # a path to write to
        return os.path.join(tmp_path, "out.h5")

    @fixture
    def grid1_path(self, tmp_path, binary_history, star_history, profile):
        # a path to a psygrid file for testing
        path = os.path.join(tmp_path, "grid1.h5")
        PSyGrid = totest.PSyGrid()
        PSyGrid.filepath = path
        PSyGrid.generate_config()
        mesa_dir1 = os.path.join(tmp_path, "m1_1.0_m2_1.0_initial_period_in"\
                                           +"_days_0.0_initial_z_0.01_idx_0")
        mesa_dir2 = os.path.join(tmp_path, "m1_1.0_m2_1.0_initial_period_in"\
                                           +"_days_{}_initial_z_0.01_idx_1".\
                                           format(\
                                            binary_history['period_days'][0]))
        with h5py.File(path, "w") as hdf5_file:
            hdf5_file.create_group("/grid/run0/")
            hdf5_file.create_dataset("/grid/run1/binary_history",\
                                     data=binary_history)
            hdf5_file.create_dataset("/grid/run1/history1", data=star_history)
            hdf5_file.create_dataset("/grid/run1/history2", data=star_history)
            hdf5_file.create_dataset("/grid/run1/final_profile1", data=profile)
            hdf5_file.create_dataset("/grid/run1/final_profile2", data=profile)
            ini_val = np.array([(np.nan), (binary_history['period_days'][0])],\
                               dtype=[('period_days', '<f8')])
            hdf5_file.create_dataset("/grid/initial_values", data=ini_val)
            fin_val = np.array([(np.nan, "TF1"),\
                                (binary_history['period_days'][-1], "ignored_no_RLO")],\
                               dtype=[('period_days', '<f8'),\
                                      ('termination_flag_1',\
                                       totest.H5_UNICODE_DTYPE)])
            hdf5_file.create_dataset("/grid/final_values", data=fin_val)
            hdf5_file.attrs["config"] =\
             totest.json.dumps(str(dict(PSyGrid.config)))
            rel_paths = np.array([(mesa_dir1), (mesa_dir2)],\
                                 dtype=totest.H5_UNICODE_DTYPE)
            hdf5_file.create_dataset("relative_file_paths", data=rel_paths)
        return path

    @fixture
    def grid2_path(self, tmp_path, binary_history, star_history, profile):
        # a path to a psygrid file for testing
        path = os.path.join(tmp_path, "grid2.h5")
        PSyGrid = totest.PSyGrid()
        PSyGrid.filepath = path
        PSyGrid.generate_config()
        mesa_dir1 = os.path.join(tmp_path, "m1_1.0_m2_1.0_initial_period_in"\
                                           +"_days_0.0_initial_z_0.01_idx_0")
        mesa_dir2 = os.path.join(tmp_path, "m1_1.0_m2_1.0_initial_period_in"\
                                           +"_days_{}_initial_z_0.01_idx_1".\
                                           format(\
                                            binary_history['period_days'][0]))
        with h5py.File(path, "w") as hdf5_file:
            hdf5_file.create_group("/grid/run0/")
            hdf5_file.create_dataset("/grid/run1/binary_history",\
                                     data=binary_history)
            hdf5_file.create_dataset("/grid/run1/history1", data=star_history)
            hdf5_file.create_dataset("/grid/run1/history2", data=star_history)
            hdf5_file.create_dataset("/grid/run1/final_profile1", data=profile)
            hdf5_file.create_dataset("/grid/run1/final_profile2", data=profile)
            ini_val = np.array([(np.nan), (binary_history['period_days'][0])],\
                               dtype=[('period_days', '<f8')])
            hdf5_file.create_dataset("/grid/initial_values", data=ini_val)
            fin_val = np.array([(np.nan, "TF1"),\
                                (binary_history['period_days'][-1], "TF1")],\
                               dtype=[('period_days', '<f8'),\
                                      ('termination_flag_1',\
                                       totest.H5_UNICODE_DTYPE)])
            hdf5_file.create_dataset("/grid/final_values", data=fin_val)
            hdf5_file.attrs["config"] =\
             totest.json.dumps(str(dict(PSyGrid.config)))
            rel_paths = np.array([(mesa_dir1), (mesa_dir2)],\
                                 dtype=totest.H5_UNICODE_DTYPE)
            hdf5_file.create_dataset("relative_file_paths", data=rel_paths)
        return path

    @fixture
    def grid_path_bad_ini(self, tmp_path, binary_history, star_history,\
                          profile):
        # a path to a psygrid file for testing
        path = os.path.join(tmp_path, "grid_bad_ini.h5")
        PSyGrid = totest.PSyGrid()
        PSyGrid.filepath = path
        PSyGrid.generate_config()
        mesa_dir1 = os.path.join(tmp_path, "m1_1.0_m2_1.0_initial_period_in"\
                                           +"_days_0.0_initial_z_0.01")
        mesa_dir2 = os.path.join(tmp_path, "m1_1.0_m2_1.0_initial_period_in"\
                                           +"_days_{}_initial_z_0.01".format(\
                                            binary_history['period_days'][0]))
        with h5py.File(path, "w") as hdf5_file:
            hdf5_file.create_group("/grid/run0/")
            hdf5_file.create_dataset("/grid/run1/binary_history",\
                                     data=binary_history)
            hdf5_file.create_dataset("/grid/run1/history1", data=star_history)
            hdf5_file.create_dataset("/grid/run1/history2", data=star_history)
            hdf5_file.create_dataset("/grid/run1/final_profile1", data=profile)
            hdf5_file.create_dataset("/grid/run1/final_profile2", data=profile)
            ini_val = np.array([(np.nan), (binary_history['period_days'][0])],\
                               dtype=[('bad', '<f8')])
            hdf5_file.create_dataset("/grid/initial_values", data=ini_val)
            fin_val = np.array([(np.nan, "TF1"),\
                                (binary_history['period_days'][-1], "TF1")],\
                               dtype=[('period_days', '<f8'),\
                                      ('termination_flag_1',\
                                       totest.H5_UNICODE_DTYPE)])
            hdf5_file.create_dataset("/grid/final_values", data=fin_val)
            hdf5_file.attrs["config"] =\
             totest.json.dumps(str(dict(PSyGrid.config)))
            rel_paths = np.array([(mesa_dir1), (mesa_dir2)],\
                                 dtype=totest.H5_UNICODE_DTYPE)
            hdf5_file.create_dataset("relative_file_paths", data=rel_paths)
        return path

    @fixture
    def grid_path_bad_fin(self, tmp_path, binary_history, star_history,\
                          profile):
        # a path to a psygrid file for testing
        path = os.path.join(tmp_path, "grid_bad_fin.h5")
        PSyGrid = totest.PSyGrid()
        PSyGrid.filepath = path
        PSyGrid.generate_config()
        mesa_dir1 = os.path.join(tmp_path, "m1_1.0_m2_1.0_initial_period_in"\
                                           +"_days_0.0_initial_z_0.01")
        mesa_dir2 = os.path.join(tmp_path, "m1_1.0_m2_1.0_initial_period_in"\
                                           +"_days_{}_initial_z_0.01".format(\
                                            binary_history['period_days'][0]))
        with h5py.File(path, "w") as hdf5_file:
            hdf5_file.create_group("/grid/run0/")
            hdf5_file.create_dataset("/grid/run1/binary_history",\
                                     data=binary_history)
            hdf5_file.create_dataset("/grid/run1/history1", data=star_history)
            hdf5_file.create_dataset("/grid/run1/history2", data=star_history)
            hdf5_file.create_dataset("/grid/run1/final_profile1", data=profile)
            hdf5_file.create_dataset("/grid/run1/final_profile2", data=profile)
            ini_val = np.array([(np.nan), (binary_history['period_days'][0])],\
                               dtype=[('period_days', '<f8')])
            hdf5_file.create_dataset("/grid/initial_values", data=ini_val)
            fin_val = np.array([(np.nan, "TF1"),\
                                (binary_history['period_days'][-1], "TF1")],\
                               dtype=[('period_days', '<f8'),\
                                      ('termination_flag_bad',\
                                       totest.H5_UNICODE_DTYPE)])
            hdf5_file.create_dataset("/grid/final_values", data=fin_val)
            hdf5_file.attrs["config"] =\
             totest.json.dumps(str(dict(PSyGrid.config)))
            rel_paths = np.array([(mesa_dir1), (mesa_dir2)],\
                                 dtype=totest.H5_UNICODE_DTYPE)
            hdf5_file.create_dataset("relative_file_paths", data=rel_paths)
        return path

    @fixture
    def grid_path_start_RLO(self, tmp_path, binary_history, star_history,\
                            profile):
        # a path to a psygrid file for testing
        path = os.path.join(tmp_path, "grid_start_RLO.h5")
        PSyGrid = totest.PSyGrid()
        PSyGrid.filepath = path
        PSyGrid.generate_config()
        PSyGrid.config["start_at_RLO"] = True
        mesa_dir1 = os.path.join(tmp_path, "m1_1.0_m2_1.0_initial_period_in"\
                                           +"_days_0.0_initial_z_0.01_idx_0")
        mesa_dir2 = os.path.join(tmp_path, "m1_1.0_m2_1.0_initial_period_in"\
                                           +"_days_{}_initial_z_0.01_idx_1".\
                                           format(\
                                            binary_history['period_days'][0]))
        with h5py.File(path, "w") as hdf5_file:
            hdf5_file.create_group("/grid/run0/")
            hdf5_file.create_dataset("/grid/run1/binary_history",\
                                     data=binary_history)
            hdf5_file.create_dataset("/grid/run1/history1", data=star_history)
            hdf5_file.create_dataset("/grid/run1/history2", data=star_history)
            hdf5_file.create_dataset("/grid/run1/final_profile1", data=profile)
            hdf5_file.create_dataset("/grid/run1/final_profile2", data=profile)
            ini_val = np.array([(np.nan), (binary_history['period_days'][0])],\
                               dtype=[('period_days', '<f8')])
            hdf5_file.create_dataset("/grid/initial_values", data=ini_val)
            fin_val = np.array([(np.nan, "TF1"),\
                                (binary_history['period_days'][-1], "TF1")],\
                               dtype=[('period_days', '<f8'),\
                                      ('termination_flag_1',\
                                       totest.H5_UNICODE_DTYPE)])
            hdf5_file.create_dataset("/grid/final_values", data=fin_val)
            hdf5_file.attrs["config"] =\
             totest.json.dumps(str(dict(PSyGrid.config)))
            rel_paths = np.array([(mesa_dir1), (mesa_dir2)],\
                                 dtype=totest.H5_UNICODE_DTYPE)
            hdf5_file.create_dataset("relative_file_paths", data=rel_paths)
        return path

    @fixture
    def grid_path_initial_RLO(self, tmp_path, binary_history, star_history,\
                              profile):
        # a path to a psygrid file for testing
        path = os.path.join(tmp_path, "grid_initial_RLO.h5")
        PSyGrid = totest.PSyGrid()
        PSyGrid.filepath = path
        PSyGrid.generate_config()
        PSyGrid.config["description"] = "UnitTest1"
        PSyGrid.config["initial_RLO_fix"] = True
        ini_val = np.array([(np.nan, np.nan, np.nan),\
                            (1.0, 1.0, binary_history['period_days'][0]),\
                            (1.0, 1.0, 0.5), (1.0, 1.0, 0.1),\
                            (3.0, 1.0, 0.1), (4.0, 1.0, 0.1)],\
                           dtype=[('star_1_mass', '<f8'),\
                                  ('star_2_mass', '<f8'),\
                                  ('period_days', '<f8')])
        mesa_dir = [os.path.join(tmp_path, "m1_1.0_m2_1.0_initial_period_in_"\
                                           +"days_0.0_initial_z_0.01_idx_0")]
        for i in range(1,6):
            m1 = ini_val['star_1_mass'][i]
            m2 = ini_val['star_2_mass'][i]
            p = ini_val['period_days'][i]
            mesa_dir += [os.path.join(tmp_path, f"m1_{m1}_m2_{m2}_initial_"\
                                                +f"period_in_days_{p}_"\
                                                +f"initial_z_0.01_idx_{i}")]
        with h5py.File(path, "w") as hdf5_file:
            hdf5_file.create_group("/grid/run0/")
            for r in range(1,6):
                hdf5_file.create_dataset(f"/grid/run{r}/binary_history",\
                                         data=binary_history)
                hdf5_file.create_dataset(f"/grid/run{r}/history1",\
                                         data=star_history)
                hdf5_file.create_dataset(f"/grid/run{r}/history2",\
                                         data=star_history)
                hdf5_file.create_dataset(f"/grid/run{r}/final_profile1",\
                                         data=profile)
                hdf5_file.create_dataset(f"/grid/run{r}/final_profile2",\
                                         data=profile)
            hdf5_file.create_dataset("/grid/initial_values", data=ini_val)
            fin_val = np.array([(np.nan, np.nan, np.nan, "TF1", "TF2", "TF3",\
                                 "TF4", "IC"),\
                                (0.9, 0.9, binary_history['period_days'][-1],\
                                 "TF1", "TF2", "TF3", "TF4", "IC"),\
                                (1.0, 1.0, 0.5, "Terminate because of "\
                                                +"overflowing initial model",\
                                 "iRLOTF2", "iRLOTF3", "iRLOTF4", "IC"),\
                                (1.0, 1.0, 0.1, "forced_initial_RLO",\
                                 "fRLOTF2", "fRLOTF3", "fRLOTF4", "IC"),\
                                (3.0, 1.0, 0.1, "TF1", "TF2", "TF3", "TF4",\
                                 "IC"),\
                                (4.0, 1.0, 0.1, "TF1", "TF2", "TF3", "TF4",\
                                 "IC")],\
                               dtype=[('star_1_mass', '<f8'),\
                                      ('star_2_mass', '<f8'),\
                                      ('period_days', '<f8'),\
                                      ('termination_flag_1',\
                                       totest.H5_UNICODE_DTYPE),\
                                      ('termination_flag_2',\
                                       totest.H5_UNICODE_DTYPE),\
                                      ('termination_flag_3',\
                                       totest.H5_UNICODE_DTYPE),\
                                      ('termination_flag_4',\
                                       totest.H5_UNICODE_DTYPE),\
                                      ('interpolation_class',\
                                       totest.H5_UNICODE_DTYPE)])
            hdf5_file.create_dataset("/grid/final_values", data=fin_val)
            hdf5_file.attrs["config"] =\
             totest.json.dumps(str(dict(PSyGrid.config)))
            rel_paths = np.array([(d) for d in mesa_dir],\
                                 dtype=totest.H5_UNICODE_DTYPE)
            hdf5_file.create_dataset("relative_file_paths", data=rel_paths)
        return path

    # test functions
    def test_downsample_history(self, ConfigFile, star_history,\
                                binary_history):
        # missing argument
        with raises(TypeError, match="missing 4 required positional "\
                                     +"arguments: 'bh', 'h1', 'h2', and "\
                                     +"'params'"):
            totest.downsample_history()
        # bad input
        with raises(TypeError, match="'NoneType' object is not subscriptable"):
            totest.downsample_history(None, None, None, None)
        # examples: failed
        for bh in [None, binary_history]:
            for h1 in [None, star_history]:
                for h2 in [None, star_history]:
                    assert totest.downsample_history(bh, h1, h2, ConfigFile)\
                           == (bh, h1, h2)
        ConfigFile.entries["history_DS_error"] = 0.1
        assert totest.downsample_history(None, None, None, ConfigFile) ==\
               (None, None, None)
        # bad input
        tests = [(bh, h1, h2) for bh in [None,\
                                         binary_history[['period_days']]]\
                              for h1 in [None, star_history[['center_he4']]]\
                              for h2 in [None, star_history]]
        tests.remove((None, None, None))
        for (bh, h1, h2) in tests:
            with raises(TypeError, match="argument of type 'NoneType' is not "\
                                         +"iterable"):
                totest.downsample_history(bh, h1, h2, ConfigFile)
        ConfigFile.entries["history_DS_exclude"] = ['star_age']
        for (bh, h1, h2) in tests:
            if h1 is None:
                with raises(TypeError, match="'NoneType' object is not "\
                                             +"subscriptable"):
                    totest.downsample_history(bh, h1, h2, ConfigFile)
            else:
                with raises(ValueError, match="no field of name star_age"):
                    totest.downsample_history(bh, h1, h2, ConfigFile)
        ConfigFile.entries["binary"] = True
        for (bh, h1, h2) in tests:
            if bh is None:
                with raises(TypeError, match="'NoneType' object is not "\
                                             +"subscriptable"):
                    totest.downsample_history(bh, h1, h2, ConfigFile)
            else:
                with raises(ValueError, match="no field of name age"):
                    totest.downsample_history(bh, h1, h2, ConfigFile)
        with raises(IndexError, match="Unequal numbers of rows in histories."):
            totest.downsample_history(binary_history, star_history,\
                                      star_history[1:], ConfigFile)
        # examples:
        bh, h1, h2 = totest.downsample_history(binary_history, star_history,\
                                               star_history, ConfigFile)
        assert np.array_equal(bh, binary_history)
        assert np.array_equal(h1, star_history)
        assert np.array_equal(h2, star_history)

    def test_downsample_profile(self, ConfigFile, profile):
        # missing argument
        with raises(TypeError, match="missing 2 required positional "\
                                     +"arguments: 'profile' and 'params'"):
            totest.downsample_profile()
        # bad input
        with raises(TypeError, match="'NoneType' object is not subscriptable"):
            totest.downsample_profile(None, None)
        # examples: failed
        assert np.array_equal(totest.downsample_profile(profile, ConfigFile),\
                              profile)
        ConfigFile.entries["profile_DS_error"] = 0.1
        assert totest.downsample_profile(None, ConfigFile) is None
        # bad input
        with raises(TypeError, match="argument of type 'NoneType' is not "\
                                     +"iterable"):
            totest.downsample_profile(profile, ConfigFile)
        # examples
        ConfigFile.entries["profile_DS_exclude"] = ['mass']
        assert np.array_equal(totest.downsample_profile(profile, ConfigFile),\
                              profile)

    def test_join_grids(self, monkeypatch, no_path, h5_out_path, grid1_path,\
                        grid2_path, grid_path_bad_ini, grid_path_bad_fin,\
                        grid_path_start_RLO, grid_path_initial_RLO):
        def check_h5_content(path, required_data={}, required_attrs={}):
            """Function to check an hdf5 file

            Parameters
            ----------
            path : str
                Location of the hdf5 file
            required_data : dict (optional)
                A dictionary where the keys are references to objects/datasets
                in the hdf5 file and the items are the ndarrays to be checked
                against (if the item is None, only the existence of the object
                is checked)
            required_attrs : dict (optional)
                A dictionary where the keys are references to attributes in the
                hdf5 file and the items are lists of strings to be part of this
                attribute
            """
            # test that file exists
            assert os.path.isfile(path)
            with h5py.File(path, "r") as hdf5_file:
                # test attribute content
                for s, vs in required_attrs.items():
                    for v in vs:
                        assert v in hdf5_file.attrs[s]
                # test data content
                for s, v in required_data.items():
                    if v is None:
                        assert s in hdf5_file.keys()
                    else:
                        assert np.array_equal(hdf5_file[s][()], v)
        def mock_get_detected_initial_RLO(grid):
            # mocked list of initial RLO systems
            if (hasattr(grid, "config") and ("description" in grid.config)):
                if grid.config["description"] == "UnitTest1":
                    return [{"star_1_mass": 1.0, "star_2_mass": 1.0,\
                             "period_days": 1.0, "termination_flag_3": "rTF3",\
                             "termination_flag_4": "rTF4"},\
                            {"star_1_mass": 1.0, "star_2_mass": 1.0,\
                             "period_days": 0.1, "termination_flag_3": "rTF3",\
                             "termination_flag_4": "rTF4"},\
                            {"star_1_mass": 1.0, "star_2_mass": 1.0,\
                             "period_days": 2.0, "termination_flag_3": "rTF3",\
                             "termination_flag_4": "rTF4"},\
                            {"star_1_mass": 2.0, "star_2_mass": 1.0,\
                             "period_days": 1.0, "termination_flag_3": "rTF3",\
                             "termination_flag_4": "rTF4"}]
                else:
                    return []
            else:
                return []
        def mock_get_nearest_known_initial_RLO(mass1, mass2,\
                                               known_initial_RLO):
            # mocked nearest initial RLO system
            if mass1 == 3.0:
                return {"star_1_mass": mass1, "star_2_mass": mass2,\
                        "period_days": 2.0, "termination_flag_3": "rTF3"}
            elif mass1 == 4.0:
                return {"star_1_mass": mass1, "star_2_mass": mass2,\
                        "period_days": 2.0, "termination_flag_4": "rTF4"}
            else:
                return {"star_1_mass": mass1, "star_2_mass": mass2,\
                        "period_days": 2.0}
        # missing argument
        with raises(TypeError, match="missing 2 required positional "\
                                     +"arguments: 'input_paths' and "\
                                     +"'output_path'"):
            totest.join_grids()
        # bad input
        with raises(FileNotFoundError):
            totest.join_grids([no_path], h5_out_path)
        # examples: one grid
        totest.join_grids([grid2_path], h5_out_path, verbose=False)
        rd = {}
        for group in ['/', '/grid/', '/grid/run0/', '/grid/run1/']:
            for key in h5py.File(grid2_path, "r")[group].keys():
                rd[group+key] = None
        ra = {"config": ["'description': 'joined'", "'compression': 'gzip9'"]\
                        +[f"'{p}': None" for p in\
                          totest.PROPERTIES_TO_BE_NONE]}
        check_h5_content(h5_out_path, required_data=rd, required_attrs=ra)
        # examples: two grids
        totest.join_grids([grid1_path, grid2_path], h5_out_path, verbose=False)
        check_h5_content(h5_out_path, required_data=rd, required_attrs=ra)
        # examples: three grids
        totest.join_grids([grid2_path, grid1_path, grid2_path], h5_out_path, verbose=False)
        check_h5_content(h5_out_path, required_data=rd, required_attrs=ra)
        # bad input: consistent values do not match
        with h5py.File(grid2_path, "r") as hdf5_file:
            config2 = hdf5_file.attrs["config"]
        for p in totest.PROPERTIES_TO_BE_CONSISTENT:
            sep = f"'{p}': "
            parts = config2.split(sep, maxsplit=1)
            if parts[1][0] == "[": # check for lists
                parts = [parts[0]+sep] + parts[1].split("],", maxsplit=1)
            elif parts[1][0] == "{": # check for dicts
                parts = [parts[0]+sep] + parts[1].split("},", maxsplit=1)
            else:
                parts = [parts[0]+sep] + parts[1].split(",", maxsplit=1)
            if len(parts) > 2:
                newconfig = parts[0] + "'test'," + parts[2]
            else: # the last entry
                newconfig = parts[0] + "'test'" + parts[1][-2:]
            with h5py.File(grid2_path, "a") as hdf5_file:
                hdf5_file.attrs["config"] = newconfig
            with raises(ValueError, match=f"Inconsistent value for `{p}`: "\
                                          +"test !="):
                totest.join_grids([grid1_path, grid2_path], h5_out_path,\
                                  verbose=False)
        with h5py.File(grid2_path, "a") as hdf5_file: # restore config
            hdf5_file.attrs["config"] = config2
        # bad input: initial values have different type
        with raises(TypeError, match="Initial values dtype's do not match."):
            totest.join_grids([grid1_path, grid_path_bad_ini], h5_out_path,\
                              verbose=False)
        # bad input: final values have different type
        with raises(TypeError, match="Final values dtype's do not match."):
            totest.join_grids([grid1_path, grid_path_bad_fin], h5_out_path,\
                              verbose=False)
        # bad input: disallowed value (reuse config2)
        for p,v in totest.PROPERTIES_ALLOWED.items():
            sep = f"'{p}': "
            parts = config2.split(sep, maxsplit=1)
            if parts[1][0] == "[": # check for lists
                parts = [parts[0]+sep] + parts[1].split("],", maxsplit=1)
            elif parts[1][0] == "{": # check for dicts
                parts = [parts[0]+sep] + parts[1].split("},", maxsplit=1)
            else:
                parts = [parts[0]+sep] + parts[1].split(",", maxsplit=1)
            if len(parts) > 2:
                newconfig = parts[0] + "'test'," + parts[2]
            else: # the last entry
                newconfig = parts[0] + "'test'" + parts[1][-2:]
            with h5py.File(grid2_path, "a") as hdf5_file:
                hdf5_file.attrs["config"] = newconfig
            with raises(ValueError, match=f"Only grid with `{p}`={v} can be "\
                                          +"joined."):
                totest.join_grids([grid2_path], h5_out_path, verbose=False)
        with h5py.File(grid2_path, "a") as hdf5_file: # restore config
            hdf5_file.attrs["config"] = config2
        # examples: RLO grid
        totest.join_grids([grid_path_start_RLO], h5_out_path, verbose=False)
        check_h5_content(h5_out_path, required_data=rd, required_attrs=ra)
        # examples: initial RLO
        monkeypatch.setattr(totest, "get_detected_initial_RLO",\
                            mock_get_detected_initial_RLO)
        monkeypatch.setattr(totest, "get_nearest_known_initial_RLO",\
                            mock_get_nearest_known_initial_RLO)
        totest.join_grids([grid_path_initial_RLO], h5_out_path, verbose=False)
        check_h5_content(h5_out_path, required_data=rd, required_attrs=ra)
        with h5py.File(grid_path_initial_RLO, "r") as hdf5_file:
            ini_val = hdf5_file['/grid/initial_values'][()]
            fin_val = hdf5_file['/grid/final_values'][()]
        with h5py.File(h5_out_path, "r") as hdf5_file:
            new_ini_val = hdf5_file['/grid/initial_values'][()]
            new_fin_val = hdf5_file['/grid/final_values'][()]
        for col in ['star_1_mass', 'star_2_mass', 'period_days']:
            assert np.array_equal(ini_val[col], new_ini_val[col],\
                                  equal_nan=True)
            assert np.array_equal(fin_val[col], new_fin_val[col],\
                                  equal_nan=True)
        for r in [1, 4, 5]: # done initial RLO replacements
            fin_val['termination_flag_1'][r] = b'forced_initial_RLO'
            fin_val['termination_flag_2'][r] = b'forced_initial_RLO'
            fin_val['interpolation_class'][r] = b'initial_MT'
        fin_val['termination_flag_3'][4] = b'rTF3'
        fin_val['termination_flag_4'][5] = b'rTF4'
        for col in ['termination_flag_1', 'termination_flag_2',\
                    'termination_flag_3', 'termination_flag_4',\
                    'interpolation_class']:
            assert np.array_equal(fin_val[col], new_fin_val[col])
        # examples: compressions
        for c in [None, "lzf", "gzip", "gzip9"]:
            totest.join_grids([grid2_path], h5_out_path, compression=c,\
                              verbose=False)
            if c is None:
                ra["config"][1] = "'compression': None"
            else:
                ra["config"][1] = "'compression': '" + c + "'"
            check_h5_content(h5_out_path, required_data=rd, required_attrs=ra)
        for c in ["gzipT", "gzip10"]:
            with raises(ValueError, match="GZIP compression level must be "\
                                          +"0-9."):
                totest.join_grids([grid2_path], h5_out_path, compression=c,\
                                  verbose=False)
        for c in ["zip"]:
            with raises(ValueError, match=f"Unknown compression `{c}`."):
                totest.join_grids([grid2_path], h5_out_path, compression=c,\
                                  verbose=False)


class TestPSyGrid:
    @fixture
    def PSyGrid(self):
        # initialize an instance of the class with defaults
        return totest.PSyGrid()

    @fixture
    def PSyGrid_with_config(self):
        # initialize an instance of the class with default config
        PSyGrid = totest.PSyGrid()
        PSyGrid.generate_config()
        return PSyGrid

    @fixture
    def Empty_dir(self, tmp_path):
        dir_path = os.path.join(tmp_path, "empty")
        os.mkdir(dir_path)
        return dir_path

    @fixture
    def eep_files(self, tmp_path):
        # create eep files and empty directory
        path = os.path.join(tmp_path, "eeps_OK")
        os.mkdir(path)
        for i,ext in enumerate(totest.EEP_FILE_EXTENSIONS):
            with open(os.path.join(path, f"file{i}"+ext), "w") as eep_file:
                eep_file.write(f"Test{i}")
        return path

    @fixture
    def same_eep_files(self, tmp_path):
        # create eep files and empty directory
        path = os.path.join(tmp_path, "eeps_same_name")
        os.mkdir(path)
        for i,ext in enumerate(totest.EEP_FILE_EXTENSIONS):
            with open(os.path.join(path, "file"+ext), "w") as eep_file:
                eep_file.write(f"Test{i}")
        return path

    # test the PSyGrid class
    def test_init(self, PSyGrid, monkeypatch):
        def mock_load(self, filepath=None):
            return filepath
        assert isroutine(PSyGrid.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(PSyGrid, totest.PSyGrid)
        assert PSyGrid.filepath is None
        assert PSyGrid.verbose == False
        monkeypatch.setattr(totest.PSyGrid, "load", mock_load)
        test_PSyGrid = totest.PSyGrid(filepath="./test", verbose=True)
        assert test_PSyGrid.filepath == "./test"
        assert test_PSyGrid.verbose == True

    def test_reset(self, PSyGrid):
        assert isroutine(PSyGrid._reset)
        # this function gets called during __init__
        for i in range(2):
            assert PSyGrid.hdf5 is None
            assert PSyGrid.MESA_dirs == []
            assert PSyGrid.initial_values is None
            assert PSyGrid.final_values is None
            assert isinstance(PSyGrid.config, totest.ConfigFile)
            assert PSyGrid.n_runs == 0
            assert PSyGrid.eeps is None
            PSyGrid.hdf5 = "Test"
            PSyGrid.MESA_dirs = "Test"
            PSyGrid.initial_values = "Test"
            PSyGrid.final_values = "Test"
            PSyGrid.config = "Test"
            PSyGrid.n_runs = "Test"
            PSyGrid.eeps = "Test"
            PSyGrid._reset()

    def test_make_compression_args(self, PSyGrid):
        assert isroutine(PSyGrid._make_compression_args)
        # this function gets called during _reset (hence, during __init__)
        # examples no compression information yet
        del PSyGrid.compression_args
        assert hasattr(PSyGrid, "compression_args") == False
        PSyGrid._make_compression_args()
        assert hasattr(PSyGrid, "compression_args")
        assert PSyGrid.compression_args == {}
        # examples: simple compressions
        for c in [None, "lzf", "gzip"]:
            PSyGrid.config["compression"] = c
            PSyGrid._make_compression_args()
            if c is None:
                assert PSyGrid.compression_args == {}
            else:
                assert PSyGrid.compression_args["compression"] == c
        # examples: gzip with optional level
        for l in range(10):
            PSyGrid.config["compression"] = f"gzip{l}"
            PSyGrid._make_compression_args()
            assert PSyGrid.compression_args["compression"] == "gzip"
            assert PSyGrid.compression_args["compression_opts"] == l
        # bad input
        for c in ["gzipt", "gzip10", "test"]:
            PSyGrid.config["compression"] = c
            if "gzip" in c:
                m = "GZIP compression level must be 0-9."
            else:
                m = f"Unknown compression `{c}`."
            with raises(ValueError, match=m):
                PSyGrid._make_compression_args()

    def test_discover_eeps(self, PSyGrid, Empty_dir, eep_files,\
                           same_eep_files):
        assert isroutine(PSyGrid._discover_eeps)
        # examples: no eeps found
        PSyGrid.eeps = "Test"
        PSyGrid._discover_eeps(Empty_dir)
        assert PSyGrid.eeps is None
        # examples: eeps with each file extension found
        PSyGrid._discover_eeps(eep_files)
        for i in range(len(totest.EEP_FILE_EXTENSIONS)):
            eep_id = f"file{i}"
            assert eep_id in PSyGrid.eeps
            assert totest.EEP_FILE_EXTENSIONS[i] ==\
                   PSyGrid.eeps[eep_id][-len(totest.EEP_FILE_EXTENSIONS[i]):]
        # bad input: eep files with same identifier
        with raises(AssertionError):
            PSyGrid._discover_eeps(same_eep_files)

    def test_say(self, PSyGrid, capsys):
        assert isroutine(PSyGrid._say)
        for v in [True, False]:
            PSyGrid.verbose = v
            PSyGrid._say("Test")
            if v:
                assert capsys.readouterr().out == "Test\n"

    def test_generate_config(self, PSyGrid):
        assert isroutine(PSyGrid.generate_config)
        # bad input
        with raises(KeyError, match="`UnitTest` is not a valid parameter "\
                                    +"name."):
            PSyGrid.generate_config(UnitTest="Test")
        # examples: set each grid property
        for c in totest.GRIDPROPERTIES:
            # compression is special, because it is checked for allowed values
            if c == "compression":
                args = {c: ""}
            else:
                args = {c: "UnitTest"}
            PSyGrid.generate_config(**args)
            for cc in totest.GRIDPROPERTIES:
                if c == cc: # the set value
                    assert PSyGrid.config[cc] == args[c]
                else: # the default value
                    assert PSyGrid.config[cc] == totest.GRIDPROPERTIES[cc]
        # examples: set all grid properties
        args = {}
        for c in totest.GRIDPROPERTIES:
            # compression is special, because it is checked for allowed values
            if c == "compression":
                args[c] = ""
            else:
                args[c] = "UnitTest"
        PSyGrid.generate_config(**args)
        for cc in totest.GRIDPROPERTIES:
            # all have the set value
            assert PSyGrid.config[cc] == args[cc]

    def test_create(self, PSyGrid, Empty_dir, capsys):
        assert isroutine(PSyGrid.create)
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'MESA_grid_path'"):
            PSyGrid.create()
        # bad input
        with raises(ValueError, match="`warn` must be in: "):
            PSyGrid.create("./", warn="UnitTest")
        # bad input
        PSyGrid.filepath = None
        with raises(ValueError, match="The path of the HDF5 file was not "\
                                      +"defined."):
            PSyGrid.create("./")
        # examples
        MESA_PATH = str(Empty_dir)
        with raises(ValueError, match=f"No folders found in {MESA_PATH}"):
            PSyGrid.create(MESA_PATH, psygrid_path=os.path.join(MESA_PATH,\
                                                                "TestGrid.h5"))
        pass

    def test_create_psygrid(self, PSyGrid_with_config, Empty_dir, capsys):
        PSyGrid = PSyGrid_with_config
        assert isroutine(PSyGrid._create_psygrid)
        # missing argument
        with raises(TypeError, match="missing 2 required positional "\
                                     +"arguments: 'MESA_path' and 'hdf5'"):
            PSyGrid._create_psygrid()
        # bad input
        PSyGrid.filepath = None
        with raises(AssertionError):
            PSyGrid._create_psygrid("./", None)
        # examples
        PSyGrid.filepath = os.path.join(Empty_dir, "TestPSyGrid.h5")
        MESA_PATH = str(Empty_dir)
        with raises(ValueError, match=f"No folders found in {MESA_PATH}"):
            PSyGrid._create_psygrid(MESA_PATH, None)
        # examples: EEPs
        PSyGrid.config["eep"] = Empty_dir
        with raises(ValueError, match=f"No folders found in {MESA_PATH}"):
            with warns(ReplaceValueWarning, match="Selected EEPs, switching "\
                                                  +"to single-star grid."):
                PSyGrid._create_psygrid(MESA_PATH, None)
        PSyGrid.config["binary"] = False
        with raises(ValueError, match=f"No folders found in {MESA_PATH}"):
            PSyGrid._create_psygrid(MESA_PATH, None)
        pass

    def test_add_column(self, PSyGrid):
        assert isroutine(PSyGrid.add_column)
        pass

    def test_update_final_values(self, PSyGrid):
        assert isroutine(PSyGrid.update_final_values)
        pass

    def test_reload_hdf5_file(self, PSyGrid):
        assert isroutine(PSyGrid._reload_hdf5_file)
        pass

    def test_load(self, PSyGrid):
        assert isroutine(PSyGrid.load)
        pass

    def test_close(self, PSyGrid):
        assert isroutine(PSyGrid.close)
        pass

    def test_str(self, PSyGrid):
        assert isroutine(PSyGrid.__str__)
        pass

    def test_getitem(self, PSyGrid):
        assert isroutine(PSyGrid.__getitem__)
        pass

    def test_get_pandas_initial_final(self, PSyGrid):
        assert isroutine(PSyGrid.get_pandas_initial_final)
        pass

    def test_len(self, PSyGrid):
        assert isroutine(PSyGrid.__len__)
        pass

    def test_contains(self, PSyGrid):
        assert isroutine(PSyGrid.__contains__)
        pass

    def test_iter(self, PSyGrid):
        assert isroutine(PSyGrid.__iter__)
        pass

    def test_del(self, PSyGrid):
        assert isroutine(PSyGrid.__del__)
        pass

    def test_rerun(self, PSyGrid):
        assert isroutine(PSyGrid.rerun)
        pass

    def test_plot2D(self, PSyGrid):
        assert isroutine(PSyGrid.plot2D)
        pass

    def test_plot(self, PSyGrid):
        assert isroutine(PSyGrid.plot)
        pass

    def test_HR(self, PSyGrid):
        assert isroutine(PSyGrid.HR)
        pass

    def test_eq(self, PSyGrid):
        assert isroutine(PSyGrid.__eq__)
        pass


class TestPSyGridIterator:
    @fixture
    def PSyGrid(self):
        # initialize an PSyGrid for testing
        return totest.PSyGrid()

    @fixture
    def PSyGridIterator(self, PSyGrid):
        # initialize an instance of the class with defaults
        return totest.PSyGridIterator(PSyGrid)

    # test the PSyGridIterator class
    def test_init(self, PSyGridIterator):
        assert isroutine(PSyGridIterator.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        pass

    def test_next(self, PSyGridIterator):
        assert isroutine(PSyGridIterator.__next__)
        pass


class TestPSyRunView:
    @fixture
    def PSyGrid(self):
        # initialize an PSyGrid for testing
        return totest.PSyGrid()

    @fixture
    def PSyRunView(self, PSyGrid):
        # initialize an instance of the class with defaults
        return totest.PSyRunView(PSyGrid, 0)

    # test the PSyRunView class
    def test_init(self, PSyRunView):
        assert isroutine(PSyRunView.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        pass

    def test_hdf5_key(self, PSyRunView):
        assert isroutine(PSyRunView._hdf5_key)
        pass

    def test_getitem(self, PSyRunView):
        assert isroutine(PSyRunView.__getitem__)
        pass

    def test_getattr(self, PSyRunView):
        assert isroutine(PSyRunView.__getattr__)
        pass

    def test_str(self, PSyRunView, capsys):
        assert isroutine(PSyRunView.__str__)
#        with capsys.disabled():
#            print("Test")
        pass
