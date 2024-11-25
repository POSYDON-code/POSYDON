"""Unit tests of posydon/grids/psygrid.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.grids.psygrid as totest
# aliases
np = totest.np

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns
from inspect import isclass, isroutine

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
        assert dir(totest) == elements, "There might be added or removed "+\
               "objects without an update on the unit test."

    def test_instance_HDF5_MEMBER_SIZE(self):
        assert isinstance(totest.HDF5_MEMBER_SIZE, (int, float)),\
               "HDF5_MEMBER_SIZE is of type: "+\
               str(type(totest.HDF5_MEMBER_SIZE))

    def test_instance_H5_UNICODE_DTYPE(self):
        assert isinstance(totest.H5_UNICODE_DTYPE, np.dtype),\
               "H5_UNICODE_DTYPE is of type: "+\
               str(type(totest.H5_UNICODE_DTYPE))

    def test_instance_H5_REC_STR_DTYPE(self):
        assert isinstance(totest.H5_REC_STR_DTYPE, str),\
               "H5_REC_STR_DTYPE is of type: "+\
               str(type(totest.H5_REC_STR_DTYPE))

    def test_instance_VALID_KEYS(self):
        assert isinstance(totest.VALID_KEYS, list),\
               "VALID_KEYS is of type: "+str(type(totest.VALID_KEYS))

    def test_instance_WARN_VALUES(self):
        assert isinstance(totest.WARN_VALUES, list),\
               "WARN_VALUES is of type: "+str(type(totest.WARN_VALUES))

    def test_instance_N_FLAGS(self):
        assert isinstance(totest.N_FLAGS, int),\
               "N_FLAGS is of type: "+str(type(totest.N_FLAGS))

    def test_instance_N_FLAGS_SINGLE(self):
        assert isinstance(totest.N_FLAGS_SINGLE, int),\
               "N_FLAGS_SINGLE is of type: "+str(type(totest.N_FLAGS_SINGLE))

    def test_instance_TERMINATION_FLAG_COLUMNS(self):
        assert isinstance(totest.TERMINATION_FLAG_COLUMNS, list),\
               "TERMINATION_FLAG_COLUMNS is of type: "+\
               str(type(totest.TERMINATION_FLAG_COLUMNS))

    def test_instance_TERMINATION_FLAG_COLUMNS_SINGLE(self):
        assert isinstance(totest.TERMINATION_FLAG_COLUMNS_SINGLE, list),\
               "TERMINATION_FLAG_COLUMNS_SINGLE is of type: "+\
               str(type(totest.TERMINATION_FLAG_COLUMNS_SINGLE))

    def test_instance_DEFAULT_BINARY_HISTORY_COLS(self):
        assert isinstance(totest.DEFAULT_BINARY_HISTORY_COLS, list),\
               "DEFAULT_BINARY_HISTORY_COLS is of type: "+\
               str(type(totest.DEFAULT_BINARY_HISTORY_COLS))

    def test_instance_DEFAULT_STAR_HISTORY_COLS(self):
        assert isinstance(totest.DEFAULT_STAR_HISTORY_COLS, list),\
               "DEFAULT_STAR_HISTORY_COLS is of type: "+\
               str(type(totest.DEFAULT_STAR_HISTORY_COLS))

    def test_instance_DEFAULT_SINGLE_HISTORY_COLS(self):
        assert isinstance(totest.DEFAULT_SINGLE_HISTORY_COLS, list),\
               "DEFAULT_SINGLE_HISTORY_COLS is of type: "+\
               str(type(totest.DEFAULT_SINGLE_HISTORY_COLS))

    def test_instance_DEFAULT_EEP_HISTORY_COLS(self):
        assert isinstance(totest.DEFAULT_EEP_HISTORY_COLS, list),\
               "DEFAULT_EEP_HISTORY_COLS is of type: "+\
               str(type(totest.DEFAULT_EEP_HISTORY_COLS))

    def test_instance_DEFAULT_PROFILE_COLS(self):
        assert isinstance(totest.DEFAULT_PROFILE_COLS, list),\
               "DEFAULT_PROFILE_COLS is of type: "+\
               str(type(totest.DEFAULT_PROFILE_COLS))

    def test_instance_EXTRA_STAR_COLS_AT_HE_DEPLETION(self):
        assert isinstance(totest.EXTRA_STAR_COLS_AT_HE_DEPLETION, list),\
               "EXTRA_STAR_COLS_AT_HE_DEPLETION is of type: "+\
               str(type(totest.EXTRA_STAR_COLS_AT_HE_DEPLETION))

    def test_instance_DEFAULT_HISTORY_DS_EXCLUDE(self):
        assert isinstance(totest.DEFAULT_HISTORY_DS_EXCLUDE, list),\
               "DEFAULT_HISTORY_DS_EXCLUDE is of type: "+\
               str(type(totest.DEFAULT_HISTORY_DS_EXCLUDE))

    def test_instance_DEFAULT_PROFILE_DS_EXCLUDE(self):
        assert isinstance(totest.DEFAULT_PROFILE_DS_EXCLUDE, list),\
               "DEFAULT_PROFILE_DS_EXCLUDE is of type: "+\
               str(type(totest.DEFAULT_PROFILE_DS_EXCLUDE))

    def test_instance_EXTRA_COLS_DS_EXCLUDE(self):
        assert isinstance(totest.EXTRA_COLS_DS_EXCLUDE, list),\
               "EXTRA_COLS_DS_EXCLUDE is of type: "+\
               str(type(totest.EXTRA_COLS_DS_EXCLUDE))

    def test_instance_GRIDPROPERTIES(self):
        assert isinstance(totest.GRIDPROPERTIES, dict),\
               "GRIDPROPERTIES is of type: "+str(type(totest.GRIDPROPERTIES))

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
               "PROPERTIES_ALLOWED is of type: "+\
               str(type(totest.PROPERTIES_ALLOWED))

    def test_instance_PROPERTIES_TO_BE_SET(self):
        assert isinstance(totest.PROPERTIES_TO_BE_SET, list),\
               "PROPERTIES_TO_BE_SET is of type: "+\
               str(type(totest.PROPERTIES_TO_BE_SET))

    def test_instance_PROPERTIES_TO_BE_NONE(self):
        assert isinstance(totest.PROPERTIES_TO_BE_NONE, dict),\
               "PROPERTIES_TO_BE_NONE is of type: "+\
               str(type(totest.PROPERTIES_TO_BE_NONE))

    def test_instance_PROPERTIES_TO_BE_CONSISTENT(self):
        assert isinstance(totest.PROPERTIES_TO_BE_CONSISTENT, list),\
               "PROPERTIES_TO_BE_CONSISTENT is of type: "+\
               str(type(totest.PROPERTIES_TO_BE_CONSISTENT))

    def test_instance_ALL_PROPERTIES(self):
        assert isinstance(totest.ALL_PROPERTIES, list),\
               "ALL_PROPERTIES is of type: "+str(type(totest.ALL_PROPERTIES))

    def test_instance_join_grids(self):
        assert isroutine(totest.join_grids)


class TestValues:
    # check that the values fit
    def test_value_HDF5_MEMBER_SIZE(self):
        assert totest.HDF5_MEMBER_SIZE == 2147483647

    def test_value_H5_UNICODE_DTYPE(self):
        assert totest.H5_UNICODE_DTYPE == totest.h5py.string_dtype()

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
    # test functions
    def test_downsample_history(self):
        # missing argument
        with raises(TypeError, match="missing 4 required positional"+\
                    " arguments: 'bh', 'h1', 'h2', and 'params'"):
            totest.downsample_history()
        pass

    def test_downsample_profile(self):
        # missing argument
        with raises(TypeError, match="missing 2 required positional"+\
                    " arguments: 'profile' and 'params'"):
            totest.downsample_profile()
        pass

    def test_join_grids(self):
        # missing argument
        with raises(TypeError, match="missing 2 required positional"+\
                    " arguments: 'input_paths' and 'output_path'"):
            totest.join_grids()
        pass


class TestPSyGrid:
    @fixture
    def PSyGrid(self):
        # initialize an instance of the class with defaults
        return totest.PSyGrid()

    # test the PSyGrid class
    def test_init(self, PSyGrid):
        assert isroutine(PSyGrid.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        pass

    def test_reset(self, PSyGrid):
        assert isroutine(PSyGrid._reset)
        pass

    def test_make_compression_args(self, PSyGrid):
        assert isroutine(PSyGrid._make_compression_args)
        pass

    def test_discover_eeps(self, PSyGrid):
        assert isroutine(PSyGrid._discover_eeps)
        pass

    def test_say(self, PSyGrid):
        assert isroutine(PSyGrid._say)
        pass

    def test_generate_config(self, PSyGrid):
        assert isroutine(PSyGrid.generate_config)
        pass

    def test_create(self, PSyGrid):
        assert isroutine(PSyGrid.create)
        pass

    def test_create_psygrid(self, PSyGrid):
        assert isroutine(PSyGrid._create_psygrid)
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

    def test_str(self, PSyRunView):
        assert isroutine(PSyRunView.__str__)
        pass
