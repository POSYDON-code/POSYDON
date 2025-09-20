"""Unit tests of posydon/grids/post_processing.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.grids.post_processing as totest
# aliases
np = totest.np

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns
from inspect import isclass, isroutine
import h5py
import json
import os
from posydon.utils.posydonwarning import InappropriateValueWarning,\
                                         ReplaceValueWarning,\
                                         POSYDONWarning
from posydon.grids.psygrid import PSyGrid
from posydon.config import PATH_TO_POSYDON, PATH_TO_POSYDON_DATA
from posydon.unit_tests._helper_functions_for_tests.psygrid import get_PSyGrid

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = {'BinaryStar', 'CC_quantities',\
                    'CEE_parameters_from_core_abundance_thresholds',\
                    'Catch_POSYDON_Warnings',\
                    'DEFAULT_MARKERS_COLORS_LEGENDS', 'SN_MODELS', 'Pwarn',\
                    'STAR_STATES_CC', 'SingleStar', 'StepSN',\
                    'TF1_POOL_STABLE', '__authors__', '__builtins__',\
                    '__cached__', '__credits__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__',\
                    'add_post_processed_quantities',\
                    'assign_core_collapse_quantities_none',\
                    'calculate_Patton20_values_at_He_depl',\
                    'check_state_of_star', 'combine_TF12', 'copy', 'np',\
                    'post_process_grid', 'print_CC_quantities', 'tqdm',\
                    'get_SN_MODEL'}
        totest_elements = set(dir(totest))
        missing_in_test = elements - totest_elements
        assert len(missing_in_test) == 0, "There are missing objects in "\
                                          +f"{totest.__name__}: "\
                                          +f"{missing_in_test}. Please "\
                                          +"check, whether they have been "\
                                          +"removed on purpose and update "\
                                          +"this unit test."
        new_in_test = totest_elements - elements
        assert len(new_in_test) == 0, "There are new objects in "\
                                      +f"{totest.__name__}: {new_in_test}. "\
                                      +"Please check, whether they have been "\
                                      +"added on purpose and update this "\
                                      +"unit test."

    def test_instance_CC_quantities(self):
        assert isinstance(totest.CC_quantities, (list)), "CC_quantities is "\
               + "of type: " + str(type(totest.CC_quantities))

    def test_instance_assign_core_collapse_quantities_none(self):
        assert isroutine(totest.assign_core_collapse_quantities_none)

    def test_instance_print_CC_quantities(self):
        assert isroutine(totest.print_CC_quantities)

    def test_instance_post_process_grid(self):
        assert isroutine(totest.post_process_grid)

    def test_instance_add_post_processed_quantities(self):
        assert isroutine(totest.add_post_processed_quantities)


class TestValues:
    # check that the values fit
    def test_value_CC_quantities(self):
        for q in ['state', 'SN_type', 'f_fb', 'mass', 'spin',\
                  'm_disk_accreted', 'm_disk_radiated',\
                  'CO_interpolation_class', 'M4', 'mu4', 'h1_mass_ej',\
                  'he4_mass_ej']:
            assert q in totest.CC_quantities


class TestFunctions:
    @fixture
    def star(self):
    # initialize a SingleStar instance, which is a required argument
        test_star = totest.SingleStar()
        test_star.state = "TestState"
        test_star.mass = 1.0
        test_star.spin = 0.0
        test_star.SN_type = "TestSNType"
        test_star.f_fb = 0.1
        test_star.m_disk_accreted = 0.2
        test_star.m_disk_radiated = 0.3
        test_star.M4 = 0.4
        test_star.mu4 = 0.5
        test_star.h1_mass_ej = np.nan
        test_star.he4_mass_ej = None
        return test_star

    @fixture
    def star_history(self):
        # a temporary star history for testing
        return np.array([(1.0, 1.0, 0.79, 0.2, 0.0, 0.0, 0.0, 0.79, 0.2,\
                          0.0, 0.0, 0.0, -1.0, -10.0, -1.0, 0.1, 0.0, 0.0),\
                         (1.0e+2, 1.1, 0.09, 0.9, 0.0, 0.0, 0.0, 0.79, 0.2,\
                          0.0, 0.0, 0.0, -1.0, -10.0, -1.0, 0.1, 1.0, 0.0),\
                         (1.0e+3, 3.2, 0.0, 0.2, 0.4, 0.09, 0.3, 0.79, 0.2,\
                          0.0, 0.0, 0.0, -10.0, -2.0, -2.0, 0.1, 2.0, 0.0),\
                         (1.0e+4, 2.3, 0.0, 0.0, 0.2, 0.1, 0.6, 0.79, 0.2,\
                          0.0, 0.0, 0.0, -10.0, -3.1, -3.0, 11.0, 13.0, 12.0)],\
                        dtype=[('star_age', '<f8'), ('log_R', '<f8'),\
                               ('center_h1', '<f8'), ('center_he4', '<f8'),\
                               ('center_c12', '<f8'), ('center_n14', '<f8'),\
                               ('center_o16', '<f8'), ('surface_h1', '<f8'),\
                               ('surface_he4', '<f8'), ('surface_c12', '<f8'),\
                               ('surface_n14', '<f8'), ('surface_o16', '<f8'),\
                               ('log_LH', '<f8'), ('log_LHe', '<f8'),\
                               ('log_Lnuc', '<f8'), ('center_gamma', '<f8'),\
                               ('he_core_mass', '<f8'),\
                               ('co_core_mass', '<f8')])

    @fixture
    def binary_history(self):
        # a temporary binary history for testing
        return np.array([(19.1, 5.2, 1.0, 1.0, -99.1, -99.2, -98.1, -98.2,\
                          -97.1, -97.2),\
                         (18.1, 5.2, 1.1, 1.0e+2, -99.1, -99.2, -98.1, -98.2,\
                          -97.1, -97.2),\
                         (16.1, 5.2, 1.2, 1.0e+3, -99.1, -99.2, -98.1, -98.2,\
                          -97.1, -97.2),\
                         (15.1, 5.1, 1.3, 1.0e+4, -99.1, -99.2, -98.1, -98.2,\
                          -97.1, -97.2)],\
                        dtype=[('star_1_mass', '<f8'), ('star_2_mass', '<f8'),\
                               ('period_days', '<f8'), ('age', '<f8'),\
                               ('lg_mstar_dot_1', '<f8'),\
                               ('lg_mstar_dot_2', '<f8'),\
                               ('lg_system_mdot_1', '<f8'),\
                               ('lg_system_mdot_2', '<f8'),\
                               ('lg_wind_mdot_1', '<f8'),\
                               ('lg_wind_mdot_2', '<f8')])

    @fixture
    def profile(self):
        # a temporary profile for testing
        return np.array([(15.1, 10**2.3, 1.0, 1.0, 1.0, 0.79, 0.2, 0.01, 0.9,\
                          0.8, 0.2),\
                         (1.1, 1.0e+2, 1.0, 1.0, 1.0, 0.79, 0.2, 0.01, 0.5,\
                          0.4, 1.2),\
                         (0.1, 1.0e-1, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,\
                          0.0, 1.9)],\
                        dtype=[('mass', '<f8'), ('radius', '<f8'),\
                               ('logRho', '<f8'), ('omega', '<f8'),\
                               ('energy', '<f8'),\
                               ('x_mass_fraction_H', '<f8'),\
                               ('y_mass_fraction_He', '<f8'),\
                               ('z_mass_fraction_metals', '<f8'),\
                               ('neutral_fraction_H', '<f8'),\
                               ('neutral_fraction_He', '<f8'),\
                               ('avg_charge_He', '<f8')])

    @fixture
    def grid_path(self, tmp_path, binary_history, star_history, profile):
        # a path to a psygrid file for testing
        return get_PSyGrid(tmp_path, 1, binary_history, star_history, profile)

    @fixture
    def link_SN_data(self):
        # link SN data in unit tests in case it is not available otherwise
        links = {}
        # create "POSYDON_data" directory if not existing
        if not os.path.exists(PATH_TO_POSYDON_DATA):
            links[PATH_TO_POSYDON_DATA] = "remove later"
            os.mkdir(PATH_TO_POSYDON_DATA)
        # create links if needed
        for SN_engine in ["Sukhbold+16", "Patton+Sukhbold20", "Couch+2020"]:
            engine_path = os.path.join(PATH_TO_POSYDON_DATA, SN_engine)
            if not os.path.exists(engine_path):
                test_engine = os.path.join(PATH_TO_POSYDON, "posydon",\
                                           "unit_tests", "_data",\
                                           "POSYDON_data", SN_engine)
                links[engine_path] = test_engine
                os.symlink(test_engine, engine_path)
        yield
        # remove created links
        for engine_path in links:
            if engine_path != PATH_TO_POSYDON_DATA:
                os.unlink(engine_path)
        # remove created directory
        if PATH_TO_POSYDON_DATA in links:
            os.rmdir(PATH_TO_POSYDON_DATA)

    @fixture
    def grid(self):
        # initialize a PSyGrid instance, which is a required argument
        test_grid = PSyGrid()
        test_grid.MESA_dirs = ["Test"]
        return test_grid

    # test functions
    def test_assign_core_collapse_quantities_none(self):
        # missing argument
        with raises(TypeError, match="missing 2 required positional "\
                                     +"arguments: 'EXTRA_COLUMNS' and "\
                                     +"'star_i'"):
            totest.assign_core_collapse_quantities_none()
        # bad input
        with raises(TypeError, match="'NoneType' object is not subscriptable"):
            totest.assign_core_collapse_quantities_none(None, 1)
        # bad input
        test_EXTRA_COLUMNS = {}
        with raises(ValueError, match="'star_i' should be 1 or 2."):
            totest.assign_core_collapse_quantities_none(test_EXTRA_COLUMNS, 0)
        # bad input
        with raises(TypeError, match="'SN_MODEL_NAME' should be a string, a "\
                                     +"list of strings, or None."):
            totest.assign_core_collapse_quantities_none(test_EXTRA_COLUMNS, 1,\
                                                        1)
        # bad input
        with raises(KeyError, match="'S1_SN_MODEL_v2_01_state'"):
            totest.assign_core_collapse_quantities_none(test_EXTRA_COLUMNS, 1)
        # examples: all models in SN_MODELS.py
        for s in [1, 2]:
            test_EXTRA_COLUMNS = {}
            for m in totest.SN_MODELS.keys():
                for q in totest.CC_quantities:
                    test_EXTRA_COLUMNS[f'S{s}_{m}_{q}'] = []
            totest.assign_core_collapse_quantities_none(test_EXTRA_COLUMNS, s)
            for m in totest.SN_MODELS.keys():
                for q in totest.CC_quantities:
                    assert test_EXTRA_COLUMNS[f'S{s}_{m}_{q}'] == [None]
        # examples: given model
        for s in [1, 2]:
            test_EXTRA_COLUMNS = {}
            for m in ["TESTMODEL"]:
                for q in totest.CC_quantities:
                    test_EXTRA_COLUMNS[f'S{s}_{m}_{q}'] = []
            totest.assign_core_collapse_quantities_none(test_EXTRA_COLUMNS, s,\
                                                        SN_MODEL_NAME=m)
            for m in ["TESTMODEL"]:
                for q in totest.CC_quantities:
                    assert test_EXTRA_COLUMNS[f'S{s}_{m}_{q}'] == [None]

    def test_print_CC_quantities(self, star, capsys):
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'star'"):
            totest.print_CC_quantities()
        # bad input
        with warns(InappropriateValueWarning) as winfo:
            totest.print_CC_quantities(None)
        assert "Failed to print star values!\nWarning in preSN: 'NoneType' "\
               + "object has no attribute 'state'"\
               == winfo._list[0].message.args[0]
        output = capsys.readouterr().out.split("\n")
        assert len(output) == 4+1
        for h in ["mechanism", "state", "SN type", "f_fb", "mass [Msun]",\
                  "spin", "m_disk_accreted [Msun]", "m_disk_radiated [Msun]",\
                  "M4 [m/Msun]", "mu4 [(dm/Msun)/(dr/1000km)]",\
                  "h1_mass_ej [Msun]", "he4_mass_ej [Msun]"]:
            assert h in output[1]
        # bad input
        with warns(InappropriateValueWarning) as winfo:
            totest.print_CC_quantities(None, "TESTMODEL")
        assert "Failed to print star values!\nWarning in TESTMODEL: "\
               + "'NoneType' object has no attribute 'spin'"\
               == winfo._list[0].message.args[0]
        # examples: pre SN
        totest.print_CC_quantities(star)
        output = capsys.readouterr().out.split("\n")
        assert len(output) == 5+1
        for h in ["mechanism", "state", "SN type", "f_fb", "mass [Msun]",\
                  "spin", "m_disk_accreted [Msun]", "m_disk_radiated [Msun]",\
                  "M4 [m/Msun]", "mu4 [(dm/Msun)/(dr/1000km)]",\
                  "h1_mass_ej [Msun]", "he4_mass_ej [Msun]"]:
            assert h in output[1]
        for v in ["PRE SN STAR", star.state,\
                  "{:7.2f} {:12.2f}".format(star.mass, star.spin)]:
            assert v in output[3]
        # examples: SN MODEL
        totest.print_CC_quantities(star, "TESTMODEL")
        output = capsys.readouterr().out.split("\n")
        assert len(output) == 1+1
        fmt = "{:1.2f} {:13.2f} {:12.2f} {:20.2f} {:20.2f} {:20.2f} {:20.2f} "\
              +"{:20.2f} {:20.2f}"
        for v in ["TESTMODEL", star.state, star.SN_type, fmt.format(star.f_fb,\
                   star.mass, star.spin, star.m_disk_accreted,\
                   star.m_disk_radiated, star.M4, star.mu4, star.h1_mass_ej,\
                   np.nan)]:
            assert v in output[0]

    def test_post_process_grid(self, grid_path, link_SN_data, capsys,\
                               monkeypatch):
        def check_EXTRA_COLUMNS(EXTRA_COLUMNS, n_runs, keys, run_range=None):
            """Function to check the EXTRA_COLUMNS

            Parameters
            ----------
            EXTRA_COLUMNS : dict
                New columns to check through.
            n_runs : int
                Count of runs.
            keys : list of str
                List of keys in EXTRA_COLUMNS.
            run_range : None or size 2 tuple (default: None)
                The range of entries to check. If None it checks 0 till
                n_runs-1.
            """
            if run_range is None:
                run_range = (0, n_runs)
            for k in keys:
                assert k in EXTRA_COLUMNS
                assert len(EXTRA_COLUMNS[k]) == n_runs
                for i in range(run_range[0], run_range[1]):
                    if ((i == 0) and (k != "mt_history")):
                        # check missing run: all None
                        assert EXTRA_COLUMNS[k][i] is None, f"i={i}, k={k}"
                    elif (("S2_SN_MODEL" in k) and (i in [1, 2, 4, 6])):
                        # check star2: all None beside state
                        assert EXTRA_COLUMNS[k][i] is None, f"i={i}, k={k}"
                    elif (("S1_SN_MODEL" in k) and (i in [2, 4, 6])):
                        # check SN of star1 for None
                        assert EXTRA_COLUMNS[k][i] is None, f"i={i}, k={k}"
        def check_EXTRA_COLUMNS_single(EXTRA_COLUMNS, n_runs, keys):
            """Function to check the EXTRA_COLUMNS of single star girds
            
            Note: The non-core collapse runs 0 and 6 will be skipped.

            Parameters
            ----------
            EXTRA_COLUMNS : dict
                New columns to check through.
            n_runs : int
                Count of runs.
            keys : list of str
                List of keys in EXTRA_COLUMNS. All 'S2' and the 'mt_history'
                will be checked for being not there.
            """
            for k in keys:
                if (('S2' in k) or (k == 'mt_history')):
                    assert k not in EXTRA_COLUMNS
            check_EXTRA_COLUMNS(EXTRA_COLUMNS, n_runs,\
                                [k for k in keys if (('S2' not in k)\
                                                     and (k!='mt_history'))],\
                                run_range=(1,6))
        def mock_check_state_of_star(star, i=None, star_CO=False):
            if star_CO:
                raise TypeError("Testing exception.")
            else:
                return 'WD'
        def mock_collapse_star(SN, star):
            for quantity in totest.CC_quantities:
                setattr(star, quantity, True)
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'grid'"):
            totest.post_process_grid()
        # bad input
        with raises(AttributeError, match="'NoneType' object has no "\
                                          +"attribute 'MESA_dirs'"):
            totest.post_process_grid(None)
        # bad input
        with raises(ValueError, match="Index range should have dim=2!"):
            totest.post_process_grid(None, index=[])
        # bad input
        with raises(TypeError, match="The argument `index` should be None, "\
                                     +"an integer, or a list of two "\
                                     +"integers."):
            totest.post_process_grid(None, index="Test")
        # get test_PSyGrid
        try:
            test_PSyGrid = PSyGrid()
            test_PSyGrid.load(grid_path)
        except: # skip test as test on load should fail
            assert "/" in grid_path
            return
        keys = ['mt_history']
        for s in [1,2]:
            for q in ['avg_c_in_c_core_at_He_depletion',\
                      'co_core_mass_at_He_depletion', 'state',\
                      'surface_other', 'center_other']:
                keys += [f'S{s}_{q}']
            for q in ['lambda_CE', 'm_core_CE', 'r_core_CE']:
                for v in [1, 10, 30, 'pure_He_star_10']:
                    keys += [f'S{s}_{q}_{v}cent']
            for m in totest.SN_MODELS:
                for q in totest.CC_quantities:
                    if q == 'state':
                        q = 'CO_type'
                    keys += [f'S{s}_{m}_{q}']
        # examples: all
        with warns(POSYDONWarning): # warnings from SN
            MESA_dirs, EXTRA_COLUMNS = totest.post_process_grid(test_PSyGrid)
        assert MESA_dirs == test_PSyGrid.MESA_dirs
        check_EXTRA_COLUMNS(EXTRA_COLUMNS, 7, keys)
        # examples: run 2 only
        MESA_dirs, EXTRA_COLUMNS = totest.post_process_grid(test_PSyGrid,\
                                                            index=2)
        assert MESA_dirs == test_PSyGrid.MESA_dirs[2:3]
        for k in keys:
            assert k in EXTRA_COLUMNS
            assert len(EXTRA_COLUMNS[k]) == 1
        # examples: 1 and 2
        with warns(POSYDONWarning): # warnings from SN
            MESA_dirs, EXTRA_COLUMNS = totest.post_process_grid(test_PSyGrid,\
                                                                index=[1,3])
        assert MESA_dirs == test_PSyGrid.MESA_dirs[1:3]
        for k in keys:
            assert k in EXTRA_COLUMNS
            assert len(EXTRA_COLUMNS[k]) == 2
        # examples: HMS-HMS
        with warns(POSYDONWarning): # warnings from SN
            with warns(InappropriateValueWarning,\
                       match="ended with TF1=gamma_center_limit however the "\
                             +"star has center_gamma < 10. This star cannot "\
                             +"go through step_SN appending NONE compact "\
                             +"object properties!"):
                with warns(InappropriateValueWarning,\
                           match="ended with TF=max_age and IC=stable_MT. "\
                                 +"This star cannot go through step_SN "\
                                 +"appending NONE compact object properties!"):
                    MESA_dirs, EXTRA_COLUMNS = totest.post_process_grid(\
                                                test_PSyGrid, star_2_CO=False)
        assert MESA_dirs == test_PSyGrid.MESA_dirs
        check_EXTRA_COLUMNS(EXTRA_COLUMNS, 7, keys)
        # examples: less SN MODELS
        TEST_MODELS = {}
        for m,v in totest.SN_MODELS.items():
            if m != 'SN_MODEL_v2_01':
                TEST_MODELS[m] = v
        with warns(POSYDONWarning): # warnings from SN
            MESA_dirs, EXTRA_COLUMNS = totest.post_process_grid(test_PSyGrid,\
                                        SN_MODELS=TEST_MODELS)
        assert MESA_dirs == test_PSyGrid.MESA_dirs
        for k in keys:
            if 'SN_MODEL_v2_01' in k:
                assert k not in EXTRA_COLUMNS
        check_EXTRA_COLUMNS(EXTRA_COLUMNS, 7,\
                            [k for k in keys if 'SN_MODEL_v2_01' not in k])
        # examples: single
        with warns(POSYDONWarning): # warnings from SN
            MESA_dirs, EXTRA_COLUMNS = totest.post_process_grid(test_PSyGrid,\
                                        single_star=True)
        assert MESA_dirs == test_PSyGrid.MESA_dirs
        check_EXTRA_COLUMNS_single(EXTRA_COLUMNS, 7, keys)
        # examples: failing check_state_of_star
        with monkeypatch.context() as mp:
            mp.setattr(totest, "check_state_of_star", mock_check_state_of_star)
            with warns(POSYDONWarning): # warnings from SN
                MESA_dirs, EXTRA_COLUMNS = totest.post_process_grid(test_PSyGrid)
        assert MESA_dirs == test_PSyGrid.MESA_dirs
        check_EXTRA_COLUMNS(EXTRA_COLUMNS, 7, keys)
        # examples: failing collapse_star
        with monkeypatch.context() as mp:
            mp.setattr(totest.StepSN, "collapse_star", mock_collapse_star)
            with warns(InappropriateValueWarning, match="is not a string!"):
                with warns(InappropriateValueWarning, match="is not a float!"):
                    with warns(InappropriateValueWarning, match="is not a "\
                                                                +"float nor "\
                                                                +"None!"):
                        MESA_dirs, EXTRA_COLUMNS = totest.post_process_grid(\
                                                    test_PSyGrid)
        assert MESA_dirs == test_PSyGrid.MESA_dirs
        check_EXTRA_COLUMNS(EXTRA_COLUMNS, 7, keys)
        # examples: single and failing collapse_star
        with monkeypatch.context() as mp:
            mp.setattr(totest.StepSN, "collapse_star", mock_collapse_star)
            with warns(InappropriateValueWarning, match="is not a string!"):
                with warns(InappropriateValueWarning, match="is not a float!"):
                    with warns(InappropriateValueWarning, match="is not a "\
                                                                +"float nor "\
                                                                +"None!"):
                        MESA_dirs, EXTRA_COLUMNS = totest.post_process_grid(\
                                                    test_PSyGrid,\
                                                    single_star=True)
        assert MESA_dirs == test_PSyGrid.MESA_dirs
        check_EXTRA_COLUMNS_single(EXTRA_COLUMNS, 7, keys)
        # examples: verbose
        with warns(POSYDONWarning): # warnings from SN
            with warns(InappropriateValueWarning, match="Failed to print "\
                                                        +"star values!"):
                MESA_dirs, EXTRA_COLUMNS = totest.post_process_grid(\
                                            test_PSyGrid, verbose=True)
        output = capsys.readouterr().out
        assert MESA_dirs == test_PSyGrid.MESA_dirs
        check_EXTRA_COLUMNS(EXTRA_COLUMNS, 7, keys)
        assert "Error during" in output
        assert "core collapse prescrition!" in output
        assert "The error was raised by" in output
        assert "in CEE_parameters_from_core_abundance_thresholds" in output
        assert "The exception was raised by" in output
        with warns(POSYDONWarning): 
            with warns(ReplaceValueWarning, match="While accessing "
                                                   +"abundances in star"):
                MESA_dirs, EXTRA_COLUMNS = totest.post_process_grid(\
                                            test_PSyGrid,\
                                            verbose=True)
                
        assert "in check_state_of_star(star_2) with IC=" in output
        # examples: single and verbose
        with warns(POSYDONWarning): # warnings from SN
            with warns(InappropriateValueWarning, match="Failed to print "\
                                                        +"star values!"):
                MESA_dirs, EXTRA_COLUMNS = totest.post_process_grid(\
                                            test_PSyGrid, single_star=True,\
                                            verbose=True)
        output = capsys.readouterr().out
        assert MESA_dirs == test_PSyGrid.MESA_dirs
        check_EXTRA_COLUMNS_single(EXTRA_COLUMNS, 7, keys)
        assert "Error during" in output
        assert "core collapse prescrition!" in output
        assert "in CEE_parameters_from_core_abundance_thresholds" in output
        assert "The exception was raised by" in output
        with warns(POSYDONWarning): 
            with warns(ReplaceValueWarning, match="While accessing "
                                                   +"abundances in star"):
                MESA_dirs, EXTRA_COLUMNS = totest.post_process_grid(\
                                            test_PSyGrid, single_star=True,\
                                            verbose=True)

    def test_add_post_processed_quantities(self, grid, monkeypatch):
        def mock_add_column(colname, array, where="final_values",\
                            overwrite=True):
            self.newcolumns[colname] = array
        # missing argument
        with raises(TypeError, match="missing 3 required positional "\
                                     +"arguments: 'grid', "\
                                     +"'MESA_dirs_EXTRA_COLUMNS', and "\
                                     +"'EXTRA_COLUMNS'"):
            totest.add_post_processed_quantities()
        # bad input
        with raises(AttributeError, match="'NoneType' object has no "\
                                          +"attribute 'MESA_dirs'"):
            totest.add_post_processed_quantities(None, None, None)
        # bad input
        with raises(ValueError, match="EXTRA_COLUMNS do not follow the "\
                                      +"correct order of grid!"):
            totest.add_post_processed_quantities(grid, ["Unit"], None)
        # examples: nothing to add
        totest.add_post_processed_quantities(grid, ["Test"], {})
        assert grid.final_values is None
        # examples: add columns
        with monkeypatch.context() as mp:
            mp.setattr(grid, "add_column", mock_add_column)
            EXTRA_COLUMNS = {'test_state': ["unit"], 'test_type': ["test"],\
                             'test_class': ["unittest"], 'mt_history': [""],\
                             'value': [1.0], 'values': [0.0, 1.0]}
            self.newcolumns = {}
            totest.add_post_processed_quantities(grid, ["Test"], EXTRA_COLUMNS)
            for k,v in EXTRA_COLUMNS.items():
                for i in range(len(v)):
                    assert self.newcolumns[k][i] == v[i]
