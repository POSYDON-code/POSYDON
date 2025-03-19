"""Unit tests of posydon/grids/termination_flags.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.grids.termination_flags as totest
# aliases
np = totest.np
os = totest.os

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns
from inspect import isroutine
from posydon.grids.psygrid import PSyGrid
from posydon.utils.posydonwarning import InappropriateValueWarning

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['LG_MTRANSFER_RATE_THRESHOLD',\
                    'MIN_COUNT_INITIAL_RLO_BOUNDARY', 'Pwarn',\
                    'RL_RELATIVE_OVERFLOW_THRESHOLD',\
                    'STAR_HISTORY_VARIABLES', 'TF1_POOL_ERROR',\
                    'TF1_POOL_INITIAL_RLO', 'TF1_POOL_STABLE',\
                    'TF1_POOL_UNSTABLE', 'TF2_POOL_INITIAL_RLO',\
                    'TF2_POOL_NO_RLO', '__authors__', '__builtins__',\
                    '__cached__', '__doc__', '__file__', '__loader__',\
                    '__name__', '__package__', '__spec__',\
                    'check_state_from_history',\
                    'cumulative_mass_transfer_flag',\
                    'get_detected_initial_RLO', 'get_flag_from_MESA_output',\
                    'get_flags_from_MESA_run', 'get_mass_transfer_flag',\
                    'get_nearest_known_initial_RLO', 'gzip',\
                    'infer_interpolation_class', 'infer_mass_transfer_case',\
                    'infer_star_state', 'np', 'os']
        assert dir(totest) == elements, "There might be added or removed "\
                                        +"objects without an update on the "\
                                        +"unit test."

    def test_instance_STAR_HISTORY_VARIABLES(self):
        assert isinstance(totest.STAR_HISTORY_VARIABLES, list)

    def test_instance_get_flag_from_MESA_output(self):
        assert isroutine(totest.get_flag_from_MESA_output)

    def test_instance_get_mass_transfer_flag(self):
        assert isroutine(totest.get_mass_transfer_flag)

    def test_instance_check_state_from_history(self):
        assert isroutine(totest.check_state_from_history)

    def test_instance_get_flags_from_MESA_run(self):
        assert isroutine(totest.get_flags_from_MESA_run)

    def test_instance_infer_interpolation_class(self):
        assert isroutine(totest.infer_interpolation_class)

    def test_instance_get_detected_initial_RLO(self):
        assert isroutine(totest.get_detected_initial_RLO)

    def test_instance_get_nearest_known_initial_RLO(self):
        assert isroutine(totest.get_nearest_known_initial_RLO)


class TestValues:
    # check that the values fit
    def test_value_STAR_HISTORY_VARIABLES(self):
        for v in ["surface_h1", "center_h1", "center_he4", "center_c12",\
                  "log_LH", "log_LHe", "log_Lnuc"]:
            # check required values
            assert v in totest.STAR_HISTORY_VARIABLES, "missing entry"


class TestFunctions:
    @fixture
    def out_path(self, tmp_path):
        # a temporary path to out file for testing
        path = os.path.join(tmp_path, "out.txt")
        with open(path, "w") as test_file:
            test_file.write("min_timestep_limit Reached TPAGB\n")
        return path

    @fixture
    def out_path2(self, tmp_path):
        # a temporary path to out file for testing
        path = os.path.join(tmp_path, "out2.txt")
        with open(path, "w") as test_file:
            test_file.write("Terminate: Unit test 2\n")
        return path

    @fixture
    def out_path3(self, tmp_path):
        # a temporary path to out file for testing
        path = os.path.join(tmp_path, "out3.txt")
        with open(path, "w") as test_file:
            test_file.write("termination code: Unit test 3\n")
            test_file.write("termination code: min_timestep_limit\n")
        return path

    @fixture
    def out_path4(self, tmp_path):
        # a temporary path to out file for testing
        path = os.path.join(tmp_path, "out4.txt")
        with open(path, "w") as test_file:
            test_file.write("Unit test 4\n")
        return path

    @fixture
    def out_path5(self, tmp_path):
        # a temporary path to out file for testing
        path = os.path.join(tmp_path, "out5.txt")
        with open(path, "w") as test_file:
            test_file.write("")
        return path

    @fixture
    def out_path6(self, tmp_path):
        # a temporary path to out file for testing
        path = os.path.join(tmp_path, "out6.txt")
        with open(path, "w") as test_file:
            test_file.write("termination code: min_timestep_limit\n")
        os.system(f"gzip -1 {path}")
        return path

    @fixture
    def star_history(self):
        # a temporary star history for testing
        return np.array([(0.5, 0.1, 0.4, 0.5, 1.0, 1.0, 1.4),\
                         (0.5, 0.1, 0.4, 0.5, 1.0, 1.0, 1.4)],\
                        dtype=[('surface_h1', '<f8'), ('center_h1', '<f8'),\
                               ('center_he4', '<f8'), ('center_c12', '<f8'),\
                               ('log_LH', '<f8'), ('log_LHe', '<f8'),\
                               ('log_Lnuc', '<f8')])

    @fixture
    def binary_history(self):
        # a temporary binary history for testing
        return np.array([(totest.RL_RELATIVE_OVERFLOW_THRESHOLD,\
                          totest.RL_RELATIVE_OVERFLOW_THRESHOLD,\
                          totest.LG_MTRANSFER_RATE_THRESHOLD, 100.0, 0.2),\
                         (totest.RL_RELATIVE_OVERFLOW_THRESHOLD,\
                          totest.RL_RELATIVE_OVERFLOW_THRESHOLD,\
                          totest.LG_MTRANSFER_RATE_THRESHOLD, 1.0, 20.0)],\
                        dtype=[('rl_relative_overflow_1', '<f8'),\
                               ('rl_relative_overflow_2', '<f8'),\
                               ('lg_mtransfer_rate', '<f8'),\
                               ('star_1_mass', '<f8'), ('star_2_mass', '<f8')])

    @fixture
    def grid(self):
        # a temporary grid for testing
        g = PSyGrid()
        # 5x5x5 gird + late small period "00-"
        iv = np.array([(m1, m2, porb) for m1 in range(5) for m2 in range(5)\
                       for porb in range(5)] + [(0, 0, -1)],\
                      dtype=[('star_1_mass', '<f8'), ('star_2_mass', '<f8'),\
                             ('period_days', '<f8')])
        setattr(g, "initial_values", iv)
        # initial RLO on diagonal plane + on the additional point
        fv = np.array([("Terminate because of overflowing initial model"\
                        if 0<m1+m2+porb<3 else "Test",\
                        str(m1)+str(m2)+str(porb), str(m1)+str(m2)+str(porb))\
                        for m1 in range(5) for m2 in range(5) for porb in\
                        range(5)] + [("Terminate because of overflowing"+\
                                      " initial model", "00-", "00-")],\
                      dtype=[('termination_flag_1', '<U50'),\
                             ('termination_flag_3', '<U50'),\
                             ('termination_flag_4', '<U50')])
        setattr(g, "final_values", fv)
        return g

    @fixture
    def iniRLO_boundary(self):
        # a temporary iniRLO_boundary for testing: upper diagonal plane of grid
        return [{'star_1_mass': 0.0, 'star_2_mass': 0.0, 'period_days': 2.0,\
                 'termination_flag_3': '002', 'termination_flag_4': '002'},\
                {'star_1_mass': 0.0, 'star_2_mass': 1.0, 'period_days': 1.0,\
                 'termination_flag_3': '011', 'termination_flag_4': '011'},\
                {'star_1_mass': 0.0, 'star_2_mass': 2.0, 'period_days': 0.0,\
                 'termination_flag_3': '020', 'termination_flag_4': '020'},\
                {'star_1_mass': 1.0, 'star_2_mass': 0.0, 'period_days': 1.0,\
                 'termination_flag_3': '101', 'termination_flag_4': '101'},\
                {'star_1_mass': 1.0, 'star_2_mass': 1.0, 'period_days': 0.0,\
                 'termination_flag_3': '110', 'termination_flag_4': '110'},\
                {'star_1_mass': 2.0, 'star_2_mass': 0.0, 'period_days': 0.0,\
                 'termination_flag_3': '200', 'termination_flag_4': '200'}]

    # test functions
    def test_get_flag_from_MESA_output(self, out_path, out_path2, out_path3,\
                                       out_path4, out_path5, out_path6):
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'MESA_log_path'"):
            totest.get_flag_from_MESA_output()
        # examples
        assert totest.get_flag_from_MESA_output(None) ==\
               "reach cluster timelimit"
        assert totest.get_flag_from_MESA_output(out_path) ==\
               "Reached TPAGB"
        assert totest.get_flag_from_MESA_output(out_path2) ==\
               "Unit test 2"
        assert totest.get_flag_from_MESA_output(out_path3) ==\
               "Unit test 3"
        assert totest.get_flag_from_MESA_output(out_path4) ==\
               "reach cluster timelimit"
        assert totest.get_flag_from_MESA_output(out_path5) ==\
               "reach cluster timelimit"
        assert totest.get_flag_from_MESA_output(out_path6+".gz") ==\
               "min_timestep_limit"

    def test_get_mass_transfer_flag(self, star_history, binary_history):
        # missing argument
        with raises(TypeError, match="missing 3 required positional "\
                                     +"arguments: 'binary_history', "\
                                     +"'history1', and 'history2'"):
            totest.get_mass_transfer_flag()
        # bad input
        with raises(TypeError, match="'NoneType' object is not subscriptable"):
            totest.get_mass_transfer_flag(None, None, None)
        # examples: failed run
        assert totest.get_mass_transfer_flag(None, None, None,\
                mesa_flag=totest.TF1_POOL_ERROR[0]) == "None"
        # examples: initial_RLOF from flag
        assert totest.get_mass_transfer_flag(None, None, None,\
                mesa_flag=totest.TF1_POOL_INITIAL_RLO[0]) == "initial_RLOF"
        # examples: no RLO
        assert totest.get_mass_transfer_flag(binary_history, None, None) ==\
               "no_RLOF"
        # examples: case A1
        binary_history["rl_relative_overflow_1"][0] = 1\
         + totest.RL_RELATIVE_OVERFLOW_THRESHOLD
        assert totest.get_mass_transfer_flag(binary_history, star_history,\
                                             star_history) == "initial_RLOF"
        assert totest.get_mass_transfer_flag(binary_history, star_history,\
                                             star_history, start_at_RLO=True)\
               == "case_A1"
        # examples: contact
        binary_history["rl_relative_overflow_2"][0] = 1\
         + totest.RL_RELATIVE_OVERFLOW_THRESHOLD
        assert totest.get_mass_transfer_flag(binary_history, None, None) ==\
               "contact_during_MS"
        # examples: case A2
        binary_history["rl_relative_overflow_1"][0] = \
         totest.RL_RELATIVE_OVERFLOW_THRESHOLD
        assert totest.get_mass_transfer_flag(binary_history, star_history,\
                                             star_history) == "initial_RLOF"
        assert totest.get_mass_transfer_flag(binary_history, star_history,\
                                             star_history, start_at_RLO=True)\
               == "case_A2"

    def test_check_state_from_history(self, star_history, binary_history):
        # missing argument
        with raises(TypeError, match="missing 2 required positional "\
                                     +"arguments: 'history' and 'mass'"):
            totest.check_state_from_history()
        # bad input
        with raises(TypeError, match="'NoneType' object is not subscriptable"):
            totest.check_state_from_history(None, None)
        for k in totest.STAR_HISTORY_VARIABLES:
            with raises(KeyError, match=f"The data column {k} is not in the "\
                                        +"history file. It is needed for the "\
                                        +"determination of the star state."):
                totest.check_state_from_history(star_history[[col for col in\
                 totest.STAR_HISTORY_VARIABLES if col!=k]], None)
        # examples: compact object
        assert totest.check_state_from_history(None,\
                binary_history['star_1_mass']) == totest.infer_star_state(\
                 star_mass=binary_history['star_1_mass'][-1], star_CO=True)
        assert totest.check_state_from_history(None,\
                binary_history['star_1_mass'], model_index=0)\
               == totest.infer_star_state(\
                star_mass=binary_history['star_1_mass'][0], star_CO=True)
        # examples: star
        assert totest.check_state_from_history(star_history, None) ==\
               "H-rich_Core_H_burning"

    def test_get_flags_from_MESA_run(self, star_history, binary_history):
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'MESA_log_path'"):
            totest.get_flags_from_MESA_run()
        # bad input
        with raises(TypeError, match="'NoneType' object is not subscriptable"):
            totest.get_flags_from_MESA_run(None)
        # examples
        fo, fmt, fs1, fs2 = totest.get_flags_from_MESA_run(None,\
                             binary_history=binary_history,\
                             history1=star_history, history2=star_history)
        assert fo == "reach cluster timelimit"
        assert fmt == "None"
        assert fs1 == "H-rich_Core_H_burning"
        assert fs2 == "H-rich_Core_H_burning"
        fo, fmt, fs1, fs2 = totest.get_flags_from_MESA_run(None,\
                             binary_history=binary_history,\
                             history1=star_history, history2=star_history,\
                             newTF1="Test")
        assert fo == "Test"
        assert fmt == "no_RLOF"
        assert fs1 == "H-rich_Core_H_burning"
        assert fs2 == "H-rich_Core_H_burning"

    def test_infer_interpolation_class(self):
        # missing argument
        with raises(TypeError, match="missing 2 required positional "\
                                     +"arguments: 'tf1' and 'tf2'"):
            totest.infer_interpolation_class()
        # bad input
        with raises(TypeError, match="argument of type 'NoneType' is not "\
                                     +"iterable"):
            totest.infer_interpolation_class(totest.TF1_POOL_STABLE[0], None)
        # examples
        assert totest.infer_interpolation_class(None, None) == "unknown"
        for tf1 in totest.TF1_POOL_INITIAL_RLO:
            assert totest.infer_interpolation_class(tf1, None) == "initial_MT"
        for tf2 in totest.TF2_POOL_INITIAL_RLO:
            assert totest.infer_interpolation_class(None, tf2) == "initial_MT"
        for tf1 in totest.TF1_POOL_ERROR:
            assert totest.infer_interpolation_class(tf1, None) ==\
                   "not_converged"
        for tf2 in totest.TF2_POOL_NO_RLO:
            assert totest.infer_interpolation_class(None, tf2) == "no_MT"
        for tf1 in totest.TF1_POOL_STABLE:
            assert totest.infer_interpolation_class(tf1, "Test") == "stable_MT"
            assert totest.infer_interpolation_class(tf1, "case_A1/A2") ==\
                   "stable_reverse_MT"
        for tf1 in totest.TF1_POOL_UNSTABLE:
            assert totest.infer_interpolation_class(tf1, None) == "unstable_MT"

    def test_get_detected_initial_RLO(self, grid, iniRLO_boundary):
        # missing argument
        with raises(TypeError, match="missing 1 required positional "\
                                     +"argument: 'grid'"):
            totest.get_detected_initial_RLO()
        # bad input
        with raises(AttributeError, match="'NoneType' object has no "\
                                          +"attribute 'initial_values'"):
            totest.get_detected_initial_RLO(None)
        # examples: boundary at higher values of diagonal plane
        assert totest.get_detected_initial_RLO(grid) == iniRLO_boundary

    def test_get_nearest_known_initial_RLO(self, iniRLO_boundary, monkeypatch):
        # missing argument
        with raises(TypeError, match="missing 3 required positional "\
                                     +"arguments: 'mass1', 'mass2', and "\
                                     +"'known_initial_RLO'"):
            totest.get_nearest_known_initial_RLO()
        # bad input
        with raises(TypeError, match="object of type 'NoneType' has no len()"):
            totest.get_nearest_known_initial_RLO(None, None, None)
        # examples: don't apply the boundary
        with warns(InappropriateValueWarning, match="Don't apply initial RLO "\
                                                    +"boundary because of "\
                                                    +"too few data points in "\
                                                    +"there."):
            assert totest.get_nearest_known_initial_RLO(0.0, 0.0,\
                    iniRLO_boundary) == {"star_1_mass": 0.0,\
                                         "star_2_mass": 0.0,\
                                         "period_days": 0.0,}
        # examples: apply the boundary
        monkeypatch.setattr(totest, "MIN_COUNT_INITIAL_RLO_BOUNDARY", 0)
        for m1 in range(5):
            for m2 in range(5):
                porb = 2 - m1 - m2
                tf = str(m1) + str(m2) + str(porb)
                if porb>=0:
                    assert totest.get_nearest_known_initial_RLO(m1, m2,\
                            iniRLO_boundary) == {'star_1_mass': m1,\
                                                 'star_2_mass': m2,\
                                                 'period_days': porb,\
                                                 'termination_flag_3': tf,\
                                                 'termination_flag_4': tf}
