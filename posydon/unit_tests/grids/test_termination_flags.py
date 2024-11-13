"""Unit tests of posydon/grids/termination_flags.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.grids.termination_flags as totest

# import other needed code for the tests, which is not already imported in the
# module you like to test
from inspect import isroutine

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
        assert dir(totest) == elements, "There might be added or removed "+\
               "objects without an update on the unit test."

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
