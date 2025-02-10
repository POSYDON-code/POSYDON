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


# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['BinaryStar', 'CC_quantities',\
                    'CEE_parameters_from_core_abundance_thresholds',\
                    'Catch_POSYDON_Warnings',\
                    'DEFAULT_MARKERS_COLORS_LEGENDS', 'MODELS', 'Pwarn',\
                    'STAR_STATES_CC', 'SingleStar', 'StepSN',\
                    'TF1_POOL_STABLE', '__authors__', '__builtins__',\
                    '__cached__', '__credits__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__',\
                    'add_post_processed_quantities',\
                    'assign_core_collapse_quantities_none',\
                    'calculate_Patton20_values_at_He_depl',\
                    'check_state_of_star', 'combine_TF12', 'copy', 'np',\
                    'post_process_grid', 'print_CC_quantities', 'tqdm']
        assert dir(totest) == elements, "There might be added or removed "\
                                        + "objects without an update on the "\
                                        + "unit test."

    def test_instance_CC_quantities(self):
        assert isinstance(totest.CC_quantities, (list)), "CC_quantities is "\
               + "of type: " + str(type(totest.CC_quantities))
