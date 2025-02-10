"""Unit tests of posydon/popsyn/normalized_pop_mass.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

# import the module which will be tested
import posydon.popsyn.normalized_pop_mass as totest
# aliases
np = totest.np
independent_sample = totest.independent_sample
quad = totest.quad
Pwarn = totest.Pwarn

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns, approx
from inspect import isroutine, isclass

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['initial_total_underlying_mass', '__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__']
        assert dir(totest) == elements, "There might be added or removed "\
                                        + "objects without an update on the "\
                                        + "unit test."

    def test_instance_initial_total_underlying_mass(self):
        assert isroutine(totest.initial_total_underlying_mass)
        
