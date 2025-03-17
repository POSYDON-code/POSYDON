"""Unit tests of posydon/popsyn/sample_from_file.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

# import the module which will be tested
import posydon.popsyn.sample_from_file as totest
# aliases
os = totest.os
np = totest.np
pd = totest.pd

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns, approx
from inspect import isroutine, isclass
from posydon.popsyn.independent_sample import (generate_orbital_periods,
                                               generate_orbital_separations,
                                               generate_eccentricities,
                                               generate_primary_masses,
                                               generate_secondary_masses)
from posydon.utils.posydonwarning import Pwarn


# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['infer_key', 'get_samples_from_file', \
                    'get_kick_samples_from_file',  '__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__']
        assert dir(totest) == elements, "There might be added or removed "\
                                        + "objects without an update on the "\
                                        + "unit test."

    def test_instance_infer_key(self):
        assert isroutine(totest.infer_key)

    def test_instance_get_samples_from_file(self):
        assert isroutine(totest.get_samples_from_file)

    def test_instance_get_kick_samples_from_file(self):
        assert isroutine(totest.get_kick_samples_from_file)

class TestFunctions:
    
    def test_infer_key():
        
        
    def test_get_samples_from_file():
        
        
    def test_get_kick_samples_from_file():
