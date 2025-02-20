"""Unit tests of posydon/popsyn/independent_sample.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

# import the module which will be tested
import posydon.popsyn.independent_sample as totest
# aliases
np = totest.np

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns, approx
from inspect import isroutine, isclass

# define test classes collecting several test functions
class TestElements:
    
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['generate_independent_samples', 'generate_orbital_periods', \
                    'generate_orbital_separations', 'generate_eccentricities',\
                    'generate_primary_masses','generate_secondary_masses',\
                    'binary_fraction_value','__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__']
        assert dir(totest) == elements, "There might be added or removed "\
                                        + "objects without an update on the "\
                                        + "unit test."

    def test_instance_generate_independent_samples(self):
        assert isroutine(totest.generate_independent_samples)

    def test_instance_generate_orbital_periods(self):
        assert isroutine(totest.generate_orbital_periods)

    def test_instance_generate_orbital_separations(self):
        assert isroutine(totest.generate_orbital_separations)

    def test_instance_generate_eccentricities(self):
        assert isroutine(totest.generate_eccentricities)

    def test_instance_generate_primary_masses(self):
        assert isroutine(totest.generate_primary_masses)

    def test_instance_generate_secondary_masses(self):
        assert isroutine(totest.generate_secondary_masses)
        
    def test_instance_binary_fraction_value(self):
        assert isroutine(totest.binary_fraction_value)

class TestFunctions:
    
    # test functions
    def test_generate_independent_samples(self):
            
    def test_generate_orbital_periods(self):
        
    def test_generate_orbital_separations(self):
        
    def test_generate_eccentricities(self):
        
    def test_generate_primary_masses(self):
        
    def test_generate_secondary_masses(self):
        
    def test_binary_fraction_value(self):

        