"""Unit tests of posydon/popsyn/GRB.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

# import the module which will be tested
import posydon.popsyn.GRB as totest
# aliases
np = totest.np

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns, approx
from inspect import isroutine, isclass
from posydon.utils.posydonwarning import Pwarn
from posydon.utils.constants import Msun, clight

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['get_GRB_properties','__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__',
                    'Pwarn','Msun','clight','np']
        assert dir(totest) == elements, "There might be added or removed "\
                                        + "objects without an update on the "\
                                        + "unit test."

    def test_instance_get_GRB_properties(self):
        assert isroutine(totest.get_GRB_properties)
        
class TestFunctions:
    
    def test_get_GRB_properties():
       