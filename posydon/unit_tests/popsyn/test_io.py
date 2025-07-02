"""Unit tests of posydon/popsyn/io.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

# import the module which will be tested
import posydon.popsyn.io as totest
# aliases
np = totest.np
pd = totest.pd

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns, approx
from inspect import isroutine, isclass
from posydon.binary_evol.simulationproperties import SimulationProperties

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['BINARYPROPERTIES_DTYPES', 'OBJECT_FIXED_SUB_DTYPES', \
                    'STARPROPERTIES_DTYPES', 'EXTRA_BINARY_COLUMNS_DTYPES', \
                    'EXTRA_STAR_COLUMNS_DTYPES', 'SCALAR_NAMES_DTYPES', \
                    'clean_binary_history_df', 'clean_binary_oneline_df', \
                    'parse_inifile', 'simprop_kwargs_from_ini', \
                    'binarypop_kwargs_from_ini', 'create_run_script_text', \
                    'create_merge_script_text', \
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__', \
                    'ConfigParser', 'ast', 'importlib', 'os', 'errno', \
                    'pprint', 'np', 'pd', 'SimulationProperties']
        assert dir(totest) == elements, "There might be added or removed "\
                                        + "objects without an update on the "\
                                        + "unit test."

    def test_instance_GRB_selection(self):
        assert isroutine(totest.GRB_selection)

    def test_instance_chi_eff(self):
        assert isroutine(totest.chi_eff)

    def test_instance_m_chirp(self):
        assert isroutine(totest.m_chirp)

    def test_instance_mass_ratio(self):
        assert isroutine(totest.mass_ratio)

    def test_instance_BBH_selection_function(self):
        assert isroutine(totest.BBH_selection_function)

    def test_instance_DCO_detectability(self):
        assert isroutine(totest.DCO_detectability)

class TestFunctions:
    
    def test_GRB_selection():
        
        
    def test_chi_eff():
        
        
    def test_m_chirp():
        
        
    def test_mass_ratio():
        
    def test_BBH_selection_function():
        
    def test_DCO_detectability():
        
        