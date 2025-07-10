"""Unit tests of posydon/popsyn/transient_select_funcs.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

# import the module which will be tested
import posydon.popsyn.transient_select_funcs as totest
# aliases
np = totest.np
pd = totest.pd
PATH_TO_POSYDON_DATA = totest.PATH_TO_POSYDON_DATA

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns, approx
from inspect import isroutine, isclass
from posydon.utils.posydonwarning import ReplaceValueWarning
import warnings
warnings.simplefilter("always")

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['PATH_TO_PDET_GRID', 'GRB_selection', 'chi_eff', 'm_chirp', \
                    'mass_ratio', 'BBH_selection_function','DCO_detectability', \
                    '__builtins__', '__cached__', '__doc__', \
                    '__file__','__loader__', '__name__', '__package__', '__spec__', \
                    'np', 'pd', 'PATH_TO_POSYDON_DATA', \
                    'os', 'tqdm', 'warnings', 'Pwarn',]
        totest_elements = set(dir(totest))
        missing_in_test = set(elements) - totest_elements
        assert len(missing_in_test) == 0, "There are missing objects in "\
                                          +f"{totest.__name__}: "\
                                          +f"{missing_in_test}. Please "\
                                          +"check, whether they have been "\
                                          +"removed on purpose and update "\
                                          +"this unit test."
        new_in_test = totest_elements - set(elements)
        assert len(new_in_test) == 0, "There are new objects in "\
                                      +f"{totest.__name__}: {new_in_test}. "\
                                      +"Please check, whether they have been "\
                                      +"added on purpose and update this "\
                                      +"unit test."

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
    
    @fixture
    def array(self):
        return np.array([1.0,2.0,3.0])
    
    @fixture
    def nan_array(self):
        return np.array([np.nan,np.nan,np.nan])
    
    @fixture
    def wrong_array(self):
        return np.array(['1.0','2.0','3.0'])
    
    def test_GRB_selection(self):
        # bad input
        # examples
        pass 
    
    def test_chi_eff(self,array,nan_array,wrong_array):
        # missing argument
        with raises(TypeError,match="missing 6 required positional arguments"):
            totest.chi_eff()
        # bad input
        with raises(TypeError,match="ufunc 'cos' not supported for the input types"):
            totest.chi_eff(array,array,array,array,array,wrong_array)
        # undefined values
        with warns(ReplaceValueWarning,match="a_1 contains undefined values"):
            totest.chi_eff(array,array,nan_array.copy(),array,array,array)
        with warns(ReplaceValueWarning,match="a_2 contains undefined values"):
            totest.chi_eff(array,array,array,nan_array.copy(),array,array)
        with warns(ReplaceValueWarning,match="tilt_1 contains undefined values"):
            totest.chi_eff(array,array,array,array,nan_array.copy(),array)
        with warns(ReplaceValueWarning,match="tilt_2 contains undefined values"):
            totest.chi_eff(array,array,array,array,array,nan_array.copy())
        # example
        assert totest.chi_eff(array,array,array,
                              array,array,array)[0] == 0.5403023058681398
        
    def test_m_chirp(self):
        # missing argument
        with raises(TypeError,match="missing 1 required positional argument: 'm_2'"):
            totest.m_chirp(3.)
        # bad input
        with raises(TypeError,match="can't multiply sequence by non-int of type 'str'"): 
            totest.m_chirp("3.","2.")
        # examples
        tests = [(4.,2.,2.433457367572823),
                 (40.,10.,16.65106414803746)]
        for (m1,m2,mc) in tests:
            assert totest.m_chirp(m1,m2) == mc
            
    def test_mass_ratio(self):
        # missing argument
        with raises(TypeError,match="missing 1 required positional argument: 'm_2'"):
            totest.mass_ratio(3.)
        # bad input
        with raises(TypeError,match="unsupported operand type"): 
            totest.mass_ratio("3.","2.")
        # examples
        tests = [(5.,1.,0.2),
                 (1.,5.,0.2),
                 (4.,2.,0.5)]
        for (m1,m2,q) in tests:
            assert totest.mass_ratio(np.array([m1]),
                                     np.array([m2])) == q
            
    def test_BBH_selection_function(self):
        # bad input
        # examples
        pass
    def test_DCO_detectability(self):
        # bad input
        # examples
        pass
        