"""Unit tests of posydon/popsyn/synthetic_population.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

# import the module which will be tested
import posydon.popsyn.synthetic_population as totest
# aliases
np = totest.np
pd = totest.pd

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns, approx
from inspect import isroutine, isclass
import warnings
warnings.simplefilter("always")

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['parameter_array', 'DFInterface','History','Oneline',
                    'Population','PopulationIO','PopulationRunner',
                    'Rates','TransientPopulation',
                    '__authors__','__builtins__', '__cached__', '__doc__', 
                    '__file__','__loader__', '__name__', '__package__', '__spec__', 
                    'np', 'pd', 'tqdm', 'os', 'plt', 
                    'Zsun', 'binarypop_kwargs_from_ini',
                    'initial_total_underlying_mass','plot_pop',
                    'convert_metallicity_to_string','Pwarn','cosmology','const',
                    'get_shell_comoving_volume', 'get_comoving_distance_from_redshift',
                    'get_cosmic_time_from_redshift', 'redshift_from_cosmic_time_interpolator',
                    'DEFAULT_SFH_MODEL', 'get_redshift_bin_edges',
                    'get_redshift_bin_centers', 'SFR_per_met_at_z',
                    'BinaryPopulation', 'HISTORY_MIN_ITEMSIZE','ONELINE_MIN_ITEMSIZE'
                   ]
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

class TestPopulationRunner:
    
    @fixture
    def fix(self):
#         return 
        pass

    def test_evolve(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_merge_parallel_runs(self):
        # missing argument
        # bad input
        # examples
        pass
        
class TestDFInterface:
    
    @fixture
    def fix(self):
#         return 
        pass

    def test_head(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_tail(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_select(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_get_repr(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_get_html_repr(self):
        # missing argument
        # bad input
        # examples
        pass
        
class TestHistory:
    
    @fixture
    def fix(self):
#         return 
        pass

    def test_head(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_tail(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_select(self):
        # missing argument
        # bad input
        # examples
        pass
        
class TestOneline:
    
    @fixture
    def fix(self):
#         return 
        pass

        pass
    def test_head(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_tail(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_select(self):
        # missing argument
        # bad input
        # examples
        pass
        
class TestPopulation:
    
    @fixture
    def fix(self):
#         return 

        pass
    def test_calculate_underlying_mass(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_export_selection(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_formation_channels(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_calculate_formation_channels(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_columns(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_create_transient_population(self):
        # missing argument
        # bad input
        # examples
        pass
        
class TestTransientPopulation:
    
    @fixture
    def fix(self):
#         return 
        pass

    def test_population(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_columns(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_select(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_get_efficiency_over_metallicity(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_calculate_cosmic_weights(self):
        # missing argument
        # bad input
        # examples
        pass
        
class TestRates:
    
    @fixture
    def fix(self):
        # return 
        pass

    def test_weights(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_z_birth(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_z_events(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_select_rate_slice(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_calculate_intrinsic_rate_density(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_calculate_observable_population(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_observable_population(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_observable_population_names(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_intrinsic_rate_density(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_edges_metallicity_bins(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_centers_metallicity_bins(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_edges_redshift_bins(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_centers_redshift_bins(self):
        # missing argument
        # bad input
        # examples    
        pass
        