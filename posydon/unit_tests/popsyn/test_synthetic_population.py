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
from posydon.utils.constants import Zsun
from posydon.popsyn.io import binarypop_kwargs_from_ini
from posydon.popsyn.normalized_pop_mass import initial_total_underlying_mass
from posydon.utils.common_functions import convert_metallicity_to_string
from posydon.utils.posydonwarning import Pwarn
from astropy.cosmology import Planck15 as cosmology
from astropy import constants as const

from posydon.popsyn.rate_calculation import (
    get_shell_comoving_volume,
    get_comoving_distance_from_redshift,
    get_cosmic_time_from_redshift,
    redshift_from_cosmic_time_interpolator,
    DEFAULT_SFH_MODEL,
    get_redshift_bin_edges,
    get_redshift_bin_centers,
)

from posydon.popsyn.star_formation_history import SFR_per_met_at_z

from posydon.popsyn.binarypopulation import (
    BinaryPopulation,
    HISTORY_MIN_ITEMSIZE,
    ONELINE_MIN_ITEMSIZE,
)

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
    def simple_ini(self,tmp_path):
        file_path = os.path.join(tmp_path, "test.ini")
        with open(file_path, "w") as f:
            f.write("[section]\nkey=value\n")
        return file_path
    
    @fixture
    def multi_ini(self,tmp_path):
        file1 = os.path.join(tmp_path, "a.ini")
        file2 = os.path.join(tmp_path, "b.ini")
        with open(file1, "w") as f:
            f.write("[section]\nkey1=value1\n")
        with open(file2, "w") as f:
            f.write("[section]\nkey2=value2\n")
        return [file1, file2]

    @fixture
    def textfile(self,tmp_path):
        file_path = os.path.join(tmp_path, "textfile.txt")
        with open(file_path, "w") as f:
            f.write("test")
        return file_path
    
    @fixture
    def sim_ini(self,tmp_path):
        ini_content = """
        [flow]
        import = ['posydon.binary_evol.flow_chart', 'flow_chart']
        absolute_import = None

        [step_HMS_HMS]
        import = ['posydon.binary_evol.MESA.step_mesa', 'MS_MS_step']
        absolute_import = None
        interpolation_method = 'linear3c_kNN'
        save_initial_conditions = True
        verbose = False

        [extra_hooks]
        import_1 = ['posydon.binary_evol.simulationproperties', 'TimingHooks']
        absolute_import_1 = None
        kwargs_1 = {}
        import_2 = ['posydon.binary_evol.simulationproperties', 'StepNamesHooks']
        absolute_import_2 = None
        kwargs_2 = {}
        """
        file_path = os.path.join(tmp_path, "sim.ini")
        with open(file_path, "w") as f:
            f.write(ini_content)
        return file_path
    
    @fixture
    def binpop_ini_stars(self, tmp_path):
        ini_content = """
        [BinaryPopulation_options]
        use_MPI = False
        metallicity = [0.02]
        number_of_binaries = 1
        temp_directory = 'tmp'

        [BinaryStar_output]
        extra_columns = {}
        only_select_columns = []
        scalar_names = []

        [SingleStar_1_output]
        include_S1 = True
        only_select_columns = [
            'state',
            'mass',
            'log_R']
            
        [SingleStar_2_output]
        include_S2 = True
        only_select_columns = [
            'log_L',
            'lg_mdot']

        [flow]
        import = ['builtins', 'int']

        [extra_hooks]
        import_1 = ['builtins', 'int']
        absolute_import_1 = None
        kwargs_1 = {}
        """
        file_path = os.path.join(tmp_path, "binpop_stars.ini")
        with open(file_path, "w") as f:
            f.write(ini_content)
        return file_path
    

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
        