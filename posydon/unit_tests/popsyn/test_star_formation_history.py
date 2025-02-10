"""Unit tests of posydon/popsyn/star_formation_history.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

# import the module which will be tested
import posydon.popsyn.star_formation_history as totest
# aliases
os = totest.os
np = totest.np
sp = totest.sp
stats = totest.stats

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns, approx
from inspect import isroutine, isclass
from posydon.utils.data_download import PATH_TO_POSYDON_DATA
from posydon.utils.constants import age_of_universe
from posydon.utils.common_functions import (
    rejection_sampler,
    histogram_sampler,
    read_histogram_from_file,
)
from posydon.utils.constants import Zsun
from scipy.interpolate import interp1d
from astropy.cosmology import Planck15 as cosmology

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['get_formation_times', 'get_illustrisTNG_data', \
                    'star_formation_rate', 'mean_metallicity',\
                    'std_log_metallicity_dist','SFR_Z_fraction_at_given_redshift',\
                    'integrated_SFRH_over_redshift','__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__']
        assert dir(totest) == elements, "There might be added or removed "\
                                        + "objects without an update on the "\
                                        + "unit test."

    def test_instance_get_formation_times(self):
        assert isroutine(totest.GRB_selection)

    def test_instance_get_illustrisTNG_data(self):
        assert isroutine(totest.chi_eff)

    def test_instance_star_formation_rate(self):
        assert isroutine(totest.m_chirp)

    def test_instance_mean_metallicity(self):
        assert isroutine(totest.mass_ratio)

    def test_instance_std_log_metallicity_dist(self):
        assert isroutine(totest.BBH_selection_function)

    def test_instance_SFR_Z_fraction_at_given_redshift(self):
        assert isroutine(totest.DCO_detectability)

    def test_instance_integrated_SFRH_over_redshift(self):
        assert isroutine(totest.DCO_detectability)

class TestFunctions:
    
    def test_get_formation_times():
        
        
    def test_get_illustrisTNG_data():
        
        
    def test_star_formation_rate():
        
        
    def test_mass_mean_metallicity():
        
    def test_BBH_std_log_metallicity_dist():
        
    def test_SFR_Z_fraction_at_given_redshift():
        
    def test_integrated_SFRH_over_redshift():
        
        