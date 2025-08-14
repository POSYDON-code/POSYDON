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
import re
from pytest import raises, approx
from inspect import isroutine, isclass

# define test classes collecting several test functions
class TestElements:
    
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['generate_independent_samples', 'generate_orbital_periods', \
                    'generate_orbital_separations', 'generate_eccentricities',\
                    'generate_primary_masses','generate_secondary_masses',\
                    'binary_fraction_value','__authors__',\
                    'np','truncnorm','rejection_sampler',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__']
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
        # bad input
        with raises(ValueError, match="Allowed orbital schemes are separation or period."):
            totest.generate_independent_samples('test')
        # examples
        tests = [("separation",42,approx(4993.106338349307,abs=6e-12)),
                 ("period",12,approx(200.82071793763188,abs=6e-12))]
        for (s,r,o) in tests:
            orb,ecc,m1,m2 = totest.generate_independent_samples(orbital_scheme=s,
                                                                RNG = np.random.default_rng(seed=42))
            assert orb[0] == approx(o,abs=6e-12)
            assert ecc[0] == approx(0.8797477186989253,abs=6e-12)
            assert m1[0] == approx(10.607132832170066,abs=6e-12)
            assert m2[0] == approx(9.182255718237206,abs=6e-12)
            
    def test_generate_orbital_periods(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional argument: 'primary_masses'"):
            totest.generate_orbital_periods()
        # bad input
        with raises(TypeError, match="expected a sequence of integers or a single integer"):
            totest.generate_orbital_periods(np.array([1.]),
                                            number_of_binaries=1.)
        with raises(ValueError, match="high - low < 0"):
            totest.generate_orbital_periods(np.array([1.]),
                                            orbital_period_min=10.,
                                            orbital_period_max=1.)
        with raises(ValueError, match="You must provide an allowed orbital period scheme."):
            totest.generate_orbital_periods(np.array([1.]),
                                            orbital_period_scheme='test')
        # examples
        tests = [(1.0,42,approx(403.44608837021764,abs=6e-12)),
                 (1.0,12,approx(3.4380527315000666,abs=6e-12))]
        for (m,r,p) in tests:
            assert totest.generate_orbital_periods(m,RNG = np.random.default_rng(seed=r))[0] == p
        
    def test_generate_orbital_separations(self):
        # missing argument
        with raises(ValueError,match="For the `log_normal separation` scheme you must give `log_orbital_separation_mean`, `log_orbital_separation_sigma`."):
            totest.generate_orbital_separations(orbital_separation_scheme='log_normal')
        # bad input
        with raises(TypeError, match="expected a sequence of integers or a single integer"):
            totest.generate_orbital_separations(number_of_binaries=1.)
        with raises(ValueError, match="high - low < 0"):
            totest.generate_orbital_separations(orbital_separation_min=10.,
                                                orbital_separation_max=1.)
        with raises(OverflowError, match="high - low range exceeds valid bounds"):
            totest.generate_orbital_separations(orbital_separation_min=0)
        with raises(ValueError, match="You must provide an allowed orbital separation scheme."):
            totest.generate_orbital_separations(orbital_separation_scheme='test')
        # examples
        tests_normal = [(0.,1.0,42,approx(39.83711402835139,abs=6e-12)),
                        (1.0,10.,42,approx(9799.179319004,abs=6e-9))]
        # larger allowance for second test because of slightly different results between 
        # running pytest locally vs github actions workflow
        for (m,s,r,sep) in tests_normal:
            assert totest.generate_orbital_separations(orbital_separation_scheme='log_normal',
                                                log_orbital_separation_mean=m,
                                                log_orbital_separation_sigma=s,
                                                RNG = np.random.default_rng(seed=r))[0] == sep
        tests_uniform = [(1.,3.,42,approx(2.3402964885050066,abs=6e-12)),
                         (2.,10.,42,approx(6.950276115688688,abs=6e-12))]
        for (mi,ma,r,sep) in tests_uniform:
            assert totest.generate_orbital_separations(orbital_separation_min=mi,
                                                orbital_separation_max=ma,
                                                RNG = np.random.default_rng(seed=r))[0] == sep
    def test_generate_eccentricities(self):
        # bad input
        with raises(TypeError, match="expected a sequence of integers or a single integer"):
            totest.generate_eccentricities(number_of_binaries=1.)
        with raises(ValueError, match="You must provide an allowed eccentricity scheme."):
            totest.generate_eccentricities(eccentricity_scheme='test')
        # examples
        tests = [('thermal',42,approx(0.8797477186989253,abs=6e-12)),
                 ('uniform',42,approx(0.7739560485559633,abs=6e-12)),
                 ('zero',42,approx(0.,abs=6e-12))]
        for (s,r,e) in tests:
            assert totest.generate_eccentricities(eccentricity_scheme=s,
                                                   RNG = np.random.default_rng(seed=r))[0] == e
        
    def test_generate_primary_masses(self):
        # bad input
        with raises(TypeError, match="expected a sequence of integers or a single integer"):
            totest.generate_primary_masses(number_of_binaries=1.)
        with raises(ValueError, match="primary_mass_max must be larger than primary_mass_min."):
            totest.generate_primary_masses(primary_mass_min=100.,primary_mass_max=10.)
        with raises(ValueError, match="You must provide an allowed primary mass scheme."):
            totest.generate_primary_masses(primary_mass_scheme='test')
        # examples
        tests = [('Salpeter',42,approx(19.97764511120556,abs=6e-12)),
                 ('Kroupa1993',42,approx(16.52331793661949,abs=6e-12)),
                 ('Kroupa2001',42,approx(20.633412780370865,abs=6e-12))]
        for (s,r,m1) in tests:
            assert totest.generate_primary_masses(primary_mass_scheme=s,
                                                   RNG = np.random.default_rng(seed=r))[0] == m1
            
    def test_generate_secondary_masses(self):
        # missing argument
        with raises(TypeError,match="missing 1 required positional argument: 'primary_masses'"):
            totest.generate_secondary_masses()
        # bad input
        with raises(TypeError, match=re.escape("unsupported operand type(s) for /: 'float' and 'list'")):
            totest.generate_secondary_masses(primary_masses=[10.])
        with raises(TypeError, match="expected a sequence of integers or a single integer"):
            totest.generate_secondary_masses(primary_masses=np.array([10.]),
                                             number_of_binaries=1.)
        with raises(ValueError, match="secondary_mass_max must be larger than secondary_mass_min"):
            totest.generate_secondary_masses(primary_masses=np.array([100.]),
                                              secondary_mass_min=10.,
                                              secondary_mass_max=1.)
        with raises(ValueError, match="`secondary_mass_min` is larger than some primary masses"):
            totest.generate_secondary_masses(primary_masses=np.array([1.]),
                                              secondary_mass_min=10.,
                                              secondary_mass_max=100.)
        with raises(ValueError, match="You must provide an allowed secondary mass scheme."):
            totest.generate_secondary_masses(primary_masses=np.array([1.]),
                                             secondary_mass_scheme='test')
        # examples
        tests = [('flat_mass_ratio',42,approx(7.852582461281652,abs=6e-12)),
                 ('q=1',42,approx(10.,abs=6e-12))]
        for (s,r,m2) in tests:
            assert totest.generate_secondary_masses(primary_masses=np.array([10.]),
                                                    secondary_mass_scheme=s,
                                                   RNG = np.random.default_rng(seed=r))[0] == m2
            
    def test_binary_fraction_value(self):
        # missing argument
        with raises(ValueError, match="There was not a primary mass provided in the inputs."):
            totest.binary_fraction_value(binary_fraction_scheme='Moe_17')
        # bad input
        with raises(ValueError, match="You must provide an allowed binary fraction scheme."):
            totest.binary_fraction_value(binary_fraction_scheme='test')
        with raises(ValueError, match="The scheme doesn't support values of m1 less than 0.8"):
            totest.binary_fraction_value(binary_fraction_scheme='Moe_17',m1=0.2)
        with raises(ValueError, match="The primary mass provided nan is not supported by the Moe_17 scheme."):
            totest.binary_fraction_value(binary_fraction_scheme='Moe_17',m1=np.nan)
        # examples
        tests_const = [1.0,1,0.5]
        for (c) in tests_const:
            assert totest.binary_fraction_value(binary_fraction_const=c,
                                                binary_fraction_scheme='const') == c
        tests_moe = [(1,0.4),
                     (3,0.59),
                     (8,0.76),
                     (10,0.84),
                     (18,0.94)]
        for (m1,f) in tests_moe:
            assert totest.binary_fraction_value(binary_fraction_scheme='Moe_17',
                                                m1=m1) == f
        
        
        
        

        