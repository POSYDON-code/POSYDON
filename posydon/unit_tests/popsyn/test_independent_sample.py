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

from pytest import approx, raises


# define test classes collecting several test functions
class TestElements:

    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['generate_independent_samples', 'use_Moe_17_PsandQs', \
                    '_gen_Moe_17_PsandQs','generate_orbital_periods', \
                    'generate_orbital_separations', 'generate_eccentricities',\
                    'generate_primary_masses','generate_secondary_masses',\
                    'generate_binary_fraction','__authors__',\
                    'np','truncnorm','rejection_sampler',\
                    'IMFs','distributions','Moe_17_PsandQs',\
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

class TestFunctions:

    # test functions
    def test_generate_independent_samples(self):
        # bad input
        with raises(ValueError, match="Allowed orbital schemes are separation or period."):
            totest.generate_independent_samples('test')

        # separation scheme
        orb, ecc, m1, m2 = totest.generate_independent_samples(
            orbital_scheme='separation',
            RNG=np.random.default_rng(seed=42))
        assert orb[0] == approx(24650.481799781122,abs=6e-12)
        assert ecc[0] == approx(0.8350856417514098,abs=6e-12)
        assert m1[0] == approx(19.97764511120556,abs=6e-12)
        assert m2[0] == approx(8.964150262412895,abs=6e-12)
        assert isinstance(orb, np.ndarray)
        assert len(orb) == 1
        assert all(np.isfinite(m1))

        # period scheme (default)
        orb_p, ecc_p, m1_p, m2_p = totest.generate_independent_samples(
            orbital_scheme='period',
            RNG=np.random.default_rng(seed=42))
        assert orb_p[0] == approx(872.213878458193,abs=6e-12)
        assert ecc_p[0] == approx(0.7259611833901314,abs=6e-12)
        assert m1_p[0] == approx(19.97764511120556,abs=6e-12)
        assert m2_p[0] == approx(8.964150262412895,abs=6e-12)
        assert isinstance(orb_p, np.ndarray)
        assert len(orb_p) == 1
        assert all(np.isfinite(m1_p))

    def test_use_Moe_17_PsandQs(self):

        # returns True for Moe+17-PsandQs secondary_mass_scheme
        assert totest.use_Moe_17_PsandQs(secondary_mass_scheme='Moe+17-PsandQs') is True

        # returns True for Moe+17-PsandQs orbital_period_scheme with period scheme
        assert totest.use_Moe_17_PsandQs(
            orbital_scheme='period',
            orbital_period_scheme='Moe+17-PsandQs') is True

        # returns True for Moe+17-PsandQs eccentricity_scheme
        assert totest.use_Moe_17_PsandQs(eccentricity_scheme='Moe+17-PsandQs') is True
        # returns False for non-Moe schemes
        assert totest.use_Moe_17_PsandQs(
            secondary_mass_scheme='flat_mass_ratio',
            orbital_scheme='period',
            orbital_period_scheme='Sana+12_period_extended',
            eccentricity_scheme='zero') is False

    def test_generate_orbital_periods(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional argument: 'primary_masses'"):
            totest.generate_orbital_periods()

        # bad input
        with raises(ValueError, match="p_max must be greater than p_min"):
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
        # missing log_normal params
        with raises(ValueError, match="For the `log_normal separation` scheme you must give"):
            totest.generate_orbital_separations(orbital_separation_scheme='log_normal')

        # bad input: min > max (raised by LogUniform distribution class)
        with raises(ValueError, match="max must be greater than min"):
            totest.generate_orbital_separations(orbital_separation_min=10.,
                                                orbital_separation_max=1.)

        # bad input: min > max with log_normal
        with raises(ValueError, match="`orbital_separation_max` must be"):
            totest.generate_orbital_separations(
                orbital_separation_scheme='log_normal',
                log_orbital_separation_mean=1.0,
                log_orbital_separation_sigma=1.0,
                orbital_separation_min=10.,
                orbital_separation_max=1.)

        # bad scheme
        with raises(ValueError, match="You must provide an allowed orbital separation scheme."):
            totest.generate_orbital_separations(orbital_separation_scheme='test')

        # log_normal examples
        tests_normal = [(0., 1.0, 42, approx(39.83711402835139, abs=6e-12)),
                        (1.0, 10., 42, approx(9799.179319004, abs=6e-9))]
        for (m, s, r, sep) in tests_normal:
            assert totest.generate_orbital_separations(
                orbital_separation_scheme='log_normal',
                log_orbital_separation_mean=m,
                log_orbital_separation_sigma=s,
                RNG=np.random.default_rng(seed=r))[0] == sep

        # log_uniform examples
        tests_uniform = [(1., 3., 42, approx(2.3402964885050066, abs=6e-12)),
                         (2., 10., 42, approx(6.950276115688688, abs=6e-12))]
        for (mi, ma, r, sep) in tests_uniform:
            assert totest.generate_orbital_separations(
                orbital_separation_min=mi,
                orbital_separation_max=ma,
                RNG=np.random.default_rng(seed=r))[0] == sep
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
        # bad input: invalid scheme
        with raises(ValueError, match="You must provide an allowed primary mass scheme."):
            totest.generate_primary_masses(primary_mass_scheme='test')

        # bad input: min > max (raised by IMF class)
        with raises(ValueError, match="m_min must be less than m_max"):
            totest.generate_primary_masses(primary_mass_min=100., primary_mass_max=10.)

        # examples for all three schemes
        tests = [('Salpeter', 42, approx(19.97764511120556, abs=6e-12)),
                 ('Kroupa1993', 42, approx(16.52331793661949, abs=6e-12)),
                 ('Kroupa2001', 42, approx(20.633204764212334, abs=6e-12))]
        for (s, r, m1) in tests:
            assert totest.generate_primary_masses(
                primary_mass_scheme=s,
                RNG=np.random.default_rng(seed=r))[0] == m1

    def test_generate_secondary_masses(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional argument: 'primary_masses'"):
            totest.generate_secondary_masses()

        # bad input: invalid scheme
        with raises(ValueError, match="You must provide an allowed secondary mass scheme."):
            totest.generate_secondary_masses(primary_masses=np.array([10.]),
                                             secondary_mass_scheme='test')

        # bad input: secondary_mass_min > primary mass
        with raises(ValueError, match="`secondary_mass_min` is larger than some primary masses"):
            totest.generate_secondary_masses(primary_masses=np.array([1.]),
                                             secondary_mass_min=10.,
                                             secondary_mass_max=100.)

        # flat_mass_ratio example
        result = totest.generate_secondary_masses(
            primary_masses=np.array([10.]),
            secondary_mass_scheme='flat_mass_ratio',
            RNG=np.random.default_rng(seed=42))
        assert len(result) == 1
        assert result[0] > 0
        assert result[0] <= 10.0

        # q=1 example
        result_q1 = totest.generate_secondary_masses(
            primary_masses=np.array([10.]),
            secondary_mass_scheme='q=1',
            RNG=np.random.default_rng(seed=42))
        assert result_q1[0] == approx(10., abs=6e-12)

    def test_generate_binary_fraction(self):
        # missing primary mass
        with raises(ValueError, match="There was not a primary mass provided in the inputs"):
            totest.generate_binary_fraction(binary_fraction_scheme='const')

        # bad scheme (m1 must be provided before scheme check)
        with raises(ValueError, match="You must provide an allowed binary fraction scheme."):
            totest.generate_binary_fraction(binary_fraction_scheme='test',
                                            m1=np.array([10.]))

        # const scheme examples
        tests_const = [1.0, 1, 0.5]
        for c in tests_const:
            assert totest.generate_binary_fraction(
                binary_fraction_const=c,
                binary_fraction_scheme='const',
                m1=np.array([10.])) == c

        # non-array m1 input (triggers np.asarray conversion)
        assert totest.generate_binary_fraction(
            binary_fraction_const=0.7,
            binary_fraction_scheme='const',
            m1=10.) == 0.7

        # Moe+17-massdependent scheme
        tests_moe = [(np.array([1.]), 0.4),
                     (np.array([3.]), 0.59),
                     (np.array([8.]), 0.76),
                     (np.array([10.]), 0.84),
                     (np.array([18.]), 0.94)]
        for (m1, f) in tests_moe:
            result = totest.generate_binary_fraction(
                binary_fraction_scheme='Moe+17-massdependent', m1=m1)
            assert result[0] == f
