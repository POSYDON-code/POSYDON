"""Test the initial_distributions module for POSYDON."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>"
]

import numpy as np
import pytest
from pytest import fixture, mark
from scipy.integrate import quad

import posydon.popsyn.IMFs as IMFs

# file to test
from posydon.popsyn import distributions, initial_distributions
from posydon.utils.common_functions import orbital_period_from_separation


# Dummy IMF to simulate a valid IMF pdf
class DummyIMF(IMFs.IMFBase):
    def __init__(self, m_min, m_max):
        self.m_min = m_min
        self.m_max = m_max
    def imf(self, m):
        return 1.0
    def rvs(self, size=1, rng=None):
        raise NotImplementedError
    def pdf(self, m):
        # returns a constant value for testing
        return 0.5


class FlatIMF(IMFs.IMFBase):
    '''Test-only flat IMF, normalised on [m_min, m_max].'''
    def imf(self, m):
        return np.ones_like(np.asarray(m, dtype=float))
    def rvs(self, size=1, rng=None):
        raise NotImplementedError


@fixture(autouse=True)
def register_test_imfs(monkeypatch):
    '''Make the test-only IMFs available as primary_mass_scheme names.'''
    monkeypatch.setattr(IMFs, 'FlatIMF', FlatIMF, raising=False)
    monkeypatch.setattr(IMFs, 'DummyIMF', DummyIMF, raising=False)


class TestGetIMFDistribution:
    @mark.parametrize('scheme', ['NonExistentIMF', '', 'IMFBase', 'np',
                                 'Pwarn', 'quad'])
    def test_invalid_imf_raises(self, scheme):
        # Unknown names, the abstract base class, and module attributes that
        # are not IMF classes must not silently fall back to another IMF.
        kwargs = {
            'primary_mass_scheme': scheme,
            'primary_mass_min': 1,
            'primary_mass_max': 100
        }
        with pytest.raises(ValueError, match="primary_mass_scheme"):
            initial_distributions.get_IMF_distribution(kwargs)

    def test_missing_imf_scheme_raises(self):
        kwargs = {'primary_mass_min': 1, 'primary_mass_max': 100}
        with pytest.raises(ValueError, match="primary_mass_scheme"):
            initial_distributions.get_IMF_distribution(kwargs)

    def test_imf_constructor_error_propagates(self, monkeypatch):
        # An AttributeError raised inside a valid IMF must not be mistaken
        # for an unknown scheme.
        class BrokenIMF(IMFs.IMFBase):
            def __init__(self, m_min, m_max):
                raise AttributeError("broken IMF")
            def imf(self, m):
                return 1.0
            def rvs(self, size=1, rng=None):
                raise NotImplementedError

        monkeypatch.setattr(IMFs, 'BrokenIMF', BrokenIMF, raising=False)
        kwargs = {
            'primary_mass_scheme': 'BrokenIMF',
            'primary_mass_min': 1,
            'primary_mass_max': 100
        }
        with pytest.raises(AttributeError, match="broken IMF"):
            initial_distributions.get_IMF_distribution(kwargs)

    @mark.parametrize('scheme', ['Salpeter', 'Kroupa1993', 'Kroupa2001'])
    def test_sampler_imfs_supported(self, scheme):
        kwargs = {
            'primary_mass_scheme': scheme,
            'primary_mass_min': 7,
            'primary_mass_max': 140
        }
        pdf_func = initial_distributions.get_IMF_distribution(kwargs).pdf
        assert np.isclose(quad(pdf_func, 7, 140)[0], 1.0)


    def test_valid_imf_returns_dummy_pdf(self, monkeypatch):
        # Use monkeypatch to add DummyIMF to IMFs module
        monkeypatch.setattr(IMFs, 'DummyIMF', DummyIMF, raising=False)

        kwargs = {
            'primary_mass_scheme': 'DummyIMF',
            'primary_mass_min': 1,
            'primary_mass_max': 100
        }
        pdf_func = initial_distributions.get_IMF_distribution(kwargs).pdf
        m_test = np.array([1, 10, 50])
        result = pdf_func(m_test)
        # Expect DummyIMF.pdf to return 0.5 regardless of m1
        expected = 0.5
        assert np.allclose(result, expected)

    def test_imf_with_scheme_kwargs(self, monkeypatch):
        # Create a test IMF class that uses scheme kwargs
        class SchemeKwargsIMF(IMFs.IMFBase):
            def __init__(self, m_min, m_max, custom_param=None):
                self.m_min = m_min
                self.m_max = m_max
                self.custom_param = custom_param
            def imf(self, m):
                return 1.0
            def rvs(self, size=1, rng=None):
                raise NotImplementedError
            def pdf(self, m):
                # Return custom_param if set, otherwise 1
                return self.custom_param if self.custom_param is not None else 1

        # Use monkeypatch to add SchemeKwargsIMF to IMFs module
        monkeypatch.setattr(IMFs,
                            'SchemeKwargsIMF',
                            SchemeKwargsIMF,
                            raising=False)

        # Create kwargs with scheme-specific parameters
        kwargs = {
            'primary_mass_scheme': 'SchemeKwargsIMF',
            'primary_mass_min': 1,
            'primary_mass_max': 100,
            'SchemeKwargsIMF': {'custom_param': 3.5}  # Scheme-specific kwargs
        }

        # Get the PDF function
        pdf_func = initial_distributions.get_IMF_distribution(kwargs).pdf

        # Test that the custom parameter was used
        m_test = np.array([1, 10, 50])
        result = pdf_func(m_test)
        expected = kwargs['SchemeKwargsIMF']['custom_param']  # The value of custom_param
        assert np.allclose(result, expected)



class TestGetMassRatioDistribution:
    def test_flat_mass_ratio_pdf(self):
        # kwargs without q_min and q_max to trigger flat_mass_ratio branch
        kwargs = {
            'primary_mass_min': 5,
            'primary_mass_max': 100,
            'secondary_mass_scheme': 'flat_mass_ratio',
            'secondary_mass_min': 5,
            'secondary_mass_max': 50
        }
        q_pdf = initial_distributions.get_mass_ratio_distribution(kwargs).pdf
        # For m1 = 10, minimum = max(5/10, 0.05) = 0.5 and
        # maximum = min(50/10, 1) = 1, so pdf = 1/(1-0.5) = 2.
        m1_val = 10
        # Test with q in range
        result = q_pdf(7.5/m1_val, m1_val)
        assert np.allclose(result, 2)
        # Test with q below the minimum
        result = q_pdf(3.0/m1_val, m1_val)
        assert np.all(result == 0)
        # Test with q above the maximum
        result = q_pdf(11.0/m1_val, m1_val)
        assert np.all(result == 0)

    def test_alternative_mass_ratio_pdf(self):
        # kwargs with q_min and q_max to trigger the alternative branch
        kwargs = {
            'secondary_mass_scheme': 'flat_mass_ratio',
            'q_min': 0.2,
            'q_max': 0.8,
        }
        q_pdf = initial_distributions.get_mass_ratio_distribution(kwargs).pdf
        # In alternative branch, pdf returns 1/(q_max-q_min) for q in [q_min,q_max]
        result = q_pdf(0.5*(kwargs['q_min']+kwargs['q_max']), None)
        assert np.all(result == 1.0/(kwargs['q_max']-kwargs['q_min']))
        result = q_pdf(0.5*kwargs['q_min'], None)
        assert np.all(result == 0)
        result = q_pdf(0.5*(kwargs['q_max']+1.0), None)
        assert np.all(result == 0)

    @mark.parametrize('scheme', ['invalid', 'q=1', ''])
    def test_invalid_mass_ratio_scheme(self, scheme):
        # Unsupported schemes (including the sampler's 'q=1') must raise
        # instead of silently reweighting with a flat mass-ratio PDF.
        kwargs = {
            'secondary_mass_scheme': scheme
        }
        with pytest.raises(ValueError, match="secondary_mass_scheme"):
            initial_distributions.get_mass_ratio_distribution(kwargs)

    def test_power_law_mass_ratio_pdf(self):
        """Test that power_law_mass_ratio scheme returns correct PDF values."""
        kwargs = {
            'secondary_mass_scheme': 'power_law_mass_ratio',
            'mass_ratio_slope': 0.0,
            'q_min': 0.1,
            'q_max': 0.9,
        }
        q_pdf = initial_distributions.get_mass_ratio_distribution(kwargs).pdf
        # alpha=0 gives a flat distribution over (0.1, 0.9]
        result_in = q_pdf(0.5, None)
        assert result_in > 0
        result_out = q_pdf(0.05, None)
        assert result_out == 0

    def test_power_law_mass_ratio_pdf_default_bounds(self):
        """Test power_law_mass_ratio uses default q_min/q_max when absent."""
        kwargs = {
            'secondary_mass_scheme': 'power_law_mass_ratio',
            'mass_ratio_slope': 1.0,
        }
        q_pdf = initial_distributions.get_mass_ratio_distribution(kwargs).pdf
        # Default q_min=0.05, q_max=1.0; value inside range should be positive
        result = q_pdf(0.5, None)
        assert result > 0

    @mark.parametrize('given', [{'q_min': 0.1}, {'q_max': 0.9}])
    def test_flat_mass_ratio_single_q_bound_raises(self, given):
        # Only one fixed bound is ambiguous between the m1-dependent and
        # the fixed-bounds flat mass ratio.
        kwargs = {
            'primary_mass_min': 7,
            'primary_mass_max': 140,
            'secondary_mass_scheme': 'flat_mass_ratio',
            'secondary_mass_min': 0.35,
            'secondary_mass_max': 140,
            **given,
        }
        missing = 'q_max' if 'q_min' in given else 'q_min'
        with pytest.raises(KeyError, match=missing):
            initial_distributions.get_mass_ratio_distribution(kwargs)

    def test_primary_mass_min_below_secondary_mass_min_raises(self):
        # The sampler rejects primaries lighter than secondary_mass_min.
        kwargs = {
            'primary_mass_min': 0.1,
            'primary_mass_max': 140,
            'secondary_mass_scheme': 'flat_mass_ratio',
            'secondary_mass_min': 0.35,
            'secondary_mass_max': 140,
        }
        with pytest.raises(ValueError, match="secondary_mass_min"):
            initial_distributions.get_mass_ratio_distribution(kwargs)

    def test_flat_mass_ratio_pdf_outside_primary_range(self):
        # Evaluating at m1 outside the primary mass range (e.g. simulated
        # systems outside the requested population) gives zero, not inf/nan.
        kwargs = {
            'primary_mass_min': 7,
            'primary_mass_max': 140,
            'secondary_mass_scheme': 'flat_mass_ratio',
            'secondary_mass_min': 1,
            'secondary_mass_max': 140,
        }
        q_pdf = initial_distributions.get_mass_ratio_distribution(kwargs).pdf
        # m1 = 1 = secondary_mass_min: empty range (width 0)
        # m1 = 0.5: negative width
        result = q_pdf(np.array([0.5, 1.0]), np.array([1.0, 0.5]))
        assert np.all(result == 0)

    @mark.parametrize('kwargs, expected', [
        ({'secondary_mass_scheme': 'flat_mass_ratio',
          'primary_mass_min': 7, 'primary_mass_max': 140,
          'secondary_mass_min': 0.35, 'secondary_mass_max': 140},
         distributions.ConditionalFlatMassRatio),
        ({'secondary_mass_scheme': 'flat_mass_ratio',
          'q_min': 0.2, 'q_max': 0.8},
         distributions.FlatMassRatio),
        ({'secondary_mass_scheme': 'power_law_mass_ratio',
          'mass_ratio_slope': 1.0},
         distributions.PowerLawMassRatio),
    ])
    def test_scheme_to_class(self, kwargs, expected):
        q_dist = initial_distributions.get_mass_ratio_distribution(kwargs)
        assert type(q_dist) is expected

    def test_m1_dependent_parameters(self):
        q_dist = initial_distributions.get_mass_ratio_distribution({
            'secondary_mass_scheme': 'flat_mass_ratio',
            'primary_mass_min': 5, 'primary_mass_max': 100,
            'secondary_mass_min': 5, 'secondary_mass_max': 50})
        q_lower, q_upper = q_dist.support(np.array([5., 10., 100.]))
        assert np.allclose(q_lower, [1.0, 0.5, 0.05])
        assert np.allclose(q_upper, [1.0, 1.0, 0.5])

    @mark.parametrize('extra, q_min, q_max', [
        ({}, 0.05, 1.0),
        ({'q_min': 0.3}, 0.3, 1.0),
        ({'q_max': 0.6}, 0.05, 0.6),
    ])
    def test_power_law_default_bounds(self, extra, q_min, q_max):
        q_dist = initial_distributions.get_mass_ratio_distribution({
            'secondary_mass_scheme': 'power_law_mass_ratio',
            'mass_ratio_slope': 1.0, **extra})
        assert q_dist.support() == (q_min, q_max)

    @mark.parametrize('scheme', ['power_law_mass_ratio', 'flat_mass_ratio'])
    def test_fixed_bounds_q_min_greater_than_q_max(self, scheme):
        kwargs = {'secondary_mass_scheme': scheme, 'mass_ratio_slope': 0.0,
                  'q_min': 0.8, 'q_max': 0.2}
        with pytest.raises(ValueError, match="q_min must be less than q_max"):
            initial_distributions.get_mass_ratio_distribution(kwargs)

    @mark.parametrize('scheme', ['invalid', 'q=1', None])
    def test_unsupported_scheme(self, scheme):
        with pytest.raises(ValueError, match="secondary_mass_scheme"):
            initial_distributions.get_mass_ratio_distribution(
                {'secondary_mass_scheme': scheme})


class TestGetBinaryFractionDistribution:

    def test_const_binary_fraction_pdf(self):
        kwargs = {
            'binary_fraction_scheme': 'const',
            'binary_fraction_const': 0.7
        }
        pdf_func = initial_distributions.get_binary_fraction_distribution(kwargs).pdf
        # For binary stars, expected = 0.7
        result = pdf_func(True)
        assert np.allclose(result, np.asarray(0.7))
        # For single stars, expected = 1 - 0.7 = 0.3
        result = pdf_func(False)
        assert np.allclose(result, np.asarray(0.3))

    def test_returns_constant_binary_fraction(self):
        f_b_dist = initial_distributions.get_binary_fraction_distribution(
            {'binary_fraction_scheme': 'const', 'binary_fraction_const': 0.7})
        assert isinstance(f_b_dist, distributions.ConstantBinaryFraction)
        assert f_b_dist.f_b == 0.7

    @mark.parametrize('scheme', ['invalid', 'Moe+17-massdependent', None])
    def test_invalid_binary_fraction_scheme(self, scheme):
        kwargs = {
            'binary_fraction_scheme': scheme
        }
        with pytest.raises(ValueError) as excinfo:
            initial_distributions.get_binary_fraction_distribution(kwargs)

        assert "Binary fraction scheme not recognized" in str(excinfo.value)

class TestGetOrbitalDistribution:
    def test_sana12_period_extended(self):
        kwargs = {
            'orbital_scheme': 'period',
            'orbital_period_scheme': 'Sana+12_period_extended',
            'orbital_period_min': 0.35,
            'orbital_period_max': 6000,
        }
        period_pdf_func = initial_distributions.get_orbital_distribution(kwargs).pdf_per_dex
        # Test that the function is callable and returns sensible values
        result = period_pdf_func(10.0, 15.0)  # P=10 days, m1=15 Msun
        assert isinstance(result, (float, np.ndarray))
        assert np.all(result >= 0)

    def test_power_law_period(self):
        kwargs = {
            'orbital_scheme': 'period',
            'orbital_period_scheme': 'power_law',
            'orbital_period_min': 1.4,
            'orbital_period_max': 3e6,
            'power_law_slope': -0.55,
        }
        period_pdf_func = initial_distributions.get_orbital_distribution(kwargs).pdf_per_dex
        # Test that the function is callable and returns sensible values
        result = period_pdf_func(100.0, 20.0)  # P=100 days, m1=20 Msun
        assert isinstance(result, (float, np.ndarray))
        assert np.all(result >= 0)
        result2 = period_pdf_func(100.0, 20.0, 1.0)
        assert isinstance(result2, (float, np.ndarray))
        assert np.all(result2 >= 0)

    def test_log_power_law_period(self):
        kwargs = {
            'orbital_scheme': 'period',
            'orbital_period_scheme': 'log_power_law',
            'orbital_period_min': 1.4,
            'orbital_period_max': 3e6,
            'power_law_slope': -0.55,
        }
        period_pdf_func = initial_distributions.get_orbital_distribution(kwargs).pdf_per_dex
        result = period_pdf_func(100.0, 20.0, 1.0)
        assert np.all(result > 0)
        # (log10 P)^-0.55 shape per dex
        ratio = period_pdf_func(1e4, 20.0) / period_pdf_func(100.0, 20.0)
        assert np.isclose(ratio, 2.0**-0.55)

    @mark.parametrize("scheme, m1", [('Sana+12_period_extended', 10.0),
                                     ('Sana+12_period_extended', 30.0),
                                     ('power_law', 10.0),
                                     ('log_power_law', 10.0)])
    def test_period_pdf_normalised_in_logP(self, scheme, m1):
        kwargs = {
            'orbital_scheme': 'period',
            'orbital_period_scheme': scheme,
            'orbital_period_min': 1.4,
            'orbital_period_max': 3e6,
            'power_law_slope': -0.55,
        }
        period_pdf_func = initial_distributions.get_orbital_distribution(kwargs).pdf_per_dex
        integral, _ = quad(lambda logP: period_pdf_func(10**logP, m1, 1.0),
                           np.log10(1.4), np.log10(3e6))
        assert np.isclose(integral, 1.0)

    def test_log_uniform_separation(self):
        kwargs = {
            'orbital_scheme': 'separation',
            'orbital_separation_scheme': 'log_uniform',
            'orbital_separation_min': 1.0,
            'orbital_separation_max': 1000.0,
        }
        separation_pdf_func = initial_distributions.get_orbital_distribution(kwargs).pdf_per_dex
        result = separation_pdf_func(10.0, 15.0, 0.8)  # a=10 Rsun, m1=15 Msun
        assert isinstance(result, (float, np.ndarray))
        assert np.all(result >= 0)

    @mark.parametrize("m1, q", [(10.0, 1.0), (30.0, 0.5), (100.0, 0.1)])
    def test_log_uniform_separation_normalised_in_logP(self, m1, q):
        a_min, a_max = 5.0, 1e5
        kwargs = {
            'orbital_scheme': 'separation',
            'orbital_separation_scheme': 'log_uniform',
            'orbital_separation_min': a_min,
            'orbital_separation_max': a_max,
        }
        separation_pdf_func = initial_distributions.get_orbital_distribution(kwargs).pdf_per_dex
        logP_min = np.log10(orbital_period_from_separation(a_min, m1, q*m1))
        logP_max = np.log10(orbital_period_from_separation(a_max, m1, q*m1))
        integral, _ = quad(lambda logP: separation_pdf_func(10**logP, m1, q),
                           logP_min, logP_max)
        assert np.isclose(integral, 1.0)

    @mark.parametrize("kwargs, expected", [
        ({'orbital_scheme': 'period',
          'orbital_period_scheme': 'Sana+12_period_extended',
          'orbital_period_min': 0.35, 'orbital_period_max': 6000},
         distributions.Sana12Period),
        ({'orbital_scheme': 'period', 'orbital_period_scheme': 'power_law',
          'orbital_period_min': 1.4, 'orbital_period_max': 3e6,
          'power_law_slope': -0.55},
         distributions.PowerLawPeriod),
        ({'orbital_scheme': 'period', 'orbital_period_scheme': 'log_power_law',
          'orbital_period_min': 1.4, 'orbital_period_max': 3e6,
          'power_law_slope': -0.55},
         distributions.LogPowerLawPeriod),
    ])
    def test_period_scheme_to_class(self, kwargs, expected):
        assert type(initial_distributions.get_orbital_distribution(kwargs)) is expected

    @mark.parametrize("scheme, separation_class, extra", [
        ('log_uniform', distributions.LogUniform, {}),
        ('log_normal', distributions.LogNormalSeparation,
         {'log_orbital_separation_mean': 2.0,
          'log_orbital_separation_sigma': 0.5}),
    ])
    def test_separation_scheme_to_class(self, scheme, separation_class,
                                        extra):
        orbit = initial_distributions.get_orbital_distribution({
            'orbital_scheme': 'separation',
            'orbital_separation_scheme': scheme,
            'orbital_separation_min': 5.0,
            'orbital_separation_max': 1e5, **extra})
        assert isinstance(orbit, distributions.PeriodFromSeparation)
        assert type(orbit.separation) is separation_class

    @mark.parametrize("m1, q", [(10.0, 1.0), (30.0, 0.5), (100.0, 0.1)])
    def test_log_normal_separation_normalised_in_logP(self, m1, q):
        a_min, a_max = 5.0, 1e5
        kwargs = {
            'orbital_scheme': 'separation',
            'orbital_separation_scheme': 'log_normal',
            'orbital_separation_min': a_min,
            'orbital_separation_max': a_max,
            'log_orbital_separation_mean': 2.0,
            'log_orbital_separation_sigma': 0.5,
        }
        separation_pdf_func = initial_distributions.get_orbital_distribution(kwargs).pdf_per_dex
        logP_min = np.log10(orbital_period_from_separation(a_min, m1, q*m1))
        logP_max = np.log10(orbital_period_from_separation(a_max, m1, q*m1))
        integral, _ = quad(lambda logP: separation_pdf_func(10**logP, m1, q),
                           logP_min, logP_max, limit=200)
        assert np.isclose(integral, 1.0)

    def test_invalid_period_scheme(self):
        kwargs = {
            'orbital_scheme': 'period',
            'orbital_period_scheme': 'invalid_scheme',
            'orbital_period_min': 0.35,
            'orbital_period_max': 6000,
        }
        with pytest.raises(ValueError) as excinfo:
            initial_distributions.get_orbital_distribution(kwargs)

        assert "Orbital period scheme not recognized" in str(excinfo.value)

        kwargs = {
            'orbital_scheme' : 'invalid_scheme',
        }

        with pytest.raises(ValueError) as excinfo:
            initial_distributions.get_orbital_distribution(kwargs)

        assert "Orbital scheme not recognized" in str(excinfo.value)

        kwargs = {
            'orbital_scheme': 'separation',
            'orbital_separation_scheme': 'invalid_scheme',
        }
        with pytest.raises(ValueError) as excinfo:
            initial_distributions.get_orbital_distribution(kwargs)

        assert "Orbital separation scheme not recognized" in str(excinfo.value)
