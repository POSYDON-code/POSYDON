"""Test norm_pop module for POSYDON."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>"
]

import os
import warnings

import numpy as np
import pandas as pd
import pytest
from pytest import fixture, mark
from scipy.integrate import quad

import posydon.popsyn.IMFs as IMFs
from posydon.config import PATH_TO_POSYDON

# file to test
from posydon.popsyn import distributions, norm_pop
from posydon.popsyn.independent_sample import generate_independent_samples
from posydon.utils.common_functions import orbital_period_from_separation
from posydon.utils.posydonwarning import ValueWarning


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


class TestGetIMFPdf:
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
            norm_pop.get_IMF_pdf(kwargs)

    def test_missing_imf_scheme_raises(self):
        kwargs = {'primary_mass_min': 1, 'primary_mass_max': 100}
        with pytest.raises(ValueError, match="primary_mass_scheme"):
            norm_pop.get_IMF_pdf(kwargs)

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
            norm_pop.get_IMF_pdf(kwargs)

    @mark.parametrize('scheme', ['Salpeter', 'Kroupa1993', 'Kroupa2001'])
    def test_sampler_imfs_supported(self, scheme):
        kwargs = {
            'primary_mass_scheme': scheme,
            'primary_mass_min': 7,
            'primary_mass_max': 140
        }
        pdf_func = norm_pop.get_IMF_pdf(kwargs)
        assert np.isclose(quad(pdf_func, 7, 140)[0], 1.0)


    def test_valid_imf_returns_dummy_pdf(self, monkeypatch):
        # Use monkeypatch to add DummyIMF to IMFs module
        monkeypatch.setattr(IMFs, 'DummyIMF', DummyIMF, raising=False)

        kwargs = {
            'primary_mass_scheme': 'DummyIMF',
            'primary_mass_min': 1,
            'primary_mass_max': 100
        }
        pdf_func = norm_pop.get_IMF_pdf(kwargs)
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
        pdf_func = norm_pop.get_IMF_pdf(kwargs)

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
        q_pdf = norm_pop.get_mass_ratio_distribution(kwargs).pdf
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
        q_pdf = norm_pop.get_mass_ratio_distribution(kwargs).pdf
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
            norm_pop.get_mass_ratio_distribution(kwargs)

    def test_power_law_mass_ratio_pdf(self):
        """Test that power_law_mass_ratio scheme returns correct PDF values."""
        kwargs = {
            'secondary_mass_scheme': 'power_law_mass_ratio',
            'mass_ratio_slope': 0.0,
            'q_min': 0.1,
            'q_max': 0.9,
        }
        q_pdf = norm_pop.get_mass_ratio_distribution(kwargs).pdf
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
        q_pdf = norm_pop.get_mass_ratio_distribution(kwargs).pdf
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
            norm_pop.get_mass_ratio_distribution(kwargs)

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
            norm_pop.get_mass_ratio_distribution(kwargs)

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
        q_pdf = norm_pop.get_mass_ratio_distribution(kwargs).pdf
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
        q_dist = norm_pop.get_mass_ratio_distribution(kwargs)
        assert type(q_dist) is expected

    def test_m1_dependent_parameters(self):
        q_dist = norm_pop.get_mass_ratio_distribution({
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
        q_dist = norm_pop.get_mass_ratio_distribution({
            'secondary_mass_scheme': 'power_law_mass_ratio',
            'mass_ratio_slope': 1.0, **extra})
        assert q_dist.support() == (q_min, q_max)

    @mark.parametrize('scheme', ['power_law_mass_ratio', 'flat_mass_ratio'])
    def test_fixed_bounds_q_min_greater_than_q_max(self, scheme):
        kwargs = {'secondary_mass_scheme': scheme, 'mass_ratio_slope': 0.0,
                  'q_min': 0.8, 'q_max': 0.2}
        with pytest.raises(ValueError, match="q_min must be less than q_max"):
            norm_pop.get_mass_ratio_distribution(kwargs)

    @mark.parametrize('scheme', ['invalid', 'q=1', None])
    def test_unsupported_scheme(self, scheme):
        with pytest.raises(ValueError, match="secondary_mass_scheme"):
            norm_pop.get_mass_ratio_distribution(
                {'secondary_mass_scheme': scheme})


class TestGetBinaryFractionPdf:

    def test_const_binary_fraction_pdf(self):
        kwargs = {
            'binary_fraction_scheme': 'const',
            'binary_fraction_const': 0.7
        }
        pdf_func = norm_pop.get_binary_fraction_pdf(kwargs)
        # For binary stars, expected = 0.7
        result = pdf_func(True)
        assert np.allclose(result, np.asarray(0.7))
        # For single stars, expected = 1 - 0.7 = 0.3
        result = pdf_func(False)
        assert np.allclose(result, np.asarray(0.3))

    def test_invalid_binary_fraction_scheme(self):
        kwargs = {
            'binary_fraction_scheme': 'invalid'
        }
        with pytest.raises(ValueError) as excinfo:
            norm_pop.get_binary_fraction_pdf(kwargs)

        assert "Binary fraction scheme not recognized" in str(excinfo.value)

class TestGetPeriodPdf:
    def test_sana12_period_extended(self):
        kwargs = {
            'orbital_scheme': 'period',
            'orbital_period_scheme': 'Sana+12_period_extended',
            'orbital_period_min': 0.35,
            'orbital_period_max': 6000,
        }
        period_pdf_func = norm_pop.get_period_pdf(kwargs)
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
        period_pdf_func = norm_pop.get_period_pdf(kwargs)
        # Test that the function is callable and returns sensible values
        result = period_pdf_func(100.0, 20.0)  # P=100 days, m1=20 Msun
        assert isinstance(result, (float, np.ndarray))
        assert np.all(result >= 0)
        result2 = period_pdf_func(100.0, 20.0, 1.0)
        assert isinstance(result2, (float, np.ndarray))
        assert np.all(result2 >= 0)

    def test_log_uniform_separation(self):
        kwargs = {
            'orbital_scheme': 'separation',
            'orbital_separation_scheme': 'log_uniform',
            'orbital_separation_min': 1.0,
            'orbital_separation_max': 1000.0,
        }
        separation_pdf_func = norm_pop.get_period_pdf(kwargs)
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
        separation_pdf_func = norm_pop.get_period_pdf(kwargs)
        logP_min = np.log10(orbital_period_from_separation(a_min, m1, q*m1))
        logP_max = np.log10(orbital_period_from_separation(a_max, m1, q*m1))
        integral, _ = quad(lambda logP: separation_pdf_func(10**logP, m1, q),
                           logP_min, logP_max)
        assert np.isclose(integral, 1.0)

    def test_invalid_period_scheme(self):
        kwargs = {
            'orbital_scheme': 'period',
            'orbital_period_scheme': 'invalid_scheme',
            'orbital_period_min': 0.35,
            'orbital_period_max': 6000,
        }
        with pytest.raises(ValueError) as excinfo:
            norm_pop.get_period_pdf(kwargs)

        assert "Orbital period scheme not recognized" in str(excinfo.value)

        kwargs = {
            'orbital_scheme' : 'invalid_scheme',
        }

        with pytest.raises(ValueError) as excinfo:
            norm_pop.get_period_pdf(kwargs)

        assert "Orbital scheme not recognized" in str(excinfo.value)

        kwargs = {
            'orbital_scheme': 'separation',
            'orbital_separation_scheme': 'invalid_scheme',
        }
        with pytest.raises(ValueError) as excinfo:
            norm_pop.get_period_pdf(kwargs)

        assert "Orbital separation scheme not recognized" in str(excinfo.value)


class TestGetMeanMass:
    def test_q_min_greater_than_q_max_error(self):
        # Test the validation error when q_min > q_max
        params = {
            'primary_mass_scheme': 'FlatIMF',
            'primary_mass_min': 1,
            'primary_mass_max': 100,
            'secondary_mass_scheme': 'flat_mass_ratio',
            'binary_fraction_scheme': 'const',
            'binary_fraction_const': 0.5,
            'orbital_scheme': 'period',
            'orbital_period_scheme': 'Sana+12_period_extended',
            'orbital_period_min': 0.35,
            'orbital_period_max': 6000,
            'q_min': 0.8,  # Greater than q_max
            'q_max': 0.2,
        }

        with pytest.raises(ValueError) as excinfo:
            norm_pop.get_mean_mass(params)

        assert "q_min must be less than q_max" in str(excinfo.value)

    def test_empty_mass_ratio_range_error(self):
        # m1-dependent bounds: q lies in (m2_min/m1, min(m2_max/m1, 1)].
        # With secondary_mass_min >= secondary_mass_max the range is empty
        # for every m1.
        params = {
            'primary_mass_scheme': 'FlatIMF',
            'primary_mass_min': 5,
            'primary_mass_max': 10,
            'secondary_mass_min': 4,
            'secondary_mass_max': 3,
            'secondary_mass_scheme': 'flat_mass_ratio',
            'binary_fraction_scheme': 'const',
            'binary_fraction_const': 0.5,
            'orbital_scheme': 'period',
            'orbital_period_scheme': 'Sana+12_period_extended',
            'orbital_period_min': 0.35,
            'orbital_period_max': 6000,
        }

        with pytest.raises(ValueError) as excinfo:
            norm_pop.get_mean_mass(params)

        assert "m2_min must be less than m2_max" in str(excinfo.value)

    @mark.parametrize('imf, m1_min, m1_max, m2_min, m2_max, f_b', [
        # m2_min/m1_min > m2_max/m1_max, but every m1 has a valid q range;
        # the sampler accepts this setup, so it must not raise.
        ('FlatIMF', 5, 10, 4, 6, 0.5),
        # typical simulation setup (base_simulation_kwargs)
        ('Salpeter', 7, 140, 0.35, 140, 1.0),
        # m2_max inside the primary mass range
        ('Salpeter', 7, 140, 0.35, 50, 1.0),
        # default population parameters in SyntheticPopulation
        ('Kroupa2001', 0.01, 200, 0.0005, 200, 0.7),
    ])
    def test_mean_mass_flat_mass_ratio_analytic(self, imf, m1_min, m1_max,
                                                m2_min, m2_max, f_b):
        # flat_mass_ratio without q_min/q_max: q ~ U(m2_min/m1,
        # min(m2_max/m1, 1)], i.e. m2 ~ U(m2_min, min(m2_max, m1)], so
        # E[m2 | m1] = (m2_min + min(m2_max, m1)) / 2.
        params = {
            'primary_mass_scheme': imf,
            'primary_mass_min': m1_min,
            'primary_mass_max': m1_max,
            'secondary_mass_min': m2_min,
            'secondary_mass_max': m2_max,
            'secondary_mass_scheme': 'flat_mass_ratio',
            'binary_fraction_scheme': 'const',
            'binary_fraction_const': f_b,
            'orbital_scheme': 'period',
            'orbital_period_scheme': 'Sana+12_period_extended',
            'orbital_period_min': 0.35,
            'orbital_period_max': 6000,
        }
        imf_pdf = getattr(IMFs, imf)(m_min=m1_min, m_max=m1_max).pdf
        points = [m for m in (0.08, 0.5, m2_max) if m1_min < m < m1_max]
        E_m1 = quad(lambda m: m * imf_pdf(m), m1_min, m1_max,
                    points=points, limit=200)[0]
        E_m2 = quad(lambda m: imf_pdf(m) * 0.5 * (m2_min + min(m2_max, m)),
                    m1_min, m1_max, points=points, limit=200)[0]
        expected = E_m1 + f_b * E_m2

        assert np.isclose(norm_pop.get_mean_mass(params), expected, rtol=1e-3)

    @mark.parametrize('extra, q_min, q_max', [
        # no bounds given: PowerLawMassRatio defaults to q in (0.05, 1]
        ({}, 0.05, 1.0),
        # secondary_mass_min/max are not used by power_law_mass_ratio
        ({'secondary_mass_min': 3.5, 'secondary_mass_max': 20}, 0.05, 1.0),
        # explicit bounds
        ({'q_min': 0.2, 'q_max': 0.6}, 0.2, 0.6),
    ])
    def test_mean_mass_power_law_mass_ratio_analytic(self, extra,
                                                     q_min, q_max):
        alpha = 1.0
        params = {
            'primary_mass_scheme': 'Salpeter',
            'primary_mass_min': 7,
            'primary_mass_max': 140,
            'secondary_mass_scheme': 'power_law_mass_ratio',
            'mass_ratio_slope': alpha,
            'binary_fraction_scheme': 'const',
            'binary_fraction_const': 0.6,
            'orbital_scheme': 'period',
            'orbital_period_scheme': 'Sana+12_period_extended',
            'orbital_period_min': 0.35,
            'orbital_period_max': 6000,
            **extra,
        }
        imf_pdf = IMFs.Salpeter(m_min=7, m_max=140).pdf
        E_m1 = quad(lambda m: m * imf_pdf(m), 7, 140)[0]
        E_q = (quad(lambda q: q**(alpha + 1), q_min, q_max)[0]
               / quad(lambda q: q**alpha, q_min, q_max)[0])
        expected = E_m1 * (1 + 0.6 * E_q)

        assert np.isclose(norm_pop.get_mean_mass(params), expected, rtol=1e-4)

    def test_mean_mass_matches_sampler(self):
        # The normalisation must agree with what the sampler actually draws.
        params = {'number_of_binaries': 100000,
                  'primary_mass_scheme': 'Salpeter',
                  'primary_mass_min': 7,
                  'primary_mass_max': 140,
                  'secondary_mass_scheme': 'flat_mass_ratio',
                  'secondary_mass_min': 0.35,
                  'secondary_mass_max': 140,
                  'binary_fraction_scheme': 'const',
                  'binary_fraction_const': 0.7,
                  'orbital_scheme': 'period',
                  'orbital_period_scheme': 'Sana+12_period_extended',
                  'orbital_period_min': 0.35,
                  'orbital_period_max': 6e3,
                  'eccentricity_scheme': 'zero'}
        sample = pop_data(params.copy())
        sampled_mean = (sample['S1_mass_i'] + sample['S2_mass_i']).mean()
        # statistical error on the sampled mean is ~0.3%
        assert np.isclose(norm_pop.get_mean_mass(params), sampled_mean,
                          rtol=0.01)


class TestGetPdf:
    def test_single_star_pdf(self):
        # Using the test-only FlatIMF -> IMF_pdf returns a normalised constant.
        kwargs = {
            'primary_mass_scheme': 'FlatIMF',
            'primary_mass_min': 1,
            'primary_mass_max': 100,
            'secondary_mass_scheme': 'flat_mass_ratio',
            'secondary_mass_min': 1,
            'secondary_mass_max': 100,
            'binary_fraction_scheme': 'const',
            'binary_fraction_const': 0.5,
            'orbital_scheme': 'period',
            'orbital_period_scheme': 'Sana+12_period_extended',
            'orbital_period_min': 0.35,
            'orbital_period_max': 6000,
        }
        pdf_func = norm_pop.get_pdf(kwargs)
        m1_val = 10
        # For non-binary, expected = (1 - f_b)*IMF_pdf = 0.5*(1/99) = 0.5/99
        imf_val = 1.0 / (kwargs['primary_mass_max'] - kwargs['primary_mass_min'])
        expected = 0.5 * imf_val
        result = pdf_func(m1_val, binary=False)
        assert np.allclose(result, expected)

    def test_binary_pdf(self):
        # Using alternative mass ratio branch (q_min and q_max provided)
        kwargs = {
            'primary_mass_scheme': 'FlatIMF',
            'primary_mass_min': 1,
            'primary_mass_max': 100,
            'secondary_mass_scheme': 'flat_mass_ratio',
            'q_min': 0,
            'q_max': 1,
            'binary_fraction_scheme': 'const',
            'binary_fraction_const': 0.3,
            'orbital_scheme': 'period',
            'orbital_period_scheme': 'Sana+12_period_extended',
            'orbital_period_min': 0.35,
            'orbital_period_max': 6000,

        }
        pdf_func = norm_pop.get_pdf(kwargs, mass_pdf=True)
        m1_val = 10
        # For binary, expected = f_b*IMF_pdf*q_pdf = 0.3*(1/99)*1 = 0.3/99
        # for valid q (e.g. 0.5)
        imf_val = 1.0 / (kwargs['primary_mass_max'] - kwargs['primary_mass_min'])
        expected = 0.3 * imf_val * 1  # f_b * IMF_pdf * q_pdf
        result = pdf_func(m1_val, q=0.5, binary=True)
        assert np.allclose(result, expected)
        # q=0 should return 0 from mass ratio pdf.
        result_invalid = pdf_func(m1_val, q=0, binary=True)
        assert np.all(result_invalid == 0)

    def test_dummy_imf_pdf(self, monkeypatch):
        # Use monkeypatch fixture to add DummyIMF to IMFs module
        monkeypatch.setattr(IMFs, 'DummyIMF', DummyIMF, raising=False)

        kwargs = {
            'primary_mass_scheme': 'DummyIMF',
            'primary_mass_min': 1,
            'primary_mass_max': 100,
            'secondary_mass_scheme': 'flat_mass_ratio',
            'secondary_mass_min': 1,
            'secondary_mass_max': 100,
            'binary_fraction_scheme': 'const',
            'binary_fraction_const': 0.7,
            'orbital_scheme': 'period',
            'orbital_period_scheme': 'Sana+12_period_extended',
            'orbital_period_min': 0.35,
            'orbital_period_max': 6000,
        }
        pdf_func = norm_pop.get_pdf(kwargs, mass_pdf=True)
        m1_val = 10
        # For non-binary stars:
        # expected = (1 - f_b)*DummyIMF.pdf = 0.3*0.5 = 0.15
        result = pdf_func(m1_val, binary=False)
        assert np.allclose(result, 0.15)

###############################
# Test the reweighting function
###############################


@fixture
def base_simulation_kwargs():
    simulation_params = {'number_of_binaries': 100000,
                        'primary_mass_scheme': 'Kroupa2001',
                        'binary_fraction_const': 1.0,
                        'binary_fraction_scheme': 'const',
                        'primary_mass_min': 7,
                        'primary_mass_max': 140,
                        'secondary_mass_scheme': 'flat_mass_ratio',
                        'secondary_mass_min': 0.35,
                        'secondary_mass_max': 140,
                        'orbital_scheme': 'period',
                        'orbital_period_scheme': 'Sana+12_period_extended',
                        'orbital_period_min': 0.35,
                        'orbital_period_max': 6e3,
                        'eccentricity_scheme': 'zero',}
    return simulation_params

@fixture
def base_population_kwargs(base_simulation_kwargs):
    population_params = base_simulation_kwargs.copy()
    population_params['q_min'] = 0.0
    population_params['q_max'] = 1.0
    return population_params


def pop_data(kwargs):
    '''generate a population of stars'''
    RNG = np.random.default_rng(seed=42)
    kwargs['RNG'] = RNG
    sample = generate_independent_samples(**kwargs)
    pop_data = pd.DataFrame(np.array(sample).T,
                            columns=['orbital_period_i',
                                     'eccentricity_i',
                                     'S1_mass_i',
                                     'S2_mass_i',])
    if kwargs['orbital_scheme'] == 'separation':
        # sampler returns separations
        pop_data['orbital_period_i'] = orbital_period_from_separation(
            pop_data['orbital_period_i'],
            pop_data['S1_mass_i'],
            pop_data['S2_mass_i'])
    pop_data['state_i'] = 'detached'
    mask = pd.isna(pop_data['S2_mass_i'])
    pop_data['state_i'][mask] = 'initially_single_star'
    pop_data['S2_mass_i'][mask] = 0
    return pop_data

class TestReweighting():

    def test_population_larger_sample_space(self, base_simulation_kwargs):
        # no binaries
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 0.

        small_sample_primary_mass_max = 20
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = small_sample_primary_mass_max
        small_sample = pop_data(base_simulation_kwargs)

        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['primary_mass_min'] = 7.
        expanded_kwargs['primary_mass_max'] = 40

        expanded_sample = pop_data(expanded_kwargs)

        M_sim = (small_sample['S1_mass_i'].sum()
                 + small_sample['S2_mass_i'].sum())
        small_weights = norm_pop.calculate_model_weights(small_sample,
                                                         M_sim,
                                                         base_simulation_kwargs,
                                                         expanded_kwargs)
        M_sim = (expanded_sample['S1_mass_i'].sum()
                 + expanded_sample['S2_mass_i'].sum())
        expanded_weights = norm_pop.calculate_model_weights(expanded_sample,
                                                            M_sim,
                                                            expanded_kwargs,
                                                            expanded_kwargs)

        assert len(expanded_weights)  == len(expanded_sample)

        mask = expanded_sample['S1_mass_i'] <= small_sample_primary_mass_max
        selection = expanded_weights[mask]
        assert np.isclose(np.nansum(selection), np.nansum(small_weights), atol=1e-3)

    def test_population_larger_sample_space_binary(self, base_simulation_kwargs):
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 1.

        small_sample_primary_mass_max = 20
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = small_sample_primary_mass_max
        small_sample = pop_data(base_simulation_kwargs)

        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['primary_mass_min'] = 7.
        expanded_kwargs['primary_mass_max'] = 40

        expanded_sample = pop_data(expanded_kwargs)

        M_sim = (small_sample['S1_mass_i'].sum()
                 + small_sample['S2_mass_i'].sum())
        small_weights = norm_pop.calculate_model_weights(small_sample,
                                                         M_sim,
                                                         base_simulation_kwargs,
                                                         expanded_kwargs)

        M_sim = (expanded_sample['S1_mass_i'].sum()
                 + expanded_sample['S2_mass_i'].sum())
        expanded_weights = norm_pop.calculate_model_weights(expanded_sample,
                                                            M_sim,
                                                            expanded_kwargs,
                                                            expanded_kwargs)

        assert len(expanded_weights)  == len(expanded_sample)

        mask = expanded_sample['S1_mass_i'] <= small_sample_primary_mass_max
        selection = expanded_weights[mask]
        assert np.isclose(np.nansum(selection), np.nansum(small_weights), atol=1e-3)

    def test_population_larger_sample_space_mixed(self, base_simulation_kwargs):
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 0.7

        small_sample_primary_mass_max = 20
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = small_sample_primary_mass_max
        small_sample = pop_data(base_simulation_kwargs)

        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['primary_mass_min'] = 7.
        expanded_kwargs['primary_mass_max'] = 40

        expanded_sample = pop_data(expanded_kwargs)

        M_sim = (small_sample['S1_mass_i'].sum()
                 + small_sample['S2_mass_i'].sum())
        small_weights = norm_pop.calculate_model_weights(small_sample,
                                                         M_sim,
                                                         base_simulation_kwargs,
                                                         expanded_kwargs)

        M_sim = (expanded_sample['S1_mass_i'].sum()
                 + expanded_sample['S2_mass_i'].sum())
        expanded_weights = norm_pop.calculate_model_weights(expanded_sample,
                                                            M_sim,
                                                            expanded_kwargs,
                                                            expanded_kwargs)

        assert len(expanded_weights)  == len(expanded_sample)

        mask = expanded_sample['S1_mass_i'] <= small_sample_primary_mass_max
        selection = expanded_weights[mask]
        assert np.isclose(np.sum(selection), np.sum(small_weights), atol=1e-3)

    def test_population_larger_sample_space_binary_fraction_change(self,
                                                        base_simulation_kwargs):

        # no binaries
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 1.

        small_sample_primary_mass_max = 20
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = small_sample_primary_mass_max
        small_sample = pop_data(base_simulation_kwargs)


        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['primary_mass_min'] = 7.
        expanded_kwargs['primary_mass_max'] = 40
        expanded_kwargs['binary_fraction_const'] = 0.7

        expanded_sample = pop_data(expanded_kwargs)

        M_sim = (small_sample['S1_mass_i'].sum()
                 + small_sample['S2_mass_i'].sum())
        small_weights = norm_pop.calculate_model_weights(small_sample,
                                                         M_sim,
                                                         base_simulation_kwargs,
                                                         expanded_kwargs)

        M_sim = (expanded_sample['S1_mass_i'].sum()
                 + expanded_sample['S2_mass_i'].sum())
        expanded_weights = norm_pop.calculate_model_weights(expanded_sample,
                                                            M_sim,
                                                            expanded_kwargs,
                                                            expanded_kwargs)

        assert len(expanded_weights)  == len(expanded_sample)

        mask = ((expanded_sample['S1_mass_i'] <= small_sample_primary_mass_max)
                & (expanded_sample['state_i'] != 'initially_single_star'))
        selection = expanded_weights[mask]
        assert np.isclose(np.sum(selection), np.sum(small_weights), atol=1e-3)

    def test_population_single_stars_binary_fraction(self,
                                                    base_simulation_kwargs):
        # no binaries
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 0

        small_sample_primary_mass_max = 20
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = small_sample_primary_mass_max
        small_sample = pop_data(base_simulation_kwargs)

        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['primary_mass_min'] = 7.
        expanded_kwargs['primary_mass_max'] = 20
        expanded_kwargs['binary_fraction_const'] = 0.7

        expanded_sample = pop_data(expanded_kwargs)

        M_sim = (small_sample['S1_mass_i'].sum()
                 + small_sample['S2_mass_i'].sum())
        small_weights = norm_pop.calculate_model_weights(small_sample,
                                                         M_sim,
                                                         base_simulation_kwargs,
                                                         expanded_kwargs)

        M_sim = (expanded_sample['S1_mass_i'].sum()
                 + expanded_sample['S2_mass_i'].sum())
        expanded_weights = norm_pop.calculate_model_weights(expanded_sample,
                                                            M_sim,
                                                            expanded_kwargs,
                                                            expanded_kwargs)

        assert len(expanded_weights)  == len(expanded_sample)

        mask = ((expanded_sample['S1_mass_i'] <= small_sample_primary_mass_max)
                & (expanded_sample['state_i'] == 'initially_single_star'))
        selection = expanded_weights[mask]
        assert np.isclose(np.sum(selection), np.sum(small_weights), atol=1e-3)

    def test_Salpeter_IMF_single(self, base_simulation_kwargs):
        # no binaries
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 0.
        base_simulation_kwargs['primary_mass_scheme'] = 'Salpeter'

        small_sample_primary_mass_max = 20
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = small_sample_primary_mass_max
        small_sample = pop_data(base_simulation_kwargs)

        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['primary_mass_min'] = 7.
        expanded_kwargs['primary_mass_max'] = 20
        expanded_kwargs['primary_mass_scheme'] = 'Salpeter'

        expanded_sample = pop_data(expanded_kwargs)

        M_sim = (small_sample['S1_mass_i'].sum()
                 + small_sample['S2_mass_i'].sum())
        small_weights = norm_pop.calculate_model_weights(small_sample,
                                                         M_sim,
                                                         base_simulation_kwargs,
                                                         expanded_kwargs)

        M_sim = (expanded_sample['S1_mass_i'].sum()
                 + expanded_sample['S2_mass_i'].sum())
        expanded_weights = norm_pop.calculate_model_weights(expanded_sample,
                                                            M_sim,
                                                            expanded_kwargs,
                                                            expanded_kwargs)

        assert len(expanded_weights)  == len(expanded_sample)

        mask = expanded_sample['S1_mass_i'] <= small_sample_primary_mass_max
        selection = expanded_weights[mask]
        assert np.isclose(np.sum(selection), np.sum(small_weights), atol=1e-3)

    def test_changing_IMF_single(self, base_simulation_kwargs):
        # no binaries
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 0.
        base_simulation_kwargs['primary_mass_scheme'] = 'Salpeter'

        small_sample_primary_mass_max = 20
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = small_sample_primary_mass_max
        small_sample = pop_data(base_simulation_kwargs)

        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['primary_mass_min'] = 7.
        expanded_kwargs['primary_mass_max'] = 20
        expanded_kwargs['primary_mass_scheme'] = 'Kroupa2001'

        expanded_sample = pop_data(expanded_kwargs)

        M_sim = (small_sample['S1_mass_i'].sum()
                 + small_sample['S2_mass_i'].sum())
        small_weights = norm_pop.calculate_model_weights(small_sample,
                                                         M_sim,
                                                         base_simulation_kwargs,
                                                         expanded_kwargs)

        M_sim = (expanded_sample['S1_mass_i'].sum()
                 + expanded_sample['S2_mass_i'].sum())
        expanded_weights = norm_pop.calculate_model_weights(expanded_sample,
                                                            M_sim,
                                                            expanded_kwargs,
                                                            expanded_kwargs)

        assert len(expanded_weights)  == len(expanded_sample)

        mask = expanded_sample['S1_mass_i'] <= small_sample_primary_mass_max
        selection = expanded_weights[mask]
        assert np.isclose(np.sum(selection), np.sum(small_weights), atol=1e-3)

    def test_changing_IMF_binary(self, base_simulation_kwargs):
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 1
        base_simulation_kwargs['primary_mass_scheme'] = 'Salpeter'

        small_sample_primary_mass_max = 20
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = small_sample_primary_mass_max
        small_sample = pop_data(base_simulation_kwargs)

        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['primary_mass_min'] = 7.
        expanded_kwargs['primary_mass_max'] = 20
        expanded_kwargs['primary_mass_scheme'] = 'Kroupa2001'

        expanded_sample = pop_data(expanded_kwargs)

        M_sim = (small_sample['S1_mass_i'].sum()
                 + small_sample['S2_mass_i'].sum())
        small_weights = norm_pop.calculate_model_weights(small_sample,
                                                         M_sim,
                                                         base_simulation_kwargs,
                                                         expanded_kwargs)

        M_sim = (expanded_sample['S1_mass_i'].sum()
                 + expanded_sample['S2_mass_i'].sum())
        expanded_weights = norm_pop.calculate_model_weights(expanded_sample,
                                                            M_sim,
                                                            expanded_kwargs,
                                                            expanded_kwargs)

        assert len(expanded_weights)  == len(expanded_sample)

        mask = expanded_sample['S1_mass_i'] <= small_sample_primary_mass_max
        selection = expanded_weights[mask]
        assert np.isclose(np.sum(selection), np.sum(small_weights), atol=1e-3)

    def test_changing_IMF_binary_range(self, base_simulation_kwargs):
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 1
        base_simulation_kwargs['primary_mass_scheme'] = 'Salpeter'

        small_sample_primary_mass_max = 20
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = small_sample_primary_mass_max
        small_sample = pop_data(base_simulation_kwargs)

        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['primary_mass_min'] = 7.
        expanded_kwargs['primary_mass_max'] = 40
        expanded_kwargs['primary_mass_scheme'] = 'Kroupa2001'

        expanded_sample = pop_data(expanded_kwargs)

        M_sim = (small_sample['S1_mass_i'].sum()
                 + small_sample['S2_mass_i'].sum())
        small_weights = norm_pop.calculate_model_weights(small_sample,
                                                         M_sim,
                                                         base_simulation_kwargs,
                                                         expanded_kwargs)

        M_sim = (expanded_sample['S1_mass_i'].sum()
                 + expanded_sample['S2_mass_i'].sum())
        expanded_weights = norm_pop.calculate_model_weights(expanded_sample,
                                                            M_sim,
                                                            expanded_kwargs,
                                                            expanded_kwargs)

        assert len(expanded_weights)  == len(expanded_sample)

        mask = expanded_sample['S1_mass_i'] <= small_sample_primary_mass_max
        selection = expanded_weights[mask]
        assert np.isclose(np.sum(selection), np.sum(small_weights), atol=1e-3)

    def test_same_sample_space(self, base_simulation_kwargs):
        '''Same sample space the weights should be the same as the simulation'''

        base_pop_data = pop_data(base_simulation_kwargs)
        M_sim = (base_pop_data['S1_mass_i'].sum()
                 + base_pop_data['S2_mass_i'].sum())
        weights = norm_pop.calculate_model_weights(base_pop_data,
                                                   M_sim,
                                                   base_simulation_kwargs,
                                                   base_simulation_kwargs)
        assert len(weights) == len(base_pop_data)
        assert np.all(weights >= 0)
        assert np.allclose(weights, 1/M_sim)

    @mark.parametrize('pop_m1_max', [140, 250])
    def test_systems_outside_simulation_warn(self, base_simulation_kwargs,
                                             pop_m1_max):
        '''Systems with zero simulation probability get weight 0 and one
        warning with their count. With pop_m1_max=140 the population PDF is
        zero too (0/0 = nan), with 250 it is not (x/0 = inf).'''
        base_simulation_kwargs['number_of_binaries'] = 1000
        data = pop_data(base_simulation_kwargs.copy())
        # two systems above the simulated primary_mass_max of 140 Msun
        data.loc[[0, 1], 'S1_mass_i'] = 200.
        data.loc[[0, 1], 'S2_mass_i'] = 100.
        data.loc[[0, 1], 'state_i'] = 'detached'
        population_kwargs = base_simulation_kwargs.copy()
        population_kwargs['primary_mass_max'] = pop_m1_max
        M_sim = data['S1_mass_i'].sum() + data['S2_mass_i'].sum()

        with pytest.warns(ValueWarning,
                          match="2 of 1000 systems have zero probability "
                                "under the simulation parameters"):
            weights = norm_pop.calculate_model_weights(
                data, M_sim, base_simulation_kwargs, population_kwargs)

        assert np.all(np.isfinite(weights))
        assert np.all(weights[:2] == 0)
        assert np.all(weights[2:] > 0)

    def test_narrower_population_does_not_warn(self, base_simulation_kwargs):
        '''Systems outside the requested population (but inside the
        simulation) get weight 0 without the outside-simulation warning.'''
        base_simulation_kwargs['number_of_binaries'] = 1000
        data = pop_data(base_simulation_kwargs.copy())
        population_kwargs = base_simulation_kwargs.copy()
        population_kwargs['primary_mass_max'] = 20
        M_sim = data['S1_mass_i'].sum() + data['S2_mass_i'].sum()

        with warnings.catch_warnings():
            warnings.filterwarnings("error", category=ValueWarning,
                                    message=r"\d+ of \d+ systems have zero")
            weights = norm_pop.calculate_model_weights(
                data, M_sim, base_simulation_kwargs, population_kwargs)

        outside_population = data['S1_mass_i'] > 20
        assert np.any(outside_population)
        assert np.all(weights[outside_population] == 0)
        assert np.all(weights[~outside_population] > 0)


    def test_mass_ratio_extension(self,
                                  base_simulation_kwargs,
                                  base_population_kwargs,):
        '''Sample space for mass ratio extension'''

        base_pop_data = pop_data(base_simulation_kwargs)
        M_sim = (base_pop_data['S1_mass_i'].sum()
                 + base_pop_data['S2_mass_i'].sum())
        weights = norm_pop.calculate_model_weights(base_pop_data,
                                                   M_sim,
                                                   base_simulation_kwargs,
                                                   base_population_kwargs)
        assert len(weights) == len(base_pop_data)
        assert np.all(weights >= 0)


@fixture
def separation_kwargs(base_simulation_kwargs):
    kwargs = base_simulation_kwargs.copy()
    for key in ['orbital_period_scheme', 'orbital_period_min',
                'orbital_period_max']:
        kwargs.pop(key)
    kwargs['orbital_scheme'] = 'separation'
    kwargs['orbital_separation_scheme'] = 'log_uniform'
    kwargs['orbital_separation_min'] = 5.0
    kwargs['orbital_separation_max'] = 1e5
    return kwargs

def check_reweighting(sim_kwargs, pop_kwargs):
    '''Compare reweighted sim to direct sampling of the population.
    The population's orbital range must be inside the sim's for all masses.'''
    sim_data = pop_data(sim_kwargs)
    M_sim = sim_data['S1_mass_i'].sum() + sim_data['S2_mass_i'].sum()
    weights = norm_pop.calculate_model_weights(sim_data, M_sim,
                                               sim_kwargs, pop_kwargs)

    direct_data = pop_data(pop_kwargs)
    M_direct = direct_data['S1_mass_i'].sum() + direct_data['S2_mass_i'].sum()

    assert np.isclose(weights.sum(), len(direct_data) / M_direct, rtol=0.03)

    logP_sim = np.log10(sim_data['orbital_period_i'])
    logP_direct = np.log10(direct_data['orbital_period_i'])
    assert np.isclose(np.average(logP_sim, weights=weights),
                      logP_direct.mean(), atol=0.05)
    assert np.isclose(np.sqrt(np.cov(logP_sim, aweights=weights)),
                      logP_direct.std(), rtol=0.05)


class TestOrbitalSchemeReweighting():
    def test_separation_to_period(self, separation_kwargs,
                                  base_simulation_kwargs):
        # a in [5, 1e5] Rsun covers P in [1, 3000] d for all masses
        base_simulation_kwargs['orbital_period_min'] = 1.0
        base_simulation_kwargs['orbital_period_max'] = 3e3
        check_reweighting(separation_kwargs, base_simulation_kwargs)

    def test_period_to_separation(self, base_simulation_kwargs,
                                  separation_kwargs):
        # P in [0.35, 6000] d covers a in [20, 2000] Rsun for all masses
        separation_kwargs['orbital_separation_min'] = 20.0
        separation_kwargs['orbital_separation_max'] = 2e3
        check_reweighting(base_simulation_kwargs, separation_kwargs)

    def test_separation_to_separation(self, separation_kwargs):
        pop_kwargs = separation_kwargs.copy()
        pop_kwargs['orbital_separation_min'] = 20.0
        pop_kwargs['orbital_separation_max'] = 2e3
        check_reweighting(separation_kwargs, pop_kwargs)


class TestBinaryFractions():

    def test_single_star_population(self, base_simulation_kwargs):
        # no binaries
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 0
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = 20

        small_sample = pop_data(base_simulation_kwargs)

        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['binary_fraction_const'] = 0.7
        expanded_sample = pop_data(expanded_kwargs)

        M_sim = (small_sample['S1_mass_i'].sum()
                 + small_sample['S2_mass_i'].sum())
        small_weights = norm_pop.calculate_model_weights(small_sample,
                                                         M_sim,
                                                         base_simulation_kwargs,
                                                         expanded_kwargs)

        M_sim = (expanded_sample['S1_mass_i'].sum()
                 + expanded_sample['S2_mass_i'].sum())
        expanded_weights = norm_pop.calculate_model_weights(expanded_sample,
                                                            M_sim,
                                                            expanded_kwargs,
                                                            expanded_kwargs)

        assert len(expanded_weights)  == len(expanded_sample)

        # Check just single stars in the population.
        mask = (expanded_sample['state_i'] == 'initially_single_star')
        selection = expanded_weights[mask]
        assert np.isclose(np.nansum(selection), np.nansum(small_weights), atol=1e-3)

    def test_single_star_population_error(self, base_simulation_kwargs):
        # no binaries
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 0
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = 20

        small_sample = pop_data(base_simulation_kwargs)

        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['binary_fraction_const'] = 1.0

        M_sim = (small_sample['S1_mass_i'].sum()
                 + small_sample['S2_mass_i'].sum())
        with pytest.raises(ValueError) as e:
            small_weights = norm_pop.calculate_model_weights(
                                                    small_sample,
                                                    M_sim,
                                                    base_simulation_kwargs,
                                                    expanded_kwargs)
        assert 'No binaries simulated, but requested' in str(e.value)

    def test_binary_star_population_error(self, base_simulation_kwargs):
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 1.0
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = 20

        small_sample = pop_data(base_simulation_kwargs)

        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['binary_fraction_const'] = 0

        M_sim = (small_sample['S1_mass_i'].sum()
                 + small_sample['S2_mass_i'].sum())
        with pytest.raises(ValueError) as e:
            small_weights = norm_pop.calculate_model_weights(
                                                        small_sample,
                                                        M_sim,
                                                        base_simulation_kwargs,
                                                        expanded_kwargs)
        assert 'No single stars simulated, but requested' in str(e.value)


    def test_mixed_population(self, base_simulation_kwargs):
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 0.3
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = 20

        small_sample = pop_data(base_simulation_kwargs)

        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['binary_fraction_const'] = 0.5
        expanded_sample = pop_data(expanded_kwargs)

        M_sim = (small_sample['S1_mass_i'].sum()
                 + small_sample['S2_mass_i'].sum())
        small_weights = norm_pop.calculate_model_weights(small_sample,
                                                         M_sim,
                                                         base_simulation_kwargs,
                                                         expanded_kwargs)

        M_sim = (expanded_sample['S1_mass_i'].sum()
                 + expanded_sample['S2_mass_i'].sum())
        expanded_weights = norm_pop.calculate_model_weights(expanded_sample,
                                                            M_sim,
                                                            expanded_kwargs,
                                                            expanded_kwargs)

        assert len(expanded_weights)  == len(expanded_sample)

        # Check just single stars in the population; need selection on both.
        mask = ((expanded_sample['S1_mass_i'] <= 20)
                & (expanded_sample['state_i'] == 'initially_single_star'))
        selection = expanded_weights[mask]

        mask2 = ((small_sample['S1_mass_i'] <= 20)
                 & (small_sample['state_i'] == 'initially_single_star'))
        assert np.isclose(np.sum(selection),
                          np.sum(small_weights[mask2]),
                          atol=1e-3)
