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
from posydon.popsyn import norm_pop
from posydon.popsyn.independent_sample import generate_independent_samples
from posydon.popsyn.initial_distributions import (
    get_IMF_distribution,
    get_Moe_17_PsandQs,
)
from posydon.unit_tests.popsyn.test_initial_distributions import (  # noqa: F401
    DummyIMF,
    FlatIMF,
    register_test_imfs,
)
from posydon.utils.common_functions import orbital_period_from_separation
from posydon.utils.posydonwarning import InappropriateValueWarning, ValueWarning


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
    pop_data.loc[mask, 'state_i'] = 'initially_single_star'
    pop_data.loc[mask, 'S2_mass_i'] = 0
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
                          match="2 of 1000 systems lie outside the "
                                "sampled parameter space"):
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
                                    message=r"\d+ of \d+ systems lie outside")
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

    def test_log_normal_separation_to_period(self, separation_kwargs,
                                             base_simulation_kwargs):
        separation_kwargs['orbital_separation_scheme'] = 'log_normal'
        separation_kwargs['log_orbital_separation_mean'] = 2.5
        separation_kwargs['log_orbital_separation_sigma'] = 1.0
        # a in [5, 1e5] Rsun covers P in [1, 3000] d for all masses
        base_simulation_kwargs['orbital_period_min'] = 1.0
        base_simulation_kwargs['orbital_period_max'] = 3e3
        check_reweighting(separation_kwargs, base_simulation_kwargs)

    def test_period_to_log_normal_separation(self, base_simulation_kwargs,
                                             separation_kwargs):
        # P in [0.35, 6000] d covers a in [20, 2000] Rsun for all masses
        separation_kwargs['orbital_separation_scheme'] = 'log_normal'
        separation_kwargs['log_orbital_separation_mean'] = 2.5
        separation_kwargs['log_orbital_separation_sigma'] = 0.5
        separation_kwargs['orbital_separation_min'] = 20.0
        separation_kwargs['orbital_separation_max'] = 2e3
        check_reweighting(base_simulation_kwargs, separation_kwargs)

    def test_separation_to_separation(self, separation_kwargs):
        pop_kwargs = separation_kwargs.copy()
        pop_kwargs['orbital_separation_min'] = 20.0
        pop_kwargs['orbital_separation_max'] = 2e3
        check_reweighting(separation_kwargs, pop_kwargs)


class TestBinaryFractions():

    def test_unsupported_scheme_in_model_weights(self):
        # an unsupported scheme must give a clear ValueError, not a
        # KeyError for the missing binary_fraction_const
        params = {'binary_fraction_scheme': 'Moe+17-massdependent'}
        with pytest.raises(ValueError,
                           match="Binary fraction scheme not recognized"):
            norm_pop.calculate_model_weights(None, 1.0, params, params)

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


###########################################
# Moe & Di Stefano (2017) distributions
###########################################

@fixture
def moe_kwargs(base_simulation_kwargs):
    kwargs = base_simulation_kwargs.copy()
    for key in ['binary_fraction_const', 'binary_fraction_scheme',
                'secondary_mass_min', 'secondary_mass_max']:
        kwargs.pop(key)
    kwargs['primary_mass_max'] = 40.0
    kwargs['secondary_mass_scheme'] = 'Moe+17-PsandQs'
    kwargs['orbital_period_scheme'] = 'Moe+17-PsandQs'
    kwargs['orbital_period_min'] = 10**0.15
    kwargs['orbital_period_max'] = 1e4
    return kwargs


def binary_stats(data, weights=None):
    '''Binary fraction, and mean/std of log10(P) and mean q of the binaries.'''
    w = np.ones(len(data)) if weights is None else np.asarray(weights)
    binary = (data['state_i'] != 'initially_single_star').to_numpy()
    logP = np.log10(data['orbital_period_i'].to_numpy()[binary])
    q = (data['S2_mass_i'] / data['S1_mass_i']).to_numpy()[binary]
    wb = w[binary]
    return (wb.sum() / w.sum(),
            np.average(logP, weights=wb),
            np.sqrt(np.cov(logP, aweights=wb)),
            np.average(q, weights=wb))


def check_reweighting_mixed(sim_kwargs, pop_kwargs):
    '''Compare a reweighted mixed simulation to direct sampling.'''
    sim_data = pop_data(sim_kwargs)
    M_sim = sim_data['S1_mass_i'].sum() + sim_data['S2_mass_i'].sum()
    weights = norm_pop.calculate_model_weights(sim_data, M_sim,
                                               sim_kwargs, pop_kwargs)
    direct_data = pop_data(pop_kwargs)
    M_direct = direct_data['S1_mass_i'].sum() + direct_data['S2_mass_i'].sum()

    # number of systems per unit mass
    assert np.isclose(weights.sum(), len(direct_data) / M_direct, rtol=0.03)
    reweighted = binary_stats(sim_data, weights)
    direct = binary_stats(direct_data)
    assert np.isclose(reweighted[0], direct[0], atol=0.01)
    assert np.isclose(reweighted[1], direct[1], atol=0.05)
    assert np.isclose(reweighted[2], direct[2], rtol=0.05)
    assert np.isclose(reweighted[3], direct[3], atol=0.01)


class TestMoe17PsandQs():

    def test_single_star_pdf(self, moe_kwargs):
        moe = get_Moe_17_PsandQs(moe_kwargs)
        m1 = np.array([0.05, 1.0, 10.0, 30.0])
        f_b = moe.binary_fraction(m1)
        # binary_fraction_* are ignored
        kwargs = moe_kwargs.copy()
        kwargs['binary_fraction_scheme'] = 'const'
        kwargs['binary_fraction_const'] = 0.3
        pdf = norm_pop.get_pdf(kwargs)
        IMF_pdf = get_IMF_distribution(kwargs).pdf
        assert np.allclose(pdf(m1, binary=False), IMF_pdf(m1) * (1 - f_b))

    def test_pdf(self, moe_kwargs):
        moe = get_Moe_17_PsandQs(moe_kwargs)
        IMF_pdf = get_IMF_distribution(moe_kwargs).pdf
        pdf = norm_pop.get_pdf(moe_kwargs)
        m1 = np.array([8.0, 15.0, 30.0, 12.0])
        q = np.array([0.5, 0.2, 0.97, 0.0])
        P = np.array([5.0, 300.0, 2e3, np.nan])
        binary = np.array([True, True, True, False])
        result = pdf(m1, q=q, P=P, binary=binary)
        # binaries: IMF * joint (q, P) density, with f_b counted once
        assert np.allclose(result[:3],
                           IMF_pdf(m1[:3]) * moe.pdf(m1[:3], P[:3], q[:3]))
        # single stars
        assert np.isclose(result[3],
                          IMF_pdf(m1[3]) * (1 - moe.binary_fraction(m1[3])))
        # outside the period range or below q = 0.1
        assert pdf(10.0, q=0.5, P=2e4, binary=True) == 0.0
        assert pdf(10.0, q=0.05, P=10.0, binary=True) == 0.0

    def test_mass_pdf_not_available(self, moe_kwargs):
        with pytest.raises(ValueError, match="mass_pdf=False"):
            norm_pop.get_pdf(moe_kwargs, mass_pdf=True)

    def test_mean_mass(self, moe_kwargs):
        mean_mass = norm_pop.get_mean_mass(moe_kwargs)
        # larger than the mean primary mass alone
        IMF_pdf = get_IMF_distribution(moe_kwargs).pdf
        mean_m1 = quad(lambda m: m * IMF_pdf(m),
                       moe_kwargs['primary_mass_min'],
                       moe_kwargs['primary_mass_max'])[0]
        assert mean_mass > 1.2 * mean_m1
        # matches the mean mass of a sampled population
        sample = pop_data(moe_kwargs)
        mass = sample['S1_mass_i'] + sample['S2_mass_i']
        assert np.isclose(mean_mass, mass.mean(), rtol=5e-3)

    def test_mean_mass_low_masses(self, moe_kwargs):
        # across the low-mass extension, f_b = 0 at M_min, and the kinks of
        # the mass interpolation
        moe_kwargs['primary_mass_min'] = 0.08
        moe_kwargs['primary_mass_max'] = 10.0
        mean_mass = norm_pop.get_mean_mass(moe_kwargs)
        sample = pop_data(moe_kwargs)
        mass = sample['S1_mass_i'] + sample['S2_mass_i']
        assert np.isclose(mean_mass, mass.mean(), rtol=1e-2)

    def test_mean_mass_const(self, base_simulation_kwargs):
        # matches the mean mass of a sampled population with f_b < 1
        base_simulation_kwargs['binary_fraction_const'] = 0.6
        mean_mass = norm_pop.get_mean_mass(base_simulation_kwargs)
        sample = pop_data(base_simulation_kwargs)
        mass = sample['S1_mass_i'] + sample['S2_mass_i']
        assert np.isclose(mean_mass, mass.mean(), rtol=5e-3)

    def test_coverage_error(self, moe_kwargs, base_simulation_kwargs):
        # Moe has no binaries at M1 <= M_min = 0.08
        moe_kwargs['primary_mass_min'] = 0.05
        pop = base_simulation_kwargs.copy()
        pop['primary_mass_min'] = 0.05
        pop['secondary_mass_min'] = 0.01
        pop['binary_fraction_const'] = 1.0
        with pytest.raises(ValueError,
                           match="No binaries simulated, but requested for "
                                 "primary masses 0.05-0.0"):
            norm_pop.check_binary_fraction_coverage(moe_kwargs, pop)
        # all single stars requested, only binaries simulated
        sim = base_simulation_kwargs.copy()
        pop = base_simulation_kwargs.copy()
        pop['binary_fraction_const'] = 0.0
        with pytest.raises(ValueError, match="No single stars simulated"):
            norm_pop.check_binary_fraction_coverage(sim, pop)

    def test_coverage_warning(self, moe_kwargs, base_simulation_kwargs):
        # Moe single stars are not in a binary-only simulation
        with pytest.warns(InappropriateValueWarning,
                          match="No single stars simulated"):
            norm_pop.check_binary_fraction_coverage(base_simulation_kwargs,
                                                    moe_kwargs)
        # binaries below M_min requested from a Moe simulation
        moe_kwargs['primary_mass_min'] = 0.05
        pop = base_simulation_kwargs.copy()
        pop['primary_mass_min'] = 0.05
        pop['binary_fraction_const'] = 0.5
        with pytest.warns(InappropriateValueWarning,
                          match="No binaries simulated"):
            norm_pop.check_binary_fraction_coverage(moe_kwargs, pop)

    def test_coverage_ok(self, moe_kwargs, base_simulation_kwargs):
        import warnings
        pop = base_simulation_kwargs.copy()
        pop['binary_fraction_const'] = 0.5
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            norm_pop.check_binary_fraction_coverage(moe_kwargs, pop)
            # no overlap of the primary-mass ranges: nothing to check
            pop['primary_mass_min'] = 50.0
            pop['primary_mass_max'] = 100.0
            pop['binary_fraction_const'] = 1.0
            sim = base_simulation_kwargs.copy()
            sim['primary_mass_max'] = 40.0
            sim['binary_fraction_const'] = 0.0
            norm_pop.check_binary_fraction_coverage(sim, pop)

    def test_Moe_to_Moe(self, moe_kwargs):
        pop = moe_kwargs.copy()
        pop['orbital_period_min'] = 10.0
        pop['orbital_period_max'] = 3e3
        check_reweighting_mixed(moe_kwargs, pop)

    def test_const_to_Moe(self, base_simulation_kwargs, moe_kwargs):
        # q >= 0.35 / m1 and P in [0.35, 6000] d cover Moe's q >= 0.1
        sim = base_simulation_kwargs.copy()
        sim['binary_fraction_const'] = 0.5
        sim['primary_mass_max'] = 40.0
        moe_kwargs['orbital_period_max'] = 5e3
        check_reweighting_mixed(sim, moe_kwargs)

    def test_Moe_to_const(self, moe_kwargs, base_simulation_kwargs):
        pop = base_simulation_kwargs.copy()
        pop['binary_fraction_const'] = 0.5
        pop['primary_mass_max'] = 40.0
        # q >= 4 / m1 >= 0.1 and P inside the Moe period range
        pop['secondary_mass_min'] = 4.0
        pop['orbital_period_min'] = 10**0.15
        pop['orbital_period_max'] = 6e3
        check_reweighting_mixed(moe_kwargs, pop)
