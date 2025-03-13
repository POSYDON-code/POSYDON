from pytest import fixture, mark
from posydon.popsyn.synthetic_population import Population, TransientPopulation, Rates
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from posydon.popsyn.norm_pop import calculate_model_weights
from posydon.popsyn.independent_sample import generate_independent_samples
import pytest


@fixture
def base_simulation_kwargs():
    simulation_params = {'number_of_binaries': 100000,
          'primary_mass_scheme':'Kroupa2001',
          'binary_fraction_const':1.0,
          'binary_fraction_scheme':'const',
          'primary_mass_min': 7,
          'primary_mass_max': 140,
          'secondary_mass_scheme':'flat_mass_ratio',
          'secondary_mass_min':0.35,
          'secondary_mass_max':140,
          'orbital_scheme':'period',
          'orbital_period_scheme':'Sana+12_period_extended',
          'orbital_period_min':0.35,
          'orbital_period_max':6e3,
          'eccentricity_scheme':'zero',}
    return simulation_params

@fixture
def base_population_kwargs(base_simulation_kwargs):
    population_params = base_simulation_kwargs.copy()
    population_params['q_min'] = 0.0
    population_params['q_max'] = 1.0
    return population_params


@fixture
def base_population_data():
    '''Example data that has q_min=0.05, and q_max=1 weights.'''
    return np.load('example_data/base_weights.npy')

def pop_data(kwargs):
    '''generate a population of stars'''
    sample = generate_independent_samples(**kwargs)
    pop_data = pd.DataFrame(np.array(sample).T, columns=['orbital_period_i', 'eccentricity_i', 'S1_mass_i', 'S2_mass_i',])
    pop_data['state_i'] = 'detached'
    f_b = kwargs['binary_fraction_const']
    n_binaries = int(f_b * len(pop_data))
    pop_data.loc[n_binaries:, 'S2_mass_i'] = 0
    pop_data.loc[n_binaries:, 'state_i'] = 'initially_single_star'
    return pop_data

class TestReweighting():
    
    def test_population_larger_sample_space(self, base_simulation_kwargs):
        
        # no binaries
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 0.

        old_max = 20
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = old_max
        small_sample = pop_data(base_simulation_kwargs)
        
        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['primary_mass_min'] = 7.
        expanded_kwargs['primary_mass_max'] = 40
        
        expanded_sample = pop_data(expanded_kwargs)
        
        M_sim = small_sample['S1_mass_i'].sum() + small_sample['S2_mass_i'].sum()
        small_weights = calculate_model_weights(small_sample, M_sim, base_simulation_kwargs, expanded_kwargs)
        
        M_sim = expanded_sample['S1_mass_i'].sum() + expanded_sample['S2_mass_i'].sum()
        expanded_weights = calculate_model_weights(expanded_sample, M_sim, expanded_kwargs, expanded_kwargs)
        
        assert len(expanded_weights)  == len(expanded_sample)
        
        mask = expanded_sample['S1_mass_i'] <= old_max
        selection = expanded_weights[mask]
        print(np.sum(selection))
        print(np.sum(small_weights))
        assert np.isclose(np.sum(selection), np.sum(small_weights), atol=1e-3)
        
    def test_population_larger_sample_space_binary(self, base_simulation_kwargs):
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 1.

        old_max = 20
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = old_max
        small_sample = pop_data(base_simulation_kwargs)
        
        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['primary_mass_min'] = 7.
        expanded_kwargs['primary_mass_max'] = 40
        
        expanded_sample = pop_data(expanded_kwargs)
        
        M_sim = small_sample['S1_mass_i'].sum() + small_sample['S2_mass_i'].sum()
        small_weights = calculate_model_weights(small_sample, M_sim, base_simulation_kwargs, expanded_kwargs)
        
        M_sim = expanded_sample['S1_mass_i'].sum() + expanded_sample['S2_mass_i'].sum()
        expanded_weights = calculate_model_weights(expanded_sample, M_sim, expanded_kwargs, expanded_kwargs)
        
        assert len(expanded_weights)  == len(expanded_sample)
        
        mask = expanded_sample['S1_mass_i'] <= old_max
        selection = expanded_weights[mask]
        print(np.sum(selection))
        print(np.sum(small_weights))
        assert np.isclose(np.sum(selection), np.sum(small_weights), atol=1e-3)
        
    def test_population_larger_sample_space_mixed(self, base_simulation_kwargs):
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 0.7

        old_max = 20
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = old_max
        small_sample = pop_data(base_simulation_kwargs)
        
        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['primary_mass_min'] = 7.
        expanded_kwargs['primary_mass_max'] = 40
        
        expanded_sample = pop_data(expanded_kwargs)
        
        M_sim = small_sample['S1_mass_i'].sum() + small_sample['S2_mass_i'].sum()
        small_weights = calculate_model_weights(small_sample, M_sim, base_simulation_kwargs, expanded_kwargs)
        
        M_sim = expanded_sample['S1_mass_i'].sum() + expanded_sample['S2_mass_i'].sum()
        expanded_weights = calculate_model_weights(expanded_sample, M_sim, expanded_kwargs, expanded_kwargs)
        
        assert len(expanded_weights)  == len(expanded_sample)
        
        mask = expanded_sample['S1_mass_i'] <= old_max
        selection = expanded_weights[mask]
        print(np.sum(selection))
        print(np.sum(small_weights))
        assert np.isclose(np.sum(selection), np.sum(small_weights), atol=1e-3)
        
    def test_population_binary_fraction_change(self, base_simulation_kwargs):
        
        # no binaries
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 1

        old_max = 20
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = old_max
        small_sample = pop_data(base_simulation_kwargs)
        
        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['primary_mass_min'] = 7.
        expanded_kwargs['primary_mass_max'] = 40
        expanded_kwargs['binary_fraction_const'] = 0.7
        
        expanded_sample = pop_data(expanded_kwargs)
        
        M_sim = small_sample['S1_mass_i'].sum() + small_sample['S2_mass_i'].sum()
        small_weights = calculate_model_weights(small_sample, M_sim, base_simulation_kwargs, expanded_kwargs)
        
        M_sim = expanded_sample['S1_mass_i'].sum() + expanded_sample['S2_mass_i'].sum()
        expanded_weights = calculate_model_weights(expanded_sample, M_sim, expanded_kwargs, expanded_kwargs)
        
        assert len(expanded_weights)  == len(expanded_sample)
        
        mask = (expanded_sample['S1_mass_i'] <= old_max) & (expanded_sample['state_i'] != 'initially_single_star')
        selection = expanded_weights[mask]
        print(np.sum(selection))
        print(np.sum(small_weights))
        assert np.isclose(np.sum(selection), np.sum(small_weights), atol=1e-3)
        
    def test_population_single_stars_binary_fraction(self, base_simulation_kwargs): 
        # no binaries
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 0

        old_max = 20
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = old_max
        small_sample = pop_data(base_simulation_kwargs)
        
        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['primary_mass_min'] = 7.
        expanded_kwargs['primary_mass_max'] = 20
        expanded_kwargs['binary_fraction_const'] = 0.7
        
        expanded_sample = pop_data(expanded_kwargs)
        
        M_sim = small_sample['S1_mass_i'].sum() + small_sample['S2_mass_i'].sum()
        small_weights = calculate_model_weights(small_sample, M_sim, base_simulation_kwargs, expanded_kwargs)
        
        M_sim = expanded_sample['S1_mass_i'].sum() + expanded_sample['S2_mass_i'].sum()
        expanded_weights = calculate_model_weights(expanded_sample, M_sim, expanded_kwargs, expanded_kwargs)
        
        assert len(expanded_weights)  == len(expanded_sample)
        
        mask = (expanded_sample['S1_mass_i'] <= old_max) & (expanded_sample['state_i'] == 'initially_single_star')
        selection = expanded_weights[mask]
        print(np.sum(selection))
        print(np.sum(small_weights))
        assert np.isclose(np.sum(selection), np.sum(small_weights), atol=1e-3)    
        
    def test_Salpeter_IMF_single(self, base_simulation_kwargs):
         # no binaries
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 0.
        base_simulation_kwargs['primary_mass_scheme'] = 'Salpeter'

        old_max = 20
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = old_max
        small_sample = pop_data(base_simulation_kwargs)
        
        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['primary_mass_min'] = 7.
        expanded_kwargs['primary_mass_max'] = 20
        expanded_kwargs['primary_mass_scheme'] = 'Salpeter'
        
        expanded_sample = pop_data(expanded_kwargs)
        
        M_sim = small_sample['S1_mass_i'].sum() + small_sample['S2_mass_i'].sum()
        small_weights = calculate_model_weights(small_sample, M_sim, base_simulation_kwargs, expanded_kwargs)
        
        M_sim = expanded_sample['S1_mass_i'].sum() + expanded_sample['S2_mass_i'].sum()
        expanded_weights = calculate_model_weights(expanded_sample, M_sim, expanded_kwargs, expanded_kwargs)
        
        assert len(expanded_weights)  == len(expanded_sample)
        
        mask = expanded_sample['S1_mass_i'] <= old_max
        selection = expanded_weights[mask]
        print(np.sum(selection))
        print(np.sum(small_weights))
        assert np.isclose(np.sum(selection), np.sum(small_weights), atol=1e-3)
        
    def test_changing_IMF_single(self, base_simulation_kwargs):
         # no binaries
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 0.
        base_simulation_kwargs['primary_mass_scheme'] = 'Salpeter'

        old_max = 20
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = old_max
        small_sample = pop_data(base_simulation_kwargs)
        
        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['primary_mass_min'] = 7.
        expanded_kwargs['primary_mass_max'] = 20
        expanded_kwargs['primary_mass_scheme'] = 'Kroupa2001'
        
        expanded_sample = pop_data(expanded_kwargs)
        
        M_sim = small_sample['S1_mass_i'].sum() + small_sample['S2_mass_i'].sum()
        small_weights = calculate_model_weights(small_sample, M_sim, base_simulation_kwargs, expanded_kwargs)
        
        M_sim = expanded_sample['S1_mass_i'].sum() + expanded_sample['S2_mass_i'].sum()
        expanded_weights = calculate_model_weights(expanded_sample, M_sim, expanded_kwargs, expanded_kwargs)
        
        assert len(expanded_weights)  == len(expanded_sample)
        
        mask = expanded_sample['S1_mass_i'] <= old_max
        selection = expanded_weights[mask]
        print(np.sum(selection))
        print(np.sum(small_weights))
        assert np.isclose(np.sum(selection), np.sum(small_weights), atol=1e-3)
        
    def test_changing_IMF_binary(self, base_simulation_kwargs):
         # no binaries
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 1
        base_simulation_kwargs['primary_mass_scheme'] = 'Salpeter'

        old_max = 20
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = old_max
        small_sample = pop_data(base_simulation_kwargs)
        
        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['primary_mass_min'] = 7.
        expanded_kwargs['primary_mass_max'] = 20
        expanded_kwargs['primary_mass_scheme'] = 'Kroupa2001'
        
        expanded_sample = pop_data(expanded_kwargs)
        
        M_sim = small_sample['S1_mass_i'].sum() + small_sample['S2_mass_i'].sum()
        small_weights = calculate_model_weights(small_sample, M_sim, base_simulation_kwargs, expanded_kwargs)
        
        M_sim = expanded_sample['S1_mass_i'].sum() + expanded_sample['S2_mass_i'].sum()
        expanded_weights = calculate_model_weights(expanded_sample, M_sim, expanded_kwargs, expanded_kwargs)
        
        assert len(expanded_weights)  == len(expanded_sample)
        
        mask = expanded_sample['S1_mass_i'] <= old_max
        selection = expanded_weights[mask]
        print(np.sum(selection))
        print(np.sum(small_weights))
        assert np.isclose(np.sum(selection), np.sum(small_weights), atol=1e-3)

    def test_changing_IMF_binary_range(self, base_simulation_kwargs):
        # no binaries
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 1
        base_simulation_kwargs['primary_mass_scheme'] = 'Salpeter'

        old_max = 20
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = old_max
        small_sample = pop_data(base_simulation_kwargs)
        
        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['primary_mass_min'] = 7.
        expanded_kwargs['primary_mass_max'] = 40
        expanded_kwargs['primary_mass_scheme'] = 'Kroupa2001'
        
        expanded_sample = pop_data(expanded_kwargs)
        
        M_sim = small_sample['S1_mass_i'].sum() + small_sample['S2_mass_i'].sum()
        small_weights = calculate_model_weights(small_sample, M_sim, base_simulation_kwargs, expanded_kwargs)
        
        M_sim = expanded_sample['S1_mass_i'].sum() + expanded_sample['S2_mass_i'].sum()
        expanded_weights = calculate_model_weights(expanded_sample, M_sim, expanded_kwargs, expanded_kwargs)
        
        assert len(expanded_weights)  == len(expanded_sample)
        
        mask = expanded_sample['S1_mass_i'] <= old_max
        selection = expanded_weights[mask]
        print(np.sum(selection))
        print(np.sum(small_weights))
        assert np.isclose(np.sum(selection), np.sum(small_weights), atol=1e-3)
    
    def test_same_sample_space(self, base_simulation_kwargs):
        '''Same sample space the weights should be the same as the simulation'''
        
        base_pop_data = pop_data(base_simulation_kwargs)
        M_sim = base_pop_data['S1_mass_i'].sum() + base_pop_data['S2_mass_i'].sum()
        weights = calculate_model_weights(base_pop_data, M_sim, base_simulation_kwargs, base_simulation_kwargs)
        assert len(weights) == len(base_pop_data)
        assert np.all(weights >= 0)
        assert np.allclose(weights, 1/M_sim)
        
    
    def test_mass_ratio_extension(self,
                                  base_simulation_kwargs,
                                  base_population_kwargs,
                                  base_population_data):
        '''Sample space for mass ratio extension'''
        
        base_pop_data = pop_data(base_simulation_kwargs)
        M_sim = base_pop_data['S1_mass_i'].sum() + base_pop_data['S2_mass_i'].sum()
        weights = calculate_model_weights(base_pop_data, M_sim, base_simulation_kwargs, base_population_kwargs)
        assert len(weights) == len(base_pop_data)
        assert np.all(weights >= 0)
        np.isclose(base_population_data, weights)
        
        

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
        
        M_sim = small_sample['S1_mass_i'].sum() + small_sample['S2_mass_i'].sum()
        small_weights = calculate_model_weights(small_sample, M_sim, base_simulation_kwargs, expanded_kwargs)
        
        M_sim = expanded_sample['S1_mass_i'].sum() + expanded_sample['S2_mass_i'].sum()
        expanded_weights = calculate_model_weights(expanded_sample, M_sim, expanded_kwargs, expanded_kwargs)
        
        assert len(expanded_weights)  == len(expanded_sample)
        
        # Check just single stars in the population.
        mask = (expanded_sample['S1_mass_i'] <= 20) & (expanded_sample['state_i'] == 'initially_single_star')
        selection = expanded_weights[mask]
        print(np.sum(selection))
        print(np.sum(small_weights))
        assert np.isclose(np.sum(selection), np.sum(small_weights), atol=1e-3)    
        
    def test_single_star_population_error(self, base_simulation_kwargs):
        # no binaries
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 0
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = 20
        
        small_sample = pop_data(base_simulation_kwargs)
        
        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['binary_fraction_const'] = 1.0
        
        M_sim = small_sample['S1_mass_i'].sum() + small_sample['S2_mass_i'].sum()
        with pytest.raises(ValueError) as e:
            small_weights = calculate_model_weights(small_sample, M_sim, base_simulation_kwargs, expanded_kwargs)
        assert 'No binaries simulated, but requested' in str(e.value)
        
    def test_binary_star_population_error(self, base_simulation_kwargs):
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 1.0
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = 20
        
        small_sample = pop_data(base_simulation_kwargs)
        
        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['binary_fraction_const'] = 0
        
        M_sim = small_sample['S1_mass_i'].sum() + small_sample['S2_mass_i'].sum()
        with pytest.raises(ValueError) as e:
            small_weights = calculate_model_weights(small_sample, M_sim, base_simulation_kwargs, expanded_kwargs)
        assert 'No single stars simulated, but requested' in str(e.value)
        
    
    def test_mixed_population(self, base_simulation_kwargs):
        # no binaries
        base_simulation_kwargs['number_of_binaries'] =  100000
        base_simulation_kwargs['binary_fraction_const'] = 0.3
        base_simulation_kwargs['primary_mass_min'] = 7.
        base_simulation_kwargs['primary_mass_max'] = 20
        
        small_sample = pop_data(base_simulation_kwargs)
        
        expanded_kwargs = base_simulation_kwargs.copy()
        expanded_kwargs['binary_fraction_const'] = 0.5
        expanded_sample = pop_data(expanded_kwargs)
        
        M_sim = small_sample['S1_mass_i'].sum() + small_sample['S2_mass_i'].sum()
        small_weights = calculate_model_weights(small_sample, M_sim, base_simulation_kwargs, expanded_kwargs)
        
        M_sim = expanded_sample['S1_mass_i'].sum() + expanded_sample['S2_mass_i'].sum()
        expanded_weights = calculate_model_weights(expanded_sample, M_sim, expanded_kwargs, expanded_kwargs)
        
        assert len(expanded_weights)  == len(expanded_sample)
        
        # Check just single stars in the population; need selection on both.
        mask = (expanded_sample['S1_mass_i'] <= 20) & (expanded_sample['state_i'] == 'initially_single_star')
        selection = expanded_weights[mask]
        print(np.sum(selection))
        
        mask2 = (small_sample['S1_mass_i'] <= 20) & (small_sample['state_i'] == 'initially_single_star')
        print(np.sum(small_weights[mask2]))
        assert np.isclose(np.sum(selection), np.sum(small_weights[mask2]), atol=1e-3)    
    

