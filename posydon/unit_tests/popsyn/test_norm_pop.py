import numpy as np
import os
import pandas as pd
import pytest
from pytest import fixture, mark

import posydon.popsyn.IMFs as IMFs
from posydon.popsyn.independent_sample import generate_independent_samples
from posydon.config import PATH_TO_POSYDON
from posydon.utils.posydonwarning import UnsupportedModelWarning

# file to test
from posydon.popsyn import norm_pop

# Dummy IMF to simulate a valid IMF pdf
class DummyIMF:
    def __init__(self, m_min, m_max):
        self.m_min = m_min
        self.m_max = m_max
    def pdf(self, m1):
        # returns a constant value for testing
        return 0.5

class TestGetIMFPdf:
    def test_invalid_imf_returns_constant(self):
        kwargs = {
            'primary_mass_scheme': 'NonExistentIMF',
            'primary_mass_min': 1,
            'primary_mass_max': 100
        }
        pdf_func = norm_pop.get_IMF_pdf(kwargs)
        m_test = np.array([1, 10, 50])
        result = pdf_func(m_test)
        assert np.allclose(result, np.ones_like(m_test))
    
    
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
        class SchemeKwargsIMF:
            def __init__(self, m_min, m_max, custom_param=None):
                self.m_min = m_min
                self.m_max = m_max
                self.custom_param = custom_param
            def pdf(self, m1):
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
        


class TestGetMassRatioPdf:
    def test_flat_mass_ratio_pdf(self):
        # kwargs without q_min and q_max to trigger flat_mass_ratio branch
        kwargs = {
            'secondary_mass_scheme': 'flat_mass_ratio',
            'secondary_mass_min': 5,
            'secondary_mass_max': 50
        }
        q_pdf = norm_pop.get_mass_ratio_pdf(kwargs)
        # For m1 = 10, minimum = max(5/10, 0.05) = 0.5 and
        # maximum = min(50/10, 1) = 1, so pdf = 1/(1-0.5) = 2.
        m1_val = 10
        # Test with q in range
        result = q_pdf(0.75, m1_val)
        assert np.allclose(result, 2)
        # Test with q below the minimum
        result = q_pdf(0.3, m1_val)
        assert np.all(result == 0)
        # Test with q above the maximum
        result = q_pdf(1.1, m1_val)
        assert np.all(result == 0)

    def test_alternative_mass_ratio_pdf(self):
        # kwargs with q_min and q_max to trigger the alternative branch
        kwargs = {
            'secondary_mass_scheme': 'flat_mass_ratio',
            'q_min': 0.2,
            'q_max': 0.8,
        }
        q_pdf = norm_pop.get_mass_ratio_pdf(kwargs)
        # In alternative branch, pdf returns 1 for q in (0,1] 
        result = q_pdf(0.5, None)
        assert np.all(result == 1/(kwargs['q_max']-kwargs['q_min']))
        result = q_pdf(0, None)
        assert np.all(result == 0)

    def test_invalid_mass_ratio_scheme(self):
        kwargs = {
            'secondary_mass_scheme': 'invalid'
        }
        
        with pytest.warns(UnsupportedModelWarning) as warning_info:
            q_pdf = norm_pop.get_mass_ratio_pdf(kwargs)
        
        assert ("The secondary_mass_scheme is not defined use a flat mass ratio"
                " distribution in (0,1]." in str(warning_info[0].message))
        results = q_pdf(0.4)
        assert np.all(results == 1)

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

class TestGetPdf:
    def test_single_star_pdf(self):
        # Using fallback IMF (non-existent) -> IMF_pdf returns 1.
        kwargs = {
            'primary_mass_scheme': 'NonExistentIMF',
            'primary_mass_min': 1,
            'primary_mass_max': 100,
            'secondary_mass_scheme': 'flat_mass_ratio',
            'secondary_mass_min': 1,
            'secondary_mass_max': 100,
            'binary_fraction_scheme': 'const',
            'binary_fraction_const': 0.5,
            'orbital_period_scheme': 'Sana+12_period_extended',
            'orbital_period_min': 0.35,
            'orbital_period_max': 6000,
        }
        pdf_func = norm_pop.get_pdf(kwargs)
        m1_val = 10
        # For non-binary, expected = (1 - f_b)*IMF_pdf = 0.5*1 = 0.5
        result = pdf_func(m1_val, binary=False)
        assert np.allclose(result, 0.5)
    
    def test_binary_pdf(self):
        # Using alternative mass ratio branch (q_min and q_max provided)
        kwargs = {
            'primary_mass_scheme': 'NonExistentIMF',
            'primary_mass_min': 1,
            'primary_mass_max': 100,
            'secondary_mass_scheme': 'flat_mass_ratio',
            'q_min': 0,
            'q_max': 1,
            'binary_fraction_scheme': 'const',
            'binary_fraction_const': 0.3,
            'orbital_period_scheme': 'Sana+12_period_extended',
            'orbital_period_min': 0.35,
            'orbital_period_max': 6000,
            
        }
        pdf_func = norm_pop.get_pdf(kwargs, mass_pdf=True)
        m1_val = 10
        # For binary, expected = f_b*IMF_pdf*q_pdf = 0.3*1*1 = 0.3 
        # for valid q (e.g. 0.5)
        result = pdf_func(m1_val, q=0.5, binary=True)
        assert np.allclose(result, 0.3)
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

# class TestGetMeanMass:
#     def test_dummy_mean_mass(self):
#         # dummy PDF returns 1 for any input
#         #dummy_pdf = lambda *args, **kwargs: 1
#         # parameters chosen for analytic integration:
#         # m1 in [1,2], q in [0.2, 0.8]
#         params = {
#             'primary_mass_min': 1,
#             'primary_mass_max': 2,
#             'q_min': 0.2,
#             'q_max': 0.8
#         }
#         # Binary integration:
#         # I_bin = ∫[m=1 to 2] m * (∫[q=0.2 to 0.8] (1+q)dq) dm
#         # ∫[q=0.2 to 0.8] (1+q)dq = [(q + 0.5*q^2)]_0.2^0.8 
#         # = (0.8+0.32) - (0.2+0.02) = 1.12-0.22 = 0.9
#         # I_bin = 0.9 * ∫[m=1 to 2] m dm = 0.9*( (2^2-1^2)/2 ) = 0.9*1.5 = 1.35
#         # Single star integration:
#         # I_single = ∫[m=1 to 2] m dm = 1.5
#         # Total expected mean mass = 1.35+1.5 = 2.85
#         result = norm_pop.get_mean_mass(params)
#         assert np.allclose(result, 2.85, rtol=1e-2)

    # def test_dummy_mean_mass_without_q_bounds(self):
    #     # dummy PDF returns 1 for any input
    #     dummy_pdf = lambda *args, **kwargs: 1
    #     # parameters without q_min and q_max, but with 
    #     #  secondary_mass_min and secondary_mass_max:
    #     # Let secondary_mass_min = 0.5, secondary_mass_max = 1.5.
    #     # Default q_min = max(secondary_mass_min/primary_mass_min, 0) 
    #     # = max(0.5/1,0)=0.5
    #     # Default q_max = min(secondary_mass_max/primary_mass_max, 1)
    #     # = min(1.5/2,1)=0.75
    #     # Binary integration:
    #     # Inner integration for each m: ∫[q=0.5 to 0.75] (1+q)dq
    #     # = [(q + 0.5*q^2)] from 0.5 to 0.75 = (0.75+0.28125) - (0.5+0.125) 
    #     # = 1.03125 - 0.625 = 0.40625
    #     # Outer integration: ∫[m=1 to 2] m*0.40625 dm = 0.40625*( (2^2-1^2)/2 ) 
    #     # = 0.40625*1.5 = 0.609375
    #     # Single star integration: ∫[m=1 to 2] m dm = 1.5
    #     # Total expected mean mass = 0.609375 + 1.5 = 2.109375
    #     params = {
    #         'primary_mass_min': 1,
    #         'primary_mass_max': 2,
    #         'secondary_mass_min': 0.5,
    #         'secondary_mass_max': 1.5
    #     }
    #     result = norm_pop.get_mean_mass(dummy_pdf, params)
    #     assert np.allclose(result, 2.109375, rtol=1e-2)

####
# Test the reweighting function
#### 


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
    return np.load(
        os.path.join(PATH_TO_POSYDON, 
                    'posydon/unit_tests/popsyn/example_data/base_weights.npy'))

def pop_data(kwargs):
    '''generate a population of stars'''
    sample = generate_independent_samples(**kwargs)
    pop_data = pd.DataFrame(np.array(sample).T,
                            columns=['orbital_period_i',
                                     'eccentricity_i',
                                     'S1_mass_i',
                                     'S2_mass_i',])
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
        
        mask = expanded_sample['S1_mass_i'] <= old_max
        selection = expanded_weights[mask]
        print(np.sum(selection))
        print(np.sum(small_weights))
        assert np.isclose(np.sum(selection), np.sum(small_weights), atol=1e-3)
        
    def test_population_larger_sample_space_binary_fraction_change(self,
                                                        base_simulation_kwargs):
        
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
        
        mask = ((expanded_sample['S1_mass_i'] <= old_max) 
                & (expanded_sample['state_i'] != 'initially_single_star'))
        selection = expanded_weights[mask]
        print(np.sum(selection))
        print(np.sum(small_weights))
        assert np.isclose(np.sum(selection), np.sum(small_weights), atol=1e-3)
        
    def test_population_single_stars_binary_fraction(self, 
                                                    base_simulation_kwargs): 
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
        
        mask = ((expanded_sample['S1_mass_i'] <= old_max)
                & (expanded_sample['state_i'] == 'initially_single_star'))
        selection = expanded_weights[mask]
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
        
        mask = expanded_sample['S1_mass_i'] <= old_max
        selection = expanded_weights[mask]
        print(np.sum(selection))
        print(np.sum(small_weights))
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
        
    
    def test_mass_ratio_extension(self,
                                  base_simulation_kwargs,
                                  base_population_kwargs,
                                  base_population_data):
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
        mask = ((expanded_sample['S1_mass_i'] <= 20)
                & (expanded_sample['state_i'] == 'initially_single_star'))
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
        # no binaries
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
        print(np.sum(selection))
        
        mask2 = ((small_sample['S1_mass_i'] <= 20)
                 & (small_sample['state_i'] == 'initially_single_star'))
        print(np.sum(small_weights[mask2]))
        assert np.isclose(np.sum(selection),
                          np.sum(small_weights[mask2]),
                          atol=1e-3)


