import numpy as np
import pytest
from posydon.popsyn.star_formation_history import SFHBase, MadauBase
from posydon.popsyn.star_formation_history import (
    MadauDickinson14,
    MadauFragos17,
    Neijssel19,
    Fujimoto24,
    IllustrisTNG,
    Chruslinska21,
    Zavala21,
    get_SFH_model
)


class TestSFHBase:
    def test_init_attributes(self):
        """Test that the initialization sets attributes correctly."""
        # Create a concrete subclass for testing
        class ConcreteSFH(SFHBase):
            def CSFRD(self, z):
                return z
            
            def mean_metallicity(self, z):
                return z
                
            def fSFR(self, z, metallicity_bins):
                return np.ones((len(z), len(metallicity_bins)-1))
        
        model_dict = {
            "test_param": 42, 
            "Z_max": 0.03,
            "another_param": "test"
        }
        sfh = ConcreteSFH(model_dict)
        
        # Check that attributes are set correctly
        assert sfh.Z_max == 0.03
        assert sfh.test_param == 42
        assert sfh.another_param == "test"
        assert sfh.MODEL == model_dict
    
    def test_validation(self):
        """Test that Z_max > 1 raises a ValueError."""
        class ConcreteSFH(SFHBase):
            def CSFRD(self, z):
                pass
            def mean_metallicity(self, z):
                pass
            def fSFR(self, z, metallicity_bins):
                pass
                
        model_dict = {"Z_max": 1.5}
        with pytest.raises(ValueError) as excinfo:
            ConcreteSFH(model_dict)
        assert "Z_max must be in absolute units!" in str(excinfo.value)
        
        # Test with Z_min > Z_max
        model_dict = {"Z_max": 0.1, "Z_min": 0.2}
        with pytest.raises(ValueError) as excinfo:
            sfh = ConcreteSFH(model_dict)
        assert "Z_min must be smaller than Z_max!" in str(excinfo.value)
    
        # Test with Z_min < 0
        model_dict = {"Z_max": 0.1, "Z_min": -0.1}
        with pytest.raises(ValueError) as excinfo:
            sfh = ConcreteSFH(model_dict)
        assert "Z_min must be in absolute units!" in str(excinfo.value)
    
    def test_abstract_methods(self):
        """Test that abstract methods must be implemented."""
        # Create incomplete subclasses that don't implement all abstract methods
        class IncompleteSFH1(SFHBase):
            def CSFRD(self, z):
                return z
                
        class IncompleteSFH2(SFHBase):
            def CSFRD(self, z):
                return z
            def mean_metallicity(self, z):
                return z
        
        model_dict = {"Z_max": 0.03}
        with pytest.raises(TypeError):
            IncompleteSFH1(model_dict)
        
        with pytest.raises(TypeError):
            IncompleteSFH2(model_dict)
    
    def test_distribute_cdf(self):
        """Test the _distribute_cdf method."""
        class ConcreteSFH(SFHBase):
            def CSFRD(self, z):
                return z
            def mean_metallicity(self, z):
                return z
            def fSFR(self, z, metallicity_bins):
                return np.ones((len(z), len(metallicity_bins)-1))
        
        model_dict = {"Z_max": 1, "Z_min": 0.0}
        sfh = ConcreteSFH(model_dict)
        
        # Create a simple CDF functions
        cdf_func = lambda x: x
        met_edges = np.array([0.0, 0.01, 0.02, 0.03])
        
        result = sfh._distribute_cdf(cdf_func, met_edges)
        expected = np.array([0.01, 0.01, 0.98])
        np.testing.assert_allclose(result, expected)
        
        # Test with normalization
        sfh.normalise = True
        result = sfh._distribute_cdf(cdf_func, met_edges)
        np.testing.assert_allclose(np.sum(result), 1.0)
        
        # Test with different model dicts
        model_dict = {"Z_max": 1, "Z_min": 0.015}
        sfh = ConcreteSFH(model_dict)
        met_edges = np.array([0.3, 0.6, 0.9])
        result = sfh._distribute_cdf(cdf_func, met_edges)
        expected = np.array([0.585, 0.4])
        np.testing.assert_allclose(result, expected)
        
        # Test with normalise
        sfh.normalise = True
        result = sfh._distribute_cdf(cdf_func, met_edges)
        np.testing.assert_allclose(np.sum(result), 1.0)
        
        # restrict the upper bound
        model_dict = {"Z_max": 0.95, "Z_min": 0.2}
        sfh = ConcreteSFH(model_dict)
        met_edges = np.array([0.3, 0.6, 0.9])
        result = sfh._distribute_cdf(cdf_func, met_edges)
        expected = np.array([0.4, 0.35])
        np.testing.assert_allclose(result, expected)
        
        # Test with normalise
        sfh.normalise = True
        result = sfh._distribute_cdf(cdf_func, met_edges)
        np.testing.assert_allclose(np.sum(result), 1.0)
    
    
    def test_call_method(self):
        """Test the __call__ method."""
        class ConcreteSFH(SFHBase):
            def CSFRD(self, z):
                return np.array([1.0, 2.0])
            
            def mean_metallicity(self, z):
                return np.array([0.01, 0.02])
                
            def fSFR(self, z, metallicity_bins):
                # Return a simple array for testing
                return np.array([[0.3, 0.7], [0.4, 0.6]])
        
        model_dict = {"Z_max": 0.03}
        sfh = ConcreteSFH(model_dict)
        
        z = np.array([0.5, 1.0])
        met_edges = np.array([0.0, 0.01, 0.03])
        
        result = sfh(z, met_edges)
        
        # Expected: CSFRD(z)[:, np.newaxis] * fSFR(z, met_edges)
        expected = np.array([
            [1.0 * 0.3, 1.0 * 0.7],
            [2.0 * 0.4, 2.0 * 0.6]
        ])
        
        np.testing.assert_allclose(result, expected)


class TestMadauBase:
    """Test class for MadauBase"""
    
    class ConcreteMadau(MadauBase):
        """Concrete subclass of MadauBase for testing"""
        def __init__(self, MODEL):
            super().__init__(MODEL)
            self.CSFRD_params = {
                "a": 0.01,
                "b": 2.6, 
                "c": 3.2,
                "d": 6.2
            }
    
    def test_init_requires_sigma(self):
        """Test that MadauBase requires a sigma parameter"""
        model_dict = {"Z_max": 0.03}
        with pytest.raises(ValueError) as excinfo:
            self.ConcreteMadau(model_dict)
        assert "sigma not given!" in str(excinfo.value)
    
    def test_init_sets_csfrd_params_to_none(self):
        """Test that CSFRD_params is set to None initially and can be set by subclass"""
        model_dict = {"sigma": 0.5, "Z_max": 0.03}
        madau = self.ConcreteMadau(model_dict)
        assert madau.CSFRD_params is not None
        assert madau.CSFRD_params["a"] == 0.01
        assert madau.CSFRD_params["b"] == 2.6
        assert madau.CSFRD_params["c"] == 3.2
        assert madau.CSFRD_params["d"] == 6.2
    
    def test_std_log_metallicity_dist(self):
        """Test the std_log_metallicity_dist method with different sigma values"""
        # Test Bavera+20
        model_dict = {"sigma": "Bavera+20", "Z_max": 0.03}
        madau = self.ConcreteMadau(model_dict)
        assert madau.std_log_metallicity_dist() == 0.5
        
        # Test Neijssel+19
        model_dict = {"sigma": "Neijssel+19", "Z_max": 0.03}
        madau = self.ConcreteMadau(model_dict)
        assert madau.std_log_metallicity_dist() == 0.39
        
        # Test float value
        model_dict = {"sigma": 0.45, "Z_max": 0.03}
        madau = self.ConcreteMadau(model_dict)
        assert madau.std_log_metallicity_dist() == 0.45
        
        # Test unknown string
        model_dict = {"sigma": "unknown", "Z_max": 0.03}
        madau = self.ConcreteMadau(model_dict)
        with pytest.raises(ValueError) as excinfo:
            madau.std_log_metallicity_dist()
        assert "Unknown sigma choice!" in str(excinfo.value)
        
        # Test invalid type
        model_dict = {"sigma": 1, "Z_max": 0.03}
        madau = self.ConcreteMadau(model_dict)
        with pytest.raises(ValueError) as excinfo:
            madau.std_log_metallicity_dist()
        assert "Invalid sigma value" in str(excinfo.value)
    
    def test_csfrd(self):
        """Test the CSFRD method"""
        model_dict = {"sigma": 0.5}
        madau = self.ConcreteMadau(model_dict)
        
        # Test with single value
        z = 0.0
        result = madau.CSFRD(z)
        expected = 0.01 * 1**2.6 / (1 + (1/3.2)**6.2)  # a * (1+z)^b / (1 + ((1+z)/c)^d) with z=0
        np.testing.assert_allclose(result, expected)
        
        # Test with array of values
        z_array = np.array([0.0, 1.0, 2.0])
        result = madau.CSFRD(z_array)
        expected = np.array([
            0.01 * 1**2.6 / (1 + (1/3.2)**6.2),  # z=0 
            0.01 * 2**2.6 / (1 + (2/3.2)**6.2),  # z=1
            0.01 * 3**2.6 / (1 + (3/3.2)**6.2)   # z=2
        ])
        np.testing.assert_allclose(result, expected)
    
    def test_mean_metallicity(self):
        """Test the mean_metallicity method"""
        from posydon.utils.constants import Zsun
        
        model_dict = {"sigma": 0.5, "Z_max": 0.03}
        madau = self.ConcreteMadau(model_dict)
        
        # Test with single value
        z = 0.0
        result = madau.mean_metallicity(z)
        expected = 10**(0.153) * Zsun
        np.testing.assert_allclose(result, expected)
        
        # Test with array of values
        z_array = np.array([0.0, 1.0, 2.0])
        result = madau.mean_metallicity(z_array)
        expected = 10**(0.153 - 0.074 * z_array**1.34) * Zsun
        np.testing.assert_allclose(result, expected)
    
    def test_fsfr(self):
        """Test the fSFR method"""
        model_dict = {"sigma": 0.5, "Z_max": 1}
        madau = self.ConcreteMadau(model_dict)
        
        # Test with redshift array and metallicity bins
        z = np.array([0.0, 1.0])
        met_bins = np.array([0.001, 0.01, 0.02, 0.03])
        
        result = madau.fSFR(z, met_bins)
        
        # Shape check - should be (len(z), len(met_bins)-1)
        assert result.shape == (2, 3)
        
        # THIS IS A VALIDATION TEST
        from scipy.stats import norm
        mean0, mean1 = (np.log10(madau.mean_metallicity(z)) 
                  - model_dict['sigma']**2 * np.log(10) / 2)
        print(mean0, mean1)
        
        expected1 = np.array(
            # integral from 0.001 to 0.01; Z_min = lowest bin edge
            [norm.cdf(np.log10(0.01), mean0, model_dict['sigma'])
             - norm.cdf(np.log10(0.001), mean0, model_dict['sigma']),
            # integral from 0.01 to 0.02
             (norm.cdf(np.log10(0.02), mean0, model_dict['sigma'])
              - norm.cdf(np.log10(0.01), mean0, model_dict['sigma'])),
            # integral from 0.02 to 1
             (norm.cdf(np.log10(1), mean0, model_dict['sigma'])
              - norm.cdf(np.log10(0.02), mean0, model_dict['sigma']))],)
        
        expected2 = np.array(
            # integral from 0.001 to 0.01;Z_min = lowest bin edge
            [norm.cdf(np.log10(0.01), mean1, model_dict['sigma'])
             - norm.cdf(np.log10(0.001), mean1, model_dict['sigma']),
            # integral from 0.01 to 0.02
             (norm.cdf(np.log10(0.02), mean1, model_dict['sigma'])
              - norm.cdf(np.log10(0.01), mean1, model_dict['sigma'])),
            # integral from 0.02 to 1
             (norm.cdf(np.log10(1), mean1, model_dict['sigma'])
              - norm.cdf(np.log10(0.02), mean1, model_dict['sigma']))],
                      )
        expected = np.array([expected1, expected2])
        np.testing.assert_allclose(result, expected)
        
        # Change Z_min to 0 to include the rest of the lowest mets
        model_dict = {"sigma": 0.5, "Z_max": 0.3, "Z_min": 0}
        madau = self.ConcreteMadau(model_dict)
        result = madau.fSFR(z, met_bins)
        
        expected1 = np.array(
            [norm.cdf(np.log10(0.01), mean0, model_dict['sigma']),
                (norm.cdf(np.log10(0.02), mean0, model_dict['sigma'])
                - norm.cdf(np.log10(0.01), mean0, model_dict['sigma'])),
                (norm.cdf(np.log10(0.3), mean0, model_dict['sigma'])
                - norm.cdf(np.log10(0.02), mean0, model_dict['sigma']))],)
        expected2 = np.array(
            [norm.cdf(np.log10(0.01), mean1, model_dict['sigma']),
                (norm.cdf(np.log10(0.02), mean1, model_dict['sigma'])
                - norm.cdf(np.log10(0.01), mean1, model_dict['sigma'])),
                (norm.cdf(np.log10(0.3), mean1, model_dict['sigma'])
                - norm.cdf(np.log10(0.02), mean1, model_dict['sigma']))],)
        expected = np.array([expected1, expected2])
        np.testing.assert_allclose(result, expected)
        
        # Test with normalise
        model_dict = {"sigma": 0.5, "Z_max": 0.3, "Z_min": 0, "normalise": True}
        madau = self.ConcreteMadau(model_dict)
        result = madau.fSFR(z, met_bins)
        print(result)
        np.testing.assert_allclose(np.sum(result, axis=1), np.ones(2))
        
        # Test with Z_min > 0.03
        model_dict = {"sigma": 0.5,
                      "Z_max": 0.3,
                      "Z_min": 0.02,
                      "normalise": True}
        madau = self.ConcreteMadau(model_dict)
        result = madau.fSFR(z, met_bins)
        expected = np.ones(2)
        np.testing.assert_allclose(np.sum(result, axis=1), expected)

class TestIllustrisTNG:
    """Tests for the IllustrisTNG SFH model with mocked data loading."""
    
    @pytest.fixture
    def mock_illustris_data(self):
        """Create mock data for the IllustrisTNG class."""
        # Create mock data structure similar to the npz file
        num_redshifts = 10
        num_metallicities = 5
        
        mock_data = {
            "SFR": np.linspace(0.1, 1.0, num_redshifts)[::-1],  # SFR decreases with redshift
            "redshifts": np.linspace(0.0, 9.0, num_redshifts)[::-1],  # Redshifts from 0 to 9
            "mets": np.logspace(-4, -1, num_metallicities),  # Metallicities from 1e-4 to 1e-1
            "M": np.ones((num_redshifts, num_metallicities))  # Equal mass in all bins for simplicity
        }
        
        # Add some variation to mass distribution for testing mean_metallicity
        for i in range(num_redshifts):
            # Linear decrease in higher metallicities as redshift increases
            scale = 1.0 - i / num_redshifts
            mock_data["M"][i] = np.linspace(1.0, scale, num_metallicities)
            
        mock_data["M"] = np.flip(mock_data["M"], axis=0)  # Reverse the mass array
        return mock_data
    
    @pytest.fixture
    def illustris_model(self, monkeypatch, mock_illustris_data):
        """Create an IllustrisTNG model instance with mocked data."""
        # Create a function that returns the mock data
        def mock_get_illustrisTNG_data(self, verbose=False):
            return mock_illustris_data
        
        # Patch the _get_illustrisTNG_data method
        monkeypatch.setattr(IllustrisTNG, "_get_illustrisTNG_data", mock_get_illustrisTNG_data)
        
        # Create and return the model
        model_dict = {"Z_max": 0.3}
        return IllustrisTNG(model_dict)
    
    def test_init_parameters(self, illustris_model, mock_illustris_data):
        """Test that initialization sets the parameters correctly."""
        # Check that data was loaded correctly
        np.testing.assert_array_equal(illustris_model.CSFRD_data, np.flip(mock_illustris_data["SFR"]))
        np.testing.assert_array_equal(illustris_model.redshifts, np.flip(mock_illustris_data["redshifts"]))
        np.testing.assert_array_equal(illustris_model.Z, mock_illustris_data["mets"])
        np.testing.assert_array_equal(illustris_model.M, np.flip(mock_illustris_data["M"], axis=0))
        
        # Check that model parameters were set correctly
        assert illustris_model.Z_max == 0.3
    
    def test_csfrd_calculation(self, illustris_model, mock_illustris_data):
        """Test the CSFRD method."""
        # Test at specific redshifts including boundary values
        z_values = np.array([0.0, 4.5, 9.0])
        result = illustris_model.CSFRD(z_values)
        
        # Expected values come from interpolating flipped SFR data
        flipped_sfr = np.flip(mock_illustris_data["SFR"])
        flipped_redshifts = np.flip(mock_illustris_data["redshifts"])
        expected = np.interp(z_values, flipped_redshifts, flipped_sfr)
        
        np.testing.assert_allclose(result, expected)
    
    def test_mean_metallicity(self, illustris_model, mock_illustris_data):
        """Test the mean_metallicity method."""
        # Test at specific redshifts
        z_values = np.array([0.0, 4.5, 9.0])
        result = illustris_model.mean_metallicity(z_values)
        
        # Calculate expected values manually
        expected = np.zeros(len(z_values))
        flipped_redshifts = np.flip(mock_illustris_data["redshifts"])
        flipped_masses = np.flip(mock_illustris_data["M"], axis=0)
        metallicities = mock_illustris_data["mets"]
        
        # Calculate expected mean metallicities at each test redshift
        out = np.zeros_like(flipped_redshifts)
        for i, z in enumerate(out):
            weights = flipped_masses[i, :]
            if np.sum(weights) == 0:
                out[i] = np.nan
            else:
                out[i] = np.average(metallicities, weights=weights)
        Z_interp = np.interp(z_values, flipped_redshifts, out)
        np.testing.assert_allclose(result, Z_interp)
    
    def test_fsfr_calculation(self, illustris_model):
        """Test the fSFR method."""
        # Test with redshift array and metallicity bins
        z = np.array([0.0, 4.5])
        met_bins = np.array([0.001, 0.01, 0.05, 0.1])
        
        result = illustris_model.fSFR(z, met_bins)
        
        # Shape check - should be (len(z), len(met_bins)-1)
        assert result.shape == (2, 3)
        
        # Test with normalise=True
        illustris_model.normalise = True
        result = illustris_model.fSFR(z, met_bins)
        for row in result:
            if np.sum(row) > 0:
                np.testing.assert_allclose(np.sum(row), 1.0)

class TestMadauDickinson14:
    """Tests for the MadauDickinson14 SFH model"""
    
    def test_init_parameters(self):
        """Test that initialization sets the correct CSFRD parameters"""
        model_dict = {"sigma": 0.5, "Z_max": 0.03}
        madau = MadauDickinson14(model_dict)
        
        # Check that CSFRD_params were set correctly
        assert madau.CSFRD_params["a"] == 0.015
        assert madau.CSFRD_params["b"] == 2.7
        assert madau.CSFRD_params["c"] == 2.9
        assert madau.CSFRD_params["d"] == 5.6
        
        # Check that it inherits correctly from MadauBase
        assert isinstance(madau, MadauBase)
    
    def test_csfrd_calculation(self):
        """Test that CSFRD calculations match expected values"""
        model_dict = {"sigma": 0.5, "Z_max": 0.03}
        madau = MadauDickinson14(model_dict)
        
        # Test at specific redshifts
        z_values = np.array([0.0, 1.0, 2.0, 6.0])
        result = madau.CSFRD(z_values)
        
        # Calculate expected values manually
        p = madau.CSFRD_params
        expected = p["a"] * (1.0 + z_values) ** p["b"] / (1.0 + ((1.0 + z_values) / p["c"]) ** p["d"])
        
        np.testing.assert_allclose(result, expected)
    
class TestMadauFragos17:
    """Tests for the MadauFragos17 SFH model"""
    
    def test_init_parameters(self):
        """Test that initialization sets the correct CSFRD parameters"""
        model_dict = {"sigma": 0.5, "Z_max": 0.03}
        madau = MadauFragos17(model_dict)
        
        # Check that CSFRD_params were set correctly
        assert madau.CSFRD_params["a"] == 0.01
        assert madau.CSFRD_params["b"] == 2.6
        assert madau.CSFRD_params["c"] == 3.2
        assert madau.CSFRD_params["d"] == 6.2
        
        # Check that it inherits correctly from MadauBase
        assert isinstance(madau, MadauBase)
    
    def test_csfrd_calculation(self):
        """Test that CSFRD calculations match expected values"""
        model_dict = {"sigma": 0.5, "Z_max": 0.03}
        madau = MadauFragos17(model_dict)
        
        # Test at specific redshifts
        z_values = np.array([0.0, 1.0, 2.0, 6.0])
        result = madau.CSFRD(z_values)
        
        # Calculate expected values manually using the formula
        p = madau.CSFRD_params
        expected = p["a"] * (1.0 + z_values) ** p["b"] / (1.0 + ((1.0 + z_values) / p["c"]) ** p["d"])
        
        np.testing.assert_allclose(result, expected)

class TestNeijssel19:
    """Tests for the Neijssel19 SFH model"""
    
    def test_init_parameters(self):
        """Test that initialization sets the correct CSFRD parameters"""
        model_dict = {"sigma": 0.5, "Z_max": 0.03}
        neijssel = Neijssel19(model_dict)
        
        # Check that CSFRD_params were set correctly
        assert neijssel.CSFRD_params["a"] == 0.01
        assert neijssel.CSFRD_params["b"] == 2.77
        assert neijssel.CSFRD_params["c"] == 2.9
        assert neijssel.CSFRD_params["d"] == 4.7
        
        # Check that it inherits correctly from MadauBase
        assert isinstance(neijssel, MadauBase)
    
    def test_csfrd_calculation(self):
        """Test that CSFRD calculations match expected values"""
        model_dict = {"sigma": 0.5, "Z_max": 0.03}
        neijssel = Neijssel19(model_dict)
        
        # Test at specific redshifts
        z_values = np.array([0.0, 1.0, 2.0, 6.0])
        result = neijssel.CSFRD(z_values)
        
        # Calculate expected values manually using the formula
        p = neijssel.CSFRD_params
        expected = p["a"] * (1.0 + z_values) ** p["b"] / (1.0 + ((1.0 + z_values) / p["c"]) ** p["d"])
        
        np.testing.assert_allclose(result, expected)
    
    def test_mean_metallicity(self):
        """Test the overridden mean_metallicity method"""
        model_dict = {"sigma": 0.5, "Z_max": 0.03}
        neijssel = Neijssel19(model_dict)
        
        # Test at specific redshifts
        z_values = np.array([0.0, 1.0, 2.0, 6.0])
        result = neijssel.mean_metallicity(z_values)
        
        # Calculate expected values based on Neijssel19's formula
        expected = 0.035 * 10 ** (-0.23 * z_values)
        
        np.testing.assert_allclose(result, expected)
    
    def test_fsfr_with_lognormal(self):
        """Test the overridden fSFR method which uses a ln-normal distribution"""
        model_dict = {"sigma": 0.5, "Z_max": 0.3}
        neijssel = Neijssel19(model_dict)
        
        # Test with redshift array and metallicity bins
        z = np.array([0.0, 1.0])
        met_bins = np.array([0.001, 0.01, 0.02, 0.03])
        
        result = neijssel.fSFR(z, met_bins)
        
        # Shape check - should be (len(z), len(met_bins)-1)
        assert result.shape == (2, 3)
        
        # Test normalization with normalise=True
        model_dict = {"sigma": 0.5, "Z_max": 0.3, "normalise": True}
        neijssel = Neijssel19(model_dict)
        result = neijssel.fSFR(z, met_bins)
        np.testing.assert_allclose(np.sum(result, axis=1), np.ones(2))

class TestFujimoto24:
    """Tests for the Fujimoto24 SFH model"""
    
    def test_init_parameters(self):
        """Test that initialization sets the correct CSFRD parameters"""
        model_dict = {"sigma": 0.5, "Z_max": 0.03}
        fujimoto = Fujimoto24(model_dict)
        
        # Check that CSFRD_params were set correctly
        assert fujimoto.CSFRD_params["a"] == 0.010
        assert fujimoto.CSFRD_params["b"] == 2.8
        assert fujimoto.CSFRD_params["c"] == 3.3
        assert fujimoto.CSFRD_params["d"] == 6.6
        
        # Check that it inherits correctly from MadauBase
        assert isinstance(fujimoto, MadauBase)
    
    def test_csfrd_calculation(self):
        """Test that CSFRD calculations match expected values"""
        model_dict = {"sigma": 0.5, "Z_max": 0.03}
        fujimoto = Fujimoto24(model_dict)
        
        # Test at specific redshifts
        z_values = np.array([0.0, 1.0, 2.0, 6.0])
        result = fujimoto.CSFRD(z_values)
        
        # Calculate expected values manually using the formula
        p = fujimoto.CSFRD_params
        expected = p["a"] * (1.0 + z_values) ** p["b"] / (1.0 + ((1.0 + z_values) / p["c"]) ** p["d"])
        
        np.testing.assert_allclose(result, expected)

class TestChruslinska21:
    """Tests for the Chruslinska21 SFH model with mocked data loading."""
    
    @pytest.fixture
    def mock_chruslinska_data(self, monkeypatch):
        """Create mock data for the Chruslinska21 class."""
        # Create mock data for FOH bins
        FOH_bins = np.linspace(5.3, 9.7, 200)
        dFOH = FOH_bins[1] - FOH_bins[0]
        redshifts = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
        delta_T = np.array([1e9, 1e9, 1e9, 1e9, 1e9, 1e9])  # Time bin widths
        
        # Mock SFR data - decreasing with redshift, varying with metallicity
        SFR_data = np.zeros((len(redshifts), len(FOH_bins)))
        for i in range(len(redshifts)):
            # Simple pattern: peak at middle metallicity, decreasing with redshift
            peak_idx = len(FOH_bins) // 2
            SFR_data[i] = np.exp(-0.5 * ((np.arange(len(FOH_bins)) - peak_idx) / 20)**2)
            SFR_data[i] *= np.exp(-redshifts[i] / 2)  # Decrease with redshift
        
        # Create method that returns tuple of time, redshift, deltaT
        def mock_load_redshift_data(self, verbose=False):
            time = np.array([1e9, 2e9, 3e9, 4e9, 5e9, 6e9])  # Fake times
            return time, redshifts, delta_T
        
        # Create method that returns the mock SFR data
        def mock_load_raw_data(self):
            return SFR_data * 1e6 * delta_T[:, np.newaxis] 
        
        # Patch the methods
        monkeypatch.setattr(Chruslinska21, "_load_redshift_data", mock_load_redshift_data)
        monkeypatch.setattr(Chruslinska21, "_load_raw_data", mock_load_raw_data)
        
        return {
            "FOH_bins": FOH_bins,
            "dFOH": dFOH,
            "redshifts": redshifts,
            "SFR_data": SFR_data
        }
    
    @pytest.fixture
    def chruslinska_model(self, mock_chruslinska_data, monkeypatch):
        """Create a Chruslinska21 model instance with mocked data."""
        model_dict = {
            "sub_model": "test_model",
            "Z_solar_scaling": "Asplund09",
            "Z_max": 0.03,
            "select_one_met": False
        }
        return Chruslinska21(model_dict)
    
    def test_init_parameters(self):
        """Test that initialization validates required parameters."""
        # Test missing sub_model
        with pytest.raises(ValueError) as excinfo:
            Chruslinska21({"Z_solar_scaling": "Asplund09", "Z_max": 0.03, "select_one_met": False})
        assert "Sub-model not given!" in str(excinfo.value)
        
        # Test missing Z_solar_scaling
        with pytest.raises(ValueError) as excinfo:
            Chruslinska21({"sub_model": "test", "Z_max": 0.03, "select_one_met": False})
        assert "Z_solar_scaling not given!" in str(excinfo.value)
        
        # Test missing Z_max
        with pytest.raises(ValueError) as excinfo:
            Chruslinska21({"sub_model": "test", "Z_solar_scaling": "Asplund09", "select_one_met": False})
        assert "Z_max not given!" in str(excinfo.value)
    
    def test_foh_to_z_conversion(self, chruslinska_model):
        """Test the _FOH_to_Z method for all scaling options."""
        # Test Asplund09 scaling
        FOH_test = np.array([7.0, 8.0, 8.69, 9.0])
        result = chruslinska_model._FOH_to_Z(FOH_test)
        
        # Expected: 10^(log10(0.0134) + FOH - 8.69)
        Zsun = 0.0134
        FOHsun = 8.69
        expected = 10**(np.log10(Zsun) + FOH_test - FOHsun)
        np.testing.assert_allclose(result, expected)
        
        # Test other scaling options
        model_dict = {
            "sub_model": "test_model",
            "Z_solar_scaling": "AndersGrevesse89",
            "Z_max": 0.03,
            "select_one_met": False
        }
        model = Chruslinska21(model_dict)
        result = model._FOH_to_Z(FOH_test)
        expected = 10**(np.log10(0.017) + FOH_test - 8.83)
        np.testing.assert_allclose(result, expected)
        
        # Test GrevesseSauval98 scaling
        model_dict["Z_solar_scaling"] = "GrevesseSauval98"
        model = Chruslinska21(model_dict)
        result = model._FOH_to_Z(FOH_test)
        expected = 10**(np.log10(0.0201) + FOH_test - 8.93)
        np.testing.assert_allclose(result, expected)
        
        # Test Villante14 scaling
        model_dict["Z_solar_scaling"] = "Villante14"
        model = Chruslinska21(model_dict)
        result = model._FOH_to_Z(FOH_test)
        expected = 10**(np.log10(0.019) + FOH_test - 8.85)
        np.testing.assert_allclose(result, expected)
        
        # Test invalid scaling
        model_dict["Z_solar_scaling"] = "InvalidScaling"
        with pytest.raises(ValueError) as excinfo:
            model = Chruslinska21(model_dict)
        assert "Invalid Z_solar_scaling!" in str(excinfo.value)
    
    def test_mean_metallicity(self, chruslinska_model, mock_chruslinska_data):
        """Test the mean_metallicity method."""
        # Test at specific redshifts
        z_values = np.array([0.0, 2.0, 4.0])
        result = chruslinska_model.mean_metallicity(z_values)
        print(result)
        # Should be an array of the same length as z_values
        assert len(result) == len(z_values)
        # Should be the same value at all redshifts
        assert np.isclose(result[0], result[1])
        assert np.isclose(result[0], result[2])
        assert np.isclose(result[0], 0.0014903210118641882)
    
    def test_csfrd_calculation(self, chruslinska_model, mock_chruslinska_data):
        """Test the CSFRD method."""
        # Test at specific redshifts
        z_values = np.array([0.0, 2.0, 4.0])
        result = chruslinska_model.CSFRD(z_values)
        
        # Should be an array of the same length as z_values
        assert len(result) == len(z_values)
        
        # Should decrease with increasing redshift
        assert result[0] > result[1] > result[2]
    
    def test_fsfr_calculation(self, chruslinska_model):
        """Test the fSFR method."""
        # Test with redshift array and metallicity bins
        z = np.array([0.0, 2.0])
        met_bins = np.array([0.001, 0.01, 0.02, 0.03])
        
        result = chruslinska_model.fSFR(z, met_bins)
        # Shape check - should be (len(z), len(met_bins)-1)
        assert result.shape == (2, 3)
        
        # Test with normalization
        chruslinska_model.normalise = True
        result = chruslinska_model.fSFR(z, met_bins)
        for row in result:
            if np.sum(row) > 0:
                np.testing.assert_allclose(np.sum(row), 1.0)

class TestZavala21:
    """Tests for the Zavala21 SFH model with mocked data loading."""
    
    @pytest.fixture
    def mock_zavala_data(self, monkeypatch):
        """Create mock data for the Zavala21 class."""
        # Create mock data - simple decreasing function with redshift
        redshifts = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
        SFRD_min = 0.1 * np.exp(-redshifts / 3.0)  # Simple declining function
        SFRD_max = 0.2 * np.exp(-redshifts / 3.0)  # Double the min values
        
        # Create a mock _load_zavala_data method
        def mock_load_zavala(self):
            self.redshifts = redshifts
            if self.sub_model == "min":
                self.SFR_data = SFRD_min
            elif self.sub_model == "max":
                self.SFR_data = SFRD_max
            else:
                raise ValueError("Invalid sub-model!")
        
        # Patch the method
        monkeypatch.setattr(Zavala21, "_load_zavala_data", mock_load_zavala)
        
        return {
            "redshifts": redshifts,
            "SFRD_min": SFRD_min,
            "SFRD_max": SFRD_max
        }
    
    def test_init_parameters(self, mock_zavala_data):
        """Test that initialization validates and sets parameters correctly."""
        # Test missing sub_model
        with pytest.raises(ValueError) as excinfo:
            Zavala21({"Z_max": 0.03, "sigma": 0.5})
        assert "Sub-model not given!" in str(excinfo.value)
        
        # Test valid initialization with min model
        model_dict = {"sub_model": "min", "Z_max": 0.03, "sigma": 0.5}
        zavala_min = Zavala21(model_dict)
        assert zavala_min.sub_model == "min"
        assert zavala_min.Z_max == 0.03
        assert zavala_min.sigma == 0.5
        
        # Test valid initialization with max model
        model_dict = {"sub_model": "max", "Z_max": 0.03, "sigma": 0.5}
        zavala_max = Zavala21(model_dict)
        assert zavala_max.sub_model == "max"
        
        # Test invalid sub_model
        model_dict = {"sub_model": "invalid", "Z_max": 0.03, "sigma": 0.5}
        with pytest.raises(ValueError) as excinfo:
            zavala_invalid = Zavala21(model_dict)
        assert "Invalid sub-model!" in str(excinfo.value)
    
    def test_csfrd_min_model(self, mock_zavala_data):
        """Test the CSFRD method with min sub-model."""
        model_dict = {"sub_model": "min", "Z_max": 0.03, "sigma": 0.5}
        zavala = Zavala21(model_dict)
        
        # Test at specific redshifts
        z_values = np.array([0.0, 2.0, 4.0, 6.0])
        result = zavala.CSFRD(z_values)
        
        # Expected values come from interpolating the mock data
        expected = 0.1 * np.exp(-z_values / 3.0)
        np.testing.assert_allclose(result, expected)
    
    def test_csfrd_max_model(self, mock_zavala_data):
        """Test the CSFRD method with max sub-model."""
        model_dict = {"sub_model": "max", "Z_max": 0.03, "sigma": 0.5}
        zavala = Zavala21(model_dict)
        
        # Test at specific redshifts
        z_values = np.array([0.0, 2.0, 4.0, 6.0])
        result = zavala.CSFRD(z_values)
        
        # Expected values come from interpolating the mock data
        expected = 0.2 * np.exp(-z_values / 3.0)
        np.testing.assert_allclose(result, expected)
    
    def test_fsfr_calculation(self, mock_zavala_data):
        """Test the fSFR method which is inherited from MadauBase."""
        model_dict = {"sub_model": "min", "Z_max": 0.03, "sigma": 0.5}
        zavala = Zavala21(model_dict)
        
        # Test with redshift array and metallicity bins
        z = np.array([0.0, 2.0])
        met_bins = np.array([0.001, 0.01, 0.02, 0.03])
        
        result = zavala.fSFR(z, met_bins)
        
        # Shape check - should be (len(z), len(met_bins)-1)
        assert result.shape == (2, 3)
        
        # Test normalization
        zavala.normalise = True
        result = zavala.fSFR(z, met_bins)
        for row in result:
            np.testing.assert_allclose(np.sum(row), 1.0)

class TestGetSFHModel:
    """Tests for the get_SFH_model function."""
    
    def test_returns_correct_instance(self):
        """Test that get_SFH_model returns the correct instance for each model."""
        # Test for MadauDickinson14
        model_dict = {"SFR": "Madau+Dickinson14", "sigma": 0.5, "Z_max": 0.03}
        model = get_SFH_model(model_dict)
        assert isinstance(model, MadauDickinson14)
        
        # Test for MadauFragos17
        model_dict = {"SFR": "Madau+Fragos17", "sigma": 0.5, "Z_max": 0.03}
        model = get_SFH_model(model_dict)
        assert isinstance(model, MadauFragos17)
        
        # Test for Neijssel19
        model_dict = {"SFR": "Neijssel+19", "sigma": 0.5, "Z_max": 0.03}
        model = get_SFH_model(model_dict)
        assert isinstance(model, Neijssel19)
        
        # Test for Fujimoto24
        model_dict = {"SFR": "Fujimoto+24", "sigma": 0.5, "Z_max": 0.03}
        model = get_SFH_model(model_dict)
        assert isinstance(model, Fujimoto24)
    