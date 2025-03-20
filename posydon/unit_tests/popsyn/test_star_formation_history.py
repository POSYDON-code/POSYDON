import numpy as np
import pytest
from posydon.popsyn.star_formation_history import SFHBase, MadauBase
from posydon.popsyn.star_formation_history import (
    MadauDickinson14,
    MadauFragos17,
    Neijssel19,
    IllustrisTNG,
    Chruslinska21,
    Zalava21,
    get_SFH_model
)


# Replace duplicate DummySFH definitions with a single merged class
class DummySFH(MadauBase):
    def CSFRD(self, z):
        if self.MODEL.get("dummy_mode", "std") == "call":
            return np.full_like(z, 2.0, dtype=float)
        else:
            return np.ones_like(z)

    def mean_metallicity(self, z):
        return np.ones_like(z)

    def fSFR(self, z, metallicity_bins):
        n_z = len(z)
        n_bins = len(metallicity_bins) - 1
        if self.MODEL.get("dummy_mode", "std") == "call":
            # Return normalized ones for __call__ test
            raw = np.ones((n_z, n_bins))
            return raw / np.sum(raw, axis=1, keepdims=True)
        else:
            return np.ones((n_z, n_bins))

# Additional tests for the std_log_metallicity_dist function
class TestStdLogMetallicityDist:
    def test_sigma_bavera(self):
        # Test with sigma as "Bavera+20" which should return 0.5
        model = {"sigma": "Bavera+20", "Z_max": 0.03, "select_one_met": False}
        dummy = DummySFH(model)
        assert dummy.std_log_metallicity_dist() == 0.5

    def test_sigma_neijssel(self):
        # Test with sigma as "Neijssel+19" which should return 0.39
        model = {"sigma": "Neijssel+19", "Z_max": 0.03, "select_one_met": False}
        dummy = DummySFH(model)
        assert dummy.std_log_metallicity_dist() == 0.39

    def test_sigma_float(self):
        # Test with sigma as a float value
        sigma_value = 0.45
        model = {"sigma": sigma_value, "Z_max": 0.03, "select_one_met": False}
        dummy = DummySFH(model)
        assert dummy.std_log_metallicity_dist() == sigma_value

    def test_unknown_sigma_string(self):
        # Test with an invalid sigma string should raise a ValueError
        model = {"sigma": "invalid_sigma",
                 "Z_max": 0.03,
                 "select_one_met": False}
        dummy = DummySFH(model)
        with pytest.raises(ValueError) as excinfo:
            dummy.std_log_metallicity_dist()
        assert "Unknown sigma choice!" in str(excinfo.value)
        
    def test_invalid_sigma(self):
        # Test with an invalid sigma value should raise a ValueError
        model = {"sigma": int(1), "Z_max": 0.03, "select_one_met": False}
        dummy = DummySFH(model)
        with pytest.raises(ValueError) as excinfo:
            dummy.std_log_metallicity_dist()
        assert "Invalid sigma value" in str(excinfo.value)
            


class TestfSFR:
    @pytest.mark.parametrize("model_class, model_name", [
        (MadauFragos17, "Madau+Fragos17"),
        (MadauDickinson14, "Madau+Dickinson14"),
        (Neijssel19, "Neijssel+19"),
    ])
    def test_fSFR_sum_is_one_empirical(self, model_class, model_name):
        # Base MODEL parameters
        base_args = {
            "SFR": model_name,
            "sigma": 0.5,
            "Z_max": 0.03,
            "select_one_met": False,
        }
        sfh_instance = model_class(base_args)
        # Arbitrary redshift array and metallicity bins
        z = np.array([0.001, 0.5, 1.0, 2.0])
        met_bins = np.linspace(0.0001, base_args["Z_max"], 50)
        
        fSFR = sfh_instance.fSFR(z, met_bins)
        # The sum over metallicity bins (axis=1) should be approximately 1 for each redshift
        np.testing.assert_allclose(np.sum(fSFR, axis=1), np.ones(len(z)), atol=1e-6)
    
    def test_fSFR_sum_is_one_illustris(self, monkeypatch):
        def dummy_get_illustrisTNG_data(self, verbose=False):
            return {
                "SFR": np.array([1.0, 1.0, 1.0]),
                "redshifts": np.array([0.0, 1.0, 2.0]),
                "mets": np.linspace(0.001, 0.03, 10),
                "M": np.ones((3, 10)),
            }
        base_args = {
            "SFR": "IllustrisTNG",
            "sigma": 0.5,
            "Z_max": 0.03,
            "select_one_met": False,
        }
        monkeypatch.setattr(IllustrisTNG, "_get_illustrisTNG_data", dummy_get_illustrisTNG_data)
        sfh_instance = IllustrisTNG(base_args)
        z = np.array([0.5, 1.0, 2.0])
        met_bins = np.linspace(0.0001, base_args["Z_max"], 50)
        fSFR = sfh_instance.fSFR(z, met_bins)
        np.testing.assert_allclose(np.sum(fSFR, axis=1), np.ones(len(z)), atol=1e-6)

class TestCallMethod:
    def test_call_method_returns_product(self):
        # Use DummySFH in 'call' mode to test __call__
        model = {"sigma": 0.5, "Z_max": 0.03, "select_one_met": False, "dummy_mode": "call"}
        dummy = DummySFH(model)
        z = np.array([0.5, 1.0, 2.0])
        met_bins = np.linspace(0.001, model["Z_max"], 10)
        expected = 2 * dummy.fSFR(z, met_bins)
        result = dummy(z, met_bins)
        np.testing.assert_allclose(result, expected, atol=1e-6)


class TestGetSFHModel:
    @pytest.mark.parametrize("model_name, model_class", [
        ("Madau+Fragos17", MadauFragos17),
        ("Madau+Dickinson14", MadauDickinson14),
        ("Neijssel+19", Neijssel19),
    ])
    def test_get_model_empirical(self, model_name, model_class):
        base_args = {
            "SFR": model_name,
            "sigma": 0.5,
            "Z_max": 0.03,
            "select_one_met": False,
        }
        model = get_SFH_model(base_args)
        assert isinstance(model, model_class)

    def test_get_model_illustris(self, monkeypatch):
        base_args = {
            "SFR": "IllustrisTNG",
            "sigma": 0.5,
            "Z_max": 0.03,
            "select_one_met": False,
        }
        # Override _get_illustrisTNG_data to avoid file I/O during tests
        def dummy_get_illustrisTNG_data(self, verbose=False):
            return {
                "SFR": np.array([1.0]),
                "redshifts": np.array([0.0]),
                "mets": np.linspace(0.001, 0.03, 10),
                "M": np.ones((1, 10)),
            }
        monkeypatch.setattr(IllustrisTNG, "_get_illustrisTNG_data", dummy_get_illustrisTNG_data)
        model = get_SFH_model(base_args)
        assert isinstance(model, IllustrisTNG)

    def test_get_model_chruslinska(self, monkeypatch):
        base_args = {
            "SFR": "Chruslinska+21",
            "sigma": 0.5,
            "Z_max": 0.03,
            "select_one_met": False,
            "sub_model": "dummy",
            "Z_solar_scaling": "Asplund09",
        }
        # Override _load_chruslinska_data to avoid file I/O during tests
        monkeypatch.setattr(Chruslinska21, "_load_chruslinska_data", lambda self, verbose=False: None)
        model = get_SFH_model(base_args)
        assert isinstance(model, Chruslinska21)
        assert model.SFR == "Chruslinska+21"
        
    def test_get_model_zalava(self, monkeypatch):
        base_args = {
            "SFR": "Zalava+21",
            "sigma": 0.5,
            "Z_max": 0.03,
            "select_one_met": False,
            "sub_model": "dummy",
            "Z_solar_scaling": "Asplund09",
        }
        # Override _load_zalava_data to avoid file I/O during tests
        monkeypatch.setattr(Zalava21, "_load_zalava_data", lambda self, verbose=False: None)
        model = get_SFH_model(base_args)
        assert isinstance(model, Zalava21)
        assert model.SFR == "Zalava+21"
        
    def test_get_fojimoto_model(self):
        base_args = {
            "SFR": "Fujimoto+24",
            "sigma": 0.5,
            "Z_max": 0.03,
            "select_one_met": False,
        }
        model = get_SFH_model(base_args)
        assert isinstance(model, SFHBase)
        assert model.MODEL["SFR"] == "Fujimoto+24"
        assert model.SFR == "Fujimoto+24"

    def test_invalid_SFR(self):
        base_args = {
            "SFR": "InvalidSFR",
            "sigma": 0.5,
            "Z_max": 0.03,
            "select_one_met": False,
        }
        with pytest.raises(ValueError) as excinfo:
            get_SFH_model(base_args)
        assert "Invalid SFR!" in str(excinfo.value)
