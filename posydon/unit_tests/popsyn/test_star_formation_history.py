import numpy as np
import pytest
from posydon.popsyn.star_formation_history import (
    MadauDickinson14,
    MadauFragos17,
    Neijssel19,
    IllustrisTNG,
)

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
    
    def dummy_get_illustrisTNG_data(self, verbose=False):
        import numpy as np
        return {
            "SFR": np.array([1.0, 1.0, 1.0]),
            "redshifts": np.array([0.0, 1.0, 2.0]),
            "mets": np.linspace(0.001, 0.03, 10),
            "M": np.ones((3, 10)),
        }
    
    def test_fSFR_sum_is_one_illustris(self, monkeypatch):
        base_args = {
            "SFR": "IllustrisTNG",
            "sigma": 0.5,
            "Z_max": 0.03,
            "select_one_met": False,
        }
        monkeypatch.setattr(IllustrisTNG, "_get_illustrisTNG_data", self.dummy_get_illustrisTNG_data)
        sfh_instance = IllustrisTNG(base_args)
        z = np.array([0.5, 1.0, 2.0])
        met_bins = np.linspace(0.0001, base_args["Z_max"], 50)
        fSFR = sfh_instance.fSFR(z, met_bins)
        np.testing.assert_allclose(np.sum(fSFR, axis=1), np.ones(len(z)), atol=1e-6)




