import pytest
import numpy as np
from posydon.popsyn.normalized_pop_mass import initial_total_underlying_mass

@pytest.fixture
def common_kwargs():
    return {
        'primary_mass_scheme': 'Salpeter',
        'secondary_mass_scheme': 'flat_mass_ratio',
        'binary_fraction_const': 1,
        'primary_mass_min': 0.08,
        'primary_mass_max': 200.,
        'secondary_mass_min': 0.5,
        'secondary_mass_max': 200.,
    }

def test_initial_total_mass_with_fb1(common_kwargs):
    """Test initial_total_underlying_mass with specific input values."""
    kwargs = common_kwargs.copy()
    simulated_mass = 1000.0
    simulated_mass_single = 0
    simulated_mass_binary = 1000.0

    underlying_mass, f_corr_single, f_corr_binary = initial_total_underlying_mass(
        simulated_mass=simulated_mass,
        simulated_mass_single=simulated_mass_single,
        simulated_mass_binaries=simulated_mass_binary,
        f_bin=1,
        **kwargs
    )
    # Check that we get valid float values
    assert isinstance(underlying_mass, float)
    assert isinstance(f_corr_single, int)
    assert isinstance(f_corr_binary, float)

    assert underlying_mass == pytest.approx(2144.551055481399)
    assert f_corr_single == pytest.approx(0)
    assert f_corr_binary == pytest.approx(0.46629806151923225)

    # Test with f_bin = 0.7
    underlying_mass_07, f_corr_single_07, f_corr_binary_07 = initial_total_underlying_mass(
        simulated_mass=simulated_mass,
        simulated_mass_single=simulated_mass_single,
        simulated_mass_binaries=simulated_mass_binary,
        f_bin=0.7,
        **kwargs
    )

    assert underlying_mass_07 == pytest.approx(2757.2799284760854)
    assert f_corr_single_07 == pytest.approx(0)
    assert f_corr_binary_07 == pytest.approx(0.3626762700705139)

    # Test with f_bin = 0.5
    underlying_mass_05, f_corr_single_05, f_corr_binary_05 = initial_total_underlying_mass(
        simulated_mass=simulated_mass,
        simulated_mass_single=simulated_mass_single,
        simulated_mass_binaries=simulated_mass_binary,
        f_bin=0.5,
        **kwargs
    )

    assert underlying_mass_05 == pytest.approx(3574.2517591356655)
    assert f_corr_single_05 == pytest.approx(0)
    assert f_corr_binary_05 == pytest.approx(0.2797788369115394)

    # Test with f_bin = 0

    underlying_mass_0, f_corr_single_0, f_corr_binary_0 = initial_total_underlying_mass(
        simulated_mass=simulated_mass,
        simulated_mass_single=simulated_mass_single,
        simulated_mass_binaries=simulated_mass_binary,
        f_bin=0,
        **kwargs
    )
    print(underlying_mass_0, f_corr_single_0, f_corr_binary_0)
    assert np.isnan(underlying_mass_0)
    assert np.isnan(f_corr_single_0)
    assert np.isnan(f_corr_binary_0)


def test_initial_total_mass_with_fb07(common_kwargs):
    """Test initial_total_underlying_mass with binary_fraction_const = 0.7."""
    kwargs = common_kwargs.copy()
    kwargs['binary_fraction_const'] = 0.7  # Change binary fraction in original population
    
    simulated_mass = 1000.0
    simulated_mass_single = 300.0
    simulated_mass_binary = 700.0

    # Test with f_bin = 1
    underlying_mass_1, f_corr_single_1, f_corr_binary_1 = initial_total_underlying_mass(
        simulated_mass=simulated_mass,
        simulated_mass_single=simulated_mass_single,
        simulated_mass_binaries=simulated_mass_binary,
        f_bin=1,
        **kwargs
    )

    assert underlying_mass_1 == pytest.approx(1167.5889079843173)
    assert f_corr_single_1 == pytest.approx(0)
    assert f_corr_binary_1 == pytest.approx(0.5995260790961559)

    # Test with f_bin = 0.7 (matching binary_fraction_const)
    underlying_mass_07, f_corr_single_07, f_corr_binary_07 = initial_total_underlying_mass(
        simulated_mass=simulated_mass,
        simulated_mass_single=simulated_mass_single,
        simulated_mass_binaries=simulated_mass_binary,
        f_bin=0.7,
        **kwargs
    )

    assert underlying_mass_07 == pytest.approx(2144.5510554813995)
    assert f_corr_single_07 == pytest.approx(0.4662980615192322)
    assert f_corr_binary_07 == pytest.approx(0.4662980615192322)

    # Test with f_bin = 0.5
    underlying_mass_05, f_corr_single_05, f_corr_binary_05 = initial_total_underlying_mass(
        simulated_mass=simulated_mass,
        simulated_mass_single=simulated_mass_single,
        simulated_mass_binaries=simulated_mass_binary,
        f_bin=0.5,
        **kwargs
    )

    assert underlying_mass_05 == pytest.approx(2303.4066892207616)
    assert f_corr_single_05 == pytest.approx(0.8393365107346181)
    assert f_corr_binary_05 == pytest.approx(0.3597156474576936)
    # Test with f_bin = 0
    with pytest.raises(ZeroDivisionError) as excinfo:
        initial_total_underlying_mass(
            simulated_mass=simulated_mass,
            simulated_mass_single=simulated_mass_single,
            simulated_mass_binaries=simulated_mass_binary,
            f_bin=0,
            **kwargs
        )
    assert "float division by zero" in str(excinfo.value)
