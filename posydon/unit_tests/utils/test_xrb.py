"""Unit tests for :mod:`posydon.utils.xrb`."""

import numpy as np
import pytest

from posydon.utils import constants as const
from posydon.utils import xrb


def test_black_hole_radiative_efficiency_scalars_and_arrays():
    assert isinstance(xrb.black_hole_radiative_efficiency(0.0), float)
    assert xrb.black_hole_radiative_efficiency(0.0) == pytest.approx(
        0.057190958417936644
    )
    efficiencies = xrb.black_hole_radiative_efficiency(
        np.array([-0.998, 0.0, 0.9, 1.0, 1.1, np.nan])
    )
    assert efficiencies[:3] == pytest.approx(
        [0.03777362529209353, 0.057190958417936644,
         0.15575299199446402]
    )
    assert efficiencies[3] == pytest.approx(1.0 - 1.0 / np.sqrt(3.0))
    assert np.isnan(efficiencies[4:]).all()


def test_neutron_star_radiative_efficiency_and_broadcasting():
    radius = 1.25e6 / const.Rsun
    expected = (const.standard_cgrav * 1.4 * const.Msun
                / (1.25e6 * const.clight**2))
    assert xrb.neutron_star_radiative_efficiency(1.4, radius) == pytest.approx(
        expected
    )
    values = xrb.neutron_star_radiative_efficiency([1.4, -1.0], radius)
    assert values[0] == pytest.approx(expected)
    assert np.isnan(values[1])
    with pytest.raises(ValueError):
        xrb.neutron_star_radiative_efficiency([1.4, 1.5], [radius] * 3)


def test_eddington_limit_units_and_invalid_elements():
    mdot_edd, luminosity_edd = xrb.eddington_limit_from_properties(
        10.0, 0.7, 0.1
    )
    expected_luminosity = (
        4.0 * np.pi * const.standard_cgrav * 10.0 * const.Msun
        * const.clight / (0.2 * 1.7)
    )
    expected_rate = (expected_luminosity / (0.1 * const.clight**2)
                     * const.secyer / const.Msun)
    assert isinstance(mdot_edd, float)
    assert isinstance(luminosity_edd, float)
    assert mdot_edd == pytest.approx(expected_rate)
    assert luminosity_edd == pytest.approx(expected_luminosity)

    rates, luminosities = xrb.eddington_limit_from_properties(
        [10.0, -1.0, 10.0], 0.7, [0.1, 0.1, 0.0]
    )
    assert np.isfinite(rates[0]) and np.isfinite(luminosities[0])
    assert np.isnan(rates[1:]).all()
    assert np.isnan(luminosities[1:]).all()


def test_wind_velocity_schemes_and_validation():
    inputs = dict(
        donor_mass=20.0,
        donor_radius=10.0,
        donor_luminosity=1.0e5,
        donor_wind_mass_loss=1.0e-6,
        donor_surface_h1=0.7,
        donor_he_core_mass=5.0,
    )
    escape_speed = np.sqrt(
        2.0 * const.standard_cgrav * 20.0 * const.Msun
        / (10.0 * const.Rsun)
    )
    assert xrb.wind_velocity(**inputs) == pytest.approx(2.65 * escape_speed)
    hurley_factor = np.sqrt(0.5 + (20.0 - 1.4) * 6.5 / 118.6)
    assert xrb.wind_velocity(**inputs, scheme="Hurley+2002") == pytest.approx(
        hurley_factor * escape_speed
    )
    with pytest.raises(ValueError, match="scheme"):
        xrb.wind_velocity(**inputs, scheme="unknown")
    invalid = inputs | {"donor_radius": [10.0, -1.0]}
    values = xrb.wind_velocity(**invalid)
    assert np.isfinite(values[0]) and np.isnan(values[1])


def test_bondi_hoyle_circular_value_and_broadcasting():
    inputs = dict(
        accretor_mass=10.0,
        donor_mass=20.0,
        donor_wind_mass_loss=1.0e-6,
        orbital_separation=100.0,
        eccentricity=0.0,
        wind_speed=1.0e8,
    )
    separation = 100.0 * const.Rsun
    orbital_speed_squared = (
        const.standard_cgrav * 30.0 * const.Msun / separation
    )
    expected_fraction = (
        1.5 * (const.standard_cgrav * 10.0 * const.Msun)**2
        / (2.0 * separation**2 * 1.0e8**4)
        * (1.0 + orbital_speed_squared / 1.0e16)**-1.5
    )
    rate = xrb.bondi_hoyle_accretion_rate(**inputs)
    assert isinstance(rate, float)
    assert rate == pytest.approx(expected_fraction * 1.0e-6)

    values = xrb.bondi_hoyle_accretion_rate(
        [10.0, -1.0], 20.0, 1.0e-6, 100.0, 0.0, 1.0e8
    )
    assert values[0] == pytest.approx(rate)
    assert np.isnan(values[1])
    with pytest.raises(ValueError, match="alpha"):
        xrb.bondi_hoyle_accretion_rate(**inputs, alpha=0.0)


def test_bondi_hoyle_eccentric_orbit_average():
    eccentricity = 0.6
    circular = xrb.bondi_hoyle_accretion_rate(
        10.0, 20.0, 1.0e-6, 100.0, 0.0, 1.0e8
    )
    eccentric = xrb.bondi_hoyle_accretion_rate(
        10.0, 20.0, 1.0e-6, 100.0, eccentricity, 1.0e8
    )

    # Numerically average 1/r^2 over mean anomaly. With the constant-wind and
    # characteristic-orbital-speed assumptions used by the analytic model,
    # this is the only phase-dependent factor.
    eccentric_anomaly = ((np.arange(200000) + 0.5)
                         * 2.0 * np.pi / 200000)
    radius_over_a = 1.0 - eccentricity * np.cos(eccentric_anomaly)
    dmean_anomaly_de = radius_over_a
    average_inverse_radius_squared = np.mean(
        dmean_anomaly_de / radius_over_a**2
    )
    assert average_inverse_radius_squared == pytest.approx(
        1.0 / np.sqrt(1.0 - eccentricity**2), rel=2.0e-10
    )
    assert eccentric / circular == pytest.approx(
        average_inverse_radius_squared, rel=2.0e-10
    )
    assert eccentric == xrb.bondi_hoyle_accretion_rate(
        10.0, 20.0, 1.0e-6, 100.0, eccentricity, 1.0e8
    )


def test_accretion_luminosity_regimes_and_boundaries():
    mdot_edd, luminosity_edd = xrb.eddington_limit_from_properties(
        10.0, 0.7, 0.1
    )
    rates = mdot_edd * np.array([0.0, 0.5, 1.0, 1.0 + 1e-10,
                                 8.5, 8.5 + 1e-10, 1.0e4])
    luminosity, beaming, ratio = xrb.accretion_luminosity(
        rates, 10.0, 0.7, 0.1
    )
    assert ratio == pytest.approx(
        [0.0, 0.5, 1.0, 1.0 + 1e-10, 8.5, 8.5 + 1e-10, 1.0e4]
    )
    assert luminosity[0] == 0.0
    assert luminosity[1] == pytest.approx(0.5 * luminosity_edd)
    assert luminosity[2] == pytest.approx(luminosity_edd)
    assert luminosity[3] == pytest.approx(
        luminosity_edd * (1.0 + np.log(ratio[3]))
    )
    assert np.all(beaming[:6] == 1.0)
    assert luminosity[5] == pytest.approx(
        luminosity_edd * (1.0 + np.log(ratio[5]))
    )
    assert beaming[-1] == pytest.approx(3.2e-3)

    scalar = xrb.accretion_luminosity(mdot_edd, 10.0, 0.7, 0.1)
    assert all(isinstance(value, float) for value in scalar)
    invalid = xrb.accretion_luminosity(-1.0, 10.0, 0.7, 0.1)
    assert all(np.isnan(value) for value in invalid)
    with pytest.raises(ValueError, match="beaming_floor"):
        xrb.accretion_luminosity(mdot_edd, 10.0, 0.7, 0.1,
                                 beaming_floor=0.0)


def test_be_xray_luminosity_valid_and_invalid_periods():
    expected = 1.0e35 * 10.0**(4.53 - 1.5 * np.log10(10.0))
    value = xrb.be_xray_luminosity(10.0)
    assert isinstance(value, float)
    assert value == pytest.approx(expected)
    values = xrb.be_xray_luminosity([10.0, 100.0, 0.0, -1.0, np.nan])
    assert values[:2] == pytest.approx(
        [expected, 1.0e35 * 10.0**(4.53 - 3.0)]
    )
    assert np.isnan(values[2:]).all()
