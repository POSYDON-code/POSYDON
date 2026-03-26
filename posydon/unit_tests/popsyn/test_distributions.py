"""Unit tests for posydon/popsyn/distributions.py

This module contains comprehensive unit tests for the distribution classes
used in population synthesis: FlatMassRatio, Sana12Period, and PowerLawPeriod.
"""

import numpy as np
import pytest
from scipy.integrate import quad

from posydon.popsyn.distributions import (
    FlatMassRatio,
    LogNormalSeparation,
    LogUniform,
    PowerLawPeriod,
    Sana12Period,
    ThermalEccentricity,
    UniformEccentricity,
    ZeroEccentricity,
)


class TestFlatMassRatio:
    """Test class for FlatMassRatio distribution."""

    @pytest.fixture
    def default_flat_ratio(self):
        """Fixture for default FlatMassRatio instance."""
        return FlatMassRatio()

    @pytest.fixture
    def custom_flat_ratio(self):
        """Fixture for custom FlatMassRatio instance."""
        return FlatMassRatio(q_min=0.1, q_max=0.8)

    def test_initialization_default(self, default_flat_ratio):
        """Test default initialization of FlatMassRatio."""
        assert default_flat_ratio.q_min == 0.05
        assert default_flat_ratio.q_max == 1.0
        assert hasattr(default_flat_ratio, 'norm')
        assert default_flat_ratio.norm > 0

    def test_initialization_custom(self, custom_flat_ratio):
        """Test custom initialization of FlatMassRatio."""
        assert custom_flat_ratio.q_min == 0.1
        assert custom_flat_ratio.q_max == 0.8
        assert hasattr(custom_flat_ratio, 'norm')
        assert custom_flat_ratio.norm > 0

    def test_initialization_invalid_parameters(self):
        """Test that initialization raises ValueError for invalid parameters."""
        with pytest.raises(ValueError, match="q_min must be in \\[0, 1\\)"):
            FlatMassRatio(q_min=-0.1, q_max=0.5)

        with pytest.raises(ValueError, match="q_min must be in \\[0, 1\\)"):
            FlatMassRatio(q_min=1.5, q_max=2.0)

        # Test q_max not in (0, 1]
        with pytest.raises(ValueError, match="q_max must be in \\(0, 1\\]"):
            FlatMassRatio(q_min=0.1, q_max=0.0)

        with pytest.raises(ValueError, match="q_max must be in \\(0, 1\\]"):
            FlatMassRatio(q_min=0.1, q_max=-0.1)

        with pytest.raises(ValueError, match="q_max must be in \\(0, 1\\]"):
            FlatMassRatio(q_min=0.1, q_max=1.5)

        # Test q_min >= q_max
        with pytest.raises(ValueError, match="q_min must be less than q_max"):
            FlatMassRatio(q_min=0.8, q_max=0.5)

        with pytest.raises(ValueError, match="q_min must be less than q_max"):
            FlatMassRatio(q_min=0.5, q_max=0.5)

    def test_repr(self, default_flat_ratio, custom_flat_ratio):
        """Test string representation of the distribution."""
        default_repr = default_flat_ratio.__repr__()
        assert "FlatMassRatio(" in default_repr
        assert "q_min=0.05" in default_repr
        assert "q_max=1" in default_repr

        custom_repr = custom_flat_ratio.__repr__()
        assert "FlatMassRatio(" in custom_repr
        assert "q_min=0.1" in custom_repr
        assert "q_max=0.8" in custom_repr

    def test_repr_html(self, default_flat_ratio):
        """Test HTML representation for Jupyter notebooks."""
        html_str = default_flat_ratio._repr_html_()
        assert "<h3>Flat Mass Ratio Distribution</h3>" in html_str
        assert "q_min = 0.05" in html_str
        assert "q_max = 1" in html_str

    def test_flat_mass_ratio_method(self, default_flat_ratio):
        """Test the flat_mass_ratio method returns constant value."""
        # Test with single value
        result = default_flat_ratio.flat_mass_ratio(0.5)
        assert result == 1.0

        # Test with array
        q_values = np.array([0.1, 0.5, 0.9])
        result = default_flat_ratio.flat_mass_ratio(q_values)
        expected = np.ones_like(q_values)
        np.testing.assert_array_equal(result, expected)

    def test_calculate_normalization(self, custom_flat_ratio):
        """Test that normalization calculation works correctly."""
        # For a flat distribution, normalization should be 1/(q_max - q_min)
        expected_norm = 1.0 / (custom_flat_ratio.q_max - custom_flat_ratio.q_min)
        np.testing.assert_allclose(custom_flat_ratio.norm, expected_norm)

    def test_pdf_within_range(self, custom_flat_ratio):
        """Test PDF returns correct values within the mass ratio range."""
        q_values = np.linspace(custom_flat_ratio.q_min, custom_flat_ratio.q_max, 10)
        pdf_values = custom_flat_ratio.pdf(q_values)
        expected_pdf = custom_flat_ratio.norm * np.ones_like(q_values)
        expected_pdf[0] = 0.0  # q_values[0] is equal to q_min, which is outside the valid range
        np.testing.assert_allclose(pdf_values, expected_pdf)

    def test_pdf_outside_range(self, custom_flat_ratio):
        """Test PDF returns zero for mass ratios outside the range."""
        # Below range
        q_below = np.array([0.05, 0.09])
        pdf_below = custom_flat_ratio.pdf(q_below)
        np.testing.assert_array_equal(pdf_below, np.zeros_like(q_below))

        # Above range
        q_above = np.array([0.85, 0.95])
        pdf_above = custom_flat_ratio.pdf(q_above)
        np.testing.assert_array_equal(pdf_above, np.zeros_like(q_above))

    def test_pdf_mixed_range(self, custom_flat_ratio):
        """Test PDF with mixture of values inside and outside range."""
        q_mixed = np.array([0.05, 0.2, 0.5, 0.9, 1.0])
        pdf_mixed = custom_flat_ratio.pdf(q_mixed)
        expected = np.array([0.0, custom_flat_ratio.norm, custom_flat_ratio.norm, 0.0, 0.0])
        np.testing.assert_allclose(pdf_mixed, expected)

    def test_pdf_scalar_input(self, default_flat_ratio):
        """Test PDF with scalar input."""
        # Within range
        q = 0.5
        pdf_value = default_flat_ratio.pdf(q)
        expected = default_flat_ratio.norm
        assert np.isclose(pdf_value, expected)

        # Outside range
        q_outside = 1.5
        pdf_outside = default_flat_ratio.pdf(q_outside)
        assert np.isclose(pdf_outside, 0.0)

    def test_normalization_integral(self, default_flat_ratio):
        """Test that the PDF integrates to 1 over the valid range."""
        integral, _ = quad(default_flat_ratio.pdf, default_flat_ratio.q_min, default_flat_ratio.q_max)
        np.testing.assert_allclose(integral, 1.0, rtol=1e-10)

    def test_rvs(self, custom_flat_ratio):
        """Test random sampling."""
        rng = np.random.default_rng(42)
        samples = custom_flat_ratio.rvs(size=1000, rng=rng)

        assert len(samples) == 1000
        assert np.all(samples >= custom_flat_ratio.q_min)
        assert np.all(samples <= custom_flat_ratio.q_max)

    def test_rvs_without_rng(self, custom_flat_ratio):
        """Test random sampling without providing an RNG."""
        samples = custom_flat_ratio.rvs(size=100)

        assert len(samples) == 100
        assert np.all(samples >= custom_flat_ratio.q_min)
        assert np.all(samples <= custom_flat_ratio.q_max)


class TestSana12Period:
    """Test class for Sana12Period distribution."""

    @pytest.fixture
    def default_sana12(self):
        """Fixture for default Sana12Period instance."""
        return Sana12Period()

    @pytest.fixture
    def custom_sana12(self):
        """Fixture for custom Sana12Period instance."""
        return Sana12Period(p_min=1.0, p_max=1000.0)

    def test_initialization_default(self, default_sana12):
        """Test default initialization of Sana12Period."""
        assert default_sana12.p_min == 0.35
        assert default_sana12.p_max == 6e3
        assert default_sana12.mbreak == 15
        assert default_sana12.slope == 0.55
        assert hasattr(default_sana12, 'low_mass_norm')
        assert hasattr(default_sana12, 'high_mass_norm')
        assert hasattr(default_sana12, 'norm')

    def test_initialization_custom(self, custom_sana12):
        """Test custom initialization of Sana12Period."""
        assert custom_sana12.p_min == 1.0
        assert custom_sana12.p_max == 1000.0
        assert custom_sana12.mbreak == 15
        assert custom_sana12.slope == 0.55

    def test_initialization_invalid_parameters(self):
        """Test that initialization raises ValueError for invalid parameters."""
        # Test p_min <= 0
        with pytest.raises(ValueError, match="p_min must be positive"):
            Sana12Period(p_min=0.0, p_max=1000.0)

        with pytest.raises(ValueError, match="p_min must be positive"):
            Sana12Period(p_min=-1.0, p_max=1000.0)

        # Test p_max <= p_min
        with pytest.raises(ValueError, match="p_max must be greater than p_min"):
            Sana12Period(p_min=1000.0, p_max=100.0)

        with pytest.raises(ValueError, match="p_max must be greater than p_min"):
            Sana12Period(p_min=100.0, p_max=100.0)

    def test_repr(self, default_sana12, custom_sana12):
        """Test string representation of the distribution."""
        default_repr = default_sana12.__repr__()
        assert "Sana12Period(" in default_repr
        assert "p_min=0.35" in default_repr
        assert "p_max=6000.0" in default_repr

        custom_repr = custom_sana12.__repr__()
        assert "Sana12Period(" in custom_repr
        assert "p_min=1.0" in custom_repr
        assert "p_max=1000.0" in custom_repr

    def test_repr_html(self, default_sana12):
        """Test HTML representation for Jupyter notebooks."""
        html_str = default_sana12._repr_html_()
        assert "<h3>Sana12 Period Distribution</h3>" in html_str
        assert "p_min = 0.35" in html_str
        assert "p_max = 6000.0" in html_str

    def test_sana12_period_low_mass(self, default_sana12):
        """Test sana12_period method for low mass stars (m1 <= 15)."""
        logp_values = np.array([0.0, 1.0, 2.0])
        m1_low = 10.0  # Below mbreak

        result = default_sana12.sana12_period(logp_values, m1_low)
        expected = np.ones_like(logp_values)  # Flat distribution for low mass
        np.testing.assert_allclose(result, expected)

    def test_sana12_period_high_mass(self, default_sana12):
        """Test sana12_period method for high mass stars (m1 > 15)."""
        logp_values = np.array([0.0, 1.0, 2.0])
        m1_high = 20.0  # Above mbreak

        result = default_sana12.sana12_period(logp_values, m1_high)
        # Should return power law values for high mass
        assert len(result) == len(logp_values)
        assert np.all(result >= 0)

    def test_sana12_period_invalid_mass(self, default_sana12):
        """Test that invalid mass values raise ValueError."""
        logp_values = np.array([0.0, 1.0])

        # Test zero mass
        with pytest.raises(ValueError, match="Mass must be positive"):
            default_sana12.sana12_period(logp_values, 0.0)

        # Test negative mass
        with pytest.raises(ValueError, match="Mass must be positive"):
            default_sana12.sana12_period(logp_values, -1.0)

        # Test array with negative mass
        m1_invalid = np.array([10.0, -5.0])
        with pytest.raises(ValueError, match="Mass must be positive"):
            default_sana12.sana12_period(logp_values, m1_invalid)

    def test_sana12_period_array_shapes(self, default_sana12):
        """Test that array shape mismatches raise ValueError."""
        logp_values = np.array([0.0, 1.0, 2.0])
        m1_wrong_shape = np.array([10.0, 15.0])  # Different length

        with pytest.raises(ValueError, match="m1 must be a single value or an array of the same length as p"):
            default_sana12.sana12_period(logp_values, m1_wrong_shape)

    def test_pdf_within_range(self, custom_sana12):
        """Test PDF returns correct values within the period range."""
        p_values = np.linspace(custom_sana12.p_min, custom_sana12.p_max, 10)
        m1_test = 10.0

        pdf_values = custom_sana12.pdf(p_values, m1_test)
        # Should return non-zero values within range
        assert np.all(pdf_values >= 0)
        assert np.any(pdf_values > 0)  # At least some should be positive

    def test_pdf_outside_range(self, custom_sana12):
        """Test PDF returns zero for periods outside the range."""
        # Below range
        p_below = np.array([0.1, 0.5])
        m1_test = 10.0
        pdf_below = custom_sana12.pdf(p_below, m1_test)
        np.testing.assert_array_equal(pdf_below, np.zeros_like(p_below))

        # Above range
        p_above = np.array([2000.0, 5000.0])
        pdf_above = custom_sana12.pdf(p_above, m1_test)
        np.testing.assert_array_equal(pdf_above, np.zeros_like(p_above))

    def test_pdf_array_broadcasting(self, default_sana12):
        """Test PDF with different array broadcasting scenarios."""
        p_values = np.array([1.0, 10.0, 100.0])

        # Single mass value
        m1_single = 10.0
        pdf_single = default_sana12.pdf(p_values, m1_single)
        assert len(pdf_single) == len(p_values)

        # Array of same length
        m1_array = np.array([10.0, 15.0, 20.0])
        pdf_array = default_sana12.pdf(p_values, m1_array)
        assert len(pdf_array) == len(p_values)

        # Wrong array length should raise error
        m1_wrong = np.array([10.0, 15.0])
        with pytest.raises(ValueError, match="m1 must be a single value or an array of the same length as p"):
            default_sana12.pdf(p_values, m1_wrong)

    def test_calculate_normalization_integration(self, default_sana12):
        """Test that normalization ensures PDF integrates to reasonable values."""
        # Test for low mass case
        m1_low = 10.0
        norm_low = default_sana12._calculate_normalization(m1_low)
        assert norm_low > 0

        # Test for high mass case
        m1_high = 20.0
        norm_high = default_sana12._calculate_normalization(m1_high)
        assert norm_high > 0

    def test_rvs_with_m1_none(self, default_sana12):
        """Test that rvs raises ValueError when m1 is None."""
        rng = np.random.default_rng(42)

        with pytest.raises(ValueError, match="m1 \\(primary mass\\) must be provided"):
            default_sana12.rvs(size=10, m1=None, rng=rng)

    def test_rvs_with_m1_wrong_size(self, default_sana12):
        """Test that rvs raises ValueError when m1 has wrong size."""
        rng = np.random.default_rng(42)
        m1_wrong_size = np.array([10.0, 15.0])

        with pytest.raises(ValueError, match="m1 must be a single value or have size="):
            default_sana12.rvs(size=10, m1=m1_wrong_size, rng=rng)

    def test_rvs_low_mass(self, default_sana12):
        """Test random sampling for low mass stars."""
        rng = np.random.default_rng(42)
        m1 = 10.0  # Below mbreak

        samples = default_sana12.rvs(size=100, m1=m1, rng=rng)

        assert len(samples) == 100
        assert np.all(samples >= default_sana12.p_min)
        assert np.all(samples <= default_sana12.p_max)

    def test_rvs_high_mass(self, default_sana12):
        """Test random sampling for high mass stars."""
        rng = np.random.default_rng(42)
        m1 = 25.0  # Above mbreak

        samples = default_sana12.rvs(size=100, m1=m1, rng=rng)

        assert len(samples) == 100
        assert np.all(samples >= default_sana12.p_min)
        assert np.all(samples <= default_sana12.p_max)

    def test_rvs_mixed_masses(self, default_sana12):
        """Test random sampling with array of masses."""
        rng = np.random.default_rng(42)
        m1 = np.array([10.0, 15.0, 20.0, 25.0, 30.0])

        samples = default_sana12.rvs(size=5, m1=m1, rng=rng)

        assert len(samples) == 5
        assert np.all(samples >= default_sana12.p_min)
        assert np.all(samples <= default_sana12.p_max)

    def test_rvs_without_rng(self, default_sana12):
        """Test random sampling without providing an RNG."""
        m1 = 20.0

        samples = default_sana12.rvs(size=100, m1=m1)

        assert len(samples) == 100
        assert np.all(samples >= default_sana12.p_min)
        assert np.all(samples <= default_sana12.p_max)


class TestPowerLawPeriod:
    """Test class for PowerLawPeriod distribution."""

    @pytest.fixture
    def default_power_law(self):
        """Fixture for default PowerLawPeriod instance."""
        return PowerLawPeriod()

    @pytest.fixture
    def custom_power_law(self):
        """Fixture for custom PowerLawPeriod instance."""
        return PowerLawPeriod(slope=-1.0, p_min=10.0, p_max=1e5)

    def test_initialization_default(self, default_power_law):
        """Test default initialization of PowerLawPeriod."""
        assert default_power_law.slope == -0.55
        assert default_power_law.p_min == 1.4
        assert default_power_law.p_max == 3e6
        assert hasattr(default_power_law, 'norm')
        assert default_power_law.norm > 0

    def test_initialization_custom(self, custom_power_law):
        """Test custom initialization of PowerLawPeriod."""
        assert custom_power_law.slope == -1.0
        assert custom_power_law.p_min == 10.0
        assert custom_power_law.p_max == 1e5
        assert hasattr(custom_power_law, 'norm')
        assert custom_power_law.norm > 0

    def test_initialization_invalid_parameters(self):
        """Test that initialization raises ValueError for invalid parameters."""
        # Test p_min <= 0
        with pytest.raises(ValueError, match="p_min must be positive"):
            PowerLawPeriod(p_min=0.0, p_max=1000.0)

        with pytest.raises(ValueError, match="p_min must be positive"):
            PowerLawPeriod(p_min=-1.0, p_max=1000.0)

        # Test p_max <= p_min
        with pytest.raises(ValueError, match="p_max must be greater than p_min"):
            PowerLawPeriod(p_min=1000.0, p_max=100.0)

        with pytest.raises(ValueError, match="p_max must be greater than p_min"):
            PowerLawPeriod(p_min=100.0, p_max=100.0)

    def test_repr(self, default_power_law, custom_power_law):
        """Test string representation of the distribution."""
        default_repr = default_power_law.__repr__()
        assert "PowerLawPeriod(" in default_repr
        assert "slope=-0.55" in default_repr
        assert "p_min=1.4" in default_repr
        assert "p_max=3000000.0" in default_repr

        custom_repr = custom_power_law.__repr__()
        assert "PowerLawPeriod(" in custom_repr
        assert "slope=-1.0" in custom_repr
        assert "p_min=10.0" in custom_repr
        assert "p_max=100000.0" in custom_repr

    def test_repr_html(self, default_power_law):
        """Test HTML representation for Jupyter notebooks."""
        html_str = default_power_law._repr_html_()
        assert "<h3>Power Law Period Distribution</h3>" in html_str
        assert "slope = -0.55" in html_str
        assert "p_min = 1.4" in html_str
        assert "p_max = 3000000.0" in html_str

    def test_power_law_period_method(self, custom_power_law):
        """Test the power_law_period method."""
        logp_values = np.array([1.0, 2.0, 3.0])
        result = custom_power_law.power_law_period(logp_values)

        # Calculate expected values
        p_values = 10**logp_values
        expected = np.zeros_like(p_values)
        valid = (p_values >= custom_power_law.p_min) & (p_values <= custom_power_law.p_max)
        expected[valid] = p_values[valid]**custom_power_law.slope

        np.testing.assert_allclose(result, expected)

    def test_power_law_period_outside_range(self, custom_power_law):
        """Test power_law_period method with values outside valid range."""
        # Below range
        logp_below = np.log10(5.0)  # p_min is 10.0
        result_below = custom_power_law.power_law_period(logp_below)
        assert result_below == 0.0

        # Above range
        logp_above = np.log10(2e5)  # p_max is 1e5
        result_above = custom_power_law.power_law_period(logp_above)
        assert result_above == 0.0

    def test_calculate_normalization(self, custom_power_law):
        """Test that normalization calculation works correctly."""
        # Verify that normalization is positive
        assert custom_power_law.norm > 0

        # Test that integral with normalization is reasonable
        integral, _ = quad(custom_power_law.power_law_period,
                          np.log10(custom_power_law.p_min),
                          np.log10(custom_power_law.p_max))
        expected_norm = 1.0 / integral
        np.testing.assert_allclose(custom_power_law.norm, expected_norm)

    def test_pdf_within_range(self, custom_power_law):
        """Test PDF returns correct values within the period range."""
        p_values = np.logspace(np.log10(custom_power_law.p_min),
                              np.log10(custom_power_law.p_max), 10)
        pdf_values = custom_power_law.pdf(p_values)

        # All should be positive within range
        assert np.all(pdf_values > 0)

        # Calculate expected values
        logp = np.log10(p_values)
        expected = custom_power_law.power_law_period(logp) * custom_power_law.norm
        np.testing.assert_allclose(pdf_values, expected)

    def test_pdf_outside_range(self, custom_power_law):
        """Test PDF returns zero for periods outside the range."""
        # Below range
        p_below = np.array([1.0, 5.0])
        pdf_below = custom_power_law.pdf(p_below)
        np.testing.assert_array_equal(pdf_below, np.zeros_like(p_below))

        # Above range
        p_above = np.array([2e5, 5e5])
        pdf_above = custom_power_law.pdf(p_above)
        np.testing.assert_array_equal(pdf_above, np.zeros_like(p_above))

    def test_pdf_mixed_range(self, custom_power_law):
        """Test PDF with mixture of values inside and outside range."""
        p_mixed = np.array([5.0, 50.0, 500.0, 2e5])  # Some in, some out of range
        pdf_mixed = custom_power_law.pdf(p_mixed)

        # First and last should be zero (outside range)
        assert pdf_mixed[0] == 0.0
        assert pdf_mixed[3] == 0.0

        # Middle values should be positive (inside range)
        assert pdf_mixed[1] > 0.0
        assert pdf_mixed[2] > 0.0

    def test_pdf_scalar_input(self, default_power_law):
        """Test PDF with scalar input."""
        # Within range
        p = 100.0
        pdf_value = default_power_law.pdf(p)
        assert pdf_value > 0

        # Outside range (below)
        p_outside = 0.1
        pdf_outside = default_power_law.pdf(p_outside)
        assert pdf_outside == 0.0

    def test_pdf_zero_and_negative_periods(self, default_power_law):
        """Test PDF handles zero and negative periods correctly."""
        # Test with zero period
        p_zero = 0.0
        pdf_zero = default_power_law.pdf(p_zero)
        assert pdf_zero == 0.0

        # Test with negative period
        p_negative = -10.0
        pdf_negative = default_power_law.pdf(p_negative)
        assert pdf_negative == 0.0

        # Test with array containing zero and negative
        p_mixed = np.array([-5.0, 0.0, 10.0, 100.0])
        pdf_mixed = default_power_law.pdf(p_mixed)
        assert pdf_mixed[0] == 0.0  # Negative
        assert pdf_mixed[1] == 0.0  # Zero
        assert pdf_mixed[2] > 0.0   # Positive within range
        assert pdf_mixed[3] > 0.0   # Positive within range

    def test_normalization_consistency(self, default_power_law):
        """Test that normalization is applied consistently."""
        # Test that the normalization constant matches expected calculation
        integral, _ = quad(default_power_law.power_law_period,
                          np.log10(default_power_law.p_min),
                          np.log10(default_power_law.p_max))
        expected_norm = 1.0 / integral
        np.testing.assert_allclose(default_power_law.norm, expected_norm, rtol=1e-10)

        # Test that PDF values are properly normalized
        p_test = np.array([10.0, 100.0, 1000.0])
        pdf_values = default_power_law.pdf(p_test)

        # Calculate expected PDF values manually
        logp_test = np.log10(p_test)
        raw_values = default_power_law.power_law_period(logp_test)
        expected_pdf = raw_values * default_power_law.norm

        np.testing.assert_allclose(pdf_values, expected_pdf)

    def test_rvs(self, custom_power_law):
        """Test random sampling."""
        rng = np.random.default_rng(42)

        samples = custom_power_law.rvs(size=1000, rng=rng)

        assert len(samples) == 1000
        assert np.all(samples >= custom_power_law.p_min)
        assert np.all(samples <= custom_power_law.p_max)

    def test_rvs_without_rng(self, custom_power_law):
        """Test random sampling without providing an RNG."""
        samples = custom_power_law.rvs(size=100)

        assert len(samples) == 100
        assert np.all(samples >= custom_power_law.p_min)
        assert np.all(samples <= custom_power_law.p_max)


class TestDistributionComparisons:
    """Test class for comparing distributions and edge cases."""

    def test_different_slope_power_laws(self):
        """Test that different slopes give different distributions."""
        steep_slope = PowerLawPeriod(slope=-2.0, p_min=1.0, p_max=1000.0)
        shallow_slope = PowerLawPeriod(slope=-0.5, p_min=1.0, p_max=1000.0)

        p_test = np.array([10.0, 100.0])
        pdf_steep = steep_slope.pdf(p_test)
        pdf_shallow = shallow_slope.pdf(p_test)

        # They should be different
        assert not np.allclose(pdf_steep, pdf_shallow)

        # Steep slope should favor smaller periods more
        assert pdf_steep[0] / pdf_steep[1] > pdf_shallow[0] / pdf_shallow[1]

    def test_mass_dependence_sana12(self):
        """Test that Sana12 period distribution depends on mass."""
        sana12 = Sana12Period()
        p_test = np.array([1.0, 10.0, 100.0])

        # Low mass vs high mass
        m1_low = 10.0  # Below break
        m1_high = 25.0  # Above break

        pdf_low = sana12.pdf(p_test, m1_low)
        pdf_high = sana12.pdf(p_test, m1_high)

        # Should be different for different masses
        assert not np.allclose(pdf_low, pdf_high)

    def test_extreme_parameter_values(self):
        """Test distributions with extreme but valid parameter values."""
        # Very narrow range for FlatMassRatio
        narrow_flat = FlatMassRatio(q_min=0.49, q_max=0.51)
        assert narrow_flat.norm > 0

        # Very wide range for PowerLawPeriod
        wide_power = PowerLawPeriod(p_min=0.1, p_max=1e8)
        assert wide_power.norm > 0

        # Very small period range for Sana12
        narrow_sana = Sana12Period(p_min=0.1, p_max=1.0)
        assert narrow_sana.low_mass_norm > 0
        assert narrow_sana.high_mass_norm > 0

    @pytest.mark.parametrize("distribution_class,default_params,test_params", [
        (FlatMassRatio, {}, {"q_min": 0.2, "q_max": 0.9}),
        (PowerLawPeriod, {}, {"slope": -1.5, "p_min": 5.0, "p_max": 1e4}),
        (Sana12Period, {}, {"p_min": 1.0, "p_max": 500.0}),
    ])
    def test_parameter_consistency(self, distribution_class, default_params, test_params):
        """Test that distribution parameters are set correctly."""
        # Test default parameters
        dist_default = distribution_class(**default_params)

        # Test custom parameters
        dist_custom = distribution_class(**test_params)

        # Verify parameters are set correctly
        for param, value in test_params.items():
            assert getattr(dist_custom, param) == value


class TestLogUniform:
    """Test class for LogUniform distribution."""

    @pytest.fixture
    def default_log_uniform(self):
        """Fixture for default LogUniform instance."""
        return LogUniform()

    def test_initialization_default(self, default_log_uniform):
        """Test default initialization of LogUniform."""
        assert default_log_uniform.min == 5.0
        assert default_log_uniform.max == 1e5
        assert hasattr(default_log_uniform, 'norm')
        assert default_log_uniform.norm > 0

    def test_initialization_custom(self):
        """Test custom initialization of LogUniform."""
        log_uniform = LogUniform(min=10.0, max=1000.0)
        assert log_uniform.min == 10.0
        assert log_uniform.max == 1000.0
        assert hasattr(log_uniform, 'norm')
        assert log_uniform.norm > 0
        # Verify normalization
        norm = 1.0 / (np.log10(log_uniform.max) - np.log10(log_uniform.min))
        np.testing.assert_allclose(log_uniform.norm, norm)

    def test_initialization_invalid_parameters(self):
        """Test that initialization raises ValueError for invalid parameters."""
        # Test min <= 0
        with pytest.raises(ValueError, match="min must be positive"):
            LogUniform(min=0.0, max=1000.0)

        with pytest.raises(ValueError, match="min must be positive"):
            LogUniform(min=-1.0, max=1000.0)

        # Test max <= min
        with pytest.raises(ValueError, match="max must be greater than min"):
            LogUniform(min=1000.0, max=100.0)

        with pytest.raises(ValueError, match="max must be greater than min"):
            LogUniform(min=100.0, max=100.0)

    def test_repr(self):
        """Test string representation."""
        log_uniform = LogUniform(min=10.0, max=1000.0)
        rep_str = log_uniform.__repr__()
        assert "LogUniform(" in rep_str
        assert "min=10.0" in rep_str
        assert "max=1000.0" in rep_str

    def test_repr_html(self):
        """Test HTML representation."""
        log_uniform = LogUniform(min=10.0, max=1000.0)
        html_str = log_uniform._repr_html_()
        assert "<h3>Log-Uniform Distribution</h3>" in html_str
        assert "min = 10.0" in html_str
        assert "max = 1000.0" in html_str

    def test_pdf_within_range(self):
        """Test PDF within the valid range."""
        log_uniform = LogUniform(min=10.0, max=1000.0)
        x_values = np.array([10.0, 50.0, 100.0, 500.0, 1000.0])
        pdf_values = log_uniform.pdf(x_values)
        expected = log_uniform.norm / x_values
        np.testing.assert_allclose(pdf_values, expected)

    def test_pdf_outside_range(self):
        """Test PDF outside the valid range."""
        log_uniform = LogUniform(min=10.0, max=1000.0)
        # Below range
        x_below = np.array([1.0, 5.0])
        pdf_below = log_uniform.pdf(x_below)
        np.testing.assert_array_equal(pdf_below, np.zeros_like(x_below))

        # Above range
        x_above = np.array([2000.0, 5000.0])
        pdf_above = log_uniform.pdf(x_above)
        np.testing.assert_array_equal(pdf_above, np.zeros_like(x_above))

    def test_rvs(self):
        """Test random sampling."""
        log_uniform = LogUniform(min=10.0, max=1000.0)
        rng = np.random.default_rng(42)
        samples = log_uniform.rvs(size=1000, rng=rng)

        assert len(samples) == 1000
        assert np.all(samples >= log_uniform.min)
        assert np.all(samples <= log_uniform.max)

    def test_rvs_without_rng(self):
        """Test random sampling without providing an RNG."""
        log_uniform = LogUniform(min=10.0, max=1000.0)
        samples = log_uniform.rvs(size=100)

        assert len(samples) == 100
        assert np.all(samples >= log_uniform.min)
        assert np.all(samples <= log_uniform.max)


class TestThermalEccentricity:
    """Test class for ThermalEccentricity distribution."""

    @pytest.fixture
    def default_thermal(self):
        """Fixture for default ThermalEccentricity instance."""
        return ThermalEccentricity()

    @pytest.fixture
    def custom_thermal(self):
        """Fixture for custom ThermalEccentricity instance."""
        return ThermalEccentricity(e_min=0.1, e_max=0.9)

    def test_initialization_default(self, default_thermal):
        """Test default initialization."""
        assert default_thermal.e_min == 0.0
        assert default_thermal.e_max == 1.0
        assert hasattr(default_thermal, 'norm')
        assert default_thermal.norm > 0

    def test_initialization_custom(self, custom_thermal):
        """Test custom initialization."""
        assert custom_thermal.e_min == 0.1
        assert custom_thermal.e_max == 0.9
        assert hasattr(custom_thermal, 'norm')
        assert custom_thermal.norm > 0

    def test_initialization_invalid_parameters(self):
        """Test that initialization raises ValueError for invalid parameters."""
        # Test e_min not in [0, 1)
        with pytest.raises(ValueError, match="e_min must be in \\[0, 1\\)"):
            ThermalEccentricity(e_min=-0.1, e_max=0.5)

        with pytest.raises(ValueError, match="e_min must be in \\[0, 1\\)"):
            ThermalEccentricity(e_min=1.5, e_max=2.0)

        # Test e_max not in (0, 1]
        with pytest.raises(ValueError, match="e_max must be in \\(0, 1\\]"):
            ThermalEccentricity(e_min=0.1, e_max=0.0)

        with pytest.raises(ValueError, match="e_max must be in \\(0, 1\\]"):
            ThermalEccentricity(e_min=0.1, e_max=-0.1)

        with pytest.raises(ValueError, match="e_max must be in \\(0, 1\\]"):
            ThermalEccentricity(e_min=0.1, e_max=1.5)

        # Test e_min >= e_max
        with pytest.raises(ValueError, match="e_min must be less than e_max"):
            ThermalEccentricity(e_min=0.8, e_max=0.5)

        with pytest.raises(ValueError, match="e_min must be less than e_max"):
            ThermalEccentricity(e_min=0.5, e_max=0.5)

    def test_repr(self, custom_thermal):
        """Test string representation."""
        rep_str = custom_thermal.__repr__()
        assert "ThermalEccentricity(" in rep_str
        assert "e_min=0.1" in rep_str
        assert "e_max=0.9" in rep_str

    def test_repr_html(self, custom_thermal):
        """Test HTML representation."""
        html_str = custom_thermal._repr_html_()
        assert "<h3>Thermal Eccentricity Distribution</h3>" in html_str
        assert "e_min = 0.1" in html_str
        assert "e_max = 0.9" in html_str

    def test_thermal_eccentricity_method(self, default_thermal):
        """Test the thermal_eccentricity method."""
        e_values = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
        result = default_thermal.thermal_eccentricity(e_values)
        expected = 2.0 * e_values
        np.testing.assert_allclose(result, expected)

    def test_pdf_within_range(self, custom_thermal):
        """Test PDF within the valid range."""
        e_values = np.linspace(custom_thermal.e_min, custom_thermal.e_max, 10)
        pdf_values = custom_thermal.pdf(e_values)

        # All should be positive within range
        assert np.all(pdf_values > 0)

        # Check normalization
        expected = custom_thermal.thermal_eccentricity(e_values) * custom_thermal.norm
        np.testing.assert_allclose(pdf_values, expected)

    def test_pdf_outside_range(self, custom_thermal):
        """Test PDF outside the valid range."""
        # Below range
        e_below = np.array([0.0, 0.05])
        pdf_below = custom_thermal.pdf(e_below)
        np.testing.assert_array_equal(pdf_below, np.zeros_like(e_below))

        # Above range
        e_above = np.array([0.95, 1.0])
        pdf_above = custom_thermal.pdf(e_above)
        np.testing.assert_array_equal(pdf_above, np.zeros_like(e_above))

    def test_pdf_scalar_input(self, default_thermal):
        """Test PDF with scalar input."""
        e = 0.5
        pdf_value = default_thermal.pdf(e)
        expected = default_thermal.thermal_eccentricity(e) * default_thermal.norm
        np.testing.assert_allclose(pdf_value, expected)

    def test_rvs(self, custom_thermal):
        """Test random sampling."""
        rng = np.random.default_rng(42)
        samples = custom_thermal.rvs(size=1000, rng=rng)

        assert len(samples) == 1000
        assert np.all(samples >= custom_thermal.e_min)
        assert np.all(samples <= custom_thermal.e_max)

    def test_rvs_without_rng(self, custom_thermal):
        """Test random sampling without providing an RNG."""
        samples = custom_thermal.rvs(size=100)

        assert len(samples) == 100
        assert np.all(samples >= custom_thermal.e_min)
        assert np.all(samples <= custom_thermal.e_max)


class TestUniformEccentricity:
    """Test class for UniformEccentricity distribution."""

    @pytest.fixture
    def default_uniform(self):
        """Fixture for default UniformEccentricity instance."""
        return UniformEccentricity()

    @pytest.fixture
    def custom_uniform(self):
        """Fixture for custom UniformEccentricity instance."""
        return UniformEccentricity(e_min=0.1, e_max=0.9)

    def test_initialization_default(self, default_uniform):
        """Test default initialization."""
        assert default_uniform.e_min == 0.0
        assert default_uniform.e_max == 1.0
        assert hasattr(default_uniform, 'norm')
        assert default_uniform.norm == 1.0

    def test_initialization_custom(self, custom_uniform):
        """Test custom initialization."""
        assert custom_uniform.e_min == 0.1
        assert custom_uniform.e_max == 0.9
        assert hasattr(custom_uniform, 'norm')
        expected_norm = 1.0 / (0.9 - 0.1)
        np.testing.assert_allclose(custom_uniform.norm, expected_norm)

    def test_initialization_invalid_parameters(self):
        """Test that initialization raises ValueError for invalid parameters."""
        # Test e_min not in [0, 1)
        with pytest.raises(ValueError, match="e_min must be in \\[0, 1\\)"):
            UniformEccentricity(e_min=-0.1, e_max=0.5)

        with pytest.raises(ValueError, match="e_min must be in \\[0, 1\\)"):
            UniformEccentricity(e_min=1.5, e_max=2.0)

        # Test e_max not in (0, 1]
        with pytest.raises(ValueError, match="e_max must be in \\(0, 1\\]"):
            UniformEccentricity(e_min=0.1, e_max=0.0)

        with pytest.raises(ValueError, match="e_max must be in \\(0, 1\\]"):
            UniformEccentricity(e_min=0.1, e_max=-0.1)

        with pytest.raises(ValueError, match="e_max must be in \\(0, 1\\]"):
            UniformEccentricity(e_min=0.1, e_max=1.5)

        # Test e_min >= e_max
        with pytest.raises(ValueError, match="e_min must be less than e_max"):
            UniformEccentricity(e_min=0.8, e_max=0.5)

        with pytest.raises(ValueError, match="e_min must be less than e_max"):
            UniformEccentricity(e_min=0.5, e_max=0.5)

    def test_repr(self, custom_uniform):
        """Test string representation."""
        rep_str = custom_uniform.__repr__()
        assert "UniformEccentricity(" in rep_str
        assert "e_min=0.1" in rep_str
        assert "e_max=0.9" in rep_str

    def test_repr_html(self, custom_uniform):
        """Test HTML representation."""
        html_str = custom_uniform._repr_html_()
        assert "<h3>Uniform Eccentricity Distribution</h3>" in html_str
        assert "e_min = 0.1" in html_str
        assert "e_max = 0.9" in html_str

    def test_pdf_within_range(self, custom_uniform):
        """Test PDF within the valid range."""
        e_values = np.linspace(custom_uniform.e_min, custom_uniform.e_max, 10)
        pdf_values = custom_uniform.pdf(e_values)

        # All should equal the normalization constant
        expected = custom_uniform.norm * np.ones_like(e_values)
        np.testing.assert_allclose(pdf_values, expected)

    def test_pdf_outside_range(self, custom_uniform):
        """Test PDF outside the valid range."""
        # Below range
        e_below = np.array([0.0, 0.05])
        pdf_below = custom_uniform.pdf(e_below)
        np.testing.assert_array_equal(pdf_below, np.zeros_like(e_below))

        # Above range
        e_above = np.array([0.95, 1.0])
        pdf_above = custom_uniform.pdf(e_above)
        np.testing.assert_array_equal(pdf_above, np.zeros_like(e_above))

    def test_pdf_scalar_input(self, default_uniform):
        """Test PDF with scalar input."""
        e = 0.5
        pdf_value = default_uniform.pdf(e)
        np.testing.assert_allclose(pdf_value, default_uniform.norm)

    def test_rvs(self, custom_uniform):
        """Test random sampling."""
        rng = np.random.default_rng(42)
        samples = custom_uniform.rvs(size=1000, rng=rng)

        assert len(samples) == 1000
        assert np.all(samples >= custom_uniform.e_min)
        assert np.all(samples <= custom_uniform.e_max)

    def test_rvs_without_rng(self, custom_uniform):
        """Test random sampling without providing an RNG."""
        samples = custom_uniform.rvs(size=100)

        assert len(samples) == 100
        assert np.all(samples >= custom_uniform.e_min)
        assert np.all(samples <= custom_uniform.e_max)


class TestZeroEccentricity:
    """Test class for ZeroEccentricity distribution."""

    @pytest.fixture
    def zero_ecc(self):
        """Fixture for ZeroEccentricity instance."""
        return ZeroEccentricity()

    def test_initialization(self, zero_ecc):
        """Test initialization."""
        # Should have no parameters
        assert isinstance(zero_ecc, ZeroEccentricity)

    def test_repr(self, zero_ecc):
        """Test string representation."""
        rep_str = zero_ecc.__repr__()
        assert "ZeroEccentricity()" in rep_str

    def test_repr_html(self, zero_ecc):
        """Test HTML representation."""
        html_str = zero_ecc._repr_html_()
        assert "<h3>Zero Eccentricity Distribution</h3>" in html_str
        assert "e = 0 (circular orbits)" in html_str

    def test_pdf_at_zero(self, zero_ecc):
        """Test PDF at e=0."""
        pdf_value = zero_ecc.pdf(0.0)
        assert pdf_value == 1.0

    def test_pdf_away_from_zero(self, zero_ecc):
        """Test PDF for non-zero eccentricities."""
        e_values = np.array([0.1, 0.5, 0.9, 1.0])
        pdf_values = zero_ecc.pdf(e_values)
        np.testing.assert_array_equal(pdf_values, np.zeros_like(e_values))

    def test_pdf_mixed(self, zero_ecc):
        """Test PDF with mixture of zero and non-zero values."""
        e_values = np.array([0.0, 0.1, 0.0, 0.5])
        pdf_values = zero_ecc.pdf(e_values)
        expected = np.array([1.0, 0.0, 1.0, 0.0])
        np.testing.assert_array_equal(pdf_values, expected)

    def test_rvs(self, zero_ecc):
        """Test random sampling."""
        rng = np.random.default_rng(42)
        samples = zero_ecc.rvs(size=1000, rng=rng)

        # All samples should be zero
        assert len(samples) == 1000
        np.testing.assert_array_equal(samples, np.zeros(1000))

    def test_rvs_without_rng(self, zero_ecc):
        """Test random sampling without providing an RNG."""
        samples = zero_ecc.rvs(size=100)

        # All samples should be zero
        assert len(samples) == 100
        np.testing.assert_array_equal(samples, np.zeros(100))


class TestLogNormalSeparation:
    """Test class for LogNormalSeparation distribution."""

    @pytest.fixture
    def default_lognormal(self):
        """Fixture for default LogNormalSeparation instance."""
        return LogNormalSeparation()

    @pytest.fixture
    def custom_lognormal(self):
        """Fixture for custom LogNormalSeparation instance."""
        return LogNormalSeparation(mean=1.0, sigma=0.5, min=10.0, max=1e4)

    def test_initialization_default(self, default_lognormal):
        """Test default initialization."""
        assert default_lognormal.mean == 0.85
        assert default_lognormal.sigma == 0.37
        assert default_lognormal.min == 5.0
        assert default_lognormal.max == 1e5

    def test_initialization_custom(self, custom_lognormal):
        """Test custom initialization."""
        assert custom_lognormal.mean == 1.0
        assert custom_lognormal.sigma == 0.5
        assert custom_lognormal.min == 10.0
        assert custom_lognormal.max == 1e4

    def test_initialization_invalid_parameters(self):
        """Test that initialization raises ValueError for invalid parameters."""
        # Test min <= 0
        with pytest.raises(ValueError, match="min must be positive"):
            LogNormalSeparation(mean=1.0, sigma=0.5, min=0.0, max=1000.0)

        with pytest.raises(ValueError, match="min must be positive"):
            LogNormalSeparation(mean=1.0, sigma=0.5, min=-1.0, max=1000.0)

        # Test max <= min
        with pytest.raises(ValueError, match="max must be greater than min"):
            LogNormalSeparation(mean=1.0, sigma=0.5, min=1000.0, max=100.0)

        with pytest.raises(ValueError, match="max must be greater than min"):
            LogNormalSeparation(mean=1.0, sigma=0.5, min=100.0, max=100.0)

        # Test sigma <= 0
        with pytest.raises(ValueError, match="sigma must be positive"):
            LogNormalSeparation(mean=1.0, sigma=0.0, min=10.0, max=1000.0)

        with pytest.raises(ValueError, match="sigma must be positive"):
            LogNormalSeparation(mean=1.0, sigma=-0.5, min=10.0, max=1000.0)

    def test_repr(self, custom_lognormal):
        """Test string representation."""
        rep_str = custom_lognormal.__repr__()
        assert "LogNormalSeparation(" in rep_str
        assert "mean=1.0" in rep_str
        assert "sigma=0.5" in rep_str
        assert "min=10.0" in rep_str
        assert "max=10000.0" in rep_str

    def test_repr_html(self, custom_lognormal):
        """Test HTML representation."""
        html_str = custom_lognormal._repr_html_()
        assert "<h3>Log-Normal Separation Distribution</h3>" in html_str
        assert "mean (log10) = 1.0" in html_str
        assert "sigma (log10) = 0.5" in html_str

    def test_pdf_within_range(self, custom_lognormal):
        """Test PDF within the valid range."""
        a_values = np.array([10.0, 50.0, 100.0, 500.0, 1000.0])
        pdf_values = custom_lognormal.pdf(a_values)

        # All should be positive within range
        assert np.all(pdf_values > 0)

    def test_pdf_outside_range(self, custom_lognormal):
        """Test PDF outside the valid range."""
        # Below range
        a_below = np.array([1.0, 5.0])
        pdf_below = custom_lognormal.pdf(a_below)
        np.testing.assert_array_equal(pdf_below, np.zeros_like(a_below))

        # Above range
        a_above = np.array([2e4, 5e4])
        pdf_above = custom_lognormal.pdf(a_above)
        np.testing.assert_array_equal(pdf_above, np.zeros_like(a_above))

    def test_pdf_zero_and_negative(self, custom_lognormal):
        """Test PDF for zero and negative values."""
        a_invalid = np.array([0.0, -10.0])
        pdf_invalid = custom_lognormal.pdf(a_invalid)
        np.testing.assert_array_equal(pdf_invalid, np.zeros_like(a_invalid))

    def test_rvs(self, custom_lognormal):
        """Test random sampling."""
        rng = np.random.default_rng(42)
        samples = custom_lognormal.rvs(size=1000, rng=rng)

        assert len(samples) == 1000
        assert np.all(samples >= custom_lognormal.min)
        assert np.all(samples <= custom_lognormal.max)

    def test_rvs_without_rng(self, custom_lognormal):
        """Test random sampling without providing an RNG."""
        samples = custom_lognormal.rvs(size=100)

        assert len(samples) == 100
        assert np.all(samples >= custom_lognormal.min)
        assert np.all(samples <= custom_lognormal.max)


        with pytest.raises(ValueError, match="max must be greater than min"):
            LogUniform(min=100.0, max=100.0)
