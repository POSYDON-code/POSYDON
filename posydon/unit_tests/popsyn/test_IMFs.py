import pytest
import numpy as np
import posydon.popsyn.IMFs as IMFs
from scipy.integrate import quad

class TestSalpeterIMF:
    @pytest.fixture
    def default_imf(self):
        """Fixture for default SalpeterIMF instance."""
        return IMFs.Salpeter()

    @pytest.fixture
    def custom_imf(self):
        """Fixture for custom SalpeterIMF instance."""
        return IMFs.Salpeter(alpha=2.5, m_min=0.5, m_max=100.0)

    def test_initialization_default(self, default_imf):
        """Test default initialization of SalpeterIMF."""
        assert default_imf.alpha == 2.35
        assert default_imf.m_min == 0.01
        assert default_imf.m_max == 200.0
        # Verify normalization
        integral, _ = quad(default_imf.pdf, default_imf.m_min, default_imf.m_max)
        assert integral == pytest.approx(1.0, rel=1e-4)

    def test_initialization_custom(self, custom_imf):
        """Test custom initialization of SalpeterIMF."""
        assert custom_imf.alpha == 2.5
        assert custom_imf.m_min == 0.5
        assert custom_imf.m_max == 100.0
        # Verify normalization
        integral, _ = quad(custom_imf.pdf, custom_imf.m_min, custom_imf.m_max)
        assert integral == pytest.approx(1.0, rel=1e-4)

    def test_salepeter_imf_method(self, default_imf):
        """Test the salpeter_IMF method for correct values."""
        m_values = np.array([1.0, 2.0, 5.0, 10.0])
        expected = m_values ** (-default_imf.alpha)
        computed = default_imf.salpeter_IMF(m_values)
        assert np.all(computed == expected)

    def test_salepeter_imf_scalar(self, default_imf):
        """Test the salpeter_IMF method with a scalar input."""
        m = 3.0
        expected = m ** (-default_imf.alpha)
        computed = default_imf.salpeter_IMF(m)
        assert computed == expected

    def test_salepeter_imf_invalid_mass(self, default_imf):
        """Test that salpeter_IMF raises ValueError for non-positive masses."""
        with pytest.raises(ValueError, match="Mass must be positive."):
            default_imf.salpeter_IMF(0)
        with pytest.raises(ValueError, match="Mass must be positive."):
            default_imf.salpeter_IMF(-1.0)

    def test_pdf_within_range(self, default_imf):
        """Test that PDF returns correct values within the mass range."""
        m = np.linspace(default_imf.m_min, default_imf.m_max, 100)
        pdf_values = default_imf.pdf(m)
        expected_pdf = default_imf.salpeter_IMF(m) * default_imf.norm
        assert np.all(pdf_values == expected_pdf)

    def test_pdf_outside_range(self, default_imf):
        """Test that PDF returns zero for masses outside the mass range."""
        m = np.array([default_imf.m_min - 0.1, default_imf.m_max + 0.1])
        pdf_values = default_imf.pdf(m)
        assert np.all(pdf_values == 0.0)

    def test_pdf_scalar(self, default_imf):
        """Test PDF with a scalar input within range."""
        m = 10.0
        pdf_value = default_imf.pdf(m)
        expected = default_imf.salpeter_IMF(m) * default_imf.norm
        assert pdf_value == expected

    def test_pdf_scalar_outside(self, default_imf):
        """Test PDF with a scalar input outside the range."""
        m = default_imf.m_max + 10.0
        pdf_value = default_imf.pdf(m)
        assert pdf_value == 0.0

    def test_normalization(self, default_imf):
        """Ensure that the integral of the PDF over the range is approximately 1."""
        m_min = default_imf.m_min
        m_max = default_imf.m_max
        integral, error = quad(default_imf.pdf, m_min, m_max)
        assert integral == pytest.approx(1.0, rel=1e-4)

    def test_invalid_initialization(self):
        """Test that initialization raises ValueError when integral is zero."""
        with pytest.raises(ValueError, match="Normalization integral is zero"):
            # Create an IMF with alpha such that the integral becomes zero
            # For example, alpha <= 1 results in a divergent integral as m approaches zero
            IMFs.Salpeter(alpha=1.0, m_min=1e-10, m_max=1e-10)

    def test_vectorization(self, default_imf):
        """Test that the pdf method correctly handles array inputs."""
        m = np.array([0.01, 0.1, 1.0, 10.0, 100.0, 200.0])
        pdf_values = default_imf.pdf(m)
        expected = np.array([
            default_imf.salpeter_IMF(0.01) * default_imf.norm,
            default_imf.salpeter_IMF(0.1) * default_imf.norm,
            default_imf.salpeter_IMF(1.0) * default_imf.norm,
            default_imf.salpeter_IMF(10.0) * default_imf.norm,
            default_imf.salpeter_IMF(100.0) * default_imf.norm,
            default_imf.salpeter_IMF(200.0) * default_imf.norm,
        ])
        assert np.all(pdf_values == expected)

    def test_pdf_integral_custom(self, custom_imf):
        """Ensure that the integral of the PDF for a custom IMF is approximately 1."""
        integral, _ = quad(custom_imf.pdf, custom_imf.m_min, custom_imf.m_max)
        assert integral == pytest.approx(1.0, rel=1e-4)

    def test_repr(self, default_imf):
        rep_str = default_imf.__repr__()
        assert "Salpeter(" in rep_str
        assert "alpha=" in rep_str

    def test_repr_html(self, default_imf):
        html_str = default_imf._repr_html_()
        assert "<h3>Salpeter IMF</h3>" in html_str

class TestKroupa2001IMF:
    @pytest.fixture
    def default_kroupa(self):
        return IMFs.Kroupa2001()

    @pytest.fixture
    def custom_kroupa(self):
        return IMFs.Kroupa2001(alpha1=0.5, alpha2=1.7, alpha3=2.7, m1break=0.09, m2break=0.6)

    def test_initialization_default(self, default_kroupa):
        assert default_kroupa.alpha1 == 0.3
        assert default_kroupa.alpha2 == 1.3
        assert default_kroupa.alpha3 == 2.3
        assert default_kroupa.m_min == 0.01
        assert default_kroupa.m_max == 200.0
        integral, _ = quad(default_kroupa.pdf, default_kroupa.m_min, default_kroupa.m_max)
        assert integral == pytest.approx(1.0, rel=1e-4)

    def test_initialization_custom(self, custom_kroupa):
        assert custom_kroupa.alpha1 == 0.5
        assert custom_kroupa.alpha2 == 1.7
        assert custom_kroupa.alpha3 == 2.7
        assert custom_kroupa.m1break == 0.09
        assert custom_kroupa.m2break == 0.6
        integral, _ = quad(custom_kroupa.pdf, custom_kroupa.m_min, custom_kroupa.m_max)
        assert integral == pytest.approx(1.0, rel=1e-4)

    def test_pdf_within_range(self, default_kroupa):
        m = np.linspace(default_kroupa.m_min, default_kroupa.m_max, 100)
        pdf_values = default_kroupa.pdf(m)
        expected_pdf = default_kroupa.kroupa_IMF(m) * default_kroupa.norm
        assert np.all(pdf_values == expected_pdf)

    def test_pdf_outside_range(self, default_kroupa):
        m = np.array([default_kroupa.m_min - 0.1, default_kroupa.m_max + 0.1])
        pdf_values = default_kroupa.pdf(m)
        assert np.all(pdf_values == 0.0)

    def test_pdf_scalar(self, default_kroupa):
        m = 1.0
        pdf_value = default_kroupa.pdf(m)
        expected = default_kroupa.kroupa_IMF(m) * default_kroupa.norm
        assert pdf_value == expected

    def test_pdf_scalar_outside(self, default_kroupa):
        m = default_kroupa.m_max + 10.0
        pdf_value = default_kroupa.pdf(m)
        assert pdf_value == 0.0

    def test_normalization(self, default_kroupa):
        integral, _ = quad(default_kroupa.pdf, default_kroupa.m_min, default_kroupa.m_max)
        assert integral == pytest.approx(1.0, rel=1e-4)

    def test_invalid_initialization(self):
        with pytest.raises(ValueError, match="Normalization integral is zero"):
            IMFs.Kroupa2001(alpha1=1.0, alpha2=1.0, alpha3=1.0, m_min=1e-10, m_max=1e-10)

    def test_vectorization(self, default_kroupa):
        m_array = np.array([0.01, 0.1, 1.0, 10.0, 100.0, 200.0])
        pdf_values = default_kroupa.pdf(m_array)
        expected = default_kroupa.kroupa_IMF(m_array) * default_kroupa.norm
        assert np.all(pdf_values == expected)

    def test_pdf_integral_custom(self, custom_kroupa):
        integral, _ = quad(custom_kroupa.pdf, custom_kroupa.m_min, custom_kroupa.m_max)
        assert integral == pytest.approx(1.0, rel=1e-4)
        
    def test_negative_mass(self, default_kroupa):
        with pytest.raises(ValueError, match="Mass must be positive."):
            default_kroupa.kroupa_IMF(-1.0)
        with pytest.raises(ValueError, match="Mass must be positive."):
            default_kroupa.kroupa_IMF(0.0)

    def test_repr(self, default_kroupa):
        rep_str = default_kroupa.__repr__()
        assert "Kroupa2001(" in rep_str
        assert "alpha1=" in rep_str

    def test_repr_html(self, default_kroupa):
        html_str = default_kroupa._repr_html_()
        assert "<h3>Kroupa (2001) IMF</h3>" in html_str

