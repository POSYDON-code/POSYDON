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
        """Test the imf method for correct values."""
        m_values = np.array([1.0, 2.0, 5.0, 10.0])
        expected = m_values ** (-default_imf.alpha)
        computed = default_imf.imf(m_values)
        assert np.all(computed == expected)

    def test_salepeter_imf_scalar(self, default_imf):
        """Test the imf method with a scalar input."""
        m = 3.0
        expected = m ** (-default_imf.alpha)
        computed = default_imf.imf(m)
        assert computed == expected

    def test_salepeter_imf_invalid_mass(self, default_imf):
        """Test that imf raises ValueError for non-positive masses."""
        with pytest.raises(ValueError, match="Mass must be positive."):
            default_imf.imf(0)
        with pytest.raises(ValueError, match="Mass must be positive."):
            default_imf.imf(-1.0)

    def test_pdf_within_range(self, default_imf):
        """Test that PDF returns correct values within the mass range."""
        m = np.linspace(default_imf.m_min, default_imf.m_max, 100)
        pdf_values = default_imf.pdf(m)
        expected_pdf = default_imf.imf(m) * default_imf.norm
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
        expected = default_imf.imf(m) * default_imf.norm
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
            default_imf.imf(0.01) * default_imf.norm,
            default_imf.imf(0.1) * default_imf.norm,
            default_imf.imf(1.0) * default_imf.norm,
            default_imf.imf(10.0) * default_imf.norm,
            default_imf.imf(100.0) * default_imf.norm,
            default_imf.imf(200.0) * default_imf.norm,
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
        expected_pdf = default_kroupa.imf(m) * default_kroupa.norm
        assert np.all(pdf_values == expected_pdf)

    def test_pdf_outside_range(self, default_kroupa):
        m = np.array([default_kroupa.m_min - 0.1, default_kroupa.m_max + 0.1])
        pdf_values = default_kroupa.pdf(m)
        assert np.all(pdf_values == 0.0)

    def test_pdf_scalar(self, default_kroupa):
        m = 1.0
        pdf_value = default_kroupa.pdf(m)
        expected = default_kroupa.imf(m) * default_kroupa.norm
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
        expected = default_kroupa.imf(m_array) * default_kroupa.norm
        assert np.all(pdf_values == expected)

    def test_pdf_integral_custom(self, custom_kroupa):
        integral, _ = quad(custom_kroupa.pdf, custom_kroupa.m_min, custom_kroupa.m_max)
        assert integral == pytest.approx(1.0, rel=1e-4)
        
    def test_negative_mass(self, default_kroupa):
        with pytest.raises(ValueError, match="Mass must be positive."):
            default_kroupa.imf(-1.0)
        with pytest.raises(ValueError, match="Mass must be positive."):
            default_kroupa.imf(0.0)

    def test_repr(self, default_kroupa):
        rep_str = default_kroupa.__repr__()
        assert "Kroupa2001(" in rep_str
        assert "alpha1=" in rep_str

    def test_repr_html(self, default_kroupa):
        html_str = default_kroupa._repr_html_()
        assert "<h3>Kroupa (2001) IMF</h3>" in html_str

class TestChabrierIMF:
    @pytest.fixture
    def default_chabrier(self):
        """Fixture for default Chabrier2003 instance."""
        return IMFs.Chabrier2003()

    @pytest.fixture
    def custom_chabrier(self):
        """Fixture for custom Chabrier2003 instance."""
        return IMFs.Chabrier2003(m_c=0.3, sigma=0.6, alpha=2.5, m_break=1.2, m_min=0.05, m_max=150.0)

    def test_initialization_default(self, default_chabrier):
        """Test default initialization of Chabrier2003 IMF."""
        assert default_chabrier.m_c == 0.22
        assert default_chabrier.sigma == 0.57
        assert default_chabrier.alpha == 2.3
        assert default_chabrier.m_break == 1.0
        assert default_chabrier.m_min == 0.01
        assert default_chabrier.m_max == 200.0
        integral, _ = quad(default_chabrier.pdf, default_chabrier.m_min, default_chabrier.m_max)
        assert integral == pytest.approx(1.0, rel=1e-4)

    def test_initialization_custom(self, custom_chabrier):
        """Test custom initialization of Chabrier2003 IMF."""
        assert custom_chabrier.m_c == 0.3
        assert custom_chabrier.sigma == 0.6
        assert custom_chabrier.alpha == 2.5
        assert custom_chabrier.m_break == 1.2
        assert custom_chabrier.m_min == 0.05
        assert custom_chabrier.m_max == 150.0
        integral, _ = quad(custom_chabrier.pdf, custom_chabrier.m_min, custom_chabrier.m_max)
        assert integral == pytest.approx(1.0, rel=1e-4)

    def test_pdf_within_range(self, default_chabrier):
        """Test that PDF returns correct values within the mass range."""
        m = np.linspace(default_chabrier.m_min, default_chabrier.m_max, 100)
        pdf_values = default_chabrier.pdf(m)
        expected_pdf = default_chabrier.imf(m) * default_chabrier.norm
        assert np.all(pdf_values == expected_pdf)

    def test_pdf_outside_range(self, default_chabrier):
        """Test that PDF returns zero for masses outside the mass range."""
        m = np.array([default_chabrier.m_min - 0.1, default_chabrier.m_max + 0.1])
        pdf_values = default_chabrier.pdf(m)
        assert np.all(pdf_values == 0.0)

    def test_pdf_scalar(self, default_chabrier):
        """Test PDF with a scalar input within range."""
        m = 0.5  # within m_break for default instance
        pdf_value = default_chabrier.pdf(m)
        expected = default_chabrier.imf(m) * default_chabrier.norm
        assert pdf_value == expected

    def test_pdf_scalar_outside(self, default_chabrier):
        """Test PDF with a scalar input outside the mass range."""
        m = default_chabrier.m_max + 10.0
        pdf_value = default_chabrier.pdf(m)
        assert pdf_value == 0.0

    def test_normalization(self, default_chabrier):
        """Ensure that the integral of the PDF over the range is approximately 1."""
        integral, _ = quad(default_chabrier.pdf, default_chabrier.m_min, default_chabrier.m_max)
        assert integral == pytest.approx(1.0, rel=1e-4)

    def test_vectorization(self, default_chabrier):
        """Test that the pdf method correctly handles array inputs."""
        m_array = np.array([0.01, 0.1, 1.0, 10.0, 100.0, 200.0])
        pdf_values = default_chabrier.pdf(m_array)
        expected = default_chabrier.imf(m_array) * default_chabrier.norm
        assert np.all(pdf_values == expected)

    def test_invalid_initialization(self):
        """Test that initialization raises ValueError when normalization integral is zero."""
        with pytest.raises(ValueError, match="Normalization integral is zero"):
            IMFs.Chabrier2003(m_c=0.22, sigma=0.57, alpha=2.3, m_break=1.0, m_min=1e-10, m_max=1e-10)

    def test_repr(self, default_chabrier):
        rep_str = default_chabrier.__repr__()
        assert "Chabrier" in rep_str
        assert "m_c=" in rep_str

    def test_repr_html(self, default_chabrier):
        html_str = default_chabrier._repr_html_()
        assert "<h3>Chabrier IMF</h3>" in html_str

