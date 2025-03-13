import pytest
import numpy as np
from posydon.popsyn import norm_pop
import posydon.popsyn.IMFs as IMFs

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
        assert np.all(result == 1)
    
    def test_valid_imf_returns_dummy_pdf(self):
        # Monkey-patch: add DummyIMF to IMFs temporarily
        IMFs.DummyIMF = DummyIMF
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
        # cleanup monkey-patch
        del IMFs.DummyIMF
