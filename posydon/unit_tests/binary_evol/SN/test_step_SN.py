"""Unit tests of posydon/binary_evol/SN/step_SN.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>"
]

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import approx, fixture

# import the module which will be tested
import posydon.binary_evol.SN.step_SN as totest
from posydon.binary_evol.singlestar import SingleStar
from posydon.utils.constants import Zsun


class TestPISNPrescription:
    @fixture
    def SN_step(self):
        # Hendriks+23 without extra mass loss and the core collapse replaced
        # by keeping the mass after PPI
        step = totest.StepSN(PISN='Hendriks+23', PPI_extra_mass_loss=0.0)
        step.compute_m_rembar = lambda star, m_PISN: (m_PISN, None, None)
        return step

    def get_m_PISN(self, step, Z_div_Zsun):
        star = SingleStar(mass=70.0, he_core_mass=60.0, co_core_mass=50.0)
        star.metallicity = Z_div_Zsun
        return step.PISN_prescription(star)

    def test_Hendriks23_metallicity(self, SN_step):
        # star.metallicity is Z/Zsun, while the fit of Hendriks et al. 2023
        # uses the absolute Z
        Z = 0.001
        delta_M_PPI = ((0.0006 * totest.np.log10(Z) + 0.0054) * (50.0-34.8)**3
                       - 0.0013 * (50.0-34.8)**2)
        assert self.get_m_PISN(SN_step, Z/Zsun) == approx(60.0 - delta_M_PPI)

    def test_Hendriks23_metallicity_floor(self, SN_step):
        # the fit is limited to the lowest absolute Z = 1e-4 in Hendriks+23
        m_PISN_floor = self.get_m_PISN(SN_step, 1e-4/Zsun)
        assert self.get_m_PISN(SN_step, 1e-6) == approx(m_PISN_floor)
        assert self.get_m_PISN(SN_step, 2e-4/Zsun) < m_PISN_floor
