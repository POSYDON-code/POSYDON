"""Unit tests of posydon/utils/constants.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the unittest module and the module which will be tested
import unittest
import posydon.utils.constants as totest

# import other needed code for the tests
# not needed

# define test classes
class TestElements(unittest.TestCase):
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['H_weight', 'He_weight', 'Lsun', 'Lsun33', 'Msun',
                    'Msun33', 'Qconv', 'Rsun', 'Rsun11', 'SNcheck_ERR',
                    'Teffsol', 'Zsun', '__authors__', '__builtins__',
                    '__cached__', '__doc__', '__file__', '__loader__',
                    '__name__', '__package__', '__spec__', 'a2rad',
                    'age_of_universe', 'agesol', 'amu', 'asol', 'au', 'aursun',
                    'avo', 'boltz_sigma', 'boltzm', 'cgas', 'clight', 'crad',
                    'day2sec', 'dayyer', 'ev2erg', 'fine', 'hbar', 'hion',
                    'inversecm2erg', 'kerg', 'kev', 'km2cm', 'loggsol', 'lsol',
                    'ly', 'm_earth', 'm_jupiter', 'mbolsol', 'mbolsun', 'me',
                    'mev_amu', 'mev_to_ergs', 'mn', 'mp', 'msol', 'pc', 'pi',
                    'planck_h', 'qe', 'r_earth', 'r_jupiter', 'rad2a', 'rbohr',
                    'rhonuc', 'rsol', 'secyer', 'semimajor_axis_jupiter',
                    'ssol', 'standard_cgrav', 'weinfre', 'weinlam']
        self.assertListEqual(dir(totest), elements,
                             msg="There might be added or removed objects "
                                 "without an update on the unit test.")

    def test_instance_pi(self):
        self.assertIsInstance(totest.pi, (float, int))

    def test_instance_a2rad(self):
        self.assertIsInstance(totest.a2rad, (float, int))

    def test_instance_rad2a(self):
        self.assertIsInstance(totest.rad2a, (float, int))

    def test_instance_standard_cgrav(self):
        self.assertIsInstance(totest.standard_cgrav, (float, int))

    def test_instance_planck_h(self):
        self.assertIsInstance(totest.planck_h, (float, int))

    def test_instance_hbar(self):
        self.assertIsInstance(totest.hbar, (float, int))

    def test_instance_qe(self):
        self.assertIsInstance(totest.qe, (float, int))

    def test_instance_avo(self):
        self.assertIsInstance(totest.avo, (float, int))

    def test_instance_clight(self):
        self.assertIsInstance(totest.clight, (float, int))

    def test_instance_kerg(self):
        self.assertIsInstance(totest.kerg, (float, int))

    def test_instance_boltzm(self):
        self.assertIsInstance(totest.boltzm, (float, int))

    def test_instance_cgas(self):
        self.assertIsInstance(totest.cgas, (float, int))

    def test_instance_kev(self):
        self.assertIsInstance(totest.kev, (float, int))

    def test_instance_amu(self):
        self.assertIsInstance(totest.amu, (float, int))

    def test_instance_mn(self):
        self.assertIsInstance(totest.mn, (float, int))

    def test_instance_mp(self):
        self.assertIsInstance(totest.mp, (float, int))

    def test_instance_me(self):
        self.assertIsInstance(totest.me, (float, int))

    def test_instance_rbohr(self):
        self.assertIsInstance(totest.rbohr, (float, int))

    def test_instance_fine(self):
        self.assertIsInstance(totest.fine, (float, int))

    def test_instance_hion(self):
        self.assertIsInstance(totest.hion, (float, int))

    def test_instance_ev2erg(self):
        self.assertIsInstance(totest.ev2erg, (float, int))

    def test_instance_inversecm2erg(self):
        self.assertIsInstance(totest.inversecm2erg, (float, int))

    def test_instance_mev_to_ergs(self):
        self.assertIsInstance(totest.mev_to_ergs, (float, int))

    def test_instance_mev_amu(self):
        self.assertIsInstance(totest.mev_amu, (float, int))

    def test_instance_Qconv(self):
        self.assertIsInstance(totest.Qconv, (float, int))

    def test_instance_boltz_sigma(self):
        self.assertIsInstance(totest.boltz_sigma, (float, int))

    def test_instance_crad(self):
        self.assertIsInstance(totest.crad, (float, int))

    def test_instance_ssol(self):
        self.assertIsInstance(totest.ssol, (float, int))

    def test_instance_asol(self):
        self.assertIsInstance(totest.asol, (float, int))

    def test_instance_weinlam(self):
        self.assertIsInstance(totest.weinlam, (float, int))

    def test_instance_weinfre(self):
        self.assertIsInstance(totest.weinfre, (float, int))

    def test_instance_rhonuc(self):
        self.assertIsInstance(totest.rhonuc, (float, int))

    def test_instance_Zsun(self):
        self.assertIsInstance(totest.Zsun, (float, int))

    def test_instance_msol(self):
        self.assertIsInstance(totest.msol, (float, int))

    def test_instance_rsol(self):
        self.assertIsInstance(totest.rsol, (float, int))

    def test_instance_lsol(self):
        self.assertIsInstance(totest.lsol, (float, int))

    def test_instance_agesol(self):
        self.assertIsInstance(totest.agesol, (float, int))

    def test_instance_Msun(self):
        self.assertIsInstance(totest.Msun, (float, int))

    def test_instance_Rsun(self):
        self.assertIsInstance(totest.Rsun, (float, int))

    def test_instance_Lsun(self):
        self.assertIsInstance(totest.Lsun, (float, int))

    def test_instance_Msun33(self):
        self.assertIsInstance(totest.Msun33, (float, int))

    def test_instance_Rsun11(self):
        self.assertIsInstance(totest.Rsun11, (float, int))

    def test_instance_Lsun33(self):
        self.assertIsInstance(totest.Lsun33, (float, int))

    def test_instance_ly(self):
        self.assertIsInstance(totest.ly, (float, int))

    def test_instance_pc(self):
        self.assertIsInstance(totest.pc, (float, int))

    def test_instance_secyer(self):
        self.assertIsInstance(totest.secyer, (float, int))

    def test_instance_dayyer(self):
        self.assertIsInstance(totest.dayyer, (float, int))

    def test_instance_age_of_universe(self):
        self.assertIsInstance(totest.age_of_universe, (float, int))

    def test_instance_Teffsol(self):
        self.assertIsInstance(totest.Teffsol, (float, int))

    def test_instance_loggsol(self):
        self.assertIsInstance(totest.loggsol, (float, int))

    def test_instance_mbolsun(self):
        self.assertIsInstance(totest.mbolsun, (float, int))

    def test_instance_mbolsol(self):
        self.assertIsInstance(totest.mbolsol, (float, int))

    def test_instance_m_earth(self):
        self.assertIsInstance(totest.m_earth, (float, int))

    def test_instance_r_earth(self):
        self.assertIsInstance(totest.r_earth, (float, int))

    def test_instance_au(self):
        self.assertIsInstance(totest.au, (float, int))

    def test_instance_aursun(self):
        self.assertIsInstance(totest.aursun, (float, int))

    def test_instance_m_jupiter(self):
        self.assertIsInstance(totest.m_jupiter, (float, int))

    def test_instance_r_jupiter(self):
        self.assertIsInstance(totest.r_jupiter, (float, int))

    def test_instance_semimajor_axis_jupiter(self):
        self.assertIsInstance(totest.semimajor_axis_jupiter, (float, int))

    def test_instance_km2cm(self):
        self.assertIsInstance(totest.km2cm, (float, int))

    def test_instance_day2sec(self):
        self.assertIsInstance(totest.day2sec, (float, int))

    def test_instance_H_weight(self):
        self.assertIsInstance(totest.H_weight, (float, int))

    def test_instance_He_weight(self):
        self.assertIsInstance(totest.He_weight, (float, int))

    def test_instance_SNcheck_ERR(self):
        self.assertIsInstance(totest.SNcheck_ERR, (float, int))


class TestValues(unittest.TestCase):
    # check that the values fit
    # use delta of last digit times 0.6
    def test_value_pi(self):
        value = 3.1415926535897932384626433832795028841971693993751
        self.assertAlmostEqual(totest.pi, value, delta=6e-50)

    def test_value_a2rad(self):
        value = 1.7453292519943295e-2
        self.assertAlmostEqual(totest.a2rad, value, delta=6e-19)

    def test_value_rad2a(self):
        value = 5.729577951308232e+1
        self.assertAlmostEqual(totest.rad2a, value, delta=6e-15)

    def test_value_standard_cgrav(self):
        value = 6.67428e-8
        self.assertAlmostEqual(totest.standard_cgrav, value, delta=6e-14)

    def test_value_planck_h(self):
        value = 6.62606896e-27
        self.assertAlmostEqual(totest.planck_h, value, delta=6e-36)

    def test_value_hbar(self):
        value = 1.05457163e-27
        self.assertAlmostEqual(totest.hbar, value, delta=6e-36)

    def test_value_qe(self):
        value = 4.80320440e-10
        self.assertAlmostEqual(totest.qe, value, delta=6e-19)

    def test_value_avo(self):
        value = 6.02214129e+23
        self.assertAlmostEqual(totest.avo, value, delta=6e+14)

    def test_value_clight(self):
        value = 2.99792458e+10
        self.assertAlmostEqual(totest.clight, value, delta=6e+1)

    def test_value_kerg(self):
        value = 1.3806504e-16
        self.assertAlmostEqual(totest.kerg, value, delta=6e-24)

    def test_value_boltzm(self):
        value = 1.3806504e-16
        self.assertAlmostEqual(totest.kerg, value, delta=6e-24)

    def test_value_cgas(self):
        value = 8.314471780895016e+7
        self.assertAlmostEqual(totest.cgas, value, delta=6e-1)

    def test_value_kev(self):
        value = 8.617385e-5
        self.assertAlmostEqual(totest.kev, value, delta=6e-12)

    def test_value_amu(self):
        value = 1.6605389210321898e-24
        self.assertAlmostEqual(totest.amu, value, delta=6e-33)

    def test_value_mn(self):
        value = 1.6749286e-24
        self.assertAlmostEqual(totest.mn, value, delta=6e-32)

    def test_value_mp(self):
        value = 1.6726231e-24
        self.assertAlmostEqual(totest.mp, value, delta=6e-32)

    def test_value_me(self):
        value = 9.1093826e-28
        self.assertAlmostEqual(totest.me, value, delta=6e-36)

    def test_value_rbohr(self):
        value = 5.291771539809704e-9
        self.assertAlmostEqual(totest.rbohr, value, delta=6e-18)

    def test_value_fine(self):
        value = 7.297352926107705e-3
        self.assertAlmostEqual(totest.fine, value, delta=6e-11)

    def test_value_hion(self):
        value = 1.3605698140e+1
        self.assertAlmostEqual(totest.hion, value, delta=6e-10)

    def test_value_ev2erg(self):
        value = 1.602176565e-12
        self.assertAlmostEqual(totest.ev2erg, value, delta=6e-22)

    def test_value_inversecm2erg(self):
        value = 1.9864455003959037e-16
        self.assertAlmostEqual(totest.inversecm2erg, value, delta=6e-25)

    def test_value_mev_to_ergs(self):
        value = 1.6021765649999999e-6
        self.assertAlmostEqual(totest.mev_to_ergs, value, delta=6e-16)

    def test_value_mev_amu(self):
        value = 9.648533645956869e+17
        self.assertAlmostEqual(totest.mev_amu, value, delta=6e+8)

    def test_value_Qconv(self):
        value = 9.648533645956868e+17
        self.assertAlmostEqual(totest.Qconv, value, delta=6e+8)

    def test_value_boltz_sigma(self):
        value = 5.670373e-5
        self.assertAlmostEqual(totest.boltz_sigma, value, delta=6e-12)

    def test_value_crad(self):
        value = 7.565731356724124e-15
        self.assertAlmostEqual(totest.crad, value, delta=6e-22)

    def test_value_ssol(self):
        value = 5.670373e-5
        self.assertAlmostEqual(totest.ssol, value, delta=6e-12)

    def test_value_asol(self):
        value = 7.565731356724124e-15
        self.assertAlmostEqual(totest.asol, value, delta=6e-22)

    def test_value_weinlam(self):
        value = 2.897768496231288e-1
        self.assertAlmostEqual(totest.weinlam, value, delta=6e-9)

    def test_value_weinfre(self):
        value = 5.878932774535368e+10
        self.assertAlmostEqual(totest.weinfre, value, delta=6e+2)

    def test_value_rhonuc(self):
        value = 2.342e+14
        self.assertAlmostEqual(totest.rhonuc, value, delta=6e+10)

    def test_value_Zsun(self):
        value = 1.42e-2
        self.assertAlmostEqual(totest.Zsun, value, delta=6e-5)

    def test_value_msol(self):
        value = 1.9892e+33
        self.assertAlmostEqual(totest.msol, value, delta=6e+28)

    def test_value_rsol(self):
        value = 6.9598e+10
        self.assertAlmostEqual(totest.rsol, value, delta=6e+5)

    def test_value_lsol(self):
        value = 3.8418e+33
        self.assertAlmostEqual(totest.lsol, value, delta=6e+28)

    def test_value_agesol(self):
        value = 4.57e+9
        self.assertAlmostEqual(totest.agesol, value, delta=6e+6)

    def test_value_Msun(self):
        value = 1.9892e+33
        self.assertAlmostEqual(totest.Msun, value, delta=6e+28)

    def test_value_Rsun(self):
        value = 6.9598e+10
        self.assertAlmostEqual(totest.Rsun, value, delta=6e+5)

    def test_value_Lsun(self):
        value = 3.8418e+33
        self.assertAlmostEqual(totest.Lsun, value, delta=6e+28)

    def test_value_Msun33(self):
        value = 1.9892
        self.assertAlmostEqual(totest.Msun33, value, delta=6e-5)

    def test_value_Rsun11(self):
        value = 6.9598e-1
        self.assertAlmostEqual(totest.Rsun11, value, delta=6e-6)

    def test_value_Lsun33(self):
        value = 3.8418
        self.assertAlmostEqual(totest.Lsun33, value, delta=6e-5)

    def test_value_ly(self):
        value = 9.460528e+17
        self.assertAlmostEqual(totest.ly, value, delta=6e+10)

    def test_value_pc(self):
        value = 3.0856770322224e+18
        self.assertAlmostEqual(totest.pc, value, delta=6e+11)

    def test_value_secyer(self):
        value = 3.1558149984e+7
        self.assertAlmostEqual(totest.secyer, value, delta=6e-4)

    def test_value_dayyer(self):
        value = 3.6525e+2
        self.assertAlmostEqual(totest.dayyer, value, delta=6e-3)

    def test_value_age_of_universe(self):
        value = 1.38e+10
        self.assertAlmostEqual(totest.age_of_universe, value, delta=6e+7)

    def test_value_Teffsol(self):
        value = 5.7770e+3
        self.assertAlmostEqual(totest.Teffsol, value, delta=6e-2)

    def test_value_loggsol(self):
        value = 4.4378893534131256
        self.assertAlmostEqual(totest.loggsol, value, delta=6e-17)

    def test_value_mbolsun(self):
        value = 4.746
        self.assertAlmostEqual(totest.mbolsun, value, delta=6e-4)

    def test_value_mbolsol(self):
        value = 4.746
        self.assertAlmostEqual(totest.mbolsol, value, delta=6e-4)

    def test_value_m_earth(self):
        value = 5.9764e+27
        self.assertAlmostEqual(totest.m_earth, value, delta=6e+22)

    def test_value_r_earth(self):
        value = 6.37e+8
        self.assertAlmostEqual(totest.r_earth, value, delta=6e+5)

    def test_value_au(self):
        value = 1.495978921e+13
        self.assertAlmostEqual(totest.au, value, delta=6e+3)

    def test_value_aursun(self):
        value = 2.1495e+2
        self.assertAlmostEqual(totest.aursun, value, delta=6e-3)

    def test_value_m_jupiter(self):
        value = 1.8986e+30
        self.assertAlmostEqual(totest.m_jupiter, value, delta=6e+25)

    def test_value_r_jupiter(self):
        value = 6.9911e+9
        self.assertAlmostEqual(totest.r_jupiter, value, delta=6e+4)

    def test_value_semimajor_axis_jupiter(self):
        value = 7.7857e+13
        self.assertAlmostEqual(totest.semimajor_axis_jupiter, value, delta=6e+8)

    def test_value_km2cm(self):
        value = 1.0e5
        self.assertAlmostEqual(totest.km2cm, value, delta=6e-16)

    def test_value_day2sec(self):
        value = 8.64000e+4
        self.assertAlmostEqual(totest.day2sec, value, delta=6e-2)

    def test_value_H_weight(self):
        value = 1.00782503207
        self.assertAlmostEqual(totest.H_weight, value, delta=6e-12)

    def test_value_He_weight(self):
        value = 4.002603254131
        self.assertAlmostEqual(totest.He_weight, value, delta=6e-13)

    def test_value_SNcheck_ERR(self):
        value = 1e-10
        self.assertAlmostEqual(totest.SNcheck_ERR, value, delta=6e-25)


if __name__ == "__main__":
    unittest.main()
