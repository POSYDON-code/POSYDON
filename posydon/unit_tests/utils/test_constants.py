"""Unit tests of posydon/utils/constants.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.utils.constants as totest

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import approx

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = {'H_weight', 'He_weight', 'Lsun', 'Lsun33', 'Msun',\
                    'Msun33', 'Qconv', 'Rsun', 'Rsun11', 'SNcheck_ERR',\
                    'Teffsol', 'Zsun', '__authors__', '__builtins__',\
                    '__cached__', '__doc__', '__file__', '__loader__',\
                    '__name__', '__package__', '__spec__', 'a2rad',\
                    'age_of_universe', 'agesol', 'amu', 'asol', 'au',\
                    'aursun', 'avo', 'boltz_sigma', 'boltzm', 'cgas',\
                    'clight', 'crad', 'day2sec', 'dayyer', 'ev2erg', 'fine',\
                    'hbar', 'hion', 'inversecm2erg', 'kerg', 'kev', 'km2cm',\
                    'loggsol', 'lsol', 'ly', 'm_earth', 'm_jupiter',\
                    'mbolsol', 'mbolsun', 'me', 'mev_amu', 'mev_to_ergs',\
                    'mn', 'mp', 'msol', 'pc', 'pi', 'planck_h', 'qe',\
                    'r_earth', 'r_jupiter', 'rad2a', 'rbohr', 'rhonuc',\
                    'rsol', 'secyer', 'semimajor_axis_jupiter', 'ssol',\
<<<<<<< HEAD
                    'standard_cgrav', 'weinfre', 'weinlam'}
=======
                    'standard_cgrav', 'weinfre', 'weinlam', 'zams_table'}
>>>>>>> eirini_CE_fix
        totest_elements = set(dir(totest))
        missing_in_test = elements - totest_elements
        assert len(missing_in_test) == 0, "There are missing objects in "\
                                          +f"{totest.__name__}: "\
                                          +f"{missing_in_test}. Please "\
                                          +"check, whether they have been "\
                                          +"removed on purpose and update "\
                                          +"this unit test."
        new_in_test = totest_elements - elements
        assert len(new_in_test) == 0, "There are new objects in "\
                                      +f"{totest.__name__}: {new_in_test}. "\
                                      +"Please check, whether they have been "\
                                      +"added on purpose and update this "\
                                      +"unit test."

    def test_instance_pi(self):
        assert isinstance(totest.pi, (float,int)),\
               "pi is of type: " + str(type(totest.pi))

    def test_instance_a2rad(self):
        assert isinstance(totest.a2rad, (float,int)),\
               "a2rad is of type: " + str(type(totest.a2rad))

    def test_instance_rad2a(self):
        assert isinstance(totest.rad2a, (float,int)),\
               "rad2a is of type: " + str(type(totest.rad2a))

    def test_instance_standard_cgrav(self):
        assert isinstance(totest.standard_cgrav, (float,int)),\
               "standard_cgrav is of type: "\
               + str(type(totest.standard_cgrav))

    def test_instance_planck_h(self):
        assert isinstance(totest.planck_h, (float,int)),\
               "planck_h is of type: " + str(type(totest.planck_h))

    def test_instance_hbar(self):
        assert isinstance(totest.hbar, (float,int)),\
               "hbar is of type: " + str(type(totest.hbar))

    def test_instance_qe(self):
        assert isinstance(totest.qe, (float,int)),\
               "qe is of type: " + str(type(totest.qe))

    def test_instance_avo(self):
        assert isinstance(totest.avo, (float,int)),\
               "avo is of type: " + str(type(totest.avo))

    def test_instance_clight(self):
        assert isinstance(totest.clight, (float,int)),\
               "clight is of type: " + str(type(totest.clight))

    def test_instance_kerg(self):
        assert isinstance(totest.kerg, (float,int)),\
               "kerg is of type: " + str(type(totest.kerg))

    def test_instance_boltzm(self):
        assert isinstance(totest.boltzm, (float,int)),\
               "boltzm is of type: " + str(type(totest.boltzm))

    def test_instance_cgas(self):
        assert isinstance(totest.cgas, (float,int)),\
               "cgas is of type: " + str(type(totest.cgas))

    def test_instance_kev(self):
        assert isinstance(totest.kev, (float,int)),\
               "kev is of type: " + str(type(totest.kev))

    def test_instance_amu(self):
        assert isinstance(totest.amu, (float,int)),\
               "amu is of type: " + str(type(totest.amu))

    def test_instance_mn(self):
        assert isinstance(totest.mn, (float,int)),\
               "mn is of type: " + str(type(totest.mn))

    def test_instance_mp(self):
        assert isinstance(totest.mp, (float,int)),\
               "mp is of type: " + str(type(totest.mp))

    def test_instance_me(self):
        assert isinstance(totest.me, (float,int)),\
               "me is of type: " + str(type(totest.me))

    def test_instance_rbohr(self):
        assert isinstance(totest.rbohr, (float,int)),\
               "rbohr is of type: " + str(type(totest.rbohr))

    def test_instance_fine(self):
        assert isinstance(totest.fine, (float,int)),\
               "fine is of type: " + str(type(totest.fine))

    def test_instance_hion(self):
        assert isinstance(totest.hion, (float,int)),\
               "hion is of type: " + str(type(totest.hion))

    def test_instance_ev2erg(self):
        assert isinstance(totest.ev2erg, (float,int)),\
               "ev2erg is of type: " + str(type(totest.ev2erg))

    def test_instance_inversecm2erg(self):
        assert isinstance(totest.inversecm2erg, (float,int)),\
               "inversecm2erg is of type: " + str(type(totest.inversecm2erg))

    def test_instance_mev_to_ergs(self):
        assert isinstance(totest.mev_to_ergs, (float,int)),\
               "mev_to_ergs is of type: " + str(type(totest.mev_to_ergs))

    def test_instance_mev_amu(self):
        assert isinstance(totest.mev_amu, (float,int)),\
               "mev_amu is of type: " + str(type(totest.mev_amu))

    def test_instance_Qconv(self):
        assert isinstance(totest.Qconv, (float,int)),\
               "Qconv is of type: " + str(type(totest.Qconv))

    def test_instance_boltz_sigma(self):
        assert isinstance(totest.boltz_sigma, (float,int)),\
               "boltz_sigma is of type: " + str(type(totest.boltz_sigma))

    def test_instance_crad(self):
        assert isinstance(totest.crad, (float,int)),\
               "crad is of type: " + str(type(totest.crad))

    def test_instance_ssol(self):
        assert isinstance(totest.ssol, (float,int)),\
               "ssol is of type: " + str(type(totest.ssol))

    def test_instance_asol(self):
        assert isinstance(totest.asol, (float,int)),\
               "asol is of type: " + str(type(totest.asol))

    def test_instance_weinlam(self):
        assert isinstance(totest.weinlam, (float,int)),\
               "weinlam is of type: " + str(type(totest.weinlam))

    def test_instance_weinfre(self):
        assert isinstance(totest.weinfre, (float,int)),\
               "weinfre is of type: " + str(type(totest.weinfre))

    def test_instance_rhonuc(self):
        assert isinstance(totest.rhonuc, (float,int)),\
               "rhonuc is of type: " + str(type(totest.rhonuc))

    def test_instance_Zsun(self):
        assert isinstance(totest.Zsun, (float,int)),\
               "Zsun is of type: " + str(type(totest.Zsun))

    def test_instance_msol(self):
        assert isinstance(totest.msol, (float,int)),\
               "msol is of type: " + str(type(totest.msol))

    def test_instance_rsol(self):
        assert isinstance(totest.rsol, (float,int)),\
               "rsol is of type: " + str(type(totest.rsol))

    def test_instance_lsol(self):
        assert isinstance(totest.lsol, (float,int)),\
               "lsol is of type: " + str(type(totest.lsol))

    def test_instance_agesol(self):
        assert isinstance(totest.agesol, (float,int)),\
               "agesol is of type: " + str(type(totest.agesol))

    def test_instance_Msun(self):
        assert isinstance(totest.Msun, (float,int)),\
               "Msun is of type: " + str(type(totest.Msun))

    def test_instance_Rsun(self):
        assert isinstance(totest.Rsun, (float,int)),\
               "Rsun is of type: " + str(type(totest.Rsun))

    def test_instance_Lsun(self):
        assert isinstance(totest.Lsun, (float,int)),\
               "Lsun is of type: " + str(type(totest.Lsun))

    def test_instance_Msun33(self):
        assert isinstance(totest.Msun33, (float,int)),\
               "Msun33 is of type: " + str(type(totest.Msun33))

    def test_instance_Rsun11(self):
        assert isinstance(totest.Rsun11, (float,int)),\
               "Rsun11 is of type: " + str(type(totest.Rsun11))

    def test_instance_Lsun33(self):
        assert isinstance(totest.Lsun33, (float,int)),\
               "Lsun33 is of type: " + str(type(totest.Lsun33))

    def test_instance_ly(self):
        assert isinstance(totest.ly, (float,int)),\
               "ly is of type: " + str(type(totest.ly))

    def test_instance_pc(self):
        assert isinstance(totest.pc, (float,int)),\
               "pc is of type: " + str(type(totest.pc))

    def test_instance_secyer(self):
        assert isinstance(totest.secyer, (float,int)),\
               "secyer is of type: " + str(type(totest.secyer))

    def test_instance_dayyer(self):
        assert isinstance(totest.dayyer, (float,int)),\
               "dayyer is of type: " + str(type(totest.dayyer))

    def test_instance_age_of_universe(self):
        assert isinstance(totest.age_of_universe, (float,int)),\
               "age_of_universe is of type: "\
               + str(type(totest.age_of_universe))

    def test_instance_Teffsol(self):
        assert isinstance(totest.Teffsol, (float,int)),\
               "Teffsol is of type: " + str(type(totest.Teffsol))

    def test_instance_loggsol(self):
        assert isinstance(totest.loggsol, (float,int)),\
               "loggsol is of type: " + str(type(totest.loggsol))

    def test_instance_mbolsun(self):
        assert isinstance(totest.mbolsun, (float,int)),\
               "mbolsun is of type: " + str(type(totest.mbolsun))

    def test_instance_mbolsol(self):
        assert isinstance(totest.mbolsol, (float,int)),\
               "mbolsol is of type: " + str(type(totest.mbolsol))

    def test_instance_m_earth(self):
        assert isinstance(totest.m_earth, (float,int)),\
               "m_earth is of type: " + str(type(totest.m_earth))

    def test_instance_r_earth(self):
        assert isinstance(totest.r_earth, (float,int)),\
               "r_earth is of type: " + str(type(totest.r_earth))

    def test_instance_au(self):
        assert isinstance(totest.au, (float,int)),\
               "au is of type: " + str(type(totest.au))

    def test_instance_aursun(self):
        assert isinstance(totest.aursun, (float,int)),\
               "aursun is of type: " + str(type(totest.aursun))

    def test_instance_m_jupiter(self):
        assert isinstance(totest.m_jupiter, (float,int)),\
               "m_jupiter is of type: " + str(type(totest.m_jupiter))

    def test_instance_r_jupiter(self):
        assert isinstance(totest.r_jupiter, (float,int)),\
               "r_jupiter is of type: " + str(type(totest.r_jupiter))

    def test_instance_semimajor_axis_jupiter(self):
        assert isinstance(totest.semimajor_axis_jupiter, (float,int)),\
               "semimajor_axis_jupiter is of type: "\
               + str(type(totest.semimajor_axis_jupiter))

    def test_instance_km2cm(self):
        assert isinstance(totest.km2cm, (float,int)),\
               "km2cm is of type: " + str(type(totest.km2cm))

    def test_instance_day2sec(self):
        assert isinstance(totest.day2sec, (float,int)),\
               "day2sec is of type: " + str(type(totest.day2sec))

    def test_instance_H_weight(self):
        assert isinstance(totest.H_weight, (float,int)),\
               "H_weight is of type: " + str(type(totest.H_weight))

    def test_instance_He_weight(self):
        assert isinstance(totest.He_weight, (float,int)),\
               "He_weight is of type: " + str(type(totest.He_weight))

    def test_instance_SNcheck_ERR(self):
        assert isinstance(totest.SNcheck_ERR, (float,int)),\
               "SNcheck_ERR is of type: " + str(type(totest.SNcheck_ERR))


class TestValues:
    # check that the values fit
    # use delta of last digit times 0.6
    def test_value_pi(self):
        assert 3.1415926535897932384626433832795028841971693993751 ==\
               approx(totest.pi, abs=6e-50)

    def test_value_a2rad(self):
        assert 1.7453292519943295e-2 == approx(totest.a2rad, abs=6e-19)

    def test_value_rad2a(self):
        assert 5.729577951308232e+1 == approx(totest.rad2a, abs=6e-15)

    def test_value_standard_cgrav(self):
        assert 6.67428e-8 == approx(totest.standard_cgrav, abs=6e-14)

    def test_value_planck_h(self):
        assert 6.62606896e-27 == approx(totest.planck_h, abs=6e-36)

    def test_value_hbar(self):
        assert 1.05457163e-27 == approx(totest.hbar, abs=6e-36)

    def test_value_qe(self):
        assert 4.80320440e-10 == approx(totest.qe, abs=6e-19)

    def test_value_avo(self):
        assert 6.02214129e+23 == approx(totest.avo, abs=6e+14)

    def test_value_clight(self):
        assert 2.99792458e+10 == approx(totest.clight, abs=6e+1)

    def test_value_kerg(self):
        assert 1.3806504e-16 == approx(totest.kerg, abs=6e-24)

    def test_value_boltzm(self):
        assert 1.3806504e-16 == approx(totest.kerg, abs=6e-24)

    def test_value_cgas(self):
        assert 8.314471780895016e+7 == approx(totest.cgas, abs=6e-1)

    def test_value_kev(self):
        assert 8.617385e-5 == approx(totest.kev, abs=6e-12)

    def test_value_amu(self):
        assert 1.6605389210321898e-24 == approx(totest.amu, abs=6e-33)

    def test_value_mn(self):
        assert 1.6749286e-24 == approx(totest.mn, abs=6e-32)

    def test_value_mp(self):
        assert 1.6726231e-24 == approx(totest.mp, abs=6e-32)

    def test_value_me(self):
        assert 9.1093826e-28 == approx(totest.me, abs=6e-36)

    def test_value_rbohr(self):
        assert 5.291771539809704e-9 == approx(totest.rbohr, abs=6e-18)

    def test_value_fine(self):
        assert 7.297352926107705e-3 == approx(totest.fine, abs=6e-11)

    def test_value_hion(self):
        assert 1.3605698140e+1 == approx(totest.hion, abs=6e-10)

    def test_value_ev2erg(self):
        assert 1.602176565e-12 == approx(totest.ev2erg, abs=6e-22)

    def test_value_inversecm2erg(self):
        assert 1.9864455003959037e-16 ==\
               approx(totest.inversecm2erg, abs=6e-25)

    def test_value_mev_to_ergs(self):
        assert 1.6021765649999999e-6 == approx(totest.mev_to_ergs, abs=6e-16)

    def test_value_mev_amu(self):
        assert 9.648533645956869e+17 == approx(totest.mev_amu, abs=6e+8)

    def test_value_Qconv(self):
        assert 9.648533645956868e+17 == approx(totest.Qconv, abs=6e+8)

    def test_value_boltz_sigma(self):
        assert 5.670373e-5 == approx(totest.boltz_sigma, abs=6e-12)

    def test_value_crad(self):
        assert 7.565731356724124e-15 == approx(totest.crad, abs=6e-22)

    def test_value_ssol(self):
        assert 5.670373e-5 == approx(totest.ssol, abs=6e-12)

    def test_value_asol(self):
        assert 7.565731356724124e-15 == approx(totest.asol, abs=6e-22)

    def test_value_weinlam(self):
        assert 2.897768496231288e-1 == approx(totest.weinlam, abs=6e-9)

    def test_value_weinfre(self):
        assert 5.878932774535368e+10 == approx(totest.weinfre, abs=6e+2)

    def test_value_rhonuc(self):
        assert 2.342e+14 == approx(totest.rhonuc, abs=6e+10)

    def test_value_Zsun(self):
        assert 1.42e-2 == approx(totest.Zsun, abs=6e-5)

    def test_value_msol(self):
        assert 1.9892e+33 == approx(totest.msol, abs=6e+28)

    def test_value_rsol(self):
        assert 6.9598e+10 == approx(totest.rsol, abs=6e+5)

    def test_value_lsol(self):
        assert 3.8418e+33 == approx(totest.lsol, abs=6e+28)

    def test_value_agesol(self):
        assert 4.57e+9 == approx(totest.agesol, abs=6e+6)

    def test_value_Msun(self):
        assert 1.9892e+33 == approx(totest.Msun, abs=6e+28)

    def test_value_Rsun(self):
        assert 6.9598e+10 == approx(totest.Rsun, abs=6e+5)

    def test_value_Lsun(self):
        assert 3.8418e+33 == approx(totest.Lsun, abs=6e+28)

    def test_value_Msun33(self):
        assert 1.9892 == approx(totest.Msun33, abs=6e-5)

    def test_value_Rsun11(self):
        assert 6.9598e-1 == approx(totest.Rsun11, abs=6e-6)

    def test_value_Lsun33(self):
        assert 3.8418 == approx(totest.Lsun33, abs=6e-5)

    def test_value_ly(self):
        assert 9.460528e+17 == approx(totest.ly, abs=6e+10)

    def test_value_pc(self):
        assert 3.0856770322224e+18 == approx(totest.pc, abs=6e+11)

    def test_value_secyer(self):
        assert 3.1558149984e+7 == approx(totest.secyer, abs=6e-4)

    def test_value_dayyer(self):
        assert 3.6525e+2 == approx(totest.dayyer, abs=6e-3)

    def test_value_age_of_universe(self):
        assert 1.38e+10 == approx(totest.age_of_universe, abs=6e+7)

    def test_value_Teffsol(self):
        assert 5.7770e+3 == approx(totest.Teffsol, abs=6e-2)

    def test_value_loggsol(self):
        assert 4.4378893534131256 == approx(totest.loggsol, abs=6e-17)

    def test_value_mbolsun(self):
        assert 4.746 == approx(totest.mbolsun, abs=6e-4)

    def test_value_mbolsol(self):
        assert 4.746 == approx(totest.mbolsol, abs=6e-4)

    def test_value_m_earth(self):
        assert 5.9764e+27 == approx(totest.m_earth, abs=6e+22)

    def test_value_r_earth(self):
        assert 6.37e+8 == approx(totest.r_earth, abs=6e+5)

    def test_value_au(self):
        assert 1.495978921e+13 == approx(totest.au, abs=6e+3)

    def test_value_aursun(self):
        assert 2.1495e+2 == approx(totest.aursun, abs=6e-3)

    def test_value_m_jupiter(self):
        assert 1.8986e+30 == approx(totest.m_jupiter, abs=6e+25)

    def test_value_r_jupiter(self):
        assert 6.9911e+9 == approx(totest.r_jupiter, abs=6e+4)

    def test_value_semimajor_axis_jupiter(self):
        assert 7.7857e+13 == approx(totest.semimajor_axis_jupiter, abs=6e+8)

    def test_value_km2cm(self):
        assert 1.0e5 == approx(totest.km2cm, abs=6e-16)

    def test_value_day2sec(self):
        assert 8.64000e+4 == approx(totest.day2sec, abs=6e-2)

    def test_value_H_weight(self):
        assert 1.00782503207 == approx(totest.H_weight, abs=6e-12)

    def test_value_He_weight(self):
        assert 4.002603254131 == approx(totest.He_weight, abs=6e-13)

    def test_value_SNcheck_ERR(self):
        assert 1e-10 == approx(totest.SNcheck_ERR, abs=6e-31)
