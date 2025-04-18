import unittest
import os
import numpy as np
from posydon.config import PATH_TO_POSYDON
from posydon.utils import common_functions as cf
from posydon.binary_evol.binarystar import BinaryStar
from posydon.binary_evol.singlestar import SingleStar
from posydon.binary_evol.CE.step_CEE import StepCEE

# spaces are read '\\ ' instead of ' '
PATH_TO_DATA = os.path.join(
    PATH_TO_POSYDON, "posydon/tests/data/POSYDON-UNIT-TESTS/binary_evol/CE/")


class TestCommonEnvelope(unittest.TestCase):
    def test_common_envelope_1(self):
        kwargs = {'prescription': 'alpha-lambda',
                  "common_envelope_option_for_lambda" : 'default_lambda'}

        CEE = StepCEE(verbose=False, **kwargs)

        # simple binary system which will experience CEE with default_lambda
        # option on.
        # no profiles needed for this
        PROPERTIES_STAR1 = {
            'mass': 10.0,
            'log_R': np.log10(1000.0),
            'he_core_mass': 3.0,
            'he_core_radius': 0.5,
            'state': 'H-rich_Shell_H_burning',
            'metallicity' : 0.0142,
            "log_Lnuc": -1e6, # arbitrary
            "log_LHe": -1e7, # arbitrary
            "center_he4" : 0.0,
            "center_h1" : 1.0,
            "center_c12" : 0.01,
        }
        PROPERTIES_STAR2 = {
            'mass': 2.0,
            'log_R': np.log10(2.0),
            'he_core_mass': 0.0,
            'he_core_radius': 0.0,
            'state': 'H-rich_Core_H_burning',
            'metallicity' : 0.0142,
            "log_Lnuc": -1e6, # arbitrary
            "log_LHe": -1e7, # arbitrary
            "center_he4" : 0.49,
            "center_h1" : 0.49,
            "center_c12" : 0.01,
        }
        giantstar = SingleStar(**PROPERTIES_STAR1)
        compstar = SingleStar(**PROPERTIES_STAR2)

        orbital_separation_for_RLOF = 10**giantstar.log_R / cf.roche_lobe_radius(
            giantstar.mass, compstar.mass, a_orb=1)
        orbital_period_for_RLOF = cf.orbital_period_from_separation(
            orbital_separation_for_RLOF, giantstar.mass, compstar.mass)
        PROPERTIES_BINARY = {
            "binary_state": "RLO1",
            "event": "oCE1",
            "orbital_period": orbital_period_for_RLOF
        }
        binary = BinaryStar(star_1=giantstar,
                            star_2=compstar,
                            **PROPERTIES_BINARY)

        CEE(binary)
        #self.assertTrue(binary.event == 'redirect', "CEE test 1 failed")
        self.assertTrue(
            abs(binary.orbital_period - 5.056621408721529) <
            1.0, "CEE test 1 failed")
        self.assertTrue("stripped_He" in binary.star_1.state, "CEE test 1 failed")

    def test_common_envelope_2(self):
        kwargs = {'prescription': 'alpha-lambda',
        "common_envelope_option_for_lambda" : 'lambda_from_profile_gravitational'}

        CEE = StepCEE(verbose=False, **kwargs)

        # testing with loading a profile of the donor at the moment of CEE to
        # calculate the lamda CEE
        profile_donor_name = os.path.join(PATH_TO_DATA,
                                          'simple_giant_profile_for_CEE.npy')
        profile_donor = np.load(profile_donor_name)
        PROPERTIES_STAR1_withprofile = {
            'mass': 22.77,
            'log_R': np.log10(1319.0),
            'he_core_mass': 11.0,
            'he_core_radius': 0.6,
            'state': 'H-rich_Shell_H_burning',
            'metallicity' : 0.0142,
            'profile': profile_donor,
            "log_Lnuc": -1e6, # arbitrary
            "log_LHe": -1e7, # arbitrary
            "center_he4" : 0.0,
            "center_h1" : 1.0,
            "center_c12" : 0.01,
        }
        giantstar_withprofile = SingleStar(**PROPERTIES_STAR1_withprofile)
        PROPERTIES_STAR2 = {
            'mass': 2.0,
            'log_R': np.log10(2.0),
            'he_core_mass': 0.0,
            'he_core_radius': 0.0,
            'state': 'H-rich_Core_H_burning',
            'metallicity' : 0.0142,
            "log_Lnuc": -1e6, # arbitrary
            "log_LHe": -1e7, # arbitrary
            "center_he4" : 0.5,
            "center_h1" : 0.5,
            "center_c12" : 0.01,
        }
        compstar = SingleStar(**PROPERTIES_STAR2)

        orbital_separation_for_RLOF = 10**giantstar_withprofile.log_R / cf.roche_lobe_radius(
            giantstar_withprofile.mass, compstar.mass, a_orb=1)
        orbital_period_for_RLOF = cf.orbital_period_from_separation(
            orbital_separation_for_RLOF, giantstar_withprofile.mass,
            compstar.mass)

        PROPERTIES_BINARY_withprofile = {
            "binary_state": "RLO1",
            "event": "oCE1",
            "orbital_period": orbital_period_for_RLOF
        }
        binary_withprofile = BinaryStar(giantstar_withprofile, compstar,
                                        **PROPERTIES_BINARY_withprofile)
        # options: 'default_lambda', 'lambda_from_profile_gravitational',
        # 'lambda_from_profile_gravitational_plus_internal',
        # 'lambda_from_profile_gravitational_plus_internal_minus_recombination'
        #binary_withprofile.properties.common_envelope_option_for_lambda = "lambda_from_profile_gravitational"
        CEE(binary_withprofile)
        #print(binary_withprofile.event)
        self.assertTrue(binary_withprofile.state == "merged",
                        "CEE test 2 failed")

    def test_common_envelope_3(self):
        kwargs = {'prescription': 'alpha-lambda',
        "common_envelope_option_for_lambda" : 'lambda_from_profile_gravitational_plus_internal_minus_recombination'}

        CEE = StepCEE(verbose=False, **kwargs)

        # testing with loading a profile of the donor at the moment of CEE to
        # calculate the lamda CEE, taking into account also the internal energy
        # - recombination energy
        profile_donor_name = os.path.join(
            PATH_TO_DATA,
            'giant_profile_for_CEE_with_recombinationenergy_calculation.npy')
        #profile_donor = np.load(profile_donor_name, mmap_mode = "r")
        profile_donor = np.load(profile_donor_name)
        PROPERTIES_STAR1_withprofile = {
            'mass': 22.77,
            'log_R': np.log10(1319.0),
            'he_core_mass': 11.0,
            'he_core_radius': 0.6,
            'state': 'H-rich_Shell_H_burning',
            'metallicity' : 0.0142,
            'profile': profile_donor,
            "log_Lnuc": -1e6, # arbitrary
            "log_LHe": -1e7, # arbitrary
            "center_he4" : 0.0,
            "center_h1" : 1.0,
            "center_c12" : 0.01,
        }
        giantstar_withprofile = SingleStar(**PROPERTIES_STAR1_withprofile)
        PROPERTIES_STAR2 = {
            'mass': 2.0,
            'log_R': np.log10(2.0),
            'he_core_mass': 0.0,
            'he_core_radius': 0.0,
            'state': 'H-rich_Core_H_burning',
            'metallicity' : 0.0142,
            "log_Lnuc": -1e6, # arbitrary
            "log_LHe": -1e7, # arbitrary
            "center_he4" : 0.49,
            "center_h1" : 0.49,
            "center_c12" : 0.01,
        }
        compstar = SingleStar(**PROPERTIES_STAR2)

        orbital_separation_for_RLOF = 10**giantstar_withprofile.log_R / cf.roche_lobe_radius(
            giantstar_withprofile.mass, compstar.mass, a_orb=1)
        orbital_period_for_RLOF = cf.orbital_period_from_separation(
            orbital_separation_for_RLOF, giantstar_withprofile.mass,
            compstar.mass)

        PROPERTIES_BINARY_withprofile = {
            "binary_state": "RLO1",
            "event": "oCE1",
            "orbital_period": orbital_period_for_RLOF
        }
        binary_withprofile = BinaryStar(giantstar_withprofile, compstar,
                                        **PROPERTIES_BINARY_withprofile)
        # options: 'default_lambda', 'lambda_from_profile_gravitational',
        # 'lambda_from_profile_gravitational_plus_internal',
        # 'lambda_from_profile_gravitational_plus_internal_minus_recombination'
        #binary_withprofile.properties.common_envelope_option_for_lambda = "lambda_from_profile_gravitational_plus_internal_minus_recombination"
        #print(binary_withprofile.properties.common_envelope_option_for_lambda)
        CEE(binary_withprofile)
        #print(binary_withprofile.event)
        self.assertTrue(binary_withprofile.state == "merged",
                        "CEE test 3 failed")

    def test_common_envelope_4(self):
        kwargs = {'prescription': 'alpha-lambda',
        "common_envelope_option_for_lambda" : 'lambda_from_profile_gravitational'}

        CEE = StepCEE(verbose=False, **kwargs)

        profile_donor_name = os.path.join(PATH_TO_DATA,
                                          'simple_giant_profile_for_CEE.npy')
        profile_donor = np.load(profile_donor_name)
        PROPERTIES_STAR1_withprofile = {
            'mass': 22.77,
            'log_R': np.log10(1319.0),
            'he_core_mass': 11.0,
            'he_core_radius': 0.6,
            'state': 'H-rich_Shell_H_burning',
            'metallicity' : 0.0142,
            'profile': profile_donor,
            "log_Lnuc": -1e6, # arbitrary
            "log_LHe": -1e7, # arbitrary
            "center_he4" : 0.0,
            "center_h1" : 1.0,
            "center_c12" : 0.01,
        }
        giantstar_withprofile = SingleStar(**PROPERTIES_STAR1_withprofile)
        PROPERTIES_STAR2 = {
            'mass': 10.0,
            'log_R': np.log10(7.0),
            'he_core_mass': 0.0,
            'he_core_radius': 0.0,
            'state': 'H-rich_Core_H_burning',
            'metallicity' : 0.0142,
            "log_Lnuc": -1e6, # arbitrary
            "log_LHe": -1e7, # arbitrary
            "center_he4" : 0.49,
            "center_h1" : 0.49,
            "center_c12" : 0.01,
        }
        compstar = SingleStar(**PROPERTIES_STAR2)

        orbital_separation_for_RLOF = 10**giantstar_withprofile.log_R / cf.roche_lobe_radius(
            giantstar_withprofile.mass, compstar.mass, a_orb=1)
        orbital_period_for_RLOF = cf.orbital_period_from_separation(
            orbital_separation_for_RLOF, giantstar_withprofile.mass,
            compstar.mass)

        PROPERTIES_BINARY_withprofile = {
            "binary_state": "RLO1",
            "event": "oCE1",
            "orbital_period": orbital_period_for_RLOF
        }
        binary_withprofile = BinaryStar(giantstar_withprofile, compstar,
                                        **PROPERTIES_BINARY_withprofile)
        # options: 'default_lambda', 'lambda_from_profile_gravitational',
        # 'lambda_from_profile_gravitational_plus_internal',
        # 'lambda_from_profile_gravitational_plus_internal_minus_recombination'
        #binary_withprofile.properties.common_envelope_option_for_lambda = "lambda_from_profile_gravitational"
        #print(binary_withprofile.properties.common_envelope_option_for_lambda)
        CEE(binary_withprofile)
        #print(binary_withprofile.event)
        self.assertTrue(binary_withprofile.state == "merged",
                        "CEE test 4 failed")

    def test_common_envelope_5(self):
        kwargs = {'prescription': 'alpha-lambda',
        "common_envelope_option_for_lambda" : 'lambda_from_profile_gravitational_plus_internal_minus_recombination'}

        CEE = StepCEE(verbose=False, **kwargs)
        # testing with loading a profile of the donor at the moment of CEE to
        # calculate the lamda CEE, taking into account also the internal
        # energy - recombination energy
        profile_donor_name = os.path.join(
            PATH_TO_DATA,
            'giant_profile_for_CEE_with_recombinationenergy_calculation.npy')
        profile_donor = np.load(profile_donor_name)
        PROPERTIES_STAR1_withprofile = {
            'mass': 22.77,
            'log_R': np.log10(1319.0),
            'he_core_mass': 11.0,
            'he_core_radius': 0.6,
            'state': 'H-rich_Shell_H_burning',
            'metallicity' : 0.0142,
            'profile': profile_donor,
            "log_Lnuc": -1e6, # arbitrary
            "log_LHe": -1e7, # arbitrary
            "center_he4" : 0.0,
            "center_h1" : 1.0,
            "center_c12" : 0.01,
        }
        giantstar_withprofile = SingleStar(**PROPERTIES_STAR1_withprofile)
        PROPERTIES_STAR2 = {
            'mass': 10.0,
            'log_R': np.log10(7.0),
            'he_core_mass': 0.0,
            'he_core_radius': 0.0,
            'state': 'H-rich_Core_H_burning',
            'metallicity' : 0.0142,
            "log_Lnuc": -1e6, # arbitrary
            "log_LHe": -1e7, # arbitrary
            "center_he4" : 0.49,
            "center_h1" : 0.49,
            "center_c12" : 0.01,
        }
        compstar = SingleStar(**PROPERTIES_STAR2)

        orbital_separation_for_RLOF = 10**giantstar_withprofile.log_R / cf.roche_lobe_radius(
            giantstar_withprofile.mass, compstar.mass, a_orb=1)
        orbital_period_for_RLOF = cf.orbital_period_from_separation(
            orbital_separation_for_RLOF, giantstar_withprofile.mass,
            compstar.mass)

        PROPERTIES_BINARY_withprofile = {
            "binary_state": "RLO1",
            "event": "oCE1",
            "orbital_period": orbital_period_for_RLOF
        }
        binary_withprofile = BinaryStar(giantstar_withprofile, compstar,
                                        **PROPERTIES_BINARY_withprofile)
        # options: 'default_lambda', 'lambda_from_profile_gravitational',
        # 'lambda_from_profile_gravitational_plus_internal',
        # 'lambda_from_profile_gravitational_plus_internal_minus_recombination'
        #binary_withprofile.properties.common_envelope_option_for_lambda = "lambda_from_profile_gravitational_plus_internal_minus_recombination"
        #print(binary_withprofile.properties.common_envelope_option_for_lambda)
        CEE(binary_withprofile)
        #print(binary_withprofile.event)
        self.assertTrue(binary_withprofile.state == "merged",
                        "CEE test 5 failed")

    def test_common_envelope_6(self):
        kwargs = {'prescription': 'alpha-lambda',
        "common_envelope_option_for_lambda" : 'lambda_from_profile_gravitational'}

        CEE = StepCEE(verbose=False, **kwargs)

        profile_donor_name = os.path.join(PATH_TO_DATA,
                                          'simple_giant_profile_for_CEE.npy')
        #profile_donor =  np.genfromtxt(profile_donor_name, skip_header=5, names=True, dtype=None)
        profile_donor = np.load(profile_donor_name)
        PROPERTIES_STAR1_withprofile = {
            'mass': 22.77,
            'log_R': np.log10(1319.0),
            'he_core_mass': 11.0,
            'he_core_radius': 0.6,
            'state': 'H-rich_Shell_H_burning',
            'metallicity' : 0.0142,
            'profile': profile_donor,
            "log_Lnuc": -1e6, # arbitrary
            "log_LHe": -1e7, # arbitrary
            "center_he4" : 0.0,
            "center_h1" : 1.0,
            "center_c12" : 0.01,
        }
        giantstar_withprofile = SingleStar(**PROPERTIES_STAR1_withprofile)
        PROPERTIES_STAR2 = {
            'mass': 10.,
            'log_R': np.log10(0.0001),
            'he_core_mass': 0.0,
            'he_core_radius': 0.0,
            'state': 'BH'
        }
        compstar = SingleStar(**PROPERTIES_STAR2)

        orbital_separation_for_RLOF = 10**giantstar_withprofile.log_R / cf.roche_lobe_radius(
            giantstar_withprofile.mass, compstar.mass, a_orb=1)
        orbital_period_for_RLOF = cf.orbital_period_from_separation(
            orbital_separation_for_RLOF, giantstar_withprofile.mass,
            compstar.mass)

        PROPERTIES_BINARY_withprofile = {
            "binary_state": "RLO1",
            "event": "oCE1",
            "orbital_period": orbital_period_for_RLOF
        }
        binary_withprofile = BinaryStar(giantstar_withprofile, compstar,
                                        **PROPERTIES_BINARY_withprofile)
        #binary_withprofile.properties.common_envelope_option_for_lambda = "lambda_from_profile_gravitational"  # options: 'default_lambda', 'lambda_from_profile_gravitational', 'lambda_from_profile_gravitational_plus_internal', 'lambda_from_profile_gravitational_plus_internal_minus_recombination'
        #print(binary_withprofile.properties.common_envelope_option_for_lambda)
        CEE(binary_withprofile)
        #print("new state of the star that triggered CEE = ",giantstar_withprofile.state)
        #print("new mass of the star that triggered CEE = ",giantstar_withprofile.mass)
        print(binary_withprofile.event)
        #self.assertTrue(binary_withprofile.event == 'redirect',
        #                "CEE test 6 failed")
        self.assertTrue((10**giantstar_withprofile.log_R - cf.roche_lobe_radius(
            giantstar_withprofile.mass, compstar.mass,
            a_orb=cf.orbital_separation_from_period(
                binary_withprofile.orbital_period, giantstar_withprofile.mass,
                compstar.mass))), "CEE test 6 failed")

        self.assertTrue(
            (abs(binary_withprofile.orbital_period - 0.12123905531545925) <
             1.0),
            "CEE test 6 failed")

    def test_common_envelope_7(self):
        kwargs = {'prescription': 'alpha-lambda',
        "common_envelope_option_for_lambda" : 'lambda_from_profile_gravitational_plus_internal_minus_recombination'}

        CEE = StepCEE(verbose=False, **kwargs)

        # testing with loading a profile of the donor at the moment of CEE to
        # calculate the lamda CEE, taking into account also the internal
        # energy - recombination energy
        profile_donor_name = os.path.join(
            PATH_TO_DATA,
            'giant_profile_for_CEE_with_recombinationenergy_calculation.npy')
        profile_donor = np.load(profile_donor_name)
        PROPERTIES_STAR1_withprofile = {
            'mass': 22.77,
            'log_R': np.log10(1319.0),
            'he_core_mass': 11.0,
            'he_core_radius': 0.6,
            'state': 'H-rich_Shell_H_burning',
            'metallicity' : 0.0142,
            'profile': profile_donor,
            "log_Lnuc": -1e6, # arbitrary
            "log_LHe": -1e7, # arbitrary
            "center_he4" : 0.0,
            "center_h1" : 1.0,
            "center_c12" : 0.01,
        }
        giantstar_withprofile = SingleStar(**PROPERTIES_STAR1_withprofile)
        PROPERTIES_STAR2 = {
            'mass': 10.,
            'log_R': np.log10(0.0001),
            'he_core_mass': 0.0,
            'he_core_radius': 0.0,
            'state': 'BH'
        }
        compstar = SingleStar(**PROPERTIES_STAR2)

        orbital_separation_for_RLOF = 10**giantstar_withprofile.log_R / cf.roche_lobe_radius(
            giantstar_withprofile.mass, compstar.mass, a_orb=1)
        orbital_period_for_RLOF = cf.orbital_period_from_separation(
            orbital_separation_for_RLOF, giantstar_withprofile.mass,
            compstar.mass)

        PROPERTIES_BINARY_withprofile = {
            "binary_state": "RLO1",
            "event": "oCE1",
            "orbital_period": orbital_period_for_RLOF
        }
        binary_withprofile = BinaryStar(giantstar_withprofile, compstar,
                                        **PROPERTIES_BINARY_withprofile)
        # options: 'default_lambda', 'lambda_from_profile_gravitational',
        # 'lambda_from_profile_gravitational_plus_internal',
        # 'lambda_from_profile_gravitational_plus_internal_minus_recombination'
        #binary_withprofile.properties.common_envelope_option_for_lambda = "lambda_from_profile_gravitational_plus_internal_minus_recombination"
        #print(binary_withprofile.properties.common_envelope_option_for_lambda)
        CEE(binary_withprofile)
        #print(binary_withprofile.event)
        #self.assertTrue(binary_withprofile.event == 'redirect',
        #                "CEE test 7 failed")
        self.assertTrue((10**giantstar_withprofile.log_R - cf.roche_lobe_radius(
            giantstar_withprofile.mass, compstar.mass,
            a_orb=cf.orbital_separation_from_period(
                binary_withprofile.orbital_period, giantstar_withprofile.mass,
                compstar.mass))), "CEE test 7 failed")
        self.assertTrue(
            (abs(binary_withprofile.orbital_period - 0.3287114957064215) <
             1.0),
            "CEE test 7 failed")

    def test_common_envelope_8(self):
        kwargs = {'prescription': 'alpha-lambda',
        "common_envelope_option_for_lambda" : 'lambda_from_profile_gravitational_plus_internal_minus_recombination'}

        CEE = StepCEE(verbose=False, **kwargs)

        profile_donor_name = os.path.join(
            PATH_TO_DATA,
            'giant_profile_for_CEE_with_recombinationenergy_calculation.npy')
        profile_donor = np.load(profile_donor_name)
        PROPERTIES_STAR1_withprofile = {
            'mass': 22.77,
            'log_R': np.log10(1319.0),
            'he_core_mass': 11.0,
            'he_core_radius': 0.6,
            'state': 'H-rich_Shell_H_burning',
            'metallicity' : 0.0142,
            'profile': profile_donor,
            "log_Lnuc": -1e6, # arbitrary
            "log_LHe": -1e7, # arbitrary
            "center_he4" : 0.0,
            "center_h1" : 1.0,
            "center_c12" : 0.01,
        }
        giantstar_withprofile = SingleStar(**PROPERTIES_STAR1_withprofile)
        PROPERTIES_STAR2 = {
            'mass': 20.,
            'log_R': np.log10(0.0001),
            'he_core_mass': 0.0,
            'he_core_radius': 0.0,
            'state': 'BH'
        }  #radius = 10km
        compstar = SingleStar(**PROPERTIES_STAR2)

        orbital_separation_for_RLOF = 10**giantstar_withprofile.log_R / cf.roche_lobe_radius(
            giantstar_withprofile.mass, compstar.mass, a_orb=1)
        orbital_period_for_RLOF = cf.orbital_period_from_separation(
            orbital_separation_for_RLOF, giantstar_withprofile.mass,
            compstar.mass)

        PROPERTIES_BINARY_withprofile = {
            "binary_state": "RLO1",
            "event": "oCE1",
            "orbital_period": orbital_period_for_RLOF
        }
        binary_withprofile = BinaryStar(giantstar_withprofile, compstar,
                                        **PROPERTIES_BINARY_withprofile)
        # options: 'default_lambda', 'lambda_from_profile_gravitational',
        # 'lambda_from_profile_gravitational_plus_internal',
        # 'lambda_from_profile_gravitational_plus_internal_minus_recombination'
        #binary_withprofile.properties.common_envelope_option_for_lambda = "lambda_from_profile_gravitational_plus_internal_minus_recombination"
        #print(binary_withprofile.properties.common_envelope_option_for_lambda)
        CEE(binary_withprofile)
        #print(binary_withprofile.event)
        #self.assertTrue(binary_withprofile.event == 'redirect',
        #                "CEE test 8 failed")
        #self.assertTrue(binary_withprofile.star_1.state == "stripped_He_Core_He_burning",
        #                "CEE test 8 failed")
        self.assertTrue("stripped_He" in binary_withprofile.star_1.state,
                        "CEE test 8 failed")
        self.assertTrue(
            (abs(binary_withprofile.orbital_period - 0.7636524660283687) <
             1.0),
            "CEE test 8 failed")

    def test_common_envelope_9(self):
        kwargs = {'prescription': 'alpha-lambda',
        "common_envelope_option_for_lambda" : 'lambda_from_profile_gravitational_plus_internal_minus_recombination'}

        CEE = StepCEE(verbose=False, **kwargs)

        profile_donor_name = os.path.join(PATH_TO_DATA,
                                          'caseB_CEE_profile.npy')
        profile_donor = np.load(profile_donor_name)
        PROPERTIES_STAR1_withprofile = {
            'mass': 28.04,
            'log_R': np.log10(927.0),
            'he_core_mass': 11.0,
            'he_core_radius': 0.6,
            'state': 'H-rich_Shell_H_burning',
            'metallicity' : 0.0142,
            'profile': profile_donor,
            "log_Lnuc": -1e6, # arbitrary
            "log_LHe": -1e7, # arbitrary
            "center_he4" : 0.0,
            "center_h1" : 1.0,
            "center_c12" : 0.01,
        }
        giantstar_withprofile = SingleStar(**PROPERTIES_STAR1_withprofile)
        PROPERTIES_STAR2 = {
            'mass': 20.,
            'log_R': np.log10(0.0001),
            'he_core_mass': 0.0,
            'he_core_radius': 0.0,
            'state': 'BH'
        }  #radius = 10km
        compstar = SingleStar(**PROPERTIES_STAR2)

        orbital_separation_for_RLOF = 10**giantstar_withprofile.log_R / cf.roche_lobe_radius(
            giantstar_withprofile.mass, compstar.mass, a_orb=1)
        orbital_period_for_RLOF = cf.orbital_period_from_separation(
            orbital_separation_for_RLOF, giantstar_withprofile.mass,
            compstar.mass)

        PROPERTIES_BINARY_withprofile = {
            "binary_state": "RLO1",
            "event": "oCE1",
            "orbital_period": orbital_period_for_RLOF
        }
        binary_withprofile = BinaryStar(giantstar_withprofile, compstar,
                                        **PROPERTIES_BINARY_withprofile)
        # options: 'default_lambda', 'lambda_from_profile_gravitational',
        # 'lambda_from_profile_gravitational_plus_internal',
        # 'lambda_from_profile_gravitational_plus_internal_minus_recombination'
        #binary_withprofile.properties.common_envelope_option_for_lambda = "lambda_from_profile_gravitational_plus_internal_minus_recombination"
        #print(binary_withprofile.properties.common_envelope_option_for_lambda)
        CEE(binary_withprofile)
        #print(binary_withprofile.event)
        #self.assertTrue(binary_withprofile.event == 'redirect',
        #                "CEE test 9 failed event")
        self.assertTrue("stripped_He" in binary_withprofile.star_1.state,
                        "CEE test 9 failed state")
        self.assertTrue(
            (abs(binary_withprofile.orbital_period - 0.166535882054919) <
             1.0),
            "CEE test 9 failed tolerance")

    def test_common_envelope_10(self):
        kwargs = {'prescription': 'alpha-lambda',
        "common_envelope_option_for_lambda" : 'lambda_from_profile_gravitational_plus_internal_minus_recombination'}

        CEE = StepCEE(verbose=False, **kwargs)

        profile_donor_name = os.path.join(PATH_TO_DATA,
                                          'caseB_CEE_profile.npy')
        profile_donor = np.load(profile_donor_name)
        PROPERTIES_STAR1_withprofile = {
            'mass': 28.04,
            'log_R': np.log10(927.0),
            'he_core_mass': 11.0,
            'he_core_radius': 0.6,
            'state': 'H-rich_Shell_H_burning',
            'metallicity' : 0.0142,
            'profile': profile_donor,
            "log_Lnuc": -1e6, # arbitrary
            "log_LHe": -1e7, # arbitrary
            "center_he4" : 0.0,
            "center_h1" : 1.0,
            "center_c12" : 0.01,
        }
        giantstar_withprofile = SingleStar(**PROPERTIES_STAR1_withprofile)
        PROPERTIES_STAR2 = {
            'mass': 2.,
            'log_R': np.log10(1.5),
            'he_core_mass': 0.0,
            'he_core_radius': 0.0,
            'state': 'H-rich_Core_H_burning',
            'metallicity' : 0.0142,
            "log_Lnuc": -1e6, # arbitrary
            "log_LHe": -1e7, # arbitrary
            "center_he4" : 0.49,
            "center_h1" : 0.49,
            "center_c12" : 0.01,
        }  #radius = 10km
        compstar = SingleStar(**PROPERTIES_STAR2)

        orbital_separation_for_RLOF = 10**giantstar_withprofile.log_R / cf.roche_lobe_radius(
            giantstar_withprofile.mass, compstar.mass, a_orb=1)
        orbital_period_for_RLOF = cf.orbital_period_from_separation(
            orbital_separation_for_RLOF, giantstar_withprofile.mass,
            compstar.mass)

        PROPERTIES_BINARY_withprofile = {
            "binary_state": "RLO1",
            "event": "oCE1",
            "orbital_period": orbital_period_for_RLOF
        }
        binary_withprofile = BinaryStar(giantstar_withprofile, compstar,
                                        **PROPERTIES_BINARY_withprofile)
        # options: 'default_lambda', 'lambda_from_profile_gravitational',
        # 'lambda_from_profile_gravitational_plus_internal',
        # 'lambda_from_profile_gravitational_plus_internal_minus_recombination'
        #binary_withprofile.properties.common_envelope_option_for_lambda = "lambda_from_profile_gravitational_plus_internal_minus_recombination"
        #print(binary_withprofile.properties.common_envelope_option_for_lambda)
        CEE(binary_withprofile)
        #print(binary_withprofile.event)
        self.assertTrue(binary_withprofile.state == "merged",
                        "CEE test 10 failed")
