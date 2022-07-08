import unittest

from posydon.utils.common_functions import PATH_TO_POSYDON
from posydon.binary_evol.SN.step_SN import StepSN
from posydon.binary_evol.binarystar import BinaryStar
from posydon.binary_evol.singlestar import SingleStar
from posydon.binary_evol import SimulationProperties

import numpy as np
import os
import pandas as pd
from scipy.stats import maxwell
import matplotlib.cm as cm

# github action are not cloning the data submoule, data for unit testing
# are therefore stored to the unit test submodule

path_to_Sukhbold_datasets = os.path.join(
    PATH_TO_POSYDON, "posydon/tests/data/POSYDON-UNIT-TESTS/binary_evol/SN/")

class TestStepSN(unittest.TestCase):
    # TODO
    '''
        """
        Test WD formation
        """
        def test_WD_formation_RAPID(self):
            M_CO = 1.3
            SN = StepSN( **{'mechanism' : 'Fryer+12-rapid',
                            'engine' : None,
                            'PISN' : 'Marchant+19',
                            'ECSN' : 'cosmic',
                            'max_neutrino_mass_loss' : 0.,
                            'kick' : True,
                            'sigma_kick_CCSN' : 265.0,
                            'sigma_kick_ECSN' : 20.0,
                            'max_NS_mass' : 2.5,
                            'verbose' : False} )

            star_prop = {'mass': M_CO / 0.7638113015667961,
                            'co_core_mass': M_CO ,
                            'he_core_mass': M_CO / 0.7638113015667961,
                            'state':'stripped_He_Central_C_depletion',
                            'profile':None ,
                            'spin': 0.0}
            star = SingleStar(**star_prop)

            star_he_core = star.he_core_mass

            M_rembar = SN.compute_m_rembar(star , None)[0]

            self.assertEqual( SN.SN_type , 'WD')
            self.assertEqual( M_rembar , star_he_core)

        def test_WD_formation_DELAYED(self):
            M_CO = 1.3
            SN = StepSN( **{'mechanism' : 'Fryer+12-delayed',
                            'engine' : None,
                            'PISN' : 'Marchant+19',
                            'ECSN' : 'cosmic',
                            'max_neutrino_mass_loss' : 0.,
                            'kick' : True,
                            'sigma_kick_CCSN' : 265.0,
                            'sigma_kick_ECSN' : 20.0,
                            'max_NS_mass' : 2.5,
                            'verbose' : False} )

            star_prop = {'mass': M_CO / 0.7638113015667961,
                            'co_core_mass': M_CO ,
                            'he_core_mass': M_CO / 0.7638113015667961,
                            'state':'stripped_He_Central_C_depletion',
                            'profile':None ,
                            'spin': 0.0}
            star = SingleStar(**star_prop)

            star_he_core = star.he_core_mass

            M_rembar = SN.compute_m_rembar(star , None)[0]

            self.assertEqual( SN.SN_type , 'WD')
            self.assertEqual( M_rembar , star_he_core)

        def test_WD_formation_SUKHBOLDN20(self):
            M_CO = 1.3
            SN = StepSN( **{'mechanism' : 'Sukhbold+16-engine',
                            'engine' : 'N20',
                            'PISN' : 'Marchant+19',
                            'ECSN' : 'cosmic',
                            'max_neutrino_mass_loss' : 0.,
                            'kick' : True,
                            'sigma_kick_CCSN' : 265.0,
                            'sigma_kick_ECSN' : 20.0,
                            'max_NS_mass' : 2.5,
                            'verbose' : False,
                            'path_to_datasets': path_to_Sukhbold_datasets} )

            star_prop = {'mass': M_CO / 0.7638113015667961,
                            'co_core_mass': M_CO ,
                            'he_core_mass': M_CO / 0.7638113015667961,
                            'state':'stripped_He_Central_C_depletion',
                            'profile':None ,
                            'spin': 0.0}
            star = SingleStar(**star_prop)

            star_he_core = star.he_core_mass

            M_rembar = SN.compute_m_rembar(star , None)[0]

            self.assertEqual( SN.SN_type , 'WD')
            self.assertEqual( M_rembar , star_he_core)


        """
        Test ECSN formation
        """
        def test_WD_formation_RAPID(self):
            M_CO = 1.38
            SN = StepSN( **{'mechanism' : 'Fryer+12-rapid',
                            'engine' : None,
                            'PISN' : 'Marchant+19',
                            'ECSN' : 'cosmic',
                            'max_neutrino_mass_loss' : 0.,
                            'kick' : True,
                            'sigma_kick_CCSN' : 265.0,
                            'sigma_kick_ECSN' : 20.0,
                            'max_NS_mass' : 2.5,
                            'verbose' : False} )

            star_prop = {'mass': M_CO / 0.7638113015667961,
                            'co_core_mass': M_CO ,
                            'he_core_mass': M_CO / 0.7638113015667961,
                            'state':'stripped_He_Central_C_depletion',
                            'profile':None ,
                            'spin': 0.0}
            star = SingleStar(**star_prop)

            star_c_core = star.co_core_mass

            M_rembar = SN.compute_m_rembar(star , None)[0]

            self.assertEqual( SN.SN_type , 'ECSN')
            self.assertEqual( M_rembar , star_c_core)

        def test_WD_formation_DELAYED(self):
            M_CO = 1.38
            SN = StepSN( **{'mechanism' : 'Fryer+12-delayed',
                            'engine' : None,
                            'PISN' : 'Marchant+19',
                            'ECSN' : 'cosmic',
                            'max_neutrino_mass_loss' : 0.,
                            'kick' : True,
                            'sigma_kick_CCSN' : 265.0,
                            'sigma_kick_ECSN' : 20.0,
                            'max_NS_mass' : 2.5,
                            'verbose' : False} )

            star_prop = {'mass': M_CO / 0.7638113015667961,
                            'co_core_mass': M_CO ,
                            'he_core_mass': M_CO / 0.7638113015667961,
                            'state':'stripped_He_Central_C_depletion',
                            'profile':None ,
                            'spin': 0.0}
            star = SingleStar(**star_prop)

            star_c_core = star.co_core_mass

            M_rembar = SN.compute_m_rembar(star , None)[0]

            self.assertEqual( SN.SN_type , 'ECSN')
            self.assertEqual( M_rembar , star_c_core)

        def test_WD_formation_SUKHBOLDN20(self):
            M_CO = 1.38
            SN = StepSN( **{'mechanism' : 'Sukhbold+16-engine',
                            'engine' : 'N20',
                            'PISN' : 'Marchant+19',
                            'ECSN' : 'cosmic',
                            'max_neutrino_mass_loss' : 0.,
                            'kick' : True,
                            'sigma_kick_CCSN' : 265.0,
                            'sigma_kick_ECSN' : 20.0,
                            'max_NS_mass' : 2.5,
                            'verbose' : False,
                            'path_to_datasets': path_to_Sukhbold_datasets} )

            star_prop = {'mass': M_CO / 0.7638113015667961,
                            'co_core_mass': M_CO ,
                            'he_core_mass': M_CO / 0.7638113015667961,
                            'state':'stripped_He_Central_C_depletion',
                            'profile':None ,
                            'spin': 0.0}
            star = SingleStar(**star_prop)

            star_c_core = star.co_core_mass

            M_rembar = SN.compute_m_rembar(star , None)[0]

            self.assertEqual( SN.SN_type , 'ECSN')
            self.assertEqual( M_rembar , star_c_core)


        """
        Test CCSN formation
        """
        def test_CCSN_formation_RAPID(self):
            M_CO = 2.0
            SN = StepSN( **{'mechanism' : 'Fryer+12-rapid',
                            'engine' : None,
                            'PISN' : 'Marchant+19',
                            'ECSN' : 'cosmic',
                            'max_neutrino_mass_loss' : 0.,
                            'kick' : True,
                            'sigma_kick_CCSN' : 265.0,
                            'sigma_kick_ECSN' : 20.0,
                            'max_NS_mass' : 2.5,
                            'verbose' : False} )

            star_prop = {'mass': M_CO / 0.7638113015667961,
                            'co_core_mass': M_CO ,
                            'he_core_mass': M_CO / 0.7638113015667961,
                            'state':'stripped_He_Central_C_depletion',
                            'profile':None ,
                            'spin': 0.0}
            star = SingleStar(**star_prop)

            star_he_core = star.he_core_mass

            M_rembar = SN.compute_m_rembar(star , None)[0]

            self.assertEqual( SN.SN_type , 'CCSN')

        def test_CCSN_formation_DELAYED(self):
            M_CO = 2.0
            SN = StepSN( **{'mechanism' : 'Fryer+12-delayed',
                            'engine' : None,
                            'PISN' : 'Marchant+19',
                            'ECSN' : 'cosmic',
                            'max_neutrino_mass_loss' : 0.,
                            'kick' : True,
                            'sigma_kick_CCSN' : 265.0,
                            'sigma_kick_ECSN' : 20.0,
                            'max_NS_mass' : 2.5,
                            'verbose' : False} )

            star_prop = {'mass': M_CO / 0.7638113015667961,
                            'co_core_mass': M_CO ,
                            'he_core_mass': M_CO / 0.7638113015667961,
                            'state':'stripped_He_Central_C_depletion',
                            'profile':None ,
                            'spin': 0.0}
            star = SingleStar(**star_prop)

            star_he_core = star.he_core_mass

            M_rembar = SN.compute_m_rembar(star , None)[0]

            self.assertEqual( SN.SN_type , 'CCSN')

        def test_CCSN_formation_SUKHBOLDN20(self):
            M_CO = 2.0
            SN = StepSN( **{'mechanism' : 'Sukhbold+16-engine',
                            'engine' : 'N20',
                            'PISN' : 'Marchant+19',
                            'ECSN' : 'cosmic',
                            'max_neutrino_mass_loss' : 0.,
                            'kick' : True,
                            'sigma_kick_CCSN' : 265.0,
                            'sigma_kick_ECSN' : 20.0,
                            'max_NS_mass' : 2.5,
                            'verbose' : False,
                            'path_to_datasets': path_to_Sukhbold_datasets} )

            star_prop = {'mass': M_CO / 0.7638113015667961,
                            'co_core_mass': M_CO ,
                            'he_core_mass': M_CO / 0.7638113015667961,
                            'state':'stripped_He_Central_C_depletion',
                            'profile':None ,
                            'spin': 0.0}
            star = SingleStar(**star_prop)

            star_he_core = star.he_core_mass

            M_rembar = SN.compute_m_rembar(star , None)[0]

            self.assertEqual( SN.SN_type , 'CCSN')

        """
        Test PPISN
        """
        def test_remnant_mass_PPISN(self):
            M_He = 35.0

            SN = StepSN( **{'mechanism' : 'Fryer+12-rapid',
                            'engine' : None,
                            'PISN' : 'Marchant+19',
                            'ECSN' : 'cosmic',
                            'max_neutrino_mass_loss' : 0.,
                            'kick' : True,
                            'sigma_kick_CCSN' : 265.0,
                            'sigma_kick_ECSN' : 20.0,
                            'max_NS_mass' : 2.5,
                            'verbose' : False} )

            star_prop = {'mass': M_He,
                            'co_core_mass': M_He * 0.7638113015667961 ,
                            'he_core_mass': M_He ,
                            'state':'stripped_He_Central_C_depletion',
                            'profile':None ,
                            'spin': 0.0}

            star = SingleStar(**star_prop)

            m_PISN = SN.PISN_prescription(star)

            SN.compute_m_rembar(star , m_PISN)[0]

            self.assertTrue( m_PISN > 0.0 )
            self.assertTrue( m_PISN <= 50.0 )
            self.assertEqual(SN.SN_type , 'PPISN')

        """
        Test PISN
        """
        def test_remnant_mass_PPISN(self):
            M_He = 70.0

            SN = StepSN( **{'mechanism' : 'Fryer+12-rapid',
                            'engine' : None,
                            'PISN' : 'Marchant+19',
                            'ECSN' : 'cosmic',
                            'max_neutrino_mass_loss' : 0.,
                            'kick' : True,
                            'sigma_kick_CCSN' : 265.0,
                            'sigma_kick_ECSN' : 20.0,
                            'max_NS_mass' : 2.5,
                            'verbose' : False} )

            star_prop = {'mass': M_He,
                            'co_core_mass': M_He * 0.7638113015667961 ,
                            'he_core_mass': M_He ,
                            'state':'stripped_He_Central_C_depletion',
                            'profile':None ,
                            'spin': 0.0}

            star = SingleStar(**star_prop)

            m_PISN = SN.PISN_prescription(star)

            SN.compute_m_rembar(star , m_PISN)[0]

            self.assertTrue( np.isnan(m_PISN) )
            self.assertEqual(SN.SN_type , 'PISN')

        """
        Test kick distribution for ECSN
        """
        def test_kick_ECSN(self):
            SN = StepSN( **{'mechanism' : 'Fryer+12-rapid',
                            'engine' : None,
                            'PISN' : 'Marchant+19',
                            'ECSN' : 'cosmic',
                            'max_neutrino_mass_loss' : 0.,
                            'kick' : True,
                            'sigma_kick_CCSN' : 265.0,
                            'sigma_kick_ECSN' : 20.0,
                            'max_NS_mass' : 2.5,
                            'verbose' : False} )

            SN_type = np.array([])
            Vkick = np.array([])
            M_co = np.full_like(np.arange(50000)*1.0 , 1.38)

            # The He stars are created
            for m_co in M_co:
                star_prop = {'mass':m_co / 0.7638113015667961,
                        'co_core_mass':m_co,
                        'he_core_mass':m_co / 0.7638113015667961,
                        'state':'stripped_He_Central_C_depletion',
                        'profile':None ,
                        'spin': 0.0}


                star = SingleStar(**star_prop)

                # The fallback fraction is stracted, this is not a random
                # variable then is a fixed value for all explotions
                f_fb = SN.compute_m_rembar(star, None)[1]

                # We perform the collapse to extract the SN type of the
                # from the code
                SN.collapse_star(star)

                if (SN.SN_type == 'CCSN') + (SN.SN_type == 'PPISN') :
                    kick = SN.generate_kick(star , SN.sigma_kick_CCSN)
                    sigma = SN.sigma_kick_CCSN
                elif SN.SN_type == 'ECSN':
                    kick = SN.generate_kick(star , SN.sigma_kick_ECSN)
                    sigma = SN.sigma_kick_ECSN


                SN_type = np.append(SN_type , SN.SN_type)
                Vkick = np.append(Vkick , kick)

                star = None

            dist = (Vkick[SN_type == 'ECSN'] / (1.0 - f_fb))

            sigma_ECSN = np.round(np.std(dist) / np.sqrt((3*np.pi - 8)/np.pi) , 2)

            print(sigma_ECSN)

            lower = sigma_ECSN <= (SN.sigma_kick_ECSN + 2)
            upper = sigma_ECSN >= (SN.sigma_kick_ECSN - 2)

            self.assertTrue( lower )
            self.assertTrue( upper )

        """
        Test kick distribution for CCSN
        """
        def test_kick_CCSN(self):
            SN = StepSN( **{'mechanism' : 'Fryer+12-rapid',
                            'engine' : None,
                            'PISN' : 'Marchant+19',
                            'ECSN' : 'cosmic',
                            'max_neutrino_mass_loss' : 0.,
                            'kick' : True,
                            'sigma_kick_CCSN' : 265.0,
                            'sigma_kick_ECSN' : 20.0,
                            'max_NS_mass' : 2.5,
                            'verbose' : False} )

            SN_type = np.array([])
            Vkick = np.array([])
            M_co = np.full_like(np.arange(50000)*1.0 , 8.0)

            # The He stars are created
            for m_co in M_co:
                star_prop = {'mass':m_co / 0.7638113015667961,
                        'co_core_mass':m_co,
                        'he_core_mass':m_co / 0.7638113015667961,
                        'state':'stripped_He_Central_C_depletion',
                        'profile':None ,
                        'spin': 0.0}


                star = SingleStar(**star_prop)

                # The fallback fraction is stracted, this is not a random
                # variable then is a fixed value for all explotions
                f_fb = SN.compute_m_rembar(star, None)[1]

                # We perform the collapse to extract the SN type of the
                # from the code
                SN.collapse_star(star)

                if (SN.SN_type == 'CCSN') + (SN.SN_type == 'PPISN') :
                    kick = SN.generate_kick(star , SN.sigma_kick_CCSN)
                    sigma = SN.sigma_kick_CCSN
                elif SN.SN_type == 'ECSN':
                    kick = SN.generate_kick(star , SN.sigma_kick_ECSN)
                    sigma = SN.sigma_kick_ECSN


                SN_type = np.append(SN_type , SN.SN_type)
                Vkick = np.append(Vkick , kick)

                star = None

            dist = (Vkick[SN_type == 'CCSN'] / (1.0 - f_fb))

            sigma_CCSN = np.round(np.std(dist) / np.sqrt((3*np.pi - 8)/np.pi) , 2)

            print(sigma_CCSN)

            lower = sigma_CCSN <= (SN.sigma_kick_CCSN + 2)
            upper = sigma_CCSN >= (SN.sigma_kick_CCSN - 2)

            self.assertTrue( lower )
            self.assertTrue( upper )

        """
        Test generate kick for expanding orbit
        """
        def test_generate_kick(self):
            SN = StepSN( **{'mechanism' : 'Fryer+12-rapid',
                            'engine' : None,
                            'PISN' : 'Marchant+19',
                            'ECSN' : 'cosmic',
                            'max_neutrino_mass_loss' : 0.,
                            'kick' : True,
                            'sigma_kick_CCSN' : 265.0,
                            'sigma_kick_ECSN' : 20.0,
                            'max_NS_mass' : 2.5,
                            'verbose' : False} )

            fallback = []

            sep_i = []
            ecc_i = []


            sep_f = []
            ecc_f = []
            Vsys_f = []

            # Loading the test data
            def end(binary):
                binary.event = 'END'

            properties_star1 = {"mass": 16.200984100257546, "state": "BH", "profile": None}
            properties_star2 = {"mass": 5.497560636139926,
                                "state": "stripped_He_Central_C_depletion",
                                "profile": None,
                                'he_core_mass': 5.497560636139926,
                                'co_core_mass': 4.1990989449324205}

            BH = SingleStar(**properties_star1)
            He_star = SingleStar(**properties_star2)
            properties_binary = {
                'orbital_period' : 6.182118856988261,
                'eccentricity' : 0.0,
                'separation': 39.5265173131476,
                'state' : 'ZAMS',
                'event' : 'CC2',
                'V_sys' : [0, 0, 0],
                'mass_transfer_case' : None,
            }
            binary = BinaryStar(BH, He_star, **properties_binary)

            pop = [binary]


            for i in range(len(pop)):
                binary = pop[i]

                # We consider that the kicks will have the same direction
                # as the velocity of the He star at the periapsis
                binary.star_2.natal_kick_array = [None , 0., 0., 0.]

                # We save the orbital separation end eccentricity pre-supernova
                sep_i.append( binary.separation )
                ecc_i.append( binary.eccentricity )

                # We save the fallback fraction f_fb of the remnant
                fallback.append(SN.compute_m_rembar(binary.star_2, None)[1])

                # The orbital kick is applied to the three dimensional orbit
                SN.orbital_kick(binary)

                # We save the orbital separation, eccentricity and kick velocity post-supernova
                sep_f.append( binary.separation )
                ecc_f.append( binary.eccentricity )
                Vsys_f.append( binary.V_sys )

                index = [sep_f[i] < sep_i[i] for i in range(len(sep_i))]

                # See if there is any orbit post supernova that shrinked more than one meter
                smaller_orbits = np.array(np.array(sep_i)[index] - np.array(sep_f)[index]  > 10**-8)

                orbit_comparision = np.sum(smaller_orbits)

                self.assertEqual( orbit_comparision , 0.0)
    '''


if __name__ == '__main__':
    unittest.main()
