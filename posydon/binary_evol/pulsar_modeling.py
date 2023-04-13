'''
This class contains all functions and parameters for pulsar evolution, e.g. NS spin, B-field, etc.
'''


__authors__ = [
    "Camille Liotine <cliotine@u.northwestern.edu>",
    "Abhishek Chattaraj <a.chattaraj@ufl.edu>"
]

import numpy as np
from posydon.binary_evol.binarystar import BINARYPROPERTIES
from posydon.binary_evol.singlestar import STARPROPERTIES
import posydon.utils.constants as constants


"""
References
----------

[1] Chattopadhyay, D., Stevenson, S., Hurley, J. R., Rossi,
L. J., & Flynn, C. 2020, MNRAS, 494, 1587,
doi: 10.1093/mnras/staa756

[2] Ye, C. S., Kremer, K., Chatterjee, S., Rodriguez, C. L., &
Rasio, F. A. 2019, ApJ, 877, 122,
doi: 10.3847/1538-4357/ab1b21

[3] Zhang, C. M., & Kojima, Y. 2006, MNRAS, 366, 137,
doi: 10.1111/j.1365-2966.2005.09802.x

[4] Kiel, P. D., Hurley, J. R., Bailes, M., & Murray, J. R. 2008,
MNRAS, 388, 393, doi: 10.1111/j.1365-2966.2008.13402.x

"""

class Pulsar:

    def __init__(self, star):
        '''
        Construct a PulsarBinary object, which is comprised of two Pulsar objects.

        Parameters
        ------------
        star: SingleStar object corresponding to the pulsar
        '''

        def draw_NS_spin():
            '''
            Draw the initial NS spin angular frequency Omega from a uniform random distribution.
            Range is 0.1 - 2.0 seconds for pulse period P.
            Omega = 2*pi/P
            '''
            return np.random.uniform(2*np.pi/.1, 2*np.pi/2) ## units are in 1/s

        def draw_NS_Bfield():
            '''
            Draw the initial NS B-field from a uniform random distribution.
            Range is 10^11.5 - 10^13.8 Gauss.
            '''
            return np.random.uniform(3.16e11, 6.31e13) ## units are in Gauss 

        def calculate_NS_radius():
            '''
            Using NS mass-radius relationship from Kiel et al. (2008)
            The original equation takes uses solar units
            Here, we have returned the radius in km
            '''
            return 2.126*10**-5*(star.mass**(-1/3))*constants.Rsun*10**-5

        def calculate_moment_of_inertia():
           '''
           Calculate moment of intertia for the neutron star
           Using relation given in Kiel et al. (2008)
           Here, we have returned the moment of inertia in CGS units
           '''

           I_0 = star.mass*((self.radius*10**5/constants.Rsun)**2)
           I_1 = 2.42*10**(-6)*star.mass/(self.radius*10**5/constants.Rsun)
           I_2 = 2.9*10**(-12)*(star.mass/(self.radius*10**5/constants.Rsun))**2

           I = (2/7)*I_0*(1-I_1-I_2)**(-1)

           return I*constants.Msun*(constants.Rsun**2)


        self.spin = draw_NS_spin()                                 ## NS spin angular frequency
        self.Bfield = draw_NS_Bfield()                             ## NS magnetic field (G)
        self.mass = star.mass                                      ## Initial mass of the NS (M_sun)
        self.radius = calculate_NS_radius()                        ## Radius calculated using an M-R relation (km)
        self.moment_of_inertia = calculate_moment_of_inertia()     ## Moment of inertia (g cm^2)

    
    
    def detached_evolve(self, binary):
        '''
        Evolve a pulsar from start to finish.

        Parameters
        ----------
        binary: BinaryStar object
        '''
    
    def RLO_evolve(self, binary):
        '''
        Evolve a pulsar from start to finish.

        Parameters
        ----------
        binary: BinaryStar object
        '''

    def CE_evolve(self, binary):
        '''
        Evolve a pulsar from start to finish.

        Parameters
        ----------
        binary: BinaryStar object
        ''' 

        ## call RLO_evolve after CE parameters (Delta M) are set






