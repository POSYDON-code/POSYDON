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

## Test push to remote branch

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

"""

class Pulsar:

    def __init__(self, star):
        '''
        Construct a PulsarBinary object, which is comprised of two Pulsar objects.

        Parameters
        ------------
        star: SingleStar object corresponding to the pulsar
        '''
        ## NEED TO FIGURE OUT UNITS
        NS_RADIUS = 2.123e-5       ## CONSTANT value in POSYDON

        self.spin = self.draw_NS_spin()          ## NS spin angular frequency
        self.Bfield = self.draw_NS_Bfield()      ## NS magnetic field 

        self.mass = star.mass               ## mass of the NS
        self.radius = NS_RADIUS 
        self.moment_inertia = self.calc_moment_of_inertia()   ## moment of inertia, STAYS CONSTANT

    def draw_NS_spin(self):
            '''
            Draw the initial NS spin angular frequency Omega from a uniform random distribution.
            Range is 0.1 - 2.0 seconds for pulse period P.
            Omega = 2*pi/P
            '''
            return np.random.uniform(2*np.pi/.1, 2*np.pi/2) ## units are in 1/s 

    def draw_NS_Bfield(self):
            '''
            Draw the initial NS B-field from a uniform random distribution.
            Range is 10^11.5 - 10^13.8 Gauss.
            '''
            return np.random.uniform(3.16e11, 6.31e13) ## units are in Gauss
    
    def calc_moment_of_inertia(self):
            '''
            Calculate the moment of intertia for the neutron 
            '''
            M = self.mass
            R = self.radius

            return 2/7 * (1 - 2.42e-6*M/R - 2.9e-12*M**2/R**2)**-1 * M*R**2
    
    def detached_evolve(self, delta_t):
        '''
        Evolve a pulsar from start to finish.

        Parameters
        ----------
        binary: BinaryStar object
        '''
        ## constants
        alpha = 45       ## angle between the axis of rotation and magnetic axis [deg]
        tau_d = 3        ## Bfield decay timescale [s]
        mu_0 = 1         ## permeability of free space [unity in cgs]
        c = 3e10         ## speed of light [cm/s]
        B_min = 1e8      ## minimum Bfield strength at which Bfield decay ceases [G]

        omega_i = self.spin
        I = self.moment_inertia
        R = self.radius

        #detached_state_start = binary.time_history[np.where((binary.state_history == 'detached') & (star.state_history == 'NS'))]
        #delta_t = 

        omega_f = np.sqrt(1/(8*np.pi*R**6*np.sin(alpha)**2/(3*mu_0*c**3*I) * B_min**2*delta_t + 1/omega_i**2)) 
        self.spin = omega_f

    
    def RLO_evolve(self, delta_m):
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






