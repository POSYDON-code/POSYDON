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


class Pulsar:

    def __init__(self, star):
        '''
        Construct a PulsarBinary object, which is comprised of two Pulsar objects.

        Parameters
        ------------
        star: SingleStar object corresponding to the pulsar
        '''
        
        NS_RADIUS = 2.123e-5       ## CONSTANT value in POSYDON
        ## NEED TO FIGURE OUT UNITS FOR THIS

        self.spin = self.draw_NS_spin()          ## NS spin angular frequency
        self.Bfield = self.draw_NS_Bfield()      ## NS magnetic field 

        self.mass = star.mass       ## initial mass of the NS
        self.radius = NS_RADIUS 
        self.moment_inertia = None 

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
    
    def calc_moment_of_inertia(self, M, R):
            '''
            Calculate the moment of intertia for the neutron 
            '''
            

    
    
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






