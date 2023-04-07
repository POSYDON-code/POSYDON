'''
This class contains all functions and parameters for pulsar evolution, e.g. NS spin, B-field, etc.
'''


__authors__ = [
    "Camille Liotine <cliotine@u.northwestern.edu>",
    "Abhishek Chattaraj <a.chattaraj@ufl.edu>"
]


#### THIS IS A TEST #########


import numpy as np
import pandas as pd
from posydon.binary_evol.binarystar import BINARYPROPERTIES
from posydon.binary_evol.singlestar import STARPROPERTIES


class Pulsar:

    def __init__(self):
        '''
        Construct a PulsarBinary object, which is comprised of two Pulsar objects.
        '''
        self.spin = None
        self.Bfield = None
        self.moment_inertia = None

    def draw_NS_spin(self):
        '''
        Draw the initial NS spin angular frequency Omega from a uniform random distribution.
        Range is 0.1 - 2.0 seconds for pulse period P.
        Omega = 2*pi/P
        '''
        self.spin = np.random.uniform(2*np.pi/.1, 2*np.pi/2) ## units are in 1/s
    
    def draw_NS_Bfield(self):
        '''
        Draw the initial NS B-field from a uniform random distribution.
        Range is 10^11.5 - 10^13.8 Gauss.
        '''
        self.Bfield = np.random.uniform(3.16e11, 6.31e13) ## units are in Gauss  
    
    def evolve(self, binary):
        '''
        Evolve a pulsar from start to finish.

        Parameters
        ----------
        binary: BinaryStar object
        '''

        self.draw_NS_spin()
        self.draw_NS_Bfield()





