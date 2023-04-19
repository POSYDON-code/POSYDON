'''
This class contains all functions and parameters for pulsar evolution, e.g. NS spin, B-field, etc.
'''


__authors__ = [
    "Camille Liotine <cliotine@u.northwestern.edu>",
    "Abhishek Chattaraj <a.chattaraj@ufl.edu>"
]

import numpy as np
from posydon.utils import constants as const
from posydon.utils.common_functions import CO_radius


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
        Construct a Pulsar object.

        Parameters
        ------------
        star: SingleStar object corresponding to the pulsar
        '''
        
        NS_RADIUS = CO_radius(star.mass, "NS")*const.Rsun    ## POSYDON constant for NS radius [cm] 

        self.mass = star.mass*const.Msun                      ## mass of the NS [g]
        self.radius = NS_RADIUS                               ## radius of the NS [cm]
        self.moment_inertia = self.calc_moment_of_inertia()   ## moment of inertia of the NS
        ## Note: Because the NS radius is constant, moment of inertia should also be constant throughout NS evolution

        self.spin = self.draw_NS_spin()          ## NS spin angular frequency [1/s]
        self.Bfield = self.draw_NS_Bfield()      ## NS magnetic field [G]

        self.alive_state = self.is_alive()
 
        
    def draw_NS_spin(self):
        '''
        Draw the initial NS spin angular frequency Omega [1/s] from a uniform random distribution.
        Range is 0.1-2.0 sec for pulse period P.
        Omega = 2*pi/P
        '''
        return np.random.uniform(2*np.pi/.1, 2*np.pi/2) 

    def draw_NS_Bfield(self):
        '''
        Draw the initial NS B-field [G] from a uniform random distribution.
        Range is 10^11.5 - 10^13.8 Gauss.
        '''
        return np.random.uniform(3.16e11, 6.31e13)
    
    def calc_moment_of_inertia(self):
        '''
        Calculate the moment of intertia of the NS in g*cm^2
        Eq. from Kiel et al. 2008 is in solar units
        '''
        M = self.mass/const.Msun       ## mass of the NS [Msun] 
        R = self.radius/const.Rsun     ## radius of the NS [Rsun]

        return 2/7 * (1 - 2.42e-6*M/R - 2.9e-12*M**2/R**2)**-1 * M*R**2 * (const.Msun*const.Rsun**2)
    
    def detached_evolve(self, delta_t):
        '''
        Evolve a pulsar during detached evolution.

        Parameters
        ----------
        delta_t: the duration of the detached evolution phase [s]
        '''
        alpha = 45*const.a2rad          ## angle between the axis of rotation and magnetic axis [rad]
        mu_0 = 1                        ## permeability of free space [has value unity in cgs]
        B_min = 1e8                     ## minimum Bfield strength at which Bfield decay ceases [G]
        tau_d = 3*1e9*const.secyer      ## B-field decay timescale = 3 Gyr [s]

        I = self.moment_inertia
        R = self.radius
        B_i = self.Bfield

        ## evolve the NS B-field
        B_f = (B_i - B_min) * np.exp(-delta_t/tau_d) + B_min
        self.Bfield = B_f

        ## evolve the NS spin
        A = 8*np.pi*R**6*np.sin(alpha)**2/(3*mu_0*const.clight**3*I)
        omega_f = np.sqrt(1/( A*(B_min**2*delta_t - tau_d*B_min*(B_f - B_i) - tau_d/2*(B_f**2 - B_i**2)) + 1/self.spin**2)) 
        self.spin = omega_f

        ## check if pulsar crossed the death line
        self.alive_state = self.is_alive()

    
    def RLO_evolve(self, delta_t, T, delta_M):
        '''
        Evolve a pulsar during Roche Lobe overflow (RLO).

        Parameters
        ----------
        delta_t: the duration of the RLO accretion phase [s]
        T: the age of the NS when RLO begins [s]
        delta_m: the total amount of mass accreted by the pulsar during RLO [g]
        '''
        G = const.standard_cgrav     ## gravitational constant [cm^3 g^-1 s^-2]
        tau_d = 3*1e9*const.secyer   ## B-field decay timescale = 3 Gyr [s]

        M_i = self.mass              ## mass of the NS before accretion [g]
        R = self.radius              ## radius of the NS [cm]

        ## evolve the NS spin
        J_i = M_i*R**2*self.spin     ## spin angular momentum (J) of the NS before accretion
        
        omega_k = np.sqrt(G*M_i/R**3)
        delta_J = delta_M*R**2*omega_k    ## change in J due to accretion

        J_f = J_i + delta_J
        M_f = M_i + delta_M

        omega_f = J_f/(M_f*R**2)
        self.spin = omega_f

        ## evolve the NS B-field
        B_f = self.Bfield/(1 + delta_M/(1e-6*const.Msun)) * np.exp(-(T-delta_t)/tau_d)
        self.Bfield = B_f

        ## check if pulsar crossed the death line
        self.alive_state = self.is_alive()


    def CE_evolve(self, T):
        '''
        Evolve a pulsar during common envelope.

        Parameters
        ----------
        T: the age of the NS when CE begins [s]
        ''' 
        delta_M = 0.1*const.Msun  ## assume amount of mass accreted during CE = 0.1 Msun
    
        ## params needed for RLO evolve, assume CE phase is instantaneous
        delta_t = 0               

        self.RLO_evolve(delta_t, T, delta_M)


    def is_alive(self):
        '''
        Check if the pulsar has crossed the death line.
        ''' 

        P = 2*np.pi/self.spin     ## spin period of the pulsar [s]
        death_line = 0.17e12

        if (self.Bfield/P**2 < death_line): return False
        else: return True



