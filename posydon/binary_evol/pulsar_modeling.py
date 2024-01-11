
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

[4] Kiel, P. D., Hurley, J. R., Bailes, M., & Murray, J. R. 2008,
MNRAS, 388, 393, doi: 10.1111/j.1365-2966.2008.13402.x

"""

class Pulsar:

    def __init__(self, initial_mass):
        '''
        Construct a Pulsar object.

        Parameters
        ------------
        initial_mass: initial mass of the NS [Msun]
        initial_Mdot: initial rate of change in mass of the NS [Msun/yr]
        '''

        NS_RADIUS = CO_radius(initial_mass, "NS")*const.Rsun  ## POSYDON constant for NS radius [cm] 
        
        self.mass = initial_mass*const.Msun                   ## mass of the NS [g]
        self.radius = NS_RADIUS                               ## radius of the NS [cm]
        self.moment_inertia = self.calc_moment_of_inertia()   ## moment of inertia of the NS  [g*cm^2]
        ## Note: Because the NS radius is constant, moment of inertia should also be constant throughout NS evolution

        self.Mdot_edd = self.calc_NS_edd_lim(np.nan)          ## Eddington accretion rate for the NS [g/s]
        self.luminosity = self.calc_NS_luminosity()           ## radio luminosity of the pulsar [erg/s]

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
        Draw the initial NS B-field [G] from a lognormal distribution.
        Range is 10^11.5 - 10^13.8 Gauss.
        '''
        # return np.random.uniform(3.16e11, 6.31e13)
        return 10**np.random.uniform(11.5, 13.8)

    def calc_moment_of_inertia(self):
        '''
        Calculate the moment of intertia of the NS [g*cm^2]
        Eq. from Kiel et al. 2008 is in solar units
        '''
        M = self.mass/const.Msun       ## mass of the NS [Msun] 
        R = self.radius/const.Rsun     ## radius of the NS [Rsun]

        return 2/7 * (1 - 2.42e-6*M/R - 2.9e-12*M**2/R**2)**-1 * M*R**2 * (const.Msun*const.Rsun**2)
    
    def calc_NS_edd_lim(self, surface_h1):
        """
        Calculate the Eddington accretion rate for a NS [g/s]
        Adapted from eddington_limit() in utils.common_functions

        Parameters
        ----------
        surface_h1: surface hydrogen fraction of donor star
        If NaN, assume = 0.7155 (approximate surface H1 abundance for HMS stars @ Zsun)
        """
        if np.isnan(surface_h1): surface_h1 = 0.7155

        eta = const.standard_cgrav*self.mass / (self.radius * const.clight**2)
        Mdot_edd = (4*np.pi*const.standard_cgrav*self.mass) / (0.2*(1 + surface_h1)*eta*const.clight)  
        return Mdot_edd

    def calc_NS_luminosity(self):
        """
        Calculate the radio luminosity of the NS at 1400 MHz.
        units of distribution from Szary et al. 2014 are mJy*kpc^2
        """
        kpc2cm = 3.086e21           ## kpc to cm
        mJy2cgs = 1e-23*1e-3        ## milliJansky in CGS units [erg/s/cm^2/Hz]
        Hz = 1400*1e6               ## Hz conversion for radio luminosity (1400 MHz)
        mJykpcsq2ergs = mJy2cgs*kpc2cm**2*Hz

        lum_draw = np.random.lognormal(0.5, 1.0)
        luminosity = 10**lum_draw*mJykpcsq2ergs
        return luminosity
    
    def calc_NS_spindown_power(self):
        """
        Calculate the spindown power of the NS [erg/s].
        From Szary et al. 2014
        """
        P = 2*np.pi/self.spin                                ## NS spin period
        P_dot =  9.87e-48*const.secyer * self.Bfield**2/P    ## NS spindown rate

        power = 4*np.pi**2 * self.moment_inertia * P_dot / P**3
        return power      
    
    def detached_evolve(self, delta_t, tau_d):
        '''
        Evolve a pulsar during detached evolution.

        Parameters
        ----------
        delta_t: the duration of the detached evolution phase [s]
        '''
        alpha = 45*const.a2rad          ## angle between the axis of rotation and magnetic axis [rad]
        mu_0 = 1                        ## permeability of free space [has value unity in cgs]
        B_min = 1e8                     ## minimum Bfield strength at which Bfield decay ceases [G]
        tau_d *= (1e9*const.secyer)     ## B-field decay timescale [s]

        delta_t *= const.secyer

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

    
    def RLO_evolve_Ye2019(self, delta_t, tau_d, delta_M, delta_Md):
        '''
        Evolve a pulsar during Roche Lobe overflow (RLO).

        Parameters
        ----------
        delta_t: the duration of the RLO accretion phase [yr]
        delta_M: the total amount of mass accreted by the pulsar during RLO [Msun]
        '''

        G = const.standard_cgrav     ## gravitational constant [cm^3 g^-1 s^-2]
        mu_0 = 1                     ## permeability of free space [has value unity in cgs]
        tau_d *= (1e9*const.secyer)  ## B-field decay timescale [s]
        delta_Md *= const.Msun       ## magnetic field mass decay scale [g]

        delta_M *= const.Msun
        delta_t *= const.secyer

        M_i = self.mass              ## mass of the NS before accretion [g]
        B_i = self.Bfield            ## B-field of the NS before accretion [G]
        R = self.radius              ## radius of the NS [cm]
        I = self.moment_inertia      ## Moment of inertia [CGS]

        R_alfven = (2*np.pi**2/(G*mu_0**2))**(1/7) * (R**6/(self.Mdot_edd*M_i**(1/2)))**(2/7) * B_i**(4/7) ## Alfven radius
        R_mag = R_alfven/2   ## magnetic radius

        ## evolve the NS spin
        #J_i = 2/5*M_i*R**2*self.spin     ## spin angular momentum (J) of the NS before accretion
        J_i = I*self.spin
        
        omega_k = np.sqrt(G*M_i/R**3)
        delta_J = 2/5*delta_M*R_mag**2*omega_k    ## change in J due to accretion

        J_f = J_i + delta_J
        M_f = M_i + delta_M
        self.mass = M_f

        #omega_f = J_f/(2/5*M_f*R**2)
        omega_f = J_f/I
        self.spin = omega_f

        ## evolve the NS B-field
        B_f = B_i/(1 + delta_M/delta_Md) * np.exp(-delta_t/tau_d)
        self.Bfield = B_f

        ## check if pulsar has crossed the death line
        self.alive_state = self.is_alive()
    
    def RLO_evolve_COMPAS(self, delta_M, delta_Md):
        '''
        Evolve a pulsar during Roche Lobe overflow (RLO).
        This uses the prescription for B-field decay applied in COMPAS from Oslowski et al. 2011.
        Spin-down is the same for now.

        Parameters
        ----------
        delta_M: the total amount of mass accreted by the pulsar during RLO [Msun]
        '''
        G = const.standard_cgrav     ## gravitational constant [cm^3 g^-1 s^-2]
        delta_Md *= const.Msun       ## magnetic field mass decay scale [g]
        B_min = 1e8                  ## minimum Bfield strength at which Bfield decay ceases [G]
        mu_0 = 1                     ## permeability of free space [has value unity in cgs]

        delta_M *= const.Msun        ## convert Msun to g

        M_i = self.mass              ## mass of the NS before accretion [g]  
        R = self.radius              ## radius of the NS [cm]
        B_i = self.Bfield            ## B-field of the NS before accretion [G]
        I = self.moment_inertia

        R_alfven = (2*np.pi**2/(G*mu_0**2))**(1/7) * (R**6/(self.Mdot_edd*M_i**(1/2)))**(2/7) * B_i**(4/7) ## Alfven radius
        R_mag = R_alfven/2   ## magnetic radius
   
        ## evolve the NS spin
        #J_i = 2/5*M_i*R**2*self.spin     ## spin angular momentum (J) of the NS before accretion
        J_i = I*self.spin
        
        omega_k = np.sqrt(G*M_i/R_mag**3)
        delta_J = 2/5*delta_M*R_mag**2*omega_k    ## change in J due to accretion

        J_f = J_i + delta_J
        M_f = M_i + delta_M
        self.mass = M_f

        #omega_f = J_f/(2/5*M_f*R**2)
        omega_f = J_f/I
        self.spin = omega_f

        ## evolve the NS B-field
        B_f = (B_i - B_min)*np.exp(-delta_M/delta_Md) + B_min 
        self.Bfield = B_f

        ## check if pulsar crossed the death line
        self.alive_state = self.is_alive()

    def CE_evolve(self, CE_acc_prescription, acc_decay_prescription, M_comp, R_comp, delta_Md, delta_t, tau_d):
        '''
        Evolve a pulsar during common envelope, accounting for mass accretion onto the NS.
        '''   
        if CE_acc_prescription == "None": delta_M = 0
        
        elif CE_acc_prescription == "uniform":
            ## assume amount of mass accreted during CE is 0.04-0.1 Msun  
            delta_M = np.random.uniform(0.04, 0.1)  

        elif CE_acc_prescription in ["MacLeod", "MacLeod_bounded"]:
            ## use the MacLeod prescription from COMPAS paper, fit to Fig. 4 in Macleod & Ramirez-Ruiz 
            a_a = -1.1e-5; a_b = 1.5e-2; b_a = 1.2e-4; b_b = -1.5e-1

            a = a_a*M_comp + b_a
            b = a_b*M_comp + b_b
            delta_M = np.abs(a*R_comp + b)

            ## assume amount of mass accreted during CE is 0.04-0.1 Msun  
            if delta_M > 0.1: delta_M = 0.1

            if CE_acc_prescription == "MacLeod_bounded":
                if delta_M < 0.04: delta_M = 0.04

        if acc_decay_prescription == "Ye2019":
            self.RLO_evolve_Ye2019(delta_t, tau_d, delta_M, delta_Md)
        elif acc_decay_prescription == "COMPAS":
            self.RLO_evolve_COMPAS(delta_M, delta_Md)

    def is_alive(self):
        '''
        Check if the pulsar has crossed the death line.
        ''' 
        P = 2*np.pi/self.spin     ## spin period of the pulsar [s]
        P_dot = 9.87e-48*const.secyer * self.Bfield**2/P  

        death_line = 0.17e12*P**2   ## death line from Ruderman & Sutherlandt 1975

        E_max = 0.01   ## threshold radio efficiency
        L = self.luminosity
        Edot = self.calc_NS_spindown_power()
  
        if ((L/Edot) < E_max): return True
        #elif (self.Bfield < death_line): return False
        else: return False
       




