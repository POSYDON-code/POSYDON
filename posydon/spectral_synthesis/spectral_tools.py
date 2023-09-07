#This is tools and useful funtion that are using in the spectral synthesis and the creation of CMD diagrams combined. 
import numpy as np
import astropy.constants as con
import astropy.units as unt
import pandas as pd
import time 

#Loading the population data in an array of binaries. 
# Inputs: the h5 file of the population and the END time. 

Zo = 0.0142



class star():
    
    def __init__(self,**kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        self.logg = None
        self.Teff = None
        self.Fe_H = np.log10(self.metallicity)
        
        
    def get_logg(self,max_logg,min_logg): 
        self.logg = np.log10(con.G*self.mass/self.R**2/(unt.cm/unt.s**2))
        if self.logg < max_logg and self.logg > min_logg:
            return self.logg
        else:
            return None
    
    def get_Teff(self,max_Teff,min_Teff):
        self.Teff = (self.L/(4*np.pi*self.R**2*con.sigma_sb))**0.25/unt.K
        if self.Teff < max_Teff and self.Teff > min_Teff:
            return self.Teff
        else:
            return None
        
    def set_metallicity(self,new_Fe_H,max_Fe_H,min_Fe_H):
        if new_Fe_H < max_Fe_H and new_Fe_H > min_Fe_H:
            self.Fe_H = new_Fe_H
        else: 
            raise TypeError
    #This is for the stars that we can't calculate their spectra so we change their logg and thus we need a new radius for them and we set and return. 
    def new_radius(self,new_logg):
        R = np.sqrt(con.G*self.mass/10**(new_logg)/unt.cm *unt.s**2).decompose()
        return R



def population_data(population_file = None ,time = None, number_of_binaries = True,**kwargs):
    if population_file is None or time is None: 
        raise Exception('File of population or time not provided')
    history = pd.read_hdf(population_file, key='history')
    final_stars = history[(history.time == time ) & (history.event == "END")].reset_index()
    zams_stars = history[history.event == 'ZAMS']
    if number_of_binaries:
        total_binaries = len(final_stars)
    elif number_of_binaries <= len(final_stars):
        total_binaries = number_of_binaries
    else: 
        raise Exception('The total number of binaries {N} is less than the input number_of_binaries {number_of_binaries}'.format(N = total_binaries, number_of_binaries = number_of_binaries  ))

    star1_properties = {}
    star2_properties = {}
    population = [None]*total_binaries
    for i in range(total_binaries):
        star1_properties['mass'],star2_properties['mass'] =  final_stars.S1_mass[i]*con.M_sun,final_stars.S2_mass[i]*con.M_sun
        star1_properties['R'],star2_properties['R']  =  10**final_stars.S1_log_R[i]*con.R_sun, 10**final_stars.S2_log_R[i]*con.R_sun
        star1_properties['L'],star2_properties['L'] = 10**final_stars.S1_log_L[i]*con.L_sun, 10**final_stars.S2_log_L[i]*con.L_sun
        star1_properties['binary_index'],star2_properties['binary_index'] = final_stars.binary_index[i],final_stars.binary_index[i]
        star1_properties['state'],star2_properties['state'] = final_stars.S1_state[i], final_stars.S2_state[i]
        star1_properties['binary_state'],star2_properties['binary_state'] = final_stars.state[i], final_stars.state[i]
        star1_properties['star_number'],star2_properties['star_number'] = 1,2
        star1_properties['binary_number'],star2_properties['binary_number'] = i,i
        star1_properties['metallicity'],star2_properties['metallicity'] = final_stars.S1_metallicity[i]/Zo, final_stars.S2_metallicity[i]/Zo

        if "stripped" in star1_properties['state']:
            star1_properties['mass'] = zams_stars.S1_mass[final_stars.binary_index[i]]*con.M_sun
        elif "stripped" in star2_properties['state']:
            star2_properties['mass'] = zams_stars.S2_mass[final_stars.binary_index[i]]*con.M_sun

        star1 = star(**star1_properties)
        star2 = star(**star2_properties)
        population[i] = [star1,star2]
    return population


def grid_global_limits(spectral_grids):
    T_max  = 0 
    T_min = 100000
    logg_max = 0 
    logg_min = 20
    for key in spectral_grids:
        specgrid = spectral_grids[key]
        for label in specgrid.axis_labels:
            if label== 'Teff':
                T_max = max(T_max,specgrid.axis_x_max[label])
                T_min = min(T_min,specgrid.axis_x_min[label])
            elif label == 'log(g)':
                logg_max = max(logg_max,specgrid.axis_x_max[label])
                logg_min = min(logg_min,specgrid.axis_x_min[label])
    return T_max,T_min,logg_max,logg_min 

