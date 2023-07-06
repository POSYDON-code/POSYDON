"""
Spectral Synthesis code 
"""
#################################### Imports ##############################################

import numpy as np
import astropy.constants as con
import astropy.units as unt
import matplotlib.pyplot as plt
import h5py
from astropy.io import fits

import math
import traceback
import copy
from spectral_tools import array_data,star
from spectral_grids import spectral_grids


#########################################Creating the class population_spectra ###############################################


class population_spectra():
    
    def __init__(self,total_binaries,first_binary_line,last_binary_line,metallicity,scaling_factor,data_file=None):
        ######Initializing the parameters ##################
        self.total_binaries = total_binaries
        self.data_file = data_file
        self.first_binary_line = first_binary_line
        self.last_binary_line = last_binary_line
        self.scaling_factor = scaling_factor 
        self.metallicity = metallicity
        self.population = [None]*total_binaries
        self.failed_stars = 0 #int, the stars that failed durin the spectra process 
        self.missing_stars = 0 #The stars that are missing due to POSYDON, some detached cases for now. 
        #Creating lists/numpy arrays for investigation for the failled and succesful stars. 
        self.stars_fails = []
        self.stars_run   = []
        self.stars_grid_fail = []
        
        #Creating readable arrays for the stars objects. 
        mass,radius,L,state = array_data(self.total_binaries,self.last_binary_line,self.first_binary_line)
        for i in range(self.total_binaries):
            star1 = star(i,0,mass[i,0],state[i,0],radius[i,0],L[i,0],self.metallicity)
            star2 = star(i,1,mass[i,1],state[i,1],radius[i,1],L[i,1],self.metallicity)
            self.population[i] = [star1,star2]
            #self.population.append((star1,star2))
        
        #Initializing the spectral_grids object and parameters used.
        #To do put an option for changing the wavelength 
        self.grids = spectral_grids()
        self.grids.global_limits()

    


    def create_spectrum_single(self,star,**kwargs):
        #global bin_num
        #bin_num = binary_number

        if "stripped" in star.state:
            M = star.mass/con.M_sun
            M_min = self.grids.specgrid_stripped.axis_x_min['M_init']
            M_max = self.grids.specgrid_stripped.axis_x_max['M_init']
            if M < M_min or M > M_max:
                self.failed_stars +=1
                return None
            return self.grids.stripped_grid_flux(star)
        
        Fe_H = star.Fe_H
        Z_Zo = star.metallicity
        Teff = star.get_Teff(self.grids.T_max,self.grids.T_min)
        logg = copy.copy(star.get_logg(self.grids.logg_max,self.grids.logg_min))
        if Teff == None or logg == None:
            """
                if create_logg(binary_number,1) > 5 or create_Teff(binary_number,1) > max_Teff or create_Teff(binary_number,1)<min_Teff:
                    self.failed_stars +=1
                    #stars_fail.append(create_logg(binary_number,1),create_Teff(binary_number,1))
                    #stars_fails_logg.append(create_logg(binary_number,1))
                    #stars_fails_Teff.append(create_Teff(binary_number,1))
                elif create_logg(binary_number,2) > 5 or create_Teff(binary_number,2) > max_Teff or create_Teff(binary_number,2) < min_Teff:
                    self.failed_stars +=1
                    #stars_fail.append(create_logg(binary_number,2),create_Teff(binary_number,2))
                    #stars_fails_logg.append(create_logg(binary_number,2))
                    #stars_fails_Teff.append(create_Teff(binary_number,2))
                else:
            """
            self.failed_stars +=1 
            return None
        
        #For temperature higher than then lo_limits of Teff in main grids we calculate the spectra using the Ostar girds. 
        if Teff >= 30000:
            #Setting the acceptable boundaries for the ostar grid in the logg.
            logg_min = self.grids.specgrid_ostar.axis_x_min['log(g)']
            logg_max = self.grids.specgrid_ostar.axis_x_max['log(g)']
            if logg > logg_min and logg<logg_max: 
                    try:
                        Flux = self.grids.ostar_grid_flux(Teff,logg,star,Z_Zo)
                        return Flux
                    except:
                        self.failed_stars +=1
                        return None
            else:
                self.failed_stars +=1
                return None
            
        try:
            Flux = self.grids.main_grid_flux(Teff,Fe_H,logg,star)
            #stars_run_logg.append(star.get_logg())
            #stars_run_Teff.append(star.get_Teff())              
            return Flux
        except Exception as e:
            try:
                # lo limits for the secondary grid: 
                logg_min = self.grids.specgrid_secondary.axis_x_min['log(g)']
                Teff_min = self.grids.specgrid_secondary.axis_x_min['Teff']
                #Calculating all the exeption that the first grid gave with the secondary. 
                if Teff>15000.0 and logg > logg_min:
                    F = self.grids.secondary_grid_flux(Teff,logg,Z_Zo=1,star=star)
                    #stars_run_logg.append(star.get_logg())
                    #stars_run_Teff.append(star.get_Teff())
                    return F
            finally:
                # if and else statements that fix the grid voids.  
                if Teff > 20000:
                    logg = max(logg, 4.0)
                elif Teff > 12000:
                    logg = max(logg, 3.0)
                elif Teff > 8000:
                    logg = max(logg, 2.0)
                elif Teff > 6000:
                    logg = max(logg, 1.0)
                try:
                    if Teff >  self.grids.specgrid_main.axis_x_max['Teff']:
                        F = self.grids.secondary_grid_flux(Teff,Fe_H,logg,star)
                    else: 
                        F = self.grids.main_grid_flux(Teff,Fe_H,logg,star)
                    #stars_run.append(star.get_logg(),star.get_Teff())
                    #stars_run_logg.append(star.get_logg())
                    #stars_run_Teff.append(star.get_Teff())
                    return F
                #*np.exp(logg/star.get_logg(self.grids.logg_max,self.grids.logg_min))
                except Exception as e: 
                    print(e)
                    self.failed_stars +=1
                    return None



    def add_spectra(self,spectrum_1,spectrum_2):
        if spectrum_1 is not None and spectrum_2 is not None:
            return spectrum_1 + spectrum_2
        elif spectrum_1 is None:
            return spectrum_2
        elif spectrum_2 is None:
            return spectrum_1
        else:
            raise Exception("Spectrum 1 and Spectrum 2 have not allowed values.")
    
    def create_spectrum_binary(self,binary):
    # check if binary has two or any non-degenerate stars
        star1 = binary[0]
        star2 = binary[1]
        if star1.state not in ['NS','BH','Failed']:
            spectrum_1 = self.create_spectrum_single(star1)
        else:
            spectrum_1 = None

        if star2.state not in ['NS','BH','Failed']:
            spectrum_2 = self.create_spectrum_single(star2)
        else: 
            spectrum_2 = None
        return self.add_spectra(spectrum_1, spectrum_2)


    def create_spectrum_population(self):
        create_spectrum_population = np.zeros(len(self.grids.lam_c))*unt.m/unt.m*0
        failed_binaries = 0
        for binary in self.population:
            spectrum = self.create_spectrum_binary(binary)
            if spectrum is not None:
                create_spectrum_population += spectrum
        return create_spectrum_population , self.grids.lam_c
       

    def create_spectrum_subpopulation(self,state):
        create_spectrum_population = np.zeros(len(self.grids.lam_c))*unt.m/unt.m*0
        failed_binaries = 0
        for binary in self.population:
            star1 = binary[0]
            if star1.binary_state == state:
                spectrum = self.create_spectrum_binary(binary)
                if spectrum is not None:
                    self.create_spectrum_population += spectrum
        return create_spectrum_population , self.grids.lam_c

