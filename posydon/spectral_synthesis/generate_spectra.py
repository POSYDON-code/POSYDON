"""
Spectral Synthesis code 


"""


#################################### Imports ##############################################
import sys
import posydon
import os
import numpy as np
import astropy.constants as con
import astropy.units as unt
import matplotlib.pyplot as plt
import h5py
from astropy.io import fits
import pandas as pd
import math
import datetime
import matplotlib.patches as mpatches
import traceback
import copy
#Connect to MSG path.
os.environ['MSG_DIR'] = '/home/kasdaglie/blue/kasdaglie/msg-1.1.2'
MSG_DIR = os.environ['MSG_DIR']

sys.path.insert(0, os.path.join(MSG_DIR, 'python'))
import pymsg

# Load the SpecGrid:

GRID_DIR = os.path.join(MSG_DIR, 'data', 'grids')

specgrid_file_name = os.path.join(GRID_DIR, 'sg-demo.h5')

specgrid = pymsg.SpecGrid(specgrid_file_name)
#%matplotlib inline
plt.rcParams.update({'font.size': 16})
kpc = 3.08e19*unt.m
##############################################################################################################################
history1 = pd.read_hdf("/home/kasdaglie/blue/kasdaglie/populations_data/metallicity_1.h5", key='history')
online1 = pd.read_hdf("/home/kasdaglie/blue/kasdaglie/populations_data/metallicity_1.h5", key='history')
last_binary_line1 = history1[history1.event == 'END']
first_binary_line1 = history1[history1.event == 'ZAMS']
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
        self.F_empty = self.grids.F_empty

    


    def new_create_spectrum_single(self,star,binary_number,**kwargs):
        #global bin_num
        #bin_num = binary_number
        F_empty = self.F_empty*0
        new_create_spectrum_single = F_empty
        if "stripped" in star.state:
            M = star.mass/con.M_sun
            M_min = self.grids.specgrid_stripped.axis_x_min['M_init']
            M_max = self.grids.specgrid_stripped.axis_x_max['M_init']
            if M < M_min or M > M_max:
                self.failed_stars +=1
                return F_empty
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
            return F_empty
        
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
                        return F_empty
            else:
                self.failed_stars +=1
                return F_empty
            
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
                    return F_empty*0  




    def add_spectra(self,spectrum_1,spectrum_2):
        return spectrum_1 + spectrum_2
    
    def create_spectrum_binary(self,binary):
    # check if binary has two or any non-degenerate stars
        star1 = binary[0]
        star2 = binary[1]
        if star1.state not in ['NS','BH','Failed']:
            spectrum_1 = self.new_create_spectrum_single(star1,star1.binary_number)
        else:
            spectrum_1 = self.F_empty*0

        if star2.state not in ['NS','BH','Failed']:
            spectrum_2 = self.new_create_spectrum_single(star2,star2.binary_number)
        else: 
            spectrum_2 = self.F_empty*0

        return self.add_spectra(spectrum_1, spectrum_2)


    def create_spectrum_population(self):
        create_spectrum_population = self.F_empty*0
        self.failed_stars = 0
        for binary in self.population:
            create_spectrum_population += self.create_spectrum_binary(binary)
        print(len(self.population))
        print(self.failed_stars)
        return create_spectrum_population
       
    
    
def array_data(N,last_binary_line,first_binary_line):
    total_binaries = N
    mass_arr = np.zeros((total_binaries,2))
    radius_arr = np.zeros((total_binaries,2))
    L_arr = np.zeros((total_binaries,2))
    state_arr = np.zeros((total_binaries,2),dtype=object)

    for i in range(N):
        try:
            mass_arr[i][:] = np.array([ last_binary_line.S1_mass[i], last_binary_line.S2_mass[i]])
            radius_arr[i][:] = np.array([ 10**last_binary_line.S1_log_R[i], 10**last_binary_line.S2_log_R[i]])
            L_arr[i][:] = np.array([ 10**last_binary_line.S1_log_L[i], 10**last_binary_line.S2_log_L[i]])
            state_arr[i][:] = np.array([ str(last_binary_line.S1_state[i]), str(last_binary_line.S2_state[i])])
        except:
            mass_arr[i][:] = np.array(None,None)
            radius_arr[i][:] = np.array(None,None)
            L_arr[i][:] =  np.array(None,None)
            state_arr[i][:] = np.array(["Failed","Failed"])
        finally:
            if "stripped" in state_arr[i][0]:
                mass_arr[i][0] = first_binary_line.S1_mass[i]
            elif "stripped" in state_arr[i][1]:
                mass_arr[i][1] = first_binary_line.S2_mass[i]
    return mass_arr,radius_arr,L_arr,state_arr





#This might be useless for the case of using a Binary Object. 
class star():
    
    def __init__(self,binary_number,number,mass,state,R,L,metallicity):
        self.binary_number = binary_number
        self.number = number
        self.mass = mass*con.M_sun
        self.state = state
        self.R = R*con.R_sun
        self.L = L*con.L_sun
        self.logg = None
        self.Teff = None
        self.metallicity = metallicity
        self.Fe_H = np.log10(metallicity)
        
        
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

######################################################################
#Calculation of the spectra
#Creating a class that will create the grid objects and potentially calculate the fluxes as well. 
class spectral_grids():
    def __init__(self,main_grid_file="sg-CAP18-coarse.h5",secondary_grid_file="sg-BSTAR2006-medium.h5",stripped_grid_file = "sg-Gotberg18.h5",ostar_grid_file = "sg-OSTAR2002-medium.h5",lam_min=3000,lam_max=7000):
        #Assigning the grid files and creating the spec_grid ojects. 
        self.main_grid_file = main_grid_file
        self.secondary_grid_file = secondary_grid_file 
        self.stripped_grid_file = stripped_grid_file
        self.ostar_grid_file = ostar_grid_file
        self.specgrid_main = self.spec_grid(main_grid_file)
        self.specgrid_secondary = self.spec_grid(secondary_grid_file)
        self.specgrid_ostar = self.spec_grid(ostar_grid_file)
        self.specgrid_stripped = self.spec_grid(stripped_grid_file)
        self.T_max = 0  #Global_min 
        self.T_min = 0  #Global_max 
        self.logg_max = 0 
        self.logg_min = 0 
        self.lam_min = lam_min
        self.lam_max = lam_max
        self.lam = np.linspace(lam_min, lam_max, 2000)
        self.lam_c = 0.5*(self.lam[1:] + self.lam[:-1])
        self.F_empty = np.zeros(((len(self.lam_c))))*unt.m/unt.m*0

    def spec_grid(self,name):
        specgrid_file_normal = os.path.join(GRID_DIR, name)
        return pymsg.SpecGrid(specgrid_file_normal)   
    
    def global_limits(self):
        pass 
    
    #function that calculates the global min. 
    def global_limits(self):
        T_max  = 0 
        T_min = 100000
        logg_max = 0 
        logg_min = 20
        spectral_grids = [self.specgrid_main,self.specgrid_secondary,self.specgrid_ostar,self.specgrid_stripped]
        for specgrid in spectral_grids:
            for label in specgrid.axis_labels:
                if label== 'Teff':
                    T_max = max(T_max,specgrid.axis_x_max[label])
                    T_min = min(T_min,specgrid.axis_x_min[label])
                elif label == 'log(g)':
                    logg_max = max(logg_max,specgrid.axis_x_max[label])
                    logg_min = min(logg_min,specgrid.axis_x_min[label])
                
        self.T_max = T_max
        self.T_min = T_min 
        self.logg_max = logg_max
        self.logg_min = logg_min
        #print(self.T_max,self.T_min,self.logg_max,self.logg_min)
    #Function that return the flux from every grid. 
    #x is hard-coded this needs to change. 
    def main_grid_flux(self,Teff,Fe_H,logg,star):
        x = {'Teff':Teff ,'log(g)': logg,'[Fe/H]': Fe_H}
        F_not = self.specgrid_main.flux(x, self.lam)
        F_lam = np.asarray(F_not)
        return F_lam*star.R**2*kpc**-2 
    
    def secondary_grid_flux(self,Teff,logg,Z_Zo,star):
        x = {'Teff':Teff ,'log(g)': logg,'Z/Zo':Z_Zo}
        F_not = self.specgrid_secondary.flux(x, self.lam)
        F_lam = np.asarray(F_not)
        return F_lam*star.R**2*kpc**-2 

    def stripped_grid_flux(self,star):
        #star.set_metallicity(0.0142)
        Z = 0.0142
        M = star.mass/con.M_sun
        x = {'M_init':M,'Z': Z}
        F_not = self.specgrid_stripped.flux(x, self.lam)
        F_lam = np.asarray(F_not)
        return F_lam
    
    def ostar_grid_flux(self,Teff,logg,star,Z_Zo):
        x = {'Teff':Teff ,'log(g)': logg,'Z/Zo':Z_Zo}
        F_not = self.specgrid_ostar.flux(x, self.lam)
        F_lam = np.asarray(F_not)
        return F_lam*star.R**2*kpc**-2 
    
    #Maybe have a function that generally can calculate the flux. 
    def flux(self):
        pass 

#newpopulation = population_spectra(200,first_binary_line1,last_binary_line1,1,0)
#spectra1 = newpopulation.create_spectrum_population()