######################################################################
#Calculation of the spectra
#Creating a class that will create the grid objects and potentially calculate the fluxes as well. 
import numpy as np
import os 
import astropy.units as unt
import astropy.constants as con
import pymsg

kpc = 3.08e19*unt.m
MSG_DIR = os.environ['MSG_DIR']
GRID_DIR = os.path.join(MSG_DIR, 'data', 'grids')

def spec_grid(name):
    specgrid_file_normal = os.path.join(GRID_DIR, name)
    return pymsg.SpecGrid(specgrid_file_normal)


########################################################################
class spectral_grids():
    def __init__(self,main_grid_file="sg-CAP18-coarse.h5",secondary_grid_file="sg-BSTAR2006-medium.h5",stripped_grid_file = "sg-Gotberg18.h5",ostar_grid_file = "sg-OSTAR2002-medium.h5",lam_min=3000,lam_max=7000,):
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
    