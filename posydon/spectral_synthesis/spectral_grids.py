######################################################################
#Calculation of the spectra
#Creating a class that will create the grid objects and potentially calculate the fluxes as well. 
import numpy as np
import os 
import astropy.units as unt
import astropy.constants as con
import pymsg
import datetime
from posydon.spectral_synthesis.default_options import default_kwargs
import copy
kpc = 3.08e19*unt.m
MSG_DIR = os.environ['MSG_DIR']
GRID_DIR = os.path.join(MSG_DIR, 'data', 'grids')

def spec_grid(name):
    specgrid_file_normal = os.path.join(GRID_DIR, name)
    return pymsg.SpecGrid(specgrid_file_normal)


########################################################################
class spectral_grids():
    def __init__(self,**kwargs):
        #Assigning the grid files and creating the spec_grid ojects. 
        self.kwargs = default_kwargs.copy()
        for key, arg in kwargs.items():
            self.kwargs[key] = arg
        
        main_grid_file = self.kwargs.get('main_grid_file')
        self.specgrid_main = spec_grid(main_grid_file)

        secondary_grid_file = self.kwargs.get('secondary_grid_file')
        self.specgrid_secondary = spec_grid(secondary_grid_file)

        ostar_grid_file = self.kwargs.get('ostar_grid_file')
        self.specgrid_ostar = spec_grid(ostar_grid_file)

        stripped_grid_file = self.kwargs.get('stripped_grid_file')
        self.specgrid_stripped = spec_grid(stripped_grid_file)
        self.T_max = 0  #Global_min 
        self.T_min = 0  #Global_max 
        self.logg_max = 0 
        self.logg_min = 0 

        lam_min = self.kwargs.get('lam_min')
        lam_max = self.kwargs.get('lam_max')
        self.lam = np.linspace(lam_min,lam_max, 2000)
        self.lam_c = 0.5*(self.lam[1:] + self.lam[:-1])
        self.F_empty = np.zeros(((len(self.lam_c))))*unt.m/unt.m*0
    
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
        

    def main_grid_flux(self,**kwargs):
        x = copy.copy(kwargs)
        for key, arg in kwargs.items():
            if key not in self.specgrid_main.axis_labels:
                x.pop(key)
        Flux = np.asarray(self.specgrid_main.flux(x, self.lam)) 
        return Flux
    
    def secondary_grid_flux(self,**kwargs):
        x = copy.copy(kwargs)
        for key, arg in kwargs.items():
            if key not in self.specgrid_secondary.axis_labels:
                x.pop(key)
        Flux = np.asarray(self.specgrid_secondary.flux(x, self.lam))
        return Flux

    def stripped_grid_flux(self,star):
        Z = 0.0142
        M = star.mass/con.M_sun
        x = {'M_init':M,'Z': Z}
        Flux = np.asarray(self.specgrid_stripped.flux(x, self.lam))
        return Flux
    
    def ostar_grid_flux(self,**kwargs):
        x = copy.copy(kwargs)
        for key, arg in kwargs.items():
            if key not in self.specgrid_ostar.axis_labels:
                x.pop(key)
        Flux = np.asarray(self.specgrid_ostar.flux(x, self.lam))
        return Flux
    