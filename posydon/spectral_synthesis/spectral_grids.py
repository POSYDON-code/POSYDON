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
PASS_DIR = os.path.join(MSG_DIR, 'data', 'passbands')

grid_keys = [
    'main_grid',
    'secondary_grid',
    'ostar_grid',
    'stripped_grid',
]



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
        
        self.filters = None
        #Create a dictonary named grids that has stored the spec_grid objects for each keys. 
        self.spectral_grids = self.grid_constractor(**self.kwargs)
        self.photgrids = self.photgrid_constractor(**self.kwargs)
        cache = self.kwargs.get('cache_limit')
        
        #assigning cache limits:
        for key in self.spectral_grids:
            self.spectral_grids[key].cache_limit = cache

        self.stripped_grid_file = self.kwargs.get('stripped_grid_file')
        self.specgrid_stripped = spec_grid(self.stripped_grid_file)
        self.specgrid_stripped.cache_limit = cache
        
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

        for key in self.spectral_grids:
            specgrid = self.spectral_grids[key]
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
        
    
    #Function that creates the specgrid objects of all the spectral libraries and then storing them in a list. 
    def grid_constractor(self,**kwargs):
        grids = {}
        for key in grid_keys:
            for file in kwargs:
                if key in str(file):
                    grids[key] = spec_grid(kwargs.get(file))
        return grids
    
    def grid_flux(self,name,**kwargs):
        if name not in grid_keys:
            raise Exception('There is no grid with that name, please refer to the the grid dictonary')
        x = copy.copy(kwargs)
        specgrid = self.spectral_grids[name]
        for key, arg in kwargs.items():
            if key not in specgrid.axis_labels:
                x.pop(key)
        Flux = np.asarray(specgrid.flux(x, self.lam)) 
        return Flux  

    #Function that return the flux from every grid. 

    def stripped_grid_flux(self,star):
        Z = 0.0142
        M = star.mass/con.M_sun
        x = {'M_init':M,'Z': Z}
        Flux = np.asarray(self.specgrid_stripped.flux(x, self.lam))
        return Flux
    

    #Constractor that makes a dictionary of photogrids. 
    def photgrid_constractor(self,**kwargs):
        if self.filters is None: 
            filters = ['U', 'B', 'V']
        else: 
            filters = self.filters
        photgrids = {}
        for key in grid_keys:
            for file in kwargs:
                if key in str(file):
                    photgrid= {}
                    for filter in filters: 
                        passband_file_name = os.path.join(PASS_DIR, f'pb-Generic-Johnson.{filter}-Vega.h5')
                        spectral_file = os.path.join(GRID_DIR, kwargs.get(file))
                        photgrid[filter] = pymsg.PhotGrid(spectral_file,passband_file_name)
                    photgrids[key]= photgrid
        return photgrids

    
    def photogrid_flux(self,name,scale,**kwargs):
        if self.filters is None: 
            filters = ['U', 'B', 'V']
        else: 
            filters = self.filters

        if name not in grid_keys:
            raise Exception('There is no grid with that name, please refer to the the grid dictonary')
        x = copy.copy(kwargs)
        photgrid = self.photgrids[name]
        specgrid = self.spectral_grids[name]
        F ={}
        for key, arg in kwargs.items():
            if key not in specgrid.axis_labels:
                x.pop(key)
        for filter in filters:
           F[filter] = photgrid[filter].flux(x)*scale
        return F

####TO DO:
# Write a function that will make sure the wavelength coverage doesn't exceed the library's coverage. 
# Set a default switch for C3K and CAP depending on the wavelength 