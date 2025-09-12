"""A class with functionalities that are setting up initializing the grids with MSG tools."""

__authors__ = [
    "Eirini Kasdagli <kasdaglie@ufl.edu>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
]


import os
import copy
from functools import reduce
import numpy as np
import astropy.units as unt
import pymsg
from posydon.spectral_synthesis.default_options import default_grid_kwargs
from posydon.spectral_synthesis.spectral_tools import grid_global_limits
from posydon.spectral_synthesis.libs.lib_tools import get_nearest_neighbor

kpc = 3.08e19*unt.m
MSG_DIR = os.environ['MSG_DIR']
GRID_DIR = os.path.join(MSG_DIR, 'data', 'grids')
PASS_DIR = os.path.join(MSG_DIR, 'data', 'passbands')

GRID_KEYS = [
    'main_grid',
    'secondary_grid',
    'bstar_grid',
    'ostar_grid',
    'stripped_grid',
    'WR_grid',
    'WNL_grid',
    'WNE_grid' ,
    'WC_grid'
]


def spec_grid(name):
    """Place docstring here."""
    specgrid_filepath = os.path.join(GRID_DIR, name)
    return pymsg.SpecGrid(specgrid_filepath)


########################################################################
class spectral_grids():
    """Constructs and initializes all the grids.
    """

    def __init__(self, **kwargs):
        """Class constructor."""
        # Assign the grid files and creating the spec_grid ojects.
        self.kwargs = default_grid_kwargs.copy()
        for key, arg in kwargs.items():
            self.kwargs[key] = arg

        # Specify which filters we should calculate
        self.filters = ['U', 'B', 'V']

        # Create a dictonary that stores the spec_grid objects for each key.
        self.spectral_grids = self.grid_constructor(**self.kwargs)
        self.photgrids = self.photgrid_constructor(**self.kwargs)
        cache = self.kwargs.get('cache_limit')

        # Assign cache limits:
        for key in self.spectral_grids:
            self.spectral_grids[key].cache_limit = cache

        # Load up stripped star grid
        self.stripped_grid_file = self.kwargs.get('stripped_grid')
        self.specgrid_stripped = spec_grid(self.stripped_grid_file)
        self.specgrid_stripped.cache_limit = cache

        # Getting the global limits for the grids
        self.T_max, self.T_min, self.logg_max, self.logg_min = \
            grid_global_limits(self.spectral_grids)
        #Getting the wavelength range
        self.lam,self.lam_c = self.wavelength_range(**self.kwargs)
        
    def grid_constructor(self, **kwargs):
        """Create the dictionary of MSG SpecGrid objects.

        Create the MSG SpecGrid objects from the list of files provided
        in kwargs. Load these up and store them in a dictionary.
        """
        grids = {}
        for key, arg in kwargs.items():
            print(key,arg)
            if key in GRID_KEYS:

                grids[key] = spec_grid(str(arg))
                #TODO
                #grids[key].points = find_grid_points(**kwargs)
        return grids

    def wavelength_range(self,**kwargs):
        """Generates the waveelength related variables 
        and checks if the input wavelength range is available.


        Raises:
            ValueError: If the wavelength range is not offered by the specific libraries

        Returns:
            lam: numpy array
        """
        grids = self.spectral_grids
        #The min wavelength is the the max min of all of the libraries
        lam_min_grids = reduce(lambda x,y: x if x > y else y, map(lambda x: grids[x].lam_min,grids))
        #The max wavelength is the the min max of all of the libraries
        lam_max_grids = reduce(lambda x,y: x if x < y else y, map(lambda x: grids[x].lam_max,grids))
        min_key = min(grids, key=lambda k: grids[k].lam_max)
        print(f"Grid with min lam_max: key = {min_key}, lam_max = {grids[min_key].lam_max}")
        lam_min = kwargs.get('lam_min', 3000)
        lam_max = kwargs.get('lam_max', 7000)
        lam_res = kwargs.get('lam_res', 2000)
        if lam_min < lam_min_grids or lam_max > lam_max_grids:
            raise ValueError('The wavelength range chosen is out of the',
                            ' range being offer by the libraries',
                             'The available wavelength range of this collection of libraries is [',
                             {lam_min_grids},r'$\AA$',{lam_max_grids},r']$\AA$')
        lam = np.linspace(lam_min, lam_max, lam_res)
        return lam,0.5*(lam[1:] + lam[:-1])

    def grid_flux(self, name, **kwargs):
        """Returns the flux of a star."""
        if name not in GRID_KEYS:
            raise ValueError(f'There is no grid with that name,{name}.' +
                             ' Please refer to the the grid dictonary')
        x = copy.copy(kwargs)
        specgrid = self.spectral_grids[name]
        for key in kwargs:
            if key not in specgrid.axis_labels:
                x.pop(key)
        Flux = np.asarray(specgrid.flux(x, self.lam))
        distance_factor = 1
        if name == "stripped_grid":
            distance_factor = (kpc)**2
            return Flux*distance_factor
        if name in ['WR_grid','WNE_grid','WNL_grid','WC_grid']:
            distance_factor = (kpc/100)**2
            return (Flux)*distance_factor
        #if 'log' in self.kwargs.get(name):
        #    print(name,'log')
            
        return Flux*distance_factor

    def NN_grid_flux(self, name, **kwargs):
        #We need to find the files_name:
        file_name = self.kwargs[name]
        if file_name is None:
            raise AttributeError("This grid doesn't correspond to a file name")
        x = copy.copy(kwargs)
        new_x = get_nearest_neighbor(file_name,x)
        Flux = self.grid_flux(name,**new_x)
        #Normalizing the new flux to correspond to the correct star temperature
        temp_norm_factor = (x['Teff']/new_x['Teff'])**4 
        return temp_norm_factor*Flux

    def photgrid_constructor(self, **kwargs):
        """Construct a dictionary of photogrids."""
        photgrids = {}
        for key in kwargs:
            if key in GRID_KEYS:
                photgrid = {}
                for filter in self.filters:
                    # Generate filename and check if file exists
                    filename = f'pb-Generic-Johnson.{filter}-Vega.h5'
                    passband_file_name = os.path.join(PASS_DIR, filename)
                    if not os.path.isfile(passband_file_name):
                        raise ValueError('We do not support the ' +
                                         str(filter) + ' filter.')
                    spectral_file = os.path.join(GRID_DIR, kwargs.get(key))
                    photgrid[filter] = pymsg.PhotGrid(spectral_file,
                                                      passband_file_name)
                photgrids[key] = photgrid

        return photgrids

    def photogrid_flux(self, name, scale, **kwargs):
        """Write a docstring here."""
        if name not in GRID_KEYS:
            raise ValueError('There is no grid with that name' +
                             'please refer to the the grid dictonary')
        x = copy.copy(kwargs)
        photgrid = self.photgrids[name]
        specgrid = self.spectral_grids[name]
        F = {}
        for key in kwargs:
            if key not in specgrid.axis_labels:
                x.pop(key)
        distance_factor = 1
        if name == "stripped_grid":
            distance_factor = (kpc)**2
        if name == 'WR_grid':
            distance_factor = (kpc/100)**2
        for filter in self.filters:
            F[filter] = photgrid[filter].flux(x)*distance_factor

        return F
    