######################################################################
#Calculation of the spectra
#Creating a class that will create the grid objects and potentially calculate the fluxes as well.
import numpy as np
import os
import astropy.units as unt
import astropy.constants as con
import pymsg
import datetime
from functools import reduce 
from posydon.spectral_synthesis.default_options import default_grid_kwargs
from posydon.spectral_synthesis.spectral_tools import grid_global_limits
import copy
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
]


def spec_grid(name):
    """Place docstring here."""
    specgrid_filepath = os.path.join(GRID_DIR, name)
    return pymsg.SpecGrid(specgrid_filepath)


########################################################################
class spectral_grids():
    """Write a docstring here."""

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
        self.lam_c = self.wavelength_range(**self.kwargs)

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
        return grids

    def wavelength_range(self,**kwargs):
        grids = self.spectral_grids
        #The min wavelength is the the max min of all of the libraries
        lam_min_grids = reduce(lambda x,y: x if x > y else y, map(lambda x: grids[x].lam_min,grids))
        #The max wavelength is the the min max of all of the libraries
        lam_max_grids = reduce(lambda x,y: x if x < y else y, map(lambda x: grids[x].lam_max,grids))
        lam_min = kwargs.get('lam_min', 3000)
        lam_max = kwargs.get('lam_max', 7000)
        lam_res = kwargs.get('lam_res', 2000)
        if lam_min < lam_min_grids or lam_max > lam_max_grids: 
            raise ValueError(f'The wavelength range chosen is out of the range being offer by the libraries',
                             'The available wavelength range of this collection of libraries is [{lam_min}$\AA$,{lam_max}]$\AA$')
        lam = np.linspace(lam_min, lam_max, lam_res)
        return 0.5*(lam[1:] + lam[:-1])
    
    def grid_flux(self, name, **kwargs):
        """Returns the flux of a star."""
        if name not in GRID_KEYS:
            raise ValueError('There is no grid with that name' +
                             'please refer to the the grid dictonary')
        x = copy.copy(kwargs)
        specgrid = self.spectral_grids[name]
        for key, arg in kwargs.items():
            if key not in specgrid.axis_labels:
                x.pop(key)
        Flux = np.asarray(specgrid.flux(x, self.lam))
        return Flux

    def stripped_grid_flux(self, star):
        """Return the flux for a stripped star."""
        Z = 0.0142
        M = star.mass/con.M_sun
        x = {'M_init': M, 'Z': Z}
        Flux = np.asarray(self.specgrid_stripped.flux(x, self.lam))
        return Flux

    def photgrid_constructor(self, **kwargs):
        """Construct a dictionary of photogrids."""
        photgrids = {}
        for key, arg in kwargs.items():
            if key in GRID_KEYS:
                photgrid = {}
                for filter in self.filters:
                    # Generate filename and check if file exists
                    filename = f'pb-Generic-Johnson.{filter}-Vega.h5'
                    passband_file_name = os.path.join(PASS_DIR, filename)
                    if not os.path.isfile(passband_file_name):
                        raise ValueError('We do not support the ' +
                                         str(filter) + ' filter.')
                        continue

                    
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
        for key, arg in kwargs.items():
            if key not in specgrid.axis_labels:
                x.pop(key)
        for filter in self.filters:
            F[filter] = photgrid[filter].flux(x)*scale
        return F

####TO DO:
# Write a function that will make sure the wavelength coverage doesn't exceed the library's coverage.
# Set a default switch for C3K and CAP depending on the wavelength
