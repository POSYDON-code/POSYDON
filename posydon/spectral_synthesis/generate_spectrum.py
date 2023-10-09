"""
Generate the spectrum of a single star. 
"""
import sys
import posydon
import os
import numpy as np
import astropy.constants as con
import astropy.units as unt
import matplotlib.pyplot as plt
import h5py
from copy import copy
from astropy.io import fits
import pandas as pd
import math
import datetime
import matplotlib.patches as mpatches
import traceback

def check_boundaries(grids,grid_name,**kwargs):
    """Checks if the stellar parameters in inside the boundaries of the spectral grids.

    It takes an input based on the grid check need to be made and returns a failled_grid message 
    if parameters are outside the boundaries. 
    """
    x = copy(kwargs)
    'First we check the global limits'
    if grid_name == "global":
        if x['Teff'] < grids.T_min or x['Teff'] > grids.T_max:
            return 'failled_grid'
        elif x['log(g)'] < grids.logg_min or x['log(g)'] > grids.logg_max:
            return 'failled_grid'
        else: 
            return True
    elif grid_name == 'stripped_grid': 
        if x['M'] < grids.spectral_grids[grid_name].axis_x_min['M_init'] or x['M'] > grids.spectral_grids[grid_name].axis_x_max['M_init']:
            return 'failled_grid'
        else: 
            return 'stripped_grid'
    else: 
        if x['Teff'] < grids.spectral_grids[grid_name].axis_x_min['Teff'] or x['Teff'] > grids.spectral_grids[grid_name].axis_x_max['Teff']:
            return 'failled_grid'
        elif x['log(g)'] < grids.spectral_grids[grid_name].axis_x_min['log(g)'] or x['log(g)'] > grids.spectral_grids[grid_name].axis_x_max['log(g)']:
            return 'failled_grid'
        else: 
            return grid_name

def point_the_grid(grids,x,ostar_temp_cut_off,label,**kwargs):
    if isinstance(check_boundaries(grids,'global',**x),str):
        return check_boundaries(grids,'global',**x)
    
    #First check for stripped stars 
    if "stripped" in x['state']:
        return check_boundaries(grids,'stripped_grid',**x)
    #Second check for ostar stars.
    if x['Teff'] > ostar_temp_cut_off:
        if label is not None: 
            return 'failed_grid'
        return check_boundaries(grids,'ostar_grid',**x)
    
    #Now we are checking at the normal grid and the number of failed trails.
    if label is None:
        check = check_boundaries(grids,'main_grid',**x)
        return check
    elif label == 'failled_attempt_1':
        return check_boundaries(grids,'secondary_grid',**x)
    elif label == 'failled_attempt_2':
        return 'failed_grid'
    else:
        raise ValueError(f'The label {label} is not recognized!')
    #TODO include a BSTAR cutoff. 
    

def generate_spectrum(grids,star,i,scale,**kwargs):
    #First we check if the star is a CO. (for future we can add WD spectra)\
    if star[f'{i}_state'] in ['massless_remnant','BH','WD','NS']:
         return None,star[f'{i}_state']
    ostar_temp_cut_off=27000
    Fe_H = np.log(star['Z/Zo'])
    Z_Zo = star['Z/Zo']
    Teff = copy(star[f'{i}_Teff'])
    logg = copy(star[f'{i}_log_g'])
    M = copy(star[f'{i}_mass'])
    state = copy(star[f'{i}_state'])
    R = 10**copy(star[f'{i}_log_R'])
    #TODO have a star label that it's going to be None in the begining"
    #label = star[f'{i}_label']
    x = {'Teff':Teff ,'log(g)': logg,'[Fe/H]': Fe_H,'Z/Zo':Z_Zo,'M_init':M,'state':state} 
    label = None
    label = point_the_grid(grids,x,ostar_temp_cut_off,label,**kwargs)
    count = 1
    while count <3:
        if label == 'failled_grid':
            return None,state
        else: 
            try:
                Flux = grids.grid_flux(label,**x)*R**2*scale**-2
                print(type(Flux))
                return Flux.value,star['state']
            except LookupError:
                label = f'failled_attempt_{count}'
        label = point_the_grid(grids,x,ostar_temp_cut_off,label,**kwargs)
        count =+ 1 
    if label == 'failed_grid':
        return None,state
    else: 
        raise ValueError(f'The label:{label} is not "failed_grid" after all the possible checks')
    #TODO add this as in the star attributes

