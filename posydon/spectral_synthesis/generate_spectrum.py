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
from astropy.io import fits
import pandas as pd
import math
import datetime
import matplotlib.patches as mpatches
import traceback

def check_boundaries(grids,**kwargs):
    x = copy(kwargs)
    'First we check the global limits'
    if x['Teff'] < grids.T_min or x['Teff'] > grids.T_max:
        return 'failled'
    if x['logg'] < grids.logg_min or x['logg'] > grids.logg_max:
        return 'failled'
def point_the_grid(grids,x,ostar_temp_cut_off,**kwargs):
    #Returns the gird label that needs to be used. 
    
    if "stripped" in x['state']:
        M_min =  grids.spectral_grids['stripped_grid'].axis_x_min['M_init']
        M_max = grids.spectral_grids['stripped_grid'].axis_x_max['M_init']
        if M < M_min or M > M_max:
            label = 'stripped_grid'
        else:
            label = 'failled_grid'
        return label
    if x['Teff'] > ostar_temp_cut_off:
        logg_min = grids.spectral_grids['ostar_grid'].axis_x_min['log(g)']
        logg_max = grids.spectral_grids['ostar_grid'].axis_x_max['log(g)']
        if x['logg'] > logg_min and x['logg']<logg_max:
            label = 'ostar_grid'
        else:
            label = 'failled_grid'
    

def generate_spectrum(grids,star,i,scale,**kwargs):
    #First we check if the star is a CO. (for future we can add WD spectra)
    if star['{x}_state'] in ['massless_remnant','BH','WD','NS']:
         return None
    ostar_temp_cut_off=27000
    Fe_H = np.log(star['Z/Zo'])
    Z_Zo = star['Z/Zo']
    Teff = copy(star[f'{i}_Teff'])
    logg = copy(star[f'{i}_log_g'])
    M = star[f'{i}_mass']
    state = star[f'{i}_state']
    #label = star[f'{i}_label']
    x = {'Teff':Teff ,'log(g)': logg,'[Fe/H]': Fe_H,'Z/Zo':Z_Zo,'M_init':M,'state':state}
    star[f'{i}_label'] = None 
    #TODO create a function that checks the global boundaries as well as the boundaries of the grids. 
    #check_boundaries
    label = point_the_grid(grids,x,ostar_temp_cut_off,label,**kwargs)
    while count <= 3:
        if label == 'failled_grid':
            return None 
        else: 
            try:
                Flux = grids.grid_flux(label,**x)
                return Flux*star.R**2*scale**-2,star['state']
            except LookupError:
                label = f'failled_attempt_{count}'
        label = point_the_grid(grids,x,ostar_temp_cut_off,label,**kwargs)
        count =+ 1 

    #TODO Create a function call check_bounds()
    #This function will check if the star is or isn't inside the Teff and logg boundaries of the library collection

        #For temperature higher than the desired cut_off we calculate the spectra using the Ostar girds.
            #Setting the acceptable boundaries for the ostar grid in the logg.

