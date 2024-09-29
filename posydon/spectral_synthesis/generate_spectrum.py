"""Generates the spectrum of a single star 

These functions are determining, based on the star's properties which is the appropiate 
grid to be used for calculated the flux.
If not such a grid is found or the interpolation fails, it returns a failed_grid message.
"""


__authors__ = [
    "Eirini Kasdagli <kasdaglie@ufl.edu>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
]

from copy import copy
import numpy as np
from astropy import constants as con
import pandas as pd
from posydon.spectral_synthesis.spectral_tools import smooth_flux_negatives
from posydon.spectral_synthesis.spectral_grids import GRID_KEYS 
#Constants
Lo = 3.828e33 #Solar Luminosity erg/s
Zo = 0.0142 #Solar metallicity

def check_boundaries(grids,grid_name,**kwargs):
    """Checks if the stellar parameters in inside the boundaries of the spectral grids.

    It takes an input based on the grid check need to be made and returns a failed_grid message 
    if parameters are outside the boundaries. 
    """
    x = copy(kwargs)
    #First we check the global limits
    if grid_name == "global":
        if x['Teff'] < grids.T_min or x['Teff'] > grids.T_max:
            return 'failed_grid',x
        elif x['log(g)'] < grids.logg_min or x['log(g)'] > grids.logg_max:
            return 'failed_grid',x
        else:
            return True
    grid = grids.spectral_grids[grid_name]
    axis_labels = grid.axis_labels
    for axis_label in axis_labels:
        if x[axis_label] < grid.axis_x_min[axis_label]:
            if abs(x[axis_label] - grid.axis_x_min[axis_label])/(grid.axis_x_min[axis_label]) < 0.1:
                x[axis_label] = grid.axis_x_min[axis_label]
            else:
                return 'failed_grid',x 
        elif x[axis_label] > grid.axis_x_max[axis_label]:
            if abs(x[axis_label] - grid.axis_x_max[axis_label])/(grid.axis_x_max[axis_label]) < 0.1:
                x[axis_label] = grid.axis_x_max[axis_label]
            else: 
                return 'failed_grid',x 
    #The star is within the grid limits.
    return grid_name,x 

def point_the_grid(grids,x,label,**kwargs):
    """Assigning the write label that would point to the spectra grid needed to used

    Args:
        grids: grid object
            Instance of the grids created in generate_spectrum
        x: dictionary
            Includes all the stellar parameters that are needed
        ostar_temp_cut_off: float
            Temperature at which we change to the OSTAR grid
        label: string
            Points to the correct grid needed to be used for the star 
            or None if the star hasn't been assessed before. 

    Raises:
        ValueError: If a label is not recognized and it's not one of the allowed options.

    Returns:
        string: Returns the name of the grid needed to be used. 
        Failed grid if it can't be matched to a grid.
    """
    ostar_temp_cut_off = kwargs.get('ostar_temp_cut_off',28000)
    bstar_temp_cut_off = kwargs.get('bstar_temp_cut_off',15000)
    #First check for stripped stars because their temp can be a lot
    # higher than the Teff of the Ostar grid limit
    if "stripped" in x['state']:
        if (label is None):
            new_label,x = check_boundaries(grids,'stripped_grid',**x)
            if new_label == 'failed_grid':
                if x['Teff'] > ostar_temp_cut_off:
                    return check_boundaries(grids,'ostar_grid',**x)
                elif x['Teff'] > bstar_temp_cut_off:
                    return check_boundaries(grids,'bstar_grid',**x)
            else:
                return new_label,x
        if label ==  'failed_attempt_1':
            return check_boundaries(grids,'stripped_grid',**x)

    if x['state'] == "WNE_star":
        if label is None or label == 'failed_attempt_1':
            new_label,x = check_boundaries(grids,'WNE_grid',**x)
            if new_label == 'failed_grid':
                new_label,x =  check_boundaries(grids,'stripped_grid',**x)
                if new_label == 'failed_grid':
                    if x['Teff'] > ostar_temp_cut_off:
                        return check_boundaries(grids,'ostar_grid',**x)
                    elif x['Teff'] > bstar_temp_cut_off:
                        return check_boundaries(grids,'bstar_grid',**x)
                    else:
                        return new_label,x
                else:
                    return new_label,x
            else:
                return new_label,x
    
    if x['state'] == "WNL_star":
        if label is None or label == 'failed_attempt_1':
            new_label,x = check_boundaries(grids,'WNL_grid',**x)
            if new_label == 'failed_grid':
                new_label,x =  check_boundaries(grids,'stripped_grid',**x)
                if new_label == 'failed_grid':
                    if x['Teff'] > ostar_temp_cut_off:
                        return check_boundaries(grids,'ostar_grid',**x)
                    elif x['Teff'] > bstar_temp_cut_off:
                        return check_boundaries(grids,'bstar_grid',**x)
                    else:
                        return new_label,x
                else:
                    return new_label,x
            else:
                return new_label,x
    #if isinstance(check_boundaries(grids,'global',**x),str):
    #    return check_boundaries(grids,'global',**x)

    #Second check for ostar stars.
    if x['Teff'] > ostar_temp_cut_off:
        if label is not None:
            if (label == 'failed_attempt_1'):
                return check_boundaries(grids,'ostar_grid',**x)
            return 'failed_grid',x
        return check_boundaries(grids,'ostar_grid',**x)
    if x['Teff'] > bstar_temp_cut_off:
        if label is not None or (label == 'failed_attempt_1') :
            return 'failed_grid',x
        return check_boundaries(grids,'bstar_grid',**x)
    #Now we are checking at the normal grid and the number of failed trails.
    if label is None or label == 'failed_grid':
        return check_boundaries(grids,'main_grid',**x)
    elif label == 'failed_attempt_1':
        return check_boundaries(grids,'secondary_grid',**x)
    elif label == 'failed_attempt_2':
        return 'failed_grid',x
    else:
        raise ValueError(f'The label {label} is not recognized!. The state of the star is {x}')

def generate_spectrum(grids,star,i,**kwargs):
    """Generates the spectrum of star. 

    Args:
        grids: grid object
            Instance of the grids created in generate_spectrum
        star: tuple
        i: string
            S1 or S2 

    Raises:
        ValueError: If the final label is not "failed_grid" or not recognised

    Returns:
        Flax: float 
        state: string
        label: string 
    """
    #First we check if the star is a CO. (for future we can add WD spectra)\
    if star[f'{i}_state'] in ['massless_remnant','BH','WD','NS']:
        return None,star[f'{i}_state'],None
    if star[f'{i}_surface_h1'] <= 0.6:
        #Check if the stars is false labeled as H_rich 
        rename_star_state(star,i)
    
    Fe_H = np.log(star['Z/Zo'])
    Z_Zo = star['Z/Zo']
    #Z= star['Z/Zo']*Zo
    Teff = copy(star[f'{i}_Teff'])
    logg = copy(star[f'{i}_log_g'])
    state = copy(star[f'{i}_state'])
    R = 10**copy(star[f'{i}_log_R'])*con.R_sun
    L = 10**copy(star[f'{i}_log_L'])
    surface_h1 = max(copy(star[f'{i}_surface_h1']),0.01)
    x = {'Teff':Teff ,
         'log(g)': logg,
         '[Fe/H]': Fe_H,
         'Z/Zo':Z_Zo,
         'state':state,
         'surface_h1' : surface_h1,
         '[alpha/Fe]':0.0}
    if state in ['WR_star','WNE_star','WNL_star']:
        x['R_t'] = star[f'{i}_Rt']
    label = None
    label,x = point_the_grid(grids,x,label,**kwargs)
    count = 1
    if label == 'failed_grid':
        return None,state,label
    while count <3:
        if 'failed_attempt' not in label:
            try:
                if label == "stripped_grid":
                    Flux = grids.grid_flux(label,**x)*4*np.pi*1e4/Lo
                elif label in ['WR_grid','WNE_grid','WNL_grid']:
                    Flux = grids.grid_flux(label,**x)*4*np.pi*1e4/Lo *(L/10**5.3)
                    #Replace the negative values for WR
                    Flux.value[Flux.value < 0] = 1e-99
                else:
                    Flux = grids.grid_flux(label,**x)*R**2*4*np.pi*1e4/Lo
                if np.min(Flux) < 0: 
                    Flux = smooth_flux_negatives(grids.lam_c,Flux.value)
                    return Flux,star['state'],label
                return Flux.value,star['state'],label
            except LookupError:
                try:
                    x = rescale_log_g(grids,label,**x)
                except Exception as e:
                    print('Under the exception',e)
                label = f'failed_attempt_{count}'
        label,x = point_the_grid(grids,x,label,**kwargs)
        count += 1
        if label == 'failed_grid':
            return None,state,label
    raise ValueError(f'The label:{label} is not "failed_grid" after all the possible checks. The star is {x}')

def regenerate_spectrum(grids,star,i,**kwargs):
    label = star[f'{i}_grid_status']
    if  pd.isna(label):
        return None,star[f'{i}_state'],label
    if label == "failed_grid":
        return None,star[f'{i}_state'],label
    Fe_H = np.log(star['Z/Zo'])
    Z_Zo = star['Z/Zo']
    #Z= star['Z/Zo']*Zo
    Teff = copy(star[f'{i}_Teff'])
    logg = copy(star[f'{i}_log_g'])
    state = copy(star[f'{i}_state'])
    R = 10**copy(star[f'{i}_log_R'])*con.R_sun
    L = 10**copy(star[f'{i}_log_L'])
    surface_h1 = max(copy(star[f'{i}_surface_h1']),0.01)
  
    x = {'Teff':Teff ,
         'log(g)': logg,
         '[Fe/H]': Fe_H,
         'Z/Zo':Z_Zo,
         'state':state,
         'surface_h1' : surface_h1,
         '[alpha/Fe]':0.0}
    #If the star is WR we need to add R_t quantity in x 

    if label == "stripped_grid":
        Flux = grids.grid_flux(label,**x)*4*np.pi*1e4/Lo
    elif label == 'WR_grid':
        if (x['Teff'] > 100000) & (((x['Teff'] -100000)/100000) < 0.1):
            x['Teff'] = 100000
        elif (x['Teff'] > 100000)  & (((x['Teff'] -100000)/100000) > 0.1): 
            return None,star[f'{i}_state'],label
        x['R_t'] = calculated_Rt(star,i)
        Flux = grids.grid_flux(label,**x)*4*np.pi*1e4/Lo *(L/10**5.3)
    else:
        Flux = grids.grid_flux(label,**x)*R**2*4*np.pi*1e4/Lo
    return Flux.value,star['state'],label

def rename_star_state(star,i):
    """Rename the star states of WR stars and 
    stripped He-stars that look like a H-rich star in the POSYDON_data

    Args:
        star: object

    """
    xH_surf = copy(star[f'{i}_surface_h1'])
    T = copy(star[f'{i}_Teff'])
    lg_M_dot = copy(star[f'{i}_lg_mdot'])
    k_e = 0.2*(1 + xH_surf)
    logg = copy(star[f'{i}_log_g'])
    Gamma = k_e * con.sigma_sb*T**4/(con.c * 10**logg)
    if lg_M_dot < -6:
        star[f'{i}_state'] = 'stripped_He_star'
    else:
        star[f'{i}_Rt'] = calculated_Rt(star,i)
        if xH_surf < 0.1: 
            star[f'{i}_state'] = 'WNE_star'
        else:
            star[f'{i}_state'] = 'WNL_star' 
        

def calculated_Rt(star,i):
    M_dot = 10**copy(star[f'{i}_lg_mdot'])
    v_terminal = 1000 #km/s
    D_max = 4 
    R  = 10**copy(star[f'{i}_log_R'])
    Rt = R*((v_terminal/2500)/(np.sqrt(D_max)*M_dot/1e-4))**(2/3)
    return Rt

def find_nearest_neighbor_stripped(**x):
    #We want a function that will find the nearest neighbor in case of a star falling
    #in the grey arays of the grids. We only do nn only in when the star has a temperature that is within the 
    #range of the grid.
    possible_loggs = np.array([4.0,4.3,4.5,4.8,5.0,5.2,5.5,5.7,6.0])  
    logg = copy(x['log(g)'])
    new_logg = possible_loggs[possible_loggs > logg][0]
    x['log(g)'] = new_logg
    return x

def rescale_log_g(grids,label,**x):
    dx = {}
    old_x = copy(x)
    if label not in GRID_KEYS:
        raise ValueError(f"The label {label} doesn't correspond to any grid")
    grid = grids.spectral_grids[label]

    for key in x:
        if key not in grid.axis_labels:
            old_x.pop(key)

    for axis_label in grid.axis_labels:
        dx[axis_label] = 0.0
    if label == 'WR_grid':
        dx['R_t'] =  grid.axis_x_max['R_t']
        new_x = grid.adjust_x(old_x, dx)
        x['R_t'] = new_x['R_t']
    else:
        dx['log(g)'] = grid.axis_x_max['log(g)']
        new_x = grid.adjust_x(old_x, dx)
        x['log(g)'] = new_x['log(g)']
    return x 

def close_to_boundaries(label,values,**x):
    for value in values:   
        if abs(x[label] - value)/(value) < 0.1:
            x[label] = value
    return x[label]

def generate_photgrid_flux(grids,star,i,**kwargs):
    """Generates the spectrum of star. 

    Args:
        grids: grid object
            Instance of the grids created in generate_spectrum
        star: tuple
        i: string
            S1 or S2 

    Raises:
        ValueError: If the final label is not "failed_grid" or not recognised

    Returns:
        Flax: float 
        state: string
        label: string 
    """
    #First we check if the star is a CO. (for future we can add WD spectra)\
    if star[f'{i}_state'] in ['massless_remnant','BH','WD','NS']:
        return None,star[f'{i}_state'],None
    if star[f'{i}_Teff'] >= 40000:
        #Check if the stars is false labeled as H_rich 
        rename_star_state(star,i)
    
    Fe_H = np.log(star['Z/Zo'])
    Z_Zo = star['Z/Zo']
    #Z= star['Z/Zo']*Zo
    Teff = copy(star[f'{i}_Teff'])
    logg = copy(star[f'{i}_log_g'])
    state = copy(star[f'{i}_state'])
    R = 10**copy(star[f'{i}_log_R'])*con.R_sun
    L = 10**copy(star[f'{i}_log_L'])
    surface_h1 = max(copy(star[f'{i}_surface_h1']),0.01)
    x = {'Teff':Teff ,
         'log(g)': logg,
         '[Fe/H]': Fe_H,
         'Z/Zo':Z_Zo,
         'state':state,
         'surface_h1' : surface_h1,
         '[alpha/Fe]':0.0}
    if state == 'WR_star':
        x['R_t'] = star[f'{i}_Rt']
    label = None
    label = point_the_grid(grids,x,label,**kwargs)
    count = 1
    if label == 'failed_grid':
        return None,state,label
    while count <= 4:
        try:
            if label == "stripped_grid":
                Flux = grids.photogrid_flux(label,**x)
            elif label == 'WR_grid':
                Flux = grids.grid_flux(label,**x)(L/10**5.3)
            else:
                Flux = grids.grid_flux(label,**x)
            return Flux.value,star['state'],label
        except LookupError:
            label = f'failed_attempt_{count}'
        label = point_the_grid(grids,x,label,**kwargs)
        count += 1
        if label == 'failed_grid':
            return None,state,label
    raise ValueError(f'The label:{label} is not "failed_grid" after all the possible checks')
