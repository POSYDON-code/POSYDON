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
            return 'failed_grid'
        elif x['log(g)'] < grids.logg_min or x['log(g)'] > grids.logg_max:
            return 'failed_grid'
        else:
            return True
    grid = grids.spectral_grids[grid_name]
    if grid_name == 'stripped_grid':
        if x['Teff'] < grid.axis_x_min['Teff'] or x['Teff'] > grid.axis_x_max['Teff']:
            return 'failed_grid'
        elif x['log(g)'] < grid.axis_x_min['log(g)'] or x['log(g)'] > grid.axis_x_max['log(g)']:
            return 'failed_grid'
        else:
            return 'stripped_grid'
    if grid_name =='WR_grid':
        if x['Teff'] < grid.axis_x_min['Teff'] or x['Teff'] > grid.axis_x_max['Teff']:
            return 'failed_grid'
        elif x['R_t'] < grid.axis_x_min['R_t'] or x['R_t'] > grid.axis_x_max['R_t']:
            return 'failed_grid'
        else: 
            return 'WR_grid'
    else:
        if x['Teff'] < grid.axis_x_min['Teff'] or x['Teff'] > grid.axis_x_max['Teff']:
            return 'failed_grid'
        elif x['log(g)'] < grid.axis_x_min['log(g)'] or x['log(g)'] > grid.axis_x_max['log(g)']:
            return 'failed_grid'
        else:
            return grid_name
    

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
        if label is not None:
            return 'failed_grid'
        return check_boundaries(grids,'stripped_grid',**x)
    
    if x['state'] == "WR_star":
        if label is not None:
            return 'failed_grid'
        return check_boundaries(grids,'WR_grid',**x)

    if isinstance(check_boundaries(grids,'global',**x),str):
        return check_boundaries(grids,'global',**x)

    #Second check for ostar stars.
    if x['Teff'] > ostar_temp_cut_off:
        if label is not None:
            return 'failed_grid'
        return check_boundaries(grids,'ostar_grid',**x)
    if x['Teff'] > bstar_temp_cut_off:
        if label is not None:
            return 'failed_grid'
        return check_boundaries(grids,'bstar_grid',**x)
    #Now we are checking at the normal grid and the number of failed trails.
    if label is None:
        check = check_boundaries(grids,'main_grid',**x)
        return check
    elif label == 'failed_attempt_1':
        return check_boundaries(grids,'secondary_grid',**x)
    elif label == 'failed_attempt_2':
        return 'failed_grid'
    else:
        raise ValueError(f'The label {label} is not recognized!')

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
    if star[f'{i}_Teff'] >= 40000:
        #Check if the stars is false labeled as H_rich 
        rename_star_state(star,i)
    
    Fe_H = np.log(star['Z/Zo'])
    Z_Zo = star['Z/Zo']
    Z= star['Z/Zo']*Zo
    Teff = copy(star[f'{i}_Teff'])
    logg = copy(star[f'{i}_log_g'])
    M_init = copy(star[f'{i}_M_init'])
    state = copy(star[f'{i}_state'])
    R = 10**copy(star[f'{i}_log_R'])*con.R_sun
    L = 10**copy(star[f'{i}_log_L'])
    surface_h1 = max(copy(star[f'{i}_surface_h1']),0.01)
    x = {'Teff':Teff ,
         'log(g)': logg,
         '[Fe/H]': Fe_H,
         'Z/Zo':Z_Zo,
         'Z':Z,
         'M_init':M_init,
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
    while count <3:
        try:
            if label == "stripped_grid":
                Flux = grids.grid_flux(label,**x)*4*np.pi*1e4/Lo
            if label == 'WR_grid':
                Flux = grids.grid_flux(label,**x)*4*np.pi*1e4/Lo *(L/10**5.3)
            else:
                Flux = grids.grid_flux(label,**x)*R**2*4*np.pi*1e4/Lo
            return Flux.value,star['state'],label
        except LookupError:
            label = f'failed_attempt_{count}'
        label = point_the_grid(grids,x,label,**kwargs)
        count += 1
        if label == 'failed_grid':
            return None,state,label
    raise ValueError(f'The label:{label} is not "failed_grid" after all the possible checks')

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
    if xH_surf < 0.4:
        if lg_M_dot < -6:
            star[f'{i}_state'] = 'stripped_He_star'
            print('here!')
        else:
           star[f'{i}_state'] = 'WR_star'
           star[f'{i}_Rt'] = calculated_Rt(star,i)
           print('WR star!')

def calculated_Rt(star,i):
    M_dot = 10**copy(star[f'{i}_lg_mdot'])
    v_terminal = 1000 #km/s
    D_max = 4 
    R  = 10**copy(star[f'{i}_log_R'])
    Rt = R*((v_terminal/2500)/(np.sqrt(D_max)*M_dot/1e-4))**(2/3)
    return Rt