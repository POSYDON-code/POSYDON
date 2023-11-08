"""
Generate the spectrum of a single star. 
"""
from copy import copy
import numpy as np
from astropy import constants as con
Lo = 3.828e33
Zo = 0.0142

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
        if x['M_init'] < grid.axis_x_min['M_init'] or x['M_init'] > grid.axis_x_max['M_init']:
            print("In here! outside limits!")
            return 'failed_grid'
        else:
            return 'stripped_grid'
    else:
        if x['Teff'] < grid.axis_x_min['Teff'] or x['Teff'] > grid.axis_x_max['Teff']:
            return 'failed_grid'
        elif x['log(g)'] < grid.axis_x_min['log(g)'] or x['log(g)'] > grid.axis_x_max['log(g)']:
            return 'failed_grid'
        else:
            return grid_name

def point_the_grid(grids,x,ostar_temp_cut_off,bstar_temp_cut_off,label,**kwargs):
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

    #First check for stripped stars because their temp can be a lot 
    # higher than the Teff of the Ostar grid limit
    if "stripped" in x['state']:
        print("Inside here!1")
        if label is not None:
            return 'failed_grid'
        return check_boundaries(grids,'stripped_grid',**x)
    
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
    #TODO include a BSTAR cutoff.

def generate_spectrum(grids,star,i,scale,**kwargs):
    """Generates the spectrum of star. 

    Args:
        grids: grid object
            Instance of the grids created in generate_spectrum
        star: tuple
        i (_type_): _description_
        scale (_type_): _description_

    Raises:
        ValueError: _description_

    Returns:
        _type_: _description_
    """
    #First we check if the star is a CO. (for future we can add WD spectra)\
    if star[f'{i}_state'] in ['massless_remnant','BH','WD','NS']:
        return None,star[f'{i}_state'],None
    ostar_temp_cut_off=28000
    bstar_temp_cut_off = 15000
    Fe_H = np.log(star['Z/Zo'])
    Z_Zo = star['Z/Zo']
    Z= star['Z/Zo']*Zo
    Teff = copy(star[f'{i}_Teff'])
    logg = copy(star[f'{i}_log_g'])
    M = copy(star[f'{i}_mass'])
    state = copy(star[f'{i}_state'])
    R = 10**copy(star[f'{i}_log_R'])*con.R_sun
    x = {'Teff':Teff ,
         'log(g)': logg,
         '[Fe/H]': Fe_H,
         'Z/Zo':Z_Zo,
         'Z':Z,
         'M_init':M,
         'state':state,
         '[alpha/Fe]':0.0}
    label = None
    label = point_the_grid(grids,x,ostar_temp_cut_off,bstar_temp_cut_off,label,**kwargs)
    count = 1
    if label == 'failed_grid':
        return None,state,label
    while count <3:
        try:
            if label == "stripped_grid":
                Flux = grids.grid_flux(label,**x)*4*np.pi*1e4/Lo
            else:
                Flux = grids.grid_flux(label,**x)*R**2*4*np.pi*1e4/Lo

            return Flux.value,star['state'],label
        except LookupError:
                label = f'failed_attempt_{count}'
        label = point_the_grid(grids,x,ostar_temp_cut_off,bstar_temp_cut_off,label,**kwargs)
        count += 1
        if label == 'failed_grid':
            return None,state,label
        else:
            raise ValueError(f'The label:{label} is not "failed_grid" after all the possible checks')