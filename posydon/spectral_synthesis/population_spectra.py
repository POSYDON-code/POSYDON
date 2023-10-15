"""
Spectral Synthesis code
"""
#################################### Imports ##############################################

import numpy as np
import astropy.constants as con
import astropy.units as unt
import datetime
import pandas as pd
import traceback
from copy import copy

from posydon.spectral_synthesis.spectral_tools import population_data,load_posydon_population
from posydon.spectral_synthesis.spectral_grids import spectral_grids
from posydon.spectral_synthesis.default_options import default_kwargs
from posydon.spectral_synthesis.generate_spectrum import generate_spectrum
grid_keys = [
    'main_grid',
    'secondary_grid',
    'ostar_grid',
    'stripped_grid',
]

######### Creating the class population_spectra ###########################


class population_spectra():
    """Write a class docstring."""

    def __init__(self,file,output_file,**kwargs):
        """Initialize a population_spectra class instance."""
        self.kwargs = default_kwargs.copy()
        for key, arg in kwargs.items():
            self.kwargs[key] = arg

        # population_file = kwargs.get('population_file')
        self.kwargs['population_file'] = file
        self.kwargs['output_file'] = output_file
        self.save_data = self.kwargs['save_data']
        #time = kwargs.get('time')
        self.failed_stars = 0  # int, the stars that failed durin the spectra process
        #Creating lists/numpy arrays for investigation for the failled and succesful stars.
        self.stars_fails = []
        self.stars_run = []
        self.stars_grid_fail = []
        self.file = file
        self.output_file = output_file
        # Create readable arrays for the stars objects.
        time_start_pop = datetime.datetime.now()
        self.population = population_data(**self.kwargs)
        time_end_pop = datetime.datetime.now()
        print('Total time is: ', time_end_pop - time_start_pop)
        self.total_binaries = len(self.population)
        # Initialize the spectral_grids object and parameters used.
        # TODO put an option for changing the wavelength
        self.grids = spectral_grids()
        self.scaling_factor = kwargs.get('scaling_factor')
        self.grid_flux = self.grids.grid_flux

    def load_population(self):
        """Function to load up a POSYDON population."""
        self.population = load_posydon_population(self.file)

        

    def create_population_spectrum(self):
        """_summary_

        Returns:
            _type_: _description_
        """
        scale = self.scaling_factor
        load_start = datetime.datetime.now()
        self.load_population()
        load_end = datetime.datetime.now()
        print('Loading the population took',load_end - load_start,'s')
        pop_spectrum = {}

        state_list = ['disrupted', 'merged', 'detached','initially_single_star','low_mass_binary','contact','RLO1','RLO2']
        if self.save_data:
            labels_S1 = []
            labels_S2 = []
        # Create empty spectral arrays
        for state in state_list:
            pop_spectrum[state] = np.zeros(len(self.grids.lam_c))       
        for i,binary in self.population.iterrows():
            #TODO write line bellow

            spectrum_1,state_1,label1 = generate_spectrum(self.grids,binary,'S1',scale)
            spectrum_2,state_2,label2 = generate_spectrum(self.grids,binary,'S2',scale)
            if self.save_data:
                labels_S1.append(label1)
                labels_S2.append(label2)
            if spectrum_1 is not None and state_1 is not None:
                pop_spectrum[state_1] += spectrum_1
            if spectrum_2 is not None and state_2 is not None:
                pop_spectrum[state_2] += spectrum_2
        if self.save_data:
            self.save_pop_data(self.population,labels_S1,labels_S2,pop_spectrum)
        return pop_spectrum,self.grids.lam_c


    def save_pop_data(self,pop_data,labels_S1,labels_S2,pop_spectrum,file_path=None):
        """_summary_

        Args:
            pop_data (_type_): _description_
            labels_S1 (_type_): _description_
            labels_S2 (_type_): _description_
            file_path (_type_, optional): _description_. Defaults to None.
        """
        pop_data['S1_grid_status'] = labels_S1
        pop_data['S2_grid_status'] = labels_S2
        spectrum_data = pd.DataFrame.from_dict(pop_spectrum)
        spectrum_data.insert(loc = 0, column='wavelength',value =self.grids.lam_c )
        if file_path is None:
            file_path = "./"
        h5file =file_path + self.output_file
        pop_data.to_hdf(h5file,key = 'data',format = 'table')
        spectrum_data.to_hdf(h5file,key = 'flux',format = 'table')


    
    def create_spectrum_single(self,star,ostar_temp_cut_off=27000,**kwargs):
        scale = self.scaling_factor

        if "stripped" in star.state:
            M = star.mass/con.M_sun
            M_min =  self.grids.spectral_grids['stripped_grid'].axis_x_min['M_init']
            M_max = self.grids.spectral_grids['stripped_grid'].axis_x_max['M_init']
            if M < M_min or M > M_max:
                self.failed_stars +=1
                return None,None
            return self.grids.stripped_grid_flux(star)

        Fe_H = star.Fe_H
        Z_Zo = star.metallicity
        Teff = copy(star.get_Teff(self.grids.T_max,self.grids.T_min))
        logg = copy(star.get_logg(self.grids.logg_max,self.grids.logg_min))
        x = {'Teff':Teff ,'log(g)': logg,'[Fe/H]': Fe_H,'Z/Zo':Z_Zo,'[alpha/Fe]':0.0}

        if Teff == None or logg == None:
            self.failed_stars +=1
            return None,None

        #For temperature higher than the desired cut_off we calculate the spectra using the Ostar girds.
        if Teff >= ostar_temp_cut_off:
            #Setting the acceptable boundaries for the ostar grid in the logg.
            logg_min = self.grids.spectral_grids['ostar_grid'].axis_x_min['log(g)']
            logg_max = self.grids.spectral_grids['ostar_grid'].axis_x_max['log(g)']
            if logg > logg_min and logg<logg_max:
                    try:
                        Flux = self.grid_flux('ostar_grid',**x)
                        return Flux*star.R**2*scale**-2
                    except LookupError:
                        self.failed_stars +=1
                        return None
            else:
                self.failed_stars +=1
                return None
        try:
            Flux = self.grid_flux('main_grid',**x)
            return Flux*star.R**2*scale**-2
        #TODO
        #except LookupError:
        except:
            try:
                # lo limits for the secondary grid:
                #logg_min = grids.specgrid_secondary.axis_x_min['log(g)']
                #Teff_min = grids.specgrid_secondary.axis_x_min['Teff']
                #Calculating all the exeption that the first grid gave with the secondary.
                #if Teff>15000.0 and logg > logg_min:
                Flux = self.grid_flux('secondary_grid',**x)
                return Flux*star.R**2*scale**-2
            #except LookupError:
            except:
                # if and else statements that fix the grid voids.
                if Teff > 20000:
                    logg = max(logg, 4.0)
                elif Teff > 12000:
                    logg = max(logg, 3.0)
                elif Teff > 8000:
                    logg = max(logg, 2.0)
                elif Teff > 6000:
                    logg = max(logg, 1.0)
                try:
                    Flux = self.grid_flux('main_grid',**x)
                    return Flux*star.R**2*scale**-2
                #except LookupError:
                except:
                    self.failed_stars +=1
                    return None



    def add_spectra(self,spectrum_1,spectrum_2):
        if spectrum_1 is not None and spectrum_2 is not None:
            return spectrum_1 + spectrum_2
        elif spectrum_1 is None:
            return spectrum_2
        elif spectrum_2 is None:
            return spectrum_1
        else:
            raise Exception("Spectrum 1 and Spectrum 2 have not allowed values.")

    def create_spectrum_binary(self,binary):
    # check if binary has two or any non-degenerate stars
        star1 = binary[0]
        star2 = binary[1]
        if star1.state not in ['NS','BH','massless_remnant','WD']:
            spectrum_1 = self.create_spectrum_single(star1)
        else:
            spectrum_1 = None

        if star2.state not in ['NS','BH','massless_remnant','WD']:
            spectrum_2 = self.create_spectrum_single(star2)
        else:
            spectrum_2 = None
        return self.add_spectra(spectrum_1, spectrum_2)


    def create_spectrum_population(self,num_binaries = None ):
        """_summary_

        Args:
            num_binaries (_type_, optional): _description_. Defaults to None.

        Raises:
            Exception: _description_

        Returns:
            _type_: _description_
        """
        if num_binaries is None:
            num_binaries = len(self.population)
        elif num_binaries > len(self.population):
            raise Exception('The number of binaries exceeds the number of binaries in the population!')

        create_spectrum_population = np.zeros(len(self.grids.lam_c))*unt.m/unt.m*0
        failed_binaries = 0
        population = self.population
        for i in range(num_binaries):
            binary = population[i]
            spectrum = self.create_spectrum_binary(binary)
            if spectrum is not None:
                create_spectrum_population += spectrum
        return create_spectrum_population , self.grids.lam_c


    def create_spectrum_subpopulation(self,state,num_binaries = None):
        if num_binaries is None:
            num_binaries = len(self.population)
        elif num_binaries > len(self.population):
            raise Exception('The number of binaries exceeds the number of binaries in the population!')
        failed_binaries = 0
        create_spectrum_population = np.zeros(len(self.grids.lam_c))*unt.m/unt.m*0
        population = self.population

        for i in range(num_binaries):
            binary = population[i]
            star1 = binary[0]
            if star1.binary_state == state:
                spectrum = self.create_spectrum_binary(binary)
                if spectrum is not None:
                    self.create_spectrum_population += spectrum
        return create_spectrum_population , self.grids.lam_c
