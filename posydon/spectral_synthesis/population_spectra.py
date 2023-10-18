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
import os
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

    def __init__(self,**kwargs):
        """Initialize a population_spectra class instance."""
        self.kwargs = default_kwargs.copy()
        for key, arg in kwargs.items():
            self.kwargs[key] = arg
        file =  self.kwargs['population_file']
        if os.path.isfile(file):
            self.file = file
        else:
            raise FileNotFoundError
        self.output_file = self.kwargs['output_file']
        self.save_data = self.kwargs['save_data']
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
        """Creates the integrated spectrum of the population.
        It also creates a file with the outputs if the save_data is True. 

        Returns:
            pop_spectrum: dictonary of type of binaries and their corresponding spectrum.
            wavelength: numpy array
        """
        scale = self.scaling_factor
        load_start = datetime.datetime.now()
        self.load_population()
        load_end = datetime.datetime.now()
        print('Loading the population took',load_end - load_start,'s')
        pop_spectrum = {}
        state_list = ['disrupted',
                      'merged', 
                      'detached',
                      'initially_single_star',
                      'low_mass_binary',
                      'contact',
                      'RLO1',
                      'RLO2']
        if self.save_data:
            labels_S1 = []
            labels_S2 = []
        # Create empty spectral arrays
        for state in state_list:
            pop_spectrum[state] = np.zeros(len(self.grids.lam_c))

        for i,binary in self.population.iterrows():
            spectrum_1,state_1,label1 = generate_spectrum(self.grids,binary,'S1',scale)
            spectrum_2,state_2,label2 = generate_spectrum(self.grids,binary,'S2',scale)
            if self.save_data:
                labels_S1.append(label1)
                labels_S2.append(label2)
            if spectrum_1 is not None and state_1 is not None:
                print(state_1, spectrum_1)
                for key,value in pop_spectrum.items():
                    print(key)
                pop_spectrum[state_1] += spectrum_1
            if spectrum_2 is not None and state_2 is not None:
                pop_spectrum[state_2] += spectrum_2
        if self.save_data:
            self.save_pop_data(self.population,labels_S1,labels_S2,pop_spectrum)
        return pop_spectrum,self.grids.lam_c


    def save_pop_data(self,pop_data,labels_S1,labels_S2,pop_spectrum,file_path=None):
        """Saves the population data and the spectrum outputs to the file 

        Args:
            pop_data: pd array
            labels_S1: string
            labels_S2: string
            file_path: string. Defaults to None.
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