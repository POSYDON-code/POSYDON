# pylint: disable=import-error
"""
Spectral Synthesis code
"""


__authors__ = [
    "Eirini Kasdagli <kasdaglie@ufl.edu>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
]
import datetime
import os
from collections import Counter
from copy import copy

import numpy as np
import pandas as pd
from mpi4py import MPI

from posydon.spectral_synthesis.default_options import default_kwargs
from posydon.spectral_synthesis.generate_spectrum import generate_spectrum
from posydon.spectral_synthesis.spectral_grids import spectral_grids
from posydon.spectral_synthesis.spectral_tools import (
    IMF_WEIGHT,
    load_posydon_population,
)

grid_keys = [
    'main_grid',
    'secondary_grid',
    'ostar_grid',
    'stripped_grid',
    'WR_grid',
    'WNL_grid',
    'WNE_grid',
    'WC_grid'
]
state_list = [
    'disrupted',
    'merged',
    'detached',
    'initially_single_star',
    'low_mass_binary',
    'contact',
    'RLO1',
    'RLO2'
    ]
spectral_types_list = [
    'main_grid',
    'bstar_grid',
    'failed_grid',
    'no_grid',
    'secondary_grid',
    'stripped_grid',
    'WR_grid',
    'ostar_grid',
    'WNL_grid',
    'WNE_grid',
    'WC_grid'  ]

class population_spectra():
    """Creates and saves the output flux of a POSYDON population"""

    def __init__(self,**kwargs):
        """Initialize a population_spectra class instance."""
        self.kwargs = default_kwargs.copy()
        for key, arg in kwargs.items():
            self.kwargs[key] = arg

        self.verbose = self.kwargs.get('verbose', False)

        file =  self.kwargs['population_file']
        if os.path.isfile(file):
            self.file = file
        else:
            raise FileNotFoundError

        self.save_data = self.kwargs['save_data']
        if self.save_data:
            self.output_file = self.kwargs.get('output_file',self.file+'_spectra.h5')
            self.output_path = self.kwargs.get('output_path','./')

        if self.verbose:
            print("Initializing population spectra with the following options:")
            print(f"  population_file : {self.file}")
            if self.save_data:
                print(f"  output_file     : {self.output_file}")
                print(f"  output_path     : {self.output_path}")
            print(f"  save_data       : {self.save_data}")

        # Initialize the spectral_grids object and parameters used.
        self.grids = spectral_grids(**self.kwargs)
        self.population = None

    def load_population(self):
        """Function to load up a POSYDON population."""
        metallicity = self.kwargs.get('metallicity')
        self.population = load_posydon_population(self.file,metallicity)

    def create_spectrum(self):
        """ It splits up the population and combines the data to be saved.
        """
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        nprocs = comm.Get_size()
        if rank == 0:
            load_start = datetime.datetime.now()
            self.load_population()
            load_end = datetime.datetime.now()
            if self.verbose:
                print(f"  Population loaded: {len(self.population)} binaries")
                print(f"  Loading took     : {load_end - load_start}")
            pop = copy(self.population)
            # determine the size of each sub-task
            ave, res = divmod(len(pop), nprocs)
            counts = [ave + 1 if p < res else ave for p in range(nprocs)]
            # determine the starting and ending indices of each sub-task
            starts = [sum(counts[:p]) for p in range(nprocs)]
            ends = [sum(counts[:p+1]) for p in range(nprocs)]
            # converts data into a list of arrays
            pop = [pop[starts[p]:ends[p]] for p in range(nprocs)]
        else:
            pop = None
        pop = comm.scatter(pop, root=0)
        if self.save_data:
            pop_spectrum, labels =  self.create_population_spectrum(pop)
            total_pop_spectrum = comm.gather(pop_spectrum, root=0)
            labels = comm.gather(labels, root=0)
            if rank ==0 :
                self.save_pop_data(self.population,np.array(labels,dtype = object),total_pop_spectrum)
        else:
            self.create_population_spectrum(pop)

    def create_population_spectrum(self,pop):
        """Creates the integrated spectrum of the population.
        It also creates a file with the outputs if the save_data is True.

        Returns:
            pop_spectrum: dictonary of type of binaries and their corresponding spectrum.
            wavelength: numpy array
        """
        isochrones = self.kwargs.get('isochrones',False)
        #
        if pop is None:
            pop = self.population

        weights = None
        if isochrones:
            mini_file = self.kwargs.get('mini_file',False)
            weights = IMF_WEIGHT(mini_file)

        num_waves = len(self.grids.lam_c)
        pop_spectrum = {}
        labels = []
        # Create empty spectral arrays


        #Setting the specific fluxes to be saved.
        include_states = self.kwargs.get('include_states', None)
        include_spectral_types = self.kwargs.get('include_spectral_types', None)

        # Validate user-provided filters against known lists
        if include_states is not None:
            invalid = set(include_states) - set(state_list)
            if invalid:
                raise ValueError(f"Unknown states in 'include_states': {invalid}. "
                                f"Valid options are: {state_list}")

        if include_spectral_types is not None:
            invalid = set(include_spectral_types) - set(spectral_types_list)
            if invalid:
                raise ValueError(f"Unknown spectral types in 'include_spectral_types': {invalid}. "
                                f"Valid options are: {spectral_types_list}")


        states = include_states if include_states is not None else state_list
        spectral_types = include_spectral_types if include_spectral_types is not None else spectral_types_list

        #Initialize the flux arrays
        for key in spectral_types:
            pop_spectrum[key] = np.zeros(num_waves)
        for state in states:
            pop_spectrum[state] = np.zeros(num_waves)

        pop_spectrum['Total'] = np.zeros(num_waves)

        #Iterate through the whole population and calculate the spectrum of S1,S2.
        pop_dict = pop.to_dict("records")
        for i,binary in enumerate(pop_dict):
            spectrum_1,state_1,label1 = generate_spectrum(self.grids,binary,'S1',**self.kwargs)
            spectrum_2,state_2,label2 = generate_spectrum(self.grids,binary,'S2',**self.kwargs)

            #Store labels
            if label1 is None:
                label1 = 'failed_grid'
            if label2 is None:
                label2 = 'failed_grid'

            labels.append([label1,label2])
            if spectrum_1 is not None and state_1 is not None:
                if isochrones:
                    spectrum_1 = spectrum_1*weights[i]
                pop_spectrum['Total'] += spectrum_1
                if state_1 in states:
                    pop_spectrum[state_1] += spectrum_1
                if label1 in spectral_types:
                    pop_spectrum[label1] += spectrum_1
            if spectrum_2 is not None and state_2 is not None:
                pop_spectrum['Total'] += spectrum_2
                if state_2 in states:
                    pop_spectrum[state_2] += spectrum_2
                if label2 in spectral_types:
                    pop_spectrum[label2] += spectrum_2

        if self.save_data:
            return pop_spectrum,labels
        return pop_spectrum

    def save_pop_data(self,pop_data,labels,pop_spectrum):
        """Saves the population data and the spectrum outputs to the file
        Args:
            pop_data: pd array
            labels_S1: string
            labels_S2: string
            file_path: string. Defaults to None.
        """
        save_population_data = self.kwargs.get('save_population_data', True)

        if type(pop_spectrum)== list:
            # Check if the population is empty (All of the stars are CO)
            if labels.size == 0:
                pop_data['S1_grid_status'] = [None]*len(pop_data)
                pop_data['S2_grid_status'] = [None]*len(pop_data)
            else:
                labels = np.vstack(labels)
                pop_data['S1_grid_status'] = labels[:,0]
                pop_data['S2_grid_status'] = labels[:,1]
            combined_spectrum = Counter(pop_spectrum[0])
            #combined_spectrum = {}
            #for key in pop_spectrum[0]:
            #    combined_spectrum[key] = sum(d.get(key, 0.0) for d in pop_spectrum)
            if len(pop_spectrum) > 0:
                for i in range(1,len(pop_spectrum)):
                    pop_dict = Counter(pop_spectrum[i])
                    combined_spectrum.update(pop_dict)
                #combined_spectrum.update(pop_dict)
            final_dict = dict(combined_spectrum)
            spectrum_data = pd.DataFrame.from_dict(final_dict)
        else:
            if labels.size == 0:
                pop_data['S1_grid_status'] = [None]*len(pop_data)
                pop_data['S2_grid_status'] = [None]*len(pop_data)
            else:
                pop_data['S1_grid_status'] = labels[:,0]
                pop_data['S2_grid_status'] = labels[:,1]
            spectrum_data = pd.DataFrame.from_dict(pop_spectrum)

        spectrum_data.insert(loc = 0, column='wavelength',value =self.grids.lam_c )
        h5file = os.path.join(self.output_path,self.output_file)

        if save_population_data:
            pop_data.to_hdf(h5file,key = 'data',format = 'table')
        spectrum_data.to_hdf(h5file,key = 'flux',format = 'table')

        if self.verbose:
            print(f"  Data saved successfully.")
