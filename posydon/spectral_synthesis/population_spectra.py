"""
Spectral Synthesis code
"""


__authors__ = [
    "Eirini Kasdagli <kasdaglie@ufl.edu>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
]
import h5py
import os
from copy import copy
import datetime
import numpy as np
import astropy.constants as con
import astropy.units as unt
from mpi4py import MPI
import pandas as pd
import traceback
from collections import Counter
from posydon.spectral_synthesis.spectral_tools import load_posydon_population
from posydon.spectral_synthesis.spectral_grids import spectral_grids
from posydon.spectral_synthesis.default_options import default_kwargs
from posydon.spectral_synthesis.generate_spectrum import generate_spectrum
grid_keys = [
    'main_grid',
    'secondary_grid',
    'ostar_grid',
    'stripped_grid',
    'WR_grid'
]



class population_spectra():
    """Creates and saves the output flux of a POSYDON population"""

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
        #self.total_binaries = len(self.population)
        # Initialize the spectral_grids object and parameters used.
        # TODO put an option for changing the wavelength
        self.grids = spectral_grids(**self.kwargs)
        self.scaling_factor = kwargs.get('scaling_factor')
        self.grid_flux = self.grids.grid_flux

    def load_population(self):
        """Function to load up a POSYDON population."""
        self.population = load_posydon_population(self.file)

    def create_spectrum(self):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        nprocs = comm.Get_size() 
        if rank == 0:
            load_start = datetime.datetime.now()
            self.load_population()
            load_end = datetime.datetime.now()
            print('Loading the population took',load_end - load_start,'s')
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
        pop_spectrum, labels_S1 , labels_S2 =  self.create_population_spectrum(pop)
        total_pop_spectrum = comm.gather(pop_spectrum, root=0)
        labels_S1 = comm.gather(labels_S1, root=0)
        labels_S2 = comm.gather(labels_S2, root=0)
        if self.save_data and rank ==0 :
            self.save_pop_data(self.population,np.array(labels_S1,dtype = object),np.array(labels_S2,dtype = object),total_pop_spectrum)


    def create_population_spectrum(self,pop):
        """Creates the integrated spectrum of the population.
        It also creates a file with the outputs if the save_data is True. 

        Returns:
            pop_spectrum: dictonary of type of binaries and their corresponding spectrum.
            wavelength: numpy array
        """
        if pop is None:
            scale = self.scaling_factor
            load_start = datetime.datetime.now()
            self.load_population()
            load_end = datetime.datetime.now()
            print('Loading the population took',load_end - load_start,'s')
            pop = self.population
        
        
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
           
        for i,binary in pop.iterrows():
            #print(pop.shape)
            spectrum_1,state_1,label1 = generate_spectrum(self.grids,binary,'S1',**self.kwargs)
            spectrum_2,state_2,label2 = generate_spectrum(self.grids,binary,'S2',**self.kwargs)
            if self.save_data:
                labels_S1.append(label1)
                labels_S2.append(label2)
            if spectrum_1 is not None and state_1 is not None:
                pop_spectrum[state_1] += spectrum_1
            if spectrum_2 is not None and state_2 is not None:
                pop_spectrum[state_2] += spectrum_2
            pop = pop.drop(i)
        """
        if self.save_data:
            self.save_pop_data(self.population,labels_S1,labels_S2,pop_spectrum)
            """
        return pop_spectrum,labels_S1,labels_S2
        


    def save_pop_data(self,pop_data,labels_S1,labels_S2,pop_spectrum,file_path=None):
        """Saves the population data and the spectrum outputs to the file 

        Args:
            pop_data: pd array
            labels_S1: string
            labels_S2: string
            file_path: string. Defaults to None.
        """
        if type(pop_spectrum)== list:
            """
            total_labels_S1 = []
            total_labels_S2 = []
            total_labels_S1.append(labels_S1[i] for i in range(len(labels_S1)))
            total_labels_S2.append(labels_S2[i] for i in range(len(labels_S1)))
            """
            pop_data['S1_grid_status'] = np.hstack(labels_S1)
            pop_data['S2_grid_status'] = np.hstack(labels_S2)
            
            """
            combined_spectrum = Counter(dict.fromkeys({'disrupted',
                      'merged', 
                      'detached',
                      'initially_single_star',
                      'low_mass_binary',
                      'contact',
                      'RLO1',
                      'RLO2'}, np.zeros(len(self.grids.lam_c))))
            """
            combined_spectrum = Counter(pop_spectrum[0])
            if len(pop_spectrum) > 0: 
                for i in range(1,len(pop_spectrum)):
                    pop_dict = Counter(pop_spectrum[i])
                    combined_spectrum.update(pop_dict)
                #combined_spectrum.update(pop_dict)
                """
                for key,value in combined_spectrum.items():
                    print(value)
                    print(key,pop_dict[key])
                    print(combined_spectrum[key])
                    #value += copy(pop_dict[key])
                    #combined_spectrum[key] += 
                """
            final_dict = dict(combined_spectrum)
            spectrum_data = pd.DataFrame.from_dict(final_dict)
        else:
            pop_data['S1_grid_status'] = labels_S1
            pop_data['S2_grid_status'] = labels_S2
            spectrum_data = pd.DataFrame.from_dict(pop_spectrum)

        spectrum_data.insert(loc = 0, column='wavelength',value =self.grids.lam_c )

        if file_path is None:
            file_path = "./"
        h5file =file_path + self.output_file
        pop_data.to_hdf(h5file,key = 'data',format = 'table')
        spectrum_data.to_hdf(h5file,key = 'flux',format = 'table')
