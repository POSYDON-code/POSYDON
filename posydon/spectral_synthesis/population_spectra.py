# pylint: disable=import-error
"""
Spectral Synthesis code
"""


__authors__ = [
    "Eirini Kasdagli <kasdaglie@ufl.edu>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
]
import os
from copy import copy
import datetime
import numpy as np
from mpi4py import MPI
import pandas as pd
from collections import Counter
from posydon.spectral_synthesis.spectral_tools import load_posydon_population,IMF_WEIGHT
from posydon.spectral_synthesis.spectral_grids import spectral_grids
from posydon.spectral_synthesis.default_options import default_kwargs
from posydon.spectral_synthesis.generate_spectrum import generate_spectrum,regenerate_spectrum
grid_keys = [
    'main_grid',
    'secondary_grid',
    'ostar_grid',
    'stripped_grid',
    'WR_grid'
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
spectral_types = [
    'main_grid', 
    'bstar_grid', 
    'failed_grid',
    np.nan,
    'stripped_grid',
    'WR_grid',
    'ostar_grid']

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
        
        self.save_data = self.kwargs['save_data']
        if self.save_data:
            self.output_file = self.kwargs.get('output_file',self.file+'_spectra.h5')
            self.output_path = self.kwargs.get('output_path','./')
        # Initialize the spectral_grids object and parameters used.
        self.grids = spectral_grids(**self.kwargs)
        self.population = None

    def load_population(self):
        """Function to load up a POSYDON population."""
        self.population = load_posydon_population(self.file)

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
        isochromes = self.kwargs.get('isochromes',False)
        spectral_type = self.kwargs.get('spectral_type',False)
        if pop is None:
            pop = self.population
        if isochromes:
            mini_file = self.kwargs.get('mini_file',False)
            weights = IMF_WEIGHT(mini_file)
        pop_spectrum = {}
        labels = []
        # Create empty spectral arrays
        for state in state_list:
            pop_spectrum[state] = np.zeros(len(self.grids.lam_c))
        #Creates arrays for the spectral types as well if indicated.
        if spectral_type:
            for key in spectral_types:
                pop_spectrum[key] = np.zeros(len(self.grids.lam_c))
        #Iterate through the whole population and calculate the spectrum of S1,S2.
        for i,binary in pop.iterrows():
            spectrum_1,state_1,label1 = generate_spectrum(self.grids,binary,'S1',**self.kwargs)
            spectrum_2,state_2,label2 = generate_spectrum(self.grids,binary,'S2',**self.kwargs)

            if spectrum_1 is not None and state_1 is not None:
                if isochromes:
                    spectrum_1 = spectrum_1*weights[i]
                pop_spectrum[state_1] += spectrum_1
                if spectral_type:
                    pop_spectrum[label1] += spectrum_1
            if spectrum_2 is not None and state_2 is not None:
                pop_spectrum[state_2] += spectrum_2
                if spectral_type:
                    pop_spectrum[label2] += spectrum_2
            pop = pop.drop(i)
            labels.append([label1,label2])

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
        file_path = self.output_file

        if type(pop_spectrum)== list:
            pop_data['S1_grid_status'] = np.hstack(labels)[:,0]
            pop_data['S2_grid_status'] = np.hstack(labels)[:,1]
            combined_spectrum = Counter(pop_spectrum[0])
            if len(pop_spectrum) > 0: 
                for i in range(1,len(pop_spectrum)):
                    pop_dict = Counter(pop_spectrum[i])
                    combined_spectrum.update(pop_dict)
                #combined_spectrum.update(pop_dict)
            final_dict = dict(combined_spectrum)
            spectrum_data = pd.DataFrame.from_dict(final_dict)
        else:
            pop_data['S1_grid_status'] = labels[:,0]
            pop_data['S2_grid_status'] = labels[:,1]
            spectrum_data = pd.DataFrame.from_dict(pop_spectrum)

        spectrum_data.insert(loc = 0, column='wavelength',value =self.grids.lam_c )

        if file_path is None:
            file_path = "./"
        h5file =file_path + self.output_file
        pop_data.to_hdf(h5file,key = 'data',format = 'table')
        spectrum_data.to_hdf(h5file,key = 'flux',format = 'table')


class  sub_population_spectra(population_spectra):
    """
    Creates the integrated spectrum of an output file.
    It also creates a file with the outputs if the save_data is True. 

    Returns:
    pop_spectrum: dictonary of type of binaries and their corresponding spectrum.
    wavelength: numpy array
    """

    def load_population(self):
        """Function to load up a spectral output population."""
        self.population = pd.read_hdf(self.file, key = 'data')
    
    def create_spectrum(self):
        print('before MPI')
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
            super().save_pop_data(self.population,np.array(labels_S1,dtype = object),np.array(labels_S2,dtype = object),total_pop_spectrum)

    def create_population_spectrum(self,pop,spectral_type = True):
        """Creates the integrated spectrum of the population.
        It also creates a file with the outputs if the save_data is True. 

        Returns:
            pop_spectrum: dictonary of type of binaries and their corresponding spectrum.
            wavelength: numpy array
        """
        if pop is None:
            load_start = datetime.datetime.now()
            self.load_population()
            load_end = datetime.datetime.now()
            print('Loading the population took',load_end - load_start,'s')
            pop = self.population     
        pop_spectrum = {}
        if self.save_data:
            labels_S1 = []
            labels_S2 = []
        # Create empty spectral arrays

        if spectral_type:
            for key in spectral_types:
                pop_spectrum[key] = np.zeros(len(self.grids.lam_c))
        else:
            for state in state_list:
                pop_spectrum[state] = np.zeros(len(self.grids.lam_c))

        for i,binary in pop.iterrows():
            spectrum_1,state_1,label1 = regenerate_spectrum(self.grids,binary,'S1',**self.kwargs)
            spectrum_2,state_2,label2 = regenerate_spectrum(self.grids,binary,'S2',**self.kwargs)
            if self.save_data:
                labels_S1.append(label1)
                labels_S2.append(label2)
            if spectrum_1 is not None and state_1 is not None:
                if spectral_type:
                    pop_spectrum[label1] += spectrum_1
                else:
                    pop_spectrum[state_1] += spectrum_1
            if spectrum_2 is not None and state_2 is not None:
                if spectral_type:
                    pop_spectrum[label2] += spectrum_2
                else:
                    pop_spectrum[state_2] += spectrum_2
               
            pop = pop.drop(i)
        return pop_spectrum,labels_S1,labels_S2

    
