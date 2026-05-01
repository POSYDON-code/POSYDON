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
spectral_types = [
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


class population_cmd():

    def __init__(self, file, **kwargs):
        ###### Initializing the parameters ##################
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
            self.output_file = self.kwargs.get('output_file',self.file+'_cmd.h5')
            self.output_path = self.kwargs.get('output_path','./')

        # Initializing the spectral_grids object and parameters used.
        # To do put an option for changing the wavelength
        self.grids = spectral_grids()
        self.scaling_factor = kwargs.get('scaling_factor')
        self.filters = self.grids.filters
        self.photgrids = self.grids.photgrids
        self.population = None
    # Making a V-(B-V) diagram! Working on a similar logic that has been used for the spectra!

    def load_population(self):
        """Function to load up a POSYDON population."""
        metallicity = self.kwargs.get('metallicity')
        self.population = load_posydon_population(self.file,metallicity)

    def calc_colours(self):
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
            self.create_population_colours(pop)


    def create_population_colours(self,pop):
        """Creates the integrated spectrum of the population.
        It also creates a file with the outputs if the save_data is True. 

        Returns:
            pop_spectrum: dictonary of type of binaries and their corresponding spectrum.
            wavelength: numpy array
        """
        isochromes = self.kwargs.get('isochromes',False)
        spectral_type = self.kwargs.get('spectral_type',False)
        #
        if pop is None:
            pop = self.population
        
        weights = None
        if isochromes:
            mini_file = self.kwargs.get('mini_file',False)
            weights = IMF_WEIGHT(mini_file)
        
        num_waves = len(self.grids.lam_c)
        labels = []
        colours = []
        # Create empty spectral arrays
        for state in state_list:
            pop_spectrum[state] = np.zeros(num_waves)
        #Creates arrays for the spectral types as well if indicated.
        if spectral_type:
            for key in spectral_types:
                pop_spectrum[key] = np.zeros(num_waves)
        #Iterate through the whole population and calculate the spectrum of S1,S2.
        pop_dict = pop.to_dict("records")
        for i,binary in enumerate(pop_dict):
            F_obs_1,state_1,label1 = generate_colour(self.grids,binary,'S1',**self.kwargs)
            F_obs_1,state_2,label2 = generate_colour(self.grids,binary,'S2',**self.kwargs)
            
            #Store labels
            if label1 is None:
                label1 = 'failed_grid'
            if label2 is None:
                label2 = 'failed_grid'


            if F_obs_1 is not None: 
                colour_1 = colour_mag(self,F_obs_1)
                
            if F_obs_2 is not None: 
                colour_2 = colour_mag(self,F_obs_2)

            
            labels.append([label1,label2])

            colours.append([colour_1,colour_2])
            

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
        if type(pop_spectrum)== list:
            # Check if the population is empty (All of the stars are CO)
            if labels.size == 0:
                pop_data['S1_grid_status'] = [None]*len(pop_data)
                pop_data['S2_grid_status'] = [None]*len(pop_data)
                pop_data['S1_filter'] = [None]*len(pop_data)
                pop_data['S2_filter'] = [None]*len(pop_data)
            else:
                labels = np.vstack(labels)
                colours = np.vstack(colours)
                pop_data['S1_grid_status'] = labels[:,0]
                pop_data['S2_grid_status'] = labels[:,1]
                pop_data['S1_filter'] = colours[:,0]
                pop_data['S2_filter'] = colours[:,1]

            final_dict = dict(combined_spectrum)
        else:
            if labels.size == 0:
                pop_data['S1_grid_status'] = [None]*len(pop_data)
                pop_data['S2_grid_status'] = [None]*len(pop_data)
            else:
                pop_data['S1_grid_status'] = labels[:,0]
                pop_data['S2_grid_status'] = labels[:,1]
                pop_data['S1_filter'] = colours[:,0]
                pop_data['S2_filter'] = colours[:,1]


    
        h5file =self.output_path + self.output_file
        pop_data.to_hdf(h5file,key = 'data',format = 'table')



    def colour_mag(self, F_obs):
        
        mags = {}
        if F_obs is None:
            return None

        for filter in F_obs:
            mags[filter] = -2.5*np.log10(F_obs[filter])
        return mags

    def population_mag(self, num_binaries=None):
        if num_binaries is None:
            num_binaries = len(self.population)
        elif num_binaries > len(self.population):
            raise Exception(
                'The number of binaries exceeds the number of binaries in the population!')
        V = []
        B_V = []
        L = []
        population = self.population\

        for i in range(num_binaries):
            binary = population[i]
            for j in (0, 1):
                newstar = binary[j]
                magnitude = self.colour_mag(newstar)
                if magnitude is not None:
                    V.append(magnitude['V'].value)
                    B_V.append(magnitude['B'].value-magnitude['V'].value)
                    L.append(newstar.L/con.L_sun)
        return B_V, V, L

    