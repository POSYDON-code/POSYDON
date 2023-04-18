"""
Spectral Synthesis code 


"""


#################################### Imports ##############################################
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
#Connect to MSG path.
#os.environ['MSG_DIR'] = '/home/kasdaglie/blue/kasdaglie/msg-1.1.2'
MSG_DIR = os.environ['MSG_DIR']

sys.path.insert(0, os.path.join(MSG_DIR, 'python'))
import pymsg

# Load the SpecGrid:

GRID_DIR = os.path.join(MSG_DIR, 'data', 'grids')

specgrid_file_name = os.path.join(GRID_DIR, 'sg-demo.h5')

specgrid = pymsg.SpecGrid(specgrid_file_name)
%matplotlib inline
plt.rcParams.update({'font.size': 16})





#########################################Creating the class population_spectra ###############################################


class population_spectra():
    
    def __init__(self,total_binaries,first_binary_line,last_binary_line,data_file=None):
        ######Initializing the parameters ##################
        self.total_binaries = total_binaries
        self.data_file = data_file
        self.first_binary_line = first_binary_line
        self.last_binary_line = last_binary_line
        self.population = []
        
        mass,radius,L,state = array_data(self.total_binaries,self.last_binary_line,self.first_binary_line)
        
        for i in range(self.total_binaries):
            star1 = newstar(i,0,mass[i,0],state[i,0],radius[i,0],L[i,0])
            star2 = newstar(i,1,mass[i,1],state[i,1],radius[i,1],L[i,1])
            self.population.append((star1,star2))
        
        
    def create_spectrum_population(self):
        create_spectrum_population = F_empty*0
        failed_binaries = 0
        for binary in self.population:
            create_spectrum_population += self.create_spectrum_binary(binary)
        print(failedstar)
        return create_spectrum_population
    
    
    def create_spectrum_binary(self,binary):
    # check if binary has two non-degenerate stars
        star1 = binary[0]
        star2 = binary[1]
        if star1.state not in ['NS','BH','Failed']:
            spectrum_1 = new_create_spectrum_single(star1,star1.binary_number)
        else:
            spectrum_1 = F_empty*0

        if star2.state not in ['NS','BH','Failed']:
            spectrum_2 = new_create_spectrum_single(star2,star2.binary_number)
        else: 
            spectrum_2 =F_empty*0

        return new_add_spectra(spectrum_1, spectrum_2)

    
    
    
    def run_spectra():
        pass 



