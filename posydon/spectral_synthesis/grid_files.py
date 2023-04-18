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
GRID_DIR = os.path.join(MSG_DIR, 'data', 'grids')
###########################################################################################

#name is string: 
#The available grids so far are: 
# 'sg-CAP18-coarse.h5' , 'sg-BSTAR2006-medium.h5','sg-Gotberg18.h5','sg-OSTAR2002-medium.h5'

def specgrid(name):
    specgrid_file_normal = os.path.join(GRID_DIR, name)
    return specgrid = pymsg.SpecGrid(specgrid_file_normal)


