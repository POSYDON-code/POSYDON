import numpy as np
import posydon.utils.constants as c
import matplotlib.pyplot as plt
import sys
import os

#Connect to MSG path.



MSG_DIR = os.environ['MSG_DIR']

sys.path.insert(0, os.path.join(MSG_DIR, 'python'))
import pymsg

# Load the SpecGrid:

GRID_DIR = os.path.join(MSG_DIR, 'data', 'grids')

specgrid_file_name = os.path.join(GRID_DIR, 'sg-demo.h5')

specgrid = pymsg.SpecGrid(specgrid_file_name)


#Create binary spectra.

def create_spectrum_binary(binary, **kwargs):

    # check if binary has two non-degenerate stars
    if binary.star_1.state not in ['NS', 'BH']:
        spectrum_1 = create_spectrum_single(binary.star_1, kwargs)

    if binary.star_2.state in ['NS', 'BH']:
        spectrum_2 = create_spectrum_single(binary.star_2, kwargs)


    if spectrum_1 is None and spectrum_2 is None:
        return None
    elif spectrum_1 is None:
        return spectrum_2
    elif spectrum_2 is None:
        return spectrum_1
    else:
        return add_spectra(spectrum_1, spectrum_2)


def create_spectrum_single(star, **kwargs):

    mass = star.mass
    radius = 10**star.log_R
    luminosity = 10**star.log_L

    # Check units here!!
    log_g = c.standard_cgrav*mass/radius**2

    # Calculate temperature from luminosity and radius
    teff = (luminosity / (4*np.pi) / c.boltz_sigma / radius**2)**(1/4)

    # Define your variables:
    x = {'Teff': Teff,'log(g)': logg}
    # Set up the wavelength abscissa

    lam_min = 3000.
    lam_max = 7000.

    lam = np.linspace(lam_min, lam_max, 501)

    lam_c = 0.5*(lam[1:] + lam[:-1])

    #


    # Call MSG to obtain spectrum


    return spectrum



def add_spectrum(spectra_1, spectra_2):

    # Combine spectra

    return combined_spectrum
