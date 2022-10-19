import numpy as np
import pymsg
import posydon.utils.constants as c


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

    # Call MSG to obtain spectrum


    return spectrum



def add_spectrum(spectra_1, spectra_2):

    # Combine spectra

    return combined_spectrum
