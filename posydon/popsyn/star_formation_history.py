"""Implements the selection of different star-formation history scenarios."""


__authors__ = [
    'Simone Bavera <Simone.Bavera@unige.ch>',
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Devina Misra <devina.misra@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
]

import os
import numpy as np
import scipy as sp
from scipy import stats
from posydon.utils.data_download import PATH_TO_POSYDON_DATA
from posydon.utils.constants import age_of_universe
from posydon.utils.common_functions import (
    rejection_sampler, histogram_sampler, read_histogram_from_file)
from posydon.utils.constants import Zsun
from scipy.interpolate import interp1d
from astropy.cosmology import Planck15 as cosmology


SFH_SCENARIOS = ["burst", "constant", "custom_linear", "custom_log10",
                 "custom_linear_histogram", "custom_log10_histogram"]


def get_formation_times(N_binaries, star_formation='constant', **kwargs):
    """Get formation times of binaries in a population based on a SFH scenario.

    Parameters
    ----------
    N_binaries : int
        Number of formation ages to produce.
    star_formation : str, {constant, burst}
        Constant - random formation times from a uniform distribution.
        Burst - all stars are born at the same time.
    burst_time : float, 0 (years)
        Sets birth time in years.
    min_time : float, 0 (years)
        If constant SF, sets minimum of random sampling.
    max_time : float, age_of_universe (years)
        If constant SF, sets maximum of random sampling.
    RNG : <class, np.random.Generator>
        Random generator instance.

    Returns
    -------
    array
        The formation times array.
    """
    RNG = kwargs.get('RNG', np.random.default_rng())

    scenario = star_formation.lower()

    if scenario == 'burst':
        burst_time = kwargs.get('burst_time', 0.0)
        return np.ones(N_binaries)*burst_time

    max_time_default = kwargs.get('max_simulation_time', age_of_universe)
    max_time = kwargs.get('max_time', max_time_default)

    if scenario == 'constant':
        min_time = kwargs.get('min_time', 0.0)
        return RNG.uniform(size=N_binaries, low=min_time, high=max_time)

    if scenario in ["custom_linear", "custom_log10"]:
        custom_ages_file = kwargs.get('custom_ages_file')
        x, y = np.loadtxt(custom_ages_file, unpack=True)
        current_binary_ages = rejection_sampler(x, y, N_binaries)
        if "log10" in scenario:
            current_binary_ages = 10.0 ** current_binary_ages
        return max_time - current_binary_ages

    if scenario in ["custom_linear_histogram", "custom_log10_histogram"]:
        custom_ages_file = kwargs.get('custom_ages_file')
        x, y = read_histogram_from_file(custom_ages_file)
        current_binary_ages = histogram_sampler(x, y)
        if "log10" in scenario:
            current_binary_ages = 10.0 ** current_binary_ages
        return max_time - current_binary_ages

    raise ValueError(
        "Unknown star formation scenario '{}' given. Valid options: {}".
        format(star_formation, ",".join(SFH_SCENARIOS)))

def get_illustrisTNG_data(verbose=False):
    """Load IllustrisTNG SFR dataset."""
    if verbose:
        print('Loading IllustrisTNG data...')
    return np.load(os.path.join(PATH_TO_POSYDON_DATA, 'SFR/IllustrisTNG.npz'))

def star_formation_rate(SFR, z):
    """Star formation rate in M_sun yr^-1 Mpc^-3.

    Parameters
    ----------
    SFR : string
        Star formation rate assumption:
        - Madau+Fragos17 see arXiv:1606.07887
        - Madau+Dickinson14 see arXiv:1403.0007
        - Neijssel+19 see arXiv:1906.08136
        - IllustrisTNG see see arXiv:1707.03395
    z : double
        Cosmological redshift.

    Returns
    -------
    double
        The total mass of stars in M_sun formed per comoving volume Mpc^-3
        per year.
    """
    if SFR == "Madau+Fragos17":
        return 0.01 * (1. + z) ** 2.6 / (1. + ((1. + z) / 3.2) ** 6.2)  # M_sun yr^-1 Mpc^-3
    elif SFR == "Madau+Dickinson14":
        return 0.015 * (1. + z) ** 2.7 / (1. + ((1. + z) / 2.9) ** 5.6)  # M_sun yr^-1 Mpc^-3
    elif SFR == "Neijssel+19":
        return 0.01 * (1. + z) ** 2.77 / (1. + ((1. + z) / 2.9) ** 4.7) # M_sun yr^-1 Mpc^-3
    elif SFR=='IllustrisTNG':
        illustris_data = get_illustrisTNG_data()
        SFR = illustris_data['SFR'] # M_sun yr^-1 Mpc^-3
        redshifts = illustris_data['redshifts']
        SFR_interp = interp1d(redshifts, SFR)
        return SFR_interp(z)
    else:
        raise ValueError('Invalid SFR!')

def mean_metallicity(SFR, z):
    """Empiric mean metallicity function.

    Parameters
    ----------
    SFR : string
        Star formation rate assumption:
        - Madau+Fragos17 see arXiv:1606.07887
        - Madau+Dickinson14 see arXiv:1403.0007
        - Neijssel+19 see arXiv:1906.08136
    z : double
        Cosmological redshift.

    Returns
    -------
    double
        Mean metallicty of the universe at the given redhist.

    """

    if SFR == "Madau+Fragos17" or SFR == "Madau+Dickinson14":
        return 10 ** (0.153 - 0.074 * z ** 1.34) * Zsun
    elif SFR == "Neijssel+19":
        return 0.035*10**(-0.23*z)
    else:
        raise ValueError('Invalid SFR!')

def std_log_metallicity_dist(sigma):
    """Standard deviation of the log-metallicity distribution.

    Returns
    -------
    double
        Standard deviation of the adopted distribution.

    """
    if isinstance(sigma, str):
        if sigma == 'Bavera+20':
            return 0.5
        elif sigma == "Neijssel+19":
            return 0.39
        else:
            raise ValueError('Uknown sigma choice!')
    elif isinstance(sigma, float):
        return sigma
    else:
        raise ValueError(f'Invalid sigma value {sigma}!')


def SFR_Z_fraction_at_given_redshift(z,SFR, sigma, metallicity_bins, Z_max, select_one_met):
    ''''Fraction of the SFR at a given redshift z in a given metallicity bin as in Eq. (B.8) of Bavera et al. (2020).
    
    Parameters
    ----------
    z : double
        Cosmological redshift.
        
    SFR : string
        Star formation rate assumption:
        - Madau+Fragos17 see arXiv:1606.07887
        - Madau+Dickinson14 see arXiv:1403.0007
        - IllustrisTNG see see arXiv:1707.03395
        - Neijssel+19 see arXiv:1906.08136
    sigma : double / string
        Standard deviation of the log-metallicity distribution.
        If string, it can be 'Bavera+20' or 'Neijssel+19'.
    metallicity_bins : array
        Metallicity bins edges in absolute metallicity.
    Z_max : double
        Maximum metallicity in absolute metallicity.
    select_one_met : bool
        If True, the function returns the fraction of the SFR in the given metallicity bin.
        If False, the function returns the fraction of the SFR in the given metallicity bin and the fraction of the SFR in the metallicity bin
        '''
        
    if SFR == "Madau+Fragos17" or SFR=="Madau+Dickinson14":
        sigma = std_log_metallicity_dist(sigma)
        mu = np.log10(mean_metallicity(SFR, z)) - sigma**2*np.log(10)/2.
        # renormalisation constant. We can use mu[0], since we integrate over the whole metallicity range
        norm = stats.norm.cdf(np.log10(Z_max), mu[0], sigma)
        fSFR = np.empty((len(z), len(metallicity_bins)-1))
        fSFR[:,:] = np.array([(stats.norm.cdf(np.log10(metallicity_bins[1:]), m, sigma)/norm
                        - stats.norm.cdf(np.log10(metallicity_bins[:-1]), m, sigma)/norm) for m in mu])
        if not select_one_met:
            fSFR[:,0] = stats.norm.cdf(np.log10(metallicity_bins[1]), mu, sigma)/norm
            fSFR[:,-1] = norm - stats.norm.cdf(np.log(metallicity_bins[-1]), mu, sigma)/norm


    elif SFR == "Neijssel+19":
        # assume a truncated ln-normal distribution of metallicities
        sigma = std_log_metallicity_dist(sigma)
        mu = np.log(mean_metallicity(SFR, z))-sigma**2/2.
        # renormalisation constant
        norm = stats.norm.cdf(np.log(Z_max), mu[0], sigma)
        fSFR = np.empty((len(z), len(metallicity_bins)-1))
        fSFR[:,:] = np.array([(stats.norm.cdf(np.log(metallicity_bins[1:]), m, sigma)/norm
                        - stats.norm.cdf(np.log(metallicity_bins[:-1]), m, sigma)/norm) for m in mu])
        if not select_one_met:
            fSFR[:,0] = stats.norm.cdf(np.log(metallicity_bins[1]), mu, sigma)/norm
            fSFR[:,-1] = norm - stats.norm.cdf(np.log(metallicity_bins[-1]), mu, sigma)/norm


    elif SFR == 'IllustrisTNG':
        # numerically itegrate the IlluystrisTNG SFR(z,Z)
        illustris_data = get_illustrisTNG_data()
        redshifts = illustris_data['redshifts']
        Z = illustris_data['mets']
        M = illustris_data['M'] # Msun
        # only use data within the metallicity bounds (no lower bound)
        Z_max_mask =  Z <= Z_max
        redshift_indices = np.array([np.where(redshifts <= i)[0][0] for i in z])
        Z_dist = M[:,Z_max_mask][redshift_indices]

        fSFR = np.zeros((len(z), len(metallicity_bins)-1))
        for i in range(len(z)):
            if Z_dist[i].sum() == 0.:
                continue
            else:
                Z_dist_cdf = np.cumsum(Z_dist[i])/Z_dist[i].sum()
                Z_dist_cdf_interp = interp1d(np.log10(Z[Z_max_mask]), Z_dist_cdf, fill_value="extrapolate")
                fSFR[i,:] = (Z_dist_cdf_interp(np.log10(metallicity_bins[1:]))
                                    - Z_dist_cdf_interp(np.log10(metallicity_bins[:-1])))
                
                if not select_one_met:
                    #left edge
                    fSFR[i,0] = Z_dist_cdf_interp(np.log10(metallicity_bins[1]))
                    fSFR[i,-1] = 1 - Z_dist_cdf_interp(np.log10(metallicity_bins[-1]))
                
    return fSFR

def fractional_SFR_at_given_redshift(z, SFR_at_z, SFR, sigma, bins, binplace, Z_max, select_one_met=False):
    """Integrated SFR over deltaZ at a given z as in Eq. (B.9) of Bavera et al. (2020).

    Parameters
    ----------
    SFR : string
        Star formation rate assumption:
        - Madau+Fragos17 see arXiv:1606.07887
        - Madau+Dickinson14 see arXiv:1403.0007
        - IllustrisTNG see see arXiv:1707.03395
        - Neijssel+19 see arXiv:1906.08136
    Z : double
        Metallicity.

    Returns
    -------
    double
        The total mass of stars formed per comoving volume at a given redshift z
        for a given metallicity range deltaZ.

    """
    # Z_left_edge == bins[binplace-1]
    # Z_right_edge == bins[binplace]
    if SFR == "Madau+Fragos17" or SFR=="Madau+Dickinson14":
        # assume a truncated log10-normal distribution of metallicities
        # see Eq. (B.2) in Bavera et al. (2020)
        sigma = std_log_metallicity_dist(sigma)
        mu = np.log10(mean_metallicity(SFR, z)) - sigma**2*np.log(10)/2.
        # renormalisation constant
        norm = stats.norm.cdf(np.log10(Z_max), mu[0], sigma)
        fSFR = SFR_at_z*(stats.norm.cdf(np.log10(bins[binplace]), mu, sigma)/norm
                        - stats.norm.cdf(np.log10(bins[binplace-1]), mu, sigma)/norm)
        if not select_one_met:
            #left edge
            fSFR[binplace == 1] = SFR_at_z[0]*stats.norm.cdf(np.log10(bins[1]), mu[0], sigma)/norm
    elif SFR == "Neijssel+19":
        # assume a truncated ln-normal distribution of metallicities
        sigma = std_log_metallicity_dist(sigma)
        mu = np.log(mean_metallicity(SFR, z))-sigma**2/2.
        # renormalisation constant
        norm = stats.norm.cdf(np.log(Z_max), mu[0], sigma)
        fSFR = SFR_at_z*(stats.norm.cdf(np.log(bins[binplace]), mu, sigma)/norm
                        - stats.norm.cdf(np.log(bins[binplace-1]), mu, sigma)/norm)
        if not select_one_met:
            #left edge
            fSFR[binplace == 1] = SFR_at_z[0]*stats.norm.cdf(np.log(bins[1]), mu[0], sigma)/norm
    elif SFR == 'IllustrisTNG':
        # numerically itegrate the IlluystrisTNG SFR(z,Z)
        illustris_data = get_illustrisTNG_data()
        redshifts = illustris_data['redshifts']
        Z = illustris_data['mets']
        M = illustris_data['M'] # Msun
        # only use data within the metallicity bounds (no lower bound)
        valid_met_idxs = np.where(Z <= Z_max)[0]
        # get the index of the correct redshift in the data
        redz_idx = np.where(redshifts <= z[0])[0][0]
        # take values of the data at this redshift between our metallicity bounds
        Z_dist = M[redz_idx, valid_met_idxs]
        if Z_dist.sum() == 0.:
            fSFR = np.zeros(len(z))
        else:
            Z_dist_cdf = np.cumsum(Z_dist)/Z_dist.sum()
            Z_dist_cdf_interp = interp1d(np.log10(Z[valid_met_idxs]), Z_dist_cdf, fill_value="extrapolate")
            fSFR = SFR_at_z*(Z_dist_cdf_interp(np.log10(bins[binplace]))
                            - Z_dist_cdf_interp(np.log10(bins[binplace-1])))
            if not select_one_met:
                #left edge
                fSFR[binplace == 1] = SFR_at_z[0]*Z_dist_cdf_interp(np.log10(bins[1]))
    else:
        raise ValueError('Invalid SFR!')

    return fSFR


def integrated_SFRH_over_redshift(SFR, sigma, Z, Z_max):
    """Integrated SFR history over z as in Eq. (B.10) of Bavera et al. (2020).

    Parameters
    ----------
    SFR : string
        Star formation rate assumption:
        - Madau+Fragos17 see arXiv:1606.07887
        - Madau+Dickinson14 see arXiv:1403.0007
        - Neijssel+19 see arXiv:1906.08136
    Z : double
        Metallicity.

    Returns
    -------
    double
        The total mass of stars formed per comoving volume at a given
        metallicity Z.

    """
    def E(z, Omega_m=cosmology.Om0):
        Omega_L = 1.- Omega_m
        return (Omega_m*(1.+z)**3 + Omega_L)**(1./2.)
    def f(z,Z):
        if SFR == "Madau+Fragos17" or SFR == "Madau+Dickinson14":
            sigma = std_log_metallicity_dist(sigma)
            mu = np.log10(mean_metallicity(SFR, z))-sigma**2*np.log(10)/2.
            H_0 = cosmology.H0.to('1/yr').value # yr
            # put a cutoff on metallicity at Z_max
            norm = stats.norm.cdf(np.log10(Z_max), mu, sigma)
            return star_formation_rate(SFR, z)*stats.norm.pdf(np.log10(Z), mu, sigma)/norm*(H_0*(1.+z)*E(z))**(-1)
        elif SFR == "Neijssel+19":
            sigma = std_log_metallicity_dist(sigma)
            mu = np.log10(mean_metallicity(SFR, z))-sigma**2/2.
            H_0 = cosmology.H0.to('1/yr').value # yr
            return star_formation_rate(SFR, z)*stats.norm.pdf(np.log(Z), mu, sigma)*(H_0*(1.+z)*E(z))**(-1)
        else:
            raise ValueError('Invalid SFR!')

    return sp.integrate.quad(f, 1e-10, np.inf, args=(Z,))[0] # M_sun yr^-1 Mpc^-3
