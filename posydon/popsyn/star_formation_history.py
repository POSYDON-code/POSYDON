"""Implements the selection of different star-formation history scenarios."""


__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Devina Misra <devina.misra@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Max Briel <max.briel@unige.ch>",
]

import os
import numpy as np
import scipy as sp
from scipy import stats
from posydon.config import PATH_TO_POSYDON_DATA
from posydon.utils.constants import age_of_universe
from posydon.utils.common_functions import (
    rejection_sampler,
    histogram_sampler,
    read_histogram_from_file,
)
from posydon.utils.constants import Zsun
from posydon.utils.interpolators import interp1d
from astropy.cosmology import Planck15 as cosmology
from abc import ABC, abstractmethod

SFH_SCENARIOS = [
    "burst",
    "constant",
    "custom_linear",
    "custom_log10",
    "custom_linear_histogram",
    "custom_log10_histogram",
]


class SFHBase(ABC):
    '''Abstract class for star formation history models'''    
    def __init__(self, MODEL):
        self.MODEL = MODEL
        # Automatically attach all model parameters as attributes
        for key, value in MODEL.items():
            setattr(self, key, value)

    @abstractmethod
    def CSFRD(self, z):
        """Compute the cosmic star formation rate density."""
        pass

    @abstractmethod
    def mean_metallicity(self, z):
        """Return the mean metallicity at redshift z."""
        pass
        
    @abstractmethod
    def fSFR(self, z, metallicity_bins):
        """Return the fractional SFR as a function of redshift and metallicity bins."""
        pass

    def std_log_metallicity_dist(self):
        sigma = self.sigma
        if isinstance(sigma, str):
            if sigma == "Bavera+20":
                return 0.5
            elif sigma == "Neijssel+19":
                return 0.39
            else:
                raise ValueError("Unknown sigma choice!")
        elif isinstance(sigma, float):
            return sigma
        else:
            raise ValueError(f"Invalid sigma value {sigma}!")

    def __call__(self, z, met_bins):
        '''Return the star formation history at a given redshift and metallicity bins
        
        Parameters
        ----------
        z : float or array-like
            Cosmological redshift.
        met_bins : array
            Metallicity bins edges in absolute metallicity.
        
        Returns
        -------
        array
            Star formation history per metallicity bin at the given redshift(s).
        '''
        return self.CSFRD(z)[:, np.newaxis] * self.fSFR(z, met_bins)
    

class MadauBase(SFHBase):
    """
    Base class for Madau-style star-formation history implementations.
    This class implements common methods for CSFRD, mean metallicity,
    and fractional SFR based on the chosen Madau parameterisation.
    The specific parameters for CSFRD must be provided by subclasses.
    """

    def CSFRD(self, z):
        '''The cosmic star formation rate density at a given redshift.
        
        Follows the Madau & Dickinson (2014) cosmic star formation rate density formula.
        
        Parameters
        ----------
        z : float or np.array
            Cosmological redshift.
        
        Returns
        -------
        float or array
            The cosmic star formation rate density at the given redshift.
        '''
        p = self.CSFRD_params
        return p["a"] * (1.0 + z) ** p["b"] / (1.0 + ((1.0 + z) / p["c"]) ** p["d"])

    def mean_metallicity(self, z):
        '''The mean metallicity at a given redshift
        
        Follows Madau & Fragos (2017) mean metallicity evolution
        
        Parameters
        ----------
        z : float or np.array
            Cosmological redshift.
        
        Returns
        -------
        float or array
            The mean metallicity at the given redshift.
        '''
        return 10 ** (0.153 - 0.074 * z ** 1.34) * Zsun

    def fSFR(self, z, metallicity_bins):
        '''Fraction of the SFR at a given redshift z in a given metallicity bin as described in Bavera et al. (2020).
        
        Parameters
        ----------
        z : np.array
            Cosmological redshift.
        metallicity_bins : array
            Metallicity bins edges in absolute metallicity.
        
        Returns
        -------
        array
            Fraction of the SFR in the given metallicity bin at the given redshift.
        '''
        sigma = self.std_log_metallicity_dist()
        # Compute mu; if z is an array, mu will be an array.
        mu = np.log10(self.mean_metallicity(z)) - sigma ** 2 * np.log(10) / 2.0
        # Ensure mu is an array for consistency
        mu_array = np.atleast_1d(mu)
        
        
        # Use the first value of mu for normalisation
        norm = stats.norm.cdf(np.log10(self.Z_max), mu_array, sigma)
        
        fSFR = np.empty((len(mu_array), len(metallicity_bins) - 1))
        
        fSFR[:, :] = np.array(
            [
                (
                    stats.norm.cdf(np.log10(metallicity_bins[1:]), m, sigma)
                    - stats.norm.cdf(np.log10(metallicity_bins[:-1]), m, sigma)
                )
                for m in mu_array
            ]
        )  / norm[:, np.newaxis]
        
        if not self.select_one_met:
            fSFR[:, 0] = stats.norm.cdf(np.log10(metallicity_bins[1]), mu_array, sigma)/norm
            fSFR[:, -1] = 1 - (stats.norm.cdf(np.log10(metallicity_bins[-2]), mu_array, sigma)/norm)
            
        return fSFR

class MadauDickinson14(MadauBase):
    '''Madau & Dickinson (2014) star formation history model using the 
    mean metallicity evolution of Madau & Fragos (2017).
    
    Madau & Dickinson (2014), ARA&A, 52, 415
    https://ui.adsabs.harvard.edu/abs/2014ARA%26A..52..415M/abstract
    '''

    def __init__(self, MODEL):
        super().__init__(MODEL)
        # Parameters for Madau+Dickinson14 CSFRD
        self.CSFRD_params = {
            "a": 0.015,
            "b": 2.7,
            "c": 2.9,
            "d": 5.6,
        }

class MadauFragos17(MadauBase):
    '''The Madau & Fragos (2017) star formation history model with the 
    mean metallicity evolution of Madau & Fragos (2017).
    
    Madau & Fragos (2017), ApJ, 840, 39
    http://adsabs.harvard.edu/abs/2017ApJ...840...39M
    '''

    def __init__(self, MODEL):
        super().__init__(MODEL)
        # Parameters for Madau+Fragos17 CSFRD
        self.CSFRD_params = {
            "a": 0.01,
            "b": 2.6,
            "c": 3.2,
            "d": 6.2,
        }
  
class Neijssel19(MadauBase):
    '''The Neijssel et al. (2019) star formation history model, which fits 
    the Madau & Dickinson (2014) cosmic star formation rate density formula
    with the BBH merger rate and uses a truncated log-normal distribution for
    the mean metallicity distribution.    
    The mean metallicity evolution follows the Langer and Normal parameterisation
    also fitted to the BBH merger rate.
    
    Neijssel et al. (2019), MNRAS, 490, 3740
    http://adsabs.harvard.edu/abs/2019MNRAS.490.3740N
    '''
    
    
    
    def __init__(self, MODEL):
        super().__init__(MODEL)
        # Parameters for Neijssel+19 CSFRD
        self.CSFRD_params = {
            "a": 0.01,
            "b": 2.77,
            "c": 2.9,
            "d": 4.7,
        }
    
    # overwrite mean_metallicity method of MadauBase    
    def mean_metallicity(self, z):
        return 0.035 * 10 ** (-0.23 * z)
    
    # overwrite std_log_metallicity_dist method of MadauBase
    # TODO: rewrite such that sigma is just changed for the Neijssel+19 case
    def fSFR(self, z, metallicity_bins):
        # assume a truncated ln-normal distribution of metallicities
        sigma = self.std_log_metallicity_dist()
        mu = np.log(self.mean_metallicity(z)) - sigma**2 / 2.0
        # renormalisation constant
        norm = stats.norm.cdf(np.log(self.Z_max), mu[0], sigma)
        fSFR = np.empty((len(z), len(metallicity_bins) - 1))
        fSFR[:, :] = np.array(
            [
                (
                    stats.norm.cdf(np.log(metallicity_bins[1:]), m, sigma) / norm
                    - stats.norm.cdf(np.log(metallicity_bins[:-1]), m, sigma) / norm
                )
                for m in mu
            ]
        )
        if not self.select_one_met:
            fSFR[:, 0] = stats.norm.cdf(np.log(metallicity_bins[1]), mu, sigma) / norm
            fSFR[:,-1] = norm - stats.norm.cdf(np.log(metallicity_bins[-2]), mu, sigma)/norm
        return fSFR
    
class IllustrisTNG(SFHBase):
    '''The IllustrisTNG star formation history model.
    
    Uses the TNG100-1 model from the IllustrisTNG simulation.
    
    https://www.tng-project.org/
    '''
    
    def __init__(self, MODEL):
        super().__init__(MODEL)
        # load the TNG data
        illustris_data = self._get_illustrisTNG_data()
        self.SFR = illustris_data["SFR"]
        self.redshifts = illustris_data["redshifts"]
        self.Z = illustris_data["mets"]
        self.M = illustris_data["M"]  # Msun
        
    def _get_illustrisTNG_data(self, verbose=False):
        """Load IllustrisTNG SFR dataset."""
        if verbose:
            print("Loading IllustrisTNG data...")
        return np.load(os.path.join(PATH_TO_POSYDON_DATA, "SFR/IllustrisTNG.npz"))
    
    def CSFRD(self, z):
        SFR_interp = interp1d(self.redshifts, self.SFR)
        return SFR_interp(z)
        
    def mean_metallicity(self, z):
        out = np.zeros_like(self.redshifts)
        for i in range(len(out)):
            if np.sum(self.M[i, :]) == 0:
                out[i] = 0
            else:
                out[i] = np.average(self.Z, weights=self.M[i, :])
        Z_interp = interp1d(self.redshifts, out)
        return Z_interp(z)
                
    def fSFR(self, z, metallicity_bins):
        # only use data within the metallicity bounds (no lower bound)
        Z_max_mask = self.Z <= self.Z_max
        redshift_indices = np.array([np.where(self.redshifts <= i)[0][0] for i in z])
        Z_dist = self.M[:, Z_max_mask][redshift_indices]
        fSFR = np.zeros((len(z), len(metallicity_bins) - 1))
        
        for i in range(len(z)):
            if Z_dist[i].sum() == 0.0:
                continue
            else:
                # Add a final point to the CDF and metallicities to ensure normalisation to 1
                Z_dist_cdf = np.cumsum(Z_dist[i]) / Z_dist[i].sum()
                Z_dist_cdf = np.append(Z_dist_cdf, 1)
                Z_x_values = np.append(np.log10(self.Z[Z_max_mask]), 0)
                Z_dist_cdf_interp = interp1d(Z_x_values, Z_dist_cdf)

                fSFR[i, :] = (Z_dist_cdf_interp(np.log10(metallicity_bins[1:])) -
                              Z_dist_cdf_interp(np.log10(metallicity_bins[:-1])))

                if not self.select_one_met:
                    if len(metallicity_bins) == 2:
                        fSFR[i, 0] = 1
                    else:
                        fSFR[i, 0] = Z_dist_cdf_interp(np.log10(metallicity_bins[1]))
                        fSFR[i, -1] = 1 - Z_dist_cdf_interp(np.log10(metallicity_bins[-2]))
                        
        return fSFR
    
class Chruslinska21(SFHBase):
    '''The Chruślińska+21 star formation history model.
    
    Chruślińska et al. (2021), MNRAS, 508, 4994
    https://ui.adsabs.harvard.edu/abs/2021MNRAS.508.4994C/abstract
    
    Data source: 
    https://ftp.science.ru.nl/astro/mchruslinska/Chruslinska_et_al_2021/
    
    
    '''
    def __init__(self, MODEL):
        '''Initialise the Chruslinska+21 model
        
        Parameters
        ----------
        MODEL : dict
            Model parameters. Chruslinska+21 requires the following parameters:
            - sub_model : str
                The sub-model to use. This is the name of the file containing the data.
            - Z_solar_scaling : str
                The scaling of the solar metallicity. Options are:
                - Asplund09
                - AndersGrevesse89
                - GrevesseSauval98
                - Villante14
        '''
        if "sub_model" not in MODEL:
            raise ValueError("Sub-model not given!")
        if 'Z_solar_scaling' not in MODEL:
            raise ValueError("Z_solar_scaling not given!")
        
        super().__init__(MODEL)
        self._load_chruslinska_data()
        
    def _load_chruslinska_data(self, verbose=False):
        '''load the data from the Chruslinska+21 models
        Transforms the data to the format used in the classes.
        
        Parameters
        ----------
        verbose : bool, optional
            Print information about the data loading.
        
        '''  
        # oxygen to hydrogen abundance ratio ( FOH == 12 + log(O/H) )
        # as used in the calculations - do not change
        # This is the metallicity bin edges used in the Chruslinska+21 calculations
        FOH_min, FOH_max = 5.3, 9.7
        self.FOH_bins = np.linspace(FOH_min,FOH_max, 200)
        self.dFOH=self.FOH_bins[1]-self.FOH_bins[0]
        # I need to use the Z_solar_scaling parameter to convert the FOH bins to absolute metallicity
        # I will use the solar metallicity as the reference point
        self.Z = self._FOH_to_Z(self.FOH_bins)
        
        self._data_folder = os.path.join(PATH_TO_POSYDON_DATA, "SFR/Chruslinska+21")
        _, self.redshifts, delta_T = self._load_redshift_data(verbose)
        M = self._load_raw_data()
        self.SFR = np.array( [M[ii]/(1e6*delta_T[ii]) for ii in range(len(delta_T))])/self.dFOH

    def _FOH_to_Z(self, FOH):
        # scalings from Chruslinksa+21
        if self.Z_solar_scaling == 'Asplund09':
            Zsun, FOHsun = [0.0134, 8.69]
        elif self.Z_solar_scaling == 'AndersGrevesse89':
            Zsun,FOHsun = [0.017, 8.83]
        elif self.Z_solar_scaling == 'GrevesseSauval98':
            Zsun,FOHsun = [0.0201, 8.93]
        elif self.Z_solar_scaling == 'Villante14':
            Zsun,FOHsun = [0.019, 8.85]
        else:
            raise ValueError("Invalid Z_solar_scaling!")
        logZ = np.log10(Zsun) + FOH - FOHsun
        ZZ=10**logZ
        return ZZ
        
    def mean_metallicity(self, z):
        '''Calculate the mean metallicity at a given redshift
        
        Parameters
        ----------
        z : float or array-like
            Cosmological redshift.
            
        Returns
        -------
        float or array-like
            The mean metallicity at the given redshift(s).
        ''' 
        mean_over_redshift = np.zeros_like(self.redshifts)
        for i in range(len(mean_over_redshift)):
            if np.sum(self.SFR[i]) == 0:
                mean_over_redshift[i] = 0
            else:
                mean_over_redshift[i] = np.average(self.Z, weights=self.SFR[i,:]*self.dFOH)
        
        Z_interp = interp1d(self.redshifts, mean_over_redshift)
        return Z_interp(z)
    
    def fSFR(self, z, metallicity_bins):
        '''Calculate the fractional SFR as a function of redshift and metallicity bins
        
        Parameters
        ----------
        z : float or array-like
            Cosmological redshift.
        metallicity_bins : array
            Metallicity bins edges in absolute metallicity.
        
        Returns
        -------
        array
            Fraction of the SFR in the given metallicity bin at the given redshift.
        '''
        # only use data within the metallicity bounds (no lower bound)
        Z_max_mask = self.Z <= self.Z_max
        redshift_indices = np.array([np.where(self.redshifts <= i)[0][0] for i in z])
        Z_dist = self.SFR[:, Z_max_mask][redshift_indices]
        fSFR = np.zeros((len(z), len(metallicity_bins) - 1))
        
        for i in range(len(z)):
            if Z_dist[i].sum() == 0.0:
                continue
            else:
                # Add a final point to the CDF and metallicities to ensure normalisation to 1
                Z_dist_cdf = np.cumsum(Z_dist[i]) / Z_dist[i].sum()
                Z_dist_cdf = np.append(Z_dist_cdf, 1)
                Z_x_values = np.append(np.log10(self.Z[Z_max_mask]), 0)
                Z_dist_cdf_interp = interp1d(Z_x_values, Z_dist_cdf)

                fSFR[i, :] = (Z_dist_cdf_interp(np.log10(metallicity_bins[1:])) -
                              Z_dist_cdf_interp(np.log10(metallicity_bins[:-1])))

                if not self.select_one_met:
                    if len(metallicity_bins) == 2:
                        fSFR[i, 0] = 1
                    else:
                        fSFR[i, 0] = Z_dist_cdf_interp(np.log10(metallicity_bins[1]))
                        fSFR[i, -1] = 1 - Z_dist_cdf_interp(np.log10(metallicity_bins[-2]))
                        
        return fSFR
        
    def _load_redshift_data(self, verbose=False):
        '''Load the redshift data from a Chruslinsk+21 model file.
        
        Returns
        -------
        time : array
            the center of the time bins 
        redshift : array
            the redshifts corresponding to the time bins
        delt : array
            the width of the time bins
        '''
        if verbose:
            print("Loading redshift data...")
        
        time, redshift, delt = np.loadtxt(
            os.path.join(self._data_folder, 'Time_redshift_deltaT.dat'), unpack=True)
        return time, redshift, delt
             
    def _load_raw_data(self):
        '''Read the sub-model data from the file
        
        The data structure is as follows:
        - mass per unit (comoving) volume formed in each z (row) - FOH (column) bin
        
        Returns
        -------
        array
            Mass formed per unit volume in each redshift and FOH bin
        '''
        input_file = os.path.join(self._data_folder, f'{self.sub_model}.dat')
        data = np.loadtxt(input_file)
        return data

    def CSFRD(self, z):
        '''Interpolate the cosmic star formation rate density at the given redshift(s)
        
        Parameters
        ----------
        z : float or array-like
            Cosmological redshift.
        
        Returns
        -------
        float or array-like
            The cosmic star formation rate density at the given redshift(s).
        '''
        SFR_interp = interp1d(self.redshifts, np.sum(self.SFR*self.dFOH, axis=1))
        return SFR_interp(z)
    
        
def get_SFH_model(MODEL):
    '''Return the appropriate SFH model based on the given parameters
    
    Parameters
    ----------
    MODEL : dict
        Model parameters.
    
    Returns
    -------
    SFHBase
        The SFH model instance.
    '''
    if MODEL["SFR"] == "Madau+Fragos17":
        return MadauFragos17(MODEL)
    elif MODEL['SFR'] == "Madau+Dickinson14":
        return MadauDickinson14(MODEL)
    elif MODEL['SFR'] == "Neijssel+19":
        return Neijssel19(MODEL)
    elif MODEL['SFR'] == "IllustrisTNG":
        return IllustrisTNG(MODEL)
    elif MODEL['SFR'] == "Chruslinska+21":
        return Chruslinska21(MODEL)
    else:
        raise ValueError("Invalid SFR!")

def SFR_per_Z_at_z(z, met_bins, MODEL):
    """Calculate the SFR per metallicity bin at a given redshift(s)
    
    Parameters
    ----------
    z : float or array-like
        Cosmological redshift.
    met_bins : array
        Metallicity bins edges in absolute metallicity.
    MODEL : dict
        Model parameters.
    
    Returns
    -------
    SFH : 2D array
        Star formation history per metallicity bin at the given redshift(s).
    
    """
    SFH = get_SFH_model(MODEL)
    return SFH(z, met_bins)


def get_formation_times(N_binaries, star_formation="constant", **kwargs):
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
    RNG = kwargs.get("RNG", np.random.default_rng())

    scenario = star_formation.lower()

    if scenario == "burst":
        burst_time = kwargs.get("burst_time", 0.0)
        return np.ones(N_binaries) * burst_time

    max_time_default = kwargs.get("max_simulation_time", age_of_universe)
    max_time = kwargs.get("max_time", max_time_default)

    if scenario == "constant":
        min_time = kwargs.get("min_time", 0.0)
        return RNG.uniform(size=N_binaries, low=min_time, high=max_time)

    if scenario in ["custom_linear", "custom_log10"]:
        custom_ages_file = kwargs.get("custom_ages_file")
        x, y = np.loadtxt(custom_ages_file, unpack=True)
        current_binary_ages = rejection_sampler(x, y, N_binaries)
        if "log10" in scenario:
            current_binary_ages = 10.0**current_binary_ages
        return max_time - current_binary_ages

    if scenario in ["custom_linear_histogram", "custom_log10_histogram"]:
        custom_ages_file = kwargs.get("custom_ages_file")
        x, y = read_histogram_from_file(custom_ages_file)
        current_binary_ages = histogram_sampler(x, y)
        if "log10" in scenario:
            current_binary_ages = 10.0**current_binary_ages
        return max_time - current_binary_ages

    raise ValueError(
        "Unknown star formation scenario '{}' given. Valid options: {}".format(
            star_formation, ",".join(SFH_SCENARIOS)
        )
    )

#### OLD CODE BELOW ####



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
        return (
            0.01 * (1.0 + z) ** 2.6 / (1.0 + ((1.0 + z) / 3.2) ** 6.2)
        )  # M_sun yr^-1 Mpc^-3
    elif SFR == "Madau+Dickinson14":
        return (
            0.015 * (1.0 + z) ** 2.7 / (1.0 + ((1.0 + z) / 2.9) ** 5.6)
        )  # M_sun yr^-1 Mpc^-3
    elif SFR == "Neijssel+19":
        return (
            0.01 * (1.0 + z) ** 2.77 / (1.0 + ((1.0 + z) / 2.9) ** 4.7)
        )  # M_sun yr^-1 Mpc^-3
    elif SFR == "IllustrisTNG":
        illustris_data = get_illustrisTNG_data()
        SFR = illustris_data["SFR"]  # M_sun yr^-1 Mpc^-3
        redshifts = illustris_data["redshifts"]
        SFR_interp = interp1d(redshifts, SFR)
        return SFR_interp(z)
    else:
        raise ValueError("Invalid SFR!")

def get_illustrisTNG_data(verbose=False):
    """Load IllustrisTNG SFR dataset."""
    if verbose:
        print("Loading IllustrisTNG data...")
    return np.load(os.path.join(PATH_TO_POSYDON_DATA, "SFR/IllustrisTNG.npz"))

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
        return 10 ** (0.153 - 0.074 * z**1.34) * Zsun
    elif SFR == "Neijssel+19":
        return 0.035 * 10 ** (-0.23 * z)
    else:
        raise ValueError("Invalid SFR!")

def std_log_metallicity_dist(sigma):
    """Standard deviation of the log-metallicity distribution.

    Returns
    -------
    double
        Standard deviation of the adopted distribution.

    """
    if isinstance(sigma, str):
        if sigma == "Bavera+20":
            return 0.5
        elif sigma == "Neijssel+19":
            return 0.39
        else:
            raise ValueError("Uknown sigma choice!")
    elif isinstance(sigma, float):
        return sigma
    else:
        raise ValueError(f"Invalid sigma value {sigma}!")


def SFR_Z_fraction_at_given_redshift(
    z, SFR, sigma, metallicity_bins, Z_max, select_one_met
):
    """'Fraction of the SFR at a given redshift z in a given metallicity bin as in Eq. (B.8) of Bavera et al. (2020).

    Parameters
    ----------
    z : np.array
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

    Returns
    -------
    array
        Fraction of the SFR in the given metallicity bin at the given redshift.
        In absolute metallicity.
    """

    if SFR == "Madau+Fragos17" or SFR == "Madau+Dickinson14":
        sigma = std_log_metallicity_dist(sigma)
        mu = np.log10(mean_metallicity(SFR, z)) - sigma**2 * np.log(10) / 2.0
        # renormalisation constant. We can use mu[0], since we integrate over the whole metallicity range
        norm = stats.norm.cdf(np.log10(Z_max), mu[0], sigma)
        fSFR = np.empty((len(z), len(metallicity_bins) - 1))
        fSFR[:, :] = np.array(
            [(stats.norm.cdf(np.log10(metallicity_bins[1:]), m, sigma) / norm
                    - stats.norm.cdf(np.log10(metallicity_bins[:-1]), m, sigma) / norm
             ) for m in mu ]
        )
        if not select_one_met:
            fSFR[:, 0] = stats.norm.cdf(np.log10(metallicity_bins[1]), mu, sigma) / norm
            fSFR[:,-1] = norm - stats.norm.cdf(np.log10(metallicity_bins[-2]), mu, sigma)/norm

    elif SFR == "Neijssel+19":
        # assume a truncated ln-normal distribution of metallicities
        sigma = std_log_metallicity_dist(sigma)
        mu = np.log(mean_metallicity(SFR, z)) - sigma**2 / 2.0
        # renormalisation constant
        norm = stats.norm.cdf(np.log(Z_max), mu[0], sigma)
        fSFR = np.empty((len(z), len(metallicity_bins) - 1))
        fSFR[:, :] = np.array(
            [
                (
                    stats.norm.cdf(np.log(metallicity_bins[1:]), m, sigma) / norm
                    - stats.norm.cdf(np.log(metallicity_bins[:-1]), m, sigma) / norm
                )
                for m in mu
            ]
        )
        if not select_one_met:
            fSFR[:, 0] = stats.norm.cdf(np.log(metallicity_bins[1]), mu, sigma) / norm
            fSFR[:,-1] = norm - stats.norm.cdf(np.log(metallicity_bins[-2]), mu, sigma)/norm

    elif SFR == "IllustrisTNG":
        # numerically itegrate the IlluystrisTNG SFR(z,Z)
        illustris_data = get_illustrisTNG_data()
        redshifts = illustris_data["redshifts"]
        Z = illustris_data["mets"]
        M = illustris_data["M"]  # Msun
        # only use data within the metallicity bounds (no lower bound)
        Z_max_mask = Z <= Z_max
        redshift_indices = np.array([np.where(redshifts <= i)[0][0] for i in z])
        Z_dist = M[:, Z_max_mask][redshift_indices]
        fSFR = np.zeros((len(z), len(metallicity_bins) - 1))
        
        for i in range(len(z)):
            if Z_dist[i].sum() == 0.0:
                continue
            else:
                # We add a final point to the CDF and metallicities to ensure normalisation to 1
                Z_dist_cdf = np.cumsum(Z_dist[i]) / Z_dist[i].sum()
                Z_dist_cdf = np.append(Z_dist_cdf, 1)
                Z_x_values = np.append(np.log10(Z[Z_max_mask]), 0)
                Z_dist_cdf_interp = interp1d(Z_x_values, Z_dist_cdf)

                fSFR[i, :] = (Z_dist_cdf_interp(np.log10(metallicity_bins[1:]))
                              - Z_dist_cdf_interp(np.log10(metallicity_bins[:-1])))

                if not select_one_met:
                    # add the fraction of the SFR in the first and last bin
                    # or the only bin without selecting one metallicity
                    if len(metallicity_bins) == 2:
                        fSFR[i, 0] = 1
                    else:
                        fSFR[i, 0] = Z_dist_cdf_interp(np.log10(metallicity_bins[1]))
                        fSFR[i, -1] = 1 - Z_dist_cdf_interp(np.log10(metallicity_bins[-2]))
    else:
        raise ValueError("Invalid SFR!")

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
        Omega_L = 1.0 - Omega_m
        return (Omega_m * (1.0 + z) ** 3 + Omega_L) ** (1.0 / 2.0)

    def f(z, Z):
        if SFR == "Madau+Fragos17" or SFR == "Madau+Dickinson14":
            sigma = std_log_metallicity_dist(sigma)
            mu = np.log10(mean_metallicity(SFR, z)) - sigma**2 * np.log(10) / 2.0
            H_0 = cosmology.H0.to("1/yr").value  # yr
            # put a cutoff on metallicity at Z_max
            norm = stats.norm.cdf(np.log10(Z_max), mu, sigma)
            return (
                star_formation_rate(SFR, z)
                * stats.norm.pdf(np.log10(Z), mu, sigma)
                / norm
                * (H_0 * (1.0 + z) * E(z)) ** (-1)
            )
        elif SFR == "Neijssel+19":
            sigma = std_log_metallicity_dist(sigma)
            mu = np.log10(mean_metallicity(SFR, z)) - sigma**2 / 2.0
            H_0 = cosmology.H0.to("1/yr").value  # yr
            return (
                star_formation_rate(SFR, z)
                * stats.norm.pdf(np.log(Z), mu, sigma)
                * (H_0 * (1.0 + z) * E(z)) ** (-1)
            )
        else:
            raise ValueError("Invalid SFR!")

    return sp.integrate.quad(f, 1e-10, np.inf, args=(Z,))[0]  # M_sun yr^-1 Mpc^-3
