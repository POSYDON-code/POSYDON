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
import pandas as pd
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
        '''Initialise the SFH model
        
        Adds the model parameters as attributes.
        
        Parameters
        ----------
        MODEL : dict
            Model parameters.
        '''
        self.MODEL = MODEL
        # Automatically attach all model parameters as attributes
        for key, value in MODEL.items():
            setattr(self, key, value)

    @abstractmethod
    def CSFRD(self, z):
        """Compute the cosmic star formation rate density.
        
        This is an abstract method that must be implemented by subclasses.
        The implementation should calculate and return the cosmic star formation
        rate density at the given redshift(s).
        
        Parameters
        ----------
        z : float or array-like
            Cosmological redshift.
        
        Returns
        -------
        float or array-like
            The cosmic star formation rate density at the given redshift(s).
        """
        pass

    @abstractmethod
    def mean_metallicity(self, z):
        """Return the mean metallicity at redshift z.

        This is an abstract method that must be implemented by subclasses.
        The implementation should calculate and return the mean metallicity
        at the given redshift(s).
        
        Parameters
        ----------
        z : float or array-like
            Cosmological redshift.
        
        Returns
        -------
        float or array-like
            The mean metallicity at the given redshift(s).        
        """
        pass
        
    @abstractmethod
    def fSFR(self, z, metallicity_bins):
        """Compute the star formation rate fraction (fSFR) at a given redshift 
        using the specified metallicity bins.

        This is an abstract method that must be implemented by subclasses.
        The implementation should calculate and return the fractional SFR per
        metallicity bins at the provided redshift (z).

        Parameters
        ---------
        z : float or array-like 
            The redshift(s) at which to compute the star formation rate.
        metallicity_bins : list or array-like
            The metallicity bin boundaries or labels used in the computation to 
            account for different metallicity contributions.

        Returns
        -------
        float or array-like
            The calculated star formation rate at the given redshift(s) and
            metallicity bins in Msun/yr.

        Raises:
            NotImplementedError: If the subclass does not implement this method.
        """
        pass


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
    def __init__(self, MODEL):
        '''Initialise the MadauBase class
        
        Parameters
        ----------
        MODEL : dict
            Model parameters. MadauBase requires the following parameters:
            - sigma : float or str
                The standard deviation of the log-normal metallicity distribution.
                Options are:
                - Bavera+20
                - Neijssel+19
                - float
            - Z_max : float
                The maximum metallicity in absolute units.
            - select_one_met : bool
                If True, the SFR is calculated for a single metallicity bin.
            
        '''
        if "sigma" not in MODEL:
            raise ValueError("sigma not given!")
        if "Z_max" not in MODEL:
            raise ValueError("Z_max not given!")
        if "select_one_met" not in MODEL:
            raise ValueError("select_one_met not given!")
        super().__init__(MODEL)
        self.CSFRD_params = None
    
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
    
    def std_log_metallicity_dist(self):
        '''Return the standard deviation of the log-normal metallicity distribution
        
        Either recognised the strings "Bavera+20" (sigma=0.5) 
        or "Neijssel+19" (sigma=0.39) or a float value.
        
        Returns
        -------
        float
            The standard deviation of the log-normal metallicity distribution.
        '''
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
        
        # Use mu for normalisation
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
        '''Initialise the Madau & Dickinson (2014) SFH model with the
        metallicity evolution of Madau & Fragos (2017).
        
        Parameters
        ----------
        MODEL : dict
            Model parameters. Madau+14 requires the following parameters:
            - sigma : float or str
                The standard deviation of the log-normal metallicity distribution.
                Options are:
                - Bavera+20
                - Neijssel+19
                - float
            - Z_max : float
                The maximum metallicity in absolute units.
            - select_one_met : bool
                If True, the SFR is calculated for a single metallicity bin.
        '''
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
    metallicity evolution of Madau & Fragos (2017).
    
    Madau & Fragos (2017), ApJ, 840, 39
    http://adsabs.harvard.edu/abs/2017ApJ...840...39M
    '''

    def __init__(self, MODEL):
        '''Initialise the Madau+17 model
        
        Parameters
        ----------
        MODEL : dict
            Model parameters. Madau+17 requires the following parameters:
            - sigma : float or str
                The standard deviation of the log-normal metallicity distribution.
                Options are:
                - Bavera+20
                - Neijssel+19
                - float
            - Z_max : float
                The maximum metallicity in absolute units.
            - select_one_met : bool
                If True, the SFR is calculated for a single metallicity bin.    
        '''
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
        '''Initialise the Neijssel+19 model
        
        Parameters
        ----------
        MODEL : dict
            Model parameters. Neijssel+19 requires the following parameters:
            - sigma : float or str
                The standard deviation of the log-normal metallicity distribution.
                Options are:
                - Bavera+20
                - Neijssel+19
                - float
            - Z_max : float
                The maximum metallicity in absolute units.
            - select_one_met : bool
                If True, the SFR is calculated for a single metallicity bin.    
        '''
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
        '''Calculate the mean metallicity at a given redshift
        
        Overwrites the mean_metallicity method of MadauBase class.
        
        Parameters
        ----------
        z : float or array-like
            Cosmological redshift.
            
        Returns
        -------
        float or array-like
            The mean metallicity at the given redshift(s).    
        '''
        return 0.035 * 10 ** (-0.23 * z)
    
    # TODO: rewrite such that sigma is just changed for the Neijssel+19 case
    def fSFR(self, z, metallicity_bins):
        '''Fraction of the SFR at a given redshift z in a given metallicity bin
        as described in Neijssel et al. (2019).
        
        Overwrites the fSFR method of MadauBase class.
        
        Parameters
        ----------
        z : np.array
            Cosmological redshift.
        metallicity_bins : array
            Metallicity bins edges in absolute metallicity.
        
        Returns
        -------
        array
            Fraction of the SFR in the given metallicity bins at the 
            given redshift.
        '''
        # assume a truncated ln-normal distribution of metallicities
        sigma = self.std_log_metallicity_dist()
        mu = np.log(self.mean_metallicity(z)) - sigma**2 / 2.0
        # renormalisation constant
        norm = stats.norm.cdf(np.log(self.Z_max), mu, sigma)
        fSFR = np.empty((len(z), len(metallicity_bins) - 1))
        fSFR[:, :] = np.array(
            [
                (
                    stats.norm.cdf(np.log(metallicity_bins[1:]), m, sigma) 
                    - stats.norm.cdf(np.log(metallicity_bins[:-1]), m, sigma)
                )
                for m in mu
            ]
        ) / norm[:, np.newaxis]
        if not self.select_one_met:
            fSFR[:, 0] = stats.norm.cdf(np.log(metallicity_bins[1]), mu, sigma) / norm
            fSFR[:,-1] = 1 - stats.norm.cdf(np.log(metallicity_bins[-2]), mu, sigma)/norm
        return fSFR
    
class IllustrisTNG(SFHBase):
    '''The IllustrisTNG star formation history model.
    
    Uses the TNG100-1 model from the IllustrisTNG simulation.
    
    https://www.tng-project.org/
    '''
    
    def __init__(self, MODEL):
        '''Initialise the IllustrisTNG model
        
        Parameters
        ----------
        MODEL : dict
            Model parameters. IllustrisTNG requires the following parameters:
            - Z_max : float
                The maximum metallicity in absolute units.
            - select_one_met : bool
                If True, the SFR is calculated for a single metallicity bin.
        '''
        super().__init__(MODEL)
        # load the TNG data
        illustris_data = self._get_illustrisTNG_data()
        self.CSFRD_data = illustris_data["SFR"]
        self.redshifts = illustris_data["redshifts"]
        self.Z = illustris_data["mets"]
        self.M = illustris_data["M"]  # Msun
        
    def _get_illustrisTNG_data(self, verbose=False):
        '''Load IllustrisTNG SFR dataset into the class.
        
        Parameters
        ----------
        verbose : bool, optional
            Print information about the data loading.
        '''
        if verbose:
            print("Loading IllustrisTNG data...")
        return np.load(os.path.join(PATH_TO_POSYDON_DATA, "SFR/IllustrisTNG.npz"))
    
    def CSFRD(self, z):
        '''The cosmic star formation rate density at a given redshift.
        
        Parameters
        ----------
        z : float or np.array
            Cosmological redshift.
        
        Returns
        -------
        float or array
            The cosmic star formation rate density at the given redshift(s).
        '''
        SFR_interp = interp1d(self.redshifts, self.CSFRD_data)
        return SFR_interp(z)
        
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
        out = np.zeros_like(self.redshifts)
        for i in range(len(out)):
            if np.sum(self.M[i, :]) == 0:
                out[i] = np.nan
            else:
                out[i] = np.average(self.Z, weights=self.M[i, :])
        Z_interp = interp1d(self.redshifts, out)
        return Z_interp(z)
                
    def fSFR(self, z, metallicity_bins):
        '''Calculate the fractional SFR as a function of redshift and 
        metallicity bins.
        
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
            - Z_max : float
        '''
        if "sub_model" not in MODEL:
            raise ValueError("Sub-model not given!")
        if 'Z_solar_scaling' not in MODEL:
            raise ValueError("Z_solar_scaling not given!")
        
        super().__init__(MODEL)
        self._load_chruslinska_data()
        
    def _load_chruslinska_data(self, verbose=False):
        '''Load the data from the Chruslinska+21 models.
        Transforms the data to the format used in the classes.
        
        Parameters
        ----------
        verbose : bool, optional
            Print information about the data loading.
        '''  
        # oxygen to hydrogen abundance ratio ( FOH == 12 + log(O/H) )
        # as used in the calculations - do not change
        # This is the metallicity bin edges used in the Chruslinska+21 calculations
        FOH_min = 5.3
        FOH_max = 9.7
        self.FOH_bins = np.linspace(FOH_min, FOH_max, 200)
        self.dFOH = self.FOH_bins[1] - self.FOH_bins[0]
        # I need to use the Z_solar_scaling parameter to convert the FOH bins to absolute metallicity
        # I will use the solar metallicity as the reference point
        self.Z = self._FOH_to_Z(self.FOH_bins)
        
        self._data_folder = os.path.join(PATH_TO_POSYDON_DATA, "SFR/Chruslinska+21")
        _, self.redshifts, delta_T = self._load_redshift_data(verbose)
        M = self._load_raw_data()
        self.SFR_data = np.array( [M[ii]/(1e6*delta_T[ii]) for ii in range(len(delta_T))])/self.dFOH

    def _FOH_to_Z(self, FOH):
        '''Convert the oxygen to hydrogen abundance ratio to absolute metallicity
        
        Parameters
        ----------
        FOH : float or array-like
            The oxygen to hydrogen abundance ratio.
        
        Returns
        -------
        float or array-like
            The absolute metallicity.    
        '''
        # scalings from Chruslinksa+21
        if self.Z_solar_scaling == 'Asplund09':
            Zsun = 0.0134
            FOHsun = 8.69
        elif self.Z_solar_scaling == 'AndersGrevesse89':
            Zsun = 0.017
            FOHsun = 8.83
        elif self.Z_solar_scaling == 'GrevesseSauval98':
            Zsun = 0.0201
            FOHsun = 8.93
        elif self.Z_solar_scaling == 'Villante14':
            Zsun = 0.019
            FOHsun = 8.85
        else:
            raise ValueError("Invalid Z_solar_scaling!"+
                             "Options are: Asplund09, AndersGrevesse89,"+
                             "GrevesseSauval98, Villante14")
        logZ = np.log10(Zsun) + FOH - FOHsun
        return 10**logZ
        
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
            if np.sum(self.SFR_data[i]) == 0:
                mean_over_redshift[i] = 0
            else:
                mean_over_redshift[i] = np.average(self.Z, weights=self.SFR_data[i,:]*self.dFOH)
        
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
        Z_dist = self.SFR_data[:, Z_max_mask][redshift_indices]
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
        SFR_interp = interp1d(self.redshifts, np.sum(self.SFR_data*self.dFOH, axis=1))
        return SFR_interp(z)
    
class Fujimoto24(MadauBase):
    '''The Fujimoto et al. (2024) star formation history model.
    mean metallicity evolution of Madau & Fragos (2017).
    
    Fujimoto et al. (2024), ApJ SS, 275, 2, 36, 59
    https://ui.adsabs.harvard.edu/abs/2024ApJS..275...36F/abstract
    '''
    def __init__(self, MODEL):
        '''Initialise the Fujimoto+24 model
        
        Parameters
        ----------
        MODEL : dict
            Model parameters. Fujimoto+24 requires the following parameters:
            - sigma : float or str
                The standard deviation of the log-normal metallicity distribution.
                Options are:
                - Bavera+20
                - Neijssel+19
                - float
            - Z_max : float
                The maximum metallicity in absolute units.
            - select_one_met : bool
                If True, the SFR is calculated for a single metallicity bin.
        '''
        super().__init__(MODEL)
        # Parameters for Fujimoto+24 CSFRD
        self.CSFRD_params = {
            "a": 0.010,
            "b": 2.8,
            "c": 3.3,
            "d": 6.6,
        }   
    
class Zalava21(MadauBase):
    
    def __init__(self, MODEL):
        '''Initialise the Zalava+21 model
        
        Requires the following parameters:
        - sub_model : str
            Either min or max
        '''
        if 'sub_model' not in MODEL:
            raise ValueError("Sub-model not given!")
        
        super().__init__(MODEL)
        self._load_zalava_data()
        
    def _load_zalava_data(self):
        '''Load the data from the Zalava+21 models
        Transforms the data to the format used in the classes.
        
        '''  
        data_file = os.path.join(PATH_TO_POSYDON_DATA, "SFR/Zalava+21.txt")
        tmp_data = pd.read_csv(data_file, names=['redshift', 'SFRD_min', 'SFRD_max'], skiprows=1, sep='\s+')
        self.redshifts = tmp_data['redshift'].values
        if self.sub_model == 'min':
            self.SFR_data = tmp_data['SFRD_min'].values
        elif self.sub_model == 'max':
            self.SFR_data = tmp_data['SFRD_max'].values
        else:
            raise ValueError("Invalid sub-model!")

    # overwrite the CSFRD method of MadauBase        
    def CSFRD(self, z):
        SFR_interp = interp1d(self.redshifts, self.SFR_data)
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
    elif MODEL['SFR'] == 'Fujimoto+24':
        return Fujimoto24(MODEL)
    elif MODEL['SFR'] == "Neijssel+19":
        return Neijssel19(MODEL)
    elif MODEL['SFR'] == "IllustrisTNG":
        return IllustrisTNG(MODEL)
    elif MODEL['SFR'] == "Chruslinska+21":
        return Chruslinska21(MODEL)
    elif MODEL['SFR'] == "Zalava+21":
        return Zalava21(MODEL)
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