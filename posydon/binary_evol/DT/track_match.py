__authors__ = [
    "Devina Misra <devina.misra@unige.ch>",
    "Zepei Xing <Zepei.Xing@unige.ch>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Nam Tran <tranhn03@gmail.com>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
    "Camille Liotine <cliotine@u.northwestern.edu>",
    "Seth Gossage <seth.gossage@northwestern.edu>"
]

import os
import time
import numpy as np
import pandas as pd
from scipy.optimize import root
from scipy.optimize import minimize
from scipy.interpolate import PchipInterpolator

import posydon.utils.constants as const
from posydon.config import PATH_TO_POSYDON_DATA
from posydon.utils.posydonwarning import Pwarn
from posydon.interpolation.data_scaling import DataScaler
from posydon.interpolation.interpolation import GRIDInterpolator
from posydon.utils.interpolators import PchipInterpolator2
from posydon.utils.posydonerror import (NumericalError, MatchingError)
from posydon.utils.common_functions import (convert_metallicity_to_string,
                                            set_binary_to_failed)


LIST_ACCEPTABLE_STATES_FOR_HMS = ["H-rich_Core_H_burning",
                                  "accreted_He_Core_H_burning"]

LIST_ACCEPTABLE_STATES_FOR_postMS = [
    "H-rich_Shell_H_burning",
    "H-rich_Core_He_burning",
    "H-rich_Central_He_depleted",
    "H-rich_Core_C_burning",
    "H-rich_Central_C_depletion",
    "H-rich_non_burning",
    "accreted_He_non_burning"]

LIST_ACCEPTABLE_STATES_FOR_HeStar = [
    'accreted_He_Core_He_burning',
    'stripped_He_Core_He_burning',
    'stripped_He_Shell_He_burning',     # includes stars burning C in core
    'stripped_He_Central_He_depleted',  # includes stars burning C in core
    'stripped_He_Central_C_depletion',
    'stripped_He_non_burning'
    ]

MATCHING_WITH_RELATIVE_DIFFERENCE = ["center_he4"]

DEFAULT_TRANSLATION = {
    "time": "time",
    "orbital_period": "porb",
    "eccentricity": "ecc",
    "separation": "sep",
    "state": None,
    "event": None,
    "rl_relative_overflow_1": "rl_relative_overflow_1",
    "rl_relative_overflow_2": "rl_relative_overflow_2",
    "lg_mtransfer_rate": "lg_mtransfer_rate",
    "V_sys": None,
    "mass": "mass",
    "log_R": "log_R",
    "R": "R",
    "lg_mdot": "mdot",
    "log_L": "log_L",
    "lg_wind_mdot": "mdot",
    "lg_system_mdot": "lg_mdot",
    "he_core_mass": "he_core_mass",
    "he_core_radius": "he_core_radius",
    "c_core_mass": "c_core_mass",
    "c_core_radius": "c_core_radius",
    "o_core_mass": "o_core_mass",
    "o_core_radius": "o_core_radius",
    "center_h1": "center_h1",
    "center_he4": "center_he4",
    "center_c12": "center_c12",
    "center_o16": "center_o16",
    "center_n14": "center_n14",
    "surface_h1": "surface_h1",
    "surface_he4": "surface_he4",
    "surface_c12": "surface_c12",
    "surface_n14": "surface_n14",
    "surface_o16": "surface_o16",
    "center_gamma": "center_gamma",
    "log_LH": "log_LH",
    "log_LHe": "log_LHe",
    "log_LZ": "log_LZ",
    "log_Lnuc": "log_Lnuc",
    "c12_c12": "c12_c12",
    "avg_c_in_c_core": "avg_c_in_c_core",
    "surf_avg_omega_div_omega_crit": "surf_avg_omega_div_omega_crit",
    "surf_avg_omega": "omega",
    "total_moment_of_inertia": "inertia",
    "log_total_angular_momentum": "log_total_angular_momentum",
    "profile": None,
    "metallicity": None,
    "spin": "spin_parameter",
    "conv_env_top_mass": "conv_env_top_mass",
    "conv_env_bot_mass": "conv_env_bot_mass",
    "conv_env_top_radius": "conv_env_top_radius",
    "conv_env_bot_radius": "conv_env_bot_radius",
    "conv_env_turnover_time_g": "conv_env_turnover_time_g",
    "conv_env_turnover_time_l_b": "conv_env_turnover_time_l_b",
    "conv_env_turnover_time_l_t": "conv_env_turnover_time_l_t",
    "envelope_binding_energy": "envelope_binding_energy",
    "mass_conv_reg_fortides": "mass_conv_reg_fortides",
    "thickness_conv_reg_fortides": "thickness_conv_reg_fortides",
    "radius_conv_reg_fortides": "radius_conv_reg_fortides",
    "lambda_CE_1cent": "lambda_CE_1cent",
    "lambda_CE_10cent": "lambda_CE_10cent",
    "lambda_CE_30cent": "lambda_CE_30cent",
    "co_core_mass": "co_core_mass",
    "co_core_radius": "co_core_radius",
    "lambda_CE_pure_He_star_10cent": "lambda_CE_pure_He_star_10cent",
    "trap_radius": "trap_radius",
    "acc_radius": "acc_radius",
    "t_sync_rad_1": "t_sync_rad_1",
    "t_sync_conv_1": "t_sync_conv_1",
    "t_sync_rad_2": "t_sync_rad_2",
    "t_sync_conv_2": "t_sync_conv_2",
    "mass_transfer_case": None,
    "nearest_neighbour_distance": None,
}

DEFAULT_TRANSLATED_KEYS = (
    'age',
    'mass',
    'mdot',
    'inertia',
    'conv_mx1_top_r',
    'conv_mx1_bot_r',
    'surface_h1',
    'center_h1',
    'mass_conv_reg_fortides',
    'thickness_conv_reg_fortides',
    'radius_conv_reg_fortides',
    'log_Teff',
    'surface_he3',
    'surface_he4',
    'center_he4',
    'avg_c_in_c_core',
    'log_LH',
    'log_LHe',
    'log_LZ',
    'log_Lnuc',
    'c12_c12',
    'center_c12',
    'he_core_mass',
    'log_L',
    'log_R',
    'c_core_mass',
    'o_core_mass',
    'co_core_mass',
    'c_core_radius',
    'o_core_radius',
    'co_core_radius',
    'spin_parameter',
    'log_total_angular_momentum',
    'center_n14',
    'center_o16',
    'surface_n14',
    'surface_o16',
    'conv_env_top_mass',
    'conv_env_bot_mass',
    'conv_env_top_radius',
    'conv_env_bot_radius',
    'conv_env_turnover_time_g',
    'conv_env_turnover_time_l_b',
    'conv_env_turnover_time_l_t',
    'envelope_binding_energy',
    'lambda_CE_1cent',
    'lambda_CE_10cent',
    'lambda_CE_30cent',
    'lambda_CE_pure_He_star_10cent',
    'center_gamma'
)

KEYS_POSITIVE = (
    'mass_conv_reg_fortides',
    'thickness_conv_reg_fortides',
    'radius_conv_reg_fortides'
)

DEFAULT_PROFILE_KEYS = (
    'radius',
    'mass',
    'logRho',
    'energy',
    'x_mass_fraction_H',
    'y_mass_fraction_He',
    'z_mass_fraction_metals',
    'neutral_fraction_H',
    'neutral_fraction_He',
    'avg_charge_He'
)

class track_matcher:
    """
        This class contains the functionality to match binary star components to
    single star models for detached evolution. Typically, this is so that the 
    binary star can continue evolving in a detached (non-interacting) state.

        Several matching methods may be used, and metrics (physical quantities) 
    that are used to determine the quality of the match may be customized. By 
    default, the metrics are as follows:

      Initial stellar matching metrics
      ========================================
        
      MS:          total mass,                      center X, radius, He core mass
      post-MS:     total mass,                      center Y, radius, He core mass
      stripped He: He core mass (total mass proxy), center Y, radius

        If an initial match can not be found, alternative sets of metrics are 
    used. For example, stars after mass transfer could swell up so that log_R is 
    not appropriate for matching. Lists for HMS and postMS stars drop radius as a 
    matching metric in these alternatives:
    
      Alternative stellar matching metrics
      ========================================
        
      MS:          total mass,                      center X, He core mass
      post-MS:     total mass,                      center Y, He core mass
      stripped He: He core mass (total mass proxy), center Y, radius
    
    The available matching methods are:

        "root":     This method invokes a modified Powell's method to 
                    minimize the difference between stellar mass and 
                    central hydrogen abundance (if the star is H-rich), 
                    or He core mass (if the star is He-rich). The function
                    used is scipy.optimize.root with the 'hybr' method.

        "minimize": This method minimizes the Euclidean distance between 
                    several matching metrics (the same mentioned above). 
                    This is the default matching method.
        
        Stellar tracks in either the H- or He-rich single star grids will 
    be searched until a match is found within an acceptable tolerance. If 
    no match can be found, the binary star will be marked as a failed 
    binary as its evolution can progress no further.

        If an acceptable match is found, then we were able to find a star 
    within the single star grids that suitably represents the given star(s) 
    at their current point in evolution. This class will return interpolator 
    objects that may be used to calculate quantities along the matched stellar 
    track so that the binary's stars may be evolved further, such as through 
    detached evolution.
            
             
    Parameters
    ----------
    path : str
        Path to the directory that contains POSYDON data HDF5 files.
    verbose : bool
        True if we want to print stuff.

    Attributes
    ----------
    KEYS : list[str]
           Contains valid keywords which are used to extract quantities from 
           the grids.

    matching_method : str
        Method to find the best match between a star from a previous step and a
        point in a single MIST-like stellar track. Options "root" (which tries
        to find a root of two matching quantities, and it is possible to not
        achieve it) or "minimize" (minimizes the sum of squares of differences
        of various quantities between the previous step and the track).
           
    grid : GRIDInterpolator object
           Object to interpolate between the time-series (i.e., along the 
           evolutionary track) in the h5 grid.

    initial_mass : list[float]
            Contains the initial masses of the stars in the grid.

    Note
    ----
    A matching between the properties of the star, and the h5 tracks are
    required. In the "root" solver matching_method, if the root solver fails
    then the evolution will immediately end and the binary state will be
    tagged with "Root solver failed". In the "minimize" matching_method, we
    minimize the sum of squares of differences of various quantities between
    the previous step and the h5 track.

    """

    def __init__(
            self,
            KEYS,
            grid_name_Hrich,
            grid_name_strippedHe,
            path=PATH_TO_POSYDON_DATA,
            metallicity=None,
            matching_method="minimize",
            initial_mass=None,
            rootm=None,
            verbose=False,
            list_for_matching_HMS=None,
            list_for_matching_postMS=None,
            list_for_matching_HeStar=None,
    ):

        # MESA history column names used as matching metrics
        # TODO: should this be singlestar.STARPROPERTIES? An 
        #       error is thrown when (possibly user defined) 
        #       matching metrics don't exist in this array. 
        #       That's not very flexible...
        self.root_keys = np.array(
            [
                "age",
                "mass",
                "he_core_mass",
                "center_h1",
                "center_he4",
                "surface_he4",
                "surface_h1",
                "log_R",
                "center_c12"
            ]
        )

        # ==================================================================================

        self.metallicity = metallicity #convert_metallicity_to_string(metallicity)
        self.matching_method = matching_method

        self.initial_mass = initial_mass
        self.rootm = rootm
        self.verbose = verbose

        self.list_for_matching_HMS = list_for_matching_HMS
        self.list_for_matching_postMS = list_for_matching_postMS
        self.list_for_matching_HeStar = list_for_matching_HeStar

        # mapping a combination of (key, htrack, method) to a pre-trained
        # DataScaler instance, created the first time it is requested
        self.stored_scalers = {}

        # these are the KEYS read from POSYDON h5 grid files (after translating
        # them to the appropriate columns)
        self.KEYS = KEYS #DEFAULT_TRANSLATED_KEYS
        self.KEYS_POSITIVE = KEYS_POSITIVE

        self.profile_keys = DEFAULT_PROFILE_KEYS

        # keys for the final value interpolation
        self.final_keys = (
            'avg_c_in_c_core_at_He_depletion',
            'co_core_mass_at_He_depletion',
            'm_core_CE_1cent',
            'm_core_CE_10cent',
            'm_core_CE_30cent',
            'm_core_CE_pure_He_star_10cent',
            'r_core_CE_1cent',
            'r_core_CE_10cent',
            'r_core_CE_30cent',
            'r_core_CE_pure_He_star_10cent'
        )

        # keys for the star profile interpolation
        self.profile_keys = DEFAULT_PROFILE_KEYS

        # should grids just get passed to this?
        if grid_name_Hrich is None:
            grid_name_Hrich = os.path.join('single_HMS', self.metallicity+'_Zsun.h5')
        self.grid_Hrich = GRIDInterpolator(os.path.join(path, grid_name_Hrich))

        if grid_name_strippedHe is None:
            grid_name_strippedHe = os.path.join('single_HeMS', self.metallicity+'_Zsun.h5')
        self.grid_strippedHe = GRIDInterpolator(os.path.join(path, grid_name_strippedHe))

        # ==================================================================================

        # Initialize the matching lists:

        # min/max ranges of initial masses for each grid
        m_min_H = np.min(self.grid_Hrich.grid_mass)
        m_max_H = np.max(self.grid_Hrich.grid_mass)
        m_min_He = np.min(self.grid_strippedHe.grid_mass)
        m_max_He = np.max(self.grid_strippedHe.grid_mass)

        # Stellar parameter matching metrics
        # ========================================
        # At first, we try to match based on...

        if self.list_for_matching_HMS is None:
            self.list_for_matching_HMS = [
                ["mass", "center_h1", "log_R", "he_core_mass"],
                [20.0, 1.0, 2.0, 10.0],
                ["log_min_max", "min_max", "min_max", "min_max"],
                [m_min_H, m_max_H], [0, None]
            ]
        if self.list_for_matching_postMS is None:
            self.list_for_matching_postMS = [
                ["mass", "center_he4", "log_R", "he_core_mass"],
                [20.0, 1.0, 2.0, 10.0],
                ["log_min_max", "min_max", "min_max", "min_max"],
                [m_min_H, m_max_H], [0, None]
            ]
        if self.list_for_matching_HeStar is None:
            self.list_for_matching_HeStar = [
                ["he_core_mass", "center_he4", "log_R"],
                [10.0, 1.0, 2.0],
                ["min_max", "min_max", "min_max"],
                [m_min_He, m_max_He], [0, None]
            ]

        # Stellar parameter lists for alternative matching metrics
        # These are used in the event that an initial match can not be found
        self.list_for_matching_HMS_alternative = [
            ["mass", "center_h1", "he_core_mass"],
            [20.0, 1.0, 10.0],
            ["log_min_max", "min_max", "min_max"],
            [m_min_H, m_max_H], [0, None]
        ]
        self.list_for_matching_postMS_alternative = [
            ["mass", "center_h1", "he_core_mass"],
            [20.0, 1.0, 10.0],
            ["log_min_max", "min_max", "min_max"],
            [m_min_H, m_max_H], [0, None]
        ]
        self.list_for_matching_HeStar_alternative = [
            ["he_core_mass", "center_he4", "log_R"],
            [10.0, 1.0, 2.0],
            ["min_max", "min_max", "min_max"],
            [m_min_He, m_max_He], [0, None]
        ]

    def square_difference(self, x, htrack, attr_names, attr_vals, attr_scalers):
        """
            Compute the square difference between values along a single star 
        evolution track and a given set of values. To find a 'good match', 
        between the given values and the single star track, we ultimately seek 
        to minimize these differences. This is the function used by the 
        minimization algoritms when using matching_method=`minimize`. 

        Parameters
        ----------
        x : list[float]
            Contains a given initial mass `m0` and time `t`. These are used to get 
            the square difference between a track of initial mass `m0` at time `t` 
            and the given star values from prior evolution. 
        
        htrack : bool
            Set True to search the single star H-rich grids, or False to search 
            the He-rich grids.

        attr_names : list[str]
            Names of SingleStar object attributes that this will calculate the 
            square differences of.

        attr_vals : list[float]
            Associated values of provided SingleStar object attribute names. 
            These values should correspond to the prior point in evolution from 
            which we are making the stellar track match.

        scalers : list[DataScaler object]
            DataScaler objects that have been trained previously. These are used 
            to scale given star attributes to the range (0, 1).

        Returns
        -------
        result : float
            The square difference between a single star track with initial mass and 
            age taken from `x` and the given attributes of a SingleStar object.

        """

        result = 0.0
        for attr_name, attr_val, attr_scaler in zip(attr_names, attr_vals, attr_scalers):
            
            scaled_grid_track_val = attr_scaler.transform(self.get_track_val(attr_name, 
                                                                             htrack, *x))
            
            scaled_star_val = attr_scaler.transform(attr_val)

            if attr_name in MATCHING_WITH_RELATIVE_DIFFERENCE:
                result += ((scaled_grid_track_val - scaled_star_val)
                           / scaled_star_val) ** 2
            else:
                result += (scaled_grid_track_val - scaled_star_val) ** 2

        return result
    
    def get_root0(self, attr_names, attr_vals, htrack, rescale_facs=None):
        """
            Get the stellar evolution track in the single star grid with values 
        closest to the requested ones. This calculates the difference in stellar 
        properties between a given track and those in the single star grids. It 
        then returns the mass of and age along the track where the minimum 
        difference occurs.

        Parameters
        ----------
        attr_names : list[str]
            Contains the keys of the requested specific quantities that will be
            matched in the single star track.

        attr_vals : list[float]
            Contains the latest values (from a previous POSYDON step) of the 
            quantities of "keys" in the POSYDON SingleStar object. This should 
            be the same length as "keys".

        htrack : bool
            Set True to search the single star H-rich grids, or False to search 
            the He-rich grids.
            
        rescale_facs : list[float]
            Contains normalization factors to be divided for rescaling attribute 
            values. This should be the same length as attr_names. 

        Returns
        -------
        m0 : float
            Mass (in solar units) of the matched model. This is NaN if no match is 
            found.
    
        t : float
            Age (in years) of the matched model. This is NaN if no match is 
            found.

        """

        # set which grid to search based on htrack condition
        grid = self.grid_Hrich if htrack else self.grid_strippedHe

        # initial masses within grid (defined but never used? used in scale())
        self.initial_mass = grid.grid_mass

        max_track_length = 0
        # search across all initial masses and get max track length
        for mass in grid.grid_mass:
            max_track_length = max(max_track_length, len(grid.get("age", mass)))

        # intialize root matrix (DIM = [N(Mi), N(max_track_length), N(root_keys)])
        self.rootm = np.inf * np.ones((len(grid.grid_mass),
                                    max_track_length, len(self.root_keys)))
        
        # for each mass, get matching metrics and store in matrix
        for i, mass in enumerate(grid.grid_mass):
            for j, key in enumerate(self.root_keys):
                track = grid.get(key, mass)
                self.rootm[i, : len(track), j] = track

        # rescaling factors
        if rescale_facs is None:
            rescale_facs = np.ones_like(attr_names)
        else:
            rescale_facs = np.asanyarray(rescale_facs)

        star_attr_vals = np.asanyarray(attr_vals)
        # indices where given attribute names match available matching metrics
        idx = np.argmax(np.asanyarray(attr_names)[:, None] == self.root_keys, axis=1)
        # Slice out just the matching metric data for all stellar tracks
        # grid_attr_vals now has shape (N(Mi), N(max_track_len), N(matching_metrics))
        grid_attr_vals = self.rootm[:, :, idx]
        # For all stellar tracks in grid:
        # Take difference btwn. track and given star values, rescale, and take (Frobenius) norm
        grid_diff = np.linalg.norm((grid_attr_vals - star_attr_vals[None, None, :]) / rescale_facs[None, None, :], axis=-1)

        # indices where minimum difference occurs (i.e., of closest matching track)
        # This contains the indices at position of (min diff mass, min diff age)
        min_diff_inds = np.unravel_index(grid_diff.argmin(), grid_attr_vals.shape[:-1])

        mindiff_mass_i = min_diff_inds[0]
        mindiff_age_i = min_diff_inds[1]

        # time series and initial mass corresp. to track w/ minimum difference
        m0 = grid.grid_mass[mindiff_mass_i]
        t = self.rootm[mindiff_mass_i][mindiff_age_i][np.argmax("age" == self.root_keys)]

        return m0, t
    
    def get_track_val(self, key, htrack, m0, t):
        """
        
            Return a single value of a stellar property from the 
        interpolated time-series along a requested stellar track 
        of mass `m0` at an age of `t`.

        Parameters
        ----------
        key : str
            Keyword of the desired quantity.

        m0 : float
            The initial mass of the desired stellar track.

        t : float
            The desired age along the stellar track.

        Returns
        -------
        val : float
            The value of the desired quantity from a stellar track of
            initial mass `m0` at the time `t`.

        """
        # htrack as a boolean determines whether H or He grid is used
        if htrack:
            grid = self.grid_Hrich
        else:
            grid = self.grid_strippedHe
        try:
            x = grid.get("age", m0)
            y = grid.get(key, m0)
        except ValueError:
            return np.array(t) * np.nan
        try:
            val = np.interp(t, x, y, left=1e99, right=1e99)
        except ValueError:
            i_bad = [None]
            while len(i_bad):
                i_bad = np.where(np.diff(x) <= 0)[0]
                x = np.delete(x, i_bad)
                y = np.delete(y, i_bad)
            val = np.interp(t, x, y)

        return val
    
    def scale(self, attr_name, htrack, scaler_method):
        """
        
            Normalize quantities in the single star grids to (0,1).

        Parameters
        ----------
        attr_name : str
            Keyword of the requested quantity.

        htrack : bool
            A boolean that specifies whether the star would be found in the 
            hydrogen rich single star grid or not (in which case it is
            matched to the helium rich single star grid).

        scaler_method : str
            Scaling method in the DataScaler class. See 
            posydon.interpolation.data_scaling.DataScaler().fit() for more 
            details.

        Returns
        -------
        scaler : DataScaler object
            This is a DataScaler object, trained to rescale the requested  
            attribute to the range (0, 1).

        """
        # TODO: why this self.grid? Why not local variable. Should this affect
        # the whole detached_step instance?

        # collect all options for the scaler
        scaler_options = (attr_name, htrack, scaler_method)

        # find if the scaler has already been fitted and return it if so...
        scaler = self.stored_scalers.get(scaler_options, None)
        if scaler is not None:
            return scaler

        # ...if not, fit a new scaler, and store it for later use
        grid = self.grid_Hrich if htrack else self.grid_strippedHe
        self.initial_mass = grid.grid_mass
        all_attributes = []

        for mass in self.initial_mass:
            for i in grid.get(attr_name, mass):
                all_attributes.append(i)

        all_attributes = np.array(all_attributes)
        scaler = DataScaler()
        scaler.fit(all_attributes, method=scaler_method, lower=0.0, upper=1.0)
        self.stored_scalers[scaler_options] = scaler

        return scaler

    def match_to_single_star(self, star, htrack):
        """
            Get the track in the grid that matches the time and mass of a star,
        that has typically undergone prior binary star evolution. A match is made 
        according to a given algorithm that minimizes the difference amongst several 
        physical properties, e.g., mass, central hydrogen abundance, radius, and 
        core helium mass, depending on the type of star being matched. However, these 
        properties may also be customized by the user.

        Parameters
        ----------
        star : SingleStar object
            A single star object that contains the star's properties.
        
        htrack: bool
            A boolean that specifies whether the star would be found in the 
            hydrogen rich single star grid or not (in which case it is
            matched to the helium rich single star grid).

        Returns
        -------
        m0: float
            Mass (in solar units) of the matched model
        
        t0: float
            Age (in years) of the matched model

        htrack: bool
            This has the same meaning as the given htrack, but the value
            may change during the course of matching. In the event that 
            a match can not be found, an He or post-MS star may be 
            alternatively matched to the H or He grid, in spite of the 
            star's actual state, as a last ditch effort to find a match.

        Warns
        -----
        EvolutionWarning
            If attempting to match an He-star with an H-rich grid or post-MS star 
            with a stripped-He grid. This can happen if an initial matching attempt 
            fails; alternative grids are checked for a match in such cases.

        Raises
        ------
        AttributeError
            If a matching parameter does not exist in the single star 
            grid options
        
        NumericalError
            If SciPy numerical differentiation occured outside boundary 
            while matching to single star track.

        MatchingError
            If the initial mass to be matched is out of the single star grid 
            range, thereby preventing a possible match.

        """
        def get_attr_values(attr_names, star):

            """
                Given a set of attribute names, gets attribute values 
            from a given SingleStar object. These values correspond 
            to the last evolution step of the SingleStar object.

            Parameters
            ----------
            attr_names: list[str]
                This list contains strings that are the names of 
                data columns (attributes) to be used as matching 
                metrics.

            star: SingleStar object
                This is a SingleStar object, typically representing 
                a member of a binary star system.

            Returns
            -------
            attr_values: list[float]
                This is a list of the values associated with the 
                provided attribute names at the last evolution 
                step experienced by `star`.

            """
            
            attr_values = []

            for attr_name in attr_names:
                attr_values.append(getattr(star, attr_name))

            return attr_values
        
        def sq_diff_function(x):

            """
                This is a wrapper that calls a function to calculate the 
            square difference between specified star attributes for matching.
            When matching_method=`minimize`, this is the function that is 
            used by the minimization algorithm.

            Parameters
            ----------
            x: list[float]
                A list containing the initial mass and age of a stellar track, 
                which will be used to get values along a track and calculate the 
                square difference between those values and the given star 
                values that we are trying to find a match to.

            """

            return self.square_difference(x, htrack=htrack, 
                                          attr_names=match_attr_names,
                                          attr_vals=match_attr_vals,
                                          attr_scalers=scalers)
        
        def get_match_attr_props(list_for_matching):

            """
               This unpacks a list_for_matching which has the following 
            structure:

                list_for_matching = [[matching attr. names], [rescale_factors],
                                    [scaling method], [mass_bnds], [age_bnds]]
            
               This also trains DataScaler objects for each provided matching 
            attribute, such that the values are scaled to the range (0, 1).

            Parameters
            ----------
            list_for_matching: list
                This is a list that contains sublists of types, str, float, 
                str, float, float. The first sublist holds the names for 
                star attributes that will be used to quantify a match. The 
                second holds rescaling factors for those attributes. The 
                third holds DataScaler scaling methods (see 
                posydon.interpolation.data_scaling.DataScaler.fit() for more). 
                The fourth and fifth hold bounds for the minimization algorithms 
                on the stellar mass and age.

            Returns
            -------
            match_attr_names: list[str]
                The names of the attributes that will be used for matching. These 
                should be SingleStar STARPROPERTIES keys.
            
            rescale_facs: list[float]
                Contains normalization factors to be used for rescaling match attribtue 
                values. This should be the same length as match_attr_names.
                

            bnds: list[list[float], list[float]]
                Lower and upper bounds on initial stellar mass and age to be used 
                in minimization algorithms.
            
            scalers: list[DataScaler object]
                DataScaler objects for each attribute, trained to rescale quantities 
                to the range (0, 1).
                
            """

            match_attr_names = list_for_matching[0]
            rescale_facs = list_for_matching[1]
            scaler_methods = list_for_matching[2]
            bnds = list_for_matching[3:]

            if self.verbose:
                print("Matching parameters and their normalizations:\n", 
                      match_attr_names, rescale_facs)
            
            # get (or train and get) scalers for attributes
            # attributes are scaled to range (0, 1)
            scalers = []
            for attr_name, scaler_method in zip(match_attr_names, scaler_methods):
                scaler_of_attribute = scale(attr_name, htrack, scaler_method)
                scalers.append(scaler_of_attribute)

            return match_attr_names, rescale_facs, bnds, scalers

            
        # setting whether to match to track from H-rich or stripped He star grid
        if htrack:
            self.grid = self.grid_Hrich
        else:
            self.grid = self.grid_strippedHe

        get_root0 = self.get_root0
        get_track_val = self.get_track_val
        matching_method = self.matching_method
        scale = self.scale

        initial_track_vals = None
        
        matching_tolerance = 1e-2
        matching_tolerance_hard = 1e-1

        if self.verbose:
            print(f"\nMatching process started in detached step for {star.state} star "
                  f"with matching method = {matching_method}")

        # matching via root method (ultimately modified Powell's method)
        if matching_method == "root":

            # if the star can be considered an HMS star
            if star.state in LIST_ACCEPTABLE_STATES_FOR_HMS:

                # get initial guess for matching via search for closest track
                x0 = get_root0(["center_h1", "mass"],
                               [star.center_h1, star.mass],
                               htrack, rescale_facs=[0.7, 300])
                
                # Using modified Powell method to find match starting from time series from above (in x0)
                # using difference of center X and stellar mass to find match
                sol = root(
                            lambda x: [
                                        get_track_val("center_h1", htrack, *x) - star.center_h1,
                                        get_track_val("mass", htrack, *x) - star.mass
                                      ],
                            x0, method="hybr"
                          )
                
            # or if not an HMS star
            else:

                # same as above but using He core mass and total mass
                x0 = get_root0(
                               ["he_core_mass", "mass"],
                               [star.he_core_mass, star.mass],
                               htrack,
                               rescale_facs=[11, 300]
                              )
                
                sol = root(
                            lambda x: [
                                        get_track_val("he_core_mass", htrack, *x) - star.he_core_mass,
                                        get_track_val("mass", htrack, *x) - star.mass],
                            x0, method="hybr"
                          )
                
            # if optimizer failed for some reason set solution as NaN
            if not sol.success or sol.x[1] < 0:
                initial_track_vals = (np.nan, np.nan)
            else:
                initial_track_vals = sol.x

        # 1st attempt: match using minimize method (ultimately Newton or Powell method)
        elif matching_method == "minimize":
            
            # set matching metrics based on star state
            if star.state in LIST_ACCEPTABLE_STATES_FOR_HMS:
                list_for_matching = self.list_for_matching_HMS
            elif star.state in LIST_ACCEPTABLE_STATES_FOR_postMS:
                list_for_matching = self.list_for_matching_postMS
            elif star.state in LIST_ACCEPTABLE_STATES_FOR_HeStar:
                list_for_matching = self.list_for_matching_HeStar
            
            match_attr_names, rescale_facs, bnds, scalers = get_match_attr_props(list_for_matching)
            
            # MESA labels are what we will match with
            for param_name in match_attr_names:
                if param_name not in self.root_keys:
                    raise AttributeError(f"Expected matching parameter {param_name} not "+\
                                         f"added in root_keys list: {self.root_keys}")
                
            match_attr_vals = get_attr_values(match_attr_names, star)

            # Get closest matching point along track single star grids. x0 contains the 
            # corresponding initial mass and age
            x0 = get_root0(match_attr_names, match_attr_vals, htrack, rescale_facs=rescale_facs)

            # Match using Newton method w/ x0 time series as guess
            # Minimizing Euclidean distance (1st w/ Newton's method)
            min_method = "TNC"
            sol = minimize(sq_diff_function, x0, method=min_method, bounds=bnds)

            # save initial matching solution as best solution so far
            best_sol = sol

            ## Alternative matching attempts if default matching fails!
            # 2nd attempt: use a different minimization method
            if (np.abs(best_sol.fun) > matching_tolerance or not best_sol.success):

                # try alternative minimization method, Powell's
                min_method = "Powell"

                if self.verbose:
                    if (np.abs(best_sol.fun) > matching_tolerance):
                        print (f"\n\nInitial matching attempt FAILED:"
                               f"\nSolution exceeds tolerance: {np.abs(best_sol.fun)} > {matching_tolerance}")
                    if (not best_sol.success):
                        print (f"Initial matching attempt FAILED:"
                               f"\nOptimizer failed (sol.success = {best_sol.success})"
                               f"\nOptimizer termination reason: {best_sol.message}")                        

                    print("\nAlternative matching started (1st attempt)")
                    print(f"(Trying alternative minimization method: {min_method})")
                
                # minimize w/ modified Powell's method
                sol = minimize(sq_diff_function, x0, method=min_method)

                ## if alternative matching has a better solution, make it the new best solution
                if (np.abs(sol.fun) < np.abs(best_sol.fun) and sol.success):
                    best_sol = sol

            # if 2nd fails, 3rd attempt: use alternative matching parameters
            if (np.abs(best_sol.fun) > matching_tolerance or not best_sol.success):

                # back to Newton's method
                min_method = "TNC" 
                   
                if self.verbose:
                    if (np.abs(best_sol.fun) > matching_tolerance):
                        print (f"Alternative matching (1st attempt) FAILED:"
                               f"\nSolution exceeds tolerance: {np.abs(best_sol.fun)} > {matching_tolerance}")
                    if (not best_sol.success):
                        print (f"Alternative matching (1st attempt) FAILED:"
                               f"\nOptimizer failed (sol.success = {best_sol.success})"
                               f"\nOptimizer termination reason: {best_sol.message}")

                    print("\nAlternative matching started (2nd attempt)")
                    print("(Now trying to match with alternative parameters)")     
                      
                # set alternative matching metrics based on star state
                if star.state in LIST_ACCEPTABLE_STATES_FOR_HMS:
                    list_for_matching = self.list_for_matching_HMS_alternative
                elif star.state in LIST_ACCEPTABLE_STATES_FOR_postMS:
                    list_for_matching = (self.list_for_matching_postMS_alternative)
                elif star.state in LIST_ACCEPTABLE_STATES_FOR_HeStar:
                    list_for_matching = (self.list_for_matching_HeStar_alternative)                
                    
                match_attr_names, rescale_facs, bnds, scalers = get_match_attr_props(list_for_matching)
                match_attr_vals = get_attr_values(match_attr_names, star)

                # Minimize w/ Newton's method using alternative metrics
                x0 = get_root0(match_attr_names, match_attr_vals, htrack, rescale_facs=rescale_facs)
                sol = minimize(sq_diff_function, x0, method=min_method, bounds=bnds)

                if (np.abs(sol.fun) < np.abs(best_sol.fun) and sol.success):
                    best_sol = sol

            # 4th attempt: match an He-star with an H-rich grid, or vice versa (not applicable for HMS stars)
            if (np.abs(best_sol.fun) > matching_tolerance or not best_sol.success):

                # Using Newton's minimization method for 4th attempt
                min_method = "TNC"

                if self.verbose:
                    if (np.abs(best_sol.fun) > matching_tolerance):
                        print (f"Alternative matching (2nd attempt) FAILED:"
                               f"\nSolution exceeds tolerance: {np.abs(best_sol.fun)} > {matching_tolerance}")
                    if (not best_sol.success):
                        print (f"Alternative matching (2nd attempt) FAILED:"
                               f"\nOptimizer failed (sol.success = {best_sol.success})"
                               f"\nOptimizer termination reason: {best_sol.message}")  
                
                # if post-MS or stripped He star
                if (star.state in LIST_ACCEPTABLE_STATES_FOR_HeStar
                    or star.state in LIST_ACCEPTABLE_STATES_FOR_postMS):

                    if self.verbose:
                        print("\nAlternative matching started (3rd attempt)")
                        print("(Now trying to match He-star or post-MS star to a different grid)")

                    Pwarn("Attempting to match an He-star with an H-rich grid or post-MS star with a"
                          " stripped-He grid", "EvolutionWarning")
                        
                    if star.state in LIST_ACCEPTABLE_STATES_FOR_HeStar:
                        new_htrack = True
                        list_for_matching = self.list_for_matching_HeStar

                    elif star.state in LIST_ACCEPTABLE_STATES_FOR_postMS:
                        new_htrack = False
                        list_for_matching = self.list_for_matching_postMS

                    match_attr_names, rescale_facs, bnds, scalers = get_match_attr_props(list_for_matching)

                    for param_name in match_attr_names:
                        if param_name not in self.root_keys:
                            raise AttributeError(f"Expected matching parameter {param_name} not "+\
                                                f"added in root_keys list: {self.root_keys}")
                        
                    match_attr_vals = get_attr_values(match_attr_names, star)
                    x0 = get_root0(match_attr_names, match_attr_vals, new_htrack, rescale_facs=rescale_facs)

                    try:
                        # minimize Euclidean dist. w/ Newton's method (TNC)
                        sol = minimize(sq_diff_function, x0, method=min_method, bounds=bnds)

                        if (np.abs(sol.fun) < np.abs(best_sol.fun) and sol.success):
                            best_sol = sol
                            htrack = new_htrack

                        if self.verbose:
                            print (f"Alternative matching (3rd attempt) completed:"
                                   f"\nBest solution: {np.abs(best_sol.fun)} (tol = {matching_tolerance})"
                                   f"\nsol.success = {best_sol.success}")
                    except:
                        raise NumericalError("SciPy numerical differentiation occured outside boundary "
                                             "while matching to single star track")

            # if matching is still not successful, set result to NaN:
            if (np.abs(best_sol.fun) > matching_tolerance_hard or not best_sol.success):
                if self.verbose:
                    print("\nFinal matching result with relaxed tolerance: FAILED")#,
                          #np.abs(best_sol.fun), ">", matching_tolerance_hard)
                    if (np.abs(best_sol.fun) > matching_tolerance_hard):
                        print ("\nSolution exceeds hard tolerance: "+\
                               f"{np.abs(best_sol.fun)} > {matching_tolerance_hard}")
                    if (not best_sol.success):
                        print (f"\nOptimizer failed, sol.success = {best_sol.success}"
                               f"\nOptimizer termination reason: {best_sol.message}") 

                initial_track_vals = (np.nan, np.nan)
                
            # or else we found a solution
            else:
                if self.verbose:
                    print("\nFinal matching result: SUCCESS. "
                         f"\nBest solution within hard tolerance: "
                         f"{np.abs(best_sol.fun):.8f}", "<", matching_tolerance_hard)
                    
                initial_track_vals = best_sol.x

        if self.verbose:
            # successful match
            if not np.isnan(initial_track_vals[0]):
                val_names = ["  ", "mass", "log_R", "center_h1", "surface_h1", "he_core_mass", "center_he4", "surface_he4", 
                             "center_c12"]

                initial_vals = [
                                "initial values",
                                f'{star.mass:.3f}', 
                                f'{star.log_R:.3f}', 
                                f'{star.center_h1:.3f}', 
                                f'{star.surface_h1:.4f}',
                                f'{star.he_core_mass:.3f}',
                                f'{star.center_he4:.4f}', 
                                f'{star.surface_he4:.4f}',  
                                f'{star.center_c12:.4f}'
                               ]

                matched_vals = [
                                "matched values",
                                f'{self.get_track_val("mass", htrack, *best_sol.x):.3f}', 
                                f'{self.get_track_val("log_R", htrack, *best_sol.x):.3f}',
                                f'{self.get_track_val("center_h1", htrack, *best_sol.x):.3f}',
                                f'{self.get_track_val("surface_h1", htrack, *best_sol.x):.4f}',
                                f'{self.get_track_val("he_core_mass", htrack, *best_sol.x):.3f}',
                                f'{self.get_track_val("center_he4", htrack, *best_sol.x):.4f}',
                                f'{self.get_track_val("surface_he4", htrack, *best_sol.x):.4f}',    
                                f'{self.get_track_val("center_c12", htrack, *best_sol.x):.4f}'
                               ]
                
                percent_diff = [
                                "% difference",
                                f'{np.nan_to_num(100.0*(self.get_track_val("mass", htrack, *best_sol.x) - star.mass)/star.mass, 0):.1e}%', 
                                f'{np.nan_to_num(100.0*(self.get_track_val("log_R", htrack, *best_sol.x) - star.log_R)/star.log_R, 0):.1e}%',
                                f'{np.nan_to_num(100.0*(self.get_track_val("center_h1", htrack, *best_sol.x) - star.center_h1)/star.center_h1, 0):.1e}%',
                                f'{np.nan_to_num(100.0*(self.get_track_val("surface_h1", htrack, *best_sol.x) - star.surface_h1)/star.surface_h1, 0):.1e}%',
                                f'{np.nan_to_num(100.0*(self.get_track_val("he_core_mass", htrack, *best_sol.x) - star.he_core_mass)/star.he_core_mass, 0):.1e}%',
                                f'{np.nan_to_num(100.0*(self.get_track_val("center_he4", htrack, *best_sol.x) - star.center_he4)/star.center_he4, 0):.1e}%',
                                f'{np.nan_to_num(100.0*(self.get_track_val("surface_he4", htrack, *best_sol.x) - star.surface_he4)/star.surface_he4, 0):.1e}%',    
                                f'{np.nan_to_num(100.0*(self.get_track_val("center_c12", htrack, *best_sol.x) - star.center_c12)/star.center_c12, 0):.1e}%'
                               ]

                output_table = [val_names, initial_vals, matched_vals, percent_diff]

                print("\nMatching completed for", star.state, "star!\n")
                for row in output_table:
                    print("{:>14}  {:>9}  {:>9}  {:>9}  {:>10}  {:>12}  {:>10}  {:>11}  {:>10}".format(*row))                

            # failed match
            else:
                print(
                    "Matching was unsuccessful for star with properties: \n"
                    f'mass = {star.mass:.3f}, ',
                    f'log_R = {star.log_R:.3f}, ',
                    f'center_he4 = {star.center_he4:.4f}, ',
                    f'surface_he4 = {star.surface_he4:.4f}, ',
                    f'surface_h1 = {star.surface_h1:.4f}, ',
                    f'he_core_mass = {star.he_core_mass:.3f}, ',
                    f'center_c12 = {star.center_c12:.4f}'
                )

        m0 = initial_track_vals[0]
        t0 = initial_track_vals[1]

        return m0, t0, htrack
    

    def get_star_match_data(self, binary, star, copy_prev_m0=None, copy_prev_t0=None):
                """
                    Match a given component of a binary (i.e., a star) to a 
                single star model. This then creates and returns interpolator
                objects that may be used to calculate properties of the star
                as a function of time.

                    In the case of a compact object, radius, mdot, and Idot are 
                set to zero. One may use another star, e.g., the companion of 
                the compact object to provide an initial mass and age.

                Parameters
                ----------
                binary: BinaryStar object
                    A binary star object, containing the binary system's properties.

                star: SingleStar object
                    A single star object that contains the star's properties.

                copy_prev_m0: float
                    A mass value that may be copied from another star in the case
                    where the target star is a compact object

                copy_prev_t0: float
                    An age value that may be copied from another star in the case
                    where the target star is a compact object

                Return
                -------
                interp1d: dict
                    A dictionary of scipy.interpolate._cubic.PchipInterpolator 
                    objects used to calculate star properties corresponding to the 
                    matched single star track.

                """

                KEYS = self.KEYS
                KEYS_POSITIVE = self.KEYS_POSITIVE

                htrack = star.htrack
                co = star.co

                with np.errstate(all="ignore"):
                    # get the initial m0, t0 track
                    if binary.event == 'ZAMS' or binary.event == 'redirect_from_ZAMS':
                        # ZAMS stars in wide (non-mass exchaging binaries) that are
                        # directed to detached step at birth
                        m0, t0 = star.mass, 0
                    elif co:
                        m0, t0 = copy_prev_m0, copy_prev_t0
                    else:
                        t_before_matching = time.time()
                        # matching to single star grids
                        m0, t0, htrack = self.match_to_single_star(star, htrack)
                        t_after_matching = time.time()
                        
                        if self.verbose:
                            print(f"Matching duration: {t_after_matching-t_before_matching:.6g} sec\n")
                
                if pd.isna(m0) or pd.isna(t0):
                    return None, None, None
                
                if htrack:
                    self.grid = self.grid_Hrich
                else:
                    self.grid = self.grid_strippedHe
                
                # check if m0 is in the grid bounds
                if m0 < self.grid.grid_mass.min() or m0 > self.grid.grid_mass.max():
                    set_binary_to_failed(binary)
                    raise MatchingError(f"The mass {m0} is out of the single star grid range and "
                                        "cannot be matched to a track.")

                # get/interpolate track values for requested mass m0
                get_track = self.grid.get

                max_time = binary.properties.max_simulation_time
                assert max_time > 0.0, "max_time is non-positive"

                age = get_track("age", m0)
                t_max = age.max()  # max timelength of the track
                interp1d = dict()
                kvalue = dict()
                for key in KEYS[1:]:
                    kvalue[key] = get_track(key, m0)
                try:
                    for key in KEYS[1:]:
                        if key in KEYS_POSITIVE:
                            positive = True
                            interp1d[key] = PchipInterpolator2(age, kvalue[key], positive=positive)
                        else:
                            interp1d[key] = PchipInterpolator2(age, kvalue[key])
                except ValueError:
                    i_bad = [None]
                    while len(i_bad) != 0:
                        i_bad = np.where(np.diff(age) <= 0)[0]
                        age = np.delete(age, i_bad)
                        for key in KEYS[1:]:
                            kvalue[key] = np.delete(kvalue[key], i_bad)

                    for key in KEYS[1:]:
                        if key in KEYS_POSITIVE:
                            positive = True
                            interp1d[key] = PchipInterpolator2(age, kvalue[key], positive=positive)
                        else:
                            interp1d[key] = PchipInterpolator2(age, kvalue[key])

                interp1d["inertia"] = PchipInterpolator(
                    age, kvalue["inertia"] / (const.msol * const.rsol**2))
                interp1d["Idot"] = interp1d["inertia"].derivative()

                interp1d["conv_env_turnover_time_l_b"] = PchipInterpolator2(
                    age, kvalue['conv_env_turnover_time_l_b'] / const.secyer)

                interp1d["L"] = PchipInterpolator(age, 10 ** kvalue["log_L"])
                interp1d["R"] = PchipInterpolator(age, 10 ** kvalue["log_R"])
                interp1d["t_max"] = t_max
                interp1d["max_time"] = max_time
                interp1d["t0"] = t0
                interp1d["m0"] = m0

                if co:
                    kvalue["mass"] = np.zeros_like(kvalue["mass"]) + star.mass
                    kvalue["R"] = np.zeros_like(kvalue["log_R"])
                    kvalue["mdot"] = np.zeros_like(kvalue["mdot"])
                    interp1d["mass"] = PchipInterpolator(age, kvalue["mass"])
                    interp1d["R"] = PchipInterpolator(age, kvalue["R"])
                    interp1d["mdot"] = PchipInterpolator(age, kvalue["mdot"])
                    interp1d["Idot"] = PchipInterpolator(age, kvalue["mdot"])

                return interp1d, m0, t0

    def update_star_properties(self, star, htrack, m0, t0):

        """
            This updates a SingleStar object (`star`) with the 
        values from a single star track that has initial mass `m0` 
        and age `t0`. This can be used after matching finds the 
        closest `m0`, `t0` to update the SingleStar object with the 
        values of the best matching single star track. This will 
        not update the `log_total_angular_momentum` or 
        `surf_avg_omega` because the single star track is non-rotating 
        and those quantities are poorly defined.

        Parameters
        ----------
        star: SingleStar object
            A single star object that contains the star's properties.
        
        htrack: bool
            A boolean that specifies whether the star would be found in the 
            hydrogen rich single star grid or not (in which case it is
            matched to the helium rich single star grid).

        m0: float
            Initial stellar mass (in solar units) of the single star track 
            that we will grab values from and update `star` with.
        
        t0: float
            Initial age (in years) of the single star track 
            that we will grab values from and update `star` with.

        """

        for key in self.KEYS:
            
            # skip updating rotation rate quantities because
            # they're 0 or not defined in non-rotating tracks
            if key == "log_total_angular_momentum":
                continue
            elif key == "surf_avg_omega":
                continue
            else:
                new_val = self.get_track_val(key, htrack, m0, t0)
                
            setattr(star, key, new_val)
           

    def get_star_final_values(self, star, htrack, m0): 
        """
            This updates the final values of a SingleStar object,
        given an initial stellar mass `m0`, typically found from 
        matching to a single star track.

        Parameters
        ----------
        star: SingleStar object
            A single star object that contains the star's properties.
        
        htrack: bool
            A boolean that specifies whether the star would be found in the 
            hydrogen rich single star grid or not (in which case it is
            matched to the helium rich single star grid).

        m0: float
            Initial stellar mass (in solar units) of the single star track 
            that we will grab values from and update `star` with.

        """

        grid = self.grid_Hrich if htrack else self.grid_strippedHe
        get_final_values = grid.get_final_values

        for key in self.final_keys:
            setattr(star, key, get_final_values('S1_%s' % (key), m0))

    def get_star_profile(self, star, htrack, m0):
        """
            This updates the stellar profile of a SingleStar object,
        given an initial stellar mass `m0`, typically found from 
        matching to a single star track. The profile of the SingleStar 
        object is updated to become the profile of the (matched) single 
        star track.

        Parameters
        ----------
        star: SingleStar object
            A single star object that contains the star's properties.
        
        htrack: bool
            A boolean that specifies whether the star would be found in the 
            hydrogen rich single star grid or not (in which case it is
            matched to the helium rich single star grid).

        m0: float
            Initial stellar mass (in solar units) of the single star track 
            that we will grab values from and update `star` with.

        """

        grid = self.grid_Hrich if htrack else self.grid_strippedHe
        get_profile = grid.get_profile
        profile_new = np.array(get_profile('mass', m0)[1])

        for i in self.profile_keys:
            profile_new[i] = get_profile(i, m0)[0]
        profile_new['omega'] = star.surf_avg_omega

        star.profile = profile_new
