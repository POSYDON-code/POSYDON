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

class track_matcher:


    def __init__(
            self,
            KEYS,
            KEYS_POSITIVE,
            DEFAULT_PROFILE_KEYS,
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
        #self.dt = dt
        #self.n_o_steps_history = n_o_steps_history
        self.matching_method = matching_method
        #self.do_wind_loss = do_wind_loss
        #self.do_tides = do_tides
        #self.do_gravitational_radiation = do_gravitational_radiation
        #self.do_magnetic_braking = do_magnetic_braking
        #self.magnetic_braking_mode = magnetic_braking_mode
        #self.do_stellar_evolution_and_spin_from_winds = (
        #    do_stellar_evolution_and_spin_from_winds
        #)
        #self.RLO_orbit_at_orbit_with_same_am = RLO_orbit_at_orbit_with_same_am

        self.initial_mass = initial_mass
        self.rootm = rootm
        self.verbose = verbose

        self.list_for_matching_HMS = list_for_matching_HMS
        self.list_for_matching_postMS = list_for_matching_postMS
        self.list_for_matching_HeStar = list_for_matching_HeStar

        # mapping a combination of (key, htrack, method) to a pre-trained
        # DataScaler instance, created the first time it is requested
        self.stored_scalers = {}

        #if self.verbose:
        #    print(
        #        dt,
        #        n_o_steps_history,
        #        matching_method,
        #        do_wind_loss,
        #        do_tides,
        #        do_gravitational_radiation,
        #        do_magnetic_braking,
        #        magnetic_braking_mode,
        #        do_stellar_evolution_and_spin_from_winds)

        #self.translate = DEFAULT_TRANSLATION

        # these are the KEYS read from POSYDON h5 grid files (after translating
        # them to the appropriate columns)
        self.KEYS = KEYS #DEFAULT_TRANSLATED_KEYS
        self.KEYS_POSITIVE = KEYS_POSITIVE #(
            #'mass_conv_reg_fortides',
            #'thickness_conv_reg_fortides',
            #'radius_conv_reg_fortides'
        #)

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
        #self.grid_Hrich = grid_Hrich

        if grid_name_strippedHe is None:
            grid_name_strippedHe = os.path.join('single_HeMS', self.metallicity+'_Zsun.h5')
        self.grid_strippedHe = GRIDInterpolator(os.path.join(path, grid_name_strippedHe))
        #self.grid_strippedHe = grid_strippedHe

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
        #
        #     MS:          total mass,                      center X, radius, He core mass
        #     post-MS:     total mass,                      center Y, radius, He core mass
        #     stripped He: He core mass (total mass proxy), center Y, radius

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

        # e.g., stars after mass transfer could swell up so that log_R
        # is not appropriate for matching. HMS and postMS lists drop radius
        # as a matching metric in these alternatives

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

    def square_difference(self, x, htrack, mesa_labels, posydon_attributes, colscalers, scales):
        """Compute the square distance used for scaling."""
        result = 0.0
        for mesa_label, posy_attr, colscaler, scale_of_mesa_label in zip(
                 mesa_labels, posydon_attributes, colscalers, scales):
            single_track_value = scale_of_mesa_label.transform(
                self.get_track_val(mesa_label, htrack, *x))
            posydon_value = scale_of_mesa_label.transform(posy_attr)
            if mesa_label in MATCHING_WITH_RELATIVE_DIFFERENCE:
                result += ((single_track_value - posydon_value)
                           / posydon_value) ** 2
            else:
                result += (single_track_value - posydon_value) ** 2

        return result
    
    def get_root0(self, keys, x, htrack, rs=None):
        """
            Get the track in the grid with values closest to the requested ones.

            Parameters
            ----------
            keys : list of str
                   Contains the keys of the required specific quantities that will be
                   matched in the MIST-like track.

            x :    list of floats, of same length as "keys"
                   Contains the latest values (from a previous POSYDON step) of the
                   quantities of "keys" in the POSYDON SingleStar object.

            htrack : bool
                     Set True to search the H-rich grids, or False to search the 
                     stripped He-star grids.
               
            rs :    list of floats, same length as "keys"
                    Contains normalization factors to be divided for rescaling
                    x values.

            Returns
            -------
            list of 2 float values
                Contains the associated initial mass (in solar units) and the time
                (in years) such that the time-series of the `keys` at that time has
                the closest values to `x`. These will become m0, t0 for the later
                integration during the detached binary evolution.
                If there is no match then NaNs will be returned instead.

        """

        # set which grid to search based on htrack condition
        grid = self.grid_Hrich if htrack else self.grid_strippedHe

        # initial masses within grid (defined but never used?)
        self.initial_mass = grid.grid_mass

        # this will track the max track length
        max_track_length = 0

        # search across all initial masses and get max track length
        for mass in grid.grid_mass:
            max_track_length = max(max_track_length, len(grid.get("age", mass)))

        # intialize root matrix (DIM = [N(Mi), N(max_track_length), N(matching_metrics)])
        self.rootm = np.inf * np.ones((len(grid.grid_mass),
                                    max_track_length, len(self.root_keys)))
        
        # for each mass, get matching metrics and store in matrix
        for i, mass in enumerate(grid.grid_mass):
            for j, key in enumerate(self.root_keys):

                # get track metric (key) for given initial mass
                track = grid.get(key, mass)

                # store metrics in grid matrix
                self.rootm[i, : len(track), j] = track

        # rescaling factors
        if rs is None:
            rs = np.ones_like(keys)
        else:
            rs = np.asanyarray(rs)

        # we're matching to the values stored in x
        x = np.asanyarray(x)
        # single out the specified matching metrics supplied in keys
        idx = np.argmax(np.asanyarray(keys)[:, None] == self.root_keys, axis=1)
        # get all masses from grid, considering only supplied keys
        X = self.rootm[:, :, idx]
        # take difference, rescale, and take (Frobenius) norm
        d = np.linalg.norm((X - x[None, None, :]) / rs[None, None, :], axis=-1)

        # indices where minimum difference occurs (i.e., of closest matching track)
        idx = np.unravel_index(d.argmin(), X.shape[:-1])

        # time series and initial mass corresp. to track w/ minimum difference
        t = self.rootm[idx][np.argmax("age" == self.root_keys)]
        m0 = grid.grid_mass[idx[0]]

        return m0, t
    
    def get_track_val(self, key, htrack, m0, t):
        """Return a single value from the interpolated time-series.

        Parameters
        ----------
        key : str
            Keyword of the required quantity.
        m0 : float
            The associated initial mass of the required quantity.
        t  : float
            The required time in the time-series.

        Returns
        -------
        float
            The value of the quantity `key` from a MIST-like track of
            initial mass `m0` at the time `t0`.

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
    
    def scale(self, key, htrack, method):
        """Nomarlize quantities in the single star grids to (0,1).

        Parameters
        ----------
        key : str
            Keyword of the required quantity.
        method : str
            Scalling method in the data normalization class

        Returns
        -------
        class
            Data normalization class

        """
        # TODO: why this self.grid? Why not local variable. Should this affect
        # the whole detached_step instance?

        # collect all options for the scaler
        scaler_options = (key, htrack, method)

        # find if the scaler has already been fitted and return it if so...
        scaler = self.stored_scalers.get(scaler_options, None)
        if scaler is not None:
            return scaler

        # ... if not, fit a new scaler, and store it for later use
        grid = self.grid_Hrich if htrack else self.grid_strippedHe
        self.initial_mass = grid.grid_mass
        all_attributes = []
        for mass in self.initial_mass:
            for i in grid.get(key, mass):
                all_attributes.append(i)
        all_attributes = np.array(all_attributes)
        scaler = DataScaler()
        scaler.fit(all_attributes, method=method, lower=0.0, upper=1.0)
        self.stored_scalers[scaler_options] = scaler
        return scaler

    def match_to_single_star(self, star, htrack):
        """Get the track in the grid that matches the time and mass of a star.

        For "root" matching_method, the properties that are matched is always
        the mass of the secondary star.
        If the secondary has the state `MS` then
        the center hydrogen abundance will also be matched
        otherwise the mass of helium-core will be matched.

        Parameters
        ----------
        star : SingleStar
            The star which properties are required
            to be matched with the single MIST-like grid.

        Returns
        -------
        list of 2 float values
            Contains the associated mass (in solar units) and the time (in years)
            such that the time-series in the grid matches
            the properties of the secondary.

        """
        def get_posydon_attributes(list_for_matching, star):
            
            list_of_attributes = []

            for attr in list_for_matching:
                list_of_attributes.append(getattr(star, attr))
            
            return list_of_attributes
        
        def sq_diff_function(x):

            return self.square_difference(
                                          x, htrack=htrack, mesa_labels=MESA_labels,
                                          posydon_attributes=posydon_attributes,
                                          colscalers=colscalers, scales=scales
                                         )
        
        def get_MESA_labels(list_for_matching):

            MESA_labels = list_for_matching[0]
            
            rs = list_for_matching[1]
            colscalers = list_for_matching[2]
            bnds = []
            for i in range(3, len(list_for_matching)):
                bnds.append(list_for_matching[i])

            if self.verbose:
                print("Matching parameters and their normalizations:\n", MESA_labels, rs)
            
            scales = []
            for MESA_label, colscaler in zip(MESA_labels, colscalers):

                scale_of_attribute = scale(MESA_label, htrack, colscaler)
                scales.append(scale_of_attribute)
            
            return MESA_labels, rs, colscalers, bnds, scales

            
        # setting whether to match to track from H-rich or stripped He star grid
        if htrack:
            self.grid = self.grid_Hrich
        else:
            self.grid = self.grid_strippedHe

        get_root0 = self.get_root0
        get_track_val = self.get_track_val
        matching_method = self.matching_method
        scale = self.scale

        initials = None
        
        tolerance_matching_integration = 1e-2
        tolerance_matching_integration_hard = 1e-1

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
                               htrack, rs=[0.7, 300])
                
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
                               rs=[11, 300]
                              )
                
                sol = root(
                            lambda x: [
                                        get_track_val("he_core_mass", htrack, *x) - star.he_core_mass,
                                        get_track_val("mass", htrack, *x) - star.mass],
                            x0, method="hybr"
                          )
                
            # if opptimizer failed for some reason set solution as NaN
            if not sol.success or sol.x[1] < 0:
                initials = (np.nan, np.nan)
            else:
                initials = sol.x

        # 1st attempt: match using minimize method (ultimately Newton minimization method)
        elif matching_method == "minimize":
            
            # set matching metrics based on star state
            if star.state in LIST_ACCEPTABLE_STATES_FOR_HMS:
                list_for_matching = self.list_for_matching_HMS
            elif star.state in LIST_ACCEPTABLE_STATES_FOR_postMS:
                list_for_matching = self.list_for_matching_postMS
            elif star.state in LIST_ACCEPTABLE_STATES_FOR_HeStar:
                list_for_matching = self.list_for_matching_HeStar
            
            MESA_labels, rs, colscalers, bnds, scales = get_MESA_labels(list_for_matching)
            
            for param_name in MESA_labels:
                if param_name not in self.root_keys:
                    raise AttributeError(f"Expected matching parameter {param_name} not added in the single star grid options.")
                
            posydon_attributes = get_posydon_attributes(MESA_labels, star)

            # get closest time series from grids
            x0 = get_root0(MESA_labels, posydon_attributes, htrack, rs=rs)

            # Match using Newton method w/ x0 time series as guess
            # Minimizing Euclidean distance
            sol = minimize(sq_diff_function, x0, method="TNC", bounds=bnds)

            # save initial matching solution as best solution so far
            best_sol = sol

            ## Alternative matching attempts if default matching fails!
            # 2nd attempt: use a different minimization method
            if (np.abs(best_sol.fun) > tolerance_matching_integration or not best_sol.success):

                if self.verbose:
                    print (f"Initial matching attempt was unsuccessful:"
                           f"\n tolerance {np.abs(sol.fun)} > {tolerance_matching_integration}, sol.success = {sol.success}")

                    print("\nAlternative matching started (1st attempt)")
                    print("(Now trying an alternative minimization method)")
                
                # minimize w/ modified Powell's method
                sol = minimize(sq_diff_function, x0, method="Powell")

                ## if alternative matching has a better solution, make it the new best solution
                if (np.abs(sol.fun) < np.abs(best_sol.fun) and sol.success):
                    best_sol = sol

            # if 2nd fails, 3rd attempt: use alternative matching parameters
            if (np.abs(best_sol.fun) > tolerance_matching_integration or not best_sol.success):      
                   
                if self.verbose:
                    print (f"Alternative matching (1st attempt) was unsuccessful:"
                           f"\n tolerance {np.abs(sol.fun)} > {tolerance_matching_integration}, sol.success = {sol.success}")

                    print("\nAlternative matching started (2nd attempt)")
                    print("(Now trying to match with alternative parameters)")     
                      
                # set alternative matching metrics based on star state
                if star.state in LIST_ACCEPTABLE_STATES_FOR_HMS:
                    list_for_matching = self.list_for_matching_HMS_alternative
                elif star.state in LIST_ACCEPTABLE_STATES_FOR_postMS:
                    list_for_matching = (self.list_for_matching_postMS_alternative)
                elif star.state in LIST_ACCEPTABLE_STATES_FOR_HeStar:
                    list_for_matching = (self.list_for_matching_HeStar_alternative)                
                    
                MESA_labels, rs, colscalers, bnds, scales = get_MESA_labels(list_for_matching)
                posydon_attributes = get_posydon_attributes(MESA_labels, star)

                # Minimize w/ Newton's method using alternative metrics
                x0 = get_root0(MESA_labels, posydon_attributes, htrack, rs=rs)
                sol = minimize(sq_diff_function, x0, method="TNC", bounds=bnds)

                if (np.abs(sol.fun) < np.abs(best_sol.fun) and sol.success):
                    best_sol = sol

            # 4th attempt: match an He-star with an H-rich grid, or vice versa (not applicable for HMS stars)
            if (np.abs(best_sol.fun) > tolerance_matching_integration or not best_sol.success):

                if self.verbose:
                        print (f"Alternative matching (2nd attempt) was unsuccessful:"
                               f"\n tolerance {np.abs(sol.fun)} > {tolerance_matching_integration}, sol.success = {sol.success}")    
                
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

                    MESA_labels, rs, colscalers, bnds, scales = get_MESA_labels(list_for_matching)

                    for param_name in MESA_labels:
                        if param_name not in self.root_keys:
                            raise AttributeError(f"Expected matching parameter {param_name} not added "
                                                  "in the single star grid options.")
                        
                    posydon_attributes = get_posydon_attributes(MESA_labels, star)
                    x0 = get_root0(MESA_labels, posydon_attributes, new_htrack, rs=rs)

                    try:
                        # minimize w/ Euclidean diff. and Newton's method
                        sol = minimize(sq_diff_function, x0, method="TNC", bounds=bnds)

                        if (np.abs(sol.fun) < np.abs(best_sol.fun) and sol.success):
                            best_sol = sol
                            htrack = new_htrack

                        if self.verbose:
                            print (f"Alternative matching (3rd attempt) completed:"
                                   f"\n tolerance {np.abs(sol.fun)} > {tolerance_matching_integration}, "
                                   f"sol.success = {sol.success}")
                    except:
                        raise NumericalError("SciPy numerical differentiation occured outside boundary "
                                             "while matching to single star track")

            # if matching is still not successful, set result to NaN:
            if (np.abs(best_sol.fun) > tolerance_matching_integration_hard or not best_sol.success):
                if self.verbose:
                    print("\nFinal matching result is NOT successful with best tolerance ",
                          np.abs(best_sol.fun), ">", tolerance_matching_integration_hard)
                initials = (np.nan, np.nan)
                
            # or else we found a solution
            else:
                if self.verbose:
                    print("\nFinal matching result is considered successful with best tolerance "
                        f'{np.abs(best_sol.fun):.8f}', "<", tolerance_matching_integration_hard)
                initials = best_sol.x

        if self.verbose:
            # successful match
            if not np.isnan(initials[0]):
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
                    "Matching completed unsuccessfully for star with properties: \n"
                    f'mass = {star.mass:.3f}, ',
                    f'log_R = {star.log_R:.3f}, ',
                    f'center_he4 = {star.center_he4:.4f}, ',
                    f'surface_he4 = {star.surface_he4:.4f}, ',
                    f'surface_h1 = {star.surface_h1:.4f}, ',
                    f'he_core_mass = {star.he_core_mass:.3f}, ',
                    f'center_c12 = {star.center_c12:.4f}'
                )

        return initials[0], initials[1], htrack
    

    def get_star_data(self, binary, star1, star2, htrack, co, copy_prev_m0=None, copy_prev_t0=None):
                """Get and interpolate the properties of stars.

                The data of a compact object can be stored as a copy of its
                companion for convenience except its mass, radius, mdot, and Idot
                are set to be zero.

                Parameters
                ----------
                htrack : bool
                    htrack of star1. False if star 1 is a stripped He star
                co: bool
                    co of star2. True if star 2 is a compact object
                Return
                -------
                interp1d
                    Contains the properties of star1 if co is false,
                    if co is true, star2 is a compact object,
                    return the properties of star2

                """

                KEYS = self.KEYS
                KEYS_POSITIVE = self.KEYS_POSITIVE

                with np.errstate(all="ignore"):
                    # get the initial m0, t0 track
                    if binary.event == 'ZAMS' or binary.event == 'redirect_from_ZAMS':
                        # ZAMS stars in wide (non-mass exchaging binaries) that are
                        # directed to detached step at birth
                        m0, t0 = star1.mass, 0
                    elif co:
                        m0, t0 = copy_prev_m0, copy_prev_t0
                    else:
                        t_before_matching = time.time()
                        # matching to single star grids
                        m0, t0, htrack = self.match_to_single_star(star1, htrack)
                        t_after_matching = time.time()
                        
                        if self.verbose:
                            print(f"Matching duration: {t_after_matching-t_before_matching:.6g} sec\n")
                
                if pd.isna(m0) or pd.isna(t0):
                    return None, None, None
                
                if htrack:
                    self.grid = self.grid_Hrich
                else:
                    self.grid = self.grid_strippedHe
                
                # check if m0 is in the grid
                if m0 < self.grid.grid_mass.min() or m0 > self.grid.grid_mass.max():
                    set_binary_to_failed(binary)
                    raise MatchingError(f"The mass {m0} is out of the single star grid range and "
                                        "cannot be matched to a track.")

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
                    kvalue["mass"] = np.zeros_like(kvalue["mass"]) + star2.mass
                    kvalue["R"] = np.zeros_like(kvalue["log_R"])
                    kvalue["mdot"] = np.zeros_like(kvalue["mdot"])
                    interp1d["mass"] = PchipInterpolator(age, kvalue["mass"])
                    interp1d["R"] = PchipInterpolator(age, kvalue["R"])
                    interp1d["mdot"] = PchipInterpolator(age, kvalue["mdot"])
                    interp1d["Idot"] = PchipInterpolator(age, kvalue["mdot"])

                return interp1d, m0, t0

    def get_star_final_values(self, star, htrack, m0):  

        grid = self.grid_Hrich if htrack else self.grid_strippedHe
        get_final_values = grid.get_final_values

        for key in self.final_keys:
            setattr(star, key, get_final_values('S1_%s' % (key), m0))

    def get_star_profile(self, star, htrack, m0):

        grid = self.grid_Hrich if htrack else self.grid_strippedHe
        get_profile = grid.get_profile
        profile_new = np.array(get_profile('mass', m0)[1])

        for i in self.profile_keys:
            profile_new[i] = get_profile(i, m0)[0]
        profile_new['omega'] = star.surf_avg_omega

        star.profile = profile_new