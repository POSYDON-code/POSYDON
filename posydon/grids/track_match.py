import numpy as np
from scipy.optimize import root
from scipy.optimize import minimize

from posydon.utils.posydonwarning import Pwarn
from posydon.interpolation.data_scaling import DataScaler
from posydon.utils.posydonerror import (NumericalError)

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

class dummy:


    def __init__(
            self,
            grid_Hrich,
            grid_strippedHe,
            metallicity=None,
            matching_method="minimize",
            initial_mass=None,
            rootm=None,
            verbose=False,
            list_for_matching_HMS=None,
            list_for_matching_postMS=None,
            list_for_matching_HeStar=None
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

        #self.metallicity = convert_metallicity_to_string(metallicity)
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
        #self.KEYS = DEFAULT_TRANSLATED_KEYS
        #self.KEYS_POSITIVE = (
        #    'mass_conv_reg_fortides',
        #    'thickness_conv_reg_fortides',
        #    'radius_conv_reg_fortides'
        #)

        # keys for the final value interpolation
        #self.final_keys = (
        #    'avg_c_in_c_core_at_He_depletion',
        #    'co_core_mass_at_He_depletion',
        #    'm_core_CE_1cent',
        #    'm_core_CE_10cent',
        #    'm_core_CE_30cent',
        #    'm_core_CE_pure_He_star_10cent',
        #    'r_core_CE_1cent',
        #    'r_core_CE_10cent',
        #    'r_core_CE_30cent',
        #    'r_core_CE_pure_He_star_10cent'
        #)

        # keys for the star profile interpolation
        #self.profile_keys = DEFAULT_PROFILE_KEYS

        # should grids just get passed to this?
        #? if grid_name_Hrich is None:
        #?     grid_name_Hrich = os.path.join('single_HMS', self.metallicity+'_Zsun.h5')
        #? self.grid_Hrich = GRIDInterpolator(os.path.join(path, grid_name_Hrich))
        self.grid_Hrich = grid_Hrich

        #? if grid_name_strippedHe is None:
        #?    grid_name_strippedHe = os.path.join('single_HeMS', self.metallicity+'_Zsun.h5')
        #? self.grid_strippedHe = GRIDInterpolator(os.path.join(path, grid_name_strippedHe))
        self.grid_strippedHe = grid_strippedHe

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
            Contains the associated (in solar units) and the time (in years)
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


            ## Alternative matching attempts if default matching fails!
            # 2nd attempt: use a different minimization method
            if (np.abs(sol.fun) > tolerance_matching_integration or not sol.success):
                
                if self.verbose:
                    print("\nAlternative matching started (1st retry) "
                          "because previous attempt was unsuccessful:\n",
                          f"tolerance {np.abs(sol.fun)} > {tolerance_matching_integration}",
                          f"or sol.success = {sol.success}")
                    print("(Now trying an alternative minimization method)")
                
                # minimize w/ modified Powell's method
                sol = minimize(sq_diff_function, x0, method="Powell")

            # if 2nd fails, 3rd attempt: use alternative matching parameters
            if (np.abs(sol.fun) > tolerance_matching_integration or not sol.success):    
                   
                if self.verbose:
                    print("\nAlternative matching started (2nd retry) "
                          "because previous attempt was unsuccessful:\n",
                          f"tolerance {np.abs(sol.fun)} > {tolerance_matching_integration}",
                          f"or sol.success = {sol.success}")
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

            # 4th attempt: match an He-star with an H-rich grid, or vice versa (not applicable for HMS stars)
            if (np.abs(sol.fun) > tolerance_matching_integration or not sol.success):
                
                # if post-MS or stripped He star
                if (star.state in LIST_ACCEPTABLE_STATES_FOR_HeStar
                    or star.state in LIST_ACCEPTABLE_STATES_FOR_postMS):

                    Pwarn("Attempting to match an He-star with an H-rich grid or post-MS star with a"
                          " stripped-He grid", "EvolutionWarning")

                    if self.verbose:
                        print("\nAlternative matching started (3rd retry) "
                              "because previous attempt was unsuccessful:\n", 
                              f"tolerance {np.abs(sol.fun)} > {tolerance_matching_integration}",
                              f"or sol.success = {sol.success}")
                        print("(Now trying to match star to a different grid)")
                        
                    if star.state in LIST_ACCEPTABLE_STATES_FOR_HeStar:
                        
                        htrack = True
                        list_for_matching = self.list_for_matching_HeStar

                    elif star.state in LIST_ACCEPTABLE_STATES_FOR_postMS:
                        
                        htrack = False
                        list_for_matching = self.list_for_matching_postMS

                    MESA_labels, rs, colscalers, bnds, scales = get_MESA_labels(list_for_matching)

                    for param_name in MESA_labels:
                        if param_name not in self.root_keys:
                            raise AttributeError(f"Expected matching parameter {param_name} not added "
                                                  "in the single star grid options.")
                        
                    posydon_attributes = get_posydon_attributes(MESA_labels, star)
                    x0 = get_root0(MESA_labels, posydon_attributes, htrack, rs=rs)

                    try:
                        # minimize w/ Euclidean diff. and Newton's method
                        sol = minimize(sq_diff_function, x0, method="TNC", bounds=bnds)
                    except:
                        raise NumericalError("SciPy numerical differentiation occured outside boundary "
                                             "while matching to single star track")

            # if matching is still not successful, set result to NaN:
            if (np.abs(sol.fun) > tolerance_matching_integration_hard or not sol.success):
                if self.verbose:
                    print("\nMatching result is NOT successful, with tolerance ",
                          np.abs(sol.fun), ">", tolerance_matching_integration_hard)
                initials = (np.nan, np.nan)
                
            # or else we found a solution
            elif np.abs(sol.fun) < tolerance_matching_integration_hard:
                if self.verbose:
                    print("\nMatching result is considered successful, with tolerance "
                        f'{np.abs(sol.fun):.8f}', "<", tolerance_matching_integration_hard)
                initials = sol.x

        if self.verbose:
            # successful match
            if not np.isnan(initials[0]):
                print(
                    "Matching completed for", star.state, "star!\n"
                    f"Matched to track with intial mass m0 = {initials[0]:.3f} [Msun]"
                    f" at time t0 = {initials[1]/1e6:.3f} [Myrs] \n",
                    "and m(t0), log10(R(t0), center_he(t0), surface_he4(t0), "
                    "surface_h1(t0), he_core_mass(t0), center_c12(t0) = \n",
                    f'{self.get_track_val("mass", htrack, *sol.x):.3f}',
                    f'{self.get_track_val("log_R", htrack, *sol.x):.3f}',
                    f'{self.get_track_val("center_he4", htrack, *sol.x):.4f}',
                    f'{self.get_track_val("surface_he4", htrack, *sol.x):.4f}',
                    f'{self.get_track_val("surface_h1", htrack, *sol.x):.4f}',
                    f'{self.get_track_val("he_core_mass", htrack, *sol.x):.3f}',
                    f'{self.get_track_val("center_c12", htrack, *sol.x):.4f}\n',
                    "The same values of the original star at the end of the previous "
                    "step were: \n",
                    f'{star.mass:.3f}',
                    f'{star.log_R:.3f}',
                    f'{star.center_he4:.4f}',
                    f'{star.surface_he4:.4f}',
                    f'{star.surface_h1:.4f}',
                    f'{star.he_core_mass:.3f}',
                    f'{star.center_c12:.4f}'
                )
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