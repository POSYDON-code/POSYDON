import numpy as np
from scipy.optimize import root
from scipy.optimize import minimize

from posydon.utils.posydonwarning import Pwarn

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

    def __init__(self):

        # for the matching
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

        # Initialize the matching lists:
        m_min_H = np.min(self.grid_Hrich.grid_mass)
        m_max_H = np.max(self.grid_Hrich.grid_mass)
        m_min_He = np.min(self.grid_strippedHe.grid_mass)
        m_max_He = np.max(self.grid_strippedHe.grid_mass)
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

        # lists of alternative matching

        # e.g., stars after mass transfer could swell up so that log_R
        # is not appropriate for matching

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
            x : list of floats, of same length as "keys"
                Contains the latest values (from a previous POSYDON step) of the
                quantities of "keys" in the POSYDON SingleStar object.
            rs : list of floats, same length as "keys"
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
        grid = self.grid_Hrich if htrack else self.grid_strippedHe
        self.initial_mass = grid.grid_mass
        n = 0

        for mass in grid.grid_mass:
            n = max(n, len(grid.get("age", mass)))

        self.rootm = np.inf * np.ones((len(grid.grid_mass),
                                    n, len(self.root_keys)))
        
        for i, mass in enumerate(grid.grid_mass):
            for j, key in enumerate(self.root_keys):
                track = grid.get(key, mass)
                self.rootm[i, : len(track), j] = track
        if rs is None:
            rs = np.ones_like(keys)
        else:
            rs = np.asanyarray(rs)

        x = np.asanyarray(x)
        idx = np.argmax(np.asanyarray(keys)[:, None] == self.root_keys, axis=1)
        X = self.rootm[:, :, idx]
        d = np.linalg.norm((X - x[None, None, :]) / rs[None, None, :], axis=-1)
        idx = np.unravel_index(d.argmin(), X.shape[:-1])
        t = self.rootm[idx][np.argmax("age" == self.root_keys)]
        m0 = grid.grid_mass[idx[0]]

        return m0, t

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
                colscalers=colscalers, scales=scales)
        
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

        if matching_method == "root":

            if star.state in LIST_ACCEPTABLE_STATES_FOR_HMS:
                x0 = get_root0(["center_h1", "mass"],
                               [star.center_h1, star.mass],
                               htrack, rs=[0.7, 300])
                sol = root(
                    lambda x: [
                        get_track_val("center_h1", htrack, *x) - star.center_h1,
                        get_track_val("mass", htrack, *x) - star.mass],
                    x0, method="hybr")
            else:
                x0 = get_root0(
                    ["he_core_mass", "mass"],
                    [star.he_core_mass, star.mass],
                    htrack,
                    rs=[11, 300])
                
                sol = root(
                    lambda x: [
                        get_track_val("he_core_mass", htrack, *x) - star.he_core_mass,
                        get_track_val("mass", htrack, *x) - star.mass],
                    x0, method="hybr")
                
            if not sol.success or sol.x[1] < 0:
                initials = (np.nan, np.nan)
            else:
                initials = sol.x

        elif matching_method == "minimize":
            
            if star.state in LIST_ACCEPTABLE_STATES_FOR_HMS:
                list_for_matching = self.list_for_matching_HMS
            elif star.state in LIST_ACCEPTABLE_STATES_FOR_postMS:
                list_for_matching = self.list_for_matching_postMS
            elif star.state in LIST_ACCEPTABLE_STATES_FOR_HeStar:
                list_for_matching = self.list_for_matching_HeStar
            
            MESA_labels, rs, colscalers, bnds, scales = get_MESA_labels(list_for_matching)
            
            for i in MESA_labels:
                if i not in self.root_keys:
                    raise AttributeError(f"Expected matching parameter {i} not added in the single star grid options.")
                
            posydon_attributes = get_posydon_attributes(MESA_labels, star)

            x0 = get_root0(MESA_labels, posydon_attributes, htrack, rs=rs)
            sol = minimize(sq_diff_function, x0, method="TNC", bounds=bnds)


            ## Alternative matching attempts if default matching fails!
            # 1st attempt: use a different minimization method
            if (np.abs(sol.fun) > tolerance_matching_integration or not sol.success):
                
                if self.verbose:
                    print("\nAlternative matching started (1st attempt) "
                          "because previous attempt was unsuccessful:\n",
                          f"tolerance {np.abs(sol.fun)} > {tolerance_matching_integration}",
                          f"or sol.success = {sol.success}")
                    print("(Now trying an alternative minimization method)")
                    
                sol = minimize(sq_diff_function, x0, method="Powell")

            # 2nd attempt: use alternative matching parameters
            if (np.abs(sol.fun) > tolerance_matching_integration or not sol.success):    
                   
                if self.verbose:
                    print("\nAlternative matching started (2nd attempt) "
                          "because previous attempt was unsuccessful:\n",
                          f"tolerance {np.abs(sol.fun)} > {tolerance_matching_integration}",
                          f"or sol.success = {sol.success}")
                    print("(Now trying to match with alternative parameters)")     
                      
                if star.state in LIST_ACCEPTABLE_STATES_FOR_HMS:
                    list_for_matching = self.list_for_matching_HMS_alternative
                elif star.state in LIST_ACCEPTABLE_STATES_FOR_postMS:
                    list_for_matching = (self.list_for_matching_postMS_alternative)
                elif star.state in LIST_ACCEPTABLE_STATES_FOR_HeStar:
                    list_for_matching = (self.list_for_matching_HeStar_alternative)                
                    
                MESA_labels, rs, colscalers, bnds, scales = get_MESA_labels(list_for_matching)
                posydon_attributes = get_posydon_attributes(MESA_labels, star)

                x0 = get_root0(MESA_labels, posydon_attributes, htrack, rs=rs)
                sol = minimize(sq_diff_function, x0, method="TNC", bounds=bnds)

            # 3rd attempt: match an He-star with an H-rich grid, or vice versa (not applicable for HMS stars)
            if (np.abs(sol.fun) > tolerance_matching_integration or not sol.success):
                
                if (star.state in LIST_ACCEPTABLE_STATES_FOR_HeStar
                    or star.state in LIST_ACCEPTABLE_STATES_FOR_postMS):

                    Pwarn("Attempting to match an He-star with an H-rich grid or post-MS star with a"
                          " stripped-He grid", "EvolutionWarning")

                    if self.verbose:
                        print("\nAlternative matching started (3rd attempt) "
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

                    for i in MESA_labels:
                        if i not in self.root_keys:
                            raise AttributeError(f"Expected matching parameter {i} not added "
                                                 "in the single star grid options.")
                        
                    posydon_attributes = get_posydon_attributes(MESA_labels, star)
                    x0 = get_root0(MESA_labels, posydon_attributes, htrack, rs=rs)

                    try:
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
                
            elif np.abs(sol.fun) < tolerance_matching_integration_hard:
                if self.verbose:
                    print("\nMatching result is considered successful, with tolerance "
                        f'{np.abs(sol.fun):.8f}', "<", tolerance_matching_integration_hard)
                initials = sol.x

        if self.verbose:
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