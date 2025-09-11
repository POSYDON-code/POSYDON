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
    "Seth Gossage <seth.gossage@northwestern.edu>",
    "Max Briel <max.briel@gmail.com>",
]

import os
import time
import numpy as np
import pandas as pd
import types
import scipy
from scipy.optimize import root
from scipy.optimize import minimize
from scipy.interpolate import PchipInterpolator

import posydon.utils.constants as const
from posydon.config import PATH_TO_POSYDON_DATA
from posydon.utils.posydonwarning import Pwarn
from posydon.interpolation.data_scaling import DataScaler
from posydon.interpolation.interpolation import GRIDInterpolator
from posydon.utils.interpolators import PchipInterpolator2
from posydon.utils.posydonerror import (NumericalError, MatchingError,
                                        POSYDONError)
from posydon.utils.common_functions import (convert_metallicity_to_string,
                                            set_binary_to_failed)

from posydon.binary_evol.DT.key_library import (DEFAULT_TRANSLATED_KEYS,
                                                KEYS_POSITIVE,
                                                DEFAULT_PROFILE_KEYS)

from posydon.binary_evol.flow_chart import (STAR_STATES_CO,
                                            STAR_STATES_H_RICH,
                                            STAR_STATES_FOR_HMS_MATCHING,
                                            STAR_STATES_FOR_postMS_MATCHING,
                                            STAR_STATES_FOR_Hestar_MATCHING)

MATCHING_WITH_RELATIVE_DIFFERENCE = ["center_he4"]


val_names = [" ", "mass", "log_R", "center_h1", "surface_h1",
                "he_core_mass", "center_he4", "surface_he4", 
                "center_c12"]
str_fmts = ["{:>14}", "{:>9}","{:>9}", 
            "{:>9}",  "{:>10}",  "{:>12}",  
            "{:>10}",  "{:>11}",  "{:>10}"]
row_str = " ".join(str_fmts)
DIVIDER_STR = "_"*len(row_str.format(*[""]*len(str_fmts)))
# MAJOR.MINOR version of imported scipy package
SCIPY_VER = float('.'.join(scipy.__version__.split('.')[:2]))

class TrackMatcher:
    """
        This class contains the functionality to match binary star components 
    to single star models for detached evolution. Typically, this is so that 
    the binary star can continue evolving in a detached (non-interacting) 
    state.

        Several matching methods may be used, and metrics (physical 
    quantities) that are used to determine the quality of the match may be 
    customized. By default, the metrics are as follows:

      Initial stellar matching metrics
      ========================================
        
      MS:          mass,                      center X, radius, He core mass
      post-MS:     mass,                      center Y, radius, He core mass
      stripped He: He core mass (mass proxy), center Y, radius

        If an initial match can not be found, alternative sets of metrics are 
    used. For example, stars after mass transfer could swell up so that log_R 
    is not appropriate for matching. Lists for HMS and postMS stars drop 
    radius as a matching metric in these alternatives:
    
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

    Attributes
    ----------
    KEYS : list[str]
        Contains valid keywords which are used to extract quantities from 
        the grids.

    KEYS_POSITIVE : list[str]
        Keys in this list are forced to be positive or else 0 by the 
        posydon.utils.PchipInterpolator2 class following interpolation 
        of the associated quantity.

    path : str
        Path to the directory that contains POSYDON data HDF5 files. Defaults 
        to the PATH_TO_POSYDON_DATA environment variable.

    metallicity : str
        The metallicity of the grid. This should be one of the eight 
        supported metallicities, stored as a string (e.g.,
        1e+00 as "1e+00_Zsun").

    matching_method : str
        Method to find the best match between a star from a previous step and a 
        point in a single star evolution track. Options: 
        
            "root": Tries to find a root of two matching quantities. It is 
                    possible to not find one, causing the evolution to fail.

            "minimize": Minimizes the sum of squares of differences of 
                        various quantities between the previous evolution step 
                        and a stellar evolution track. 

    root_keys : numpy.ndarray
        An array of keys corresponding to possible matching metrics. These 
        keys should exist as MESA history data column names. In practice, 
        we match using only a subset of these.
           
    rootm : numpy.ndarray
        A 3D matrix to hold roots with dimensions 
            
            DIM = [N(initial_masses), 
                   N(max_evolution_track_length), 
                   N(matching_metrics)]

        Structured to hold the matching metrics along the entire evolution 
        track of each stellar evolution track of a given initial mass in 
        a single star grid. Assigned after loading a grid and before 
        storing matching metrics.

    grid_name_Hrich : str
        Name of the single star H-rich grid h5 file, 
        including its parent directory. This is set to 
        (w/ Z = 1e+00, for example):

            grid_name_Hrich = 'single_HMS/1e+00_Zsun.h5'  

        by default if not specified.

    grid_name_strippedHe : str
        Name of the single star He-rich grid h5 file. This is 
        set to (w/ Z = 1e+00, for example):

            grid_name_strippedHe = 'single_HeMS/1e+00_Zsun.h5'
        
        by default if not specified.
        
    grid_Hrich : GRIDInterpolator object
        Object to interpolate between the time-series (i.e., along the 
        evolutionary track) in the H-rich single star h5 grid.

    grid_strippedHe : GRIDInterpolator object
        Object to interpolate between the time-series (i.e., along the 
        evolutionary track) in the He-rich single star h5 grid.

    initial_mass : list[float]
            Contains the initial masses of the stars in the single star 
            grid in which we are searching for a match. Assigned after 
            loading a grid to match to.

    list_for_matching_HMS : list
        A list of mixed type that specifies properties of the matching 
        process for HMS stars. This list has the following structure: 
        
            list_for_matching = [[matching attr. names], [rescale_factors],
                                 [scaling method], [mass_bnds], [age_bnds]]

    list_for_matching_postMS : list
        A list of mixed type that specifies properties of the matching 
        process for postMS stars. This list has the following structure: 
        
            list_for_matching = [[matching attr. names], [rescale_factors],
                                 [scaling method], [mass_bnds], [age_bnds]]

    list_for_matching_HeStar : list
        A list of mixed type that specifies properties of the matching 
        process for He stars. This list has the following structure: 
        
            list_for_matching = [[matching attr. names], [rescale_factors],
                                 [scaling method], [mass_bnds], [age_bnds]]

    stored_scalers : dict
        Mapping a combination of (key, htrack, scaling_method) to a 
        pre-trained DataScaler instance.

    final_keys : tuple
        Contains keys for final value interpolation.

    profile_keys : tuple
        Contains keys for profile interpolation.

    matching_tolerance : float
        When using the "minimize" matching method, a computed square 
        Euclidean distance between the pre-match and post-match values 
        less than this must be achieved for a successful match.

    matching_tolerance_hard : float
        When using the "minimize" matching method, a computed square 
        Euclidean distance between the pre-match and post-match values 
        less than at most this must be achieved for a successful match. 
        This tolerance is checked after all else fails, as a last attempt 
        to find a solution.

    record_matching : bool
        True if we want to append matched quantities to the binary history.

    verbose : bool
        True if we want to print stuff.

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
            grid_name_Hrich,
            grid_name_strippedHe,
            path=PATH_TO_POSYDON_DATA,
            metallicity=None,
            matching_method="minimize",
            matching_tolerance=1e-2,
            matching_tolerance_hard=1e-1,
            list_for_matching_HMS=None,
            list_for_matching_postMS=None,
            list_for_matching_HeStar=None,
            record_matching=False,
            verbose=False
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

        # =====================================================================

        self.metallicity = convert_metallicity_to_string(metallicity)
        self.matching_method = matching_method
        self.matching_tolerance = matching_tolerance # DEFAULT: 1e-2
        self.matching_tolerance_hard = matching_tolerance_hard # DEFAULT: 1e-1

        self.initial_mass = None
        self.rootm = None
        self.verbose = verbose

        self.list_for_matching_HMS = list_for_matching_HMS
        self.list_for_matching_postMS = list_for_matching_postMS
        self.list_for_matching_HeStar = list_for_matching_HeStar

        # mapping a combination of (key, htrack, method) to a pre-trained
        # DataScaler instance, created the first time it is requested
        self.stored_scalers = {}

        # these are the KEYS read from POSYDON h5 grid files (after translating
        # them to the appropriate columns)
        self.KEYS = DEFAULT_TRANSLATED_KEYS #KEYS #DEFAULT_TRANSLATED_KEYS
        self.KEYS_POSITIVE = KEYS_POSITIVE

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
            grid_name_Hrich = os.path.join('single_HMS', 
                                           self.metallicity+'_Zsun.h5')
        grid_path_Hrich = os.path.join(path, grid_name_Hrich)
        self.grid_Hrich = GRIDInterpolator(grid_path_Hrich)

        if grid_name_strippedHe is None:
            grid_name_strippedHe = os.path.join('single_HeMS', 
                                                self.metallicity+'_Zsun.h5')
        grid_path_strippedHe = os.path.join(path, grid_name_strippedHe)
        self.grid_strippedHe = GRIDInterpolator(grid_path_strippedHe)

        # =====================================================================

        # Initialize the matching lists:

        # min/max ranges of initial masses for each grid
        m_min_H = np.min(self.grid_Hrich.grid_mass)
        m_max_H = np.max(self.grid_Hrich.grid_mass)
        m_min_He = np.min(self.grid_strippedHe.grid_mass)
        m_max_He = np.max(self.grid_strippedHe.grid_mass)

        # min/max ranges of ages for each grid
        t_min_H = 0.0
        t_max_H = None
        t_min_He = 0.0
        t_max_He = None

        # Stellar parameter matching metrics
        # ========================================
        # At first, we try to match based on...

        if self.list_for_matching_HMS is None:
            self.list_for_matching_HMS = [
                ["mass", "center_h1", "log_R", "he_core_mass"],
                [20.0, 1.0, 2.0, 10.0],
                ["log_min_max", "min_max", "min_max", "min_max"],
                [m_min_H, m_max_H], [t_min_H, t_max_H]
            ]
        if self.list_for_matching_postMS is None:
            self.list_for_matching_postMS = [
                ["mass", "center_he4", "log_R", "he_core_mass"],
                [20.0, 1.0, 2.0, 10.0],
                ["log_min_max", "min_max", "min_max", "min_max"],
                [m_min_H, m_max_H], [t_min_H, t_max_H]
            ]
        if self.list_for_matching_HeStar is None:
            self.list_for_matching_HeStar = [
                ["he_core_mass", "center_he4", "log_R"],
                [10.0, 1.0, 2.0],
                ["min_max", "min_max", "min_max"],
                [m_min_He, m_max_He], [t_min_He, t_max_He]
            ]

        # Stellar parameter lists for alternative matching metrics
        # These are used in the event that an initial match can not be found
        self.list_for_matching_HMS_alternative = [
            ["mass", "center_h1", "he_core_mass"],
            [20.0, 1.0, 10.0],
            ["log_min_max", "min_max", "min_max"],
            [m_min_H, m_max_H], [t_min_H, t_max_H]
        ]
        self.list_for_matching_postMS_alternative = [
            ["mass", "center_h1", "he_core_mass"],
            [20.0, 1.0, 10.0],
            ["log_min_max", "min_max", "min_max"],
            [m_min_H, m_max_H], [t_min_H, t_max_H]
        ]
        self.list_for_matching_HeStar_alternative = [
            ["he_core_mass", "center_he4", "log_R"],
            [10.0, 1.0, 2.0],
            ["min_max", "min_max", "min_max"],
            [m_min_He, m_max_He], [t_min_He, t_max_He]
        ]

        self.record_matching = record_matching
    
    def get_root0(self, attr_names, attr_vals, htrack, rescale_facs=None):
        """
            Get the stellar evolution track in the single star grid with values 
        closest to the requested ones. This calculates the difference in 
        stellar properties between a given track and those in the single star 
        grids. It then returns the mass of and age along the track where the 
        minimum difference occurs.

        Parameters
        ----------
        attr_names : list[str]
            Contains the keys of the requested specific quantities that will 
            be matched in the single star track.

        attr_vals : list[float]
            Contains the latest values (from a previous POSYDON step) of the 
            quantities of "keys" in the POSYDON SingleStar object. This should 
            be the same length as "keys".

        htrack : bool
            Set True to search the single star H-rich grids, or False to 
            search the He-rich grids.
            
        rescale_facs : list[float]
            Contains normalization factors to be divided for rescaling 
            attribute values. This should be the same length as attr_names. 

        Returns
        -------
        m0 : float
            Initial mass (in solar units) of the minimum diff model for initial 
            guess.
    
        t0 : float
            Age (in years) of the minimum diff model for initial guess.

        """

        # set which grid to search based on htrack condition
        grid = self.grid_Hrich if htrack else self.grid_strippedHe

        # initial masses within grid (defined but never used? used in scale())
        self.initial_mass = grid.grid_mass

        # search across all initial masses and get max track length
        max_track_length = 0
        for mass in grid.grid_mass:
            track_length = len(grid.get("age", mass))
            max_track_length = max(max_track_length, track_length)

        # intialize root matrix 
        # (DIM = [N(Mi), N(max_track_length), N(root_keys)])
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
        idx = np.argmax(np.asanyarray(attr_names)[:, None] == self.root_keys,
                        axis=1)

        # Slice out just the matching metric data for all stellar tracks
        # grid_attr_vals now has shape 
        # (N(Mi), N(max_track_len), N(matching_metrics))
        grid_attr_vals = self.rootm[:, :, idx]

        # For all stellar tracks in grid:
        # Take difference btwn. grid track and given star values...
        grid_diff = grid_attr_vals - star_attr_vals[None, None, :]
        # rescale... 
        grid_diff /= rescale_facs[None, None, :]
        # and take (Frobenius) norm.
        grid_diff = np.linalg.norm(grid_diff, axis=-1)

        # indices where minimum difference occurs (i.e., of 
        # closest matching track). This contains the indices at 
        # position of (min diff mass, min diff age)
        min_diff_inds = np.unravel_index(grid_diff.argmin(), 
                                         grid_attr_vals.shape[:-1])

        mass_i = min_diff_inds[0]
        age_i = min_diff_inds[1]

        # time and initial mass corresp. to track w/ minimum difference
        m0 = grid.grid_mass[mass_i]
        t0 = self.rootm[mass_i][age_i][np.argmax("age" == self.root_keys)]

        return m0, t0
    
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

    def match_to_single_star(self, star):
        """
            Get the track in the grid that matches the time and mass of a star,
        that has typically undergone prior binary star evolution. A match is 
        made according to a given algorithm that minimizes the difference 
        amongst several physical properties, e.g., mass, central hydrogen 
        abundance, radius, and core helium mass, depending on the type of star 
        being matched. However, these properties may also be customized by the 
        user.

        Parameters
        ----------
        star : SingleStar object
            A single star object that contains the star's properties.

        Returns
        -------
        match_m0 : float
            Initial mass (in solar units) of the matched model.
        
        match_t0 : float
            Age (in years) of the matched model.

        Warns
        -----
        EvolutionWarning
            If attempting to match an He-star with an H-rich grid or post-MS 
            star with a stripped-He grid. This can happen if an initial 
            matching attempt fails; alternative grids are checked for a match 
            in such cases.

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

        # setting whether to match to track from H-rich 
        # or stripped He star grid
        if star.htrack:
            self.grid = self.grid_Hrich
        else:
            self.grid = self.grid_strippedHe

        get_root0 = self.get_root0
        get_track_val = self.get_track_val

        match_vals = [None, None]

        if self.verbose:
            print(f"\nMatching process started in detached step for "
                  f"{star.state} star with matching "
                  f"method = {self.matching_method}")
            if self.matching_method == "minimize":
                  print(f"matching_tolerance = {self.matching_tolerance}\n"
                        "matching_tolerance_hard = "
                        f"{self.matching_tolerance_hard}")

        # matching via root method (ultimately modified Powell's method)
        if self.matching_method == "root":
            match_vals, best_sol = self.match_through_root(star, 
                                                            get_root0, 
                                                            get_track_val)

        # matching using minimize method (Newton or Powell method)
        elif self.matching_method == "minimize":
            match_vals, best_sol = self.match_through_minimize(star, 
                                                                get_root0, 
                                                                get_track_val)

        if self.verbose:
            # successful match
            if all(pd.notna(match_vals)):

                getv = lambda n: get_track_val(n, star.htrack, *best_sol.x)
                sol_vals = [getv(vn) for vn in val_names[1::]]
                init_vals = [getattr(star, vn) for vn in val_names[1::]]

                init_val_str = ["initial values"]+\
                               [f"{v:.3e}" for v in init_vals]

                sol_val_str = ["matched values"]+\
                              [f"{v:.3e}" for v in sol_vals]

                pcntd_fn = lambda s, i: np.nan_to_num(100.0*(s - i)/i, 0)
                zip_vals = zip(sol_vals, init_vals)
                pcntd_vals = [pcntd_fn(sval, ival) for sval, ival in zip_vals]
                percent_diff_str = ["% difference"]+\
                                   [f"{pcntd:.1e}%" for pcntd in pcntd_vals]

                output_table = [val_names, init_val_str, sol_val_str, 
                                percent_diff_str]

                print("\nMatching completed for", star.state, "star!\n")
                for row in output_table:
                    print(row_str.format(*row))

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

            # done with matching attempts
            print(DIVIDER_STR)

        match_m0 = match_vals[0]
        match_t0 = match_vals[1]

        return match_m0, match_t0

    def match_through_root(self, star, get_root0, get_track_val):
        """Match the star through a root method.
        
        Parameters
        ----------
        star : SingleStar object
            This is a SingleStar object, typically representing 
            a member of a binary star system that we are trying 
            to match to a single star track.

        get_root0 : function
            Function that is used to get an initial guess for a 
            closely matching stellar track.

        get_track_val : function
            Function that returns a single value of a stellar property
            from the interpolated time-series along a requested stellar
            track of mass `m0` at an age of `t`.

        Returns
        -------
        match_vals : array
            Mass and age of the star found through matching. These serve 
            as starting points for subsequent (post-match) evolution.

        best_sol : OptimizeResult object
            The OptimizeResult object that contains attributes of the 
            closest matching stellar track. Produced by SciPy's `root()`.
                        
        """

        htrack = star.htrack

        # if the star can be considered an HMS star
        if star.state in STAR_STATES_FOR_HMS_MATCHING:

            # get initial guess for matching via search for closest track
            x0 = get_root0(["center_h1", "mass"],
                            [star.center_h1, star.mass],
                            htrack, rescale_facs=[0.7, 300])

            # Using modified Powell method to find match starting 
            # from time series from above (in x0) using difference 
            # of center X and stellar mass to find match
            sol = root(
                        lambda x: [
                                    get_track_val("center_h1",
                                                    htrack,
                                                    *x) - star.center_h1,
                                    get_track_val("mass",
                                                    htrack,
                                                    *x) - star.mass
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
                                    get_track_val("he_core_mass", 
                                                    htrack, 
                                                    *x) - star.he_core_mass,
                                    get_track_val("mass", 
                                                    htrack, 
                                                    *x) - star.mass],
                        x0, method="hybr"
                        )

        # if optimizer failed for some reason set solution as NaN
        if not sol.success or sol.x[1] < 0:
            match_vals = np.array([np.nan, np.nan])
        else:
            match_vals = sol.x

        return match_vals, sol

    def match_through_minimize(self, star, get_root0, get_track_val):
        """
            Match the star through a minimization method. An initial 
        match is found using the `get_root0()` function. Then, those
        initial values are used in SciPy's minimize method to find 
        the closest matching point in evolution, based on several 
        physical matching metrics.
        
        Parameters
        ----------
        star : SingleStar object
            This is a SingleStar object, typically representing 
            a member of a binary star system that we are trying 
            to match to a single star track.

        get_root0 : function
            Function that is used to get an initial guess for a 
            closely matching stellar track.

        get_track_val : function
            Function that returns a single value of a stellar property
            from the interpolated time-series along a requested stellar
            track of mass `m0` at an age of `t`.

        Returns
        -------
        match_vals : array
            Mass and age of the star found through matching. These serve 
            as starting points for subsequent (post-match) evolution.

        best_sol : OptimizeResult object
            The OptimizeResult object that contains attributes of the 
            closest matching stellar track. Produced by SciPy's `minimize()`.
                        
        """
        # START SUBFUNCTION DEFIITIONS
        def get_attr_values(attr_names):

            """
                Given a set of attribute names, gets attribute values 
            from a given SingleStar object. These values correspond 
            to the last evolution step of the SingleStar object.

            Parameters
            ----------
            attr_names : list[str]
                This list contains strings that are the names of 
                data columns (attributes) to be used as matching 
                metrics.

            Returns
            -------
            attr_values : list[float]
                This is a list of the values associated with the 
                provided attribute names at the last evolution 
                step experienced by `star`.

            """

            attr_values = []

            for attr_name in attr_names:
                attr_values.append(getattr(star, attr_name))

            return attr_values

        def square_difference(x, new_htrack, attr_names, attr_vals, 
                              attr_scalers):
            """
                Compute the square difference between values along a single 
            star evolution track and a given set of values. To find a 'good 
            match', between the given values and the single star track, we 
            ultimately seek to minimize these differences. This is the function 
            used by the minimization algoritms when using matching_method = 
            `minimize`. The minimize function passes values of mass and age 
            that are used to search either the H-rich or He-rich single star 
            grid tracks for a match.

            Parameters
            ----------
            x : list[float]
                Contains a mass and time. These are used to get the square 
                difference between a track of that mass at that time and 
                the given star values from prior evolution. 
            
            new_htrack : bool
                Set True to search the single star H-rich grids, or False to 
                search the He-rich grids. This is determined primarily 
                through the current `match_type`.

            attr_names : list[str]
                Names of SingleStar object attributes that this will calculate 
                the square differences of.

            attr_vals : list[float]
                Associated values of provided SingleStar object attribute names. 
                These values should correspond to the prior point in evolution 
                to which we are finding the stellar track match.

            scalers : list[DataScaler object]
                DataScaler objects that have been trained previously. These are 
                used to scale given star attributes to the range (0, 1).

            Returns
            -------
            result : float
                The square difference between a single star track with mass and 
                age taken from `x` and the given attributes of a SingleStar 
                object.

            """

            result = 0.0
            zip_attr_props = zip(attr_names, attr_vals, attr_scalers)
            for attr_name, star_val, attr_scaler in zip_attr_props:
                
                # get interpolated attribute value from grid at 
                # proposed x (m0, t) values
                # TODO: can this be sped up? w/o this call, 
                #       execution time lowers by a factor of 100
                grid_track_val = get_track_val(attr_name, new_htrack, *x)

                # scale values
                scaled_grid_track_val = attr_scaler.transform(grid_track_val)
                scaled_star_val = attr_scaler.transform(star_val)

                if attr_name in MATCHING_WITH_RELATIVE_DIFFERENCE:
                    result += ((scaled_grid_track_val - scaled_star_val)
                            / scaled_star_val) ** 2
                else:
                    result += (scaled_grid_track_val - scaled_star_val) ** 2

            return result

        def get_attr_props(new_htrack, list_for_matching):

            """
               This unpacks a list_for_matching which has the following 
            structure:

                list_for_matching = [[matching attr. names], [rescale_factors],
                                    [scaling method], [mass_bnds], [age_bnds]]
            
               This also trains DataScaler objects for each provided matching 
            attribute, such that the values are scaled to the range (0, 1).

            Parameters
            ----------
            new_htrack : bool
                Set True to search the single star H-rich grids, or False to 
                search the He-rich grids. This is determined primarily 
                through the current `match_type`.

            list_for_matching : list
                This is a list that contains sublists of types, str, float, 
                str, float, float. The first sublist holds the names for 
                star attributes that will be used to quantify a match. The 
                second holds rescaling factors for those attributes. The 
                third holds DataScaler scaling methods (see 
                posydon.interpolation.data_scaling.DataScaler.fit() for more). 
                The fourth and fifth hold bounds for the minimization 
                algorithms on the stellar mass and age.

            Returns
            -------
            match_attr_names : list[str]
                The names of the attributes that will be used for matching. 
                These should be SingleStar STARPROPERTIES keys.
            
            rescale_facs : list[float]
                Contains normalization factors to be used for rescaling match 
                attribtue values. This should be the same length as 
                match_attr_names.

            bnds : list[list[float], list[float]]
                Lower and upper bounds on initial stellar mass and age to be 
                used in minimization algorithms.
            
            scalers : list[DataScaler object]
                DataScaler objects for each attribute, trained to rescale 
                quantities to the range (0, 1).
                
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
            for attr_name, method in zip(match_attr_names, scaler_methods):
                # check that attributes are allowed as matching attributes
                if attr_name not in self.root_keys:
                    raise AttributeError("Expected matching attribute "
                                         f"{attr_name} not "
                                         "added in root_keys list: "
                                         f"{self.root_keys}")
                # create attribute scalers
                scaler_of_attribute = self.scale(attr_name, new_htrack, method)
                scalers.append(scaler_of_attribute)

            return match_attr_names, rescale_facs, bnds, scalers

        def get_match_attrs(match_type="default", match_ok=True):

            """
                Gather the attribute names and values for matching
            a given star, plus the attribute bounds and scalings.
                
            Parameters
            ----------
            match_type : str
                A string that sets which matching list to use when obtaining 
                attributes. The defined options are:

                    "default" :     This gets the set default list of 
                                    attributes, which are dependent on the 
                                    star's evolutionary state.

                    "alt" :         This selects alternate lists of attributes 
                                    that are meant to ease the match finding 
                                    process, ignoring certain parameters like 
                                    stellar radius.

                    "alt_evolved" : This switches the grids that are searched 
                                    for a match for evolved stars, meaning 
                                    post-MS and stripped He stars as a last 
                                    ditch effort to find a match. 

            Returns
            -------
            match_attrs : tuple(list[str], list[float])
                A tuple that contains the match_attr_names, which are the 
                data column names of the match attributes and the 
                match_attr_vals, which are the associated values.

            scls_bnds : tuple(float, float, DataScaler)
                A tuple that contains the associated rescaling factors,
                boundaries, and DataScaler objects of the matching 
                attributes.
                    
            """

            # set matching metrics based on star state
            new_htrack = star.htrack
            if match_type == "default":
                if star.state in STAR_STATES_FOR_HMS_MATCHING:
                    match_list = self.list_for_matching_HMS
                elif star.state in STAR_STATES_FOR_postMS_MATCHING:
                    match_list = self.list_for_matching_postMS
                elif star.state in STAR_STATES_FOR_Hestar_MATCHING:
                    match_list = self.list_for_matching_HeStar
                else:
                    raise ValueError(f"{star.state} invalid for matching step")
                    
            elif match_type == "alt":
                if self.verbose:
                    print("(Now trying to match with alternative parameters)")

                # set alternative matching metrics based on star state
                if star.state in STAR_STATES_FOR_HMS_MATCHING:
                    match_list = self.list_for_matching_HMS_alternative
                elif star.state in STAR_STATES_FOR_postMS_MATCHING:
                    match_list = (self.list_for_matching_postMS_alternative)
                elif star.state in STAR_STATES_FOR_Hestar_MATCHING:
                    match_list = (self.list_for_matching_HeStar_alternative)
                else:
                    raise ValueError(f"{star.state} invalid for matching step")
                    
            elif match_type == "evolved_alt":
                Pwarn("Attempting to match an He-star with an H-rich grid or " \
                      "post-MS star with a stripped-He grid","EvolutionWarning")

                if star.state in STAR_STATES_FOR_HMS_MATCHING:
                    if self.verbose:
                        print(f"We cannot use match_type={match_type} " \
                                "since star is on the MS. Skipping...")
                    match_ok = False
                    match_attrs = (None, None, new_htrack, match_ok)
                    scls_bnds = (None, None, None)
                    return match_attrs, scls_bnds

                # if he star, try matching to H-rich grid
                elif star.state in STAR_STATES_FOR_Hestar_MATCHING:
                    new_htrack = True
                    match_list = self.list_for_matching_HeStar

                # if post MS star, try matching to He-rich grid
                elif star.state in STAR_STATES_FOR_postMS_MATCHING:
                    new_htrack = False
                    match_list = self.list_for_matching_postMS

            else:
                raise POSYDONError("In getting match attributes, match_type "
                        f"'{match_type}' is not a recognized option. Please "
                        "check how this happened and either create a new "
                        "option in get_match_attrs() or correct it.")

            # have approriate matching list, now get necessary attr properties
            attr_names, rescl_facs, bnds, sclrs = get_attr_props(new_htrack,
                                                                 match_list)
            attr_vals = get_attr_values(attr_names)

            match_attrs = (attr_names, attr_vals, new_htrack, match_ok)
            scls_bnds = (rescl_facs, bnds, sclrs)

            return match_attrs, scls_bnds

        def get_match_params(match_type, match_ok=True):
            """
                Get initial guess for minimization and matching specs. 
            The initial guess is found via the `get_root0()` function (see 
            that function for further details). The other specs gathered 
            by this function are the function arguments needed for the 
            square_difference function matching attribute boundaries 
            used by SciPy's minimization function.
            
            Parameters
            ----------
            match_type :str
                This string sets the type of matching to be done. This can 
                be 'default', 'alt', or 'alt_evolved'. This string will 
                dictate which matching parameters to use when considering 
                which stellar track is the closest match. (See 
                get_match_attrs() for more details.)

            match_ok : bool
                This boolean tracks whether something went wrong or not. If 
                this is False, the matching is considered failed and will 
                abort.
            
            Returns
            -------
            x0 : list[float]
                A list that contains the initial stellar mass and age of a 
                stellar track. Used as an initial guess for 
                `scipy.optimize.minimize()`.

            fnc_args : tuple(bool, list[str], list[float], list[DataScaler])
                This tuple contains the arguments that will be passed to the 
                minimization function. The minimization function uses the 
                square_difference function. See that function for more details.

            bnds : list[float]
                Lower and upper bounds on initial stellar mass and age to be 
                used in minimization algorithms.

            match_ok : bool
                This boolean tracks whether something went wrong or not. If 
                this is False, the matching is considered failed and will 
                abort.
            """
            # get names, values, bounds, and scalings of 
            # attributes to be matched
            match_attrs, scls_bnds = get_match_attrs(match_type, 
                                                     match_ok=match_ok)
            
            # unpack matching properties
            match_attr_names, match_attr_vals = match_attrs[:2]
            new_htrack, match_ok = match_attrs[2:]
            rescale_facs, bnds, scalers = scls_bnds

            if not match_ok:
                fnc_args = (new_htrack, None, None, None)
                return None, fnc_args, None, match_ok

            # Get closest matching point along track single star grids. x0 
            # contains the corresponding initial guess for mass and age
            x0 = get_root0(match_attr_names, match_attr_vals, 
                           new_htrack, rescale_facs=rescale_facs)
            fnc_args = (new_htrack, match_attr_names, match_attr_vals, scalers)

            return x0, fnc_args, bnds, match_ok

        def do_minimization(method='TNC', match_type="default"):
            """
                Perform minimization procedure with a given method and 
            match type, using SciPy's minimize function.

            Parameters
            ----------
            method : str
               This sets the desired solver to be used by SciPy's minimize 
               function. The default is 'TNC', but we also use 'Powell' as 
               an alternative.

            match_type :str
                This string sets the type of matching to be done. This can 
                be 'default', 'alt', or 'alt_evolved'. This string will 
                dictate which matching parameters to use when considering 
                which stellar track is the closest match. (See 
                get_match_attrs() for more details.)

            Returns
            -------
            sol : OptimizeResult object
                This is the OptimizeResult object returned by SciPy's 
                minimize function. It contains the mass and age of the 
                closest matching stellar track, found through 
                minimization.

            new_htrack : bool
                This boolean tracks whether the star should be matched to 
                the H- or He-rich single star grid. This can change in the 
                match finding process.
            
            match_ok : bool
                This boolean tracks whether something went wrong or not. If 
                this is False, the matching is considered failed and will 
                abort.
            """

            x0, fnc_args, bounds, match_ok = get_match_params(match_type)
            new_htrack = fnc_args[0]

            if not match_ok:
                return failed_sol, new_htrack, match_ok

            # Powell's method does not work with bound in SciPy version < 1.5.x
            if method == 'Powell' and SCIPY_VER < 1.5:
                if self.verbose:
                    print(f"Ignoring bounds because method = {method} "
                          f"does not accept them in SciPy v{SCIPY_VER}.x.")
                bounds = None

            # Run minimization method
            try:
                # Minimize sq, Euclidean dist. w/ Newton's method (TNC)
                sol = minimize(square_difference, x0, args=fnc_args,
                                method=method, bounds=bounds)
            
                # guard against NaN solutions, ensuring they will fail
                sol.fun = 1e99 if np.isnan(sol.fun) else sol.fun

                if self.verbose:
                    print (f"Matching attempt completed:"
                            f"\nBest solution: {np.abs(sol.fun)} " 
                            f"(tol = {self.matching_tolerance})"
                            f"\nsol.success = {sol.success}")
            except:
                raise NumericalError("SciPy numerical differentiation "
                                        "occurred outside boundary while "
                                        "matching to single star track")

            # check for failures:
            if (np.abs(sol.fun) > self.matching_tolerance):
                if self.verbose:
                    print (f"Matching result: FAILED"
                            "\nReason: Solution exceeds tolerance "
                            f"({np.abs(sol.fun)} > {self.matching_tolerance})")
                match_ok = False
            if (not sol.success):
                if self.verbose:
                    print (f"Matching result: FAILED"
                            "\nReason: Optimizer failed (sol.success = "
                            f"{sol.success})"
                            f"\nOptimizer termination reason: {sol.message}")
                match_ok = False

            if match_ok and self.verbose:
                print("Matching result: OK")

            return sol, new_htrack, match_ok

        # END SUBFUNCTION DEFIITIONS
        
        if self.verbose:
            print(DIVIDER_STR)

        # defined sequence of matching attempt types and methods
        match_sequence = {
                            1: {"type":"default", 
                                "method":"TNC"},
                            2: {"type":"default", 
                                "method":"Powell"},
                            3: {"type":"alt", 
                                "method":"TNC"},
                            4: {"type":"evolved_alt", 
                                "method":"TNC"}
                            }
        
        # dummy
        failed_sol = types.SimpleNamespace()
        failed_sol.success = False
        failed_sol.fun = 1e99
        failed_sol.message = "This is a dummy SciPy OptimizeResult object."
        
        best_sol = failed_sol
        for attempt_num in match_sequence.keys():

            match_type = match_sequence[attempt_num]["type"]
            method = match_sequence[attempt_num]["method"]

            if self.verbose:
                print(f"\nMatching attempt {attempt_num} started...")
                print(f"match_type = {match_type}")
                print(f"method = {method}")

            sol, new_htrack, match_ok = do_minimization(method, match_type)

            better_match = np.abs(sol.fun) < np.abs(best_sol.fun)
            sol_is_better =  better_match and sol.success
            # TODO: case for nan solution
            # if a better solution is found, update it
            if sol_is_better:
                star.htrack = new_htrack
                best_sol = sol

            # if we matched successfully and solution was better, we're done
            if sol_is_better and match_ok:
                break

        # if matching is still not successful, set result to NaN:
        exceeds_tol = np.abs(best_sol.fun) > self.matching_tolerance_hard
        if (exceeds_tol or not best_sol.success):
            if self.verbose:
                print("\nFinal matching result: FAILED")
                if (np.abs(best_sol.fun) > self.matching_tolerance_hard):
                    print ("\nReason: Solution exceeds hard tolerance "+\
                           f"({np.abs(best_sol.fun)} > "
                           f"{self.matching_tolerance_hard})")
                if (not best_sol.success):
                    print ("\nReason: Optimizer failed, "
                           "sol.success = {best_sol.success}"
                           "\nOptimizer termination reason: "
                           f"{best_sol.message}") 

            match_vals = np.array([np.nan, np.nan])

        # or else we found a solution
        else:
            if self.verbose:
                print("\nFinal matching result: SUCCESS"
                        f"\nBest solution within hard tolerance: "
                        f"{np.abs(best_sol.fun):.8f}", "<", 
                        self.matching_tolerance_hard)

            match_vals = best_sol.x

        return match_vals, best_sol


    def get_star_match_data(self, binary, star,
                            copy_prev_m0=None, copy_prev_t0=None):
        """
            Match a given component of a binary (i.e., a star) to a 
        single star model. This then creates and returns interpolator
        objects that may be used to calculate properties of the star
        as a function of time.

            In the case of a compact object, radius, mdot, and Idot are 
        set to zero. One may use another star, e.g., the companion of 
        the compact object to provide a mass and age.

        Parameters
        ----------
        binary : BinaryStar object
            A binary star object, containing the binary system's properties.

        star : SingleStar object
            A single star object that contains the star's properties.

        copy_prev_m0 : float
            A mass value that may be copied from another star in the case
            where the target star is a compact object

        copy_prev_t0 : float
            An age value that may be copied from another star in the case
            where the target star is a compact object

        Returns
        -------
        match_m0 : float
            Initial mass (in solar units) of the matched model.
        
        match_t0 : float
            Age (in years) of the matched model.

        """

        with np.errstate(all="ignore"):
            # get the initial m0, t0 track
            if binary.event == 'ZAMS' or binary.event == 'redirect_from_ZAMS':
                # ZAMS stars in wide (non-mass exchaging binaries) that are
                # directed to detached step at birth
                match_m0, match_t0 = star.mass, 0
            elif star.co:
                match_m0, match_t0 = copy_prev_m0, copy_prev_t0
            else:
                t_before_matching = time.time()
                # matching to single star grids (getting mass, age of 
                # closest track)
                match_m0, match_t0 = self.match_to_single_star(star)
                t_after_matching = time.time()

                if self.verbose:
                    match_tspan = t_after_matching-t_before_matching
                    print(f"Matching duration: {match_tspan:.6g} sec\n")

        # bad result
        if pd.isna(match_m0) or pd.isna(match_t0):
            star.interp1d = None
            return None, None

        if star.htrack:
            self.grid = self.grid_Hrich
        else:
            self.grid = self.grid_strippedHe

        # check if m0 is in the grid bounds
        outside_low = match_m0 < self.grid.grid_mass.min()
        outside_high = match_m0 > self.grid.grid_mass.max()
        if outside_low or outside_high:
            set_binary_to_failed(binary)
            raise MatchingError(f"The mass {match_m0} is out of "
                                "the single star grid range and "
                                "cannot be matched to a track.")

        # get/interpolate track values for requested mass match_m0
        get_track = self.grid.get

        max_time = binary.properties.max_simulation_time
        assert max_time > 0.0, "max_time is non-positive"

        # getting track of mass match_m0's age data
        age = get_track("age", match_m0)
        # max timelength of the track
        t_max = age.max()
        interp1d = dict()
        kvalue = dict()
        for key in self.KEYS[1:]:
            kvalue[key] = get_track(key, match_m0)
        try:
            for key in self.KEYS[1:]:
                if key in self.KEYS_POSITIVE:
                    positive = True
                    interp1d[key] = PchipInterpolator2(age, kvalue[key], 
                                                       positive=positive)
                else:
                    interp1d[key] = PchipInterpolator2(age, kvalue[key])
        except ValueError:
            i_bad = [None]
            while len(i_bad) != 0:
                i_bad = np.where(np.diff(age) <= 0)[0]
                age = np.delete(age, i_bad)
                for key in self.KEYS[1:]:
                    kvalue[key] = np.delete(kvalue[key], i_bad)

            for key in self.KEYS[1:]:
                if key in self.KEYS_POSITIVE:
                    positive = True
                    interp1d[key] = PchipInterpolator2(age, kvalue[key], 
                                                       positive=positive)
                else:
                    interp1d[key] = PchipInterpolator2(age, kvalue[key])

        interp1d["inertia"] = PchipInterpolator2(age, 
                                                kvalue["inertia"] / (const.msol * const.rsol**2))

        interp1d["Idot"] = PchipInterpolator2(age, 
                                              kvalue["inertia"] / (const.msol * const.rsol**2),
                                              derivative=True)

        interp1d["conv_env_turnover_time_l_b"] = PchipInterpolator2(
            age, kvalue['conv_env_turnover_time_l_b'] / const.secyer)

        interp1d["L"] = PchipInterpolator2(age, 10 ** kvalue["log_L"])
        interp1d["R"] = PchipInterpolator2(age, 10 ** kvalue["log_R"])
        interp1d["t_max"] = t_max
        interp1d["max_time"] = max_time
        interp1d["t0"] = match_t0
        interp1d["m0"] = match_m0

        if star.co:
            kvalue["mass"] = np.zeros_like(kvalue["mass"]) + star.mass
            kvalue["R"] = np.zeros_like(kvalue["log_R"])
            kvalue["mdot"] = np.zeros_like(kvalue["mdot"])
            interp1d["mass"] = PchipInterpolator2(age, kvalue["mass"])
            interp1d["R"] = PchipInterpolator2(age, kvalue["R"])
            interp1d["mdot"] = PchipInterpolator2(age, kvalue["mdot"])
            interp1d["Idot"] = PchipInterpolator2(age, kvalue["mdot"])

        # update star with interp1d object built from matched values
        star.interp1d = interp1d

        return match_m0, match_t0

    def calc_omega(self, star):
        """
        
            Calculate the spin of a star from its (pre-match) moment
        of inertia and angular momentum (or rotation rates). This is 
        required because we match a rotating model (from the binary grids) 
        to a non-rotating model (from the single star grids).

        Parameters
        ----------
        star : SingleStar object
            Star object containing the star properties.

        Returns
        -------
        omega_in_rad_per_year: float
            The rotation rate of the star in radians per year, 
            calculated from the star's (pre-match) angular momentum
            and moment of inertia.

        Warns
        -----
        InappropriateValueWarning
            If the pre-match rotation quantities are NaN or None and can 
            not calculate a the post-match rotation rate, we setting the 
            post-match rotation rate to zero.
        """

        log_total_J = star.log_total_angular_momentum
        total_MOI = star.total_moment_of_inertia
        omega_div_omega_c = star.surf_avg_omega_div_omega_crit
        omega = star.surf_avg_omega

        if (pd.notna(log_total_J) and pd.notna(total_MOI)):

            if self.verbose:
                print("Calculating post-match omega using " \
                      "pre-match angular momentum and moment of inertia")

            # the last factor converts rad/s to rad/yr
            omega_in_rad_per_yr = (10.0 ** log_total_J 
                                     / total_MOI * const.secyer)

        else:
            # we equate the secondary's initial omega to surf_avg_omega
            # (although the critical rotation should be improved to
            # take into account radiation pressure)
            if pd.notna(omega):

                if self.verbose:
                    print("Calculating post-match omega using " \
                          "pre-match surf_avg_omega")

                omega_in_rad_per_yr = omega * const.secyer                    

            elif pd.notna(omega_div_omega_c):

                if self.verbose:
                    print("Calculating post-match omega using " \
                          "pre-match surf_avg_omega_div_omega_crit")

                if pd.notna(star.log_R):
                    numerator = const.standard_cgrav * star.mass * const.msol
                    denominator = (10.0 ** (star.log_R) * const.rsol) ** 3 
                    omega_c = np.sqrt(numerator / denominator) # rad/s
                    omega_c *= const.secyer # rad/yr
                    omega_in_rad_per_yr = omega_div_omega_c * omega_c

                else:
                    radius_interp = star.interp1d["R"](star.interp1d["t0"])
                    mass_interp = star.interp1d["mass"](star.interp1d["t0"])

                    numerator = const.standard_cgrav * mass_interp * const.msol
                    denominator = (radius_interp * const.rsol) ** 3 
                    omega_c = np.sqrt(numerator / denominator) # rad/s
                    omega_c *= const.secyer # rad/yr
                    omega_in_rad_per_yr = omega_div_omega_c * omega_c

            else:
                omega_in_rad_per_yr = 0.0
                if self.verbose:
                    print("Could not calculate post-match omega, " \
                          "pre-match values are None or NaN.")
                    print("Pre-match rotation rates:")
                    print("surf_avg_omega = ", omega)
                    print("surf_avg_omega_div_omega_crit = ", 
                          omega_div_omega_c)
                    Pwarn("Setting (post-match) rotation rate to zero.", 
                          "InappropriateValueWarning")

        if self.verbose and omega is not None and (omega_in_rad_per_yr != 0):
            print("pre-match omega [rad/yr] = ", omega * const.secyer)
            print("calculated omega [rad/yr] = ", omega_in_rad_per_yr)
            pcdiff = 100.0*(omega_in_rad_per_yr-omega * const.secyer) \
                                                / omega_in_rad_per_yr
            omega_percent_diff = np.nan_to_num(pcdiff, 0)
            print("omega [rad/yr] % difference = ", 
                    f"{omega_percent_diff:.1f}%")

        return omega_in_rad_per_yr

    def update_rotation_info(self, primary, secondary):

        """
            Once we have matched to a non-rotating single star, we need to 
        calculate what the single star's spin should be, based the star's 
        rotation before matching. This calculates a new rotation rate for the 
        matched star from the previous moment of inertia and angular 
        momentum (or rotation rate in lieu of those). This also updates the 
        stars in the binary with the newly calculated values to reflect this.

        Parameters
        ----------
        primary: SingleStar object
            A single star object, representing the primary (more evolved) star 
            in the binary and containing its properties.
        
        secondary: SingleStar object
            A single star object, representing the secondary (less evolved) 
            star in the binary and containing its properties. 

        Returns
        -------
        omega0_pri: float 
            The rotation rate of the primary star in radians per year, 
            calculated from the star's (pre-match) angular momentum
            and moment of inertia.

        omega0_sec: float 
            The rotation rate of the secondary star in radians per year, 
            calculated from the star's (pre-match) angular momentum
            and moment of inertia.

        """

        omega0_sec = self.calc_omega(secondary)

        # recalculate rotation quantities using the newly calculated 
        # omega [rad/s] here. Need to be careful about case where 
        # omega was made/is 0 to avoid div by zero errors or None * float
        surf_avg_omega = secondary.surf_avg_omega
        if pd.notna(surf_avg_omega) and surf_avg_omega > 0.0:
            new_omega = omega0_sec / const.secyer # rad/s
            old_ratio = secondary.surf_avg_omega_div_omega_crit
            old_omega = surf_avg_omega
            old_omega_c = old_omega / old_ratio
            new_ratio = new_omega / old_omega_c
            new_lgJ = np.log10(new_omega * secondary.total_moment_of_inertia)
            setattr(secondary, "surf_avg_omega_div_omega_crit", new_ratio)
            setattr(secondary, "surf_avg_omega", new_omega)
            setattr(secondary, "log_total_angular_momentum", new_lgJ)

        else:
            secondary.surf_avg_omega_div_omega_crit = 0.0
            secondary.surf_avg_omega = 0.0
            secondary.log_total_angular_momentum = 0.0

        if self.primary_normal:
            omega0_pri = self.calc_omega(primary)

            surf_avg_omega = primary.surf_avg_omega
            if pd.notna(surf_avg_omega) and surf_avg_omega > 0.0:
                new_omega = omega0_pri / const.secyer # rad/s
                old_ratio = primary.surf_avg_omega_div_omega_crit
                old_omega = surf_avg_omega
                old_omega_c = old_omega / old_ratio
                new_ratio = new_omega / old_omega_c
                new_lgJ = np.log10(new_omega * primary.total_moment_of_inertia)
                setattr(primary, "surf_avg_omega_div_omega_crit", new_ratio)
                setattr(primary, "surf_avg_omega", new_omega)
                setattr(primary, "log_total_angular_momentum", new_lgJ)

            else:
                primary.surf_avg_omega_div_omega_crit = 0.0
                primary.surf_avg_omega = 0.0
                primary.log_total_angular_momentum = 0.0
        else:
            # omega of compact objects or massless remnant 
            # (won't be used for integration)
            omega0_pri = omega0_sec

        return omega0_pri, omega0_sec

    def do_matching(self, binary, step_name="step_match"):

        """
            Perform binary to single star grid matching. This is currently
        used when transitioning to detached star evolution from binary but 
        may be used in other steps. This performs several actions:

            a. Determines which star is primary/secondary in the evolution
               and their evolutionary states. If evolvable, matching will 
               proceed.
            b. Match either one or both stars to a (non-rotating) single 
               star evolution track.
            c. Calculate the matched star's rotation using the pre-match 
               step's angular momentum-related properties.
            d. Returns the primary/secondary stars with interpolator 
               objects that may be used to calculate quantities along the 
               time series of each star for further evolution.
        
        Parameters
        ----------
        binary : BinaryStar object
            A binary star object, containing the binary system's properties.

        step_name : str
            If self.record_matching is True, then the matched quantities of 
            the star will be appended to the history. This is a string that 
            can be used as a custom label in the BinaryStar object's history, 
            meant to indicate the relevant evolution step's name. This should 
            normally match the name of the step in which the matching was 
            made, e.g., "step_detached".

        Returns
        -------
        primary_out : tuple(SingleStar, PchipInterpolator, float)
            The first element is the SingleStar object of the primary 
            (more evolved) star. The second element is the primary's 
            PchipInterpolator object, used to interpolate values along 
            this star's time series. The third is the primary's 
            calculated (post-match) rotation rate, using angular 
            momentum-related quantites from the pre-match step.
        
        secondary_out :  tuple(SingleStar, PchipInterpolator, float)
            The first element is the SingleStar object of the secondary 
            (less evolved) star. The second element is the secondary's 
            PchipInterpolator object, used to interpolate values along 
            this star's time series. The third is the secondary's 
            calculated (post-match) rotation rate, using angular 
            momentum-related quantites from the pre-match step.
        
        only_CO : bool
            A boolean indicating whether the binary system contains only a 
            single compact object (True) or not (False). As in the detached 
            step, it may be desirable to use this flag to exit an evolution 
            step, as a single compact obejct (point mass) can not be evolved 
            further.

        Raises
        ------
        ValueError
            If the `primary_normal` and `primary_not_normal` flags are 
            both determined to be False. One or the other should be True.

        MatchingError
            If the stellar matching to a single star model fails, or the 
            PchipInterpolator object returned from matching is None.

        """

        # determine star states for matching
        primary, secondary, only_CO = self.determine_star_states(binary)
        if only_CO:
            return (None, None, None), (None, None, None), only_CO

        # record which star we performed matching on for reporting purposes
        self.matched_s1 = False
        self.matched_s2 = False

        # get the matched data of binary components
        # match secondary:
        m0, t0 = self.get_star_match_data(binary, secondary)
        # record which star got matched
        if secondary == binary.star_2:
            self.matched_s2 = True
        elif secondary == binary.star_1:
            self.matched_s1 = True

        # primary is a CO or massless remnant, or else it is "normal"
        # TODO: should these be star properties? also, do we only really need one?
        has_non_existent = binary.non_existent_companion in [1,2]
        all_exist = binary.non_existent_companion == 0
        self.primary_not_normal = primary.co or has_non_existent
        self.primary_normal = not primary.co and all_exist

        if self.primary_not_normal:
            # copy the secondary star except mass which is of the primary,
            # and radius, mdot, Idot = 0
            self.get_star_match_data(binary, primary,
                                     copy_prev_m0 = m0,
                                     copy_prev_t0 = t0)
        elif self.primary_normal:
            # match primary
            self.get_star_match_data(binary, primary)

            if primary == binary.star_1:
                self.matched_s1 = True
            elif primary == binary.star_2:
                self.matched_s2 = True
        else:
            raise ValueError("During matching, the primary should either be "
                             "normal (stellar object) or "
                             "not normal (a CO or nonexistent companion).",
                            f"\nprimary.co = {primary.co}",
                            "\nnon_existent_companion = "
                            f"{binary.non_existent_companion}",
                            "\ncompanion_1_exists = "
                            f"{binary.companion_1_exists}",
                            "\ncompanion_2_exists = "
                            f"{binary.companion_2_exists}")


        if (secondary.interp1d == None) or (primary.interp1d == None):
            failed_state = binary.state
            set_binary_to_failed(binary)
            raise MatchingError("Grid matching failed for " 
                                f"{failed_state} binary.")

        # recalculate rotation quantities after matching
        omega0_pri, omega0_sec = self.update_rotation_info(primary, secondary)

        # update binary history with matched values
        # (only shown in history if record_matching = True)
        # (this gets overwritten after detached evolution)
        self.update_star_properties(secondary, secondary.htrack)
        if self.primary_normal:
            self.update_star_properties(primary, primary.htrack)

        if self.record_matching:
            # append matching information as a part of step_detached
            binary.step_names.append(step_name)
            if self.matched_s1 and self.matched_s2:
                binary.event = "Match12"
            elif self.matched_s1:
                binary.event = "Match1"
            elif self.matched_s2:
                binary.event = "Match2"

            binary.append_state()

        primary.omega0 = omega0_pri
        secondary.omega0 = omega0_sec

        return primary, secondary, only_CO

    def determine_star_states(self, binary):

        """
            Determines which star is primary (further evolved) and which is 
        secondary (less evolved). Determines whether stars should be 
        matched to the H- or He-rich grid, whether they exist, or if they
        are compact objects/massless remnants. This is used to determine
        how to match the stars.

        Parameters
        ----------
        binary: BinaryStar object
            A binary star object, containing the binary system's properties.

        Returns
        -------
        primary : SingleStar object
            A single star object, representing the primary (more evolved) star 
            in the binary and containing its properties.
        
        secondary : SingleStar object
            A single star object, representing the secondary (less evolved) 
            star in the binary and containing its properties. 
        
        only_CO : bool
            A boolean indicating whether the binary system contains only a 
            single compact object (True) or not (False). As in the detached 
            step, it may be desirable to use this flag to exit an evolution 
            step, as a single compact obejct (point mass) can not be evolved 
            further.

        Raises
        ------
        POSYDONError
            If there is no star to evolve (both are massless remnants), then 
            evolution can no continue and detached step should not have been 
            called.

        ValueError
            If the State of star 1 or 2 is not recognized.

        POSYDONError
            If the `non_existent_companion` of the binary is determined to 
            be not equal to 0 (both stars exist), 1 (only star 2 exists), 
            or 2 (only star 1 exists), something has gone wrong.
        
        """

        only_CO = False

        # update BinaryStar instance with attributes storing whether star 1/2 
        # exists or not
        binary.check_who_exists()
        if binary.non_existent_companion == -1:
            raise POSYDONError("There is no star to evolve. Who summoned me?")

        # where both stars exist. The primary is a potential compact object, or
        # the more evolved star
        s_arr = np.array([binary.star_1, binary.star_2])
        s_CO = np.array([s.state in STAR_STATES_CO for s in s_arr])
        s_H = np.array([s.state in STAR_STATES_H_RICH for s in s_arr])
        s_He = np.array([s.state in STAR_STATES_FOR_Hestar_MATCHING for s in s_arr])
        s_massless = np.array([s.state == "massless_remnant" for s in s_arr])
        s_valid = s_H | s_He | s_CO | s_massless    # states considered here
        s_htrack = s_H & ~(s_CO)   # only true if h rich and not a CO

        # check if star states are recognizable
        if any(~s_valid):
            raise ValueError(f"Star1 state: {binary.star_1.state} "
                                f"(valid: {s_valid[0]})\n"
                                f"Star2 state: {binary.star_2.state} "
                                f"(valid: {s_valid[1]})\n")

        if binary.non_existent_companion == 0: # both stars exist, detached step of a binary

            # states match, either both H stars or both He stars and not any COs
            # prevents He+CO going into here.
            if (all(s_valid) and (all(s_htrack) or all(~s_htrack)) and not (any(s_CO))):
                primary = s_arr[0]
                primary.co = s_CO[0]
                primary.htrack = s_htrack[0]

                secondary = s_arr[1]
                secondary.co = s_CO[1]
                secondary.htrack = s_htrack[1]
            # states mismatch, one is an H star and the other is an He star
            elif (all(s_valid) and not any(s_CO)):
                htrack_mask = s_htrack == True

                primary = s_arr[~htrack_mask].item()
                primary.co = s_CO[~htrack_mask].item()
                primary.htrack = s_htrack[~htrack_mask].item()

                secondary = s_arr[htrack_mask].item()
                secondary.co = s_CO[htrack_mask].item()
                secondary.htrack = s_htrack[htrack_mask].item()
            # states mismatch, one is a CO and other is an H or He star
            elif (all(s_valid) and (any(s_CO) and not all(s_CO))):
                CO_mask = s_CO == True

                primary = s_arr[CO_mask].item()
                primary.co = s_CO[CO_mask].item()
                primary.htrack = s_htrack[~CO_mask].item()

                secondary = s_arr[~CO_mask].item()
                secondary.co = s_CO[~CO_mask].item()
                secondary.htrack = s_htrack[~CO_mask].item()

            else:
                # both stars are compact objects, should redirect to step_dco
                only_CO = True
                return None, None, only_CO

        # In case a star is a massless remnant:
        # We force primary.co = True for all isolated evolution
        # where the primary does not exist (is a massless remnant)
        # and the secondary is the one evolving

        # star 1 is a massless remnant, only star 2 exists
        elif binary.non_existent_companion == 1:
            # massless remnant
            primary = s_arr[0]
            primary.co = True
            primary.htrack = s_htrack[1]

            secondary = s_arr[1]
            secondary.htrack = s_htrack[1]
            secondary.co = s_CO[1]
            if secondary.co:
                only_CO = True
                return primary, secondary, only_CO

        # star 2 is a massless remnant, only star 1 exists
        elif binary.non_existent_companion == 2:
            primary = s_arr[1]
            primary.co = True
            primary.htrack = s_htrack[0]

            secondary = s_arr[0]
            secondary.htrack = s_htrack[0]
            secondary.co = s_CO[0]
            if secondary.co:
                only_CO = True
                return primary, secondary, only_CO

        else:
            raise POSYDONError("non_existent_companion = "
                f"{binary.non_existent_companion} (should be -1, 0, 1, or 2).")

        return primary, secondary, only_CO

    def update_star_properties(self, star, htrack):

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
        star : SingleStar object
            A single star object that contains the star's properties.
        
        htrack : bool
            A boolean that specifies whether the star would be found in the 
            hydrogen rich single star grid or not (in which case it is
            matched to the helium rich single star grid).
        """

        # initial mass and age at point of closest match
        m0 = star.interp1d["m0"]
        t0 = star.interp1d["t0"]

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


    def get_star_final_values(self, star):
        """
            This updates the final values of a SingleStar object,
        given an initial stellar mass `m0`, typically found from 
        matching to a single star track.

        Parameters
        ----------
        star : SingleStar object
            A single star object that contains the star's properties.
        
        htrack : bool
            A boolean that specifies whether the star would be found in the 
            hydrogen rich single star grid or not (in which case it is
            matched to the helium rich single star grid).

        m0 : float
            Initial stellar mass (in solar units) of the single star track 
            that we will grab values from and update `star` with.

        """

        m0 = star.interp1d["m0"]
        htrack = star.htrack

        grid = self.grid_Hrich if htrack else self.grid_strippedHe
        get_final_values = grid.get_final_values

        for key in self.final_keys:
            setattr(star, key, get_final_values('S1_%s' % (key), m0))

    def get_star_profile(self, star):
        """
            This updates the stellar profile of a SingleStar object,
        given an initial stellar mass `m0`, typically found from 
        matching to a single star track. The profile of the SingleStar 
        object is updated to become the profile of the (matched) single 
        star track.

        Parameters
        ----------
        star : SingleStar object
            A single star object that contains the star's properties.
        
        htrack : bool
            A boolean that specifies whether the star would be found in the 
            hydrogen rich single star grid or not (in which case it is
            matched to the helium rich single star grid).

        m0 : float
            Initial stellar mass (in solar units) of the single star track 
            that we will grab values from and update `star` with.

        """

        m0 = star.interp1d["m0"]
        htrack = star.htrack

        grid = self.grid_Hrich if htrack else self.grid_strippedHe
        get_profile = grid.get_profile
        profile_new = np.array(get_profile('mass', m0)[1])

        for i in self.profile_keys:
            profile_new[i] = get_profile(i, m0)[0]
        profile_new['omega'] = star.surf_avg_omega

        star.profile = profile_new
