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

class TrackMatcher:
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
                        various quantities between the previous evolution step and 
                        a stellar evolution track. 

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
        a single star grid. Assigned after loading a grid and before storing 
        matching metrics.

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
        
        See

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
        Mapping a combination of (key, htrack, scaling_method) to a pre-trained
        DataScaler instance.

    final_keys : tuple
        Contains keys for final value interpolation.

    profile_keys : tuple
        Contains keys for profile interpolation.

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

        # ==================================================================================

        self.metallicity = convert_metallicity_to_string(metallicity)
        self.matching_method = matching_method

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

        self.record_matching = record_matching

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
        
        htrack : bool
            A boolean that specifies whether the star would be found in the 
            hydrogen rich single star grid or not (in which case it is
            matched to the helium rich single star grid).

        Returns
        -------
        m0 : float
            Mass (in solar units) of the matched model
        
        t0 : float
            Age (in years) of the matched model

        htrack : bool
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
            attr_names : list[str]
                This list contains strings that are the names of 
                data columns (attributes) to be used as matching 
                metrics.

            star : SingleStar object
                This is a SingleStar object, typically representing 
                a member of a binary star system.

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
        
        def sq_diff_function(x):

            """
                This is a wrapper that calls a function to calculate the 
            square difference between specified star attributes for matching.
            When matching_method=`minimize`, this is the function that is 
            used by the minimization algorithm.

            Parameters
            ----------
            x : list[float]
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
            list_for_matching : list
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
            match_attr_names : list[str]
                The names of the attributes that will be used for matching. These 
                should be SingleStar STARPROPERTIES keys.
            
            rescale_facs : list[float]
                Contains normalization factors to be used for rescaling match attribtue 
                values. This should be the same length as match_attr_names.
                

            bnds : list[list[float], list[float]]
                Lower and upper bounds on initial stellar mass and age to be used 
                in minimization algorithms.
            
            scalers : list[DataScaler object]
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
            if star.state in STAR_STATES_FOR_HMS_MATCHING:

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
            if star.state in STAR_STATES_FOR_HMS_MATCHING:
                list_for_matching = self.list_for_matching_HMS
            elif star.state in STAR_STATES_FOR_postMS_MATCHING:
                list_for_matching = self.list_for_matching_postMS
            elif star.state in STAR_STATES_FOR_Hestar_MATCHING:
                list_for_matching = self.list_for_matching_HeStar

            # begin matching attempts
            divider_str = "_______________________________________________________________________\n"
            if self.verbose:
                print(divider_str)

            print("Initial matching attempt started...")
            
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
                        print (f"\nInitial matching attempt: FAILED"
                               f"\nReason: Solution exceeds tolerance ({np.abs(best_sol.fun)} > {matching_tolerance})")
                    if (not best_sol.success):
                        print (f"Initial matching attempt: FAILED"
                               f"\nReason: Optimizer failed (sol.success = {best_sol.success})"
                               f"\nOptimizer termination reason: {best_sol.message}")                        

                    print("\nAlternative matching started (2nd attempt)...")
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
                        print (f"Alternative matching (2nd attempt): FAILED"
                               f"\nReason: Solution exceeds tolerance ({np.abs(best_sol.fun)} > {matching_tolerance})")
                    if (not best_sol.success):
                        print (f"Alternative matching (2nd attempt): FAILED"
                               f"\nReason: Optimizer failed (sol.success = {best_sol.success})"
                               f"\nOptimizer termination reason: {best_sol.message}")

                    print("\nAlternative matching started (3rd attempt)...")
                    print("(Now trying to match with alternative parameters)")     
                      
                # set alternative matching metrics based on star state
                if star.state in STAR_STATES_FOR_HMS_MATCHING:
                    list_for_matching = self.list_for_matching_HMS_alternative
                elif star.state in STAR_STATES_FOR_postMS_MATCHING:
                    list_for_matching = (self.list_for_matching_postMS_alternative)
                elif star.state in STAR_STATES_FOR_Hestar_MATCHING:
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
                        print (f"Alternative matching (3rd attempt): FAILED"
                               f"\nReason: Solution exceeds tolerance ({np.abs(best_sol.fun)} > {matching_tolerance})")
                    if (not best_sol.success):
                        print (f"Alternative matching (3rd attempt): FAILED"
                               f"\nReason: Optimizer failed (sol.success = {best_sol.success})"
                               f"\nOptimizer termination reason: {best_sol.message}")  
                
                # if post-MS or stripped He star
                if (star.state in STAR_STATES_FOR_Hestar_MATCHING
                    or star.state in STAR_STATES_FOR_postMS_MATCHING):

                    if self.verbose:
                        print("\nAlternative matching started (4th attempt)...")
                        print("(Now trying to match He-star or post-MS star to a different grid)")

                    Pwarn("Attempting to match an He-star with an H-rich grid or post-MS star with a"
                          " stripped-He grid", "EvolutionWarning")
                        
                    if star.state in STAR_STATES_FOR_Hestar_MATCHING:
                        new_htrack = True
                        list_for_matching = self.list_for_matching_HeStar

                    elif star.state in STAR_STATES_FOR_postMS_MATCHING:
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
                            print (f"Alternative matching (4th attempt) completed:"
                                   f"\nBest solution: {np.abs(best_sol.fun)} (tol = {matching_tolerance})"
                                   f"\nsol.success = {best_sol.success}")
                    except:
                        raise NumericalError("SciPy numerical differentiation occured outside boundary "
                                             "while matching to single star track")

            # if matching is still not successful, set result to NaN:
            if (np.abs(best_sol.fun) > matching_tolerance_hard or not best_sol.success):
                if self.verbose:
                    print("\nFinal matching result: FAILED")#,
                          #np.abs(best_sol.fun), ">", matching_tolerance_hard)
                    if (np.abs(best_sol.fun) > matching_tolerance_hard):
                        print ("\nReason: Solution exceeds hard tolerance "+\
                               f"({np.abs(best_sol.fun)} > {matching_tolerance_hard})")
                    if (not best_sol.success):
                        print (f"\nReason: Optimizer failed, sol.success = {best_sol.success}"
                               f"\nOptimizer termination reason: {best_sol.message}") 

                initial_track_vals = (np.nan, np.nan)
                
            # or else we found a solution
            else:
                if self.verbose:
                    print("\nFinal matching result: SUCCESS"
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

            # done with matching attempts
            print(divider_str)

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

                Return
                -------
                interp1d : dict
                    A dictionary of scipy.interpolate._cubic.PchipInterpolator 
                    objects used to calculate star properties corresponding to the 
                    matched single star track.

                """

                htrack = star.htrack

                with np.errstate(all="ignore"):
                    # get the initial m0, t0 track
                    if binary.event == 'ZAMS' or binary.event == 'redirect_from_ZAMS':
                        # ZAMS stars in wide (non-mass exchaging binaries) that are
                        # directed to detached step at birth
                        m0, t0 = star.mass, 0
                    elif star.co:
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
                for key in self.KEYS[1:]:
                    kvalue[key] = get_track(key, m0)
                try:
                    for key in self.KEYS[1:]:
                        if key in self.KEYS_POSITIVE:
                            positive = True
                            interp1d[key] = PchipInterpolator2(age, kvalue[key], positive=positive)
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

                if star.co:
                    kvalue["mass"] = np.zeros_like(kvalue["mass"]) + star.mass
                    kvalue["R"] = np.zeros_like(kvalue["log_R"])
                    kvalue["mdot"] = np.zeros_like(kvalue["mdot"])
                    interp1d["mass"] = PchipInterpolator(age, kvalue["mass"])
                    interp1d["R"] = PchipInterpolator(age, kvalue["R"])
                    interp1d["mdot"] = PchipInterpolator(age, kvalue["mdot"])
                    interp1d["Idot"] = PchipInterpolator(age, kvalue["mdot"])

                return interp1d, m0, t0

    def calc_omega(self, star, prematch_rotation, interp1d):
        """
        
            Calculate the spin of a star from its (pre-match) moment
        of inertia and angular momentum (or rotation rates). This is 
        required because we match a rotating model (from the binary grids) 
        to a non-rotating model (from the single star grids).

        Parameters
        ----------
        star : SingleStar object
            Star object containing the star properties.

        prematch_rotation: dict 
            A dictionary containing the pre-match values of 
            angular momentum, moment of inertia, omega/omega_crit, 
            and omega

        interp1d: dict
            A dictionary of scipy.interpolate._cubic.PchipInterpolator 
            objects. Here, used to calculate radius and mass of the star 
            in the event that star.log_R is NaN.

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

        log_total_angular_momentum = prematch_rotation['log_total_angular_momentum']
        total_moment_of_inertia = prematch_rotation['total_moment_of_inertia']
        surf_avg_omega_div_omega_crit = prematch_rotation['surf_avg_omega_div_omega_crit']
        surf_avg_omega = prematch_rotation['surf_avg_omega']

        if (log_total_angular_momentum is not None
                and total_moment_of_inertia is not None
                and pd.notna(log_total_angular_momentum)
                and pd.notna(total_moment_of_inertia)):
            
            if self.verbose:
                print("Calculating post-match omega using angular momentum and moment of inertia")

            # the last factor converts rad/s to rad/yr
            omega_in_rad_per_year = (10.0 ** log_total_angular_momentum
                                        / total_moment_of_inertia * const.secyer)  
            
        else:
            # we equate the secondary's initial omega to surf_avg_omega
            # (although the critical rotation should be improved to
            # take into account radiation pressure)
            if pd.notna(surf_avg_omega):

                if self.verbose:
                    print("Calculating post-match omega using surf_avg_omega")

                omega_in_rad_per_year = surf_avg_omega * const.secyer                    
                    
            elif pd.notna(surf_avg_omega_div_omega_crit):
                
                if self.verbose:
                    print("Calculating post-match omega using surf_avg_omega_div_omega_crit")

                if pd.notna(star.log_R):
                    omega_in_rad_per_year = (surf_avg_omega_div_omega_crit * np.sqrt(
                                                const.standard_cgrav * star.mass * const.msol
                                                / ((10.0 ** (star.log_R) * const.rsol) ** 3)) * const.secyer)
                    
                else:
                    
                    radius_to_be_used = interp1d["R"](interp1d["t0"])
                    mass_to_be_used = interp1d["mass"](interp1d["t0"])

                    omega_in_rad_per_year = (surf_avg_omega_div_omega_crit * np.sqrt(
                                                const.standard_cgrav * mass_to_be_used * const.msol
                                                / ((radius_to_be_used * const.rsol) ** 3)) * const.secyer)
                
            else:
                omega_in_rad_per_year = 0.0
                if self.verbose:
                    print("Could not calculate post-match omega, pre-match values are None or NaN.")
                    print("Pre-match rotation rates:")
                    print("surf_avg_omega = ", surf_avg_omega)
                    print("surf_avg_omega_div_omega_crit = ", surf_avg_omega_div_omega_crit)
                    Pwarn("Setting (post-match) rotation rate to zero.", "InappropriateValueWarning")
                    
                    
        if self.verbose:
            print("pre-match omega_in_rad_per_year = ", surf_avg_omega * const.secyer)
            print("calculated omega_in_rad_per_year = ", omega_in_rad_per_year)
            omega_percent_diff = np.nan_to_num(100.0*(omega_in_rad_per_year-surf_avg_omega*const.secyer) \
                                                / omega_in_rad_per_year, 0)
            print("omega_in_rad_per_year % difference = ", 
                    f"{omega_percent_diff:.1f}%")

        return omega_in_rad_per_year
        
    def update_rotation_info(self, primary, secondary, interp1d_pri, interp1d_sec):

        """
            Once we have matched to a non-rotating single star, we need to calculate
        what the single star's spin should be, based the star's rotation before matching. 
        This calculates a new rotation rate for the matched star from the previous moment 
        of intertia and angular momentum (or rotation rate in lieu of those). This also 
        updates the stars in the binary with the newly calculated values to reflect this.

        Parameters
        ----------
        primary: SingleStar object
            A single star object, representing the primary (more evolved) star 
            in the binary and containing its properties.
        
        secondary: SingleStar object
            A single star object, representing the secondary (less evolved) star 
            in the binary and containing its properties. 

        Returns
        -------
        omega0_pri: PchipInterpolator 
            The rotation rate of the primary star in radians per year, 
            calculated from the star's (pre-match) angular momentum
            and moment of inertia.

        omega0_sec: PchipInterpolator 
            The rotation rate of the secondary star in radians per year, 
            calculated from the star's (pre-match) angular momentum
            and moment of inertia.

        """

        prematch_rotation_sec = {"log_total_angular_momentum": secondary.log_total_angular_momentum,
                                "total_moment_of_inertia": secondary.total_moment_of_inertia,
                                "surf_avg_omega_div_omega_crit": secondary.surf_avg_omega_div_omega_crit,
                                "surf_avg_omega": secondary.surf_avg_omega}
        
        prematch_rotation_pri = {"log_total_angular_momentum": primary.log_total_angular_momentum,
                                "total_moment_of_inertia": primary.total_moment_of_inertia,
                                "surf_avg_omega_div_omega_crit": primary.surf_avg_omega_div_omega_crit,
                                "surf_avg_omega": primary.surf_avg_omega}

        omega0_sec = self.calc_omega(secondary, prematch_rotation_sec, interp1d_sec)

        # recalculate rotation quantities using the newly calculated omega
        secondary.surf_avg_omega_div_omega_crit = (omega0_sec / const.secyer / secondary.surf_avg_omega) \
                                                * secondary.surf_avg_omega_div_omega_crit
        secondary.surf_avg_omega = omega0_sec / const.secyer
        secondary.log_total_angular_momentum = np.log10((omega0_sec / const.secyer) \
                                                        * secondary.total_moment_of_inertia)

        if self.primary_normal:
            omega0_pri = self.calc_omega(primary, prematch_rotation_pri, interp1d_pri)

            primary.surf_avg_omega_div_omega_crit = (omega0_pri / const.secyer / primary.surf_avg_omega) \
                                            * primary.surf_avg_omega_div_omega_crit
            primary.surf_avg_omega = omega0_pri / const.secyer
            primary.log_total_angular_momentum = np.log10((omega0_pri / const.secyer) \
                                                                * primary.total_moment_of_inertia)
        else:
            # omega of compact objects or massless remnant 
            # (won't be used for integration)
            omega0_pri = omega0_sec
        
        return omega0_pri, omega0_sec

    def do_matching(self, binary):
        
        """
            Perform binary to single star grid matching. This is currently
        used when transitioning to detached star evolution from binary.
        
        Parameters
        ----------
        binary : BinaryStar object
            A binary star object, containing the binary system's properties.

        primary : SingleStar object
            A single star object, representing the primary (more evolved) star 
            in the binary and containing its properties.
        
        secondary : SingleStar object
            A single star object, representing the secondary (less evolved) star 
            in the binary and containing its properties. 

        Returns
        -------
        interp1d_pri : PchipInterpolator object
            Interpolator object used to interpolate star properties and
            simulate detached evolution for the primary star.

        interp1d_sec : PchipInterpolator object
            Interpolator object used to interpolate star properties and
            simulate detached evolution for the secondary star.

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
        interp1d_sec, m0, t0 = self.get_star_match_data(binary, secondary)
        # record which star got matched
        if secondary == binary.star_2:
            self.matched_s2 = True
        elif secondary == binary.star_1:
            self.matched_s1 = True

        # primary is a CO or massless remnant, or else it is "normal"
        # TODO: should these be star properties? also, do we only really need one?
        self.primary_not_normal = (primary.co) or (binary.non_existent_companion in [1,2])
        self.primary_normal = (not primary.co) and binary.non_existent_companion == 0

        if self.primary_not_normal:
            # copy the secondary star except mass which is of the primary,
            # and radius, mdot, Idot = 0
            interp1d_pri = self.get_star_match_data(binary, primary, 
                                                            copy_prev_m0 = m0, 
                                                            copy_prev_t0 = t0)[0]
        elif self.primary_normal:
            # match primary
            interp1d_pri = self.get_star_match_data(binary, primary)[0]

            if primary == binary.star_1:
                self.matched_s1 = True
            elif primary == binary.star_2:
                self.matched_s2 = True
        else:
            raise ValueError("During matching, the primary should either be normal (stellar object) or ",
                            "not normal (a CO or nonexistent companion).",
                            f"\nprimary.co = {primary.co}",
                            f"\nnon_existent_companion = {binary.non_existent_companion}",
                            f"\ncompanion_1_exists = {binary.companion_1_exists}",
                            f"\ncompanion_2_exists = {binary.companion_2_exists}")


        if interp1d_sec is None or interp1d_pri is None:
            failed_state = binary.state
            set_binary_to_failed(binary)
            raise MatchingError(f"Grid matching failed for {failed_state} binary.")   

        # recalculate rotation quantities after matching
        omega0_pri, omega0_sec = self.update_rotation_info(primary, secondary, 
                                                           interp1d_pri, interp1d_sec)
        
        # update binary history with matched values (only shown in history if record_match = True)
        # (this gets overwritten after detached evolution)
        self.update_star_properties(secondary, secondary.htrack, 
                                    interp1d_sec["m0"], interp1d_sec["t0"])
        if self.primary_normal:
            self.update_star_properties(primary, primary.htrack, 
                                        interp1d_pri["m0"], interp1d_pri["t0"]) 

        if self.record_matching:
            # append matching information as a part of step_detached
            binary.step_names.append("step_detached")
            if self.matched_s1 and self.matched_s2:
                binary.event = "Match1,2"
            elif self.matched_s1:
                binary.event = "Match1"
            elif self.matched_s2:
                binary.event = "Match2"

            binary.append_state()

        primary_out = (primary, interp1d_pri, omega0_pri)
        secondary_out = (secondary, interp1d_sec, omega0_sec)

        return primary_out, secondary_out, only_CO
    
    def determine_star_states(self, binary):

        """
            Determines which star is primary (further evolved) and which is 
        secondary (less evolved). Determines whether stars should be 
        matched to the H- or He-rich grid, whether they exist, or if they
        are compact objects/massless remnants. THis is used to determine
        how to match the stars and also treat their detached evolution.

        Parameters
        ----------
        binary: BinaryStar object
            A binary star object, containing the binary system's properties.

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

        # update BinaryStar instance with attributes storing whether star 1/2 exists or not            
        binary.check_who_exists()
        if binary.non_existent_companion == -1:
                raise POSYDONError("There is no star to evolve. Who summoned me?")

        if binary.non_existent_companion == 0: # both stars exist, detached step of a binary 
            # where both stars exist. The primary is a potential compact object, or 
            # the more evolved star

            # star 1 is a CO, star 2 is H-rich
            if (binary.star_1.state in STAR_STATES_CO and
                binary.star_2.state in STAR_STATES_H_RICH):

                secondary = binary.star_2
                secondary.htrack = True
                secondary.co = False

                primary = binary.star_1
                primary.htrack = secondary.htrack
                primary.co = True

            # star 1 is a CO, star 2 is an He star
            elif (binary.star_1.state in STAR_STATES_CO and
                binary.star_2.state in STAR_STATES_FOR_Hestar_MATCHING):
                
                secondary = binary.star_2
                secondary.htrack = False
                secondary.co = False

                primary = binary.star_1
                primary.htrack = secondary.htrack
                primary.co = True

            # star 1 is H-rich, star 2 is a CO
            elif (binary.star_1.state in STAR_STATES_H_RICH and
                binary.star_2.state in STAR_STATES_CO):

                secondary = binary.star_1
                secondary.htrack = True
                secondary.co = False

                primary = binary.star_2
                primary.htrack = secondary.htrack
                primary.co = True

            # star 1 is an He star, star 2 is a CO
            elif (binary.star_1.state in STAR_STATES_FOR_Hestar_MATCHING and
                binary.star_2.state in STAR_STATES_CO):

                secondary = binary.star_1
                secondary.htrack = False
                secondary.co = False

                primary = binary.star_2
                primary.htrack = secondary.htrack
                primary.co = True
            
            # star 1 is H-rich, star 2 is H-rich
            elif (binary.star_1.state in STAR_STATES_H_RICH and
                binary.star_2.state in STAR_STATES_H_RICH):
                
                secondary = binary.star_2
                secondary.htrack = True
                secondary.co = False

                primary = binary.star_1
                primary.htrack = True
                primary.co = False

            # star 1 is an He star, star 2 is H-rich
            elif (binary.star_1.state in STAR_STATES_FOR_Hestar_MATCHING and
                binary.star_2.state in STAR_STATES_H_RICH):
                
                secondary = binary.star_2
                secondary.htrack = True
                secondary.co = False

                primary = binary.star_1
                primary.htrack = False
                primary.co = False

            # star 1 is H-rich, star 2 is an He star
            elif (binary.star_1.state in STAR_STATES_H_RICH and
                binary.star_2.state in STAR_STATES_FOR_Hestar_MATCHING):

                secondary = binary.star_1
                secondary.htrack = True
                secondary.co = False

                primary = binary.star_2
                primary.htrack = False
                primary.co = False

            # star 1 is an He star, star 2 is an He star
            elif (binary.star_1.state in STAR_STATES_FOR_Hestar_MATCHING and
                binary.star_2.state in STAR_STATES_FOR_Hestar_MATCHING):
                
                secondary = binary.star_2
                secondary.htrack = False
                secondary.co = False

                primary = binary.star_1
                primary.htrack = False
                primary.co = False

            else:
                raise ValueError(f"State {binary.star_1.state} is not recognized!")

        # In case a star is a massless remnant:
        # We force primary.co = True for all isolated evolution
        # where the primary does not exist (is a massless remnant) 
        # and the secondary is the one evolving

        # star 1 is a massless remnant, only star 2 exists
        elif binary.non_existent_companion == 1:
            primary = binary.star_1
            primary.co = True
            primary.htrack = False
            secondary = binary.star_2
            secondary.co = False

            # only H-rich star left
            if (binary.star_2.state in STAR_STATES_H_RICH):
                secondary.htrack = True
            # only He star left
            elif (binary.star_2.state in STAR_STATES_FOR_Hestar_MATCHING):
                secondary.htrack = False
            # only a compact object left
            elif (binary.star_2.state in STAR_STATES_CO):
                only_CO = True
                return primary, secondary, only_CO
            else:
                raise ValueError(f"State {binary.star_2.state} is not recognized!")

        # star 2 is a massless remnant, only star 1 exists
        elif binary.non_existent_companion == 2:
            primary = binary.star_2
            primary.co = True
            primary.htrack = False
            secondary = binary.star_1
            secondary.co = False

            if (binary.star_1.state in STAR_STATES_H_RICH):
                secondary.htrack = True
            elif (binary.star_1.state in STAR_STATES_FOR_Hestar_MATCHING):
                secondary.htrack = False
            elif (binary.star_1.state in STAR_STATES_CO):
                only_CO = True
                return primary, secondary, only_CO
            else:
                raise ValueError(f"State {binary.star_1.state} is not recognized!")
        else:
            raise POSYDONError(f"non_existent_companion = {binary.non_existent_companion} (should be -1, 0, 1, or 2).")

        return primary, secondary, only_CO

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
        star : SingleStar object
            A single star object that contains the star's properties.
        
        htrack : bool
            A boolean that specifies whether the star would be found in the 
            hydrogen rich single star grid or not (in which case it is
            matched to the helium rich single star grid).

        m0 : float
            Initial stellar mass (in solar units) of the single star track 
            that we will grab values from and update `star` with.
        
        t0 : float
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

        grid = self.grid_Hrich if htrack else self.grid_strippedHe
        get_profile = grid.get_profile
        profile_new = np.array(get_profile('mass', m0)[1])

        for i in self.profile_keys:
            profile_new[i] = get_profile(i, m0)[0]
        profile_new['omega'] = star.surf_avg_omega

        star.profile = profile_new
