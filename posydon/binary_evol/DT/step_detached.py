"""Detached evolution step."""


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
]

import os
import numpy as np
import pandas as pd
import time
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

from posydon.config import PATH_TO_POSYDON_DATA
from posydon.binary_evol.binarystar import BINARYPROPERTIES
from posydon.binary_evol.singlestar import STARPROPERTIES
from posydon.interpolation.interpolation import GRIDInterpolator
#from posydon.interpolation.data_scaling import DataScaler
from posydon.utils.common_functions import (bondi_hoyle,
                                            orbital_period_from_separation,
                                            roche_lobe_radius,
                                            check_state_of_star,
                                            convert_metallicity_to_string,
                                            set_binary_to_failed)
from posydon.utils.interpolators import PchipInterpolator2
from posydon.binary_evol.flow_chart import (STAR_STATES_CC, 
                                            STAR_STATES_CO, 
                                            STAR_STATES_H_RICH_EVOLVABLE,
                                            STAR_STATES_HE_RICH_EVOLVABLE)
import posydon.utils.constants as const
from posydon.utils.posydonerror import (NumericalError, MatchingError,
                                        POSYDONError, FlowError,
                                        ClassificationError)
from posydon.utils.posydonwarning import Pwarn

from posydon.grids.track_match import track_matcher

LIST_ACCEPTABLE_STATES_FOR_HMS = ["H-rich_Core_H_burning",
                                  "accreted_He_Core_H_burning"] # REMOVE

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

STAR_STATES_H_RICH = [
    'H-rich_Core_H_burning',
    'H-rich_Core_He_burning',
    'H-rich_Shell_H_burning',
    'H-rich_Central_He_depleted',
    'H-rich_Shell_He_burning',
    'H-rich_Core_C_burning',
    'H-rich_Central_C_depletion',
    'H-rich_non_burning',
    'accreted_He_Core_H_burning',
    'accreted_He_non_burning'
]


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

MATCHING_WITH_RELATIVE_DIFFERENCE = ["center_he4"]




class detached_step:
    """Evolve a detached binary.

    The binary will be evolved until Roche-lobe overflow, core-collapse or
    maximum simulation time, using the standard equations that govern the
    orbital evolution.

    Parameters
    ----------
    path : str
        Path to the directory that contains a HDF5 grid.
    dt : float
        The timestep size, in years, to be appended to the history of the
        binary. None means only the final step.
        Note: do not select very small timesteps cause it may mess with the
        solving of the ODE.
    n_o_steps_history: int
        Alternatively, we can define the number of timesteps to be appended to
        the history of the binary. None means only the final step. If both `dt`
        and `n_o_steps_history` are different than None, `dt` has priority.
    matching_method: str
        Method to find the best match between a star from a previous step and a
        point in a single MIST-like stellar track. Options "root" (which tries
        to find a root of two matching quantities, and it is possible to not
        achieve it) or "minimize" (minimizes the sum of squares of differences
        of various quantities between the previous step and the track).
    verbose : Boolean
        True if we want to print stuff.
    do_wind_loss: Boolean
        If True, take into account change of separation due to mass loss from
        the star.
    do_tides: Booleans
        If True, take into account change of separation, eccentricity and star
        spin due to tidal forces.
    do_gravitational_radiation: Boolean
        If True, take into account change of separation and eccentricity due to
        gravitational wave radiation.
    do_magnetic_braking: Boolean
        If True, take into account change of star spin due to magnetic braking.
    magnetic_braking_mode: String
        A string corresponding to the desired magnetic braking prescription.
            -- RVJ83: Rappaport, Verbunt, & Joss 1983
            -- M15: Matt et al. 2015
            -- G18: Garraffo et al. 2018
            -- CARB: Van & Ivanova 2019
    do_stellar_evolution_and_spin_from_winds: Boolean
        If True, take into account change of star spin due to change of its
        moment of inertia during its evolution and due to spin angular momentum
        loss due to winds.

    Attributes
    ----------
    KEYS : list of str
           Contains valid keywords which is used
           to extract quantities from the grid.
    grid : GRIDInterpolator
           Object to interpolate between the time-series in
           the h5 grid.
    initial_mass : list of float
            Contains the initial masses of the stars in the grid.

    Note
    ----
    A matching between the properties of the star, and the h5 tracks are
    required. In the "root" solver matching_method, if the root solver fails
    then the evolution will immediately end, and the binary state will be
    tagged with "Root solver failed". In the "minimize" matching_method, we
    minimize the sum of squares of differences of various quantities between
    the previous step and the h5 track.

    Warns
    -----
    UserWarning
        If the call cannot determine the primary or secondary in the binary.

    Raises
    ------
    Exception
        If the ode-solver fails to solve the differential equation
        that governs the orbital evolution.
    """

    def __init__(
            self,
            grid_name_Hrich=None,
            grid_name_strippedHe=None,
            metallicity=None,
            path=PATH_TO_POSYDON_DATA,
            dt=None,
            n_o_steps_history=None,
            matching_method="minimize",
            initial_mass=None,
            rootm=None,
            verbose=False,
            do_wind_loss=True,
            do_tides=True,
            do_gravitational_radiation=True,
            do_magnetic_braking=True,
            magnetic_braking_mode="RVJ83",
            do_stellar_evolution_and_spin_from_winds=True,
            RLO_orbit_at_orbit_with_same_am=False,
            list_for_matching_HMS=None,
            list_for_matching_postMS=None,
            list_for_matching_HeStar=None
    ):
        """Initialize the step. See class documentation for details."""
        self.metallicity = convert_metallicity_to_string(metallicity)
        self.dt = dt
        self.n_o_steps_history = n_o_steps_history
        self.matching_method = matching_method
        self.do_wind_loss = do_wind_loss
        self.do_tides = do_tides
        self.do_gravitational_radiation = do_gravitational_radiation
        self.do_magnetic_braking = do_magnetic_braking
        self.magnetic_braking_mode = magnetic_braking_mode
        self.do_stellar_evolution_and_spin_from_winds = (
            do_stellar_evolution_and_spin_from_winds
        )
        self.RLO_orbit_at_orbit_with_same_am = RLO_orbit_at_orbit_with_same_am
        self.initial_mass = initial_mass
        self.rootm = rootm
        self.verbose = verbose
        self.list_for_matching_HMS = list_for_matching_HMS
        self.list_for_matching_postMS = list_for_matching_postMS
        self.list_for_matching_HeStar = list_for_matching_HeStar

        # mapping a combination of (key, htrack, method) to a pre-trained
        # DataScaler instance, created the first time it is requested
        self.stored_scalers = {}

        if self.verbose:
            print(
                dt,
                n_o_steps_history,
                matching_method,
                do_wind_loss,
                do_tides,
                do_gravitational_radiation,
                do_magnetic_braking,
                magnetic_braking_mode,
                do_stellar_evolution_and_spin_from_winds)

        self.translate = DEFAULT_TRANSLATION

        # these are the KEYS read from POSYDON h5 grid files (after translating
        # them to the appropriate columns)
        self.KEYS = DEFAULT_TRANSLATED_KEYS
        self.KEYS_POSITIVE = (
            'mass_conv_reg_fortides',
            'thickness_conv_reg_fortides',
            'radius_conv_reg_fortides'
        )

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

        if grid_name_Hrich is None:
            grid_name_Hrich = os.path.join('single_HMS', self.metallicity+'_Zsun.h5')
        self.grid_Hrich = GRIDInterpolator(os.path.join(path, grid_name_Hrich))

        if grid_name_strippedHe is None:
            grid_name_strippedHe = os.path.join('single_HeMS', self.metallicity+'_Zsun.h5')
        self.grid_strippedHe = GRIDInterpolator(os.path.join(path, grid_name_strippedHe))

    def __repr__(self):
        """Return the type of evolution type."""
        return "Detached Step."

    def __call__(self, binary):
        """Evolve the binary until RLO or compact object formation."""
        
        @event(True, 1)
        def ev_rlo1(t, y):
            """Difference between radius and Roche lobe at a given time.

            Used to check if there is RLOF mass transfer during the detached
            binary evolution interpolation.

            Parameters
            ----------
            t : float
                Time of the evolution, in years.
            y : tuple of floats
                [separation, eccentricity] at that time. Separation should be
                in solar radii.

            Returns
            -------
            float
                Difference between stellar radius and Roche lobe radius in
                solar radii.

            """
            pri_mass = interp1d_pri["mass"](t - t_offset_pri)
            sec_mass = interp1d_sec["mass"](t - t_offset_sec)

            sep = y[0]
            ecc = y[1]
           
            RL = roche_lobe_radius(sec_mass, pri_mass, (1 - ecc) * sep)
            
            # 95% filling of the RL is enough to assume beginning of RLO,
            # as we do in CO-HMS_RLO grid
            return interp1d_sec["R"](t - t_offset_sec) - 0.95*RL

        @event(True, 1)
        def ev_rlo2(t, y):
            """Difference between radius and Roche lobe at a given time.

            Used to check if there is RLOF mass transfer during the detached
            binary evolution interpolation.

            Parameters
            ----------
            t : float
                Time of the evolution, in years
            y : tuple of floats
                [separation, eccentricity] at that time. Separation should be
                in solar radii.

            Returns
            -------
            float
                Difference between stellar radius and Roche lobe radius in
                solar radii.

            """
            pri_mass = interp1d_pri["mass"](t - t_offset_pri)
            sec_mass = interp1d_sec["mass"](t - t_offset_sec)
            
            sep = y[0]
            ecc = y[1]
           
            RL = roche_lobe_radius(pri_mass, sec_mass, (1 - ecc) * sep)
            
            return interp1d_pri["R"](t - t_offset_pri) - 0.95*RL

        @event(True, 1)
        def ev_rel_rlo1(t, y):
            """Relative difference between radius and Roche lobe.

            Used to check if there is RLOF mass transfer during the detached
            binary evolution interpolation.

            Parameters
            ----------
            t : float
                Time of the evolution, in years.
            y : tuple of floats
                [separation, eccentricity] at that time. Separation should be
                in solar radii.

            Returns
            -------
            float
                Relative difference between stellar radius and Roche lobe
                radius.

            """
            pri_mass = interp1d_pri["mass"](t - t_offset_pri)
            sec_mass = interp1d_sec["mass"](t - t_offset_sec)
            
            sep = y[0]
            ecc = y[1]
           
            RL = roche_lobe_radius(sec_mass, pri_mass, (1 - ecc) * sep)

            return (interp1d_sec["R"](t - t_offset_sec) - RL) / RL

        @event(True, 1)
        def ev_rel_rlo2(t, y):
            """Relative difference between radius and Roche lobe.

            Used to check if there is RLOF mass transfer during the detached
            binary evolution interpolation.

            Parameters
            ----------
            t : float
                Time of the evolution, in years.
            y : tuple of floats
                [separation, eccentricity] at that time. Separation should be
                in solar radii.

            Returns
            -------
            float
                Relative difference between stellar radius and Roche lobe
                radius.

            """
            pri_mass = interp1d_pri["mass"](t - t_offset_pri)
            sec_mass = interp1d_sec["mass"](t - t_offset_sec)

            sep = y[0]
            ecc = y[1]
           
            RL = roche_lobe_radius(pri_mass, sec_mass, (1 - ecc) * sep)

            return (interp1d_pri["R"](t - t_offset_pri) - RL) / RL

        @event(True, -1)
        def ev_max_time1(t, y):
            return t_max_sec + t_offset_sec - t

        @event(True, -1)
        def ev_max_time2(t, y):
            return t_max_pri + t_offset_pri - t
        
        def get_omega(star, is_secondary=True):
            """Calculate the spin of a star.

            Parameters
            ----------
            star : SingleStar object

            Returns
            -------
            float
                spin of the star in radians per year
            """

            if (star.log_total_angular_momentum is not None
                    and star.total_moment_of_inertia is not None
                    and not np.isnan(star.log_total_angular_momentum)
                    and not np.isnan(star.total_moment_of_inertia)):
                
                # the last factor converts rad/s to rad/yr
                omega_in_rad_per_year = (
                        10.0 ** star.log_total_angular_momentum
                        / star.total_moment_of_inertia * const.secyer)  
                
                if self.verbose:
                    print("calculating initial omega using angular momentum and moment of inertia")
            else:
                # we equate the secondary's initial omega to surf_avg_omega
                # (although the critical rotation should be improved to
                # take into account radiation pressure)

                if (star.surf_avg_omega is not None and not np.isnan(star.surf_avg_omega)):

                    if self.verbose:
                        print("calculating initial omega using surf_avg_omega")

                    omega_in_rad_per_year = star.surf_avg_omega * const.secyer                    
                        
                elif (star.surf_avg_omega_div_omega_crit is not None
                        and not np.isnan(star.surf_avg_omega_div_omega_crit)):
                    
                    if (star.log_R is not None and not np.isnan(star.log_R)):
                        omega_in_rad_per_year = (
                            star.surf_avg_omega_div_omega_crit * np.sqrt(
                                const.standard_cgrav * star.mass * const.msol
                                / ((10.0 ** (star.log_R) * const.rsol) ** 3)) * const.secyer)
                        
                    else:
                        if is_secondary:
                            radius_to_be_used = interp1d_sec["R"](interp1d_sec["t0"])
                            mass_to_be_used = interp1d_sec["mass"](interp1d_sec["t0"])
                        else:
                            radius_to_be_used = interp1d_pri["R"](interp1d_pri["t0"])
                            mass_to_be_used = interp1d_pri["mass"](interp1d_pri["t0"])
                        
                        if self.verbose:
                            print("calculating initial omega using surf_avg_omega_div_omega_crit")

                        omega_in_rad_per_year = (
                            star.surf_avg_omega_div_omega_crit * np.sqrt(
                                const.standard_cgrav * mass_to_be_used * const.msol
                                / ((radius_to_be_used * const.rsol) ** 3)) * const.secyer)
                   
                else:
                    omega_in_rad_per_year = 0.0
                    if self.verbose:
                        print("could not calculate initial omega, setting to zero")
                        
            if self.verbose:
                print("calculated omega_in_rad_per_year: ", omega_in_rad_per_year)

            return omega_in_rad_per_year
            

        KEYS = self.KEYS
        KEYS_POSITIVE = self.KEYS_POSITIVE

        # creating a track matching object
        self.track_matcher = track_matcher(KEYS, KEYS_POSITIVE)

        binary_sim_prop = getattr(binary, "properties")   ## simulation properties of the binary
        all_step_names = getattr(binary_sim_prop, "all_step_names")

        companion_1_exists = (binary.star_1 is not None
                              and binary.star_1.state != "massless_remnant")
        companion_2_exists = (binary.star_2 is not None
                              and binary.star_2.state != "massless_remnant")

        if companion_1_exists:
            if companion_2_exists:                  # to evolve a binary star
                self.non_existent_companion = 0
            else:                                   # star1 is a single star
                self.non_existent_companion = 2
        else:
            if companion_2_exists:                  # star2 is a single star
                self.non_existent_companion = 1
            else:                                   # no star in the system
                raise POSYDONError("There is no star to evolve. Who summoned me?")

        if self.non_existent_companion == 0: #no isolated evolution, detached step of an actual binary
            # the primary in a real binary is potential compact object, or the more evolved star
            if (binary.star_1.state in STAR_STATES_CO
                    and binary.star_2.state in STAR_STATES_H_RICH):
                primary = binary.star_1
                secondary = binary.star_2
                secondary.htrack = True
                primary.htrack = secondary.htrack
                primary.co = True

            elif (binary.star_1.state in STAR_STATES_CO
                    and binary.star_2.state in LIST_ACCEPTABLE_STATES_FOR_HeStar):
                primary = binary.star_1
                secondary = binary.star_2
                secondary.htrack = False
                primary.htrack = secondary.htrack
                primary.co = True

            elif (binary.star_2.state in STAR_STATES_CO
                    and binary.star_1.state in STAR_STATES_H_RICH):
                primary = binary.star_2
                secondary = binary.star_1
                secondary.htrack = True
                primary.htrack = secondary.htrack
                primary.co = True

            elif (binary.star_2.state in STAR_STATES_CO
                    and binary.star_1.state in LIST_ACCEPTABLE_STATES_FOR_HeStar):
                primary = binary.star_2
                secondary = binary.star_1
                secondary.htrack = False
                primary.htrack = secondary.htrack
                primary.co = True
            elif (binary.star_1.state in STAR_STATES_H_RICH
                    and binary.star_2.state in STAR_STATES_H_RICH):
                primary = binary.star_1
                secondary = binary.star_2
                secondary.htrack = True
                primary.htrack = True
                primary.co = False
            elif (binary.star_1.state in LIST_ACCEPTABLE_STATES_FOR_HeStar
                    and binary.star_2.state in STAR_STATES_H_RICH):
                primary = binary.star_1
                secondary = binary.star_2
                secondary.htrack = True
                primary.htrack = False
                primary.co = False
            elif (binary.star_2.state in LIST_ACCEPTABLE_STATES_FOR_HeStar
                    and binary.star_1.state in STAR_STATES_H_RICH):
                primary = binary.star_2
                secondary = binary.star_1
                secondary.htrack = True
                primary.htrack = False
                primary.co = False
            elif (binary.star_1.state in LIST_ACCEPTABLE_STATES_FOR_HeStar
                    and binary.star_2.state
                    in LIST_ACCEPTABLE_STATES_FOR_HeStar):
                primary = binary.star_1
                secondary = binary.star_2
                secondary.htrack = False
                primary.htrack = False
                primary.co = False
            else:
                raise ValueError("States not recognized!")

        # star 1 is a massless remnant, only star 2 exists
        elif self.non_existent_companion == 1:
            # we force primary.co=True for all isolated evolution,
            # where the secondary is the one evolving one
            primary = binary.star_1
            primary.co = True
            primary.htrack = False
            secondary = binary.star_2
            if (binary.star_2.state in STAR_STATES_H_RICH):
                secondary.htrack = True
            elif (binary.star_2.state in LIST_ACCEPTABLE_STATES_FOR_HeStar):
                secondary.htrack = False
            elif (binary.star_2.state in STAR_STATES_CO):
                # only a compact object left
                return
            else:
                raise ValueError("State not recognized!")

        # star 2 is a massless remnant, only star 1 exists
        elif self.non_existent_companion == 2:
            primary = binary.star_2
            primary.co = True
            primary.htrack = False
            secondary = binary.star_1
            if (binary.star_1.state in STAR_STATES_H_RICH):
                secondary.htrack = True
            elif (binary.star_1.state in LIST_ACCEPTABLE_STATES_FOR_HeStar):
                secondary.htrack = False
            elif (binary.star_1.state in STAR_STATES_CO):
                return
            else:
                raise ValueError("State not recognized!")
        else:
            raise POSYDONError("Non existent companion has not a recognized value!")

        # get the matched data of two stars, respectively
        interp1d_sec, m0, t0 = track_matcher.get_star_data(binary, secondary, primary, secondary.htrack, co=False)

        primary_not_normal = (primary.co) or (self.non_existent_companion in [1,2])
        primary_normal = (not primary.co) and self.non_existent_companion == 0

        if primary_not_normal:
            # copy the secondary star except mass which is of the primary,
            # and radius, mdot, Idot = 0
            interp1d_pri = track_matcher.get_star_data(
                binary, secondary, primary, secondary.htrack, co=True,
                copy_prev_m0=m0, copy_prev_t0=t0)[0]
        elif primary_normal:
            interp1d_pri = track_matcher.get_star_data(
                binary, primary, secondary, primary.htrack, False)[0]
        else:
            raise ValueError("During matching, the primary should either be normal (stellar object) or ",
                             "not normal (CO, nonexistent companion).")


        if interp1d_sec is None or interp1d_pri is None:
            failed_state = binary.state
            set_binary_to_failed(binary)
            raise MatchingError(f"Grid matching failed for {failed_state} binary.")    
                  
                   
        t0_sec = interp1d_sec["t0"]
        t0_pri = interp1d_pri["t0"]
        m01 = interp1d_sec["m0"]
        m02 = interp1d_pri["m0"]
        t_max_sec = interp1d_sec["t_max"]
        t_max_pri = interp1d_pri["t_max"]
        t_offset_sec = binary.time - t0_sec
        t_offset_pri = binary.time - t0_pri
        max_time = interp1d_sec["max_time"]


        if (ev_rlo1(binary.time, [binary.separation, binary.eccentricity]) >= 0
                or ev_rlo2(binary.time, [binary.separation, binary.eccentricity]) >= 0):
            binary.state = "initial_RLOF"
            return
            
        else:
            if not (max_time - binary.time > 0.0):
                raise ValueError("max_time is lower than the current time. "
                                "Evolution of the detached binary will go to "
                                "lower times.")
            with np.errstate(all="ignore"):
                omega_in_rad_per_year_sec = get_omega(secondary)
                if primary_not_normal:
                    # omega of compact objects or masslessremnant won't be used for intergration
                    omega_in_rad_per_year_pri = omega_in_rad_per_year_sec
                elif not primary.co:
                    omega_in_rad_per_year_pri = get_omega(primary,is_secondary = False)

                t_before_ODEsolution = time.time()
                try:
                    s = solve_ivp(
                        lambda t, y: diffeq(
                            t,
                            y,
                            interp1d_sec["R"](t - t_offset_sec),
                            interp1d_sec["L"](t - t_offset_sec),
                            *[
                                interp1d_sec[key](t - t_offset_sec)
                                for key in KEYS[1:11]
                            ],
                            interp1d_sec["Idot"](t - t_offset_sec),
                            interp1d_sec["conv_env_turnover_time_l_b"](
                                t - t_offset_sec),
                            interp1d_pri["R"](t - t_offset_pri),
                            interp1d_pri["L"](t - t_offset_pri),
                            *[
                                interp1d_pri[key](t - t_offset_pri)
                                for key in KEYS[1:11]
                            ],
                            interp1d_pri["Idot"](t - t_offset_pri),
                            interp1d_pri["conv_env_turnover_time_l_b"](
                                t - t_offset_pri),
                            self.do_wind_loss,
                            self.do_tides,
                            self.do_gravitational_radiation,
                            self.do_magnetic_braking,
                            self.magnetic_braking_mode,
                            self.do_stellar_evolution_and_spin_from_winds
                            # ,self.verbose
                        ),
                        events=[ev_rlo1, ev_rlo2, ev_max_time1, ev_max_time2],
                        method="Radau",
                        t_span=(binary.time, max_time),
                        y0=[
                            binary.separation,
                            binary.eccentricity,
                            omega_in_rad_per_year_sec,
                            omega_in_rad_per_year_pri,
                        ],
                        dense_output=True,
                        # vectorized=True
                    )
                except Exception:
                    s = solve_ivp(
                        lambda t, y: diffeq(
                            t,
                            y,
                            interp1d_sec["R"](t - t_offset_sec),
                            interp1d_sec["L"](t - t_offset_sec),
                            *[
                                interp1d_sec[key](t - t_offset_sec)
                                for key in KEYS[1:11]
                            ],
                            interp1d_sec["Idot"](t - t_offset_sec),
                            interp1d_sec["conv_env_turnover_time_l_b"](
                                t - t_offset_sec),
                            interp1d_pri["R"](t - t_offset_pri),
                            interp1d_pri["L"](t - t_offset_pri),
                            *[
                                interp1d_pri[key](t - t_offset_pri)
                                for key in KEYS[1:11]
                            ],
                            interp1d_pri["Idot"](t - t_offset_pri),
                            interp1d_pri["conv_env_turnover_time_l_b"](
                                t - t_offset_pri),
                            self.do_wind_loss,
                            self.do_tides,
                            self.do_gravitational_radiation,
                            self.do_magnetic_braking,
                            self.magnetic_braking_mode,
                            self.do_stellar_evolution_and_spin_from_winds
                            # ,self.verbose
                        ),
                        events=[ev_rlo1, ev_rlo2, ev_max_time1, ev_max_time2],
                        method="RK45",
                        t_span=(binary.time, max_time),
                        y0=[
                            binary.separation,
                            binary.eccentricity,
                            omega_in_rad_per_year_sec,
                            omega_in_rad_per_year_pri,
                        ],
                        dense_output=True,
                        # vectorized=True
                    )

            t_after_ODEsolution = time.time()

            if self.verbose:
                print(f"\nODE solver duration: {t_after_ODEsolution-t_before_ODEsolution:.6g}")
                print("solution of ODE", s)

            if s.status == -1:
                failed_state = binary.state
                set_binary_to_failed(binary)
                raise NumericalError(f"Integration failed for {failed_state} binary.")               

            if self.dt is not None and self.dt > 0:
                t = np.arange(binary.time, s.t[-1] + self.dt/2.0, self.dt)[1:]
                if t[-1] < s.t[-1]:
                    t = np.hstack([t, s.t[-1]])
            elif (self.n_o_steps_history is not None
                    and self.n_o_steps_history > 0):
                t_step = (s.t[-1] - binary.time) / self.n_o_steps_history
                t = np.arange(binary.time, s.t[-1] + t_step / 2.0, t_step)[1:]
                if t[-1] < s.t[-1]:
                    t = np.hstack([t, s.t[-1]])
            else:  # self.dt is None and self.n_o_steps_history is None
                t = np.array([s.t[-1]])

            sep_interp, ecc_interp, omega_interp_sec, omega_interp_pri = s.sol(t)
            mass_interp_sec = interp1d_sec[self.translate["mass"]]
            mass_interp_pri = interp1d_pri[self.translate["mass"]]

            interp1d_sec["sep"] = sep_interp
            interp1d_sec["ecc"] = ecc_interp    
            interp1d_sec["omega"] = omega_interp_sec
            interp1d_pri["omega"] = omega_interp_pri

            interp1d_sec["porb"] = orbital_period_from_separation(
                sep_interp, mass_interp_sec(t - t_offset_sec),
                mass_interp_pri(t - t_offset_pri))
            interp1d_pri["porb"] = orbital_period_from_separation(
                sep_interp, mass_interp_pri(t - t_offset_pri),
                mass_interp_sec(t - t_offset_sec))

            interp1d_sec["time"] = t

            ## UPDATE STAR AND BINARY PROPERTIES WITH INTERPOLATED VALUES
            for obj, prop in zip([secondary, primary, binary], 
                                 [STARPROPERTIES, STARPROPERTIES, BINARYPROPERTIES]):
                
                for key in prop:
                    if key in ["event",
                               "mass_transfer_case",
                               "nearest_neighbour_distance",
                               "state", "metallicity", "V_sys"]:
                        current = getattr(obj, key)
                        # For star objects, the state is calculated further below
                        history = [current] * len(t[:-1])

                    elif (key in ["surf_avg_omega_div_omega_crit"] and obj == secondary):
                        # replace the actual surf_avg_w with the effective omega,
                        # which takes into account the whole star
                        # key  = 'effective_omega' # in rad/sec
                        # current = s.y[2][-1] / 3.1558149984e7
                        # history_of_attribute = s.y[2][:-1] / 3.1558149984e7
                        omega_crit_current_sec = np.sqrt(const.standard_cgrav
                            * interp1d_sec[self.translate["mass"]](t[-1] - t_offset_sec).item() * const.msol
                            / (interp1d_sec[self.translate["R"]](t[-1] - t_offset_sec).item() * const.rsol)**3)
                        
                        omega_crit_hist_sec = np.sqrt(const.standard_cgrav
                            * interp1d_sec[self.translate["mass"]](t[:-1] - t_offset_sec) * const.msol
                            / (interp1d_sec[self.translate["R"]](t[:-1] - t_offset_sec) * const.rsol)**3)

                        current = (interp1d_sec["omega"][-1] / const.secyer / omega_crit_current_sec)
                        history = (interp1d_sec["omega"][:-1] / const.secyer / omega_crit_hist_sec)
                        
                    elif (key in ["surf_avg_omega_div_omega_crit"] and obj == primary):
                        if primary.co:
                            current = None
                            history = [current] * len(t[:-1])
                        elif not primary.co:
                            # TODO: change `item()` to 0
                            omega_crit_current_pri = np.sqrt(const.standard_cgrav
                                * interp1d_pri[self.translate["mass"]](t[-1] - t_offset_pri).item() * const.msol
                                / (interp1d_pri[self.translate["R"]](t[-1] - t_offset_pri).item() * const.rsol)**3)

                            omega_crit_hist_pri = np.sqrt(const.standard_cgrav
                                * interp1d_pri[self.translate["mass"]](t[:-1] - t_offset_pri) * const.msol
                                / (interp1d_pri[self.translate["R"]](t[:-1] - t_offset_pri) * const.rsol)**3)

                            current = (interp1d_pri["omega"][-1] / const.secyer / omega_crit_current_pri)
                            history = (interp1d_pri["omega"][:-1] / const.secyer / omega_crit_hist_pri)
                            
                    elif (key in ["surf_avg_omega"] and obj == secondary):
                        current = interp1d_sec["omega"][-1] / const.secyer
                        history = interp1d_sec["omega"][:-1] / const.secyer

                    elif (key in ["surf_avg_omega"] and obj == primary):
                        if primary.co:
                            current = None
                            history = [current] * len(t[:-1])
                        else:
                            current = interp1d_pri["omega"][-1] / const.secyer
                            history = interp1d_pri["omega"][:-1] / const.secyer

                    elif (key in ["rl_relative_overflow_1"] and obj == binary):
                        if binary.star_1.state in ("BH", "NS", "WD","massless_remnant"):
                            current = None
                            history = [current] * len(t[:-1])

                        elif secondary == binary.star_1:
                            current = ev_rel_rlo1(t[-1], [interp1d_sec["sep"][-1], interp1d_sec["ecc"][-1]])
                            history = ev_rel_rlo1(t[:-1], [interp1d_sec["sep"][:-1], interp1d_sec["ecc"][:-1]])

                        elif secondary == binary.star_2:
                            current = ev_rel_rlo2(t[-1], [interp1d_sec["sep"][-1], interp1d_sec["ecc"][-1]])
                            history = ev_rel_rlo2(t[:-1], [interp1d_sec["sep"][:-1], interp1d_sec["ecc"][:-1]])
                            
                    elif (key in ["rl_relative_overflow_2"] and obj == binary):
                        if binary.star_2.state in ("BH", "NS", "WD","massless_remnant"):
                            current = None
                            history = [current] * len(t[:-1])

                        elif secondary == binary.star_2:
                            current = ev_rel_rlo1(t[-1], [interp1d_sec["sep"][-1], interp1d_sec["ecc"][-1]])
                            history = ev_rel_rlo1(t[:-1], [interp1d_sec["sep"][:-1], interp1d_sec["ecc"][:-1]])

                        elif secondary == binary.star_1:
                            current = ev_rel_rlo2(t[-1], [interp1d_sec["sep"][-1], interp1d_sec["ecc"][-1]])
                            history = ev_rel_rlo2(t[:-1], [interp1d_sec["sep"][:-1], interp1d_sec["ecc"][:-1]])
                            
                    elif key in ["separation", "orbital_period", "eccentricity", "time"]:
                        current = interp1d_sec[self.translate[key]][-1].item()
                        history = interp1d_sec[self.translate[key]][:-1]

                    elif (key in ["total_moment_of_inertia"] and obj == secondary):
                        current = interp1d_sec[self.translate[key]](
                            t[-1] - t_offset_sec).item() * (const.msol * const.rsol**2)
                        history = interp1d_sec[self.translate[key]](
                            t[:-1] - t_offset_sec) * (const.msol * const.rsol**2)
                        
                    elif (key in ["total_moment_of_inertia"] and obj == primary):
                        if primary.co:
                            current = getattr(obj, key)
                            history = [current] * len(t[:-1])
                        else:
                            current = interp1d_pri[self.translate[key]](
                                t[-1] - t_offset_pri).item() * (const.msol * const.rsol**2)
                            history = interp1d_pri[self.translate[key]](
                                t[:-1] - t_offset_pri) * (const.msol * const.rsol**2)
                            
                    elif (key in ["log_total_angular_momentum"] and obj == secondary):

                        current_omega = interp1d_sec["omega"][-1]

                        ## add a warning catch if the current omega has an invalid value
                        ## (otherwise python will throw an insuppressible warning when taking the log)
                        if interp1d_sec["omega"][-1] <=0:
                            Pwarn("Trying to compute log angular momentum for object with no spin", 
                                  "InappropriateValueWarning")
                            current_omega = np.nan
                                               
                        current = np.log10(
                            (current_omega / const.secyer)
                                * (interp1d_sec[
                                    self.translate["total_moment_of_inertia"]](t[-1] - t_offset_sec).item() * 
                                    (const.msol * const.rsol ** 2)))                                             
                        history = np.log10(
                            (interp1d_sec["omega"][:-1] / const.secyer)
                            * (interp1d_sec[self.translate["total_moment_of_inertia"]](
                                    t[:-1] - t_offset_sec) * (const.msol * const.rsol**2)))
                        
                    elif (key in ["log_total_angular_momentum"] and obj == primary):
                        if primary.co:
                            current = getattr(obj, key)
                            history = [current] * len(t[:-1])
                        else:
                            current = np.log10(
                                (interp1d_pri["omega"][-1] / const.secyer)
                                * (interp1d_pri[self.translate["total_moment_of_inertia"]](
                                        t[-1] - t_offset_pri).item() * (const.msol * const.rsol**2)))
                            history = np.log10(
                                (interp1d_pri["omega"][:-1] / const.secyer)
                                * (interp1d_pri[self.translate["total_moment_of_inertia"]](
                                        t[:-1] - t_offset_pri) * (const.msol * const.rsol**2)))
                            
                    elif (key in ["spin"] and obj == secondary):
                        current = (const.clight
                            * (interp1d_sec["omega"][-1] / const.secyer)
                            * interp1d_sec[self.translate["total_moment_of_inertia"]](
                                t[-1] - t_offset_sec).item() * (const.msol * const.rsol**2)
                            / (const.standard_cgrav * (interp1d_sec[self.translate["mass"]](
                                    t[-1] - t_offset_sec).item() * const.msol)**2))
                        history = (const.clight
                            * (interp1d_sec["omega"][:-1] / const.secyer)
                            * interp1d_sec[self.translate["total_moment_of_inertia"]](
                                    t[:-1] - t_offset_sec)* (const.msol * const.rsol**2)
                            / (const.standard_cgrav * (interp1d_sec[self.translate["mass"]](
                                    t[:-1] - t_offset_sec) * const.msol)**2))
                        
                    elif (key in ["spin"] and obj == primary):
                        if primary.co:
                            current = getattr(obj, key)
                            history = [current] * len(t[:-1])
                        else:
                            current = (const.clight
                                * (interp1d_pri["omega"][-1] / const.secyer)
                                * interp1d_pri[self.translate["total_moment_of_inertia"]](
                                        t[-1] - t_offset_pri).item() * (const.msol * const.rsol**2)
                                / (const.standard_cgrav * (interp1d_pri[self.translate["mass"]](
                                        t[-1] - t_offset_pri).item() * const.msol)**2))
                            history = (const.clight 
                                * (interp1d_pri["omega"][:-1] / const.secyer)
                                * interp1d_pri[self.translate["total_moment_of_inertia"]](
                                        t[:-1] - t_offset_pri) * (const.msol * const.rsol**2)
                                / (const.standard_cgrav * (interp1d_pri[self.translate["mass"]](
                                        t[:-1] - t_offset_pri) * const.msol)**2))
                            
                    elif (key in ["lg_mdot", "lg_wind_mdot"] and obj == secondary):
                        # in detached step, lg_mdot = lg_wind_mdot
                        if interp1d_sec[self.translate[key]](t[-1] - t_offset_sec) == 0:
                            current = -98.99
                        else:
                            current = np.log10(np.abs(interp1d_sec[self.translate[key]](
                                t[-1] - t_offset_sec))).item()
                            
                        history = np.ones_like(t[:-1])

                        for i in range(len(t)-1):
                            if interp1d_sec[self.translate[key]](t[i] - t_offset_sec) == 0:
                                history[i] = -98.99
                            else:
                                history[i] = np.log10(np.abs(interp1d_sec[self.translate[key]](
                                        t[i] - t_offset_sec)))

                    elif (key in ["lg_mdot", "lg_wind_mdot"] and obj == primary):
                        if primary.co:
                            current = None
                            history = [current] * len(t[:-1])
                        else:
                            if interp1d_sec[self.translate[key]](t[-1] - t_offset_sec) == 0:
                                current = -98.99
                            else:
                                current = np.log10(np.abs(interp1d_sec[self.translate[key]](
                                        t[-1] - t_offset_sec))).item()
                            history = np.ones_like(t[:-1])
                            for i in range(len(t)-1):
                                if (interp1d_sec[self.translate[key]](t[i] - t_offset_sec) == 0):
                                    history[i] = -98.99
                                else:
                                    history[i] = np.log10(np.abs(interp1d_sec[self.translate[key]](
                                            t[i] - t_offset_sec)))
                                    
                    elif (self.translate[key] in interp1d_sec and obj == secondary):
                        current = interp1d_sec[self.translate[key]](t[-1] - t_offset_sec).item()
                        history = interp1d_sec[self.translate[key]](t[:-1] - t_offset_sec)
                        
                    elif (self.translate[key] in interp1d_pri and obj == primary):
                        if primary.co:
                            current = getattr(obj, key)
                            history = [current] * len(t[:-1])
                        else:
                            current = interp1d_pri[self.translate[key]](t[-1] - t_offset_pri).item()
                            history = interp1d_pri[self.translate[key]](t[:-1] - t_offset_pri)
                            
                    elif key in ["profile"]:
                        current = None
                        history = [current] * len(t[:-1])

                    else:
                        current = np.nan
                        history = np.ones_like(t[:-1]) * current

                    setattr(obj, key, current)
                    getattr(obj, key + "_history").extend(history)

            secondary.state = check_state_of_star(secondary, star_CO=False)

            for timestep in range(-len(t[:-1]), 0):
                secondary.state_history[timestep] = check_state_of_star(secondary, i=timestep, star_CO=False)


            if primary.state == "massless_remnant":
                pass

            elif primary.co:
                mdot_acc = np.atleast_1d(bondi_hoyle(
                    binary, primary, secondary, slice(-len(t), None),
                    wind_disk_criteria=True, scheme='Kudritzki+2000'))
                primary.lg_mdot = np.log10(mdot_acc.item(-1))
                primary.lg_mdot_history[len(primary.lg_mdot_history) - len(t) + 1:] = np.log10(mdot_acc[:-1])
            else:
                primary.state = check_state_of_star(primary, star_CO=False)
                for timestep in range(-len(t[:-1]), 0):

                    primary.state_history[timestep] = check_state_of_star(primary, i=timestep, star_CO=False)

            ## CHECK IF THE BINARY IS IN RLO
            if s.t_events[0] or s.t_events[1]:  

                if self.RLO_orbit_at_orbit_with_same_am:
                    # final circular orbit conserves angular momentum
                    # compared to the eccentric orbit
                    binary.separation *= (1 - s.y[1][-1]**2)
                    binary.orbital_period *= (1 - s.y[1][-1]**2) ** 1.5
                else:
                    # final circular orbit is at periastron of the ecc. orbit
                    binary.separation *= (1 - s.y[1][-1])
                    binary.orbital_period *= (1 - s.y[1][-1]) ** 1.5

                assert np.abs(
                    binary.orbital_period - orbital_period_from_separation(
                        binary.separation, secondary.mass, primary.mass)
                ) / binary.orbital_period < 10 ** (-2)
                binary.eccentricity = 0

                if s.t_events[0]:
                    if secondary == binary.star_1:
                        binary.state = "RLO1"
                        binary.event = "oRLO1"
                    else:
                        binary.state = "RLO2"
                        binary.event = "oRLO2"

                elif s.t_events[1]:
                    if secondary == binary.star_1:
                        binary.state = "RLO2"
                        binary.event = "oRLO2"
                    else:
                        binary.state = "RLO1"
                        binary.event = "oRLO1"
                
                if ('step_HMS_HMS_RLO' not in all_step_names):
                    if ((binary.star_1.state in STAR_STATES_HE_RICH_EVOLVABLE 
                         and binary.star_2.state in STAR_STATES_H_RICH_EVOLVABLE)
                    or (binary.star_1.state in STAR_STATES_H_RICH_EVOLVABLE
                         and binary.star_2.state in STAR_STATES_HE_RICH_EVOLVABLE)):
                        set_binary_to_failed(binary)
                        raise FlowError("Evolution of H-rich/He-rich stars in RLO onto H-rich/He-rich stars after " 
                                    "HMS-HMS not yet supported.") 
                
                    elif (binary.star_1.state in STAR_STATES_H_RICH_EVOLVABLE
                         and binary.star_2.state in STAR_STATES_H_RICH_EVOLVABLE):
                        set_binary_to_failed(binary)
                        raise ClassificationError("Binary is in the detached step but has stable RLO with two HMS stars - "
                                              "should it have undergone CE (was its HMS-HMS interpolation class unstable MT?)") 
                    

            ## CHECK IF STARS WILL UNDERGO CC
            elif s.t_events[2]:
                # reached t_max of track. End of life (possible collapse) of secondary
                if secondary == binary.star_1:
                    binary.event = "CC1"
                else:
                    binary.event = "CC2"

                track_matcher.get_star_final_values(secondary, secondary.htrack, m01)
                track_matcher.get_star_profile(secondary, secondary.htrack, m01)

                if not primary.co and primary.state in STAR_STATES_CC:
                    # simultaneous core-collapse of the other star as well
                    primary_time = t_max_pri + t_offset_pri - t[-1]
                    secondary_time = t_max_sec + t_offset_sec - t[-1]

                    if primary_time == secondary_time:
                        # we manually check if s.t_events[3] should also be happening simultaneously
                        track_matcher.get_star_final_values(primary, primary.htrack, m02)
                        track_matcher.get_star_profile(primary, primary.htrack, m02)

                    if primary.mass != secondary.mass:
                        raise POSYDONError(
                            "Both stars are found to be ready for collapse "
                            "(i.e. end of their life) during the detached "
                            "step, but do not have the same mass")
                    
            elif s.t_events[3]:
                # reached t_max of track. End of life (possible collapse) of primary
                if secondary == binary.star_1:
                    binary.event = "CC2"
                else:
                    binary.event = "CC1"

                track_matcher.get_star_final_values(primary, primary.htrack, m02)
                track_matcher.get_star_profile(primary, primary.htrack, m02)
                
            else:  # Reached max_time asked.
                if binary.properties.max_simulation_time - binary.time < 0.0:
                    binary.event = "MaxTime_exceeded"
                else:
                    binary.event = "maxtime"


def event(terminal, direction=0):
    """Return a helper function to set attributes for solve_ivp events."""
    def dec(f):
        f.terminal = True
        f.direction = direction
        return f
    return dec


def diffeq(
        t,
        y,
        R_sec,
        L_sec,
        M_sec,
        Mdot_sec,
        I_sec,
        # he_core_mass,
        # mass_conv_core,
        # conv_mx1_top,
        # conv_mx1_bot,
        conv_mx1_top_r_sec,
        conv_mx1_bot_r_sec,
        surface_h1_sec,
        center_h1_sec,
        M_env_sec,
        DR_env_sec,
        Renv_middle_sec,
        Idot_sec,
        tau_conv_sec,
        R_pri,
        L_pri,
        M_pri,
        Mdot_pri,
        I_pri,
        # he_core_mass,
        # mass_conv_core,
        # conv_mx1_top,
        # conv_mx1_bot,
        conv_mx1_top_r_pri,
        conv_mx1_bot_r_pri,
        surface_h1_pri,
        center_h1_pri,
        M_env_pri,
        DR_env_pri,
        Renv_middle_pri,
        Idot_pri,
        tau_conv_pri,
        do_wind_loss=True,
        do_tides=True,
        do_gravitational_radiation=True,
        do_magnetic_braking=True,
        magnetic_braking_mode="RVJ83",
        do_stellar_evolution_and_spin_from_winds=True,
        verbose=False,
):
    """Diff. equation describing the orbital evolution of a detached binary.

    The equation handles wind mass-loss [1]_, tidal [2]_, gravational [3]_
    effects and magnetic braking [4]_, [5]_, [6]_, [7]_, [8]_. It also handles
    the change of the secondary's stellar spin due to its change of moment of
    intertia and due to mass-loss from its spinning surface. It is assumed that
    the mass loss is fully non-conservative. Magnetic braking is fully applied
    to secondary stars with mass less than 1.3 Msun and fully off for stars
    with mass larger then 1.5 Msun. The effect of magnetic braking falls
    linearly for stars with mass between 1.3 Msun and 1.5 Msun.

    TODO: exaplin new features (e.g., double COs)

    Parameters
    ----------
    t : float
        The age of the system in years
    y : list of float
        Contains the separation, eccentricity and angular velocity, in Rsolar,
        dimensionless and rad/year units, respectively.
    M_pri : float
        Mass of the primary in Msolar units.
    M_sec : float
        Mass of the secondary in Msolar units.
    Mdot : float
        Rate of change of mass of the star in Msolar/year units.
        (Negative for wind mass loss.)
    R : float
        Radius of the star in Rsolar units.
    I : float
        Moment of inertia of the star in Msolar*Rsolar^2.
    tau_conv: float
        Convective turnover time of the star, calculated @
        0.5*pressure_scale_height above the bottom of the outer convection
        zone in yr.
    L : float
        Luminosity of the star in solar units.
    #mass_conv_core : float
    #    Convective core mass of the secondary in Msolar units.
    conv_mx1_top_r : float
        Coordinate of top convective mixing zone coordinate in Rsolar.
    conv_mx1_bot_r : float
        Coordinate of bottom convective mixing zone coordinate in Rsolar.
    surface_h1 : float
        surface mass Hydrogen abundance
    center_h1 : float
        center mass Hydrogen abundance
    M_env : float
        mass of the dominant convective region for tides above the core,
        in Msolar.
    DR_env : float
        thickness of the dominant convective region for tides above the core,
        in Rsolar.
    Renv_middle : float
        position of the dominant convective region for tides above the core,
        in Rsolar.
    Idot : float
        Rate of change of the moment of inertia of the star in
        Msolar*Rsolar^2 per year.
    do_wind_loss: Boolean
        If True, take into account change of separation due to mass loss from
        the secondary. Default: True.
    do_tides: Booleans
       If True, take into account change of separation, eccentricity and
       secondary spin due to tidal forces. Default: True.
    do_gravitational_radiation: Boolean
       If True, take into account change of separation and eccentricity due to
       gravitational wave radiation. Default: True
    do_magnetic_braking: Boolean
        If True, take into account change of star spin due to magnetic braking.
        Default: True.
    magnetic_braking_mode: String
        A string corresponding to the desired magnetic braking prescription.
            - RVJ83 : Rappaport, Verbunt, & Joss 1983 [4]_
            - M15 : Matt et al. 2015 [5]_
            - G18 : Garraffo et al. 2018 [6]_
            - CARB : Van & Ivanova 2019 [7]_
    do_stellar_evolution_and_spin_from_winds: Boolean
        If True, take into account change of star spin due to change of its
        moment of inertia during its evolution and due to spin angular momentum
        loss due to winds. Default: True.
    verbose : Boolean
        If we want to print stuff. Default: False.

    Returns
    -------
    list of float
        Contains the change of
        the separation, eccentricity and angular velocity, in Rsolar,
        dimensionless and rad/year units, respectively.

    References
    ----------
    .. [1] Tauris, T. M., & van den Heuvel, E. 2006,
           Compact stellar X-ray sources, 1, 623
    .. [2] Hut, P. 1981, A&A, 99, 126
    .. [3] Junker, W., & Schafer, G. 1992, MNRAS, 254, 146
    .. [4] Rappaport, S., Joss, P. C., & Verbunt, F. 1983, ApJ, 275, 713
    .. [5] Matt et al. 2015, ApJ, 799, L23
    .. [6] Garraffo et al. 2018, ApJ, 862, 90
    .. [7] Van & Ivanova 2019, ApJ, 886, L31
    .. [8] Gossage et al. 2021, ApJ, 912, 65

    """
    y[0] = np.max([y[0], 0])  # We limit separation to non-negative values
    a = y[0]
    y[1] = np.max([y[1], 0])  # We limit eccentricity to non-negative values
    e = y[1]
    if e > 0 and e < 10.0 ** (-3):
        # we force a negligible eccentricity to become 0
        # for computational stability
        e = 0.0
        if verbose and verbose != 1:
            print("negligible eccentricity became 0 for "
                  "computational stability")
    y[2] = np.max([y[2], 0])  # We limit omega spin to non-negative values
    Omega_sec = y[2]  # in rad/yr
    y[3] = np.max([y[3], 0])
    Omega_pri = y[3]

    da = 0.0
    de = 0.0
    dOmega_sec = 0.0
    dOmega_pri = 0.0
    #  Mass Loss
    if do_wind_loss:
        q1 = M_sec / M_pri
        k11 = (1 / (1 + q1)) * (Mdot_sec / M_sec)
        k21 = Mdot_sec / M_sec
        k31 = Mdot_sec / (M_pri + M_sec)
        # This is simplified to da_mt = -a * Mdot/(M+Macc), for only (negative)
        # wind Mdot from star M.
        da_mt_sec = a * (2 * k11 - 2 * k21 + k31)

        q2 = M_pri / M_sec
        k12 = (1 / (1 + q2)) * (Mdot_pri / M_pri)
        k22 = Mdot_pri / M_pri
        k32 = Mdot_pri / (M_pri + M_sec)
        da_mt_pri = a * (
                2 * k12 - 2 * k22 + k32
        )
        if verbose and verbose != 1:
            print("da_mt = ", da_mt_sec, da_mt_pri)

        da = da + da_mt_sec + da_mt_pri

    #  Tidal forces
    if do_tides:
        q1 = M_pri / M_sec
        q2 = M_sec / M_pri
        # P_orb in years. From 3rd Kepler's law, transforming separation from
        # Rsolar to AU, to avoid using constants
        # TODO: we aleady have this function!
        P_orb = np.sqrt((a / const.aursun) ** 3 / (M_pri + M_sec))
        n = 2.0 * const.pi / P_orb  # mean orbital ang. vel. in rad/year
        f1 = (
                1
                + (31 / 2) * e ** 2
                + (255 / 8) * e ** 4
                + (185 / 16) * e ** 6
                + (25 / 64) * e ** 8
        )
        f2 = 1 + (15 / 2) * e ** 2 + (45 / 8) * e ** 4 + (5 / 16) * e ** 6
        f3 = 1 + (15 / 4) * e ** 2 + (15 / 8) * e ** 4 + (5 / 64) * e ** 6
        f4 = 1 + (3 / 2) * e ** 2 + (1 / 8) * e ** 4
        f5 = 1 + 3 * e ** 2 + (3 / 8) * e ** 4

        # equilibrium timecale
        if ((M_env_sec != 0.0 and not np.isnan(M_env_sec))
                and (DR_env_sec != 0.0 and not np.isnan(DR_env_sec)) and (
                    Renv_middle_sec != 0.0 and not np.isnan(Renv_middle_sec))):
            # eq. (31) of Hurley et al. 2002, generalized for convective layers
            # not on surface too
            tau_conv_sec = 0.431 * ((M_env_sec * DR_env_sec * Renv_middle_sec
                                     / (3 * L_sec)) ** (1.0 / 3.0))
        else:
            if verbose and verbose != 1:
                print("something wrong with M_env/DR_env/Renv_middle",
                      M_env_sec, DR_env_sec, Renv_middle_sec)
            tau_conv_sec = 1.0e99
        if ((M_env_pri != 0.0 and not np.isnan(M_env_pri))
                and (DR_env_pri != 0.0 and not np.isnan(DR_env_pri)) and (
                    Renv_middle_pri != 0.0 and not np.isnan(Renv_middle_pri))):
            # eq. (31) of Hurley et al. 2002, generalized for convective layers
            # not on surface too
            tau_conv_pri = 0.431 * ((M_env_pri * DR_env_pri * Renv_middle_pri
                                     / (3 * L_pri)) ** (1.0/3.0))
        else:
            if verbose and verbose != 1:
                print("something wrong with M_env/DR_env/Renv_middle",
                      M_env_pri, DR_env_pri, Renv_middle_pri)
            tau_conv_pri = 1.0e99
        P_spin_sec = 2 * np.pi / Omega_sec
        P_spin_pri = 2 * np.pi / Omega_pri
        P_tid_sec = np.abs(1 / (1 / P_orb - 1 / P_spin_sec))
        P_tid_pri = np.abs(1 / (1 / P_orb - 1 / P_spin_pri))
        f_conv_sec = np.min([1, (P_tid_sec / (2 * tau_conv_sec)) ** 2])
        f_conv_pri = np.min([1, (P_tid_pri / (2 * tau_conv_pri)) ** 2])
        F_tid = 1.  # not 50 as before
        kT_conv_sec = (
                (2. / 21) * (f_conv_sec / tau_conv_sec) * (M_env_sec / M_sec)
        )  # eq. (30) of Hurley et al. 2002
        kT_conv_pri = (
                (2. / 21) * (f_conv_pri / tau_conv_pri) * (M_env_pri / M_pri)
        )
        if kT_conv_sec is None or not np.isfinite(kT_conv_sec):
            kT_conv_sec = 0.0
        if kT_conv_pri is None or not np.isfinite(kT_conv_pri):
            kT_conv_pri = 0.0

            if verbose:
                print("kT_conv_sec is", kT_conv_sec, ", set to 0.")
                print("kT_conv_pri is", kT_conv_pri, ", set to 0.")
        # this is the 1/timescale of all d/dt calculted below in yr^-1
        if verbose and verbose != 1:
            print(
                "Equilibrium tides in deep convective envelope",
                M_env_sec,
                DR_env_sec,
                Renv_middle_sec,
                R_sec,
                M_sec,
                M_env_pri,
                DR_env_pri,
                Renv_middle_pri,
                R_pri,
                M_pri
            )
            print("convective tiimescales and efficiencies",
                  tau_conv_sec, P_orb, P_spin_sec, P_tid_sec,
                  f_conv_sec,
                  F_tid,
                  tau_conv_pri, P_orb, P_spin_pri, P_tid_pri,
                  f_conv_pri,
                  F_tid,
                  )

        # dynamical timecale
        F_tid = 1
        # E2 = 1.592e-9*M**(2.84) # eq. (43) of Hurley et al. 2002. Deprecated
        R_conv_sec = conv_mx1_top_r_sec - conv_mx1_bot_r_sec
        R_conv_pri = conv_mx1_top_r_pri - conv_mx1_bot_r_pri  # convective core
        # R_conv = conv_mx1_top_r  # convective core
        if (R_conv_sec > R_sec or R_conv_sec <= 0.0
                or conv_mx1_bot_r_sec / R_sec > 0.1):
            # R_conv = 0.5*R
            # if verbose:
            #     print(
            #         "R_conv of the convective core is not behaving well or "
            #         "we are not calculating the convective core, we make it "
            #         "equal to half of Rstar",
            #         R_conv,
            #         R,
            #         conv_mx1_top_r,
            #         conv_mx1_bot_r,
            #     )
            # we switch to Zahn+1975 calculation of E2
            E21 = 1.592e-9 * M_sec ** (2.84)
        else:
            if R_sec <= 0:
                E21 = 0
            elif surface_h1_sec > 0.4:
                E21 = 10.0 ** (-0.42) * (R_conv_sec / R_sec) ** (
                    7.5
                )  # eq. (9) of Qin et al. 2018, 616, A28
            elif surface_h1_sec <= 0.4:  # "HeStar":
                E21 = 10.0 ** (-0.93) * (R_conv_sec / R_sec) ** (
                    6.7
                )  # eq. (9) of Qin et al. 2018, 616, A28
            else:  # in principle we should not go here
                E21 = 1.592e-9 * M_sec ** (
                    2.84
                )  # eq. (43) of Hurley et al. 2002 from Zahn+1975 Depreciated
        # kT = 1.9782e4 * np.sqrt(M * R**2 / a**5) * (1 + q)**(5. / 6) * E2
        # eq. (42) of Hurley et al. 2002. Depreciated

        if (R_conv_pri > R_pri or R_conv_pri <= 0.0
                or conv_mx1_bot_r_pri / R_pri > 0.1):
            E22 = 1.592e-9 * M_pri ** (2.84)
            if verbose and verbose != 1:
                print(
                    "R_conv of the convective core is not behaving well or we "
                    "are not calculating the convective core, we switch to "
                    "Zahn+1975 calculation of E2",
                    R_conv_sec,
                    R_sec,
                    conv_mx1_top_r_sec,
                    conv_mx1_bot_r_sec,
                    E21,
                    R_conv_pri,
                    R_pri,
                    conv_mx1_top_r_pri,
                    conv_mx1_bot_r_pri,
                    E22
                )
        else:
            if R_pri <= 0:
                E22 = 0
            elif surface_h1_pri > 0.4:
                E22 = 10.0 ** (-0.42) * (R_conv_pri / R_pri) ** (
                    7.5
                )  # eq. (9) of Qin et al. 2018, 616, A28
            elif surface_h1_pri <= 0.4:  # "HeStar":
                E22 = 10.0 ** (-0.93) * (R_conv_pri / R_pri) ** (
                    6.7
                )  # eq. (9) of Qin et al. 2018, 616, A28
            else:  # in principle we should not go here
                E22 = 1.592e-9 * M_pri ** (
                    2.84
                )
        kT_rad_sec = (
            np.sqrt(const.standard_cgrav * (M_sec * const.msol)
                    * (R_sec * const.rsol)**2 / (a * const.rsol)**5)
            * (1 + q1) ** (5.0 / 6)
            * E21
            * const.secyer)
        kT_rad_pri = (
            np.sqrt(const.standard_cgrav * (M_pri * const.msol)
                    * (R_pri * const.rsol)**2 / (a * const.rsol)**5)
            * (1 + q2) ** (5.0 / 6)
            * E22
            * const.secyer)
        # this is the 1/timescale of all d/dt calculted below in yr^-1
        if verbose and verbose != 1:
            print(
                "Dynamical tides in radiative envelope",
                conv_mx1_top_r_sec,
                conv_mx1_bot_r_sec,
                R_conv_sec,
                E21,
                conv_mx1_top_r_pri,
                conv_mx1_bot_r_pri,
                R_conv_pri,
                E22,
                F_tid
            )
        kT_sec = max(kT_conv_sec, kT_rad_sec)
        kT_pri = max(kT_conv_pri, kT_rad_pri)
        if verbose and verbose != 1:
            print("kT_conv/rad of tides is ", kT_conv_sec, kT_rad_sec,
                  kT_conv_pri, kT_rad_pri, "in 1/yr, and we picked the ",
                  kT_sec, kT_pri)

        da_tides_sec = (
                -6
                * F_tid
                * kT_sec
                * q1
                * (1 + q1)
                * (R_sec / a) ** 8
                * (a / (1 - e ** 2) ** (15 / 2))
                * (f1 - (1 - e ** 2) ** (3 / 2) * f2 * Omega_sec / n)
        )  # eq. (9) of Hut 1981, 99, 126

        da_tides_pri = (
                -6
                * F_tid
                * kT_pri
                * q2
                * (1 + q2)
                * (R_pri / a) ** 8
                * (a / (1 - e ** 2) ** (15 / 2))
                * (f1 - (1 - e ** 2) ** (3 / 2) * f2 * Omega_pri / n)
        )
        de_tides_sec = (
                -27
                * F_tid
                * kT_sec
                * q1
                * (1 + q1)
                * (R_sec / a) ** 8
                * (e / (1 - e ** 2) ** (13 / 2))
                * (f3 - (11 / 18) * (1 - e ** 2) ** (3 / 2) * f4 * Omega_sec/n)
        )  # eq. (10) of Hut 1981, 99, 126

        de_tides_pri = (
                -27
                * F_tid
                * kT_pri
                * q2
                * (1 + q2)
                * (R_pri / a) ** 8
                * (e / (1 - e ** 2) ** (13 / 2))
                * (f3 - (11 / 18) * (1 - e ** 2) ** (3 / 2) * f4 * Omega_pri/n)
        )

        dOmega_tides_sec = (
                (3 * F_tid * kT_sec * q1 ** 2 * M_sec * R_sec ** 2 / I_sec)
                * (R_sec / a) ** 6
                * n
                / (1 - e ** 2) ** 6
                * (f2 - (1 - e ** 2) ** (3 / 2) * f5 * Omega_sec / n)
        )  # eq. (11) of Hut 1981, 99, 126
        dOmega_tides_pri = (
                (3 * F_tid * kT_pri * q2 ** 2 * M_pri * R_pri ** 2 / I_pri)
                * (R_pri / a) ** 6
                * n
                / (1 - e ** 2) ** 6
                * (f2 - (1 - e ** 2) ** (3 / 2) * f5 * Omega_pri / n)
        )
        if verbose:
            print("da,de,dOmega_tides = ",
                  da_tides_sec, de_tides_sec, dOmega_tides_sec,
                  da_tides_pri, de_tides_pri, dOmega_tides_pri)
        da = da + da_tides_sec + da_tides_pri
        de = de + de_tides_sec + de_tides_pri
        dOmega_sec = dOmega_sec + dOmega_tides_sec
        dOmega_pri = dOmega_pri + dOmega_tides_pri

    #  Gravitional radiation
    if do_gravitational_radiation:
        v = (M_pri * M_sec / (M_pri + M_sec) ** 2)
        da_gr = (
            (-2 * const.clight / 15) * (v / ((1 - e ** 2) ** (9 / 2)))
            * (const.standard_cgrav * (M_pri + M_sec) * const.msol
               / (a * const.rsol * const.clight ** 2)) ** 3
            * ((96 + 292 * e ** 2 + 37 * e ** 4) * (1 - e ** 2)
               - (1 / 28 * const.standard_cgrav * (M_pri + M_sec) * const.msol
               / (a * const.rsol * const.clight ** 2))
               * ((14008 + 4707 * v)+(80124 + 21560 * v) * e ** 2
               + (17325 + 10458) * e ** 4 - 0.5 * (5501 - 1036 * v) * e ** 6))
            ) * const.secyer / const.rsol
        # eq. (35) of Junker et al. 1992, 254, 146
        de_gr = (
            (-1 / 15) * ((v * const.clight ** 3) / (
                const.standard_cgrav * (M_pri + M_sec) * const.msol))
            * (const.standard_cgrav * (M_pri + M_sec) * const.msol / (
                a * const.rsol * const.clight ** 2)) ** 4
            * (e / (1 - e ** 2) ** (7 / 2))
            * ((304 + 121 * e ** 2) * (1 - e ** 2)
               - (1 / 56 * const.standard_cgrav * (M_pri + M_sec) * const.msol
               / (a * const.rsol * const.clight ** 2))
               * (8 * (16705 + 4676 * v) + 12 * (9082 + 2807 * v) * e ** 2
               - (25211 - 3388 * v) * e ** 4))
            ) * const.secyer
        # eq. (36) of Junker et al. 1992, 254, 146
        if verbose:
            print("da,de_gr = ", da_gr, de_gr)
        da = da + da_gr
        de = de + de_gr

    #  Magnetic braking
    if do_magnetic_braking:
        # domega_mb / dt = torque_mb / I is calculated below.
        # All results are in units of [yr^-2], i.e., the amount of change
        # in Omega over 1 year.

        if magnetic_braking_mode == "RVJ83":
            # Torque from Rappaport, Verbunt, and Joss 1983, ApJ, 275, 713
            # The torque is eq.36 of Rapport+1983, with γ = 4. Torque units
            # converted from cgs units to [Msol], [Rsol], [yr] as all stellar
            # parameters are given in units of [Msol], [Rsol], [yr] and so that
            # dOmega_mb/dt is in units of [yr^-2].
            dOmega_mb_sec = (
                -3.8e-30 * (const.rsol**2 / const.secyer)
                * M_sec
                * R_sec**4
                * Omega_sec**3
                / I_sec
                * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
            )
            dOmega_mb_pri = (
                -3.8e-30 * (const.rsol**2 / const.secyer)
                * M_pri
                * R_pri**4
                * Omega_pri**3
                / I_pri
                * np.clip((1.5 - M_pri) / (1.5 - 1.3), 0, 1)
            )
            # Converting units:
            # The constant 3.8e-30 from Rappaport+1983 has units of [cm^-2 s]
            # which need to be converted...
            #
            # -3.8e-30 [cm^-2 s] * (const.rsol**2/const.secyer) -> [Rsol^-2 yr]
            # * M [Msol]
            # * R ** 4 [Rsol^4]
            # * Omega ** 3 [yr^-3]
            # / I [Msol Rsol^2 ]
            #
            # Thus, dOmega/dt comes out to [yr^-2]

        elif magnetic_braking_mode == "M15":

            # Torque prescription from Matt et al. 2015, ApJ, 799, L23
            # Constants:
            # [erg] or [g cm^2 s^-2] -> [Msol Rsol^2 yr^-2]
            K = 1.4e30 * const.secyer**2 / (const.msol * const.rsol**2)
            # m = 0.22
            # p = 2.6
            # Above constants were calibrated as in
            # Gossage et al. 2021, ApJ, 912, 65

            # TODO: I am not sure which constants are used from each reference

            # Below, constants are otherwise as assumed as in
            # Matt et al. 2015, ApJ, 799, L23
            omega_sol = 2.6e-6 * const.secyer   # [s^-1] -> [yr^-1]
            # solar rossby = 2
            # solar convective turnover time = 12.9 days
            # Rossby number saturation threshold = 0.14
            chi = 2.0 / 0.14
            tau_conv_sol = 12.9 / 365.25        # 12.9 [days] -> [yr]

            Prot_pri = 2 * np.pi / Omega_pri    # [yr]
            Rossby_number_pri = Prot_pri / tau_conv_pri
            Prot_sec = 2 * np.pi / Omega_sec    # [yr]
            Rossby_number_sec = Prot_sec / tau_conv_sec

            # critical rotation rate in rad/yr
            Omega_crit_pri = np.sqrt(
                const.standard_cgrav * M_pri * const.msol
                / ((R_pri * const.rsol) ** 3)) * const.secyer
            Omega_crit_sec = np.sqrt(
                const.standard_cgrav * M_sec * const.msol
                / ((R_sec * const.rsol) ** 3)) * const.secyer

            # omega/omega_c
            wdivwc_pri = Omega_pri / Omega_crit_pri
            wdivwc_sec = Omega_sec / Omega_crit_sec

            gamma_pri = (1 + (wdivwc_pri / 0.072)**2)**0.5
            T0_pri = K * R_pri**3.1 * M_pri**0.5 * gamma_pri**(-2 * 0.22)
            gamma_sec = (1 + (wdivwc_sec / 0.072)**2)**0.5
            T0_sec = K * R_sec**3.1 * M_sec**0.5 * gamma_sec**(-2 * 0.22)

            if (Rossby_number_sec < 0.14):
                dOmega_mb_sec = (
                    T0_sec * (chi**2.6) * (Omega_sec / omega_sol) / I_sec
                    * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
                )
            else:
                dOmega_mb_sec = (
                    T0_sec * ((tau_conv_sec/tau_conv_sol)**2.6)
                           * ((Omega_sec/omega_sol)**(2.6 + 1)) / I_sec
                    * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
                )

            if (Rossby_number_pri < 0.14):
                dOmega_mb_pri = (
                    T0_pri * (chi**2.6) * (Omega_pri / 2.6e-6) / I_pri
                    * np.clip((1.5 - M_pri) / (1.5 - 1.3), 0, 1)
                )
            else:
                dOmega_mb_pri = (
                    T0_pri * ((tau_conv_pri/tau_conv_sol)**2.6)
                           * ((Omega_pri/omega_sol)**(2.6 + 1)) / I_pri
                    * np.clip((1.5 - M_pri) / (1.5 - 1.3), 0, 1)
                )

        elif magnetic_braking_mode == "G18":

            # Torque prescription from Garraffo et al. 2018, ApJ, 862, 90
            # a = 0.03
            # b = 0.5
            # [g cm^2] -> [Msol Rsol^2]
            c = 3e41 / (const.msol * const.rsol**2)
            # Above are as calibrated in Gossage et al. 2021, ApJ, 912, 65

            Prot_pri = 2 * np.pi / Omega_pri            # [yr]
            Rossby_number_pri = Prot_pri / tau_conv_pri
            Prot_sec = 2 * np.pi / Omega_sec            # [yr]
            Rossby_number_sec = Prot_sec / tau_conv_sec

            n_pri = (0.03 / Rossby_number_pri) + 0.5 * Rossby_number_pri + 1.0
            n_sec = (0.03 / Rossby_number_sec) + 0.5 * Rossby_number_sec + 1.0

            Qn_pri = 4.05 * np.exp(-1.4 * n_pri)
            Qn_sec = 4.05 * np.exp(-1.4 * n_sec)

            dOmega_mb_sec = (
                    c * Omega_sec**3 * tau_conv_sec * Qn_sec / I_sec
                    * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
            )

            dOmega_mb_pri = (
                    c * Omega_pri**3 * tau_conv_pri * Qn_pri / I_pri
                    * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
            )

        elif magnetic_braking_mode == "CARB":

            # Torque prescription from Van & Ivanova 2019, ApJ, 886, L31
            # Based on files hosted on Zenodo:
            #         https://zenodo.org/record/3647683#.Y_TfedLMKUk,
            # with units converted from [cm], [g], [s] to [Rsol], [Msol], [yr]

            # Constants as assumed in Van & Ivanova 2019, ApJ, 886, L31
            omega_sol = 3e-6 * const.secyer         # [s^-1] -> [yr^-1]
            tau_conv_sol = 2.8e6 / const.secyer     # [s] -> yr
            K2 = 0.07**2

            tau_ratio_sec = tau_conv_sec / tau_conv_sol
            tau_ratio_pri = tau_conv_pri / tau_conv_sol
            rot_ratio_sec = Omega_sec / omega_sol
            rot_ratio_pri = Omega_pri / omega_sol

            # below in units of [Rsol yr^-1]^2
            v_esc2_sec = ((2 * const.standard_cgrav * M_sec / R_sec)
                          * (const.msol * const.secyer**2 / const.rsol**3))
            v_esc2_pri = ((2 * const.standard_cgrav * M_pri / R_pri)
                          * (const.msol * const.secyer**2 / const.rsol**3))
            v_mod2_sec = v_esc2_sec + (2 * Omega_sec**2 * R_sec**2) / K2
            v_mod2_pri = v_esc2_pri + (2 * Omega_pri**2 * R_pri**2) / K2

            # Van & Ivanova 2019, MNRAS 483, 5595 replace the magnetic field
            # with Omega * tau_conv phenomenology. Thus, the ratios
            # (rot_ratio_* and tau_ratio_*) inherently have units of Gauss
            # [cm^-0.5 g^0.5 s^-1] that needs to be converted to [Rsol],
            # [Msol], [yr]. VI2019 assume the solar magnetic field strength is
            # on average 1 Gauss.
            if (abs(Mdot_sec) > 0):
                R_alfven_div_R3_sec = (
                    R_sec**4 * rot_ratio_sec**4 * tau_ratio_sec**4
                    / (Mdot_sec**2 * v_mod2_sec)
                    * (const.rsol**2 * const.secyer / const.msol**2))
            else:
                R_alfven_div_R3_sec = 0.0

            if (abs(Mdot_pri) > 0):
                R_alfven_div_R3_pri = (
                    R_pri**4 * rot_ratio_pri**4 * tau_ratio_pri**4
                    / (Mdot_pri**2 * v_mod2_pri)
                    * (const.rsol**2 * const.secyer / const.msol**2))
            else:
                R_alfven_div_R3_pri = 0.0

            # Alfven radius in [Rsol]
            R_alfven_sec = R_sec * R_alfven_div_R3_sec**(1./3.)
            R_alfven_pri = R_pri * R_alfven_div_R3_pri**(1./3.)

            dOmega_mb_sec = (
                  (2./3.) * Omega_sec * Mdot_sec * R_alfven_sec**2 / I_sec
                  * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
            )

            dOmega_mb_pri = (
                  (2./3.) * Omega_pri * Mdot_pri * R_alfven_pri**2 / I_pri
                  * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
            )

        else:
            Pwarn("WARNING: Magnetic braking is not being calculated in the "
                  "detached step. The given magnetic_braking_mode string \"",
                  magnetic_braking_mode, "\" does not match the available "
                  "built-in cases. To enable magnetic braking, please set "
                  "magnetc_braking_mode to one of the following strings: "
                  "\"RVJ83\" for Rappaport, Verbunt, & Joss 1983"
                  "\"G18\" for Garraffo et al. 2018"
                  "\"M15\" for Matt et al. 2015"
                  "\"CARB\" for Van & Ivanova 2019", "UnsupportedModelWarning")

        if verbose and verbose != 1:
            print("magnetic_braking_mode = ", magnetic_braking_mode)
            print("dOmega_mb = ", dOmega_mb_sec, dOmega_mb_pri)
            dOmega_sec = dOmega_sec + dOmega_mb_sec
            dOmega_pri = dOmega_pri + dOmega_mb_pri

    if do_stellar_evolution_and_spin_from_winds:
        # Due to the secondary's own evolution, we have:
        # domega_spin/dt = d(Jspin/I)/dt = dJspin/dt * 1/I + Jspin*d(1/I)/dt.
        # These are the two terms calculated below.
        dOmega_spin_wind_sec = (
                2.0 / 3.0 * R_sec ** 2 * Omega_sec * Mdot_sec / I_sec
        )
        # jshell*Mdot/I : specific angular momentum of a thin sperical shell
        # * mdot  / moment of inertia
        # Omega is in rad/yr here, and R, M in solar (Mdot solar/yr).
        dOmega_deformation_sec = np.min(
            [-Omega_sec * Idot_sec / I_sec, 100]
        )
        # This is term of Jspin*d(1/I)/dt term of the domega_spin/dt. #
        # We limit its increase due to contraction to 100 [(rad/yr)/yr]
        # increase, otherwise the integrator will fail
        # (usually when we have WD formation).
        dOmega_spin_wind_pri = (
                2.0 / 3.0 * R_pri ** 2 * Omega_pri * Mdot_pri / I_pri
        )
        # jshell*Mdot/I : specific angular momentum of a thin sperical shell
        # * mdot  / moment of inertia
        # Omega is in rad/yr here, and R, M in solar (Mdot solar/yr).
        dOmega_deformation_pri = np.min(
            [-Omega_pri * Idot_pri / I_pri, 100]
        )
        if verbose:
            print(
                "dOmega_spin_wind , dOmega_deformation = ",
                dOmega_spin_wind_sec,
                dOmega_deformation_sec,
                dOmega_spin_wind_pri,
                dOmega_deformation_pri,
            )
        dOmega_sec = dOmega_sec + dOmega_spin_wind_sec
        dOmega_pri = dOmega_pri + dOmega_spin_wind_pri
        dOmega_sec = dOmega_sec + dOmega_deformation_sec
        dOmega_pri = dOmega_pri + dOmega_deformation_pri
    if verbose:
        print("a[Ro],e,Omega[rad/yr] have been =", a, e, Omega_sec, Omega_pri)
        print("da,de,dOmega (all) in 1yr is = ",
              da, de, dOmega_sec, dOmega_pri)

    return [
        da,
        de,
        dOmega_sec,
        dOmega_pri,
    ]
