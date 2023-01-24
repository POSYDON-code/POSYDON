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
]


import os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize
from scipy.optimize import root

from posydon.utils.data_download import PATH_TO_POSYDON_DATA
from posydon.binary_evol.binarystar import BINARYPROPERTIES
from posydon.binary_evol.singlestar import STARPROPERTIES
from posydon.interpolation import GRIDInterpolator
from posydon.interpolation.data_scaling import DataScaler
from posydon.utils.common_functions import (
    bondi_hoyle,
    orbital_period_from_separation,
    orbital_separation_from_period,
    roche_lobe_radius,
    check_state_of_star,
    PchipInterpolator2
)
from posydon.binary_evol.flow_chart import (STAR_STATES_CC)
import posydon.utils.constants as const
from posydon.binary_evol.step_detached import detached_step
from posydon.binary_evol.step_isolated import isolated_step


STAR_STATES_CO = ['BH', 'NS', 'WD']

STAR_STATES_ALL = [
    'WD',
    'NS',
    'BH',
    'H-rich_Core_H_burning',
    'H-rich_Core_He_burning',
    'H-rich_Shell_H_burning',
    'H-rich_Central_He_depleted',
    'H-rich_Shell_He_burning',
    'H-rich_Core_C_burning',
    'H-rich_Central_C_depletion',
    'H-rich_non_burning',
    'stripped_He_Core_He_burning',
    'stripped_He_Central_He_depleted',
    'stripped_He_Central_C_depletion',
    'stripped_He_non_burning'
]

STAR_STATES_CO = ['BH', 'NS', 'WD']

STAR_STATES_NOT_CO = STAR_STATES_ALL.copy()
[STAR_STATES_NOT_CO.remove(x) for x in STAR_STATES_CO]

STAR_STATES_H_RICH = STAR_STATES_NOT_CO.copy()
[STAR_STATES_H_RICH.remove(x) for x in ['stripped_He_Core_He_burning',
                                        'stripped_He_Central_He_depleted',
                                        'stripped_He_Central_C_depletion',
                                        'stripped_He_non_burning']]

STAR_STATES_HE_RICH = STAR_STATES_NOT_CO.copy()
[STAR_STATES_HE_RICH.remove(x) for x in ['H-rich_Core_H_burning',
                                         'H-rich_Core_He_burning',
                                         'H-rich_Shell_H_burning',
                                         'H-rich_Central_He_depleted',
                                         'H-rich_Shell_He_burning',
                                         'H-rich_Core_C_burning',
                                         'H-rich_Central_C_depletion']]

LIST_ACCEPTABLE_STATES_FOR_HMS = ["H-rich_Core_H_burning"]
LIST_ACCEPTABLE_STATES_FOR_HeMS = ["stripped_He_Core_He_burning"]

LIST_ACCEPTABLE_STATES_FOR_POSTMS = STAR_STATES_H_RICH.copy()
[LIST_ACCEPTABLE_STATES_FOR_POSTMS.remove(x) for x in LIST_ACCEPTABLE_STATES_FOR_HMS]

LIST_ACCEPTABLE_STATES_FOR_POSTHeMS = STAR_STATES_HE_RICH.copy()
[LIST_ACCEPTABLE_STATES_FOR_POSTHeMS.remove(x) for x in LIST_ACCEPTABLE_STATES_FOR_HeMS]

#TODO: Ability to change the psygrid of the sinlge star tracks for different evolution (mass acretors for disrupted, merger tracks for merged?)
# For his check grid_name_Hrich, grid_name_strippedHe in detached step
# But still how I point to different grids for disrupted and merged, if I initialize isolated_step once...?!!
# POSSIBLE ANSWER: MERGEING step should be a different prior step before isolataed, where the pointing to th grid also chages

class isolated_step(isolated_step):
    """Evolve an isolated star (a single star, a merger product, a runaway star, etc.)

    The star will be matched in the beginning of the step and will be evolved
    until core-collapse or maximum simulation time,
    based on a grid of single star HDF5 grid.

    """

    def __init__(self,
        grid_name_Hrich=None,
        grid_name_stripedHe=None,
        path=PATH_TO_POSYDON_DATA,
        dt=None,
        n_o_steps_history=None,
        do_wind_loss=False,
        do_tides=False,
        do_gravitational_radiation=False,
        do_magnetic_braking=False,
        *args, **kwargs):

        super().__init__(
        grid_name_Hrich=grid_name_Hrich,
        grid_name_stripedHe=grid_name_stripedHe,
        path=path,
        dt=dt,
        n_o_steps_history=n_o_steps_history,
        do_wind_loss=do_wind_loss,
        do_tides=do_tides,
        do_gravitational_radiation=do_gravitational_radiation,
        do_magnetic_braking=do_magnetic_braking,
        *args,
        **kwargs)

    def merged_star_properties(star_base,comp, position_of_star_base):
        """
        Make assumptions about the core/total mass of the star of a merged product.

        Similar to the table of merging in BSE

        star_base: Single Star
            is our base star that engulfs its companions. The merged star will have this star as a base
        comp: Single Star
            is the star that is engulfed
        position_of_star_base: int
            position in the binary of the star_base (1 if binary.star_1 etc)
        """


        '''
        MERGED_STAR= {}
            'mass': 10.0,
            'log_R': np.log10(1000.0),
            'he_core_mass': 3.0,
            'he_core_radius': 0.5,
            'state': 'H-rich_Shell_H_burning',
            'metallicity' : 0.0142,
            "log_Lnuc": -1e6, # arbitrary
            "log_LHe": -1e7, # arbitrary
            "center_he4" : 0.0,
            "center_h1" : 1.0,
            "center_c12" : 0.01,
        }
        '''
        #by default the stellar attributes that keep the same value from the
        merged_star = star_base


        '''
        for key in BINARYPROPERTIES:
            # the binary attributes that are changed in the CE step
            if key not in ["state", "event", 'V_sys', 'mass_transfer_case',
                           'nearest_neighbour_distance']:
                setattr(binary, key, np.nan)    # the rest become np.nan
            if key == 'mass_transfer_case':
                setattr(binary, key, 'None')
        stars = [donor_star, comp_star]
        # for now we just add all initial mass to the (merger) star_1
        masses = [donor_star.mass + comp_star.mass, np.nan]
        star_states = [donor_star.state, comp_star.state]
        for star, star_state, mass in zip(stars, star_states, masses):
            star.mass = mass
            star.state = star_state
            for key in STARPROPERTIES:
                # the binary attributes that are changed in the CE step
                if key not in ["mass", "state"]:
                    setattr(star, key, np.nan)
        '''


        s1 = star_base.state
        s1 = comp.state

        for s1 in LIST_ACCEPTABLE_STATES_FOR_HMS:
            for s2 in LIST_ACCEPTABLE_STATES_FOR_HMS:
                #these stellar attributes change value
                merged_star.mass = star_base.mass + comp.mass
                center_h1_merged = (star_base.mass * star_base.center_h1 + comp.mass * comp.center_h1) / (star_base.mass + comp.mass) #mass weighted
                merged_star.center_h1 = center_h1_merged
                center_he4_merged = (star_base.mass * star_base.center_he4 + comp.mass * comp.center_he4) / (star_base.mass + comp.mass) #mass weighted
                merged_star.center_he4 = center_he4_merged
                center_c12_merged = (star_base.mass * star_base.center_c12 + comp.mass * comp.center_c12) / (star_base.mass + comp.mass) #mass weighted
                merged_star.center_c12 = center_c12_merged

                #TODO log_total_angular_momentum after the merger

                #TODO: should I check if the abundaces above end up in ~1 (?)

                #TODO: WE DO NOT HAVE LOGR FOR THE MATCHING

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surface","surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(binary, key, np.nan)
                        if key in [ "c12_c12", "center_gamma",
                                   "avg_c_in_c_core", "total_moment_of_inertia", "spin", "envelope_binding_energy"]:
                            setattr(binary, key, np.nan)


        for s1 in STAR_STATE_POST_MS:
            for s2 in LIST_ACCEPTABLE_STATES_FOR_HMS:

                merged_star.mass = star_base.mass + comp.mass
                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surface","surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(binary, key, np.nan)
                        if key in [ "c12_c12", "center_gamma",
                                   "avg_c_in_c_core", "total_moment_of_inertia", "spin"]:
                            setattr(binary, key, np.nan)




    STAR_STATES_CO
    STAR_STATES_NOT_CO

    STAR_STATES_H_RICH
    [STAR_STATES_HE_RICH

    LIST_ACCEPTABLE_STATES_FOR_HMS
    LIST_ACCEPTABLE_STATES_FOR_HeMS

    LIST_ACCEPTABLE_STATES_FOR_POSTMS
    [LIST_ACCEPTABLE_STATES_FOR_POSTMS

    LIST_ACCEPTABLE_STATES_FOR_POSTHeMS
    [LIST_ACCEPTABLE_STATES_FOR_POSTHeMS



        if position_of_star_base == 1:
            binary.star_1 = merged_star
            binary.star_2 = None
        elif position_of_star_base == 2:
            binary.star_1 = None
            binary.star_2 = merged_star
        else:
            raise ValueError("position_of_star_base != 1 or 2")


        for key in STARPROPERTIES:
            # the binary attributes that are changed in the CE step
            if key not in ["separation", "orbital_period",
                           "eccentricity", "state", "event"]:
                # the binary attributes that keep the same value from the
                # previous step
                if key not in ["time", "V_sys", "mass_transfer_case",
                               "nearest_neighbour_distance"]:
                    setattr(binary, key, np.nan)  # the rest become np.nan


        return

    def __call__(self,binary):

        if binary.state == "merged":
            if binary.event == 'oMerging1':
                merged_star_properties(binary.star_1,binary.star_2,1)
            elif binary.event == 'oMerging2':
                merged_star_properties(binary.star_2,binary.star_1,2)
            else:
                raise ValueError("binary.state='merged' but binary.event != 'oMerging1/2'")
        else:
            raise ValueError("step_merging initated but binary.state != 'merged'")

        binary.event == None 

        super().__call__(binary)
