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

STAR_STATE_POST_MS = [
    "H-rich_Core_H_burning",
    "H-rich_Shell_H_burning",
    "H-rich_Core_He_burning",
    "H-rich_Central_He_depleted",
    "H-rich_Central_C_depletion",
    "H-rich_non_burning"
]

STAR_STATE_POST_HeMS = [
    'stripped_He_Core_He_burning',
    'stripped_He_Central_He_depleted',
    'stripped_He_Central_C_depletion',
    'stripped_He_non_burning'
]






class isolated_step(detached_step):
    """Evolve an isolated star (a single star, a merger product, a runaway star, etc.)

    The star will be matched in the beginning of the step and will be evolved
    until core-collapse or maximum simulation time,
    based on a grid of single star HDF5 grid.

    """

    def __init__(self,
        do_wind_loss=False,
        do_tides=False,
        do_gravitational_radiation=False,
        do_magnetic_braking=False,
        *args, **kwargs):

        super().__init__(do_wind_loss=do_wind_loss,
        do_tides=do_tides,
        do_gravitational_radiation=do_gravitational_radiation,
        do_magnetic_braking=do_magnetic_braking,
        *args,
        **kwargs)

    def merged_star_properties(star1,star2):
        """
        Make assumptions about the core/total mass of the star of a merged product.

        Similar to the table of merging in BSE
        """
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

        merged_star = SingleStar()


        for s1 in LIST_ACCEPTABLE_STATES_FOR_HMS:
            for s2 in LIST_ACCEPTABLE_STATES_FOR_HMS:
                merged_star.mass = star1.mass + star2.mass
                center_he4_merged = (star1.mass * star1.center_he4 + star2.mass * star2.center_he4) / (star1.mass + star2.mass) #mass weighted
                merged_star.center_he4 = center_he4_merged
                merged_star.metallicity = star1.metallicity


                            STAR_STATES_CO

                            STAR_STATES_H_RICH

                            STAR_STATES_HE_RICH
                            -----------
                            LIST_ACCEPTABLE_STATES_FOR_HMS

                            STAR_STATE_POST_MS

                            STAR_STATE_POST_HeMS

        binary.star_1 = merged_star
        binary.star_2 = None
        return

    def initialize_isolated_binary_orbit():
        """
        and isolated star is treated as a extremely far away binary for the purpose of keeping the same code structure
        put period at extreme, and initiate detached step with one star (and one non-evolving compact object),
         with no orbital changes apart from spin change due to winds and deformation
        """
        binary = self.binary
        binary.orbital_period = 10.**99
        binary.eccentricity = 0
        binary.separation = orbital_separation_from_period()

        if binary.state == "single_star":
            binary.star_2 = None
        elif binary.state == "disrupted":
            #find which star is a CO, the other will be evolved in isolation
            if binary.star_1 in STAR_STATES_CO:
                binary.star_1 = None
            elif binary.star_2 in STAR_STATES_CO:
                binary.star_2 = None
        elif binary.state == "merged":

            merged_star_properties(binary.star_1,binary.star_2)

    super().__call__(binary)
