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
LIST_ACCEPTABLE_STATES_FOR_HeMS = ["stripped_He_Core_He_burning"]

LIST_ACCEPTABLE_STATES_FOR_POSTMS = STAR_STATES_H_RICH.copy()
[LIST_ACCEPTABLE_STATES_FOR_POSTMS.remove(x) for x in LIST_ACCEPTABLE_STATES_FOR_HMS]

LIST_ACCEPTABLE_STATES_FOR_POSTHeMS = STAR_STATES_HE_RICH.copy()
[LIST_ACCEPTABLE_STATES_FOR_POSTHeMS.remove(x) for x in LIST_ACCEPTABLE_STATES_FOR_HeMS]

#TODO: Ability to change the psygrid of the sinlge star tracks for different evolution (mass acretors for disrupted, merger tracks for merged?)
# For his check grid_name_Hrich, grid_name_strippedHe in detached step
# But still how I point to different grids for disrupted and merged, if I initialize isolated_step once...?!!
# POSSIBLE ANSWER: MERGEING step should be a different prior step before isolataed, where the pointing to th grid also chages

class isolated_step(detached_step):
    """Evolve an isolated star (a single star, a merger product, a runaway star, etc.)

    The star will be matched in the beginning of the step and will be evolved
    until core-collapse or maximum simulation time,
    based on a grid of single star HDF5 grid.

    """

    def __init__(self,
        grid_name_Hrich=None,
        grid_name_stripedHe=None,
        path=PATH_TO_POSYDON_DATA,
        #dt=None,
        #n_o_steps_history=None,
        do_wind_loss=False,
        do_tides=False,
        do_gravitational_radiation=False,
        do_magnetic_braking=False,
        *args, **kwargs):

        super().__init__(
        grid_name_Hrich=grid_name_Hrich,
        grid_name_stripedHe=grid_name_stripedHe,
        path=path,
        #dt=dt,
        #n_o_steps_history=n_o_steps_history,
        do_wind_loss=do_wind_loss,
        do_tides=do_tides,
        do_gravitational_radiation=do_gravitational_radiation,
        do_magnetic_braking=do_magnetic_braking,
        *args,
        **kwargs)



     def __call__(self, binary):

         initialize_isolated_binary_orbit()

         if binary.star_1 == None or binary.star_2 == None: # already one star became None in step_merging
             continue
         elif binary.state == "initially_single_star":
             binary.star_2 = None
         elif binary.state == "disrupted":
             #find which star is a CO, the other will be evolved in isolation
             if binary.star_1 in STAR_STATES_CO:
                 binary.star_1 = None
             elif binary.star_2 in STAR_STATES_CO:
                 binary.star_2 = None
        '''
        elif binary.state == "merged":
            if binary.event == 'oMerging1':
                merged_star_properties(binary.star_1,binary.star_2,1)
            elif binary.event == 'oMerging2':
                merged_star_properties(binary.star_2,binary.star_1,2)
            else:
                raise ValueError("binary.state='merged' but binary.event != 'oMerging1/2'")
        '''

         super().__call__(binary)


    def initialize_isolated_binary_orbit():
        """
        and isolated star is treated as a extremely far away binary for the purpose of keeping the same code structure
        put period at extreme, and initiate detached step with one star (and one non-evolving compact object),
         with no orbital changes apart from spin change due to winds and deformation
        """
        binary = self.binary
        binary.orbital_period = 10.**99
        binary.eccentricity = 0
        binary.separation = orbital_separation_from_period(binary.orbital_period, 1.,1.)
