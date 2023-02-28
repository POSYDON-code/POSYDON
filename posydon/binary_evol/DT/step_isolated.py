"""Isolated evolution step."""


__authors__ = [
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>"
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
from posydon.binary_evol.DT.step_detached import detached_step


class isolated_step(detached_step):
    """Evolve an isolated star (a single star, a merger product, a runaway star, etc.)

    The star will be matched in the beginning of the step and will be evolved
    until core-collapse or maximum simulation time,
    based on a grid of single star HDF5 grid.

    """

    def __init__(self,
        grid_name_Hrich=None,
        grid_name_strippedHe=None,
        path=PATH_TO_POSYDON_DATA,
        #dt=None,
        #n_o_steps_history=None,
        do_wind_loss=False,
        do_tides=False,
        do_gravitational_radiation=False,
        do_magnetic_braking=False,
        *args, **kwargs):
        print('Init isolated')
        super().__init__(
        grid_name_Hrich=grid_name_Hrich,
        grid_name_strippedHe=grid_name_strippedHe,
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
        print('before initialize orbit')
        initialize_isolated_binary_orbit(binary)
        print('after initialize orbit')
        if binary.star_1 == None or binary.star_2 == None: # already one star became None in step_merged or step_initially_single
            pass
        elif binary.state == "disrupted":
            pass
        else:
            raise ValueError("In isolated step one of the two stars should be None or the the binary.state=='disrupted' ")
        print('before iso super')
        super().__call__(binary)
        

         # TODO maybe stuff after the call of the detached step


    def initialize_isolated_binary_orbit(binary):
        """
        and isolated star is treated as a extremely far away binary for the purpose of keeping the same code structure
        put period at extreme, and initiate detached step with one star (and one non-evolving compact object),
         with no orbital changes apart from spin change due to winds and deformation
        """
        binary.orbital_period = 10.**99
        print("Isolated initialize")
        binary.eccentricity = 0.0
        binary.separation = orbital_separation_from_period(binary.orbital_period, 1.,1.)
