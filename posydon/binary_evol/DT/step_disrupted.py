"""Merging and isolated evolution step."""


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
from posydon.binary_evol.DT.step_isolated import isolated_step

import warnings

from posydon.binary_evol.flow_chart import (STAR_STATES_ALL,
    STAR_STATES_CO,
    STAR_STATES_H_RICH,
    STAR_STATES_HE_RICH,
    STAR_STATES_NOT_CO
    )



LIST_ACCEPTABLE_STATES_FOR_HMS = ["H-rich_Core_H_burning"]
LIST_ACCEPTABLE_STATES_FOR_HeMS = ["stripped_He_Core_He_burning"]

LIST_ACCEPTABLE_STATES_FOR_POSTMS = STAR_STATES_H_RICH.copy()
[LIST_ACCEPTABLE_STATES_FOR_POSTMS.remove(x) for x in LIST_ACCEPTABLE_STATES_FOR_HMS]

LIST_ACCEPTABLE_STATES_FOR_POSTHeMS = STAR_STATES_HE_RICH.copy()
[LIST_ACCEPTABLE_STATES_FOR_POSTHeMS.remove(x) for x in LIST_ACCEPTABLE_STATES_FOR_HeMS]

class disrupted_step(isolated_step):
    """
    Prepare a runaway star to do an an isolated_step)
    """

    def __init__(self,
        grid_name_Hrich=None,
        grid_name_strippedHe=None,
        path=PATH_TO_POSYDON_DATA,
        *args, **kwargs):
        super().__init__(
        grid_name_Hrich=grid_name_Hrich,
        grid_name_strippedHe=grid_name_strippedHe,
        *args,
        **kwargs)

    def __call__(self,binary):

        '''
        if binary.state == "disrupted":
            #find which star is a CO, the other will be evolved in isolation
            if binary.star_1 in STAR_STATES_CO:   ## TODO KEEP IT AS CORE-COLLAPSE
                binary.star_1 = None
            elif binary.star_2 in STAR_STATES_CO:
                binary.star_2 = None
        '''
        print('disrupted super')
        super().__call__(binary)
