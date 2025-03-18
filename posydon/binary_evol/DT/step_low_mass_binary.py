


__authors__ = [
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
from posydon.interpolation.data_scaling import DataScaler
from posydon.utils.common_functions import (
    bondi_hoyle,
    orbital_period_from_separation,
    orbital_separation_from_period,
    roche_lobe_radius,
    check_state_of_star
)
from posydon.binary_evol.flow_chart import (STAR_STATES_CC)
import posydon.utils.constants as const
from posydon.binary_evol.DT.step_detached import detached_step

class LowMassBinaryStep(detached_step):
    
    """ 
    Evolving the low mass binaries as detached systems. 
    Low mass for now is consindered a binary with m1 < 5Mo. 
    This step is directed in the detached to go and do an detached evolution. 
    Both of the stars will go through matching to a single star grid. 
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
        print('Inside init!')



    def __call__(self, binary):
        """
            and isolated star is treated as a extremely far away binary for the purpose of keeping the same code structure
            put period at extreme, and initiate detached step with one star (and one non-evolving compact object),
            with no orbital changes apart from spin change due to winds and deformation

        """

        if binary.state == 'low_mass_binary':
            pass
        else:
            raise ValueError("System should be in a detached state!")

        super().__call__(binary)


  


