"""Isolated evolution step."""


__authors__ = [
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>"
]

import numpy as np

from posydon.utils.data_download import PATH_TO_POSYDON_DATA
from posydon.utils.common_functions import orbital_separation_from_period
from posydon.binary_evol.DT.step_detached import detached_step
from posydon.utils.posydonerror import FlowError


class IsolatedStep(detached_step):
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
        """
            and isolated star is treated as a extremely far away binary for the purpose of keeping the same code structure
            put period at extreme, and initiate detached step with one star (and one non-evolving compact object),
            with no orbital changes apart from spin change due to winds and deformation

        """

        self.initialize_isolated_binary_orbit(binary)

        if binary.state == 'initially_single_star' or binary.state == 'merged':
            pass
            '''
            if binary.star_1.state.state == 'massless_remnant' or binary.star_2.state == 'massless_remnant':
                pass
            else:
                raise FlowError("In merged or initially single stars, step one of the two stars should be 'massless_remnant' ")
            '''
        elif binary.state == "disrupted":
            pass
        else:
            raise FlowError("In isolated step binary.state=='disrupted' or 'initially_single_star' or 'merged' ")

        super().__call__(binary)

        self.re_erase_isolated_binary_orbit(binary)


    def initialize_isolated_binary_orbit(self,binary):
        # I give values to the orbital parameters so that the detached step will not complain
        binary.orbital_period = 10.**99
        binary.eccentricity = 0.0
        binary.separation = orbital_separation_from_period(binary.orbital_period, binary.star_1.mass, binary.star_2.mass)

    def re_erase_isolated_binary_orbit(self,binary):
        # I give values to the orbital parameters so that the detached step will not complain
        binary.orbital_period = np.nan
        binary.eccentricity = np.nan
        binary.separation = np.nan
