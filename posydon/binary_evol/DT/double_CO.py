"""Detached evolution for double compact-object binaries."""


__authors__ = [
    "Zepei Xing <Zepei.Xing@unige.ch>",
    "Devina Misra <devina.misra@unige.ch>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
]


import numpy as np
from scipy.integrate import solve_ivp

import posydon.utils.constants as constants
from posydon.binary_evol.binarystar import BINARYPROPERTIES
from posydon.binary_evol.DT.step_detached import detached_evolution, detached_step
from posydon.binary_evol.singlestar import STARPROPERTIES
from posydon.utils.common_functions import (
    CO_radius,
    orbital_period_from_separation,
    set_binary_to_failed,
)
from posydon.utils.posydonerror import NumericalError


def event(terminal, direction=0):
    """Return a helper function to set attributes for solve_ivp events."""
    def dec(f):
        f.terminal = terminal
        f.direction = direction
        return f
    return dec

class DoubleCO(detached_step):
    """Evolve a double compact-object binary due to gravitational radiation.

    The binary will be evolved until the two compact objects come into contact
    or until maximum simulation time, based on the quadrupole approximation
    of gravitational radiation.

    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        # Reassign detached_evolution to handle
        # special DCO evolution
        self.evo = double_CO_evolution(**self.evo_kwargs)

    def __call__(self, binary):

        super().__call__(binary)

        binary.separation = self.res.y[0][-1] * 100_000 / constants.Rsun
        binary.orbital_period = orbital_period_from_separation(
            binary.separation, binary.star_1.mass, binary.star_2.mass
        )
        binary.V_sys = binary.V_sys_history[-1]

        #if self.res.status == -1:
        #    set_binary_to_failed(binary)
        #    raise NumericalError(f"Integration failed for {binary.state} DCO "
        #                         f"({binary.star_1.state}, {binary.star_2.state}): "
        #                         f"{self.res.message}")
        # contact event triggered
        if self.res.status == 1:
            binary.time = self.res.t[-1]
            binary.eccentricity = 0.0
            binary.state = "contact"
            binary.event = "CO_contact"
        # max time reached
        elif self.res.status != -1:
            binary.time = self.max_time
            binary.eccentricity = self.res.y[1][-1]
            binary.state = "detached"
            binary.event = "maxtime"

    def solve_ODEs(self, binary, primary, secondary):

        self.max_time = binary.properties.max_simulation_time
        time_sol = []
        sol = []
        t0 = binary.time
        n = 0
        a = binary.separation * constants.Rsun / 100_000
        e = binary.eccentricity
        status = -1

        while status == -1 and n < 6:
            try:
                res = solve_ivp(self.evo,
                                events=self.evo.ev_contact,
                                method="BDF",
                            t_span=(0, self.max_time - t0),
                            y0=[a, e,
                                secondary.omega0, primary.omega0],
                            rtol=1e-10,
                            atol=1e-10,
                            dense_output=True)

            except Exception as e:
                set_binary_to_failed(binary)
                raise NumericalError(f"SciPy encountered termination edge case while solving GR equations: {e}")

            t0 += res.t[-1]
            status = res.status
            n += 1
            a = res.y[0][-1]
            e = res.y[1][-1]
            time_sol.append(res.t)
            sol.append(res)

        return res



class double_CO_evolution(detached_evolution):

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        # For DCO system, only gravitational radiation is considered
        self.do_magnetic_braking = False
        self.do_tides = False
        self.do_wind_loss = False
        self.do_stellar_evolution_and_spin_from_winds = False
        self.do_gravitational_radiation = True


    def set_stars(self, primary, secondary, t0=0.0):

        super().set_stars(primary, secondary, t0)

        # calculate CO radii in kilometers (km rather than Rsol for better scaling in ODE solutions)
        self.r1 = CO_radius(self.primary.mass, self.primary.state) * constants.Rsun / 100_000 # [km]
        self.r2 = CO_radius(self.secondary.mass, self.secondary.state) * constants.Rsun / 100_000 # [km]

    @event(True, -1)
    def ev_contact(self, t, y):
        # stop when binary separation = r1+r2 [km]
        return y[0] - (self.r1 + self.r2)

    def gravitational_radiation(self):
        """Calculate the change in orbital separation and eccentricity due to
        gravitational wave radiation from two point masses. Calculated as
        according to [1].

            [1] Peters, P. C. (1964). Gravitational Radiation and the
                Motion of Two Point Masses. Physical Review, 136(4B),
                B1224–B1232
        """
        g = constants.standard_cgrav
        c = constants.clight

        a = self.a
        e = self.e

        m1 = self.primary.latest["mass"] * constants.msol
        m2 = self.secondary.latest["mass"] * constants.msol

        a = a * 100_000

        da_gr = ((-64 / 5)
                * (g ** 3 * (m1 + m2) * m1 * m2)
                / (1 - e ** 2) ** (7 / 2) / a ** 3 / c ** 5
                * (1 + (73 / 24) * e ** 2 + (37 / 96) * e ** 4)
                * (constants.secyer / 100_000))
        de_gr = ((-304 / 15)
                * e * (g ** 3 * (m1 + m2) * m1 * m2)
                / (1 - e ** 2) ** (5 / 2) / a ** 4 / c ** 5
                * (1 + (121 / 304) * e ** 2)
                * constants.secyer)

        self.da += da_gr
        self.de += de_gr
