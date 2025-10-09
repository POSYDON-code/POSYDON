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
from posydon.binary_evol.singlestar import STARPROPERTIES
from posydon.utils.common_functions import (
    CO_radius,
    orbital_period_from_separation,
    set_binary_to_failed,
)
from posydon.utils.posydonerror import NumericalError
from posydon.binary_evol.DT.step_detached import detached_step, detached_evolution

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
        self.evo = double_CO_evolution(**self.evo_kwargs)

    def __call__(self, binary):
        super().__call__(binary)

        binary.state = "detached"
        binary.time = self.max_time
        binary.separation = self.res.y[0][-1] * 100000 / constants.Rsun
        binary.eccentricity = self.res.y[1][-1]
        binary.orbital_period = orbital_period_from_separation(
            binary.separation, binary.star_1.mass, binary.star_2.mass
        )
        binary.V_sys = binary.V_sys_history[-1]
        binary.event = "maxtime"

    def solve_ODEs(self, binary, primary, secondary):

        self.max_time = binary.properties.max_simulation_time

        try:
            res = solve_ivp(self.evo, 
                            events=self.evo.ev_contact,
                            method="BDF", 
                        t_span=(0, self.max_time - binary.time),
                        y0=[binary.separation * constants.Rsun / 100000, 
                            binary.eccentricity,
                            secondary.omega0, primary.omega0],
                        rtol=1e-10,
                        atol=1e-10,
                        dense_output=True)
        except Exception as e:
            set_binary_to_failed(binary)
            raise NumericalError(f"SciPy encountered termination edge case while solving GR equations: {e}")

        
        #res.y[0] = res.y[0] * 100000 / constants.Rsun  # convert back to Rsun
        #print("!!!", res.y[0][-1])

        #binary.state = "detached"
        #binary.time = self.max_time
        #binary.separation = res.y[0][-1] * 100000 / constants.Rsun
        #binary.eccentricity = res.y[1][-1]
        #binary.orbital_period = orbital_period_from_separation(
        #    binary.separation, binary.star_1.mass, binary.star_2.mass
        #)
        #binary.V_sys = binary.V_sys_history[-1]
        #binary.event = "maxtime"
            
        return res
    


class double_CO_evolution(detached_evolution):

    def __init__(self, **kwargs):

        super().__init__(**kwargs)
        self.do_magnetic_braking = False
        self.do_tides = False
        self.do_wind_loss = False
        self.do_stellar_evolution_and_spin_from_winds = False
        self.do_gravitational_radiation = True

    #def __call__(self, t, y):
    #    print("CALL")
    #    result = super().__call__(t, y)

    #    return [self.da, self.de]

    def set_stars(self, primary, secondary, t0=0.0):
        """Sets memory references for primary and secondary star associated with
        this evolution. It is expected that primary/secondary have interp1d 
        objects already, as required for detached evolution.

        Parameters
        ----------
        primary : SingleStar object
            A single star object, representing the primary (more evolved) star 
            in the binary and containing its properties.
        
        secondary : SingleStar object
            A single star object, representing the secondary (less evolved) star 
            in the binary and containing its properties.

        t0 : float
            The time at the start of detached evolution. Typically should be the 
            binary.time prior to detached evolution.

        """

        super().set_stars(primary, secondary, t0)
        self.r1 = CO_radius(self.primary.mass, self.primary.state) * constants.Rsun / 100_000 # [km]
        self.r2 = CO_radius(self.secondary.mass, self.secondary.state) * constants.Rsun / 100_000 # [km]

    @event(True, -1)
    def ev_contact(self, t, y):
        # stop when binary separation = r1+r2 [km]
        return y[0] - (self.r1 + self.r2)
    
    def gravitational_radiation(self):
        """TODO: add description and reference for the equations."""
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

