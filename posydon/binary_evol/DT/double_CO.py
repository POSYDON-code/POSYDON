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

        binary.V_sys = binary.V_sys_history[-1]

        #if self.res.status == -1:
        #    set_binary_to_failed(binary)
        #    raise NumericalError(f"Integration failed for {binary.state} DCO "
        #                         f"({binary.star_1.state}, {binary.star_2.state}): "
        #                         f"{self.res.message}")
        # contact event triggered
        if self.res.status == 1:
            binary.time = self.res.real_time[-1]
            if len(self.res.real_time) > 1:
                for k, real_time in enumerate(self.res.real_time[::-1]):
                    binary.time_history[-(k+1)] = real_time
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
                                t_span=(0, self.max_time-t0),
                                y0=[a, e,
                                    secondary.omega0, primary.omega0],
                                rtol=1e-10,
                                atol=1e-10,
                                dense_output=True)
            except Exception as e:
                set_binary_to_failed(binary)
                raise NumericalError(f"SciPy encountered termination edge case while solving GR equations: {e}")

            time_sol.append(t0)
            t0 += res.t[-1]
            status = res.status
            n += 1
            a = res.y[0][-1]
            e = res.y[1][-1]
            sol.append(res)


        class CombinedSolution:
            pass
        # combine results from multiple solve_ivp calls if needed

        output_solution = CombinedSolution()
        output_solution.t = np.concatenate([s.t for s in sol])
        output_solution.y = np.hstack([s.y for s in sol])
        output_solution.status = sol[-1].status
        output_solution.message = sol[-1].message
        output_solution.t_events = sol[-1].t_events
        output_solution.y_events = sol[-1].y_events
        output_solution.success = sol[-1].success
        output_solution.time_sols = time_sol
        output_solution.time_arrs = [s.t for s in sol]

        # dynamically create a combined sol method that can interpolate across the combined solution
        def combined_sol(t):

            t = np.array(t)
            solutions = []
            for time in t:
                appended = False
                for s in sol[::-1]:
                    if 0 <= time <= s.t[-1]:
                        solutions.append(s.sol(time))
                        appended = True
                        break
                if not appended:
                    raise ValueError(f"Time {t} is out of bounds for the combined solution.")

            solutions = np.array(solutions).T

            # convert separation back from [km] to [Rsun]
            solutions[0] *= 100_000 / constants.Rsun

            return solutions

        output_solution.sol = combined_sol

        return output_solution

    def get_time_after_evo(self, binary):
        """
            After detached evolution, this uses the ODESolver result
        to determine what the current time is.

        Parameters
        ----------
        res : ODESolver object
            This is the ODESolver object produced by SciPy's
            solve_ivp function that contains calculated values
            of the stars evolution through the detached step.

        binary: BinaryStar object
            A binary star object, containing the binary system's properties.

        Returns
        -------
        t : float or array[float]
            This is the time elapsed as a result of detached
            evolution in years. This is a float unless the
            user specifies a timestep (see `n_o_steps_history`
            or `dt`) to use via the simulation properties ini
            file, in which case it is an array. For DCO objects,
            this needs to be on the time scale of the inerpolants
            from solve_ivp, so does not represent the real binary
            time. This is corrected later in the __call__ method.

        """

        t = np.array([])
        real_time = np.array([])
        if self.dt is not None and self.dt > 0:

            for k, time_arr in enumerate(self.res.time_arrs):
                t_chunk = np.arange(time_arr[0], time_arr[-1] + self.dt/2.0, self.dt)[1:]

                if len(t_chunk) == 0:
                    continue

                # ensure last element matches (rounding can mess things up for small numbers)
                if t_chunk[-1] != time_arr[-1]:
                    t_chunk[-1] = time_arr[-1]

                t = np.hstack([t, t_chunk])
                real_time = np.hstack([real_time, t_chunk + self.res.time_sols[k]])

            if t[-1] != self.res.t[-1]:
                t[-1] = self.res.t[-1]
                real_time[-1] = self.res.t[-1] + self.res.time_sols[-1]

        elif (self.n_o_steps_history is not None and self.n_o_steps_history > 0):

            for k, time_arr in enumerate(self.res.time_arrs):
                t_step = (time_arr[-1] - time_arr[0]) / self.n_o_steps_history
                t_chunk = np.arange(time_arr[0], time_arr[-1] + t_step/2.0, t_step)[1:]

                if t_chunk[-1] != time_arr[-1]:
                    t_chunk[-1] = time_arr[-1]

                t = np.hstack([t, t_chunk])
                real_time = np.hstack([real_time, t_chunk + self.res.time_sols[k]])

        else:  # self.dt is None and self.n_o_steps_history is None
            t = np.array([self.res.t[-1]])
            real_time = t + self.res.time_sols[-1]

        # store this for later to convert time back to real time
        self.res.real_time = real_time

        return t


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

        # Sanitize output to prevent propagation of inf/NaN
        if not np.isfinite(da_gr):
            da_gr = 0.0
        if not np.isfinite(de_gr):
            de_gr = 0.0

        self.da += da_gr
        self.de += de_gr
