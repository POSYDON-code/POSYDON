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
from posydon.utils.common_functions import orbital_period_from_separation
from posydon.utils.common_functions import CO_radius
from posydon.utils.common_functions import set_binary_to_failed
from posydon.utils.posydonerror import NumericalError

class DoubleCO:
    """The double compact-object step class."""

    def __init__(self, n_o_steps_interval=None):
        """Initialize a DoubleCO instance."""
        self.n_o_steps_interval = n_o_steps_interval

    def __call__(self, binary):
        """Apply the double CO step on a BinaryStar object."""
        solt = []
        sol = []
        t_final = []
        a_final = []
        e_final = []

        self.m1 = binary.star_1.mass
        self.m2 = binary.star_2.mass
        self.eccentricity = binary.eccentricity
        self.separation = binary.separation * constants.Rsun / 100000   # in km

        self.state1 = binary.star_1.state
        self.state2 = binary.star_2.state

        t_inspiral = binary.time
        max_time = binary.properties.max_simulation_time
        assert (
                max_time - binary.time > 0.0
        ), "max_time is loXer than the current time."
        r1 = CO_radius(self.m1, self.state1) * constants.Rsun / 100000  # in km
        r2 = CO_radius(self.m2, self.state2) * constants.Rsun / 100000  # in km

        @event(True, -1)
        def ev_contact(t, y):
            # stop when binary separation = r1+r2 km
            return y[0] - (r1 + r2)

        status = -1
        n = 0
        a = self.separation
        e = self.eccentricity
        # solve the equations at most 6 times
        while(status == -1 and n < 6):
            try: 
                s = solve_ivp(
                lambda t, y: gr(t, y, self.m1, self.m2,),
                t_span=[0, max_time - t_inspiral],
                y0=[a, e],
                method='BDF',
                events=ev_contact,
                rtol=1e-10,
                atol=1e-10,
                dense_output=True)
            except:
                set_binary_to_failed(binary)
                raise NumericalError("SciPy encountered termination edge case while solving GR equations")
            t_inspiral += s.t[-1]
            status = s.status
            n += 1
            a = s.y[0][-1]
            e = s.y[1][-1]
            solt.append(s.t)
            sol.append(s)
        if self.n_o_steps_interval is not None:
            for i in range(len(solt)):
                # solt[-1] = [0]
                t_step = (solt[i][-1] - solt[i][0]) / self.n_o_steps_interval
                t = np.arange(
                    solt[i][0], solt[i][-1] + t_step / 2.0, t_step)[1:]
                a = sol[i].sol(t)[0]
                e = sol[i].sol(t)[1]
                if i == 0:
                    t = t + binary.time
                if i == 1:
                    t = t + binary.time + solt[0][-1]
                elif i == 2:
                    t = t + binary.time + solt[0][-1] + solt[1][-1]
                for k in range(len(a)):
                    a_final.append(a[k])
                    e_final.append(e[k])
                    t_final.append(t[k])
        else:
            a_final = [s.y[0][-1]]
            e_final = [s.y[1][-1]]
            t_final = [t_inspiral]

        if s.status == -1:
            failed_state = binary.state
            set_binary_to_failed(binary)
            raise NumericalError(f"Integration failed for {failed_state} DCO ({self.state1}, {self.state2}): ", s.message)

        elif s.status == 1:

            if self.n_o_steps_interval is not None:
                p_history = []
                a_history = []
                for i in range(len(a_final)-1):
                    a_history.append(a_final[i] * 100000 / constants.Rsun)
                for i in range(len(a_history)):
                    p_history.append(orbital_period_from_separation(
                        a_history[i], binary.star_1.mass, binary.star_2.mass
                    ))
                for i in range(len(a_history)):
                    binary.time_history.append(t_final[i])
                    binary.separation_history.append(a_history[i])
                    binary.eccentricity_history.append(e_final[i])
                    binary.orbital_period_history.append(p_history[i])

                for key in BINARYPROPERTIES:
                    if key not in ["eccentricity", "time",
                                   "orbital_period", "separation"]:
                        current = getattr(binary, key)
                        history = [current] * len(a_history)
                        getattr(binary, key + "_history").extend(history)
                for key in STARPROPERTIES:
                    current1 = getattr(binary.star_1, key)
                    history1 = [current1] * len(a_history)
                    getattr(binary.star_1, key + "_history").extend(history1)
                    current2 = getattr(binary.star_2, key)
                    history2 = [current2] * len(a_history)
                    getattr(binary.star_2, key + "_history").extend(history2)
            binary.state = "contact"
            binary.separation = s.y[0][-1] * 100000 / constants.Rsun
            binary.time = t_inspiral
            binary.eccentricity = 0.0
            binary.orbital_period = orbital_period_from_separation(
                binary.separation, binary.star_1.mass, binary.star_2.mass)
            binary.V_sys = binary.V_sys_history[-1]
            binary.event = "CO_contact"

            for key in BINARYPROPERTIES:
                if key not in ["eccentricity", "time",
                               "orbital_period", "separation"]:
                    current = getattr(binary, key)
                    setattr(binary, key, current)
            for key in STARPROPERTIES:
                current1 = getattr(binary.star_1, key)
                current2 = getattr(binary.star_2, key)
                setattr(binary.star_1, key, current1)
                setattr(binary.star_2, key, current2)

        else:
            if self.n_o_steps_interval is not None:
                p_history = []
                a_history = []
                for i in range(len(a_final)-1):
                    a_history.append(a_final[i] * 100000 / constants.Rsun)
                for i in range(len(a_history)):
                    p_history.append(orbital_period_from_separation(
                        a_history[i], binary.star_1.mass, binary.star_2.mass
                    ))
                for i in range(len(a_history)):
                    binary.time_history.append(t_final[i])
                    binary.separation_history.append(a_history[i])
                    binary.eccentricity_history.append(e_final[i])
                    binary.orbital_period_history.append(p_history[i])

                for key in BINARYPROPERTIES:
                    if key not in [
                        "eccentricity",
                        "time",
                        "orbital_period",
                        "separation"
                    ]:
                        current = getattr(binary, key)
                        history = [current] * len(a_history)
                        getattr(binary, key + "_history").extend(history)
                for key in STARPROPERTIES:
                    current1 = getattr(binary.star_1, key)
                    history1 = [current1] * len(a_history)
                    getattr(binary.star_1, key + "_history").extend(history1)
                    current2 = getattr(binary.star_2, key)
                    history2 = [current2] * len(a_history)
                    getattr(binary.star_2, key + "_history").extend(history2)
            binary.state = "detached"
            binary.time = max_time
            binary.separation = s.y[0][-1] * 100000 / constants.Rsun
            binary.eccentricity = s.y[1][-1]
            binary.orbital_period = orbital_period_from_separation(
                binary.separation, binary.star_1.mass, binary.star_2.mass
            )
            binary.V_sys = binary.V_sys_history[-1]
            binary.event = "maxtime"

            for key in BINARYPROPERTIES:
                if key not in [
                    "eccentricity",
                    "time",
                    "orbital_period",
                    "separation"
                ]:
                    current = getattr(binary, key)
                    setattr(binary, key, current)
            for key in STARPROPERTIES:
                current1 = getattr(binary.star_1, key)
                current2 = getattr(binary.star_2, key)
                setattr(binary.star_1, key, current1)
                setattr(binary.star_2, key, current2)


def gr(t, y, M_acc, M):
    """TODO: add description and reference for the equations."""
    g = constants.standard_cgrav
    c = constants.clight

    y[0] = np.max(y[0], 0)
    a = y[0]
    y[1] = np.max(y[1], 0)
    e = y[1]

    da = 0
    de = 0
    M_acc = M_acc * constants.msol
    M = M * constants.msol
    a = a * 100000
    da_gr = ((-64 / 5)
             * (g ** 3 * (M_acc + M) * M_acc * M)
             / (1 - e ** 2) ** (7 / 2) / a ** 3 / c ** 5
             * (1 + (73 / 24) * e ** 2 + (37 / 96) * e ** 4)
             * (constants.secyer / 100000))
    de_gr = ((-304 / 15)
             * e * (g ** 3 * (M_acc + M) * M_acc * M)
             / (1 - e ** 2) ** (5 / 2) / a ** 4 / c ** 5
             * (1 + (121 / 304) * e ** 2)
             * constants.secyer)

    da = da_gr
    de = de_gr

    return [da, de]


def event(terminal, direction=0):
    """Return a helper function to set attributes for `solve_ivp` events."""
    def dec(f):
        f.terminal = terminal
        f.direction = direction
        return f

    return dec
