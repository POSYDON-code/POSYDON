"""Detached evolution for double compact-object binaries.
"""


__authors__ = [
    "Zepei Xing <Zepei.Xing@unige.ch>",
    "Devina Misra <devina.misra@unige.ch>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
]


import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

import posydon.utils.constants as constants
from posydon.binary_evol.DT.step_detached import detached_evolution, detached_step
from posydon.utils.common_functions import (
    CO_radius,
    orbital_period_from_separation,
    set_binary_to_failed,
)
from posydon.utils.posydonerror import NumericalError


def _event(terminal, direction=0):
    """Return a decorator that sets solve_ivp event attributes."""
    def dec(f):
        f.terminal = terminal
        f.direction = direction
        return f
    return dec


class _SBrentqDenseOutput:
    """Dense output for the s = -ln(alpha) formulation of the Peters equations.

    Inverts the monotone tau(s) mapping via ``brentq`` root-finding on the
    native ODE dense output, then converts back to physical variables
    [separation_Rsun, eccentricity, ω_sec, ω_pri].

    Parameters
    ----------
    sol : OdeSolution
        Dense output from solve_ivp in s-space.
    a0_Rsun : float
        Initial separation [Rsun].
    t_scale : float
        Characteristic GW timescale [yr].
    t0_phys : float
        Physical time at integration start [yr].
    s_lo : float
        Lower bound of s integration range.
    s_hi : float
        Upper bound of s integration range.
    """

    def __init__(self, sol, a0_Rsun, t_scale, t0_phys, s_lo, s_hi):
        self.sol = sol
        self.a0_Rsun = a0_Rsun
        self.t_scale = t_scale
        self.t0_phys = t0_phys
        self.s_lo = float(s_lo)
        self.s_hi = float(s_hi)

    def __call__(self, t_phys):
        t_phys = np.atleast_1d(np.asarray(t_phys, dtype=float))
        tau_target = (t_phys - self.t0_phys) / self.t_scale

        n = len(tau_target)
        result = np.empty((4, n))
        for i in range(n):
            s_star = brentq(
                lambda s: self.sol(s)[1] - tau_target[i],
                self.s_lo, self.s_hi,
                xtol=1e-14, rtol=1e-14,
            )
            raw = self.sol(s_star)          # [l, tau, omega_sec, omega_pri]
            alpha = np.exp(-s_star)
            result[0, i] = alpha * self.a0_Rsun   # separation (Rsun)
            result[1, i] = np.exp(raw[0])          # eccentricity = exp(l)
            result[2, i] = raw[2]                  # omega_sec
            result[3, i] = raw[3]                  # omega_pri

        if result.shape[1] == 1:
            return result[:, 0]
        return result


class DoubleCO(detached_step):
    """Evolve a double compact-object binary due to gravitational radiation.

    The binary will be evolved until the two compact objects come into contact
    or until maximum simulation time, based on the quadrupole approximation
    of gravitational radiation.  The integration is performed in dimensionless
    variables to guarantee numerical stability.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # Use the DCO-specific evolution object (disables tides, winds, etc.)
        self.evo = double_CO_evolution(**self.evo_kwargs)

    def __call__(self, binary):
        super().__call__(binary)

        binary.V_sys = binary.V_sys_history[-1]

        if self.res.status == 1:
            binary.eccentricity = 0.0
            binary.state = "contact"
            binary.event = "CO_contact"

        elif self.res.status != -1:
            binary.time = self.max_time
            binary.eccentricity = self.res.y[1][-1]
            binary.state = "detached"
            binary.event = "maxtime"

    def solve_ODEs(self, binary, primary, secondary):
        """Solve the Peters (1964) GW inspiral.

        We do two transformations to the standard Peters equations:
        1. We remove the dimensionality by using the dimensionless separation \alpha = a/a0 and
           dimensionless time \tau = t / t0, where a0 is the initial separation
           and t0 is the characteristic GW timescale.
        2. We then substitute s = -ln(\alpha) to eliminate the singularity in the Peters equations,
           which removes the alpha^-3 and alpha^-4 prefactors that cause numerical issues for tight binaries.
        """
        self.max_time = binary.properties.max_simulation_time
        t0_phys = binary.time                                     # [yr]

        # --- physical scales ------------------------------------------------
        a0_Rsun = binary.separation                               # [Rsun]
        a0_km = a0_Rsun * constants.Rsun / 1e5                    # [km]
        e0 = binary.eccentricity
        a0_cgs = a0_km * 1e5                                      # [cm]

        m1_cgs = primary.mass * constants.msol                     # [g]
        m2_cgs = secondary.mass * constants.msol                   # [g]
        G = constants.standard_cgrav                               # [cm^3 g^-1 s^02]
        c = constants.clight                                       # [cm s^-1]

        # Characteristic circular inspiral timescale t0 = a0^4 / (beta) [yr],
        # where beta = (64/5) * G^3 * m1 * m2 * M / c^5.
        t_scale = ((5.0 / 64.0) * c**5 * a0_cgs**4
                   / (G**3 * (m1_cgs + m2_cgs) * m1_cgs * m2_cgs)
                   / constants.secyer)

        # dimensionless contact threshold  \alpha_contact = (r1 + r2) / a0
        r1_km = (CO_radius(primary.mass, primary.state)
                 * constants.Rsun / 1e5)                           # [km]
        r2_km = (CO_radius(secondary.mass, secondary.state)
                 * constants.Rsun / 1e5)                           # [km]

        # Unitless separation at contact.
        alpha_contact = (r1_km + r2_km) / a0_km

        # dimensionless max-time limit
        tau_max = (self.max_time - t0_phys) / t_scale

        s_contact = -np.log(alpha_contact)

        # --- max-time event: tau reaches tau_max before contact ---------------
        @_event(True, +1)
        def ev_maxtime(s, y):
            return y[1] - tau_max

        # --- integrate -------------------------------------------------------
        # State vector: [l, tau, secondary.omega, primary.omega]
        #   where l = ln(e) and tau = t/t0
        # Independent variable: s = −ln(a/a0), from 0 to s_contact
        l0 = (np.log(e0) if e0 > 0
              else np.log(np.finfo(float).tiny))
        try:
            res = solve_ivp(self.evo.rhs,
                            events=ev_maxtime,
                            method="RK45",
                            t_span=(0.0, s_contact),
                            y0=[l0, 0.0,
                                secondary.omega0, primary.omega0],
                            rtol=1e-10,
                            atol=1e-10,
                            dense_output=True)
        except Exception as exc:
            set_binary_to_failed(binary)
            raise NumericalError(
                "SciPy encountered an error while solving "
                f"GR equations (s-formulation): {exc}"
            )

        # --- map solution back to physical units ----------------------------
        s_vals = res.t.copy()
        alpha_vals = np.exp(-s_vals)
        l_vals = res.y[0].copy()
        ecc_vals = np.exp(l_vals)
        tau_vals = res.y[1].copy()

        # Wrap dense output to convert from ln-space back to physical variables
        res.sol = _SBrentqDenseOutput(res.sol, a0_Rsun, t_scale, t0_phys,
                                      s_vals[0], s_vals[-1])

        # Map solver status to the convention expected by __call__:
        #   status  0 -> reached s_contact -> contact  -> report as 1
        #   status  1 -> maxtime event     -> no merge -> report as 0
        #   status −1 -> solver failure                -> keep as −1
        if res.status == 0:
            res.status = 1   # contact
        elif res.status == 1:
            res.status = 0   # maxtime

        # Rearrange state vector in ln-space
        # from [l, tau, secondary.omega, primary.omega]
        # to   [sep_Rsun, ecc, secondary.omega, primary.omega]
        res.y[0] = alpha_vals * a0_Rsun                            # sep (Rsun)
        res.y[1] = ecc_vals                                        # ecc

        res.t = tau_vals * t_scale + t0_phys                       # tau -> yr

        # Convert events from s-space to physical units
        for i in range(len(res.t_events)):
            if len(res.t_events[i]) > 0:
                s_ev = res.t_events[i].copy()
                alpha_ev = np.exp(-s_ev)
                l_ev = res.y_events[i][:, 0].copy()
                ecc_ev = np.exp(l_ev)
                tau_ev = res.y_events[i][:, 1].copy()
                res.y_events[i][:, 0] = alpha_ev * a0_Rsun         # sep (Rsun)
                res.y_events[i][:, 1] = ecc_ev                     # ecc
                res.t_events[i] = tau_ev * t_scale + t0_phys       # yr

        return res


class double_CO_evolution(detached_evolution):
    """Evolution object for double compact-object binaries.

    Only gravitational radiation is active; tides, winds, and magnetic
    braking are disabled.  The actual ODE right-hand side is defined as a
    dimensionless local function inside ``DoubleCO.solve_ODEs``; this class
    exists so that the ``detached_step`` pipeline (track matching, property
    updates) has a valid ``evo`` object.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.do_magnetic_braking = False
        self.do_tides = False
        self.do_wind_loss = False
        self.do_stellar_evolution_and_spin_from_winds = False
        self.do_gravitational_radiation = True

    def __call__(self, s, y):
        """Differential equation describing the orbital evolution of a double
        compact object binary.

        Parameters
        ----------
        s : float
            Independent variable, s = -ln(a/a0).
        y : array_like
            State vector at s, [l, tau, secondary.omega, primary.omega],
            where l = ln(e) and tau = t / t0.
        """

        if self.do_gravitational_radiation:
            return self.rhs(s, y)


    def _g(e2):
        """Peters (1964) f(e): denominator function for the da/dt equation.

        Defined as f(e) = 1 + (73/24)e^2 + (37/96)e^4 in Peters (1964), Eq. 5.11.
        Named _g here to match the GW-integration convention.
        """
        return 1.0 + (73.0 / 24.0) * e2 + (37.0 / 96.0) * e2 * e2


    def _f(e2):
        """Peters (1964) g(e): numerator function for the de/dt equation.

        Defined as g(e) = 1 + (121/304)e^2 in Peters (1964), Eq. 5.13.
        Named _f here to match the GW-integration convention.
        """
        return 1.0 + (121.0 / 304.0) * e2

    def rhs(self, s, y):
        """Right-hand side of the ODE system in ln-space.

        Parameters
        ----------
        s : float
            Independent variable, s = -ln(a/a0).
        y : array_like
            State vector at s, [l, tau, secondary.omega, primary.omega],
            where l = ln(e).
        """
        l = y[0]
        e = np.exp(l)
        e2 = e * e

        if e2 >= 1.0:
            return [0.0, 0.0, 0.0, 0.0]

        one_minus_e2 = 1.0 - e2
        G = self._g(e2)
        F = self._f(e2)

        dl_ds = -(19.0 / 12.0) * one_minus_e2 * F / G
        dtau_ds = np.exp(-4.0 * s) * one_minus_e2**3.5 / G

        # primary.omega0 and secondary.omega0 are constant for DCO (no tides / MB)
        return [dl_ds, dtau_ds, 0.0, 0.0]
