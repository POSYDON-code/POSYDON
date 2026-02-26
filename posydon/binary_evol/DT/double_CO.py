"""Detached evolution for double compact-object binaries.

The ODEs are non-dimensionalized using the initial orbital separation a₀
and the characteristic gravitational-wave inspiral timescale

    t₀ = (5/64) c⁵ a₀⁴ / (G³ M m₁ m₂)

so that the Peters (1964) equations become unit-free with O(1) coefficients:

    dα/dτ = −1 / [α³ (1−e²)^(7/2)] × (1 + 73/24 e² + 37/96 e⁴)
    de/dτ = −(19/12) e / [α⁴ (1−e²)^(5/2)] × (1 + 121/304 e²)

where α = a/a₀ and τ = t/t₀.  This avoids the "step size smaller than
machine epsilon" failures that occur when solve_ivp works with physical
units spanning many orders of magnitude.
"""


__authors__ = [
    "Zepei Xing <Zepei.Xing@unige.ch>",
    "Devina Misra <devina.misra@unige.ch>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
]


import numpy as np
from scipy.integrate import solve_ivp

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


class _ScaledDenseOutput:
    """Wrapper around OdeSolution that maps dimensionless ↔ physical units.

    Parameters
    ----------
    sol : OdeSolution
        The dense output from solve_ivp in dimensionless variables.
    a0_km : float
        Initial separation scale factor [km].
    t_scale : float
        Characteristic timescale [yr].
    t0_phys : float
        Physical time at the start of integration [yr].
    """

    def __init__(self, sol, a0_km, t_scale, t0_phys):
        self.sol = sol
        self.a0_km = a0_km
        self.t_scale = t_scale
        self.t0_phys = t0_phys

    def __call__(self, t_phys):
        tau = (np.asarray(t_phys) - self.t0_phys) / self.t_scale
        y = self.sol(tau)
        y[0] = y[0] * self.a0_km          # α → km
        return y


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

        # Convert final separation from km → Rsun
        binary.separation = self.res.y[0][-1] * 1e5 / constants.Rsun
        binary.orbital_period = orbital_period_from_separation(
            binary.separation, binary.star_1.mass, binary.star_2.mass
        )
        binary.V_sys = binary.V_sys_history[-1]

        # contact event triggered
        if self.res.status == 1:
            binary.eccentricity = 0.0
            binary.state = "contact"
            binary.event = "CO_contact"
        # max time reached (status == 0)
        elif self.res.status != -1:
            binary.time = self.max_time
            binary.eccentricity = self.res.y[1][-1]
            binary.state = "detached"
            binary.event = "maxtime"

    # ------------------------------------------------------------------
    #  Dimensionless ODE integration
    # ------------------------------------------------------------------
    def solve_ODEs(self, binary, primary, secondary):
        """Solve the dimensionless Peters (1964) equations for GW inspiral.

        The equations are scaled so that all dependent variables are O(1),
        eliminating the "step size < machine epsilon" problem that plagues
        the dimensional formulation.  A single call to ``solve_ivp`` is
        sufficient.
        """
        self.max_time = binary.properties.max_simulation_time
        t0_phys = binary.time                                     # [yr]

        # --- physical scales ------------------------------------------------
        a0_km = binary.separation * constants.Rsun / 1e5          # [km]
        e0 = binary.eccentricity
        a0_cgs = a0_km * 1e5                                      # [cm]

        m1_cgs = primary.mass * constants.msol                     # [g]
        m2_cgs = secondary.mass * constants.msol                   # [g]
        G = constants.standard_cgrav                               # [cm³ g⁻¹ s⁻²]
        c = constants.clight                                       # [cm s⁻¹]

        # characteristic GW inspiral timescale  [yr]
        t_scale = ((5.0 / 64.0) * c**5 * a0_cgs**4
                   / (G**3 * (m1_cgs + m2_cgs) * m1_cgs * m2_cgs)
                   / constants.secyer)

        # dimensionless contact threshold  α_contact = (r1 + r2) / a0
        r1_km = (CO_radius(primary.mass, primary.state)
                 * constants.Rsun / 1e5)                           # [km]
        r2_km = (CO_radius(secondary.mass, secondary.state)
                 * constants.Rsun / 1e5)                           # [km]
        alpha_contact = (r1_km + r2_km) / a0_km

        # dimensionless integration limits
        tau_max = (self.max_time - t0_phys) / t_scale

        # --- dimensionless RHS (Peters 1964) --------------------------------
        def rhs(tau, y):
            alpha, ecc = y[0], y[1]
            e2 = ecc * ecc
            one_minus_e2 = 1.0 - e2

            dalpha = (-1.0
                      / (alpha**3 * one_minus_e2**3.5)
                      * (1.0 + (73.0 / 24.0) * e2
                         + (37.0 / 96.0) * e2 * e2))

            decc = (-(19.0 / 12.0) * ecc
                    / (alpha**4 * one_minus_e2**2.5)
                    * (1.0 + (121.0 / 304.0) * e2))

            # ω_sec and ω_pri are constant for DCO (no tides / MB)
            return [dalpha, decc, 0.0, 0.0]

        # --- contact event: α hits α_contact from above --------------------
        @_event(True, -1)
        def ev_contact(tau, y):
            return y[0] - alpha_contact

        # --- integrate -------------------------------------------------------
        try:
            res = solve_ivp(rhs,
                            events=ev_contact,
                            method="BDF",
                            t_span=(0.0, tau_max),
                            y0=[1.0, e0,
                                secondary.omega0, primary.omega0],
                            rtol=1e-10,
                            atol=1e-10,
                            dense_output=True)
        except Exception as exc:
            set_binary_to_failed(binary)
            raise NumericalError(
                "SciPy encountered an error while solving dimensionless "
                f"GR equations: {exc}"
            )

        # --- map solution back to physical units ----------------------------
        # Wrap dense output *before* mutating res.t (need original tau range)
        res.sol = _ScaledDenseOutput(res.sol, a0_km, t_scale, t0_phys)

        res.y[0] = res.y[0] * a0_km                               # α → km
        res.t = res.t * t_scale + t0_phys                          # τ → yr

        for i in range(len(res.t_events)):
            if len(res.t_events[i]) > 0:
                res.t_events[i] = res.t_events[i] * t_scale + t0_phys
        for i in range(len(res.y_events)):
            if len(res.y_events[i]) > 0:
                res.y_events[i][:, 0] *= a0_km

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
