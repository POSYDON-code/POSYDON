"""Unit tests of the configurable ODE solver in
posydon/binary_evol/DT/step_detached.py (detached_step.solve_ODEs)."""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>"
]

from types import SimpleNamespace
from unittest import mock

import numpy as np
from pytest import fixture, raises
from scipy.integrate import solve_ivp

import posydon.binary_evol.DT.step_detached as totest


def _never(t, y):
    return 1.0


_never.terminal = True
_never.direction = 0


class ToyEvo:
    """Stand-in for detached_evolution: a 4-variable ODE with events."""

    def __init__(self, stiffness=1e-5):
        self.k = stiffness
        self.ev_rlo1 = self.ev_rlo2 = _never
        self.ev_max_time1 = self.ev_max_time2 = _never

    def __call__(self, t, y):
        return [-1e-9 * y[0], -1e-9 * y[1],
                -self.k * (y[2] - 1.0 - np.cos(t * 1e-6)), -1e-8 * y[3]]


def _setup(stiffness=1e-5, verbose=False, **kwargs):
    step = totest.detached_step(**kwargs)
    step.verbose = verbose        # set after __init__ (no track_matcher here)
    step.evo = ToyEvo(stiffness)
    step.max_time = 1e6
    binary = SimpleNamespace(time=0.0, separation=30.0, eccentricity=0.1)
    primary = SimpleNamespace(omega0=2.0)
    secondary = SimpleNamespace(omega0=3.0)
    return step, binary, primary, secondary


class TestODEOptions:

    def test_defaults(self):
        step = totest.detached_step()
        assert step.ode_method == "Radau"
        assert step.ode_rtol == 1e-3 and step.ode_atol == 1e-6
        assert step.ode_fallback_method == "RK45"
        assert step.ode_max_nfev is None

    def test_default_identical_to_previous_call(self):
        # previously: solve_ivp(..., method="Radau") with scipy's default
        # tolerances; the new explicit defaults must give the same numbers
        step, b, p, s = _setup()
        res = step.solve_ODEs(b, p, s)
        old = solve_ivp(step.evo, events=[_never] * 4, method="Radau",
                        t_span=(0.0, 1e6), y0=[30.0, 0.1, 3.0, 2.0],
                        dense_output=True)
        assert np.array_equal(res.t, old.t)
        assert np.array_equal(res.y, old.y)
        assert res.nfev == old.nfev

    def test_method_and_tolerances_passed(self):
        step, b, p, s = _setup(ode_method="RK45", ode_rtol=1e-4,
                               ode_atol=1e-7)
        with mock.patch.object(totest, "solve_ivp", wraps=solve_ivp) as m:
            res = step.solve_ODEs(b, p, s)
        assert m.call_count == 1
        kw = m.call_args.kwargs
        assert kw["method"] == "RK45"
        assert kw["rtol"] == 1e-4 and kw["atol"] == 1e-7
        assert res.status in (0, 1)

    def test_evaluation_limit_triggers_fallback(self, capsys):
        step, b, p, s = _setup(stiffness=1e2, ode_method="RK45",
                               ode_fallback_method="Radau",
                               ode_max_nfev=200, verbose=True)
        with mock.patch.object(totest, "solve_ivp", wraps=solve_ivp) as m:
            res = step.solve_ODEs(b, p, s)
        methods = [c.kwargs["method"] for c in m.call_args_list]
        assert methods == ["RK45", "Radau"]
        assert res.status in (0, 1)
        assert "RK45 failed" in capsys.readouterr().out

    def test_limit_not_reached(self):
        step, b, p, s = _setup(ode_method="RK45", ode_max_nfev=10**5)
        with mock.patch.object(totest, "solve_ivp", wraps=solve_ivp) as m:
            step.solve_ODEs(b, p, s)
        assert m.call_count == 1

    def test_limit_evaluations_helper(self):
        f = totest._limit_evaluations(lambda t, y: y, 2)
        assert f(0, 1) == 1 and f(0, 2) == 2
        with raises(totest.ODEEvaluationLimitError):
            f(0, 3)
