"""Unit tests for direct-collapse kick suppression in StepSN."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import numpy as np

from posydon.binary_evol.SN.step_SN import BINARYPROPERTIES, StepSN


class FakeStar:
    """Minimal stand-in for a component star."""

    def __init__(self, state, f_fb, mass, he_core_mass=20.0):
        self.state = state
        self.f_fb = f_fb
        self.mass = mass
        self.he_core_mass = he_core_mass
        self.co_core_mass = mass
        self.SN_type = "CCSN"
        self.state_history = ["stripped_He_Core_He_depleted"]
        self.mass_history = [30.0]
        self.natal_kick_velocity = None
        self.natal_kick_azimuthal_angle = None
        self.natal_kick_polar_angle = None
        self.natal_kick_mean_anomaly = None
        self.spin_orbit_tilt_first_SN = None
        self.spin_orbit_tilt_second_SN = None


class FakeCompanion(FakeStar):
    """Companion that must never be treated as the collapsing star."""

    def __init__(self, mass=15.0):
        super().__init__("H-rich_Core_He_depleted", f_fb=0.0, mass=mass)


class FakeBinary:
    """Minimal binary driven through a single CC1 event."""

    def __init__(self, star_1, star_2):
        self.star_1 = star_1
        self.star_2 = star_2
        self.event = "CC1"
        self.state = "detached"
        self.eccentricity = 0.0
        self.separation = 1.0e4
        self.time_history = [1.0]
        self.orbital_period = None
        self.mass_transfer_case = "None"
        self.first_SN_already_occurred = False
        self.V_sys = None
        for key in BINARYPROPERTIES:
            if not hasattr(self, key):
                setattr(self, key, None)


def make_step(no_kick_direct_collapse):
    return StepSN(
        kick=True,
        kick_normalisation="one_over_mass",
        kick_prescription="maxwellian",
        sigma_kick_CCSN_BH=265.0,
        no_kick_direct_collapse=no_kick_direct_collapse,
        RNG=np.random.default_rng(1234),
        verbose=False,
    )


def test_default_no_kick_direct_collapse_false():
    assert StepSN().no_kick_direct_collapse is False


def test_direct_collapse_no_kick_when_enabled():
    step = make_step(no_kick_direct_collapse=True)
    binary = FakeBinary(FakeStar(state="BH", f_fb=1.0, mass=10.0),
                        FakeCompanion())
    step.orbital_kick(binary)
    assert binary.star_1.natal_kick_velocity == 0.0


def test_direct_collapse_kick_when_disabled():
    step = make_step(no_kick_direct_collapse=False)
    binary = FakeBinary(FakeStar(state="BH", f_fb=1.0, mass=10.0),
                        FakeCompanion())
    step.orbital_kick(binary)
    assert binary.star_1.natal_kick_velocity > 0.0


def test_fallback_bh_keeps_kick_when_enabled():
    step = make_step(no_kick_direct_collapse=True)
    binary = FakeBinary(FakeStar(state="BH", f_fb=0.99, mass=10.0),
                        FakeCompanion())
    step.orbital_kick(binary)
    assert binary.star_1.natal_kick_velocity > 0.0
