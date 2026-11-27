"""Unit tests of posydon/binary_evol/SN/step_SN.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
    "David R. Aguilera-Dena <daguilera@astro.unam.mx>",
]

# import other needed code for the tests, which is not already imported in the
# module you like to test
import numpy as np
from pytest import approx, fixture, mark, raises, warns

# import the module which will be tested
import posydon.binary_evol.SN.step_SN as totest
import posydon.utils.posydonwarning as pw
from posydon.binary_evol.binarystar import BinaryStar
from posydon.binary_evol.singlestar import SingleStar
from posydon.binary_evol.SN.engine_data import Sukhbold16_table
from posydon.binary_evol.SN.explodability.CCSN_explosion_criteria import (
    explodes_Ertl16,
    get_CCSN_prescription,
)
from posydon.binary_evol.SN.remnant_mass.remnant_mass_prescriptions import (
    PISN_mass_Marchant,
    get_remnant_mass_prescription,
)
from posydon.utils.common_functions import (
    CO_radius,
    orbital_period_from_separation,
)
from posydon.utils.constants import Zsun
from posydon.utils.posydonerror import ModelError


def make_star(m_CO, state="stripped_He_Core_C_depleted", **kwargs):
    """A star at core collapse with a CO core of m_CO Msun.

    The He core is 1.25 m_CO + 0.4 Msun; a stripped star is its He core, a
    H-rich star has a H envelope on top.
    """
    stripped = state.startswith("stripped")
    m_He = 1.25 * m_CO + 0.4
    properties = dict(
        state=state, mass=m_He if stripped else 1.6 * m_He + 2.0,
        he_core_mass=m_He, co_core_mass=m_CO, metallicity=1.0,
        log_R=0.0 if stripped else 3.0, center_he4=0.0,
        surface_h1=0.0 if stripped else 0.7,
        log_total_angular_momentum=50.0, spin=0.3, profile=None)
    properties.update(kwargs)
    return SingleStar(**properties)


def make_binary(m_CO, event="CC1", companion=("H-rich_Core_H_burning", 10.0),
                separation=100.0, eccentricity=0.0, **kwargs):
    """A detached binary in which one star is about to collapse."""
    dying_star = make_star(m_CO, **kwargs)
    other_star = SingleStar(state=companion[0], mass=companion[1],
                            metallicity=1.0, log_R=0.5, spin=0.1,
                            profile=None)
    if event == "CC2":
        star_1, star_2 = other_star, dying_star
    else:
        star_1, star_2 = dying_star, other_star
    binary = BinaryStar(star_1=star_1, star_2=star_2, state="detached",
                        event=event, time=1.0e7, separation=separation,
                        eccentricity=eccentricity, V_sys=[0.0, 0.0, 0.0])
    binary.orbital_period = orbital_period_from_separation(
        separation, star_1.mass, star_2.mass)
    return binary


def make_step(**kwargs):
    """A StepSN computing every collapse, with a seeded RNG."""
    settings = dict(use_interp_values=False, RNG=np.random.default_rng(42))
    settings.update(kwargs)
    return totest.StepSN(**settings)


def m_grav(step, m_rembar):
    """Gravitational mass of a NS/BH with baryonic mass m_rembar."""
    return max((20.0 / 3.0) * (np.sqrt(1.0 + 0.3 * m_rembar) - 1.0),
               m_rembar - step.max_neutrino_mass_loss)


class TestPISNPrescription:
    @fixture
    def SN_step(self):
        # Hendriks+23 without extra mass loss, and direct collapse of the He
        # core, so the remnant is the mass left after the pulses
        return totest.StepSN(PISN='Hendriks+23', PPI_extra_mass_loss=0.0,
                             mechanism='direct', use_interp_values=False)

    def get_m_PISN(self, step, Z_div_Zsun):
        star = SingleStar(mass=70.0, he_core_mass=60.0, co_core_mass=50.0)
        star.metallicity = Z_div_Zsun
        return step._get_remnant_mass._mass_after_pulses(star, None)

    def test_Hendriks23_metallicity(self, SN_step):
        # star.metallicity is Z/Zsun, while the fit of Hendriks et al. 2023
        # uses the absolute Z
        Z = 0.001
        delta_M_PPI = ((0.0006 * totest.np.log10(Z) + 0.0054) * (50.0-34.8)**3
                       - 0.0013 * (50.0-34.8)**2)
        assert self.get_m_PISN(SN_step, Z/Zsun) == approx(60.0 - delta_M_PPI)

    def test_Hendriks23_metallicity_floor(self, SN_step):
        # the fit is limited to the lowest absolute Z = 1e-4 in Hendriks+23
        m_PISN_floor = self.get_m_PISN(SN_step, 1e-4/Zsun)
        assert self.get_m_PISN(SN_step, 1e-6) == approx(m_PISN_floor)
        assert self.get_m_PISN(SN_step, 2e-4/Zsun) < m_PISN_floor


class TestSetup:
    def test_defaults(self):
        step = totest.StepSN()
        assert step.mechanism == "Fryer+12-delayed"
        assert isinstance(step.RNG, np.random.Generator)
        # the default settings are those of the first SN model
        assert step._interp_collapse.SN_MODEL_NAME == "SN_MODEL_v2_01"

    def test_invalid_argument(self):
        with raises(ValueError, match="Invalid arguments passed to StepSN"):
            totest.StepSN(not_an_option=1.0)

    def test_no_matching_SN_model(self):
        # use_interp_values=True needs an SN model with the same settings
        with raises(ValueError, match="no SN model"):
            totest.StepSN(PISN_CO_shift=123.0)

    @mark.parametrize("kwargs, message", [
        ({"mechanism": "unknown"}, "given to mechanism"),
        ({"mechanism": "Sukhbold+16-engine", "engine": "bad"},
         "Invalid engine"),
        ({"mechanism": "Patton&Sukhbold20-engine", "engine": "bad"},
         "Invalid engine"),
        ({"mechanism": "Couch+20-engine", "engine": "bad"},
         "Invalid turbulence strength"),
        ({"ECSN": "bad"}, "given to ECSN prescription"),
        ({"PISN": "bad"}, "given to PISN prescription"),
        ({"kick_prescription": "bad"}, "Invalid kick_prescription"),
        ({"kick_normalisation": "bad"}, "Invalid kick_normalisation"),
    ])
    def test_invalid_options(self, kwargs, message):
        with raises(ValueError, match=message):
            make_step(**kwargs)

    def test_deprecated_direct_he_core(self):
        with warns(pw.DeprecationWarning, match="direct_he_core"):
            step = make_step(mechanism="direct_he_core")
        assert step.mechanism == "direct"
        assert step.conserve_hydrogen_envelope is False
        with raises(ValueError, match="cannot be combined"):
            make_step(mechanism="direct_he_core",
                      conserve_hydrogen_envelope=True)

    def test_unavailable_Maltsev25(self):
        with warns(pw.UnsupportedModelWarning, match="not available"):
            step = make_step(mechanism="Maltsev+25-engine", engine="")
        assert step.mechanism == "Patton&Sukhbold20-engine"
        assert step.engine == "N20"
        # the factories do not build it either
        for get_prescription in (get_CCSN_prescription,
                                 get_remnant_mass_prescription):
            with raises(NotImplementedError, match="not available"):
                get_prescription(mechanism="Maltsev+25-engine")

    @mark.parametrize("option", ["asym_ej", "linear"])
    def test_deprecated_kick_normalisation(self, option):
        with warns(pw.DeprecationWarning, match="kick_prescription"):
            step = make_step(kick_normalisation=option)
        assert step.kick_prescription == option
        assert step.kick_normalisation == "one"


class TestSNType:
    @mark.parametrize("ECSN, PISN, m_CO, SN_type", [
        ("Tauris+15", "Hendriks+23", 1.0, "WD"),
        ("Tauris+15", "Hendriks+23", 1.4, "ECSN"),
        ("Tauris+15", "Hendriks+23", 1.5, "CCSN"),
        ("Podsiadlowski+04", "Hendriks+23", 0.5, "WD"),     # He core 1.025
        ("Podsiadlowski+04", "Hendriks+23", 1.2, "ECSN"),    # He core 1.9
        ("Podsiadlowski+04", "Hendriks+23", 2.0, "CCSN"),    # He core 2.9
        ("No_ECSN", "Hendriks+23", 1.0, "WD"),
        ("No_ECSN", "Hendriks+23", 1.4, "CCSN"),
        ("Tauris+15", "Hendriks+23", 50.0, "PPISN"),
        ("Tauris+15", "Hendriks+23", 120.0, "CCSN"),         # above the gap
        ("Tauris+15", "Marchant+19", 40.0, "PPISN"),         # He core 50.4
        ("Tauris+15", "Marchant+19", 60.0, "PISN"),          # He core 75.4
        ("Tauris+15", "Marchant+19", 100.0, "CCSN"),         # above the gap
        ("Tauris+15", 45.0, 40.0, "PPISN"),                  # He core 50.4
        ("Tauris+15", 45.0, 30.0, "CCSN"),                   # He core 37.9
        ("Tauris+15", None, 60.0, "CCSN"),
    ])
    def test_collapse_outcome(self, ECSN, PISN, m_CO, SN_type):
        step = make_step(ECSN=ECSN, PISN=PISN)
        assert step.collapse_outcome(make_star(m_CO)) == SN_type

    @mark.parametrize("PISN, m_CO, explosion", [
        ("Hendriks+23", 5.0, None),         # Fryer+12: the mass decides
        ("Hendriks+23", 1.0, None),         # WD
        ("Hendriks+23", 120.0, "BH"),       # above the gap
        ("Marchant+19", 100.0, "BH"),       # above the gap
        (45.0, 100.0, None),                # no gap with a constant cap
    ])
    def test_explosion_outcome(self, PISN, m_CO, explosion):
        step = make_step(PISN=PISN)
        star = make_star(m_CO)
        SN_type = step.collapse_outcome(star)
        assert step.explosion_outcome(star, SN_type) == explosion

    @mark.parametrize("mechanism, engine, m_CO, explosion", [
        ("Sukhbold+16-engine", "N20", 2.0, "NS"),
        ("Sukhbold+16-engine", "N20", 20.0, "BH"),
        ("Couch+20-engine", "1.0", 5.0, "BH"),
        ("Couch+20-engine", "1.3", 20.0, "NS"),     # all models explode
    ])
    def test_explosion_criteria_from_tables(self, mechanism, engine, m_CO,
                                            explosion):
        step = make_step(mechanism=mechanism, engine=engine)
        assert step.explosion_outcome(make_star(m_CO), "CCSN") == explosion

    @mark.parametrize("M_CO, X_C", [(2.0, 0.3), (3.0, 0.3), (5.0, 0.2),
                                    (7.0, 0.1), (9.0, 0.3), (12.0, 0.3)])
    def test_explosion_criterion_Patton20(self, M_CO, X_C):
        step = make_step(mechanism="Patton&Sukhbold20-engine", engine="N20")
        star = make_star(5.0, avg_c_in_c_core_at_He_depletion=X_C,
                         co_core_mass_at_He_depletion=M_CO)
        explosion = step.explosion_outcome(star, "CCSN")
        # M4 and mu4 of the matching model are recorded on the star
        assert isinstance(star.M4, float) and isinstance(star.mu4, float)
        if M_CO <= 2.5:
            assert explosion == "NS"
        elif M_CO >= 10.0:
            assert explosion == "BH"
        else:
            explodes = explodes_Ertl16(star.M4, star.mu4, 0.194, 0.0580)
            assert explosion == ("NS" if explodes else "BH")

    def test_Patton20_without_CO_core_at_He_depletion(self):
        step = make_step(mechanism="Patton&Sukhbold20-engine", engine="N20")
        with raises(ModelError, match="no CO core at He depletion"):
            step.explosion_outcome(make_star(5.0), "CCSN")
        # unless it is approximated from the final CO core
        step = make_step(mechanism="Patton&Sukhbold20-engine", engine="N20",
                         approx_at_he_depletion=True)
        assert step.explosion_outcome(make_star(2.0), "CCSN") == "NS"


class TestRemnant:
    def test_WD(self):
        # WDs keep their radius
        step = make_step(use_core_masses=True)
        star = make_star(1.0)
        assert step.collapse_star(star) is None
        assert (star.SN_type, star.state, star.f_fb) == ("WD", "WD", 1.0)
        assert star.mass == approx(1.0)
        assert star.co_core_mass == approx(1.0)
        assert star.log_R == approx(np.log10(CO_radius(1.0, "WD")))

    @mark.parametrize("ECSN, m_rembar", [("Tauris+15", 1.4),
                                         ("Podsiadlowski+04", 1.38)])
    def test_ECSN(self, ECSN, m_rembar):
        # ECSNe also work without profile and with use_core_masses=False
        step = make_step(ECSN=ECSN, use_core_masses=False)
        star = make_star(1.4 if ECSN == "Tauris+15" else 1.2)
        assert step.collapse_star(star) is None
        assert (star.SN_type, star.state, star.f_fb) == ("ECSN", "NS", 0.0)
        assert star.mass == approx(m_grav(step, m_rembar))
        assert star.spin == 0.0

    def test_Fryer12_delayed_NS(self):
        step = make_step(mechanism="Fryer+12-delayed")
        star = make_star(2.0)           # He core 2.9
        assert step.collapse_star(star) is None
        # m_proto = 1.2 and m_fb = 0.2 for m_CO < 2.5
        assert (star.SN_type, star.state) == ("CCSN", "NS")
        assert star.f_fb == approx(0.2 / (2.9 - 1.2))
        assert star.mass == approx(m_grav(step, 1.4))
        assert star.log_R == approx(np.log10(CO_radius(star.mass, "NS")))

    @mark.parametrize("mechanism", ["Fryer+12-delayed", "Fryer+12-rapid"])
    @mark.parametrize("m_CO", [6.0, 7.0, 11.0, 15.0])
    def test_Fryer12_full_fallback(self, mechanism, m_CO):
        step = make_step(mechanism=mechanism)
        star = make_star(m_CO)
        m_He = star.he_core_mass
        assert step.collapse_star(star) is None
        if m_CO >= 11.0 or (mechanism == "Fryer+12-rapid" and m_CO < 7.0):
            assert star.f_fb == approx(1.0)
            assert star.mass == approx(m_grav(step, m_He))
        assert star.state == "BH"
        assert star.metallicity == 1.0
        assert star.co_core_mass == m_CO

    @mark.parametrize("state, conserve_hydrogen_envelope", [
        ("stripped_He_Core_C_depleted", False),
        ("H-rich_Core_C_depleted", False),
        ("H-rich_Core_C_depleted", True),
    ])
    def test_direct(self, state, conserve_hydrogen_envelope):
        step = make_step(mechanism="direct",
                         conserve_hydrogen_envelope=conserve_hydrogen_envelope)
        star = make_star(5.0, state=state)
        m_collapsing = (star.mass if conserve_hydrogen_envelope
                        else star.he_core_mass)
        assert step.collapse_star(star) is None
        assert (star.SN_type, star.state, star.f_fb) == ("CCSN", "BH", 1.0)
        assert star.mass == approx(m_grav(step, m_collapsing))

    @mark.parametrize("mechanism, engine", [
        ("Fryer+12-delayed", ""), ("Sukhbold+16-engine", "N20"),
        ("Couch+20-engine", "1.3"), ("Patton&Sukhbold20-engine", "N20")])
    def test_direct_collapse_above_the_gap(self, mechanism, engine):
        step = make_step(mechanism=mechanism, engine=engine)
        star = make_star(120.0)
        m_He = star.he_core_mass
        assert step.collapse_star(star) is None
        assert (star.SN_type, star.state, star.f_fb) == ("CCSN", "BH", 1.0)
        assert star.mass == approx(m_grav(step, m_He))

    def test_Sukhbold16_NS(self):
        step = make_step(mechanism="Sukhbold+16-engine", engine="N20")
        star = make_star(2.0)
        table = Sukhbold16_table("N20")
        NS = table[table["stellar_state"] == 13]
        m_rembar = np.interp(star.he_core_mass, NS["He_c_mass"],
                             NS["Rem_mass"])
        assert step.collapse_star(star) is None
        assert (star.SN_type, star.state, star.f_fb) == ("CCSN", "NS", 0.0)
        assert star.mass == approx(m_grav(step, m_rembar))

    def test_Couch20_NS(self):
        step = make_step(mechanism="Couch+20-engine", engine="1.3")
        star = make_star(5.0)
        assert step.collapse_star(star) is None
        assert (star.SN_type, star.state, star.f_fb) == ("CCSN", "NS", 0.0)
        assert star.mass == approx(m_grav(step, 1.4))

    @mark.parametrize("M_CO, X_C", [(2.0, 0.3), (3.0, 0.3), (5.0, 0.2)])
    def test_Patton20_NS(self, M_CO, X_C):
        step = make_step(mechanism="Patton&Sukhbold20-engine", engine="N20")
        star = make_star(5.0, avg_c_in_c_core_at_He_depletion=X_C,
                         co_core_mass_at_He_depletion=M_CO)
        assert step.collapse_star(star) is None
        assert (star.SN_type, star.state, star.f_fb) == ("CCSN", "NS", 0.0)
        # M4 of the matching model, or 1.25 Msun below the models
        m_rembar = 1.25 if M_CO <= 2.5 else star.M4
        assert star.mass == approx(m_grav(step, m_rembar))

    def test_PPISN_Marchant19(self):
        step = make_step(PISN="Marchant+19")
        star = make_star(40.0)          # He core 50.4
        m_PPI = PISN_mass_Marchant(star)
        assert m_PPI < star.he_core_mass
        assert step.collapse_star(star) is None
        assert (star.SN_type, star.state) == ("PPISN", "BH")
        assert star.mass == approx(m_grav(step, m_PPI))

    def test_PPISN_constant_cap(self):
        step = make_step(PISN=45.0)
        star = make_star(40.0)          # He core 50.4
        assert step.collapse_star(star) is None
        assert (star.SN_type, star.state) == ("PPISN", "BH")
        assert star.mass == approx(m_grav(step, 45.0))

    def test_PPISN_Hendriks23(self):
        # without extra mass loss the pulses remove mass, ...
        step = make_step(PISN="Hendriks+23", PPI_extra_mass_loss=0.0)
        star = make_star(50.0)          # He core 62.9
        assert step.collapse_star(star) is None
        assert (star.SN_type, star.state) == ("PPISN", "BH")
        assert star.mass < m_grav(step, 62.9)
        # ... while the default extra mass loss of -20 Msun leaves no pulses
        step = make_step(PISN="Hendriks+23")
        star = make_star(50.0)
        assert step.collapse_star(star) is None
        assert (star.SN_type, star.state) == ("CCSN", "BH")
        assert star.mass == approx(m_grav(step, 62.9))

    @mark.parametrize("PISN, m_CO", [("Marchant+19", 60.0),
                                     ("Hendriks+23", 100.0)])
    def test_PISN(self, PISN, m_CO):
        step = make_step(PISN=PISN)
        star = make_star(m_CO)
        assert step.collapse_star(star) is None
        assert (star.SN_type, star.state, star.mass) == (
            "PISN", "massless_remnant", 0.0)

    def test_spin(self):
        # use_core_masses: BH spin from the angular momentum, capped at 1
        step = make_step(mechanism="direct")
        star = make_star(10.0, log_total_angular_momentum=60.0)
        assert step.collapse_star(star) is None
        assert star.spin == 1.0
        # allow_spin_None: no spin
        step = make_step(mechanism="direct", use_core_masses=False,
                         allow_spin_None=True)
        star = make_star(10.0)
        assert step.collapse_star(star) is None
        assert star.state == "BH" and star.spin is None

    def test_NS_or_BH_decided_by_the_mass(self):
        step = make_step()
        with warns(pw.ReplaceValueWarning, match="makes it a BH"):
            assert step._update_params._NS_or_BH("NS", 3.0) == "BH"
        assert step._update_params._NS_or_BH(None, step.max_NS_mass) == "BH"
        assert step._update_params._NS_or_BH(None, 2.0) == "NS"

    def test_failed_collapse(self):
        # no way to set the spin
        step = make_step(use_core_masses=False)
        assert step.collapse_star(make_star(5.0)) == "FAILED core collapse!"
        binary = make_binary(5.0)
        with raises(ModelError, match="FAILED core collapse!"):
            step(binary)
        assert binary.state == "ERR"

    def test_star_cannot_collapse(self):
        step = make_step()
        star = make_star(5.0, state="H-rich_Core_H_burning")
        assert "cannot collapse" in step.collapse_star(star)

    @mark.parametrize("use_interp_values", [False, True])
    @mark.parametrize("m_CO", [None, np.nan])
    def test_CO_core_not_set(self, use_interp_values, m_CO):
        step = make_step(use_interp_values=use_interp_values)
        star = make_star(5.0, co_core_mass=m_CO)
        star.SN_MODEL_v2_01 = TestInterpolatedPath.outcome()
        with raises(ValueError, match="The CO core mass is not correct"):
            step.collapse_star(star)


class TestInterpolatedPath:
    @fixture
    def step(self):
        # the default settings match SN_MODEL_v2_01
        return totest.StepSN(use_interp_values=True,
                             RNG=np.random.default_rng(42))

    @staticmethod
    def outcome(**kwargs):
        values = dict(state="BH", SN_type="CCSN", f_fb=1.0, mass=9.0,
                      spin=0.2, m_disk_accreted=0.0, m_disk_radiated=0.0,
                      M4=None, mu4=None, h1_mass_ej=0.0, he4_mass_ej=0.0)
        values.update(kwargs)
        return values

    def test_used(self, step):
        star = make_star(10.0)
        star.SN_MODEL_v2_01 = self.outcome()
        assert step.collapse_star(star) is None
        assert (star.SN_type, star.state, star.mass, star.spin) == (
            "CCSN", "BH", 9.0, 0.2)
        assert star.log_R == approx(np.log10(CO_radius(9.0, "BH")))
        # metallicity and CO core are kept, as in the computed path
        assert star.metallicity == 1.0
        assert star.co_core_mass == 10.0

    @mark.parametrize("SN_type, state, mass", [
        ("PPISN", "BH", 40.0), ("CCSN", "BH", 60.0),
        ("PISN", "massless_remnant", 0.0)])
    def test_used_for_PPISN_candidates(self, step, SN_type, state, mass):
        star = make_star(50.0)              # PPISN candidate
        star.SN_MODEL_v2_01 = self.outcome(SN_type=SN_type, state=state,
                                           mass=mass)
        assert step.collapse_star(star) is None
        assert (star.SN_type, star.state, star.mass) == (SN_type, state,
                                                         mass)

    def test_WD(self, step):
        star = make_star(1.0)
        star.SN_MODEL_v2_01 = self.outcome(SN_type="WD", state="WD",
                                           mass=0.9, spin=0.0)
        assert step.collapse_star(star) is None
        assert (star.state, star.mass) == ("WD", 0.9)
        assert star.co_core_mass == 0.9 and star.he_core_mass == 0.9

    def test_ECSN_computed(self, step):
        # ECSNe are decided by the ECSN option of the step
        star = make_star(1.4)
        star.SN_MODEL_v2_01 = self.outcome(SN_type="CCSN", state="NS",
                                           mass=1.5)
        assert step.collapse_star(star) is None
        assert (star.SN_type, star.state) == ("ECSN", "NS")
        assert star.mass == approx(m_grav(step, 1.4))

    @mark.parametrize("values, warning, message", [
        (None, pw.InterpolationWarning, "not interpolated"),
        ({"SN_type": "ECSN", "state": "NS", "mass": 1.3},
         pw.ApproximationWarning, "does not match the star"),
        ({"SN_type": "CCSN", "state": "WD"},
         pw.ApproximationWarning, "does not match the predicted CO"),
        ({"mass": np.nan}, pw.ApproximationWarning,
         "no interpolated remnant mass"),
    ])
    def test_computed_instead(self, step, values, warning, message):
        star = make_star(5.0)
        star.SN_MODEL_v2_01 = (None if values is None
                               else self.outcome(**values))
        with warns(warning, match=message):
            assert step.collapse_star(star) is None
        # the same as without interpolated values
        computed = make_star(5.0)
        assert make_step().collapse_star(computed) is None
        assert (star.SN_type, star.state, star.f_fb) == (
            computed.SN_type, computed.state, computed.f_fb)
        assert star.mass == approx(computed.mass)
        # the stored outcome is not changed
        if values is not None:
            assert star.SN_MODEL_v2_01 == self.outcome(**values)


class TestKicks:
    def test_kick(self):
        binary = make_binary(5.0)
        make_step()(binary)
        star = binary.star_1
        assert star.natal_kick_velocity > 0.0
        for angle in (star.natal_kick_azimuthal_angle,
                      star.natal_kick_polar_angle,
                      star.natal_kick_mean_anomaly):
            assert np.isfinite(angle)

    def test_no_kick(self):
        # the SN mass loss still changes the orbit
        binary = make_binary(5.0)
        make_step(kick=False)(binary)
        assert binary.star_1.natal_kick_velocity == 0.0
        assert np.isfinite(binary.star_1.natal_kick_polar_angle)
        assert binary.state == "detached"
        assert binary.separation > 100.0
        assert binary.eccentricity > 0.0

    def test_preset_values(self):
        binary = make_binary(5.0)
        star = binary.star_1
        star.natal_kick_velocity = 100.0
        star.natal_kick_azimuthal_angle = 1.0
        star.natal_kick_polar_angle = 0.5
        star.natal_kick_mean_anomaly = 2.0
        make_step()(binary)
        assert (star.natal_kick_velocity, star.natal_kick_azimuthal_angle,
                star.natal_kick_polar_angle,
                star.natal_kick_mean_anomaly) == (100.0, 1.0, 0.5, 2.0)

    def test_NaN_presets_are_drawn(self):
        presets = dict(natal_kick_velocity=np.nan,
                       natal_kick_azimuthal_angle=np.nan,
                       natal_kick_polar_angle=np.nan,
                       natal_kick_mean_anomaly=np.nan)
        binary_NaN, binary_unset = make_binary(5.0, **presets), make_binary(5.0)
        make_step()(binary_NaN)
        make_step()(binary_unset)
        for key in presets:
            assert getattr(binary_NaN.star_1, key) == approx(
                getattr(binary_unset.star_1, key))

    def test_WD(self):
        binary = make_binary(1.0)
        make_step()(binary)
        assert binary.star_1.state == "WD"
        assert binary.star_1.natal_kick_velocity == 0.0
        # gradual mass loss: wider orbit, same eccentricity
        assert binary.separation > 100.0
        assert binary.eccentricity == 0.0

    @mark.parametrize("kwargs", [
        {"mechanism": "direct", "no_kick_direct_collapse": True},
        {"mechanism": "direct", "sigma_kick_CCSN_BH": None},
        {"kick_normalisation": "zero"},
    ])
    def test_zero_kick(self, kwargs):
        binary = make_binary(10.0)
        make_step(**kwargs)(binary)
        assert binary.star_1.natal_kick_velocity == 0.0

    @mark.parametrize("prescription", ["maxwellian", "log_normal",
                                       "asym_ej", "linear"])
    def test_kick_prescriptions(self, prescription):
        binary = make_binary(2.0)
        make_step(kick_prescription=prescription)(binary)
        assert np.isfinite(binary.star_1.natal_kick_velocity)


class TestStep:
    def test_invalid_event(self):
        binary = make_binary(5.0, event=None)
        with raises(ValueError, match="Invalid call of step_SN"):
            make_step()(binary)

    def test_CC2(self):
        binary = make_binary(5.0, event="CC2")
        make_step()(binary)
        assert binary.star_2.state == "NS"
        # star_1 does not collapse
        assert binary.star_1.mass == 10.0
        assert binary.event is None

    def test_disrupted_by_PISN(self):
        binary = make_binary(60.0)
        make_step(PISN="Marchant+19")(binary)
        assert binary.star_1.state == "massless_remnant"
        assert binary.state == "disrupted"
        assert np.isnan(binary.inspiral_time)
        assert np.isnan(binary.star_1.spin_orbit_tilt_first_SN)

    def test_double_compact_object(self):
        binary = make_binary(2.0, companion=("NS", 1.4))
        make_step(kick=False)(binary)
        assert binary.state == "detached"
        assert np.isfinite(binary.inspiral_time)
        assert np.isfinite(binary.star_1.spin_orbit_tilt_first_SN)
        assert binary.first_SN_already_occurred

    @mark.parametrize("states, event", [
        (("NS", "stripped_He_Core_C_depleted"), "CC2"),
        (("H-rich_Core_C_depleted", "BH"), "CC1"),
        (("NS", "H-rich_Core_H_burning"), None),
    ])
    def test_companion_collapse_follows(self, states, event):
        binary = make_binary(2.0, event=None)
        binary.star_1.state, binary.star_2.state = states
        make_step().finalise_step(binary)
        assert binary.event == event
