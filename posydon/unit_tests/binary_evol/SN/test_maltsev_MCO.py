"""Unit tests for the Maltsev+25 M_CO (rapid) supernova prescription."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import numpy as np
from pytest import approx, fixture, raises

from posydon.binary_evol.SN.maltsev_MCO import (
    Maltsev25_MCO_corecollapse,
)
from posydon.binary_evol.SN.step_SN import StepSN


class FakeStar:
    """Minimal stand-in for a collapsing SingleStar."""

    def __init__(self, co_core_mass, metallicity, he_core_mass=5.0, mass=20.0):
        self.co_core_mass = co_core_mass
        self.metallicity = metallicity  # Z / Z_sun
        self.he_core_mass = he_core_mass
        self.mass = mass
        self.SN_categorisation = None


@fixture
def engine():
    return Maltsev25_MCO_corecollapse(RNG=np.random.default_rng(2025))


def test_registered_mechanism():
    s = StepSN(mechanism="Maltsev+25-MCO-rapid", verbose=False)
    assert s.Maltsev25_MCO_rapid in s.mechanisms
    assert hasattr(s, "Maltsev25_MCO_engine")


def test_boundaries_at_Zsun(engine):
    # single @ Z_sun: M1=6.6, M2=7.2, M3=13.0
    M1, M2, M3 = engine.get_boundaries("single", 1.0)
    assert M1 == approx(6.6)
    assert M2 == approx(7.2)
    assert M3 == approx(13.0)
    # case_A @ Z_sun
    M1, M2, M3 = engine.get_boundaries("case_A", 1.0)
    assert (M1, M2, M3) == approx((7.4, 8.4, 15.4))
    # case_B @ Z_sun (Be + Bl merged)
    M1, M2, M3 = engine.get_boundaries("case_B", 1.0)
    assert (M1, M2, M3) == approx((7.7, 8.3, 15.2))
    # case_C @ Z_sun
    M1, M2, M3 = engine.get_boundaries("case_C", 1.0)
    assert (M1, M2, M3) == approx((6.6, 7.1, 13.2))


def test_explodability_deterministic(engine):
    # M_CO < M1 -> successful
    assert engine.explodability(5.0, 1.0, "single") is True
    # M1 <= M_CO <= M2 -> failed (direct BH)
    assert engine.explodability(7.0, 1.0, "single") is False
    # M2 < M_CO < M3 -> successful
    assert engine.explodability(10.0, 1.0, "single") is True
    # M_CO >= M3 -> failed (direct BH)
    assert engine.explodability(14.0, 1.0, "single") is False


def test_explodability_Z_dependence(engine):
    # At Z/Z_sun = 0.1 the single boundaries shift to 6.1, 6.6, 12.9
    M1, M2, M3 = engine.get_boundaries("single", 0.1)
    assert (M1, M2, M3) == approx((6.1, 6.6, 12.9))
    # A M_CO that explodes at solar may fail at sub-solar (window shifts down)
    assert engine.explodability(6.4, 1.0, "single") is True
    assert engine.explodability(6.4, 0.1, "single") is False


def test_categorisation_guaranteed_NS_window(engine):
    # Inside the guaranteed-NS window (9.0, 10.2) for single @ Z_sun
    assert engine.categorisation(9.5, 1.0, "single") == "NS"
    assert engine.categorisation(10.0, 1.0, "single") == "NS"
    # The window is open: exactly at the lower boundary the outcome is not
    # guaranteed to be an NS (it follows the stochastic split).
    NS1, NS2 = engine.get_NS_window("single", 1.0)
    assert NS1 < 9.5 < NS2
    assert 9.0 == NS1  # boundary value is excluded from the guaranteed window


def test_categorisation_probabilistic_split():
    # Use a fixed RNG and count the NS / fallback_BH split (model A: 15%)
    rng = np.random.default_rng(7)
    eng = Maltsev25_MCO_corecollapse(RNG=rng, fallback_model="A")
    counts = {"NS": 0, "fallback_BH": 0}
    # M_CO=12.0 (in (M2,M3) but outside NS window) for many draws
    for _ in range(20000):
        counts[eng.categorisation(12.0, 1.0, "single")] += 1
    frac_fb = counts["fallback_BH"] / sum(counts.values())
    assert frac_fb == approx(0.15, abs=0.02)


def test_low_MCO_successful_SN_is_always_NS():
    # Regression: a star with M_CO < M1 is a successful SN that forms an NS
    # only (the paper treats fallback-BH as statistically insignificant there,
    # line 571). The stochastic NS/fallback-BH split must NOT apply, so this
    # must hold for every RNG draw.
    rng = np.random.default_rng(12345)
    eng = Maltsev25_MCO_corecollapse(RNG=rng, fallback_model="A")
    # single @ Z_sun: M1 = 6.6, so 6.457 < M1 -> successful, must be NS
    for _ in range(5000):
        assert eng.categorisation(6.457, 1.0, "single") == "NS"
    # same via __call__
    for _ in range(5000):
        star = FakeStar(co_core_mass=6.457, metallicity=1.0)
        _, _, state = eng(star, "single")
        assert state == "NS"
        assert star.SN_categorisation == "NS"


def test_fallback_model_B():
    eng = Maltsev25_MCO_corecollapse(
        RNG=np.random.default_rng(3), fallback_model="B")
    assert eng.fallback_model == "B"
    counts = {"NS": 0, "fallback_BH": 0}
    for _ in range(20000):
        counts[eng.categorisation(12.0, 1.0, "single")] += 1
    frac_fb = counts["fallback_BH"] / sum(counts.values())
    assert frac_fb == approx(0.10, abs=0.02)


def test_call_direct_BH(engine):
    star = FakeStar(co_core_mass=7.0, metallicity=1.0)
    m_rem, f_fb, state = engine(star, "single", conserve_hydrogen_envelope=False)
    assert state == "BH"
    assert f_fb == 1.0
    assert star.SN_categorisation == "direct_BH"
    assert m_rem == approx(star.he_core_mass)


def test_call_NS(engine):
    star = FakeStar(co_core_mass=5.0, metallicity=1.0)
    m_rem, f_fb, state = engine(star, "single", conserve_hydrogen_envelope=False)
    assert state == "NS"
    assert f_fb == 0.0
    assert star.SN_categorisation == "NS"
    assert m_rem == approx(engine.NS_mass)


def test_call_fallback_BH(engine):
    # Force the fallback-BH outcome deterministically via a stub RNG that always
    # draws below the 0.15 fallback probability (model A).
    class _StubRNG:
        def uniform(self, a, b):
            return 0.01
    eng = Maltsev25_MCO_corecollapse(RNG=_StubRNG(), fallback_model="A")
    star = FakeStar(co_core_mass=12.0, metallicity=1.0)
    m_rem, f_fb, state = eng(star, "single", conserve_hydrogen_envelope=False)
    assert state == "BH"
    assert star.SN_categorisation == "fallback_BH"
    assert f_fb == approx(eng.fallback_fraction)
    assert m_rem == approx(star.he_core_mass)


def test_custom_NS_mass_model():
    def my_ns_mass(star):
        return 2.0
    eng = Maltsev25_MCO_corecollapse(NS_mass_model=my_ns_mass,
                                     RNG=np.random.default_rng(0))
    star = FakeStar(co_core_mass=5.0, metallicity=1.0)
    m_rem, _, _ = eng(star, "single")
    assert m_rem == approx(2.0)


def test_out_of_range_Z_extrapolates(engine):
    # Z/Z_sun > 1 is outside the calibrated range (0.1, 1); the boundaries
    # are linearly extrapolated, so the engine completes without warning.
    star = FakeStar(co_core_mass=5.0, metallicity=2.0)  # > Z_sun
    m_rem, f_fb, state = engine(star, "single")
    assert state in ("NS", "BH")


def test_missing_MCO_raises():
    eng = Maltsev25_MCO_corecollapse(RNG=np.random.default_rng(0))
    star = FakeStar(co_core_mass=np.nan, metallicity=1.0)
    from posydon.utils.posydonerror import ModelError
    with raises(ModelError):
        eng(star, "single")


def test_resolve_mt_class_reads_star_mt_class(engine):
    # The MT class of the first mass-transfer episode (whichever star is the
    # donor) is stored on the star during step_MESA; the engine returns it
    # as-is.
    star = FakeStar(co_core_mass=2.0, metallicity=1.0)
    star.mt_class = "case_B"
    assert engine._resolve_mt_class(star) == "case_B"
    star.mt_class = "case_A"
    assert engine._resolve_mt_class(star) == "case_A"


def test_resolve_mt_class_defaults_to_single(engine):
    # No class is stored on the star (e.g. single-star path, or an
    # interpolation trained without first_mt_case), or the stored value is not
    # a valid MT class (e.g. the verbatim termination flag token 'no_RLOF'),
    # so the recipe falls back to 'single' rather than a guessed class.
    star = FakeStar(co_core_mass=2.0, metallicity=1.0)
    assert not hasattr(star, "mt_class")
    assert engine._resolve_mt_class(star) == "single"
    # no collapsing star (single-star path) -> single
    assert engine._resolve_mt_class(None) == "single"
    # invalid / retained-verbatim flags normalize to 'single'
    for mt_class in ["no_RLOF", "initial_RLOF", "not_converged", "None"]:
        star.mt_class = mt_class
        assert engine._resolve_mt_class(star) == "single"
    # an explicit fallback class is honoured when no valid class is stored
    assert engine._resolve_mt_class(star, "case_C") == "case_C"


def test_resolve_mt_class_first_episode_overall(engine):
    # The recipe uses the class of the first mass-transfer episode, whichever
    # star was the donor. step_MESA stores that class on both stars, so either
    # collapsing star resolves to the same class.
    from posydon.utils.common_functions import mt_class_from_cumulative
    star_1 = FakeStar(co_core_mass=2.0, metallicity=1.0)
    star_2 = FakeStar(co_core_mass=2.0, metallicity=1.0)
    star_1.mt_class = star_2.mt_class = mt_class_from_cumulative("case_A1/B2")
    assert star_1.mt_class == "case_A"
    assert engine._resolve_mt_class(star_1) == "case_A"
    assert engine._resolve_mt_class(star_2) == "case_A"
    # ...even when the collapsing star was not the first donor
    star_1.mt_class = star_2.mt_class = mt_class_from_cumulative("case_B2/A1")
    assert engine._resolve_mt_class(star_2) == "case_B"


def test_resolve_mt_class_latest_grid_wins(engine):
    # Each step_MESA run overwrites star.mt_class with the class resolved from
    # that grid's first_mt_case, so the most recent grid wins (e.g. a
    # CO_HMS_RLO / CO_HeMS run takes precedence over the earlier HMS_HMS one).
    from posydon.utils.common_functions import mt_class_from_cumulative
    star = FakeStar(co_core_mass=2.0, metallicity=1.0)
    # HMS_HMS grid first (case A donation)...
    star.mt_class = mt_class_from_cumulative("case_A1/B1")
    assert star.mt_class == "case_A"
    # ...then a later CO_HMS_RLO grid overwrites with its own class. The class
    # is the first episode overall, so the second grid's case_B wins.
    star.mt_class = mt_class_from_cumulative("case_B2/A1")
    assert engine._resolve_mt_class(star) == "case_B"
    # the latest grid always wins, even if it reports no RLO at all
    star.mt_class = mt_class_from_cumulative("no_RLO")
    assert star.mt_class == "no_RLO"
    assert engine._resolve_mt_class(star) == "single"


def test_call_auto_resolves_star_mt_class(engine):
    # M_CO = 7.0 sits in the 'single' direct-collapse window (M1-M2 = 6.6-7.2)
    # but below the case_A M1 = 7.4, so the star's MT class (set by step_MESA)
    # fully determines the outcome.
    star = FakeStar(co_core_mass=7.0, metallicity=1.0)
    star.mt_class = "case_A"
    _, _, state = engine(star, conserve_hydrogen_envelope=False)
    assert state == "NS"
    # an invalid/verbatim stored class falls back to 'single' boundaries
    star.mt_class = "no_RLOF"
    _, _, state = engine(star, conserve_hydrogen_envelope=False)
    assert state == "BH"
    # without a stored class the explicit argument is the fallback
    delattr(star, "mt_class")
    _, _, state = engine(star, "case_A", conserve_hydrogen_envelope=False)
    assert state == "NS"
