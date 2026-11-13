"""Unit tests for the first MT case bookkeeping of step_MESA."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

from types import SimpleNamespace

from pytest import mark

from posydon.binary_evol.MESA.step_mesa import (
    FIRST_MT_CASES,
    _set_first_mt_case,
    _update_first_mt_cases,
)
from posydon.binary_evol.SN.maltsev_MCO import Maltsev25_MCO_corecollapse


def new_star():
    # BinaryStar initialises the attribute to None
    return SimpleNamespace(first_mt_case=None)


def mt_class(star):
    return Maltsev25_MCO_corecollapse()._resolve_mt_class(star)


class TestSetFirstMTCase:

    @mark.parametrize("mt_case", FIRST_MT_CASES)
    def test_records_MT_cases(self, mt_case):
        star = new_star()
        _set_first_mt_case(star, mt_case)
        assert star.first_mt_case == mt_case

    @mark.parametrize("mt_case", [None, "no_RLOF", "no_RLO", "initial_RLOF",
                                  "not_converged", "case_nonburning",
                                  "case_undetermined_MT", "case_A2",
                                  "?contact_during_MS"])
    def test_ignores_non_MT_cases(self, mt_case):
        star = new_star()
        _set_first_mt_case(star, mt_case)
        assert star.first_mt_case is None

    def test_keeps_earlier_case(self):
        star = new_star()
        _set_first_mt_case(star, "case_B")
        for mt_case in ["case_BB", "case_C", "no_RLOF", None]:
            _set_first_mt_case(star, mt_case)
        assert star.first_mt_case == "case_B"

    def test_star_without_attribute(self):
        star = SimpleNamespace()
        _set_first_mt_case(star, "case_A")
        assert star.first_mt_case == "case_A"


class TestSequences:
    """Evolutionary sequences through several grids (grid order per step)."""

    def test_accretor_of_disrupted_binary_is_single(self):
        # HMS-HMS case A from star 1; SN1 disrupts the binary and star 2
        # collapses without entering another grid
        s1, s2 = new_star(), new_star()
        _update_first_mt_cases([s1, s2], [False, False], "case_A1")
        assert mt_class(s1) == "case_A"
        assert s2.first_mt_case is None
        assert mt_class(s2) == "single"

    def test_reverse_MT_gives_each_star_its_own_case(self):
        s1, s2 = new_star(), new_star()
        _update_first_mt_cases([s1, s2], [False, False], "case_A1/B2")
        assert mt_class(s1) == "case_A"
        assert mt_class(s2) == "case_B"

    def test_compact_companion_is_not_updated(self):
        star, co = new_star(), new_star()
        _update_first_mt_cases([star, co], [False, True], "case_B1")
        assert star.first_mt_case == "case_B"
        assert co.first_mt_case is None

    @mark.parametrize("tf2_CO_HeMS", ["no_RLOF", "initial_RLOF", "case_BB1",
                                      "case_BA1/BB1", "not_converged"])
    def test_CE_stripped_star_stays_case_B(self, tf2_CO_HeMS):
        # star 2 accretes in HMS-HMS (case_A1), donates in CO-HMS_RLO
        # (unstable case B -> CE) and then evolves in CO-HeMS
        s1, s2 = new_star(), new_star()
        _update_first_mt_cases([s1, s2], [False, False], "case_A1")
        assert s2.first_mt_case is None
        # in the CO grids the non-compact star is star 1 of the grid
        _update_first_mt_cases([s2, s1], [False, True], "case_B1")
        _update_first_mt_cases([s2, s1], [False, True], tf2_CO_HeMS)
        assert mt_class(s2) == "case_B"
        assert mt_class(s1) == "case_A"

    def test_first_donor_episode_in_a_later_grid(self):
        # an accretor in HMS-HMS that donates for the first time in CO-HMS_RLO
        s1, s2 = new_star(), new_star()
        _update_first_mt_cases([s1, s2], [False, False], "case_B1")
        _update_first_mt_cases([s2, s1], [False, True], "case_C1")
        assert mt_class(s2) == "case_C"
