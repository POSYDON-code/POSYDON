"""Unit tests of posydon/popsyn/Moes_distributions.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>",
    "Max Briel <max.briel@gmail.com>",
]

# import the module which will be tested
import posydon.popsyn.Moes_distributions as totest

# aliases
np = totest.np

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import approx, fixture, raises, warns

from posydon.utils.posydonwarning import InappropriateValueWarning


def _integrate_pdf(model, M1, n_lp=801, n_q=901, q_power=0):
    """Integrate q**q_power * model.pdf over log10(P) and q on a grid."""
    lp = np.linspace(model.lp_min, model.lp_max, n_lp)
    q = np.linspace(totest.Q_MIN, 1.0, n_q)
    L, Q = np.meshgrid(lp, q, indexing='ij')
    dens = model.pdf(M1, 10**L, Q) * Q**q_power
    return np.trapz(np.trapz(dens, q, axis=1), lp)


# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['Moe_17_PsandQs', '__authors__',
                    '__builtins__', '__cached__', '__doc__', '__file__',
                    '__loader__', '__name__', '__package__', '__spec__',
                    'np', 'Pwarn',
                    'ALPHA', 'DLOGP', 'Q_MIN', 'Q_SPLIT', 'Q_TWIN', 'E_MIN',
                    'M1_MIN', 'M1_MAX', 'M_STAR_MIN',
                    '_wlog', '_mass_interp', '_powint', '_powinv',
                    'f_twin', 'gamma_largeq', 'gamma_smallq', 'eta', 'e_max',
                    'is_circular', 'f_logP_q03', 'q_components', 'q_factor',
                    'prob_inner_binary', 'f_logP_q01',
                    'inner_binary_density']
        # created when other tests trigger a warning in this module
        totest_elements = set(dir(totest)) - {'__warningregistry__'}
        missing_in_test = set(elements) - totest_elements
        assert len(missing_in_test) == 0, "There are missing objects in "\
                                          +f"{totest.__name__}: "\
                                          +f"{missing_in_test}. Please "\
                                          +"check, whether they have been "\
                                          +"removed on purpose and update "\
                                          +"this unit test."
        new_in_test = totest_elements - set(elements)
        assert len(new_in_test) == 0, "There are new objects in "\
                                      +f"{totest.__name__}: {new_in_test}. "\
                                      +"Please check, whether they have been "\
                                      +"added on purpose and update this "\
                                      +"unit test."


class TestHelpers:

    def test_wlog(self):
        assert totest._wlog(1.2, 1.2, 3.5) == approx(0.0)
        assert totest._wlog(3.5, 1.2, 3.5) == approx(1.0)
        assert totest._wlog(np.sqrt(1.2 * 3.5), 1.2, 3.5) == approx(0.5)
        # clipped outside [lo, hi]
        assert totest._wlog(0.5, 1.2, 3.5) == approx(0.0)
        assert totest._wlog(10.0, 1.2, 3.5) == approx(1.0)

    def test_mass_interp(self):
        args = (0.0, 1.0, 3.0, 1.2, 3.5, 6.0)
        assert totest._mass_interp(1.0, *args) == approx(0.0)
        assert totest._mass_interp(3.5, *args) == approx(1.0)
        assert totest._mass_interp(np.sqrt(3.5 * 6.0), *args) == approx(2.0)
        assert totest._mass_interp(20.0, *args) == approx(3.0)

    def test_powint(self):
        assert totest._powint(0.0, 1.0, 0.0) == approx(1.0)
        assert totest._powint(1.0, 2.0, 1.0) == approx(1.5)
        # g = -1: logarithm
        assert totest._powint(1.0, np.e, -1.0) == approx(1.0)

    def test_powinv(self):
        u = np.linspace(0.0, 1.0, 11)
        for g in [-2.5, -1.0, 0.0, 0.7]:
            x = totest._powinv(u, 0.1, 0.3, g)
            assert x[0] == approx(0.1)
            assert x[-1] == approx(0.3)
            # inverse of the CDF
            cdf = totest._powint(0.1, x, g) / totest._powint(0.1, 0.3, g)
            assert cdf == approx(u)


class TestPaperEquations:

    def test_f_twin(self):
        # Eqns. 5-7
        assert totest.f_twin(1.0, 0.5) == approx(0.3)
        assert totest.f_twin(1.0, 4.5) == approx(0.3 * (1 - 3.5 / 6.0))
        assert totest.f_twin(1.0, 7.5) == approx(0.0)
        assert totest.f_twin(10.0, 1.25) == approx(0.15 * 0.5)
        assert totest.f_twin(10.0, 2.0) == approx(0.0)

    def test_gamma_largeq(self):
        # Eqns. 9-11
        assert totest.gamma_largeq(1.0, 3.0) == approx(-0.5)
        assert totest.gamma_largeq(1.0, 6.0) == approx(-0.8)
        assert totest.gamma_largeq(3.5, 0.5) == approx(-0.5)
        assert totest.gamma_largeq(3.5, 3.0) == approx(-0.9)
        assert totest.gamma_largeq(3.5, 5.5) == approx(-1.6)
        assert totest.gamma_largeq(3.5, 7.0) == approx(-2.0)
        assert totest.gamma_largeq(10.0, 0.5) == approx(-0.5)
        assert totest.gamma_largeq(10.0, 1.5) == approx(-0.95)
        assert totest.gamma_largeq(10.0, 3.4) == approx(-1.82)
        assert totest.gamma_largeq(10.0, 7.0) == approx(-2.0)
        # interpolation in log10(M1): average at the geometric mid point
        mid = np.sqrt(1.2 * 3.5)
        assert totest.gamma_largeq(mid, 3.0) == approx(0.5 * (-0.5 - 0.9))

    def test_gamma_smallq(self):
        # Eqns. 13-15
        assert totest.gamma_smallq(1.0, 4.0) == approx(0.3)
        assert totest.gamma_smallq(3.5, 1.0) == approx(0.2)
        assert totest.gamma_smallq(3.5, 4.0) == approx(-0.25)
        assert totest.gamma_smallq(3.5, 6.5) == approx(-0.9)
        assert totest.gamma_smallq(10.0, 0.5) == approx(0.1)
        assert totest.gamma_smallq(10.0, 2.0) == approx(-0.05)
        assert totest.gamma_smallq(10.0, 4.0) == approx(-0.7)
        assert totest.gamma_smallq(10.0, 6.0) == approx(-1.5)

    def test_eta(self):
        # Eqns. 17-18, held constant below log10(P) = 0.7
        assert totest.eta(1.0, 1.5) == approx(-0.1)
        assert totest.eta(10.0, 1.5) == approx(0.7)
        assert totest.eta(1.0, 0.3) == approx(-2.9)
        assert totest.eta(10.0, 0.3) == approx(-0.1)
        assert totest.eta(np.sqrt(3.0 * 7.0), 1.5) == approx(0.3)

    def test_e_max(self):
        # Eqn. 3
        assert totest.e_max(np.log10(16.0)) == approx(0.75)
        assert totest.e_max(np.log10(2.0)) == approx(0.0)
        assert totest.e_max(0.0) == approx(0.0)

    def test_is_circular(self):
        assert totest.is_circular(np.log10(1.5))
        assert totest.is_circular(np.log10(2.0))
        assert totest.is_circular(np.log10(2.0001))
        assert not totest.is_circular(np.log10(16.0))

    def test_f_logP_q03(self):
        # Eqns. 20-23 at M1 = 1 Msun (log10(M1) = 0)
        assert totest.f_logP_q03(1.0, 0.5) == approx(0.020)
        assert totest.f_logP_q03(1.0, 2.7) == approx(0.039)
        assert totest.f_logP_q03(1.0, 5.5) == approx(0.078)
        assert totest.f_logP_q03(1.0, 6.5) == approx(0.078 * np.exp(-0.3))
        # continuity at the break points
        for lp in [1.0, 2.0, 3.4, 5.5]:
            lo = totest.f_logP_q03(12.0, lp - 1e-9)
            hi = totest.f_logP_q03(12.0, lp + 1e-9)
            assert lo == approx(hi, rel=1e-6)
        # M1 = 10 Msun at short periods
        assert totest.f_logP_q03(10.0, 0.5) == approx(0.13)

    def test_q_components(self):
        for m, lp in [(1.0, 0.5), (3.5, 4.0), (12.0, 1.2), (28.0, 6.0)]:
            w_s, w_l, w_t, gs, gl, c = totest.q_components(m, lp)
            # paper normalisation: integral over 0.3 < q < 1 is 1
            assert w_l + w_t == approx(1.0)
            assert w_t == approx(totest.f_twin(m, lp))
            assert gl == approx(totest.gamma_largeq(m, lp))
            assert gs == approx(totest.gamma_smallq(m, lp))
            # continuity of the density at q = 0.3
            assert c * 0.3**(gl - gs) * 0.3**gs == approx(c * 0.3**gl)
            assert totest.q_factor(m, lp) == approx(1.0 + w_s)

    def test_prob_inner_binary(self):
        assert totest.prob_inner_binary(10.0, 1.0) == approx(1.0)
        assert totest.prob_inner_binary(10.0, 2.5) == approx(0.89)
        assert totest.prob_inner_binary(40.0, 8.0) == approx(0.0)

    def test_frequencies(self):
        m, lp = 7.0, 3.0
        f01 = totest.f_logP_q01(m, lp)
        assert f01 == approx(totest.f_logP_q03(m, lp)
                             * totest.q_factor(m, lp))
        assert totest.inner_binary_density(m, lp) == approx(
            f01 * totest.prob_inner_binary(m, lp))


class TestMoe17PsandQs:

    @fixture(scope="class")
    def model(self):
        return totest.Moe_17_PsandQs(RNG=np.random.default_rng(seed=42))

    @fixture(scope="class")
    def model_short(self):
        """POSYDON default period range, binary fraction capped at 1."""
        with warns(InappropriateValueWarning, match="exceeds 1"):
            return totest.Moe_17_PsandQs(orbital_period_min=0.75,
                                         orbital_period_max=6000.0)

    # __init__

    def test_init_bad_input(self):
        with raises(ValueError, match="orbital_period_min must be positive"):
            totest.Moe_17_PsandQs(orbital_period_min=0.0)
        with raises(ValueError, match="orbital_period_max must be finite"):
            totest.Moe_17_PsandQs(orbital_period_min=10.0,
                                  orbital_period_max=5.0)
        with raises(ValueError, match="M_min must be in"):
            totest.Moe_17_PsandQs(M_min=1.0)

    def test_init(self, model):
        assert model.orbital_period_min == approx(10**0.15)
        assert model.orbital_period_max == approx(1e8)
        assert model.M_min == approx(0.08)
        assert isinstance(model.RNG, np.random.Generator)
        # extra kwargs are ignored
        other = totest.Moe_17_PsandQs(primary_mass_min=7.0)
        assert isinstance(other.RNG, np.random.Generator)

    def test_init_former_crash(self):
        # the grid version crashed for orbital_period_min = 1.0
        with warns(InappropriateValueWarning):
            model = totest.Moe_17_PsandQs(orbital_period_min=1.0)
        assert model.lp_min == approx(0.0)

    def test_repr(self, model):
        assert "1e+08" in repr(model)
        assert "orbital_period_max" in model._repr_html_()

    # binary fraction

    def test_binary_fraction(self, model):
        # Sect. 9.4 and Moe's code: ~41% at 1 Msun, ~96% at 28 Msun
        fb = model.binary_fraction([1.0, 28.0])
        assert fb[0] == approx(0.407, abs=2e-3)
        assert fb[1] == approx(0.955, abs=2e-3)
        # 0.8 and 40 Msun statistics outside the fitted range
        assert model.binary_fraction(100.0) == approx(
            model.binary_fraction(40.0))
        assert model.binary_fraction(0.8) == approx(0.4035, abs=2e-3)
        # low-mass extension: zero at M_min
        assert model.binary_fraction(0.08) == approx(0.0)
        assert 0.0 < model.binary_fraction(0.3) < model.binary_fraction(0.8)

    def test_binary_fraction_capped(self, model_short):
        assert model_short.binary_fraction(40.0) == approx(1.0)
        assert model_short._fbin_raw(40.0) > 1.0

    # PDF

    def test_pdf_normalisation(self, model, model_short):
        for M1 in [0.3, 1.0, 7.0, 40.0]:
            assert _integrate_pdf(model, M1) == approx(
                model.binary_fraction(M1), abs=1e-3)
        # capped: integrates to 1, not to the raw frequency
        assert _integrate_pdf(model_short, 40.0) == approx(1.0, abs=1e-3)

    def test_mean_q(self, model, model_short):
        for m, M1 in [(model, 0.3), (model, 1.0), (model, 7.0),
                      (model, 100.0), (model_short, 40.0)]:
            expected = (_integrate_pdf(m, M1, n_q=3601, q_power=1)
                        / _integrate_pdf(m, M1, n_q=3601))
            assert m.mean_q(M1) == approx(expected, abs=1e-3)
        # vectorised; 0 where there are no binaries
        mq = model.mean_q(np.array([0.05, 1.0, 7.0]))
        assert mq.shape == (3,)
        assert mq[0] == 0.0
        assert mq[1:] == approx([model.mean_q(1.0), model.mean_q(7.0)])
        # the mean of q drawn by the sampler
        M1 = np.full(200000, 10.0)
        M2, P, e = model(M1, RNG=np.random.default_rng(1))
        q = M2[np.isfinite(M2)] / 10.0
        assert q.mean() == approx(model.mean_q(10.0), abs=3e-3)

    def test_logP_pdf_support(self, model):
        assert model.logP_pdf(10.0, 1.0) == approx(0.0)
        assert model.logP_pdf(10.0, 1e9) == approx(0.0)
        assert model.logP_pdf(10.0, 100.0) > 0.0

    def test_q_pdf(self, model):
        q = np.linspace(0.1, 1.0, 9001)
        for M1, P in [(1.0, 5.0), (12.0, 1e3), (0.3, 1e4)]:
            dens = model.q_pdf(M1, P, q)
            assert np.trapz(dens, q) == approx(1.0, abs=1e-3)
        # zero outside the support, including below M_min / M1
        assert model.q_pdf(10.0, 100.0, 0.05) == approx(0.0)
        assert model.q_pdf(0.3, 100.0, 0.2) == approx(0.0)
        assert model.q_pdf(0.3, 100.0, 0.3) > 0.0
        # twin excess starts at q = 0.95 for short periods
        jump = model.q_pdf(1.0, 5.0, 0.95) - model.q_pdf(1.0, 5.0, 0.9499)
        assert jump > 0.0

    def test_e_pdf(self, model):
        e = np.geomspace(totest.E_MIN, 0.999, 100001)
        for M1, P in [(1.0, 5.0), (12.0, 30.0), (28.0, 1e5)]:
            assert np.trapz(model.e_pdf(M1, P, e), e) == approx(1.0,
                                                                abs=1e-3)
        # circular and beyond e_max
        assert model.e_pdf(10.0, 1.5, 0.1) == approx(0.0)
        assert model.e_pdf(10.0, 16.0, 0.8) == approx(0.0)

    def test_pdf_with_e(self, model):
        q = 0.5
        # circular: the e factor is left out
        assert model.pdf(10.0, 1.5, q, e=0.0) == approx(
            model.pdf(10.0, 1.5, q))
        # eccentric: product with e_pdf
        assert model.pdf(10.0, 100.0, q, e=0.3) == approx(
            model.pdf(10.0, 100.0, q) * model.e_pdf(10.0, 100.0, 0.3))

    # sampler

    def test_call_bad_input(self, model):
        with raises(ValueError, match="must all be finite"):
            model(np.array([np.nan]))
        with raises(ValueError, match="must be positive"):
            model(np.array([0.0, 10.0]))

    def test_call_output(self, model):
        M1 = np.full(20000, 10.0)
        M2, P, e = model(M1, RNG=np.random.default_rng(seed=1))
        assert M2.shape == P.shape == e.shape == M1.shape
        single = np.isnan(M2)
        assert np.array_equal(single, np.isnan(P))
        assert np.array_equal(single, np.isnan(e))
        # binary fraction within 4 sigma
        fb = model.binary_fraction(10.0)
        sigma = np.sqrt(fb * (1 - fb) / M1.size)
        assert np.mean(~single) == approx(fb, abs=4 * sigma)
        q, P, e = M2[~single] / 10.0, P[~single], e[~single]
        assert np.all((q >= 0.1) & (q <= 1.0))
        assert np.all((P >= model.orbital_period_min)
                      & (P <= model.orbital_period_max))
        assert np.all(e <= totest.e_max(np.log10(P)))
        assert np.all(e[P <= 2.0] == 0.0)

    def test_call_matches_pdf(self, model):
        """Histogram of log10(P) agrees with logP_pdf."""
        N = 200000
        _, P, _ = model(np.full(N, 7.0), RNG=np.random.default_rng(seed=2))
        edges = np.linspace(model.lp_min, model.lp_max, 21)
        counts, _ = np.histogram(np.log10(P[~np.isnan(P)]), bins=edges)
        expected = []
        for a, b in zip(edges[:-1], edges[1:]):
            lp = np.linspace(a, b, 201)
            expected.append(np.trapz(model.logP_pdf(7.0, 10**lp), lp) * N)
        expected = np.array(expected)
        ok = expected > 50
        chi2 = np.sum((counts[ok] - expected[ok])**2 / expected[ok])
        assert chi2 / ok.sum() < 2.0

    def test_call_seed(self, model):
        M1 = np.full(10, 10.0)
        a = model(M1, RNG=np.random.default_rng(seed=3))
        b = model(M1, RNG=np.random.default_rng(seed=3))
        for x, y in zip(a, b):
            assert np.array_equal(x, y, equal_nan=True)
        # default: the generator's own RNG
        assert len(model(M1)[0]) == 10

    def test_call_all_binaries(self, model):
        M2, P, e = model(np.full(5000, 10.0), all_binaries=True,
                         RNG=np.random.default_rng(seed=4))
        assert not np.any(np.isnan(M2))
        # no pile-up at the maximum period
        assert np.mean(np.isclose(P, model.orbital_period_max)) == 0.0

    def test_call_low_mass(self, model):
        M1 = 0.3
        M2, P, e = model(np.full(5000, M1), all_binaries=True,
                         RNG=np.random.default_rng(seed=5))
        assert np.all(M2 >= model.M_min)
        # primaries just above M_min are (almost) all single
        M2, _, _ = model(np.full(100, 0.0801),
                         RNG=np.random.default_rng(seed=6))
        assert np.all(np.isnan(M2))

    def test_call_below_M_min(self, model):
        # primaries at or below M_min are single stars, not an error
        M1 = np.array([0.05, 0.08, 10.0])
        M2, P, e = model(np.tile(M1, 2000), RNG=np.random.default_rng(seed=8))
        single = np.isnan(M2).reshape(-1, 3)
        assert np.all(single[:, :2])
        assert not np.all(single[:, 2])
        # all_binaries: they stay single, with a warning
        with warns(InappropriateValueWarning, match="2 of 3 primaries"):
            M2, P, e = model(M1, all_binaries=True,
                             RNG=np.random.default_rng(seed=9))
        assert np.all(np.isnan(M2[:2])) and np.isfinite(M2[2])
        assert np.all(np.isnan(P[:2])) and np.all(np.isnan(e[:2]))

    def test_call_no_inner_binary_density(self):
        # 40 Msun: prob_inner_binary = 0 above log10(P) ~ 4.2, so there is no
        # inner binary in 1e5-1e8 d; this used to loop forever
        model = totest.Moe_17_PsandQs(orbital_period_min=1e5,
                                      orbital_period_max=1e8)
        assert model.binary_fraction(40.0) == approx(0.0)
        with warns(InappropriateValueWarning, match="1 of 2 primaries"):
            M2, P, e = model(np.array([40.0, 1.0]), all_binaries=True,
                             RNG=np.random.default_rng(seed=10))
        assert np.isnan(M2[0]) and np.isnan(P[0]) and np.isnan(e[0])
        assert P[1] >= 1e5
        M2, _, _ = model(np.full(100, 40.0),
                         RNG=np.random.default_rng(seed=11))
        assert np.all(np.isnan(M2))

    def test_q_pdf_below_M_min_no_warning(self, model):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            dens = model.q_pdf(np.array([0.05, 0.08]), 10.0, 0.5)
            assert np.all(dens == 0.0)
            assert model.pdf(0.05, 10.0, 0.5, e=0.1) == 0.0

    def test_call_capped(self, model_short):
        """With the binary fraction capped, long periods are still drawn."""
        _, P, _ = model_short(np.full(20000, 40.0),
                              RNG=np.random.default_rng(seed=7))
        assert not np.any(np.isnan(P))
        assert np.log10(P.max()) > model_short.lp_max - 0.1

    def test_rejection_envelope_check(self):
        model = totest.Moe_17_PsandQs()
        model._gmax_table = model._gmax_table * 1e-3
        with raises(RuntimeError, match="Rejection envelope"):
            model(np.full(10, 10.0), all_binaries=True)
