"""Unit tests of posydon/popsyn/Moes_distributions.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

# import the module which will be tested
import posydon.popsyn.Moes_distributions as totest

# aliases
np = totest.np

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import approx, fixture, mark, raises


# define test classes collecting several test functions
class TestElements:
    def test_public_api_contract(self):
        """Assert only required public contract, not full module symbol list."""
        assert hasattr(totest, 'Moe_17_PsandQs')
        assert hasattr(totest, '__authors__')
        assert hasattr(totest, 'np')
        assert hasattr(totest, 'newton_cotes')
        assert hasattr(totest, 'quad')

        model = totest.Moe_17_PsandQs(n_M1=3, n_logP=5, n_q=5, n_e=10)
        assert callable(model)

class TestMoe17PsandQs:

    @fixture
    def small_model(self):
        """Create a Moe_17_PsandQs with small grid for fast testing."""
        return totest.Moe_17_PsandQs(
            n_M1=5, n_logP=10, n_q=10, n_e=20,
            RNG=np.random.default_rng(seed=42))

    # _idl_tabulate

    def test_idl_tabulate_constant(self, small_model):
        """Integral of f=1 from 0 to 1 should be 1."""
        x = np.linspace(0.0, 1.0, 11)
        f = np.ones_like(x)
        result = small_model._idl_tabulate(x, f)
        assert result == approx(1.0, abs=1e-10)

    def test_idl_tabulate_linear(self, small_model):
        """Integral of f=x from 0 to 1 should be 0.5."""
        x = np.linspace(0.0, 1.0, 11)
        f = x.copy()
        result = small_model._idl_tabulate(x, f)
        assert result == approx(0.5, abs=1e-10)

    def test_idl_tabulate_quadratic(self, small_model):
        """Integral of f=x^2 from 0 to 1 should be 1/3."""
        x = np.linspace(0.0, 1.0, 21)
        f = x**2
        result = small_model._idl_tabulate(x, f)
        assert result == approx(1.0 / 3.0, abs=1e-6)

    def test_idl_tabulate_single_point(self, small_model):
        """Single point: integral over zero range should be 0."""
        x = np.array([1.0])
        f = np.array([5.0])
        result = small_model._idl_tabulate(x, f)
        assert result == approx(0.0)

    # __init__

    def test_init_grid_shapes(self, small_model):
        """Verify grid dimensions match requested sizes."""
        assert small_model.numM1 == 5
        assert small_model.numlogP == 10
        assert small_model.numq == 10
        assert small_model.nume == 20
        assert small_model.M1v.shape == (5,)
        assert small_model.logPv.shape == (10,)
        assert small_model.qv.shape == (10,)
        assert small_model.ev.shape == (20,)
        assert small_model.flogP_sq.shape == (10, 5)
        assert small_model.cumqdist.shape == (10, 10, 5)
        assert small_model.cumedist.shape == (20, 10, 5)
        assert small_model.probbin.shape == (10, 5)
        assert small_model.cumPbindist.shape == (10, 5)

    def test_init_mass_range(self, small_model):
        """M1v should span 0.8 to 40 Msun."""
        assert small_model.M1v[0] == approx(0.8, abs=1e-10)
        assert small_model.M1v[-1] == approx(40.0, abs=1e-10)

    def test_init_q_range(self, small_model):
        """qv should span 0.1 to 1.0."""
        assert small_model.qv[0] == approx(0.1)
        assert small_model.qv[-1] == approx(1.0)

    def test_init_cumulative_distributions(self, small_model):
        """Cumulative distributions should end at 1.0."""
        # cumqdist should reach 1.0 at q=1.0 for each (logP, M1)
        for i in range(small_model.numM1):
            for j in range(small_model.numlogP):
                assert small_model.cumqdist[-1, j, i] == approx(1.0, abs=1e-6)
                assert small_model.cumedist[-1, j, i] == approx(1.0, abs=1e-6)

    def test_init_default_params(self):
        """Test with default parameters (expensive — just verify it constructs)."""
        # Use non-default but still small grid to confirm kwarg handling
        model = totest.Moe_17_PsandQs(
            n_M1=3, n_logP=5, n_q=5, n_e=10)
        assert model.numM1 == 3

    # __call__

    def test_call_single_mass(self, small_model):
        """Generate sample for a single primary mass."""
        M2, P, e, Z = small_model(10.0)
        assert len(M2) == 1
        assert len(P) == 1
        assert len(e) == 1
        assert len(Z) == 1
        assert P[0] > 0
        assert Z[0] > 0

    def test_call_array(self, small_model):
        """Generate samples for multiple primary masses."""
        M1 = np.array([5.0, 10.0, 20.0])
        M2, P, e, Z = small_model(M1)
        assert len(M2) == 3
        assert len(P) == 3
        assert len(e) == 3
        assert len(Z) == 3

    def test_seeded_reproducibility(self):
        """Two generators with the same seed must produce identical draws."""
        M1 = np.array([1.0, 5.0, 10.0, 20.0, 35.0])
        model_a = totest.Moe_17_PsandQs(
            n_M1=5, n_logP=10, n_q=10, n_e=20,
            RNG=np.random.default_rng(seed=1234))
        model_b = totest.Moe_17_PsandQs(
            n_M1=5, n_logP=10, n_q=10, n_e=20,
            RNG=np.random.default_rng(seed=1234))

        out_a = model_a(M1, all_binaries=False)
        out_b = model_b(M1, all_binaries=False)

        for arr_a, arr_b in zip(out_a, out_b):
            assert np.array_equal(arr_a, arr_b, equal_nan=True)

    def test_seeded_numeric_snapshot(self):
        """Deterministic numeric snapshot for regression protection."""
        model = totest.Moe_17_PsandQs(
            n_M1=5, n_logP=10, n_q=10, n_e=20,
            RNG=np.random.default_rng(seed=2024))
        M1 = np.array([1.0, 5.0, 10.0, 20.0])

        M2, P, e, Z = model(M1, all_binaries=True)

        expected_M2 = np.array([
            0.2470181820658170,
            0.5979425977891516,
            3.7627175309509338,
            2.0780559303165750,
        ])
        expected_P = np.array([
            1.0000000000000000e+08,
            1.0000000000000000e+08,
            4.3578331121073444e+02,
            3.7529701591440934e+00,
        ])
        expected_e = np.array([
            0.2992502119431766,
            0.2739546694743913,
            0.2807044989336530,
            0.0344916607326516,
        ])
        expected_Z = np.array([
            0.0095611041995118,
            0.0002810281072445,
            0.0033747794810215,
            0.0014212325106382,
        ])

        np.testing.assert_allclose(M2, expected_M2, rtol=0.0, atol=1e-14)
        np.testing.assert_allclose(P, expected_P, rtol=0.0, atol=1e-10)
        np.testing.assert_allclose(e, expected_e, rtol=0.0, atol=1e-14)
        np.testing.assert_allclose(Z, expected_Z, rtol=0.0, atol=1e-14)

    def test_call_all_binaries_true(self):
        """With all_binaries=True, no single stars should be produced."""
        model = totest.Moe_17_PsandQs(
            n_M1=5, n_logP=10, n_q=10, n_e=20,
            RNG=np.random.default_rng(seed=42))
        M1 = np.array([10.0] * 20)
        M2, P, e, Z = model(M1, all_binaries=True)
        # all_binaries=True means mybinfrac=1.0, so no NaN values
        assert not np.any(np.isnan(M2))
        assert not np.any(np.isnan(P))
        assert not np.any(np.isnan(e))
        assert np.all(np.isfinite(Z))

        q = M2 / M1
        assert np.all(q >= 0.1)
        assert np.all(q <= 1.0)
        assert np.all(P > 0.0)
        assert np.all(e >= 0.0)
        assert np.all(e < 1.0)

    def test_call_all_binaries_false(self):
        """With all_binaries=False, some single stars may be produced."""
        model = totest.Moe_17_PsandQs(
            n_M1=5, n_logP=10, n_q=10, n_e=20,
            RNG=np.random.default_rng(seed=0))
        M1 = np.array([1.0] * 50)
        M2, P, e, Z = model(M1, all_binaries=False)
        # With 50 draws at M1=1.0, expect some single stars (NaN)
        # Z is always set, never NaN
        assert not np.any(np.isnan(Z))
        assert len(M2) == 50

    def test_call_single_star_branch_consistency(self):
        """When binary fraction is forced to zero, all companions are NaN."""
        model = totest.Moe_17_PsandQs(
            n_M1=5, n_logP=10, n_q=10, n_e=20,
            RNG=np.random.default_rng(seed=42))

        M1 = np.array([0.08] * 10)
        M2, P, e, Z = model(M1, M_min=0.08, all_binaries=False)

        assert np.all(np.isnan(M2))
        assert np.all(np.isnan(P))
        assert np.all(np.isnan(e))
        assert np.all(np.isfinite(Z))

    def test_call_high_mass(self, small_model):
        """M1 > 40 Msun should adopt binary statistics of M1 = 40 Msun."""
        M2, P, e, Z = small_model(80.0)
        assert len(M2) == 1
        assert P[0] > 0

    def test_call_low_mass(self):
        """M1 < 0.8 Msun should rescale binary fraction."""
        model = totest.Moe_17_PsandQs(
            n_M1=5, n_logP=10, n_q=10, n_e=20,
            RNG=np.random.default_rng(seed=42))
        M2, P, e, Z = model(0.5, M_min=0.08, all_binaries=False)
        assert len(M2) == 1
        assert Z[0] > 0

    def test_call_low_mass_q_truncation(self):
        """M1 < 0.8 inside the binary path should truncate q distribution."""
        model = totest.Moe_17_PsandQs(
            n_M1=5, n_logP=10, n_q=10, n_e=20,
            RNG=np.random.default_rng(seed=42))
        M2, P, e, Z = model(0.5, M_min=0.08, all_binaries=True)
        assert len(M2) == 1
        # M2 = M1 * q, and q >= q_min = M_min/M1 = 0.08/0.5 = 0.16
        assert M2[0] >= 0.08 * 0.99  # M2 >= M_min (small tolerance)

    def test_call_metallicity_range(self, small_model):
        """Metallicities should be within the expected range."""
        M1 = np.array([10.0] * 100)
        _, _, _, Z = small_model(M1)
        Zsun = 0.02
        Z_min = Zsun * 10**(-2.3)
        Z_max = Zsun * 10**(0.176)
        assert all(Z >= Z_min * 0.99)  # small tolerance
        assert all(Z <= Z_max * 1.01)

    def test_cumulative_monotonicity(self, small_model):
        """Cumulative distributions should be non-decreasing."""
        tol = 1e-10
        assert np.all(np.diff(small_model.cumqdist, axis=0) >= -tol)
        assert np.all(np.diff(small_model.cumedist, axis=0) >= -tol)
        assert np.all(np.diff(small_model.cumPbindist, axis=0) >= -tol)

    @mark.parametrize(
        "M1, M_min, M_max",
        [
            (np.array([-1.0, 0.0]), 0.08, 150.0),
            (np.array([np.nan, 1.0]), 0.08, 150.0),
            (np.array([0.05, 1.0]), 0.08, 150.0),
            (np.array([1.0, 160.0]), 0.08, 150.0),
            (np.array([1.0]), 0.0, 150.0),
            (np.array([1.0]), 0.5, 0.5),
        ],
    )
    def test_invalid_input_raises(self, small_model, M1, M_min, M_max):
        """Invalid mass inputs and bounds should raise ValueError."""
        with raises(ValueError):
            small_model(M1, M_min=M_min, M_max=M_max)
