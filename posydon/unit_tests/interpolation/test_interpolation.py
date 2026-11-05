"""Unit tests of GRIDInterpolator.get in posydon/interpolation/interpolation.py
(per-mass caching must not change any returned value).
"""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>"
]

import numpy as np
from pytest import fixture, raises

import posydon.interpolation.interpolation as totest


def _reference_get(self, key, M_new):
    """GRIDInterpolator.get as it was before caching was added."""
    if M_new in self.grid_mass:
        try:
            kvalue = self.grid_data[M_new][key]
        except KeyError:
            self.load_grid(M_new)
            kvalue = self.grid_data[M_new][key]
        log_L = self.grid_data[M_new]['log_L']
        log_LH = self.grid_data[M_new]['log_LH']
        i_zams = np.argmax(np.log10(0.99) + log_L <= log_LH)
    else:
        mass_low = np.max(self.grid_mass[M_new > self.grid_mass])
        mass_high = np.min(self.grid_mass[M_new < self.grid_mass])
        weight = (M_new - mass_low) / (mass_high - mass_low)
        try:
            kvalue_low = self.grid_data[mass_low][key]
        except KeyError:
            self.load_grid(mass_low)
            kvalue_low = self.grid_data[mass_low][key]
        try:
            kvalue_high = self.grid_data[mass_high][key]
        except KeyError:
            self.load_grid(mass_high)
            kvalue_high = self.grid_data[mass_high][key]
        i_cut = min(len(kvalue_low), len(kvalue_high))
        kvalue = kvalue_low[:i_cut] + weight * (kvalue_high[:i_cut]
                                                - kvalue_low[:i_cut])
        log_L_low = self.grid_data[mass_low]['log_L']
        log_L_high = self.grid_data[mass_high]['log_L']
        log_LH_low = self.grid_data[mass_low]['log_LH']
        log_LH_high = self.grid_data[mass_high]['log_LH']
        log_L = log_L_low[:i_cut] + weight * (log_L_high[:i_cut]
                                              - log_L_low[:i_cut])
        log_LH = log_LH_low[:i_cut] + weight * (log_LH_high[:i_cut]
                                                - log_LH_low[:i_cut])
        i_zams = np.argmax(np.log10(0.99) + log_L <= log_LH)
    if key == 'age':
        kvalue = kvalue - kvalue[i_zams]
    return kvalue[i_zams:]


def _make_grid(seed=1, n_masses=25, lazy=False):
    rng = np.random.default_rng(seed)
    g = object.__new__(totest.GRIDInterpolator)
    g.grid_mass = np.sort(rng.uniform(0.5, 120.0, n_masses))
    full = {}
    for m in g.grid_mass:
        n = int(rng.integers(150, 900))
        log_L = np.cumsum(rng.normal(0.01, 0.02, n))
        d = {'log_L': log_L,
             'log_LH': log_L - np.abs(rng.normal(0.01, 0.05, n))
                       + (np.arange(n) > rng.integers(0, 25)) * 0.02,
             'age': np.cumsum(rng.uniform(1e3, 1e5, n)),
             'mass': m - np.cumsum(rng.uniform(0, 1e-3, n)),
             'center_h1': rng.uniform(size=n)}
        full[m] = d
    g.grid_data = {} if lazy else dict(full)

    def load_grid(*masses):
        for m in masses:
            g.grid_data[m] = full[m]
    g.load_grid = load_grid
    return g


class TestGRIDInterpolatorGet:

    @fixture
    def grid(self):
        return _make_grid()

    def _check(self, g, ref, masses, keys):
        for M in masses:
            for key in keys:
                try:
                    expected = _reference_get(ref, key, M)
                except ValueError:
                    with raises(ValueError):
                        g.get(key, M)
                    continue
                got = g.get(key, M)
                assert got.dtype == expected.dtype
                assert np.array_equal(got, expected)

    def test_identical_to_uncached(self, grid):
        ref = _make_grid()
        rng = np.random.default_rng(2)
        masses = (list(rng.uniform(0.2, 130.0, 150))     # includes off-grid
                  + list(grid.grid_mass[:8])               # exact grid masses
                  + [np.float64(grid.grid_mass[3])])
        keys = ['age', 'mass', 'center_h1', 'log_L', 'age']
        self._check(grid, ref, masses, keys)
        self._check(grid, ref, masses, keys)   # second pass: cache hits
        assert len(grid._get_cache) > 0

    def test_lazy_loading(self):
        g, ref = _make_grid(lazy=True), _make_grid(lazy=True)
        self._check(g, ref, [5.0, 37.3, 80.1, g.grid_mass[4]], ['mass', 'age'])

    def test_cache_is_bounded(self, grid):
        grid._GET_CACHE_SIZE = 4
        for M in np.linspace(10, 100, 20):
            grid.get('mass', M)
        assert len(grid._get_cache) <= 4

    def test_unhashable_mass(self, grid):
        ref = _make_grid()
        M = np.array(42.0)          # 0-d array is not hashable
        assert np.array_equal(grid.get('age', M), _reference_get(ref, 'age', M))
