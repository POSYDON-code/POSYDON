"""Initial binary parameters following Moe & Di Stefano (2017).

https://ui.adsabs.harvard.edu/abs/2017ApJS..230...15M/abstract

The fits of Sect. 9 of the paper are evaluated directly (Eqns. 3, 5-7, 9-11,
13-15, 17-18 and 20-23), so the sampler and its probability density function
(PDF) are exact and continuous in the primary mass M1 and the orbital period P:

* the mass ratio q and the eccentricity e are drawn from closed-form inverse
  cumulative distributions at the drawn (M1, P) of each system;
* log10(P) is drawn by rejection sampling from the companion frequency;
* ``Moe_17_PsandQs.pdf`` is the density the sampler draws from.

The earlier implementation was an adoption of the translation of Moe's IDL code
into python by Mads Sørensen (V. 0.1; 2017/02/03). Two ingredients of that code
are not in the paper and are kept here:

* ``prob_inner_binary``: the probability that a companion at period P is the
  inner binary of the system. Following Maxwell Moe, only single stars and
  (inner) binaries are generated; outer tertiaries in hierarchical triples are
  not. This approximation reproduces the overall multiplicity statistics.
* the eccentricity distribution turns over linearly between 0.8 e_max and
  e_max, and has a lower cut-off ``E_MIN`` (needed for eta < -1).

Notation: ``lp`` = log10(P / day), ``m`` = M1 / Msun.
"""


__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
    "Max Briel <max.briel@gmail.com>",
]


import numpy as np

from posydon.utils.posydonwarning import Pwarn

ALPHA = 0.018          # Eqn. 23: slope across the intermediate-period bump
DLOGP = 0.7            # Eqn. 23: half width of the bump
Q_MIN = 0.1            # smallest mass ratio of the M+D17 statistics
Q_SPLIT = 0.3          # break between the small-q and large-q power laws
Q_TWIN = 0.95          # lower edge of the twin excess
E_MIN = 1e-4           # lower cut-off of the eccentricity power law
M1_MIN = 0.8           # primary-mass range with M+D17 statistics
M1_MAX = 40.0
M_STAR_MIN = 0.08      # minimum stellar mass (low-mass extension)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _wlog(m, lo, hi):
    """Interpolation weight in log10(M1) between lo and hi, clipped to [0, 1]."""
    return np.clip(np.log10(m / lo) / np.log10(hi / lo), 0.0, 1.0)


def _mass_interp(m, v_lo, v_mid, v_hi, m_lo, m_mid, m_hi):
    """Interpolate in log10(M1) between fits at m_lo, m_mid and m_hi."""
    return np.where(m <= m_mid,
                    v_lo + (v_mid - v_lo) * _wlog(m, m_lo, m_mid),
                    v_mid + (v_hi - v_mid) * _wlog(m, m_mid, m_hi))


def _powint(a, b, g):
    """Integral of x**g from a to b (vectorised, including g = -1)."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    k = g + 1.0
    small = np.abs(k) < 1e-10
    ks = np.where(small, 1.0, k)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(small, np.log(b / a), (b**ks - a**ks) / ks)


def _powinv(u, a, b, g):
    """Inverse cumulative distribution of p(x) ~ x**g on [a, b]."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    k = g + 1.0
    small = np.abs(k) < 1e-10
    ks = np.where(small, 1.0, k)
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        x = (a**ks + u * (b**ks - a**ks))**(1.0 / ks)
        x_log = a * (b / a)**u
    return np.where(small, x_log, x)


# ---------------------------------------------------------------------------
# Paper equations (Sect. 9)
# ---------------------------------------------------------------------------
def f_twin(m, lp):
    """Excess twin fraction, Eqns. 5-7."""
    ft1 = 0.30 - 0.15 * np.log10(m)
    lptwin = np.where(m <= 6.5, 8.0 - m, 1.5)
    return np.where(lp < 1.0, ft1,
                    np.where(lp < lptwin,
                             ft1 * (1.0 - (lp - 1.0) / (lptwin - 1.0)), 0.0))


def gamma_largeq(m, lp):
    """Power-law slope across 0.3 < q < 1, Eqns. 9-11."""
    g12 = np.where(lp < 5.0, -0.5, -0.5 - 0.3 * (lp - 5.0))
    g35 = np.select([lp < 1.0, lp < 4.5, lp < 6.5],
                    [-0.5, -0.5 - 0.2 * (lp - 1.0), -1.2 - 0.4 * (lp - 4.5)],
                    -2.0)
    g6 = np.select([lp < 1.0, lp < 2.0, lp < 4.0],
                   [-0.5, -0.5 - 0.9 * (lp - 1.0), -1.4 - 0.3 * (lp - 2.0)],
                   -2.0)
    return _mass_interp(m, g12, g35, g6, 1.2, 3.5, 6.0)


def gamma_smallq(m, lp):
    """Power-law slope across 0.1 < q < 0.3, Eqns. 13-15."""
    g12 = np.full(np.broadcast(m, lp).shape, 0.3)
    g35 = np.select([lp < 2.5, lp < 5.5],
                    [0.2, 0.2 - 0.3 * (lp - 2.5)], -0.7 - 0.2 * (lp - 5.5))
    g6 = np.select([lp < 1.0, lp < 3.0, lp < 5.6],
                   [0.1, 0.1 - 0.15 * (lp - 1.0), -0.2 - 0.5 * (lp - 3.0)],
                   -1.5)
    return _mass_interp(m, g12, g35, g6, 1.2, 3.5, 6.0)


def eta(m, lp):
    """Power-law slope of the eccentricity distribution, Eqns. 17-18.

    The fits diverge at log10(P) = 0.5; as in Moe's code, eta is held at its
    log10(P) = 0.7 value below that.
    """
    x = np.maximum(lp, 0.7) - 0.5
    e3 = 0.6 - 0.7 / x
    e7 = 0.9 - 0.2 / x
    return e3 + (e7 - e3) * _wlog(m, 3.0, 7.0)


def e_max(lp):
    """Maximum eccentricity, Eqn. 3 (0 for P <= 2 d, which are circular)."""
    return np.where(lp > np.log10(2.0),
                    1.0 - (10**np.maximum(lp, np.log10(2.0)) / 2.0)
                    **(-2.0 / 3.0), 0.0)


def is_circular(lp):
    """True for P <= 2 d (paper) or when e_max is too small to exceed E_MIN."""
    return e_max(lp) <= E_MIN / 0.8


def f_logP_q03(m, lp):
    """Companion frequency per decade of P with q > 0.3, Eqns. 20-23."""
    lm = np.log10(m)
    f1 = 0.020 + 0.04 * lm + 0.07 * lm**2
    f27 = 0.039 + 0.07 * lm + 0.01 * lm**2
    f55 = 0.078 - 0.05 * lm + 0.04 * lm**2
    a, b = 2.7 - DLOGP, 2.7 + DLOGP
    return np.select(
        [lp < 1.0, lp < a, lp < b, lp < 5.5],
        [f1,
         f1 + (lp - 1.0) / (1.7 - DLOGP) * (f27 - f1 - ALPHA * DLOGP),
         f27 + ALPHA * (lp - 2.7),
         f27 + ALPHA * DLOGP + (lp - b) / (2.8 - DLOGP)
         * (f55 - f27 - ALPHA * DLOGP)],
        f55 * np.exp(-0.3 * (lp - 5.5)))


def q_components(m, lp):
    """Masses and shape parameters of the three mass-ratio components.

    The mass-ratio distribution p_q is normalised so that its integral over
    0.3 < q < 1 is 1, as in the paper. It is the sum of

    * small q: c * 0.3**(gl - gs) * q**gs on [0.1, 0.3) with mass w_small,
    * large q: c * q**gl on [0.3, 1] with mass 1 - F_twin,
    * twins: F_twin / 0.05 on [0.95, 1] with mass F_twin.

    Returns
    -------
    w_small, w_large, w_twin, gs, gl, c : ndarrays
    """
    gl = gamma_largeq(m, lp)
    gs = gamma_smallq(m, lp)
    ft = f_twin(m, lp)
    c = (1.0 - ft) / _powint(Q_SPLIT, 1.0, gl)
    w_small = c * Q_SPLIT**(gl - gs) * _powint(Q_MIN, Q_SPLIT, gs)
    return w_small, 1.0 - ft, ft, gs, gl, c


def q_factor(m, lp):
    """Ratio f_logP;q>0.1 / f_logP;q>0.3, i.e. the integral of p_q over q>0.1."""
    w_s, w_l, w_t, *_ = q_components(m, lp)
    return w_s + w_l + w_t


def prob_inner_binary(m, lp):
    """Probability that a companion is the inner binary (Moe's code).

    1 for log10(P) <= 1.5; tuned to reproduce the multiplicity statistics of
    Sect. 9.4 (41% binaries at 1 Msun, 96% at 28 Msun).
    """
    p = 1.0 - 0.11 * np.clip(lp - 1.5, 0.0, None)**1.43 * (m / 10.0)**0.56
    return np.clip(p, 0.0, 1.0)


def f_logP_q01(m, lp):
    """Companion frequency per decade of P with q > 0.1."""
    return f_logP_q03(m, lp) * q_factor(m, lp)


def inner_binary_density(m, lp):
    """Inner-binary density per decade of P (integrates to f_bin)."""
    return f_logP_q01(m, lp) * prob_inner_binary(m, lp)


# ---------------------------------------------------------------------------
# Sampler and PDF
# ---------------------------------------------------------------------------
class Moe_17_PsandQs():
    """Sampler and PDF of the Moe & Di Stefano (2017) [1]_ binary model.

    For a primary of mass M1 the model gives one inner binary companion with
    probability ``binary_fraction(M1)`` and a single star otherwise. The binary
    has a log10(P) distribution, and q and e distributions conditional on
    (M1, P). Primaries outside 0.8-40 Msun use the statistics of 0.8 or
    40 Msun; below 0.8 Msun the binary fraction is scaled down to zero at
    ``M_min`` and q is truncated at ``M_min / M1``. Primaries with
    M1 <= ``M_min`` cannot have a companion of at least ``M_min`` with q <= 1
    and are always single stars.

    References
    ----------
    .. [1] Moe, M. and Di Stefano, R., “Mind Your Ps and Qs: The Interrelation
        between Period (P) and Mass-ratio (Q) Distributions of Binary Stars”,
        <i>The Astrophysical Journal Supplement Series</i>, vol. 230, no. 2,
        Art. no. 15, IOP, 2017. doi:10.3847/1538-4365/aa6fb6.
    """

    def __init__(self, orbital_period_min=10**0.15, orbital_period_max=1e8,
                 M_min=M_STAR_MIN, n_M1_table=600, n_logP_table=4001,
                 RNG=None, **kwargs):
        """Initialise the model for a period range.

        Parameters
        ----------
        orbital_period_min, orbital_period_max : float (default: 10**0.15, 1e8)
            Period range in days. A primary whose inner companion would lie
            outside this range is a single star.
        M_min : float (default: 0.08)
            Minimum stellar mass in Msun, used only for primaries below
            0.8 Msun.
        n_M1_table : int (default: 600)
            Number of log10(M1) nodes of the binary-fraction table.
        n_logP_table : int (default: 4001)
            Number of log10(P) nodes used to integrate the binary fraction.
        RNG : numpy.random.Generator or None (default: None)
            Default random-number generator of ``__call__``.
        **kwargs : dict
            Ignored; allows passing the full population settings.
        """
        if not (np.isfinite(orbital_period_min) and orbital_period_min > 0.0):
            raise ValueError("orbital_period_min must be positive and finite.")
        if not (np.isfinite(orbital_period_max)
                and orbital_period_max > orbital_period_min):
            raise ValueError("orbital_period_max must be finite and larger "
                             "than orbital_period_min.")
        if not (0.0 < M_min < M1_MIN):
            raise ValueError(f"M_min must be in (0, {M1_MIN}) Msun.")
        self.orbital_period_min = orbital_period_min
        self.orbital_period_max = orbital_period_max
        self.lp_min = np.log10(orbital_period_min)
        self.lp_max = np.log10(orbital_period_max)
        self.M_min = M_min
        self.RNG = np.random.default_rng() if RNG is None else RNG

        # tables of the binary fraction and of the maximum period density on
        # a log10(M1) grid; the kinks of the mass interpolation are included
        kinks = [1.2, 3.0, 3.5, 6.0, 6.5, 7.0]
        m = np.unique(np.concatenate(
            [np.geomspace(M1_MIN, M1_MAX, n_M1_table), kinks]))
        lp = np.linspace(self.lp_min, self.lp_max, n_logP_table)
        g = inner_binary_density(m[:, None], lp[None, :])
        self._n_logP_table = n_logP_table
        self._lm_table = np.log10(m)
        self._fbin_raw_table = np.trapz(g, lp, axis=1)
        self._gmax_table = g.max(axis=1)
        if self._fbin_raw_table.max() > 1.0:
            m_exceed = m[np.argmax(self._fbin_raw_table > 1.0)]
            Pwarn("The integrated inner-binary frequency of Moe & Di Stefano "
                  f"(2017) exceeds 1 for M1 >= {m_exceed:.1f} Msun (maximum "
                  f"{self._fbin_raw_table.max():.3f}), because "
                  f"orbital_period_min = {orbital_period_min} d is below the "
                  "range of the fits (log10(P) = 0.15-0.2). The binary "
                  "fraction is capped at 1 and the period distribution is "
                  "rescaled there.", "InappropriateValueWarning")

    def __repr__(self):
        return ("Moe and Di Stefano 2017 distributions for periods "
                f"{self.orbital_period_min:g} - {self.orbital_period_max:g} d")

    def _repr_html_(self):
        return ("<h3>Moe and Di Stefano 2017 distributions</h3>"
                f"<p>orbital_period_min={self.orbital_period_min:g} d</p>"
                f"<p>orbital_period_max={self.orbital_period_max:g} d</p>")

    # ---- binary fraction --------------------------------------------------
    @staticmethod
    def _m_stats(M1):
        """Mass whose statistics are used (clipped to 0.8-40 Msun)."""
        return np.clip(M1, M1_MIN, M1_MAX)

    def _fbin_raw(self, M1):
        """Integrated inner-binary frequency (may exceed 1)."""
        return np.interp(np.log10(self._m_stats(M1)), self._lm_table,
                         self._fbin_raw_table)

    def _lowmass_scale(self, M1):
        """Binary-fraction scaling: 0 at M_min, 1 from 0.8 Msun upwards."""
        return np.interp(np.log10(M1), np.log10([self.M_min, M1_MIN]),
                         [0.0, 1.0])

    def binary_fraction(self, M1):
        """Fraction of primaries with an inner companion in the period range.

        Parameters
        ----------
        M1 : float or ndarray
            Primary mass(es) in Msun.

        Returns
        -------
        ndarray
            Binary fraction, between 0 and 1.
        """
        M1 = np.asarray(M1, dtype=float)
        return np.minimum(self._fbin_raw(M1), 1.0) * self._lowmass_scale(M1)

    # ---- PDF --------------------------------------------------------------
    def _q_lower(self, M1):
        return np.maximum(Q_MIN, self.M_min / M1)

    def _q_masses(self, m, lp, q_lo):
        """Masses of the three q components above q_lo (q_lo >= 0.1).

        Returns
        -------
        (small, large, twin), (w_t, gs, gl, c)
            Component masses above q_lo, and the shape parameters of
            ``q_components`` so that callers do not evaluate it twice.
        """
        _, _, w_t, gs, gl, c = q_components(m, lp)
        small = (c * Q_SPLIT**(gl - gs)
                 * _powint(np.minimum(q_lo, Q_SPLIT), Q_SPLIT, gs))
        large = c * _powint(np.clip(q_lo, Q_SPLIT, 1.0), 1.0, gl)
        twin = w_t * (1.0 - np.clip(q_lo, Q_TWIN, 1.0)) / (1.0 - Q_TWIN)
        return (small, large, twin), (w_t, gs, gl, c)

    def _q_first_moments(self, m, lp, q_lo, shape):
        """First moments in q of the three q components above q_lo.

        Counterpart of ``_q_masses``; ``shape`` is its second return value.
        """
        w_t, gs, gl, c = shape
        small = (c * Q_SPLIT**(gl - gs)
                 * _powint(np.minimum(q_lo, Q_SPLIT), Q_SPLIT, gs + 1.0))
        large = c * _powint(np.clip(q_lo, Q_SPLIT, 1.0), 1.0, gl + 1.0)
        q_t = np.clip(q_lo, Q_TWIN, 1.0)
        twin = w_t * 0.5 * (1.0 - q_t**2) / (1.0 - Q_TWIN)
        return small, large, twin

    def mean_q(self, M1):
        """Mean mass ratio of the binaries, E[q | M1, binary].

        The mean of q over ``pdf`` at fixed M1, i.e. averaged over the period
        distribution in the period range. It is 0 where the binary fraction
        is 0.

        Parameters
        ----------
        M1 : float or ndarray
            Primary mass(es) in Msun.

        Returns
        -------
        ndarray
            Mean mass ratio, of the shape of M1.
        """
        M1 = np.asarray(M1, dtype=float)
        lp = np.linspace(self.lp_min, self.lp_max, self._n_logP_table)
        m1 = M1[..., None]
        q_lo = self._q_lower(m1)
        masses, shape = self._q_masses(self._m_stats(m1), lp, q_lo)
        moments = self._q_first_moments(self._m_stats(m1), lp, q_lo, shape)
        norm = sum(masses)
        with np.errstate(invalid='ignore', divide='ignore'):
            mean_q_P = np.where(norm > 0.0, sum(moments) / norm, 0.0)
        g = self.logP_pdf(m1, 10**lp)
        weight = np.trapz(g, lp, axis=-1)
        with np.errstate(invalid='ignore', divide='ignore'):
            out = np.trapz(g * mean_q_P, lp, axis=-1) / weight
        return np.where(weight > 0.0, out, 0.0)

    def logP_pdf(self, M1, P):
        """Binary period density per decade, p(log10 P, binary | M1).

        Integrates over log10(P) to ``binary_fraction(M1)``.
        """
        M1, lp = np.broadcast_arrays(np.asarray(M1, dtype=float),
                                     np.log10(P))
        g = inner_binary_density(self._m_stats(M1), lp)
        g = g / np.maximum(self._fbin_raw(M1), 1.0) * self._lowmass_scale(M1)
        inside = (lp >= self.lp_min) & (lp <= self.lp_max)
        return np.where(inside, g, 0.0)

    def q_pdf(self, M1, P, q):
        """Mass-ratio density p(q | M1, P, binary), normalised to 1."""
        M1, lp, q = np.broadcast_arrays(np.asarray(M1, dtype=float),
                                        np.log10(P), np.asarray(q, dtype=float))
        m = self._m_stats(M1)
        q_lo = self._q_lower(M1)
        masses, (w_t, gs, gl, c) = self._q_masses(m, lp, q_lo)
        norm = sum(masses)
        # norm is 0 for M1 <= M_min (q_lo >= 1); those are masked to 0 below
        with np.errstate(invalid='ignore', divide='ignore', over='ignore'):
            dens = np.where(q < Q_SPLIT, c * Q_SPLIT**(gl - gs) * q**gs,
                            c * q**gl + np.where(q >= Q_TWIN,
                                                 w_t / (1.0 - Q_TWIN), 0.0))
            out = dens / norm
        ok = (q >= q_lo) & (q <= 1.0) & (norm > 0.0)
        return np.where(ok, out, 0.0)

    def e_pdf(self, M1, P, e):
        """Eccentricity density p(e | M1, P, binary), normalised to 1.

        Zero for circular orbits (``is_circular``), where e = 0 exactly.
        """
        M1, lp, e = np.broadcast_arrays(np.asarray(M1, dtype=float),
                                        np.log10(P), np.asarray(e, dtype=float))
        et = eta(self._m_stats(M1), lp)
        emx = e_max(lp)
        e8 = 0.8 * emx
        with np.errstate(invalid='ignore', divide='ignore', over='ignore'):
            h = e8**et
            norm = _powint(E_MIN, e8, et) + 0.5 * (emx - e8) * h
            dens = np.where(e < e8, e**et, h * (emx - e) / (emx - e8))
            out = dens / norm
        ok = (e >= E_MIN) & (e <= emx) & ~is_circular(lp)
        return np.where(ok, out, 0.0)

    def pdf(self, M1, P, q, e=None):
        """Joint binary density p(log10 P, q[, e], binary | M1).

        Integrates to ``binary_fraction(M1)``; single stars have probability
        ``1 - binary_fraction(M1)``. Circular binaries have e = 0 exactly, so
        for them the eccentricity factor is left out.

        Parameters
        ----------
        M1 : float or ndarray
            Primary mass(es) in Msun.
        P : float or ndarray
            Orbital period(s) in days.
        q : float or ndarray
            Mass ratio(s) M2 / M1.
        e : float or ndarray or None (default: None)
            Eccentricity(ies). If None, the density in (log10 P, q) is given.

        Returns
        -------
        ndarray
            Density per decade of P and per unit q (and e).
        """
        dens = self.logP_pdf(M1, P) * self.q_pdf(M1, P, q)
        if e is not None:
            circular = is_circular(np.log10(P))
            dens = dens * np.where(circular, 1.0, self.e_pdf(M1, P, e))
        return dens

    # ---- sampler ----------------------------------------------------------
    def _sample_logP(self, M1, rng):
        """Rejection sampling of log10(P) from inner_binary_density.

        Returns NaN for primaries whose inner-binary density is zero over the
        whole period range (no inner binary possible), instead of looping
        forever.
        """
        m = self._m_stats(M1)
        envelope = 1.02 * np.interp(np.log10(m), self._lm_table,
                                    self._gmax_table)
        lp_out = np.full(M1.shape, np.nan)
        todo = np.flatnonzero(envelope > 0.0)
        while todo.size > 0:
            lp = rng.uniform(self.lp_min, self.lp_max, todo.size)
            g = inner_binary_density(m[todo], lp)
            if np.any(g > envelope[todo]):
                raise RuntimeError("Rejection envelope of the Moe & Di "
                                   "Stefano period distribution too small.")
            accept = rng.uniform(0.0, envelope[todo]) < g
            lp_out[todo[accept]] = lp[accept]
            todo = todo[~accept]
        return lp_out

    def _sample_q(self, M1, lp, rng):
        """Mixture sampling of q: choose a component, invert its CDF."""
        m = self._m_stats(M1)
        q_lo = self._q_lower(M1)
        (w_s, w_l, w_t), (_, gs, gl, _) = self._q_masses(m, lp, q_lo)
        u = rng.random(M1.size) * (w_s + w_l + w_t)
        v = rng.random(M1.size)
        lo_twin = np.clip(q_lo, Q_TWIN, 1.0)
        return np.where(
            u < w_s, _powinv(v, np.minimum(q_lo, Q_SPLIT), Q_SPLIT, gs),
            np.where(u < w_s + w_l,
                     _powinv(v, np.clip(q_lo, Q_SPLIT, 1.0), 1.0, gl),
                     lo_twin + v * (1.0 - lo_twin)))

    def _sample_e(self, M1, lp, rng):
        """Mixture sampling of e: power law up to 0.8 e_max, then linear."""
        e = np.zeros(M1.size)
        ecc = ~is_circular(lp)
        et = eta(self._m_stats(M1[ecc]), lp[ecc])
        emx = e_max(lp[ecc])
        e8 = 0.8 * emx
        a_pow = _powint(E_MIN, e8, et)
        a_tri = 0.5 * (emx - e8) * e8**et
        u = rng.random(ecc.sum()) * (a_pow + a_tri)
        v = rng.random(ecc.sum())
        e[ecc] = np.where(u < a_pow, _powinv(v, E_MIN, e8, et),
                          emx - (emx - e8) * np.sqrt(v))
        return e

    def __call__(self, M1, all_binaries=False, RNG=None):
        """Draw companions for the primaries M1.

        Parameters
        ----------
        M1 : float or ndarray of float
            Mass(es) of the primary in Msun; must be positive. Primaries with
            M1 <= ``M_min`` are single stars.
        all_binaries : bool (default: False)
            If True, every primary that can have an inner companion gets one;
            the period distribution is then normalised to 1. Primaries with
            M1 <= ``M_min``, or without any inner-binary density in the
            period range, stay single (with a warning).
        RNG : numpy.random.Generator or None (default: None)
            Random-number generator; ``self.RNG`` if None.

        Returns
        -------
        M2, P, e : ndarrays of float
            Secondary mass in Msun, period in days and eccentricity; NaN for
            single stars.
        """
        rng = self.RNG if RNG is None else RNG
        M1 = np.atleast_1d(np.asarray(M1, dtype=float))
        if not np.all(np.isfinite(M1)):
            raise ValueError("Primary masses M1 must all be finite.")
        if np.any(M1 <= 0.0):
            raise ValueError("Primary masses M1 must be positive.")
        # no companion of at least M_min with q <= 1 is possible
        can_pair = M1 > self.M_min
        if all_binaries:
            is_binary = can_pair.copy()
        else:
            is_binary = rng.random(M1.size) < self.binary_fraction(M1)
        M2 = np.full(M1.size, np.nan)
        P = np.full(M1.size, np.nan)
        e = np.full(M1.size, np.nan)
        idx = np.flatnonzero(is_binary)
        lp = self._sample_logP(M1[idx], rng)
        # no inner-binary density anywhere in the period range: single star
        has_P = np.isfinite(lp)
        idx, lp = idx[has_P], lp[has_P]
        if all_binaries and idx.size < M1.size:
            Pwarn(f"{M1.size - idx.size} of {M1.size} primaries cannot have "
                  "an inner companion (M1 <= M_min or no inner-binary density "
                  f"in {self.orbital_period_min:g}-"
                  f"{self.orbital_period_max:g} d) and are returned as single "
                  "stars despite all_binaries=True.",
                  "InappropriateValueWarning")
        m_bin = M1[idx]
        M2[idx] = m_bin * self._sample_q(m_bin, lp, rng)
        P[idx] = 10**lp
        e[idx] = self._sample_e(m_bin, lp, rng)
        return M2, P, e
