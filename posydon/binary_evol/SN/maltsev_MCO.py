"""Maltsev+25 M_CO-based core-collapse supernova prescription.

This module implements the rapid binary-population-synthesis (BPS) CCSN recipe
of Maltsev et al. 2025 (arXiv:2503.23856, Sects. 3.2.1-3.2.3, Eq. 11). The
recipe predicts the *compact-object type* (NS, fallback BH or direct-collapse
BH) from the carbon-oxygen core mass ``M_CO``, the metallicity ``Z`` and the
mass-transfer (MT) history class of the progenitor (single / case_A / case_B
[Be+Bl merged] / case_C).

Crucially, the recipe **separates** two distinct questions:

1. *Explodability* (does the star explode?): a deterministic decision based on
   the three ``M_CO`` boundaries ``M1 <= M_CO <= M2`` (direct collapse) and
   ``M_CO >= M3`` (direct collapse), with successful SNe elsewhere. Implemented
   in :meth:`Maltsev25_MCO_corecollapse.explodability`.
2. *Categorisation* (NS versus fallback BH): for a successful SN, a probabilistic
   decision inside the ``(M2, M3)`` window, using a guaranteed-NS sub-window.
   Implemented in :meth:`Maltsev25_MCO_corecollapse.categorisation`.

The recipe as published returns only the type; remnant *masses* are not part of
it. The class therefore ships with a simple, **swappable** mass scheme (see the
``NS_mass_model`` argument of :meth:`__init__`).

Note
----
The recipe is calibrated for ``Z/Z_sun in (0.1, 1)``. Outside this range the
boundaries are linearly extrapolated and a warning is emitted.

"""

import numpy as np

from posydon.utils.constants import Zsun as Z_SUN
from posydon.utils.posydonerror import ModelError
from posydon.utils.posydonwarning import Pwarn

# Solar metallicity in absolute units (Asplund+09, POSYDON convention). Used to
# normalise ``star.metallicity`` (which holds the *absolute* metallicity Z during
# binary evolution, i.e. MESA ``initial_z``) to the Z/Z_sun ratio required by the
# recipe. A directly constructed star may instead hold the Z/Z_sun ratio; the
# ``__call__`` method detects this and normalises accordingly.

# Coefficients (a, b) of the linear fits
#   M_i(Z) / M_sun = a_i + b_i * log10(Z / Z_sun)
# for the three direct-collapse boundaries M1, M2, M3 (Tables 3 & 4 of the
# paper). Index 0 -> M1, 1 -> M2, 2 -> M3.
_BOUNDARIES = {
    'single': [(6.6, 0.5), (7.2, 0.6), (13.0, 0.1)],
    'case_A': [(7.4, 0.5), (8.4, 1.0), (15.4, 1.7)],
    'case_B': [(7.7, 0.7), (8.3, 0.4), (15.2, 1.2)],
    'case_C': [(6.6, 0.3), (7.1, 0.0), (13.2, 0.9)],
}

# Coefficients (a, b) of the linear fits
#   NS_k(Z) / M_sun = a_k + b_k * log10(Z / Z_sun)
# for the guaranteed-NS sub-window (NS1, NS2) inside the successful-SN region
# (Tables 5 & 6 of the paper). Index 0 -> NS1, 1 -> NS2.
_NS_WINDOW = {
    'single': [(9.0, 1.6), (10.2, -0.8)],
    'case_A': [(11.1, 0.7), (12.1, 1.0)],
    'case_B': [(9.9, 0.6), (10.3, 0.0)],
    'case_C': [(9.6, 0.7), (10.7, 1.2)],
}

# Valid MT classes recognised by the recipe.
MT_CLASSES = tuple(_BOUNDARIES.keys())

# Default probabilities of forming a fallback BH (instead of an NS) for a
# successful SN that lies outside the guaranteed-NS sub-window.
_FALLBACK_PROB = {
    'A': 0.15,  # Section 3.2.2 of the paper (model A)
    'B': 0.10,  # Appendix A.6 / averaged 10% model (model B)
}


class Maltsev25_MCO_corecollapse(object):
    """Maltsev+25 M_CO-based CCSN recipe (rapid BPS variant).

    The class exposes the explodability and categorisation steps as separate
    methods and only needs ``M_CO``, ``Z`` and the MT class to predict the
    compact-object type. Remnant masses are assigned by a simple, swappable
    scheme (see ``NS_mass_model``).

    Parameters
    ----------
    RNG : numpy.random.Generator, optional
        Random number generator used for the stochastic NS/fallback-BH split.
        If ``None``, a fresh generator is created.
    NS_mass_model : callable, optional
        Callable ``f(star) -> float`` returning the baryonic NS mass in Msun.
        If ``None``, a constant baryonic mass ``NS_mass`` is used.
    NS_mass : float
        Default baryonic NS mass (Msun) used by the built-in ``NS_mass_model``.
    fallback_fraction : float
        Fallback mass fraction ``f_fb`` assigned to a fallback BH (the collapsing
        mass is the remnant mass). Default 0.99.
    fallback_model : {'A', 'B'}
        Stochastic model for the NS/fallback-BH split outside the guaranteed-NS
        window. 'A' uses a 15% fallback probability, 'B' a uniform 10%.
    Zsun : float
        Solar metallicity (only kept for reference / potential overrides).
    verbose : bool
        Verbosity flag.

    """

    def __init__(self, RNG=None, NS_mass_model=None, NS_mass=1.4,
                 fallback_fraction=0.99, fallback_model='A', Zsun=Z_SUN,
                 verbose=False):
        """Initialize a Maltsev25_MCO_corecollapse instance."""
        if RNG is None:
            self.RNG = np.random.default_rng()
        else:
            self.RNG = RNG

        if NS_mass_model is None:
            self.NS_mass = NS_mass
            self.NS_mass_model = self._default_NS_mass
        else:
            if not callable(NS_mass_model):
                raise ValueError("NS_mass_model must be a callable f(star).")
            self.NS_mass = NS_mass
            self.NS_mass_model = NS_mass_model

        if not (0.0 <= fallback_fraction <= 1.0):
            raise ValueError("fallback_fraction must be in [0, 1].")
        self.fallback_fraction = fallback_fraction

        if fallback_model not in _FALLBACK_PROB:
            raise ValueError(
                "fallback_model must be one of %s." % list(_FALLBACK_PROB))
        self.fallback_model = fallback_model

        self.Zsun = Zsun
        self.verbose = verbose

    # ------------------------------------------------------------------
    # Built-in (swappable) NS mass scheme
    # ------------------------------------------------------------------
    def _default_NS_mass(self, star):
        """Built-in NS mass scheme: constant baryonic mass."""
        return self.NS_mass

    # ------------------------------------------------------------------
    # Boundary / window accessors
    # ------------------------------------------------------------------
    def get_boundaries(self, mt_class, Z):
        """Return the three M_CO direct-collapse boundaries for ``mt_class``.

        Parameters
        ----------
        mt_class : str
            One of 'single', 'case_A', 'case_B', 'case_C'.
        Z : float
            Metallicity ``Z/Z_sun`` of the progenitor.

        Returns
        -------
        (M1, M2, M3) : tuple of float
            Direct-collapse boundaries in Msun, evaluated at ``Z``.

        """
        if mt_class not in _BOUNDARIES:
            raise ValueError(
                "mt_class must be one of %s, got '%s'." % (MT_CLASSES, mt_class))
        x = np.log10(Z)
        return tuple(a + b * x for (a, b) in _BOUNDARIES[mt_class])

    def get_NS_window(self, mt_class, Z):
        """Return the guaranteed-NS sub-window (NS1, NS2) for ``mt_class``.

        Parameters
        ----------
        mt_class : str
            One of 'single', 'case_A', 'case_B', 'case_C'.
        Z : float
            Metallicity ``Z/Z_sun`` of the progenitor.

        Returns
        -------
        (NS1, NS2) : tuple of float
            Guaranteed-NS window boundaries in Msun, evaluated at ``Z``.

        """
        if mt_class not in _NS_WINDOW:
            raise ValueError(
                "mt_class must be one of %s, got '%s'." % (MT_CLASSES, mt_class))
        x = np.log10(Z)
        return tuple(a + b * x for (a, b) in _NS_WINDOW[mt_class])

    # ------------------------------------------------------------------
    # 1) Explodability (deterministic)
    # ------------------------------------------------------------------
    def explodability(self, M_CO, Z, mt_class):
        """Decide whether the star explodes (successful SN) or not.

        Parameters
        ----------
        M_CO : float
            Carbon-oxygen core mass in Msun.
        Z : float
            Metallicity ``Z/Z_sun``.
        mt_class : str
            One of 'single', 'case_A', 'case_B', 'case_C'.

        Returns
        -------
        bool
            ``True`` if the star undergoes a successful SN, ``False`` if it
            collapses directly to a BH (failed SN).

        """
        M1, M2, M3 = self.get_boundaries(mt_class, Z)
        if M_CO < M1:
            return True
        if M_CO <= M2:
            # M1 <= M_CO <= M2 : direct-collapse window (failed SN)
            return False
        if M_CO < M3:
            # M2 < M_CO < M3 : successful SN
            return True
        # M_CO >= M3 : direct collapse (failed SN)
        return False

    # ------------------------------------------------------------------
    # 2) Categorisation (probabilistic, only for successful SNe)
    # ------------------------------------------------------------------
    def categorisation(self, M_CO, Z, mt_class):
        """Decide the compact-object type for a *successful* SN.

        A successful SN inside the guaranteed-NS sub-window ``(NS1, NS2)``
        always forms an NS. Outside that window it forms an NS with probability
        ``1 - p_fallback`` and a fallback BH with probability ``p_fallback``
        (``p_fallback`` depends on ``fallback_model``).

        Parameters
        ----------
        M_CO : float
            Carbon-oxygen core mass in Msun.
        Z : float
            Metallicity ``Z/Z_sun``.
        mt_class : str
            One of 'single', 'case_A', 'case_B', 'case_C'.

        Returns
        -------
        str
            'NS' or 'fallback_BH'.

        """
        M1, _, _ = self.get_boundaries(mt_class, Z)
        # M_CO < M1: successful SN that forms an NS only. The paper treats the
        # fallback-BH channel as statistically insignificant in this low-M_CO
        # region (Maltsev et al. 2025, line 571), so no stochastic NS/fallback
        # split is applied here.
        if M_CO < M1:
            return 'NS'
        NS1, NS2 = self.get_NS_window(mt_class, Z)
        if NS1 < M_CO < NS2:
            return 'NS'
        p_fallback = _FALLBACK_PROB[self.fallback_model]
        if self.RNG.uniform(0.0, 1.0) <= p_fallback:
            return 'fallback_BH'
        return 'NS'

    # ------------------------------------------------------------------
    # Remnant mass / fallback assignment
    # ------------------------------------------------------------------
    def remnant_mass_and_fallback(self, outcome, star,
                                  conserve_hydrogen_envelope=False):
        """Assign the baryonic remnant mass and fallback fraction.

        Parameters
        ----------
        outcome : str
            One of 'NS', 'direct_BH', 'fallback_BH'.
        star : object
            Collapsing star (used for the NS mass model and the collapsing
            mass of BH remnants).
        conserve_hydrogen_envelope : bool
            If ``True``, the whole stellar mass collapses to the BH; otherwise
            only the He-core mass does.

        Returns
        -------
        (m_rembar, f_fb) : tuple of float
            Baryonic remnant mass (Msun) and fallback mass fraction.

        """
        if outcome == 'NS':
            return self.NS_mass_model(star), 0.0
        if outcome in ('direct_BH', 'fallback_BH'):
            if conserve_hydrogen_envelope:
                m_collapsing = star.mass
            else:
                m_collapsing = star.he_core_mass
            f_fb = 1.0 if outcome == 'direct_BH' else self.fallback_fraction
            return m_collapsing, f_fb
        raise ValueError("Unknown outcome '%s'." % outcome)

    # ------------------------------------------------------------------
    # Entry point mirroring the other engine classes
    # ------------------------------------------------------------------
    def __call__(self, star, mt_class='single',
                 conserve_hydrogen_envelope=False):
        """Compute the remnant type, mass and fallback for a collapsing star.

        Parameters
        ----------
        star : object
            Collapsing star object. Must expose ``co_core_mass`` (M_CO, Msun)
            and ``metallicity`` (Z/Z_sun).
        mt_class : str
            MT-history class: 'single', 'case_A', 'case_B' or 'case_C'.
        conserve_hydrogen_envelope : bool
            Whether to assume the hydrogen envelope is conserved in direct
            collapse to a BH.

        Returns
        -------
        m_rembar : float
            Baryonic remnant mass in Msun.
        f_fb : float
            Fallback mass fraction.
        state : str
            'NS' or 'BH'.

        """
        M_CO = star.co_core_mass
        if M_CO is None or (isinstance(M_CO, float) and np.isnan(M_CO)):
            raise ModelError(
                "The CO core mass (M_CO) is not available; the Maltsev+25-MCO "
                "recipe cannot be applied.")

        Z = star.metallicity
        if Z is None or Z <= 0.0:
            raise ModelError(
                "The metallicity Z is not available; the Maltsev+25-MCO recipe "
                "requires Z/Z_sun > 0.")
        # star.metallicity holds the relative metallicity Z/Z_sun.
        Z_rel = Z
        if not (0.1 <= Z_rel <= 1.0):
            Pwarn(
                "Z/Z_sun=%.4f is outside the (0.1, 1) calibration range of the "
                "Maltsev+25-MCO recipe; boundaries are linearly extrapolated."
                % Z_rel, "ApproximationWarning")

        explodes = self.explodability(M_CO, Z_rel, mt_class)
        if explodes:
            outcome = self.categorisation(M_CO, Z_rel, mt_class)
        else:
            outcome = 'direct_BH'

        m_rembar, f_fb = self.remnant_mass_and_fallback(
            outcome, star, conserve_hydrogen_envelope)
        state = 'BH' if outcome.endswith('BH') else 'NS'

        # preserve the separated categorisation for inspection
        star.SN_categorisation = outcome

        return m_rembar, f_fb, state
