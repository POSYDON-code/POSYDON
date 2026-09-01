"""Maltsev+25 M_CO-based core-collapse supernova prescription.

This module implements the rapid binary-population-synthesis (BPS) CCSN recipe
of Maltsev et al. 2025 (arXiv:2503.23856, Sects. 3.2.1-3.2.3, Eq. 11). The
recipe predicts the *compact-object type* (NS, fallback BH or direct-collapse
BH) from the carbon-oxygen core mass ``M_CO``, the metallicity ``Z`` and the
mass-transfer (MT) history class of the progenitor. The class is taken from
the **first** mass-transfer episode (Maltsev+25, Appendix A.5.1); see
:meth:`Maltsev25_MCO_corecollapse._resolve_mt_class`.
Crucially, the recipe **separates** two distinct questions:

1. *Explodability* (does the star explode?): a deterministic decision based on
   the three ``M_CO`` boundaries ``M1 <= M_CO <= M2`` (direct collapse) and
   ``M_CO >= M3`` (direct collapse), with successful SNe elsewhere. Implemented
   in :meth:`Maltsev25_MCO_corecollapse.explodability`.
2. *Categorisation* (NS versus fallback BH): for a successful SN, a probabilistic
   decision inside the ``(M2, M3)`` window, using a guaranteed-NS sub-window.
   Implemented in :meth:`Maltsev25_MCO_corecollapse.categorisation`.

The recipe only returns the type; remnant *masses* are not part of
it. The class default therefore contains a NS mass scheme (see the
``NS_mass_model`` argument of :meth:`__init__`).

Extrapolation modes
-------------------
The recipe is calibrated for ``Z/Z_sun in (0.1, 1)``.  Outside this range the
boundaries must be extrapolated.  Three schemes are available, following
Willcox et al. 2025 (arXiv:2510.07573, Sec. 3.1.1):

* ``'optimistic'`` (default) -- linear extrapolation in ``log10(Z)`` extending
  the trend between the calibration boundaries, continued to arbitrary low/high
  ``Z``.
* ``'pessimistic'`` -- nearest-neighbour (clamp) at both calibration
  boundaries; the boundary values at ``Z/Z_sun = 0.1`` and ``Z/Z_sun = 1``
  are held constant outside the calibrated range.
* ``'balanced'`` -- linear extrapolation in ``log10(Z)`` down to
  ``Z/Z_sun = 1/50`` (~ 0.02, inspired by I Zwicky 18), then
  nearest-neighbour below that floor.  Above ``Z/Z_sun = 1`` the behaviour
  is the same as the optimistic mode (linear extrapolation).

"""

import numpy as np

from posydon.utils.posydonerror import ModelError

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

#   M_i(Z) / M_sun = a_i + b_i * log10(Z / Z_sun)
# for the three direct-collapse boundaries M1, M2, M3 (Tables 3 & 4 of the
# paper). Index 0 -> M1, 1 -> M2, 2 -> M3.
_BOUNDARIES = {
    'single': [(6.6, 0.5), (7.2, 0.6), (13.0, 0.1)],
    'case_A': [(7.4, 0.5), (8.4, 1.0), (15.4, 1.7)],
    'case_B': [(7.7, 0.7), (8.3, 0.4), (15.2, 1.2)],
    'case_C': [(6.6, 0.3), (7.1, 0.0), (13.2, 0.9)],
}

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

# Valid extrapolation modes.
EXTRAPOLATION_MODES = ('optimistic', 'balanced', 'pessimistic')

# Lower metallicity floor for the 'balanced' extrapolation mode
# (Z/Z_sun = 1/50 ~ 0.02, inspired by I Zwicky 18).
Z_BALANCED_FLOOR = 1.0 / 50.0

# Calibration boundaries of the Maltsev+25 recipe (Z/Z_sun).
Z_CALIBRATION_LOW = 0.1
Z_CALIBRATION_HIGH = 1.0

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
    compact-object type. Remnant masses are assigned by separately by the NS_mass_model.

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
    extrapolation_mode : {'optimistic', 'balanced', 'pessimistic'}
        How to extrapolate the ``M_CO`` boundaries outside the calibrated
        metallicity range ``Z/Z_sun in [0.1, 1]``.  See the module docstring
        for details.  Default ``'optimistic'``.
    verbose : bool
        Verbosity flag.

    """

    def __init__(self, RNG=None, NS_mass_model=None, NS_mass=1.4,
                 fallback_fraction=0.99, fallback_model='A',
                 extrapolation_mode='optimistic',
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

        if extrapolation_mode not in EXTRAPOLATION_MODES:
            raise ValueError(
                "extrapolation_mode must be one of %s, got '%s'."
                % (EXTRAPOLATION_MODES, extrapolation_mode))
        self.extrapolation_mode = extrapolation_mode

        self.verbose = verbose

    # ------------------------------------------------------------------
    # Built-in (swappable) NS mass scheme
    # ------------------------------------------------------------------
    def _default_NS_mass(self, star):
        """Built-in NS mass scheme: constant baryonic mass."""
        return self.NS_mass

    # ------------------------------------------------------------------
    # MT-class resolution
    # ------------------------------------------------------------------
    def _resolve_mt_class(self, star, default='single'):
        """Return the Maltsev+25 MT class of the collapsing star.

        The MT class determines which set of ``M_CO`` boundaries is used by
        the recipe, derived from the **first** mass-transfer episode of the
        binary (never a later one). This follows Maltsev+25, Appendix A.5.1:
        for repeated episodes the first episode is sufficient -- Case A/B
        pre-SN properties are similar, and for Case BC donors the first
        (Case B) episode is preferred because Case C progenitors are closer
        to single stars, so the Case B critical values capture the binary
        interaction effects more adequately.
        Since each grid (step_MESA) overwrites the class, the most recent grid
        wins. Values that are not a valid MT class (e.g. no MT interaction
        occurred, or an interpolator without a ``first_mt_case`` key) fall
        back to the ``default`` class.

        Parameters
        ----------
        star : object or None
            The collapsing star (``None`` in degenerate cases).
        default : str
            MT class to fall back to when ``star`` carries no valid class.
            One of 'single', 'case_A', 'case_B', 'case_C'.

        Returns
        -------
        str
            One of 'single', 'case_A', 'case_B', 'case_C'.

        """
        if star is None:
            return default
        mt_class = getattr(star, 'first_mt_class', default)
        if mt_class not in MT_CLASSES:
            # nothing valid stored (e.g. 'no_RLOF', 'initial_RLOF', 'None')
            return default
        return mt_class

    # ------------------------------------------------------------------
    # Boundary / window accessors
    # ------------------------------------------------------------------
    def _require_finite_M_CO(self, M_CO):
        """Raise ``ModelError`` if the CO core mass cannot feed the recipe."""
        if M_CO is None or not np.isfinite(M_CO):
            raise ModelError(
                "Maltsev+25 MCO prescription requires a finite CO core mass; "
                "got co_core_mass = %s." % M_CO)

    def _eval_boundary(self, a, b, Z):
        """Evaluate a single linear boundary ``a + b*log10(Z)`` with the
        configured extrapolation mode.

        See the module docstring (Extrapolation modes) for the behaviour of
        each mode outside the calibrated ``Z/Z_sun in [0.1, 1]`` range.

        Parameters
        ----------
        a, b : float
            Intercept and slope of the linear boundary fit.
        Z : float
            Metallicity ``Z/Z_sun``.

        Returns
        -------
        float
            The boundary value in Msun at ``Z``.

        """
        # Within the calibrated range the value is always the straight fit.
        if Z_CALIBRATION_LOW <= Z <= Z_CALIBRATION_HIGH:
            return a + b * np.log10(Z)

        if Z < Z_CALIBRATION_LOW:
            if self.extrapolation_mode == 'optimistic':
                # linear continuation all the way down
                return a + b * np.log10(Z)
            # clamp threshold for 'balanced' (floor) or 'pessimistic' (lower
            # calibration boundary).
            if self.extrapolation_mode == 'balanced':
                return a + b * np.log10(max(Z, Z_BALANCED_FLOOR))
            # pessimistic: nearest-neighbour at the lower calibration boundary
            return a + b * np.log10(Z_CALIBRATION_LOW)

        # Z > Z_CALIBRATION_HIGH : linear extrapolation in all modes
        return a + b * np.log10(Z)

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
        if Z is None or not np.isfinite(Z) or Z <= 0.0:
            raise ModelError(
                "Maltsev+25 MCO boundaries require a positive, finite "
                "Z/Z_sun; got Z = %s." % Z)
        return tuple(self._eval_boundary(a, b, Z)
                     for (a, b) in _BOUNDARIES[mt_class])

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
        if Z is None or not np.isfinite(Z) or Z <= 0.0:
            raise ModelError(
                "Maltsev+25 MCO boundaries require a positive, finite "
                "Z/Z_sun; got Z = %s." % Z)
        return tuple(self._eval_boundary(a, b, Z)
                     for (a, b) in _NS_WINDOW[mt_class])

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
        self._require_finite_M_CO(M_CO)
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
        self._require_finite_M_CO(M_CO)
        M1, _, _ = self.get_boundaries(mt_class, Z)
        # M_CO < M1: successful SN that forms an NS only.
        if M_CO < M1:
            return 'NS'

        # Guaranteed NS window
        NS1, NS2 = self.get_NS_window(mt_class, Z)
        if NS1 < M_CO < NS2:
            return 'NS'

        # fallback BH window
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
            and ``metallicity`` (Z/Z_sun). Its ``first_mt_class`` attribute (set by
            step_MESA) is used as the MT-history class when valid.
        mt_class : str
            Fallback MT-history class, used only when ``star.first_mt_class`` is
            missing or not a valid class: 'single', 'case_A', 'case_B' or
            'case_C'. When ``star.first_mt_class`` is valid it takes precedence.
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

        # get the MT class of the collapsing star
        mt_class = self._resolve_mt_class(star, mt_class)

        # main properties for prescription
        M_CO = star.co_core_mass
        Z = star.metallicity

        # step 1: determine whether the star explodes
        explodes = self.explodability(M_CO, Z, mt_class)

        # step 2: determine the outcome of the explosion
        if explodes:
            outcome = self.categorisation(M_CO, Z, mt_class)
        else:
            outcome = 'direct_BH'

        # step 3: compute the remnant mass and fallback fraction
        m_rembar, f_fb = self.remnant_mass_and_fallback(
            outcome, star, conserve_hydrogen_envelope)

        # clean up the outcome state
        state = 'BH' if outcome.endswith('BH') else 'NS'

        # preserve the separated categorisation
        # for example for successful or failed SN determination
        star.SN_categorisation = outcome

        return m_rembar, f_fb, state
