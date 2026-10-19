"""Period and mass ratio distributions module for POSYDON.

All distributions follow the same convention for their densities:

- ``pdf(x)`` is the probability density per unit ``x`` (as in scipy.stats),
  i.e. it integrates to 1 over ``x``.
- ``pdf_per_dex(x)`` is the probability density per unit ``log10(x)``,
  evaluated at ``x``, i.e. it integrates to 1 over ``log10(x)``.

The two are related by ``pdf_per_dex(x) = pdf(x) * x * ln(10)``.
"""

__authors__ = [
    "Max Briel <max.briel@gmail.com>"
]
import numpy as np
from scipy.integrate import quad
from scipy.stats import truncnorm


def _per_dex_from_pdf(x, pdf):
    """Convert a density per unit x into a density per unit log10(x).

    Parameters
    ----------
    x : float or array_like
        Value(s) the density was evaluated at.
    pdf : float or array_like
        Probability density per unit x at x.

    Returns
    -------
    float or ndarray
        Probability density per unit log10(x) at x.
    """
    x = np.asarray(x, dtype=float)
    return np.asarray(pdf, dtype=float) * x * np.log(10)


def _pdf_from_per_dex(x, per_dex):
    """Convert a density per unit log10(x) into a density per unit x.

    Parameters
    ----------
    x : float or array_like
        Value(s) the density was evaluated at.
    per_dex : float or array_like
        Probability density per unit log10(x) at x.

    Returns
    -------
    float or ndarray
        Probability density per unit x at x (0 where x <= 0).
    """
    x = np.asarray(x, dtype=float)
    per_dex = np.asarray(per_dex, dtype=float)
    return np.divide(per_dex, x * np.log(10),
                     out=np.zeros_like(per_dex), where=x > 0)


def _power_law_norm(slope, x_min, x_max):
    """Normalisation constant of x^slope on [x_min, x_max], with x_min > 0.

    Parameters
    ----------
    slope : float
        Power law exponent.
    x_min : float
        Lower bound.
    x_max : float
        Upper bound.

    Returns
    -------
    float
        1 / integral of x^slope from x_min to x_max.
    """
    if slope == -1:
        return 1.0 / np.log(x_max / x_min)
    return (slope + 1) / (x_max**(slope + 1) - x_min**(slope + 1))


def _power_law_inverse_cdf(u, slope, x_min, x_max):
    """Inverse CDF of x^slope on [x_min, x_max], with x_min > 0.

    Parameters
    ----------
    u : float or array_like
        Uniform random number(s) in [0, 1].
    slope : float
        Power law exponent.
    x_min : float
        Lower bound.
    x_max : float
        Upper bound.

    Returns
    -------
    float or ndarray
        Value(s) x with CDF(x) = u.
    """
    u = np.asarray(u, dtype=float)
    if slope == -1:
        return x_min * (x_max / x_min)**u
    k = slope + 1
    return (x_min**k + u * (x_max**k - x_min**k))**(1.0 / k)


class FlatMassRatio:
    """Flat mass ratio distribution for binary star systems.

    A uniform distribution for mass ratios q = m2/m1 within specified bounds.
    This distribution assigns equal probability to all mass ratios within the
    given range (q_min, q_max], exclusive bottom, inclusive top.

    Parameters
    ----------
    q_min : float, optional
        Minimum mass ratio (default: 0.05). Must be in [0, 1).
    q_max : float, optional
        Maximum mass ratio (default: 1.0). Must be in (0, 1].

    Raises
    ------
    ValueError
        If q_min or q_max are not in (0, 1], or if q_min >= q_max.

    Examples
    --------
    >>> dist = FlatMassRatio(q_min=0.1, q_max=0.8)
    >>> pdf_value = dist.pdf(0.5)
    """

    def __init__(self, q_min=0.05, q_max=1):
        """Initialize the flat mass ratio distribution.

        Parameters
        ----------
        q_min : float, optional
            Minimum mass ratio (default: 0.05). Must be in [0, 1).
        q_max : float, optional
            Maximum mass ratio (default: 1.0). Must be in (0, 1].

        Raises
        ------
        ValueError
            If q_min or q_max are not in valid range, or if q_min >= q_max.
        """
        if not (0 <= q_min < 1):
            raise ValueError("q_min must be in [0, 1)")
        if not (0 < q_max <= 1):
            raise ValueError("q_max must be in (0, 1]")
        if q_min >= q_max:
            raise ValueError("q_min must be less than q_max")

        self.q_min = q_min
        self.q_max = q_max
        self.norm = self._calculate_normalization()

    def __repr__(self):
        """Return string representation of the distribution.

        Returns
        -------
        str
            String representation showing the distribution parameters.
        """
        return f"FlatMassRatio(q_min={self.q_min}, q_max={self.q_max})"

    def _repr_html_(self):
        """Return HTML representation for Jupyter notebooks.

        Returns
        -------
        str
            HTML string for rich display in notebooks.
        """
        return (f"<h3>Flat Mass Ratio Distribution</h3>"
                f"<p>q_min = {self.q_min}</p>"
                f"<p>q_max = {self.q_max}</p>")

    def _calculate_normalization(self):
        """Calculate the normalization constant for the flat mass ratio
        distribution.

        Returns
        -------
        float
            The normalization constant ensuring the PDF integrates to 1.
        """
        integral, _ = quad(self.flat_mass_ratio, self.q_min, self.q_max)
        if (integral == 0): # pragma: no cover
            raise ValueError("Normalization integral is zero. "
                             "Check mass ratio parameters.")
        return 1.0 / integral

    def flat_mass_ratio(self, q):
        """Compute the flat mass ratio distribution at a given mass.

        Parameters
        ----------
        q : float
            Mass ratio.

        Returns
        -------
        float
            Distribution value (always 1.0 for flat distribution).
        """
        return 1.0

    def pdf(self, q):
        """Probability density function of the flat mass ratio distribution.

        Parameters
        ----------
        q : float or array_like
            Mass ratio(s).

        Returns
        -------
        float or ndarray
            Probability density at mass ratio q.
        """
        q = np.asarray(q)
        valid = (q > self.q_min) & (q <= self.q_max)
        pdf_values = np.zeros_like(q, dtype=float)
        pdf_values[valid] = self.flat_mass_ratio(q[valid]) * self.norm
        return pdf_values

    def pdf_per_dex(self, q):
        """Probability density per unit log10(q) of the flat mass ratio distribution.

        Normalised in log10(q) space to 1.

        Parameters
        ----------
        q : float or array_like
            Mass ratio(s).

        Returns
        -------
        float or ndarray
            Probability density per unit log10(q) at q.
        """
        return _per_dex_from_pdf(q, self.pdf(q))

    def rvs(self, size=1, rng=None):
        """Draw random samples from the flat mass ratio distribution.

        Parameters
        ----------
        size : int, optional
            Number of samples to draw (default: 1).
        rng : numpy.random.Generator, optional
            Random number generator. If None, uses np.random.default_rng().

        Returns
        -------
        float or ndarray
            Random samples from the distribution.
        """
        if rng is None:
            rng = np.random.default_rng()

        return rng.uniform(self.q_min, self.q_max, size=size)


class PowerLawMassRatio:
    """Power law mass ratio distribution for binary star systems.

    A distribution where the PDF follows q^alpha within specified bounds
    (q_min, q_max], exclusive bottom, inclusive top.

    Parameters
    ----------
    alpha : float, optional
        Power law exponent (default: 0.0, i.e. flat). Can be any real number.
    q_min : float, optional
        Minimum mass ratio (default: 0.05). Must be in [0, 1).
    q_max : float, optional
        Maximum mass ratio (default: 1.0). Must be in (0, 1].

    Raises
    ------
    ValueError
        If q_min or q_max are not in valid range, or if q_min >= q_max.

    Examples
    --------
    >>> dist = PowerLawMassRatio(alpha=-1.0, q_min=0.1, q_max=1.0)
    >>> pdf_value = dist.pdf(0.5)
    """

    def __init__(self, alpha=0.0, q_min=0.05, q_max=1.0):
        """Initialize the power law mass ratio distribution.

        Parameters
        ----------
        alpha : float, optional
            Power law exponent (default: 0.0).
        q_min : float, optional
            Minimum mass ratio (default: 0.05). Must be in [0, 1).
        q_max : float, optional
            Maximum mass ratio (default: 1.0). Must be in (0, 1].

        Raises
        ------
        ValueError
            If q_min or q_max are not in valid range, or if q_min >= q_max.
        """
        if not (0 <= q_min < 1):
            raise ValueError("q_min must be in [0, 1)")
        if not (0 < q_max <= 1):
            raise ValueError("q_max must be in (0, 1]")
        if q_min >= q_max:
            raise ValueError("q_min must be less than q_max")
        if alpha <= -1 and q_min == 0:
            raise ValueError("q_min must be > 0 for alpha <= -1 "
                             "to avoid divergent integral")

        self.alpha = alpha
        self.q_min = q_min
        self.q_max = q_max
        self.norm = self._calculate_normalization()

    def __repr__(self):
        """Return string representation of the distribution."""
        return (f"PowerLawMassRatio(alpha={self.alpha}, "
                f"q_min={self.q_min}, q_max={self.q_max})")

    def _repr_html_(self):
        """Return HTML representation for Jupyter notebooks."""
        return (f"<h3>Power Law Mass Ratio Distribution</h3>"
                f"<p>alpha = {self.alpha}</p>"
                f"<p>q_min = {self.q_min}</p>"
                f"<p>q_max = {self.q_max}</p>")

    def _calculate_normalization(self):
        """Calculate the normalization constant for the power law mass ratio
        distribution.

        Returns
        -------
        float
            The normalization constant ensuring the PDF integrates to 1.
        """
        integral, _ = quad(self.power_law_mass_ratio, self.q_min, self.q_max)
        if integral == 0:  # pragma: no cover
            raise ValueError("Normalization integral is zero. "
                             "Check mass ratio parameters.")
        return 1.0 / integral

    def power_law_mass_ratio(self, q):
        """Compute the power law mass ratio distribution value.

        Parameters
        ----------
        q : float or array_like
            Mass ratio.

        Returns
        -------
        float or ndarray
            Distribution value q^alpha.
        """
        return np.asarray(q, dtype=float)**self.alpha

    def pdf(self, q):
        """Probability density function of the power law mass ratio distribution.

        Parameters
        ----------
        q : float or array_like
            Mass ratio(s).

        Returns
        -------
        float or ndarray
            Probability density at mass ratio q.
        """
        q = np.asarray(q)
        valid = (q > self.q_min) & (q <= self.q_max)
        pdf_values = np.zeros_like(q, dtype=float)
        pdf_values[valid] = self.power_law_mass_ratio(q[valid]) * self.norm
        return pdf_values

    def pdf_per_dex(self, q):
        """Probability density per unit log10(q) of the power law mass ratio distribution.

        Normalised in log10(q) space to 1.

        Parameters
        ----------
        q : float or array_like
            Mass ratio(s).

        Returns
        -------
        float or ndarray
            Probability density per unit log10(q) at q.
        """
        return _per_dex_from_pdf(q, self.pdf(q))

    def rvs(self, size=1, n_points=1000, rng=None):
        """Draw random samples from the power law mass ratio distribution.

        Parameters
        ----------
        size : int, optional
            Number of samples to draw (default: 1).
        rng : numpy.random.Generator, optional
            Random number generator. If None, uses np.random.default_rng().

        Returns
        -------
        ndarray
            Random samples from the distribution.
        """
        if rng is None:
            rng = np.random.default_rng()

        from posydon.utils.common_functions import inverse_sampler
        q_grid = np.linspace(self.q_min, self.q_max, n_points)
        pdf_values = self.power_law_mass_ratio(q_grid)
        return inverse_sampler(q_grid, pdf_values, size=size, rng=rng)


class Sana12Period():
    """Period distribution from Sana et al. (2012).

    Orbital period distribution for massive binary stars based on the
    observational study by Sana H., et al., 2012, Science, 337, 444.
    The distribution has different behaviors for low-mass (≤15 M☉) and
    high-mass (>15 M☉) primary stars.

    For low-mass primaries: flat distribution in log period
    For high-mass primaries: power law with slope -0.55 and a break at log P = 0.15

    Parameters
    ----------
    p_min : float, optional
        Minimum period in days (default: 0.35).
    p_max : float, optional
        Maximum period in days (default: 6000).

    Attributes
    ----------
    mbreak : float
        Mass break at 15 M☉ separating low and high mass regimes.
    slope : float
        Power law slope for high-mass regime (-0.55).

    References
    ----------
    Sana H., et al., 2012, Science, 337, 444
    Bavera S., et al., 2021, A&A, 647, A153
    (https://ui.adsabs.harvard.edu/abs/2021A%26A...647A.153B/abstract)
    """

    def __init__(self, p_min=0.35, p_max=6e3):
        """Initialize the Sana12 period distribution.

        Parameters
        ----------
        p_min : float, optional
            Minimum period in days (default: 0.35).
        p_max : float, optional
            Maximum period in days (default: 6000).

        Raises
        ------
        ValueError
            If p_min is not positive or p_max <= p_min.
        """
        if p_min <= 0:
            raise ValueError("p_min must be positive")
        if p_max <= p_min:
            raise ValueError("p_max must be greater than p_min")
        self.p_min = p_min
        self.p_max = p_max
        self.mbreak = 15
        self.slope = 0.55  # Set slope before normalization calculations
        # below this log10(P) the high-mass distribution is flat
        self.logp_break = 0.15

        # mass boundary is at 15 Msun
        self.low_mass_norm = self._calculate_normalization(self.mbreak-1)
        self.high_mass_norm = self._calculate_normalization(self.mbreak+1)
        self.norm = lambda m1: np.where(m1 <= self.mbreak,
                                        self.low_mass_norm,
                                        self.high_mass_norm)

    def __repr__(self):
        """Return string representation of the distribution."""
        return f"Sana12Period(p_min={self.p_min}, p_max={self.p_max})"

    def _repr_html_(self):
        """Return HTML representation for Jupyter notebooks."""
        return (f"<h3>Sana12 Period Distribution</h3>"
                f"<p>p_min = {self.p_min}</p>"
                f"<p>p_max = {self.p_max}</p>")

    def _calculate_normalization(self, m1):
        """Calculate the normalization constant for the Sana12 period
        distribution.

        Returns
        -------
        float
            The normalization constant ensuring the PDF integrates to 1.
        """
        logp_min, logp_max = np.log10(self.p_min), np.log10(self.p_max)
        # integrate the two pieces separately around the break
        breaks = ([self.logp_break]
                  if logp_min < self.logp_break < logp_max else None)
        integral = quad(self.sana12_period, logp_min, logp_max,
                        args=(m1,), points=breaks)[0]
        if integral == 0: # pragma: no cover
            raise ValueError("Normalization integral is zero. "
                             "Check period parameters.")
        return 1.0 / integral

    def sana12_period(self, logp, m1):
        """Compute the Sana12 period distribution at given periods.

        Parameters
        ----------
        logp : float or array
            Log10 of period(s).
        m1 : float or array
            Mass value(s). Single value or an array of the same length as p.

        Returns
        -------
        float or ndarray
            Distribution value(s).
        """
        if np.any(m1 <= 0):
            raise ValueError("Mass must be positive. "
                             "Check parameters.")

        logp = np.asarray(logp, dtype=float)
        m1 = np.asarray(m1)
        if m1.size == 1:
            m1 = np.full_like(logp, m1)
        elif m1.shape != logp.shape:
            raise ValueError("m1 must be a single value or an array of "
                             "the same length as p.")

        result = np.zeros_like(logp)

        # low-mass primaries: flat in log10(P)
        mask1 = (m1 > 0) & (m1 <= self.mbreak)
        result[mask1] = 1

        # high-mass primaries: (log10 P)^-slope above the break, and constant
        # below it, so the distribution is continuous at the break.
        # The normalisation is done in `_calculate_normalization`, which keeps
        # this well-defined for any period range (also p_max < 10**logp_break).
        mask2 = m1 > self.mbreak
        logp_high = logp[mask2]
        in_range = ((logp_high >= np.log10(self.p_min))
                    & (logp_high <= np.log10(self.p_max)))
        temp_result = np.zeros_like(logp_high)
        temp_result[in_range] = np.maximum(logp_high[in_range],
                                           self.logp_break)**(-self.slope)
        result[mask2] = temp_result

        return result


    def pdf(self, p, m1):
        '''Probability density function of the Sana12 period distribution.
        Conditional on m1.

        Normalised in P space to 1. See `pdf_per_dex` for the density
        per unit log10(P).

        Parameters
        ----------
        p : float or array_like
            Period(s).
        m1 : float or array_like
            Mass value(s). Single value or an array of the same length as p.

        Returns
        -------
        float or ndarray
            Probability density per unit P at period p.

        '''
        return _pdf_from_per_dex(p, self.pdf_per_dex(p, m1))

    def pdf_per_dex(self, p, m1):
        '''Probability density per unit log10(P) of the Sana12 period
        distribution. Conditional on m1.

        Normalised in log10(P) space to 1.

        Parameters
        ----------
        p : float or array_like
            Period(s).
        m1 : float or array_like
            Mass value(s). Single value or an array of the same length as p.

        Returns
        -------
        float or ndarray
            Probability density per unit log10(P) at period p.

        '''
        p = np.asarray(p, dtype=float)
        m1 = np.asarray(m1)
        if m1.size == 1:
            m1 = np.full_like(p, m1)
        elif m1.shape != p.shape:
            raise ValueError("m1 must be a single value or an array of "
                             "the same length as p.")

        valid = (p >= self.p_min) & (p <= self.p_max)
        pdf_values = np.zeros_like(p, dtype=float)

        logp = np.log10(p[valid])
        pdf_values[valid] = (self.sana12_period(logp, m1[valid])
                             * self.norm(m1[valid]))
        return pdf_values

    def rvs(self, size=1, m1=None, rng=None):
        """Draw random samples from the Sana12 period distribution.

        Parameters
        ----------
        size : int, optional
            Number of samples to draw (default: 1).
        m1 : float or array_like, optional
            Primary mass(es). If array, must have length equal to size.
        rng : numpy.random.Generator, optional
            Random number generator. If None, uses np.random.default_rng().

        Returns
        -------
        ndarray
            Random period samples in days.

        Raises
        ------
        ValueError
            If m1 is None or has incorrect size.
        """
        if rng is None:
            rng = np.random.default_rng()

        if m1 is None:
            raise ValueError("m1 (primary mass) must be provided for Sana12Period sampling")

        m1 = np.atleast_1d(m1)
        if m1.size == 1:
            m1 = np.full(size, m1[0])
        elif m1.size != size:
            raise ValueError(f"m1 must be a single value or have size={size}"
                                      f"\n m1 = {m1}\n m1.size = {m1.size}")

        # Import here to avoid circular dependency
        from posydon.utils.common_functions import rejection_sampler

        periods = np.zeros(size)

        # For low mass stars (m1 <= 15), use log-uniform distribution
        low_mass_mask = m1 <= self.mbreak
        n_low = np.sum(low_mass_mask)
        if n_low > 0:
            periods[low_mass_mask] = 10**rng.uniform(
                np.log10(self.p_min),
                np.log10(self.p_max),
                size=n_low
            )

        # For high mass stars (m1 > 15), use rejection sampling
        high_mass_mask = ~low_mass_mask
        n_high = np.sum(high_mass_mask)
        if n_high > 0:
            # Create PDF function for rejection sampler
            def pdf_high_mass(logp):
                return self.sana12_period(logp, self.mbreak + 1)

            periods[high_mass_mask] = 10**rejection_sampler(
                size=n_high,
                x_lim=[np.log10(self.p_min), np.log10(self.p_max)],
                pdf=pdf_high_mass,
                rng=rng
            )

        return periods


class PowerLawPeriod():
    '''Power law period distribution in log10(P), as in Sana et al. (2012).

    The density per unit log10(P) follows

        f(log10 P) ∝ (log10 P)^slope

    between p_min and p_max. With the default slope of -0.55 this is the
    Sana et al. (2012) period distribution extended to a wider period range.
    See `LinearPowerLawPeriod` for a power law in P itself.

    Normalised in log10(P) space to 1.

    Parameters
    ----------
    slope : float
        Exponent of log10(P) (default: -0.55).
    p_min : float
        Minimum period in days (default: 1.4). Must be greater than 1 day,
        since (log10 P)^slope is undefined for log10(P) <= 0.
    p_max : float
        Maximum period in days (default: 3e6).

    References
    ----------
    Sana H., et al., 2012, Science, 337, 444
    '''

    def __init__(self, slope=-0.55, p_min=1.4, p_max=3e6):
        """Initialize the power law period distribution.

        Parameters
        ----------
        slope : float, optional
            Exponent of log10(P) (default: -0.55).
        p_min : float, optional
            Minimum period in days (default: 1.4). Must be greater than 1.
        p_max : float, optional
            Maximum period in days (default: 3e6).

        Raises
        ------
        ValueError
            If p_min is not greater than 1 day or p_max <= p_min.
        """
        if p_min <= 1:
            raise ValueError("p_min must be greater than 1 day, since "
                             "(log10 P)^slope is undefined for P <= 1 day")
        if p_max <= p_min:
            raise ValueError("p_max must be greater than p_min")
        self.slope = slope
        self.p_min = p_min
        self.p_max = p_max

        self.norm = self._calculate_normalization()

    def __repr__(self):
        """Return string representation of the distribution.

        Returns
        -------
        str
            String representation showing the distribution parameters.
        """
        return f"PowerLawPeriod(slope={self.slope}, p_min={self.p_min}, p_max={self.p_max})"

    def _repr_html_(self):
        """Return HTML representation for Jupyter notebooks.

        Returns
        -------
        str
            HTML string for rich display in notebooks.
        """
        return (f"<h3>Power Law Period Distribution</h3>"
                f"<p>slope = {self.slope}</p>"
                f"<p>p_min = {self.p_min}</p>"
                f"<p>p_max = {self.p_max}</p>")

    def _calculate_normalization(self):
        """Calculate the normalization constant for the power law period
        distribution.

        Returns
        -------
        float
            The normalization constant ensuring the PDF integrates to 1
            over log10(P).
        """
        return _power_law_norm(self.slope,
                               np.log10(self.p_min),
                               np.log10(self.p_max))

    def power_law_period(self, logp):
        """Compute the unnormalised power law period distribution.

        Parameters
        ----------
        logp : float or array
            Log10 of period(s).

        Returns
        -------
        float or ndarray
            Distribution value(s), (log10 P)^slope inside [p_min, p_max]
            and 0 outside.
        """
        logp = np.asarray(logp, dtype=float)

        result = np.zeros_like(logp)
        valid = ((logp >= np.log10(self.p_min))
                 & (logp <= np.log10(self.p_max)))
        result[valid] = logp[valid]**self.slope

        return result

    def pdf(self, p):
        """Probability density function of the power law period distribution.

        Normalised in P space to 1. See `pdf_per_dex` for the density
        per unit log10(P).

        Parameters
        ----------
        p : float or array_like
            Period(s).

        Returns
        -------
        float or ndarray
            Probability density per unit P at period p.
        """
        return _pdf_from_per_dex(p, self.pdf_per_dex(p))

    def pdf_per_dex(self, p):
        """Probability density per unit log10(P) of the power law period
        distribution.

        Normalised in log10(P) space to 1.

        Parameters
        ----------
        p : float or array_like
            Period(s).

        Returns
        -------
        float or ndarray
            Probability density per unit log10(P) at period p.
        """
        p = np.asarray(p, dtype=float)
        valid = (p >= self.p_min) & (p <= self.p_max)
        pdf_values = np.zeros_like(p)

        pdf_values[valid] = (self.power_law_period(np.log10(p[valid]))
                             * self.norm)

        return pdf_values

    def rvs(self, size=1, rng=None):
        """Draw random samples from the power law period distribution.

        Uses the analytical inverse CDF in log10(P).

        Parameters
        ----------
        size : int, optional
            Number of samples to draw (default: 1).
        rng : numpy.random.Generator, optional
            Random number generator. If None, uses np.random.default_rng().

        Returns
        -------
        ndarray
            Random period samples in days.
        """
        if rng is None:
            rng = np.random.default_rng()

        logp_samples = _power_law_inverse_cdf(rng.uniform(size=size),
                                              self.slope,
                                              np.log10(self.p_min),
                                              np.log10(self.p_max))
        return 10**logp_samples


class LinearPowerLawPeriod():
    '''Power law period distribution in P.

    The density per unit P follows

        f(P) ∝ P^slope

    between p_min and p_max. A slope of -1 is flat in log10(P) (Öpik's law).
    See `PowerLawPeriod` for a power law in log10(P).

    Normalised in P space to 1.

    Parameters
    ----------
    slope : float
        Exponent of P (default: -1.0).
    p_min : float
        Minimum period in days (default: 1.4).
    p_max : float
        Maximum period in days (default: 3e6).
    '''

    def __init__(self, slope=-1.0, p_min=1.4, p_max=3e6):
        """Initialize the linear power law period distribution.

        Parameters
        ----------
        slope : float, optional
            Exponent of P (default: -1.0).
        p_min : float, optional
            Minimum period in days (default: 1.4).
        p_max : float, optional
            Maximum period in days (default: 3e6).

        Raises
        ------
        ValueError
            If p_min is not positive or p_max <= p_min.
        """
        if p_min <= 0:
            raise ValueError("p_min must be positive")
        if p_max <= p_min:
            raise ValueError("p_max must be greater than p_min")
        self.slope = slope
        self.p_min = p_min
        self.p_max = p_max

        self.norm = _power_law_norm(self.slope, self.p_min, self.p_max)

    def __repr__(self):
        """Return string representation of the distribution.

        Returns
        -------
        str
            String representation showing the distribution parameters.
        """
        return (f"LinearPowerLawPeriod(slope={self.slope}, "
                f"p_min={self.p_min}, p_max={self.p_max})")

    def _repr_html_(self):
        """Return HTML representation for Jupyter notebooks.

        Returns
        -------
        str
            HTML string for rich display in notebooks.
        """
        return (f"<h3>Linear Power Law Period Distribution</h3>"
                f"<p>slope = {self.slope}</p>"
                f"<p>p_min = {self.p_min}</p>"
                f"<p>p_max = {self.p_max}</p>")

    def linear_power_law_period(self, p):
        """Compute the unnormalised linear power law period distribution.

        Parameters
        ----------
        p : float or array_like
            Period(s) in days.

        Returns
        -------
        float or ndarray
            Distribution value(s), P^slope inside [p_min, p_max] and 0
            outside.
        """
        p = np.asarray(p, dtype=float)

        result = np.zeros_like(p)
        valid = (p >= self.p_min) & (p <= self.p_max)
        result[valid] = p[valid]**self.slope

        return result

    def pdf(self, p):
        """Probability density function of the linear power law period
        distribution.

        Normalised in P space to 1.

        Parameters
        ----------
        p : float or array_like
            Period(s).

        Returns
        -------
        float or ndarray
            Probability density per unit P at period p.
        """
        return self.linear_power_law_period(p) * self.norm

    def pdf_per_dex(self, p):
        """Probability density per unit log10(P) of the linear power law
        period distribution.

        Normalised in log10(P) space to 1.

        Parameters
        ----------
        p : float or array_like
            Period(s).

        Returns
        -------
        float or ndarray
            Probability density per unit log10(P) at period p.
        """
        return _per_dex_from_pdf(p, self.pdf(p))

    def rvs(self, size=1, rng=None):
        """Draw random samples from the linear power law period distribution.

        Uses the analytical inverse CDF in P.

        Parameters
        ----------
        size : int, optional
            Number of samples to draw (default: 1).
        rng : numpy.random.Generator, optional
            Random number generator. If None, uses np.random.default_rng().

        Returns
        -------
        ndarray
            Random period samples in days.
        """
        if rng is None:
            rng = np.random.default_rng()

        return _power_law_inverse_cdf(rng.uniform(size=size),
                                      self.slope, self.p_min, self.p_max)


class LogUniform():
    """Log-uniform distribution between specified minimum and maximum values.
    """

    def __init__(self, min=5.0, max=1e5):
        """
        Initialize the log-uniform distribution.

        Parameters
        ----------
        min : float, optional
            Minimum value
        max : float, optional
            Maximum value

        Raises
        ------
        ValueError
            If min is not positive or max <= min.
        """
        if min <= 0:
            raise ValueError("min must be positive")
        if max <= min:
            raise ValueError("max must be greater than min")
        self.min = min
        self.max = max

        self.norm = self._calculate_normalization()

    def __repr__(self):
        """Return string representation of the distribution.

        Returns
        -------
        str
            String representation showing the distribution parameters.
        """
        return f"LogUniform(min={self.min}, max={self.max})"

    def _repr_html_(self):
        """Return HTML representation for Jupyter notebooks.

        Returns
        -------
        str
            HTML string for rich display in notebooks.
        """
        return (f"<h3>Log-Uniform Distribution</h3>"
                f"<p>min = {self.min}</p>"
                f"<p>max = {self.max}</p>")

    def _calculate_normalization(self):
        """
        Calculate the normalization constant for the log-uniform distribution.

        Returns
        -------
        float
            The normalization constant ensuring the PDF integrates to 1.
        """
        return 1.0 / (np.log10(self.max) - np.log10(self.min))

    def pdf(self, x):
        """
        Probability density function of the log-uniform distribution.

        Normalised in x space to 1, as in scipy.stats.loguniform.
        See `pdf_per_dex` for the density per unit log10(x).

        Parameters
        ----------
        x : float or array_like
            Value(s) to evaluate the PDF at.

        Returns
        -------
        float or ndarray
            Probability density per unit x at x.
        """
        return _pdf_from_per_dex(x, self.pdf_per_dex(x))

    def pdf_per_dex(self, x):
        """
        Probability density per unit log10(x) of the log-uniform distribution.

        Normalised in log10(x) space to 1.

        Parameters
        ----------
        x : float or array_like
            Value(s) to evaluate the density at.

        Returns
        -------
        float or ndarray
            Probability density per unit log10(x) at x.
        """
        x = np.asarray(x)
        valid = (x > 0) & (x >= self.min) & (x <= self.max)
        pdf_values = np.zeros_like(x, dtype=float)

        pdf_values[valid] = self.norm

        return pdf_values

    def rvs(self, size=1, rng=None):
        """Draw random samples from the log-uniform distribution.

        Parameters
        ----------
        size : int, optional
            Number of samples to draw (default: 1).
        rng : numpy.random.Generator, optional
            Random number generator. If None, uses np.random.default_rng().

        Returns
        -------
        ndarray
            Random samples from the distribution.
        """
        if rng is None:
            rng = np.random.default_rng()

        # Sample uniformly in log space
        log_samples = rng.uniform(np.log10(self.min), np.log10(self.max), size=size)

        # Convert back to linear space
        return 10**log_samples


class ThermalEccentricity:
    """Thermal eccentricity distribution for binary star systems.

    The thermal distribution follows pdf(e) = 2*e, which is the
    distribution expected for binaries that have undergone significant
    dynamical interactions or thermal relaxation.

    Parameters
    ----------
    e_min : float, optional
        Minimum eccentricity (default: 0.0). Must be in [0, 1).
    e_max : float, optional
        Maximum eccentricity (default: 1.0). Must be in (0, 1].

    Raises
    ------
    ValueError
        If e_min or e_max are not in [0, 1], or if e_min >= e_max.

    Examples
    --------
    >>> dist = ThermalEccentricity()
    >>> e_samples = dist.rvs(size=1000)
    """

    def __init__(self, e_min=0.0, e_max=1.0):
        """Initialize the thermal eccentricity distribution.

        Parameters
        ----------
        e_min : float, optional
            Minimum eccentricity (default: 0.0). Must be in [0, 1).
        e_max : float, optional
            Maximum eccentricity (default: 1.0). Must be in (0, 1].

        Raises
        ------
        ValueError
            If e_min or e_max are not in [0, 1], or if e_min >= e_max.
        """
        if not (0 <= e_min < 1):
            raise ValueError("e_min must be in [0, 1)")
        if not (0 < e_max <= 1):
            raise ValueError("e_max must be in (0, 1]")
        if e_min >= e_max:
            raise ValueError("e_min must be less than e_max")

        self.e_min = e_min
        self.e_max = e_max
        self.norm = self._calculate_normalization()

    def __repr__(self):
        """Return string representation of the distribution."""
        return f"ThermalEccentricity(e_min={self.e_min}, e_max={self.e_max})"

    def _repr_html_(self):
        """Return HTML representation for Jupyter notebooks."""
        return (f"<h3>Thermal Eccentricity Distribution</h3>"
                f"<p>e_min = {self.e_min}</p>"
                f"<p>e_max = {self.e_max}</p>")

    def _calculate_normalization(self):
        """Calculate the normalization constant for the thermal eccentricity
        distribution.

        Returns
        -------
        float
            The normalization constant ensuring the PDF integrates to 1.
        """
        # Integral of 2*e from e_min to e_max is e_max^2 - e_min^2
        integral = self.e_max**2 - self.e_min**2
        if integral == 0: # pragma: no cover
            raise ValueError("Cannot normalize distribution: e_min == e_max")
        return 1.0 / integral

    def thermal_eccentricity(self, e):
        """Compute the thermal eccentricity distribution value.

        Parameters
        ----------
        e : float or array_like
            Eccentricity value(s).

        Returns
        -------
        float or ndarray
            Distribution value (2*e).
        """
        return 2.0 * np.asarray(e)

    def pdf(self, e):
        """Probability density function of the thermal eccentricity distribution.

        Parameters
        ----------
        e : float or array_like
            Eccentricity value(s).

        Returns
        -------
        float or ndarray
            Probability density at eccentricity e.
        """
        e = np.asarray(e)
        valid = (e >= self.e_min) & (e <= self.e_max)
        pdf_values = np.zeros_like(e, dtype=float)
        pdf_values[valid] = self.thermal_eccentricity(e[valid]) * self.norm
        return pdf_values

    def pdf_per_dex(self, e):
        """Probability density per unit log10(e) of the thermal eccentricity distribution.

        Normalised in log10(e) space to 1.

        Parameters
        ----------
        e : float or array_like
            Eccentricity value(s).

        Returns
        -------
        float or ndarray
            Probability density per unit log10(e) at e.
        """
        return _per_dex_from_pdf(e, self.pdf(e))

    def rvs(self, size=1, rng=None):
        """Draw random samples from the thermal eccentricity distribution.

        Uses the analytical inverse CDF: e = sqrt(u * (e_max^2 - e_min^2) + e_min^2)

        Parameters
        ----------
        size : int, optional
            Number of samples to draw (default: 1).
        rng : numpy.random.Generator, optional
            Random number generator. If None, uses np.random.default_rng().

        Returns
        -------
        ndarray
            Random eccentricity samples in [e_min, e_max].
        """
        if rng is None:
            rng = np.random.default_rng()

        # Inverse CDF: e = sqrt(u * (e_max^2 - e_min^2) + e_min^2)
        u = rng.uniform(size=size)
        return np.sqrt(u * (self.e_max**2 - self.e_min**2) + self.e_min**2)


class UniformEccentricity:
    """Uniform eccentricity distribution for binary star systems.

    A flat distribution over eccentricities between e_min and e_max.

    Parameters
    ----------
    e_min : float, optional
        Minimum eccentricity (default: 0.0). Must be in [0, 1).
    e_max : float, optional
        Maximum eccentricity (default: 1.0). Must be in (0, 1].

    Raises
    ------
    ValueError
        If e_min or e_max are not in [0, 1], or if e_min >= e_max.

    Examples
    --------
    >>> dist = UniformEccentricity()
    >>> e_samples = dist.rvs(size=1000)
    """

    def __init__(self, e_min=0.0, e_max=1.0):
        """Initialize the uniform eccentricity distribution.

        Parameters
        ----------
        e_min : float, optional
            Minimum eccentricity (default: 0.0). Must be in [0, 1).
        e_max : float, optional
            Maximum eccentricity (default: 1.0). Must be in (0, 1].

        Raises
        ------
        ValueError
            If e_min or e_max are not in [0, 1], or if e_min >= e_max.
        """
        if not (0 <= e_min < 1):
            raise ValueError("e_min must be in [0, 1)")
        if not (0 < e_max <= 1):
            raise ValueError("e_max must be in (0, 1]")
        if e_min >= e_max:
            raise ValueError("e_min must be less than e_max")

        self.e_min = e_min
        self.e_max = e_max
        self.norm = 1.0 / (e_max - e_min)

    def __repr__(self):
        """Return string representation of the distribution."""
        return f"UniformEccentricity(e_min={self.e_min}, e_max={self.e_max})"

    def _repr_html_(self):
        """Return HTML representation for Jupyter notebooks."""
        return (f"<h3>Uniform Eccentricity Distribution</h3>"
                f"<p>e_min = {self.e_min}</p>"
                f"<p>e_max = {self.e_max}</p>")

    def pdf(self, e):
        """Probability density function of the uniform eccentricity distribution.

        Parameters
        ----------
        e : float or array_like
            Eccentricity value(s).

        Returns
        -------
        float or ndarray
            Probability density at eccentricity e.
        """
        e = np.asarray(e)
        valid = (e >= self.e_min) & (e <= self.e_max)
        pdf_values = np.zeros_like(e, dtype=float)
        pdf_values[valid] = self.norm
        return pdf_values

    def pdf_per_dex(self, e):
        """Probability density per unit log10(e) of the uniform eccentricity distribution.

        Normalised in log10(e) space to 1.

        Parameters
        ----------
        e : float or array_like
            Eccentricity value(s).

        Returns
        -------
        float or ndarray
            Probability density per unit log10(e) at e.
        """
        return _per_dex_from_pdf(e, self.pdf(e))

    def rvs(self, size=1, rng=None):
        """Draw random samples from the uniform eccentricity distribution.

        Parameters
        ----------
        size : int, optional
            Number of samples to draw (default: 1).
        rng : numpy.random.Generator, optional
            Random number generator. If None, uses np.random.default_rng().

        Returns
        -------
        ndarray
            Random eccentricity samples in [e_min, e_max].
        """
        if rng is None:
            rng = np.random.default_rng()

        return rng.uniform(self.e_min, self.e_max, size=size)


class ZeroEccentricity:
    """Zero eccentricity distribution for circular binary orbits.

    All samples are exactly zero (circular orbits). The PDF is a Dirac delta
    at e=0.

    Examples
    --------
    >>> dist = ZeroEccentricity()
    >>> e_samples = dist.rvs(size=1000)  # All zeros
    """

    def __init__(self):
        """Initialize the zero eccentricity distribution."""
        pass

    def __repr__(self):
        """Return string representation of the distribution."""
        return "ZeroEccentricity()"

    def _repr_html_(self):
        """Return HTML representation for Jupyter notebooks."""
        return "<h3>Zero Eccentricity Distribution</h3><p>e = 0 (circular orbits)</p>"

    def pdf(self, e):
        """Probability density function of the zero eccentricity distribution.

        This is formally a Dirac delta at e=0. Returns 1.0 at e=0, 0.0 elsewhere.

        Parameters
        ----------
        e : float or array_like
            Eccentricity value(s).

        Returns
        -------
        float or ndarray
            1.0 where e==0, 0.0 elsewhere.
        """
        e = np.asarray(e)
        return np.where(e == 0, 1.0, 0.0)

    def pdf_per_dex(self, e):
        """Probability density per unit log10(e) of the zero eccentricity distribution.

        Parameters
        ----------
        e : float or array_like
            Eccentricity value(s).

        Returns
        -------
        float or ndarray
            Probability density per unit log10(e) at e.
            Always 0, since the point mass at e=0 has no density
            in log10(e) space.
        """
        return _per_dex_from_pdf(e, self.pdf(e))

    def rvs(self, size=1, rng=None):
        """Draw random samples from the zero eccentricity distribution.

        Parameters
        ----------
        size : int, optional
            Number of samples to draw (default: 1).
        rng : numpy.random.Generator, optional
            Random number generator (unused, included for API consistency).

        Returns
        -------
        ndarray
            Array of zeros with shape (size,).
        """
        return np.zeros(size)


class LogNormalSeparation:
    """Log-normal orbital separation distribution for binary star systems.

    Orbital separations are drawn from a log-normal distribution in log10 space,
    truncated between specified minimum and maximum values.

    Parameters
    ----------
    mean : float, optional
        Mean of the log10 distribution (default: 0.85, corresponding to ~7.08 Rsun).
    sigma : float, optional
        Standard deviation of the log10 distribution (default: 0.37).
    min : float, optional
        Minimum orbital separation in solar radii (default: 5.0).
    max : float, optional
        Maximum orbital separation in solar radii (default: 1e5).

    Raises
    ------
    ValueError
        If min is not positive, max <= min, or sigma <= 0.

    Examples
    --------
    >>> dist = LogNormalSeparation(mean=3.0, sigma=1.5, min=10, max=1e6)
    >>> separations = dist.rvs(size=1000)

    Notes
    -----
    Uses scipy.stats.truncnorm for efficient sampling from the truncated
    normal distribution in log10 space.
    """

    def __init__(self, mean=0.85, sigma=0.37, min=5.0, max=1e5):
        """Initialize the log-normal separation distribution.

        Parameters
        ----------
        mean : float, optional
            Mean of the log10 distribution (default: 0.85).
        sigma : float, optional
            Standard deviation of the log10 distribution (default: 0.37).
        min : float, optional
            Minimum orbital separation in solar radii (default: 5.0).
        max : float, optional
            Maximum orbital separation in solar radii (default: 1e5).

        Raises
        ------
        ValueError
            If min is not positive, max <= min, or sigma <= 0.
        """
        if min <= 0:
            raise ValueError("min must be positive")
        if max <= min:
            raise ValueError("max must be greater than min")
        if sigma <= 0:
            raise ValueError("sigma must be positive")

        self.mean = mean
        self.sigma = sigma
        self.min = min
        self.max = max

        # Compute truncation bounds for scipy.stats.truncnorm
        self.a_low = (np.log10(min) - mean) / sigma
        self.a_high = (np.log10(max) - mean) / sigma

    def __repr__(self):
        """Return string representation of the distribution."""
        return (f"LogNormalSeparation(mean={self.mean}, sigma={self.sigma}, "
                f"min={self.min}, max={self.max})")

    def _repr_html_(self):
        """Return HTML representation for Jupyter notebooks."""
        return (f"<h3>Log-Normal Separation Distribution</h3>"
                f"<p>mean (log10) = {self.mean}</p>"
                f"<p>sigma (log10) = {self.sigma}</p>"
                f"<p>min = {self.min} R☉</p>"
                f"<p>max = {self.max} R☉</p>")

    def pdf(self, a):
        """Probability density function of the log-normal separation distribution.

        Normalised in a space to 1. See `pdf_per_dex` for the density
        per unit log10(a).

        Parameters
        ----------
        a : float or array_like
            Orbital separation(s) in solar radii.

        Returns
        -------
        float or ndarray
            Probability density per unit a at separation a.
        """
        # pdf(a) = pdf_per_dex(a) * |d(log a)/da| = pdf_per_dex(a) / (a * ln(10))
        return _pdf_from_per_dex(a, self.pdf_per_dex(a))

    def pdf_per_dex(self, a):
        """Probability density per unit log10(a) of the log-normal separation
        distribution.

        Normalised in log10(a) space to 1.

        Parameters
        ----------
        a : float or array_like
            Orbital separation(s) in solar radii.

        Returns
        -------
        float or ndarray
            Probability density per unit log10(a) at separation a.
        """
        a = np.asarray(a)
        valid = (a > 0) & (a >= self.min) & (a <= self.max)
        pdf_values = np.zeros_like(a, dtype=float)

        # truncated normal in log10 space
        pdf_values[valid] = truncnorm.pdf(
            np.log10(a[valid]),
            self.a_low, self.a_high,
            loc=self.mean,
            scale=self.sigma
        )

        return pdf_values

    def rvs(self, size=1, rng=None):
        """Draw random samples from the log-normal separation distribution.

        Parameters
        ----------
        size : int, optional
            Number of samples to draw (default: 1).
        rng : numpy.random.Generator, optional
            Random number generator. If None, uses np.random.default_rng().

        Returns
        -------
        ndarray
            Random orbital separation samples in solar radii.
        """
        if rng is None:
            rng = np.random.default_rng()

        # Sample from truncated normal in log10 space
        log_separations = truncnorm.rvs(
            self.a_low, self.a_high,
            loc=self.mean,
            scale=self.sigma,
            size=size,
            random_state=rng
        )

        # Convert back to linear space
        return 10**log_separations
