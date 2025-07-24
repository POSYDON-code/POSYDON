from scipy.integrate import quad
import numpy as np

class FlatMassRatio:
    """Flat mass ratio distribution for binary star systems.
    
    A uniform distribution for mass ratios q = m2/m1 within specified bounds.
    This distribution assigns equal probability to all mass ratios within the
    given range [q_min, q_max].
    
    Parameters
    ----------
    q_min : float, optional
        Minimum mass ratio (default: 0.05). Must be in (0, 1].
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
            Minimum mass ratio (default: 0.05). Must be in (0, 1].
        q_max : float, optional
            Maximum mass ratio (default: 1.0). Must be in (0, 1].
            
        Raises
        ------
        ValueError
            If q_min or q_max are not in (0, 1], or if q_min >= q_max.
        """
        if not (0 < q_min <= 1):
            raise ValueError("q_min must be in (0, 1)")
        if not (0 < q_max <= 1):
            raise ValueError("q_max must be in (0, 1)")
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
        valid = (q >= self.q_min) & (q <= self.q_max)
        pdf_values = np.zeros_like(q, dtype=float)
        pdf_values[valid] = self.flat_mass_ratio(q[valid]) * self.norm
        return pdf_values
    

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
        integral =  quad(self.sana12_period,
                         np.log10(self.p_min),
                         np.log10(self.p_max),
                         args=(m1,))[0]
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

        logp = np.asarray(logp)
        m1 = np.asarray(m1)
        if m1.size == 1:
            m1 = np.full_like(logp, m1)
        elif m1.shape != logp.shape:
            raise ValueError("m1 must be a single value or an array of "
                             "the same length as p.")
        
        result = np.zeros_like(logp, dtype=float)        

        mask1 = (m1 > 0) & (m1 <= self.mbreak)
        result[mask1] = 1

        mask2 = m1 > self.mbreak
        beta = 1 - self.slope
        A = np.log10(10**0.15)**(-self.slope) * (np.log10(10**0.15) - np.log10(self.p_min))
        B = 1./beta * (np.log10(self.p_max)**beta - np.log10(10**0.15)**beta)
        C = 1./(A + B)
        logp_valid = logp[mask2]
        result[mask2] = np.where(
            (np.log10(self.p_min) <= logp_valid) & (logp_valid < 0.15),
            C * 0.15**(-self.slope),
            np.where(
                (0.15 <= logp_valid) & (logp_valid <= np.log10(self.p_max)),
                C * logp_valid**(-self.slope),
                0
            )
        )
       
        return result
            
        
    def pdf(self, p, m1):
        '''Probability density function of the Sana12 period distribution.
        Conditional on m1.
        
        Normalised in logP space to 1.
        
        Parameters
        ----------
        p : float or array_like
            Period(s).
        m1 : float or array_like
            Mass value(s). Single value or an array of the same length as p.
        
        Returns
        -------
        float or ndarray
            Probability density at period p.
        
        '''
        p = np.asarray(p)
        m1 = np.asarray(m1)
        if m1.size == 1:
            m1 = np.full_like(p, m1)
        elif m1.shape != p.shape:
            raise ValueError("m1 must be a single value or an array of "
                             "the same length as p.")
        
        valid = (p >= self.p_min) & (p <= self.p_max)
        pdf_values = np.zeros_like(p, dtype=float)
        
        logp = np.log10(p)
        pdf_values[valid] = (self.sana12_period(logp[valid], m1[valid]) 
                             * self.norm(m1[valid]))
        return pdf_values


class PowerLawPeriod():
    '''Power law period distribution with slope pi and boundaries m_min, m_max.
    
    Normalised in logP space to 1.
    
    Parameters
    ----------
    slope : float
        Slope of the power law.
    p_min : float
        Minimum period.
    p_max : float
        Maximum period.
    
    '''
    
    def __init__(self, slope=-0.55, p_min=1.4, p_max=3e6):
        """Initialize the power law period distribution.
        
        Parameters
        ----------
        slope : float, optional
            Power law slope (default: -0.55).
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
            The normalization constant ensuring the PDF integrates to 1.
        """
        integral = quad(self.power_law_period,
                        np.log10(self.p_min),
                        np.log10(self.p_max))[0]
        if integral == 0: # pragma: no cover
            raise ValueError("Normalization integral is zero. "
                             "Check period parameters.")
        return 1.0 / integral
    
    def power_law_period(self, logp):
        """Compute the power law period distribution at given periods.

        Parameters
        ----------
        logp : float or array
            Log10 of period(s).

        Returns
        -------
        float or ndarray
            Distribution value(s).
        """
        p = 10**np.asarray(logp)

        result = np.zeros_like(p, dtype=float)
        valid = (p >= self.p_min) & (p <= self.p_max)
        result[valid] = p[valid]**self.slope

        return result
        
    def pdf(self, p):
        """Probability density function of the power law period distribution.

        Parameters
        ----------
        p : float or array_like
            Period(s).

        Returns
        -------
        float or ndarray
            Probability density at period p.
        """
        p = np.asarray(p)
        valid = (p > 0) & (p >= self.p_min) & (p <= self.p_max)
        pdf_values = np.zeros_like(p, dtype=float)

        logp = np.log10(p[valid])
        pdf_values[valid] = (self.power_law_period(logp)
                             * self.norm)
        
        return pdf_values
