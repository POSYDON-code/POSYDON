from scipy.integrate import quad
import numpy as np

class flat_mass_ratio():
    
    def __init__(self, q_min=0.05, q_max=1):
        self.q_min = q_min
        self.q_max = q_max
        self.norm = self._calculate_normalization()
        
    def __repr__(self):
        return f"FlatMassRatio(q_min={self.q_min}, q_max={self.q_max})"
    
    def _repr_html_(self):
        return f"<h3>Flat Mass Ratio Distribution</h3><p>q_min = {self.q_min}</p><p>q_max = {self.q_max}</p>"
    
    def _calculate_normalization(self):
        """Calculate the normalization constant for the flat mass ratio distribution.

        Returns
        -------
        float
            The normalization constant ensuring the PDF integrates to 1.
        """
        integral, _ = quad(self.flat_mass_ratio, self.q_min, self.q_max)
        if integral == 0:
            raise ValueError("Normalization integral is zero. Check mass ratio parameters.")
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

        """
        if np.any(q < 0) or np.any(q > 1):
            raise ValueError("Mass ratio must be between 0 and 1. Check mass ratio parameters.")
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
            Probability density at mass ratio m.
        """
        q = np.asarray(q)
        valid = (q >= self.q_min) & (q <= self.q_max)
        pdf_values = np.zeros_like(q, dtype=float)
        pdf_values[valid] = self.flat_mass_ratio(q[valid]) * self.norm
        return pdf_values
    

class Sana12Period():
    '''Sana H., et al., 2012, Science, 337, 444'''
    
    def __init__(self, p_min=0.35, p_max=10**3.5):
        self.p_min = p_min
        self.p_max = p_max
        self.mbreak = 15
        
        # mass boundary is at 15 Msun
        self.low_mass_norm = self._calculate_normalization(self.mbreak-1)
        self.high_mass_norm = self._calculate_normalization(self.mbreak+1)
        self.norm = lambda m1: np.where(m1 <= self.mbreak, self.low_mass_norm, self.high_mass_norm)
        
        
    def __repr__(self):
        return f"Sana12Period(p_min={self.p_min}, p_max={self.p_max})"
    
    def _repr_html_(self):
        return f"<h3>Sana12 Period Distribution</h3><p>p_min = {self.p_min}</p><p>p_max = {self.p_max}</p>"
    
    def _calculate_normalization(self, m1):
        """Calculate the normalization constant for the Sana12 period distribution.

        Returns
        -------
        float
            The normalization constant ensuring the PDF integrates to 1.
        """
        integral =  quad(self.sana12_period, np.log10(self.p_min), np.log10(self.p_max), args=(m1))[0]
        if integral == 0:
            raise ValueError("Normalization integral is zero. Check period parameters.")
        return 1.0 / integral
    
    def sana12_period(self, logp, m1):
        """Compute the Sana12 period distribution at given periods.

        Parameters
        ----------
        logp : float or array_like
            Period(s).
        m1 : float or array_like
            Mass value(s). Can be a single value or an array of the same length as p.

        Returns
        -------
        float or ndarray
            Distribution value(s).
        """
        if np.any(m1 <= 0):
            raise ValueError("Period and mass must be positive. Check parameters.")

        logp = np.asarray(logp)
        m1 = np.asarray(m1)
        if m1.size == 1:
            m1 = np.full_like(logp, m1)
        elif m1.shape != logp.shape:
            raise ValueError("m1 must be a single value or an array of the same length as p.")
        
        result = np.zeros_like(logp, dtype=float)        

        mask1 = (m1 > 0) & (m1 <= self.mbreak)
        result[mask1] = 1

        mask2 = m1 > self.mbreak
        pi = 0.55
        beta = 1 - pi
        A = np.log10(10**0.15)**(-pi) * (np.log10(10**0.15) - np.log10(self.p_min))
        B = 1./beta * (np.log10(self.p_max)**beta - np.log10(10**0.15)**beta)
        C = 1./(A + B)
        logp_valid = logp[mask2]
        result[mask2] = np.where(
            (np.log10(self.p_min) <= logp_valid) & (logp_valid < 0.15),
            C * 0.15**(-pi),
            np.where(
                (0.15 <= logp_valid) & (logp_valid <= np.log10(self.p_max)),
                C * logp_valid**(-pi),
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
            Mass value(s). Can be a single value or an array of the same length as p.
        
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
            raise ValueError("m1 must be a single value or an array of the same length as p.")
        
        p = np.asarray(p)
        valid = (p >= self.p_min) & (p <= self.p_max)
        pdf_values = np.zeros_like(p, dtype=float)
        
        logp = np.log10(p)
        pdf_values[valid] = self.sana12_period(logp[valid], m1[valid]) * self.norm(m1[valid])
        return pdf_values

        