import numpy as np
from scipy.integrate import quad
from abc import ABC, abstractmethod

class IMFBase(ABC):
    '''Base class for Initial Mass Functions (IMFs)'''
    def __init__(self, m_min, m_max):
        self.m_min = m_min
        self.m_max = m_max
        self.norm = self._calculate_normalization()

    def _calculate_normalization(self):
        integral, _ = quad(self.imf, self.m_min, self.m_max)
        if not (integral > 0):
            raise ValueError("Normalization integral is zero. Check IMF parameters.")
        return 1.0 / integral

    def pdf(self, m):
        """
        Compute the probability density function (PDF) values for the input mass(es).

        Parameters:
            m (scalar or array-like): The mass or masses at which to evaluate the PDF [Msun].

        Returns:
            numpy.ndarray: An array of PDF values corresponding to each entry in m, with invalid
            masses (outside [self.m_min, self.m_max]) assigned a value of zero.
        """
        m = np.asarray(m)
        valid = (m >= self.m_min) & (m <= self.m_max)
        pdf_values = np.zeros_like(m, dtype=float)
        pdf_values[valid] = self.imf(m[valid]) * self.norm
        return pdf_values

    # # rejection sampling
    # def sample(self, n, num_grid=10000, max_iter=100000):
    #     # Estimate maximum pdf value over a grid
    #     x_grid = np.linspace(self.m_min, self.m_max, num_grid)
    #     f_max = np.max(self.pdf(x_grid))
    #     samples = []
    #     iterations = 0
    #     while len(samples) < n and iterations < max_iter:
    #         x_cand = np.random.uniform(self.m_min, self.m_max)
    #         f_val = self.pdf(np.array([x_cand]))[0]
    #         if np.random.uniform(0, f_max) < f_val:
    #             samples.append(x_cand)
    #         iterations += 1
    #     if len(samples) < n:
    #         raise RuntimeError("Rejection sampling failed to converge")
    #     return np.array(samples)

    # def sample_inverse(self, n, num_grid=10000):
    #     x_grid = np.linspace(self.m_min, self.m_max, num_grid)
    #     pdf_vals = self.pdf(x_grid)
    #     cdf_vals = np.cumsum(pdf_vals)
    #     cdf_vals /= cdf_vals[-1]  # normalize CDF to 1
    #     u = np.random.uniform(0, 1, n)  # uniform samples in [0,1]
    #     return np.interp(u, cdf_vals, x_grid)

    # This forces that the method is implemented in a sub-class
    @abstractmethod
    def imf(self, m):
        pass

class Salpeter(IMFBase):
    """
    Initial Mass Function based on Salpeter (1955), which is defined as:
    
        dN/dM = m^-2.35
    
    References
    ----------
    Salpeter, 1955, ApJ, 121, 161
    https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S/abstract
    
    Parameters
    ----------
    alpha : float, optional
        The power-law index of the IMF (default is 2.35).
    m_min : float, optional
        The minimum allowable mass (default is 0.01) [Msun].
    m_max : float, optional
        The maximum allowable mass (default is 200.0) [Msun].
        
    Attributes
    ----------
    alpha : float
        Power-law index used in the IMF calculation.
    m_min : float
        Minimum stellar mass for the IMF [Msun].
    m_max : float
        Maximum stellar mass for the IMF [Msun].
    """
    
    def __init__(self, alpha=2.35, m_min=0.01, m_max=200.0):
        self.alpha = alpha
        super().__init__(m_min, m_max)
        
    def __repr__(self):
        return f"Salpeter(alpha={self.alpha}, m_min={self.m_min}, m_max={self.m_max})"
    
    def _repr_html_(self):
        return f"<h3>Salpeter IMF</h3><p>alpha = {self.alpha}</p><p>m_min = {self.m_min}</p><p>m_max = {self.m_max}</p>"
    
    def imf(self, m):
        '''Computes the IMF value for a given mass or array of masses 'm'. Raises a
        ValueError if any value in 'm' is less than or equal to zero.
        
        Parameters
        ----------
        m : float or array_like
            Stellar mass or array of stellar masses [Msun].
            
        Returns
        -------
        float or ndarray
            The IMF value for the given mass or masses.
        '''
        if np.any(m <= 0):
            raise ValueError("Mass must be positive.")
        return m ** (-self.alpha)

class Kroupa2001(IMFBase):
    """Initial Mass Function based on Kroupa (2001), which is 
    defined as a broken power-law:
    
        dN/dM = m^-0.3 for 0.01 <= m < 0.08
        dN/dM = m^-1.3 for 0.08 <= m < 0.5
        dN/dM = m^-2.3 for 0.5 <= m <= 200.0
        
    References
    ----------
    Kroupa, 2001, MNRAS, 322, 231
    https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K/abstract

    Parameters
    ----------
    alpha1 : float, optional
        The power-law index for m < m1break (default is 0.3).
    alpha2 : float, optional
        The power-law index for m1break <= m < m2break (default is 1.3).
    alpha3 : float, optional
        The power-law index for m >= m2break (default is 2.3).
    m1break : float, optional
        The first mass break (default is 0.08) [Msun].
    m2break : float, optional
        The second mass break (default is 0.5) [Msun].
    m_min : float, optional
        The minimum allowable mass (default is 0.01) [Msun].
    m_max : float, optional
        The maximum allowable mass (default is 200.0) [Msun].
        
    Attributes
    ----------
    alpha1 : float
        Power-law index for m < m1break.
    alpha2 : float
        Power-law index for m1break <= m < m2break.
    alpha3 : float
        Power-law index for m >= m2break.
    m1break : float
        First mass break [Msun].
    m2break : float
        Second mass break [Msun].
    m_min : float
        Minimum stellar mass for the IMF [Msun].
    m_max : float
        Maximum stellar mass for the IMF [Msun].        
   
    """
    
    def __init__(self,
                 alpha1=0.3, alpha2=1.3, alpha3=2.3,
                 m1break=0.08, m2break=0.5,
                 m_min=0.01, m_max=200.0):
        self.alpha1 = alpha1
        self.alpha2 = alpha2
        self.alpha3 = alpha3
        self.m1break = m1break
        self.m2break = m2break
        super().__init__(m_min, m_max)
        
    def __repr__(self):
        return (f"Kroupa2001(alpha1={self.alpha1}, alpha2={self.alpha2}, alpha3={self.alpha3}, "
                f"m1break={self.m1break}, m2break={self.m2break}, m_min={self.m_min}, m_max={self.m_max})")
        
    def _repr_html_(self):
        return (f"<h3>Kroupa (2001) IMF</h3><p>alpha1 = {self.alpha1}</p><p>alpha2 = {self.alpha2}</p>"
                f"<p>alpha3 = {self.alpha3}</p><p>m1break = {self.m1break}</p><p>m2break = {self.m2break}</p>"
                f"<p>m_min = {self.m_min}</p><p>m_max = {self.m_max}</p>")
        
    def imf(self, m):
        '''Computes the IMF value for a given mass or array of masses 'm'. Raises a
        ValueError if any value in 'm' is less than or equal to zero.
        
        Parameters
        ----------
        m : float or array_like
            Stellar mass or array of stellar masses [Msun].
            
        Returns
        -------
        float or ndarray
            The IMF value for the given mass or masses.
        '''
        m = np.asarray(m)
        if np.any(m <= 0):
            raise ValueError("Mass must be positive.")
        
        mask1 = m < self.m1break
        mask2 = (m >= self.m1break) & (m < self.m2break)
        mask3 = m >= self.m2break

        # Precompute constants to ensure continuity
        const1 = (self.m1break / self.m_min) ** (-self.alpha1)
        const2 = const1 * (self.m2break / self.m1break) ** (-self.alpha2)

        out = np.empty_like(m, dtype=float)
        out[mask1] = (m[mask1] / self.m_min) ** (-self.alpha1)
        out[mask2] = const1 * (m[mask2] / self.m1break) ** (-self.alpha2)
        out[mask3] = const2 * (m[mask3] / self.m2break) ** (-self.alpha3)
        return out

class Chabrier2003(IMFBase):
    """
    Chabrier2003 Initial Mass Function (IMF), which is defined as a lognormal distribution
    for low-mass stars and a power-law distribution for high-mass stars:
    
        dN/dM = 1/(m * sqrt(2*pi) * sigma) * exp(- (log10(m) - log10(m_c))^2 / (2 * sigma^2)) for m < m_break
        dN/dM = C * (m / m_break)^(-alpha) for m >= m_break
        
    Reference
    ----------
    Chabrier, (2003). PASP, 115(809), 763-795.
    https://ui.adsabs.harvard.edu/abs/2003PASP..115..763C/abstract
        
    Parameters
    ----------
    m_c : float, optional
        The characteristic mass of the lognormal component (default is 0.22) [Msun].
    sigma : float, optional
        The dispersion (standard deviation) of the lognormal component (default is 0.57).
    alpha : float, optional
        The power-law index governing the high-mass end of the IMF (default is 2.3).
    m_break : float, optional
        The mass at which the IMF transitions from a lognormal to a power-law behavior (default is 1.0) [Msun].
    m_min : float, optional
        The minimum mass considered in the IMF (default is 0.01) [Msun].
    m_max : float, optional
        The maximum mass considered in the IMF (default is 200.0) [Msun].
    """
    
    def __init__(self, m_c=0.22, sigma=0.57, alpha=2.3, m_break=1.0, m_min=0.01, m_max=200.0):
        self.m_c = m_c
        self.sigma = sigma
        self.alpha = alpha
        self.m_break = m_break
        super().__init__(m_min, m_max)
    
    def __repr__(self):
        return (f"Chabrier(m_c={self.m_c}, sigma={self.sigma}, alpha={self.alpha}, "
                f"m_break={self.m_break}, m_min={self.m_min}, m_max={self.m_max})")

    def _repr_html_(self):
        return (f"<h3>Chabrier IMF</h3><p>m_c = {self.m_c}</p><p>sigma = {self.sigma}</p>"
                f"<p>alpha = {self.alpha}</p><p>m_break = {self.m_break}</p>"
                f"<p>m_min = {self.m_min}</p><p>m_max = {self.m_max}</p>")

    def imf(self, m):
        '''Computes the IMF value for a given mass or array of masses 'm'. Raises a
        ValueError if any value in 'm' is less than or equal to zero.
        
        Parameters
        ----------
        m : float or array_like
            Stellar mass or array of stellar masses [Msun].
            
        Returns
        -------
        float or ndarray
            The IMF value for the given mass or masses
        '''
        m = np.asarray(m)
        lognormal = (1.0 / (m * np.sqrt(2 * np.pi) * self.sigma) *
                     np.exp(- (np.log10(m) - np.log10(self.m_c))**2 / (2 * self.sigma**2)))
        C = (1.0 / (self.m_break * np.sqrt(2 * np.pi) * self.sigma) *
             np.exp(- (np.log10(self.m_break) - np.log10(self.m_c))**2 / (2 * self.sigma**2)))
        powerlaw = C * (m / self.m_break) ** (-self.alpha)
        return np.where(m < self.m_break, lognormal, powerlaw)