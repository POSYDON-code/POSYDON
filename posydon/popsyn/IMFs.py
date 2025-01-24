import numpy as np
from scipy.integrate import quad

class Salpeter:
    '''Salpeter E. E., 1955, ApJ, 121, 161'''
    
    def __init__(self, alpha=2.35, m_min=0.01, m_max=200.0):
        """Initialize the Salpeter Initial Mass Function (IMF).

        Parameters
        ----------
        alpha : float, optional
            The power-law index of the Salpeter IMF (default is 2.35).
        m_min : float, optional
            The minimum stellar mass (default is 0.01).
        m_max : float, optional
            The maximum stellar mass (default is 200.0).
        """
        self.alpha = alpha
        self.m_min = m_min
        self.m_max = m_max
        self.norm = self._calculate_normalization()
        
    def __repr__(self):
        return f"Salpeter(alpha={self.alpha}, m_min={self.m_min}, m_max={self.m_max})"
    
    def _repr_html_(self):
        return f"<h3>Salpeter IMF</h3><p>alpha = {self.alpha}</p><p>m_min = {self.m_min}</p><p>m_max = {self.m_max}</p>"
    
    def _calculate_normalization(self):
        """Calculate the normalization constant for the IMF.

        Returns
        -------
        float
            The normalization constant ensuring the PDF integrates to 1.
        """
        integral, _ = quad(self.salpeter_IMF, self.m_min, self.m_max)
        if integral == 0:
            raise ValueError("Normalization integral is zero. Check IMF parameters.")
        return 1.0 / integral

    def salpeter_IMF(self, m):
        """Compute the Salpeter IMF at a given mass.

        Parameters
        ----------
        m : float
            Stellar mass.

        Returns
        -------
        float
            Number of stars per unit mass.
        """
        if np.any(m <= 0):
            raise ValueError("Mass must be positive.")
        return m ** (-self.alpha)

    def pdf(self, m):
        """Probability density function of the Salpeter IMF.

        Parameters
        ----------
        m : float or array_like
            Stellar mass(es).

        Returns
        -------
        float or ndarray
            Probability density at mass m.
        """
        m = np.asarray(m)
        valid = (m >= self.m_min) & (m <= self.m_max)
        pdf_values = np.zeros_like(m, dtype=float)
        pdf_values[valid] = self.salpeter_IMF(m[valid]) * self.norm
        return pdf_values
    
    
class Kroupa2001:
    '''Kroupa P., 2001, MNRAS, 322, 231'''
    
    def __init__(self,
                 alpha1=0.3, alpha2=1.3, alpha3=2.3,
                 m1break=0.08, m2break=0.5,
                 m_min=0.01, m_max=200.0,):
        """Initialize the Kroupa (2001) Initial Mass Function (IMF).
        
        Parameters
        ----------
        alpha1 : float, optional
            The power-law index for m < m1break (default is 0.3).
        alpha2 : float, optional
            The power-law index for m1break <= m < m2break (default is 1.3).
        alpha3 : float, optional
            The power-law index for m >= m2break (default is 2.3).        
        m1break : float, optional
            The mass at the first break point (default is 0.08).
        m2break : float, optional
            The mass at the second break point (default is 0.5).
        m_min : float, optional
            The minimum stellar mass (default is 0.01).
        m_max : float, optional
            The maximum stellar mass (default is 200.0).
        """
        self.alpha1 = alpha1
        self.alpha2 = alpha2
        self.alpha3 = alpha3
        self.m_min = m_min
        self.m_max = m_max
        self.m1break = m1break
        self.m2break = m2break
        self.norm = self._calculate_normalization()
        
    def __repr__(self):
        return (f"Kroupa2001(alpha1={self.alpha1}, alpha2={self.alpha2}, alpha3={self.alpha3}, "
                f"m1break={self.m1break}, m2break={self.m2break}, m_min={self.m_min}, m_max={self.m_max})")
        
    def _repr_html_(self):
        return (f"<h3>Kroupa (2001) IMF</h3><p>alpha1 = {self.alpha1}</p><p>alpha2 = {self.alpha2}</p>"
                f"<p>alpha3 = {self.alpha3}</p><p>m1break = {self.m1break}</p><p>m2break = {self.m2break}</p>"
                f"<p>m_min = {self.m_min}</p><p>m_max = {self.m_max}</p>")
        
    def _calculate_normalization(self):
        """Calculate the normalization constant for the IMF.

        Returns
        -------
        float
            The normalization constant ensuring the PDF integrates to 1.
        """
        integral, _ = quad(self.kroupa_IMF, self.m_min, self.m_max)
        if integral == 0:
            raise ValueError("Normalization integral is zero. Check IMF parameters.")
        return 1.0 / integral
    
    def kroupa_IMF(self, m):
        '''Compute the Kroupa IMF at a given mass.
        
        Parameters
        ----------
        m : float
            Stellar mass.
        
        Returns
        -------
        float
            Number of stars per unit mass.
        '''
        if np.any(m <= 0):
            raise ValueError("Mass must be positive.")
        
        conditions = [
            m < self.m1break,
            (m >= self.m1break) & (m < self.m2break),
            m >= self.m2break
        ]
        choices = [
            lambda m: (m/self.m_min) ** (-self.alpha1),
            lambda m: (m/self.m1break) ** (-self.alpha2) * (self.m_min/self.m1break) ** self.alpha1, 
            lambda m: (m/self.m2break) ** (-self.alpha3) * (self.m_min/self.m1break) ** self.alpha1 * (self.m1break/self.m2break) ** self.alpha2
        ]
        imf_values = np.select(conditions, [f(m) for f in choices], default=0.0)
        
        return imf_values
        
    def pdf(self, m):
        """Probability density function of the Kroupa IMF.

        Parameters
        ----------
        m : float or array_like
            Stellar mass(es).

        Returns
        -------
        float or ndarray
            Probability density at mass m.
        """
        m = np.asarray(m)
        valid = (m >= self.m_min) & (m <= self.m_max)
        pdf_values = np.zeros_like(m, dtype=float)
        pdf_values[valid] = self.kroupa_IMF(m[valid]) * self.norm
        return pdf_values
