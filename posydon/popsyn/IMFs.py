import numpy as np
from scipy.integrate import quad
from abc import ABC, abstractmethod

class IMFBase(ABC):
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
    '''Salpeter E. E., 1955, ApJ, 121, 161'''
    
    def __init__(self, alpha=2.35, m_min=0.01, m_max=200.0):
        self.alpha = alpha
        super().__init__(m_min, m_max)
        
    def __repr__(self):
        return f"Salpeter(alpha={self.alpha}, m_min={self.m_min}, m_max={self.m_max})"
    
    def _repr_html_(self):
        return f"<h3>Salpeter IMF</h3><p>alpha = {self.alpha}</p><p>m_min = {self.m_min}</p><p>m_max = {self.m_max}</p>"
    
    def imf(self, m):
        if np.any(m <= 0):
            raise ValueError("Mass must be positive.")
        return m ** (-self.alpha)

class Kroupa2001(IMFBase):
    '''Kroupa P., 2001, MNRAS, 322, 231'''
    
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
        return np.select(conditions, [f(m) for f in choices], default=0.0)

class Chabrier2003(IMFBase):
    '''Chabrier, 2003 IMF
    https://ui.adsabs.harvard.edu/abs/2003PASP..115..763C/abstract'''
    
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
        m = np.asarray(m)
        lognormal = (1.0 / (m * np.sqrt(2 * np.pi) * self.sigma) *
                     np.exp(- (np.log10(m) - np.log10(self.m_c))**2 / (2 * self.sigma**2)))
        C = (1.0 / (self.m_break * np.sqrt(2 * np.pi) * self.sigma) *
             np.exp(- (np.log10(self.m_break) - np.log10(self.m_c))**2 / (2 * self.sigma**2)))
        powerlaw = C * (m / self.m_break) ** (-self.alpha)
        return np.where(m < self.m_break, lognormal, powerlaw)