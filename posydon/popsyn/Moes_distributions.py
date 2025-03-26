"""
Generate initial parameters following Moe & Di Stefano (2017).
This module is an adoption of the translation of Moe's IDL code into python by
Mads Sørensen (V. 0.1; 2017/02/03)

NOTE my Maxwell Moe: This version produces only the statistical distributions
of single stars, binaries, and inner binaries in hierarchical triples. Outer
tertiaries in hierarchical triples are NOT generated. Moreover, given a set of
companions, all with period P to primary mass M1, this version currently uses
an approximation to determine the fraction of those companions that are inner
binaries vs. outer triples. Nevertheless, this approximation reproduces the
overall multiplicity statistics.
"""


__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]


import numpy as np
from scipy.integrate import newton_cotes, quad


class Moe2017PsandQs():
    """Generate initial parameters following Moe & Di Stefano (2017) [1]_.
    
    References
    ----------
    .. [1] Moe, M. and Di Stefano, R., “Mind Your Ps and Qs: The Interrelation
        between Period (P) and Mass-ratio (Q) Distributions of Binary Stars”,
        <i>The Astrophysical Journal Supplement Series</i>, vol. 230, no. 2,
        Art. no. 15, IOP, 2017. doi:10.3847/1538-4365/aa6fb6.
    """

    def _idl_tabulate(self, x, f, p=5):
        """Discrete integration in chunks to avoid instabilities for
        Newton-Cotes rule at high degree.
        This function is similar to the IDL function int_tabulated except, this
        function seems to be slightly more exact in its solution. Therefore,
        relative to the IDL code, there are small numerical differences.

        Parameters
        ----------
        x : ndarray of floats
            Values of integration variable.
        f : ndarray of floats
            Function values at x, therefore need to have the same shape as x.
        p : int (default: 5)
            Chunk size. Too small and too high values may suffer from numerical
            uncertainties.

        Returns
        -------
        ret : float
            Integral value.
        """

        def integrate_newton_cotes(x, f):
            """Integrate f from x[0] to x[-1] via Newton cotes.

            Parameters
            ----------
            x : ndarray of floats
                Values of integration variable.
            f : ndarray of floats
                Function values at x, therefore need to have the same shape as
                x.

            Returns
            -------
            float
                Integral value.
            """
            if x.shape[0] < 2:
                # x[0] and x[-1] are the same, hence integral over no range
                return 0
            # average sample spacing
            dx = (x[-1] - x[0]) / (x.shape[0] - 1)
            # sample positions normalized to [0,x.shape[0]-1]
            rn = (x - x[0]) / dx
            # round edges to avoid numerical errors
            rn[0] = np.round(rn[0], 5)
            rn[-1] = np.round(rn[-1], 5)
            # weights of Newton cotes
            weights = newton_cotes(rn)[0]
            # integration
            return dx * np.dot(weights, f)

        ret = 0
        for idx in range(0, x.shape[0], p - 1):
            # integrate by chunks of size p
            ret += integrate_newton_cotes(x[idx:idx + p], f[idx:idx + p])
        return ret

    def __init__(self, n_M1=101, n_logP=158, n_q=91, n_e=100):
        """Initializing the class.

        Parameters
        ----------
        n_M1 : int (Default: 101)
            Number of supporting points in the primary's mass.
        n_logP : int (Default: 158)
            Number of supporting points in the logarythmic period.
        n_q : int (Default: 91)
            Number of supporting points in the mass ratio.
        n_e : int (Default: 100)
            Number of supporting points in the eccentricity.
        """
        self.numM1 = n_M1
        self.numlogP = n_logP
        self.numq = n_q
        self.nume = n_e
        # ranges where M+D17 has statistics corrected for selection effects:
        # 0.8 < M1/Msun < 40 with self.numM1 steps in log space
        self.M1v = 10**(np.linspace(np.log10(0.8), np.log10(40.0), self.numM1))
        # 0.15 < log10(P/day) < 8.0 with self.numlogP steps
        self.logPv = np.linspace(0.15, 8.0, self.numlogP)
        # 0.10 < q < 1.00 with self.numq steps
        self.qv = np.linspace(0.1, 1.0, self.numq)
        # 0.0001 < e < 0.9901 with self.nume steps
        # Maxwell Moe: set minimum to non-zero value to avoid numerical errors
        self.ev = 0.0001 + np.linspace(0.0, 0.99, self.nume)
        # Distribution functions (first setup gird, fill later in loop):
        # Frequency of companions with q > 0.1 per decade of orbital period.
        # Bottom panel in Fig. 36 of M+D17
        self.flogP_sq = np.zeros([self.numlogP, self.numM1])
        # Given M1 and P, the cumulative distribution of mass ratios q
        self.cumqdist = np.zeros([self.numq, self.numlogP, self.numM1])
        # Given M1 and P, the cumulative distribution of eccentricities e
        self.cumedist = np.zeros([self.nume, self.numlogP, self.numM1])
        # Given M1 and P, the probability that the companion is a member of the
        # inner binary (currently an approximation). 100% for log P < 1.5,
        # decreases with increasing P.
        self.probbin = np.zeros([self.numlogP, self.numM1])
        # Given M1, the cumulative period distribution of the inner binary
        # Normalized so that max(cumPbindist) = total binary frac. (NOT unity)
        self.cumPbindist = np.zeros([self.numlogP, self.numM1])
        # Slope alpha of period distribution across intermediate periods
        # 2.7 - DlogP < log P < 2.7 + DlogP, see Section 9.3 and Eqn. 23.
        alpha = 0.018
        DlogP = 0.7
        # Heaviside function for twins with 0.95 < q < 1.00
        Heaviside = np.zeros_like(self.qv)
        Heaviside[self.qv >= 0.95] = 1.0
        # normalize so that integral is unity
        Heaviside = Heaviside / self._idl_tabulate(self.qv, Heaviside)
        # Relevant indices with respect to mass ratio
        indlq = np.flatnonzero(self.qv >= 0.3)
        indsq = (self.qv < 0.3)
        indq0p3 = np.min(indlq)

        # Loop through primary mass
        for i in range(0, self.numM1):
            myM1 = self.M1v[i]
            # Twin fraction parameters that are dependent on M1 only;
            # section 9.1
            FtwinlogPle1 = 0.3 - 0.15 * np.log10(myM1) # Eqn. 6
            logPtwin = 8.0 - myM1                      # Eqn. 7a
            if (myM1 >= 6.5):
                logPtwin = 1.5                         # Eqn. 7b
            # Frequency of companions with q > 0.3 at different orbital periods
            # and dependent on M1 only; section 9.3
            flogPle1 = (0.020 + 0.04 * np.log10(myM1)
                        + 0.07 * (np.log10(myM1))**2.0)   #; Eqn. 20
            flogPeq2p7 = (0.039 + 0.07 * np.log10(myM1)
                          + 0.01 * (np.log10(myM1))**2.0) #; Eqn. 21
            flogPeq5p5 = (0.078 - 0.05 * np.log10(myM1)
                          + 0.04 * (np.log10(myM1))**2.0) #; Eqn. 22

            # Loop through orbital period P
            for j in range(0, self.numlogP):
                mylogP = self.logPv[j]
                # Given M1 and P, set excess twin fraction;
                # section 9.1 and Eqn. 5
                if(mylogP <= 1.0):
                    Ftwin = FtwinlogPle1
                else:
                    Ftwin = FtwinlogPle1 * (1.0-(mylogP-1)/(logPtwin-1.0))
                if(mylogP >= logPtwin):
                    Ftwin = 0.0
                # Power-law slope gamma_largeq for M1 < 1.2 Msun and various P;
                # Eqn. 9
                if(mylogP <= 5.0):
                    gl_1p2 = -0.5
                if(mylogP > 5.0):
                    gl_1p2 = -0.5 - 0.3 * (mylogP-5.0)
                # Power-law slope gamma_largeq for M1 = 3.5 Msun and various P;
                # Eqn. 10
                if(mylogP <= 1.0):
                    gl_3p5 = -0.5
                if((mylogP > 1.0) and (mylogP <= 4.5)):
                    gl_3p5 = -0.5 - 0.2 * (mylogP-1.0)
                if((mylogP > 4.5) and (mylogP <= 6.5)):
                    gl_3p5 = -1.2 - 0.4 * (mylogP-4.5)
                if(mylogP > 6.5):
                    gl_3p5 = -2.0
                # Power-law slope gamma_largeq for M1 > 6 Msun and various P;
                # Eqn. 11
                if(mylogP <= 1.0):
                    gl_6 = -0.5
                if((mylogP > 1.0) and (mylogP <= 2.0)):
                    gl_6 = -0.5 - 0.9 * (mylogP-1.0)
                if((mylogP > 2.0) and (mylogP <= 4.0)):
                    gl_6 = -1.4 - 0.3 * (mylogP-2.0)
                if(mylogP > 4.0):
                    gl_6 = -2.0
                # Given P, interpolate gamma_largeq w/ respect to M1 at myM1
                if(myM1 <= 1.2):
                    gl = gl_1p2
                if((myM1 > 1.2) and (myM1 <= 3.5)):
                    gl = np.interp(np.log10(myM1), np.log10([1.2,3.5]),
                                   [gl_1p2,gl_3p5])
                if((myM1 > 3.5) and (myM1 <= 6.0)):
                    gl = np.interp(np.log10(myM1), np.log10([3.5,6.0]),
                                   [gl_3p5,gl_6])
                if(myM1 > 6.0):
                    gl = gl_6
                # Power-law slope gamma_smallq for M1 < 1.2 Msun and all P;
                # Eqn. 13
                gs_1p2 = 0.3
                # Power-law slope gamma_smallq for M1 = 3.5 Msun and various P;
                # Eqn. 14
                if(mylogP <= 2.5):
                    gs_3p5 = 0.2
                if((mylogP > 2.5) and (mylogP <= 5.5)):
                    gs_3p5 = 0.2 - 0.3 * (mylogP-2.5)
                if(mylogP > 5.5):
                    gs_3p5 = -0.7 - 0.2 * (mylogP-5.5)
                # Power-law slope gamma_smallq for M1 > 6 Msun and various P;
                # Eqn. 15
                if(mylogP <= 1.0):
                    gs_6 = 0.1
                if((mylogP > 1.0) and (mylogP <= 3.0)):
                    gs_6 = 0.1 - 0.15 * (mylogP-1.)
                if((mylogP > 3.0) and (mylogP <= 5.6)):
                    gs_6 = -0.2 - 0.50 * (mylogP-3.)
                if(mylogP > 5.6):
                    gs_6 = -1.5
                # Given P, interpolate gamma_smallq w/ respect to M1 at myM1
                if(myM1 <= 1.2):
                    gs = gs_1p2
                if((myM1 > 1.2) and (myM1 <= 3.5)):
                    gs = np.interp(np.log10(myM1), np.log10([1.2,3.5]),
                                   [gs_1p2,gs_3p5])
                if((myM1 > 3.5) and (myM1 <= 6.0)):
                    gs = np.interp(np.log10(myM1), np.log10([3.5,6.0]),
                                   [gs_3p5,gs_6])
                if(myM1 > 6.0):
                    gs = gs_6
                # Given Ftwin, gamma_smallq, and gamma_largeq at the specified
                # M1 & P, tabulate the cumulative mass ratio distribution
                # across 0.1 < q < 1.0
                # slope across 0.3 < q < 1.0
                fq = self.qv**gl
                # normalize to 0.3 < q < 1.0
                fq = fq / self._idl_tabulate(self.qv[indlq], fq[indlq])
                # add twins
                fq = fq * (1.0-Ftwin) + Heaviside * Ftwin
                # slope across 0.1 < q < 0.3
                fq[indsq] = (fq[indq0p3] * (self.qv[indsq]/0.3)**gs)
                # cumulative distribution
                cumfq = np.cumsum(fq) - fq[0]
                # normalize cumfq(q=1.0) = 1
                cumfq = cumfq / np.max(cumfq)
                # save to grid
                self.cumqdist[:,j,i] = cumfq
                # Given M1 and P, q_factor is the ratio of all binaries
                # 0.1 < q < 1.0 to those with 0.3 < q < 1.0
                q_factor = self._idl_tabulate(self.qv, fq)
                # Given M1 & P, calculate power-law slope eta of eccentricity
                # distribution
                if(mylogP >= 0.7):
                    # For log P > 0.7 use fits in Section 9.2.
                    # Power-law slope eta for M1 < 3 Msun and log P > 0.7
                    eta_3 = 0.6 - 0.7 / (mylogP-0.5)  #; Eqn. 17
                    # Power-law slope eta for M1 > 7 Msun and log P > 0.7
                    eta_7 = 0.9 - 0.2 / (mylogP-0.5)  #; Eqn. 18
                else:
                    # For log P < 0.7, set eta to fitted values at log P = 0.7
                    eta_3 = -2.9
                    eta_7 = -0.1
                # Given P, interpolate eta with respect to M1 at myM1
                if(myM1 <= 3.0):
                    eta = eta_3
                if((myM1 > 3.0) and (myM1 <= 7.0)):
                    eta = np.interp(np.log10(myM1), np.log10([3.,7.]),
                                    [eta_3, eta_7])
                if(myM1 > 7.0):
                    eta = eta_7
                # Given eta at the specified M1 and P, tabulate eccentricity
                # distribution
                if(10**mylogP <= 2.0):
                    # For P < 2 days, assume all systems are close to circular
                    # For adopted ev (spacing and minimum value), eta = -3.2
                    # satisfies this
                    fe = self.ev**(-3.2)
                else:
                    fe = self.ev**eta
                    # maximum eccentricity for given P
                    e_max = 1.0 - (10**mylogP/2.0)**(-2.0/3.0)
                    ind = np.where(self.ev >= e_max)
                    # set distribution = 0 for e > e_max
                    fe[ind] = 0.0                         
                    # Assume e distribution has power-law slope eta for
                    # 0.0 < e / e_max < 0.8 and then linear turnover between
                    # 0.8 < e / e_max < 1.0 so that distribution is continuous
                    # at e / e_max = 0.8 and zero at e = e_max
                    ind = np.where((self.ev >= 0.8*e_max)&(self.ev <= 1.0*e_max))
                    ind_cont = np.min(ind) - 1
                    fe[ind] = np.interp(self.ev[ind], [0.8*e_max,1.0*e_max],
                                        [fe[ind_cont],0.0])
                # cumulative distribution
                cumfe = np.cumsum(fe) - fe[0]
                # normalize cumfe(e=e_max) = 1
                cumfe = cumfe / np.max(cumfe)
                # save to grid
                self.cumedist[:,j,i] = cumfe
                # Given constants alpha and DlogP and M1 dependent values
                # flogPle1, flogPeq2p7, and flogPeq5p5, calculate frequency
                # flogP of companions with q > 0.3 per decade of orbital period
                # at given P (Section 9.3 and Eqn. 23)
                if(mylogP <= 1.0):
                    flogP = flogPle1
                if((mylogP > 1.0) and (mylogP <= 2.7-DlogP)):
                    flogP = flogPle1 + ((flogPeq2p7-flogPle1-alpha*DlogP)
                                        * (mylogP-1.0)/(1.7-DlogP))
                if((mylogP > 2.7-DlogP) and (mylogP <= 2.7+DlogP)):
                    flogP = flogPeq2p7 + alpha * (mylogP-2.7)
                if((mylogP > 2.7+DlogP) and (mylogP <= 5.5)):
                    flogP = (flogPeq2p7 + alpha * DlogP
                             + (flogPeq5p5-flogPeq2p7-alpha*DlogP)
                               *(mylogP-2.7-DlogP)/(2.8-DlogP))
                if(mylogP > 5.5):
                    flogP = flogPeq5p5 * np.exp(-0.3*(mylogP-5.5))
                # Convert frequency of companions with q > 0.3 to frequency of
                # companions with q > 0.1 according to q_factor; save to grid
                self.flogP_sq[j,i] = flogP * q_factor
                # Calculate prob. that a companion to M1 with period P is the
                # inner binary. Currently this is an approximation. 100% for
                # log P < 1.5; For log P > 1.5 adopt functional form that
                # reproduces M1 dependent multiplicity statistics in
                # Section 9.4, including a 41% binary star faction (59% single
                # star fraction) for M1 = 1 Msun and 96% binary star fraction
                # (4% single star fraction) for M1 = 28 Msun
                if(mylogP <= 1.5):
                    self.probbin[j,i] = 1.0
                else:
                    self.probbin[j,i] = (1.0 - 0.11 * (mylogP-1.5)**1.43
                                                    * (myM1/10.0)**0.56)
                if(self.probbin[j,i] <= 0.):
                    self.probbin[j,i] = 0.

            # Given M1, calculate cumulative binary period distribution
            mycumPbindist = (np.cumsum(self.flogP_sq[:,i]*self.probbin[:,i])
                             - self.flogP_sq[0,i] * self.probbin[0,i])
            # Normalize so that max(cumPbindist) = total binary star fraction
            # (NOT 1)
            mycumPbindist = (mycumPbindist / np.max(mycumPbindist)
                             * self._idl_tabulate(self.logPv,
                                                  self.flogP_sq[:,i]
                                                  *self.probbin[:,i]))
            # save to grid
            self.cumPbindist[:,i] = mycumPbindist

    def __repr__(self):
        return ("Moe and Di Stefano 2017 distributions on a grid of "
                f"n_M1={self.n_M1}, n_logP={self.n_logP}, n_q={self.n_q}, and "
                f"n_e={self.n_e}")

    def _repr_html_(self):
        return ("<h3>Moe and Di Stefano 2017 distributions on a grid of</h3>"
                f"<p>n_M1={self.n_M1}</p><p>n_logP={self.n_logP}</p>"
                f"<p>n_q={self.n_q}</p><p>n_e={self.n_e}</p>")

    def __call__(self, M1, M_min=0.08, M_max=150.0, all_binaries=True):
        """Initializing the class.

        Parameters
        ----------
        M1 : float or ndarray of float
            Mass of the primary in Msun.
        M_min : float (default: 0.08)
            Minimum mass of a star in Msun.
        M_max : float (default: 150.0)
            Maximum mass of a star in Msun.
        all_binaries : bool (default: True)
            If true set binary fration to 1.

        Returns
        -------
        M2, logP, e, Z : each being of same type as M1
            The four values are:
                - the mass of the secondary in Msun
                - the log10 of the period in days
                - the eccentricity
                - the metallicity
        """
        np.random.seed(1234567890) #TODO: remove
        M1s = np.atleast_1d(M1)
        M2s = []
        logPs = []
        es = []
        Zs = []
        # Implement Monte Carlo method / random number generator to select
        # single stars and binaries from the grids of distributions
        # Find index of M1v that is closest to M1.
        # For M1 = 40 - M1_max Msun, adopt binary statistics of M1 = 40 Msun.
        # For M1 = M_min - 0.8 Msun, adopt P and e dist of M1 = 0.8Msun, scale
        # and interpolate the companion frequencies so that the binary star
        # fraction of M1 = 0.08 Msun primaries is zero, and truncate the q
        # distribution so that q > q_min = M_min/M1
        for M1 in M1s:
            indM1 = np.where(abs(M1-self.M1v) == min(abs(M1-self.M1v)))
            indM1 = indM1[0]
            # Given M1, determine cumulative binary period distribution
            mycumPbindist = (self.cumPbindist[:,indM1]).flatten
            # If M1 < 0.8 Msun, rescale to appropriate binary star fraction
            if(M1 <= 0.8):
                mycumPbindist = mycumPbindist() * np.interp(np.log10(M1),
                                                            np.log10([M_min,
                                                                      0.8]),
                                                            [0.0, 1.0])
            # Given M1, determine the binary star fraction
            if all_binaries:
                mybinfrac = 1.0
            else:
                mybinfrac = np.max(mycumPbindist())
            # Generate random number myrand between 0 and 1
            myrand = np.random.rand()
            # If random number < binary star fraction, generate a binary
            if(myrand < mybinfrac):
                # Given myrand, select P and corresponding index in logPv
                mylogP = np.interp(myrand, mycumPbindist(), self.logPv)
                indlogP = np.where(abs(mylogP-self.logPv)
                                   == min(abs(mylogP-self.logPv)))
                indlogP = indlogP[0]
                # Given M1 & P, select e from eccentricity distribution
                mye = np.interp(np.random.rand(),
                                self.cumedist[:,indlogP,indM1].flatten(), self.ev)
                # Given M1 & P, determine mass ratio distribution.
                # If M1 < 0.8 Msun, truncate q distribution and consider only mass
                # ratios q > q_min = M_min / M1
                mycumqdist = self.cumqdist[:,indlogP,indM1].flatten()
                if(M1 < 0.8):
                    q_min = M_min / M1
                    # Calculate cumulative probability at q = q_min
                    cum_qmin = np.interp(q_min, self.qv, mycumqdist)
                    # Rescale and renormalize cumulative distribution for q > q_min
                    mycumqdist = mycumqdist - cum_qmin
                    mycumqdist = mycumqdist / max(mycumqdist)
                    # Set probability = 0 where q < q_min
                    indq = np.where(self.qv <= q_min)
                    mycumqdist[indq] = 0.0
                # Given M1 & P, select q from cumulative mass ratio distribution
                myq = np.interp(np.random.rand(), mycumqdist, self.qv)
            else:
                # If instead random number > binary star fraction, generate single
                # star
                myq = 0.0
                mylogP = np.nan
                mye = np.nan
            # Get metallicity
            Zsun = 0.02
            logZ = -2.3 + (0.176-(-2.3)) * np.random.rand()
            Z = Zsun*10**logZ
            M2s.append(M1*myq)
            logPs.append(mylogP)
            es.append(mye)
            Zs.append(Z)
        return np.array(M2s), np.array(logPs), np.array(es), np.array(Zs)


