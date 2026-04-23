"""Module for converting a MESA history file to an EEP track.

Reference: Dotter, Aaron (2016), AJSS, 222, 1.
"""


__authors__ = [
    "Aaron Dotter <aaron.dotter@gmail.com>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
]


import gzip

import numpy as np
from scipy.interpolate import pchip

# suggested lists of EEPs
ZAMSlo = ['ZAMS', 'IAMS', 'TAMS', 'TRGB', 'ZACHEB', 'TACHEB']
ZAMShi = ['ZAMS', 'IAMS', 'TAMS', 'ZACHEB', 'TACHEB', 'CBURN']
HeZAMS = ['ZACHEB', 'TACHEB', 'CBURN']


class EEP:
    """Convert a MESA history file to Equivalent Evolutionary Phase track."""

    def __init__(self, filename, EEP_NAMES=ZAMShi, EEP_INTERVAL=100):
        """Load an MESA history file and construct the EEP instance."""
        self.filename = filename.strip()
        try:
            if self.filename.endswith(".gz"):
                f = gzip.open(self.filename, "rt")
            else:
                f = open(self.filename, "r")

            self.header1 = f.readline()
            self.header2 = f.readline()
            self.header3 = f.readline()
            f.close()
                # this is compatible with `.gz` files
            tr = np.genfromtxt(self.filename, names=True, skip_header=5)
            names = tr.dtype.names
        except IOError:
            print("Failed to open: ")
            print(self.filename)

        # this section attempts to find each of the EEPs for the track
        prems = self._PreMS(tr)
        zams = self._ZAMS(tr)
        iams = self._IAMS(tr, Xc=0.2, guess=zams+1)
        tams = self._TAMS(tr, guess=iams+1)
        trgb = self._TRGB(tr, guess=tams+1)
        zacheb = self._ZACHEB(tr, guess=trgb+1)
        tacheb = self._TACHEB(tr, guess=zacheb+1)
        tpagb = self._TPAGB(tr, guess=tacheb+1)
        pagb = self._PAGB(tr, guess=tpagb+1)
        wdcs = self._WDCS(tr, guess=pagb+1)
        cburn = self._CBurn(tr, guess=zacheb+1)

        # compute the distance metric along the track that is used to assign
        # secondary EEPs
        metric = self._metric_function(tr)

        eep_index = []

        # if the EEP is in the input list, and it exists (>0)
        # then add it to the official list of EEPs
        for eep in EEP_NAMES:
            if eep == 'PreMS' and prems >= 0:
                eep_index.append(prems)
            if eep == 'ZAMS' and zams >= 0:
                eep_index.append(zams)
            if eep == 'IAMS' and iams >= 0:
                eep_index.append(iams)
            if eep == 'TAMS' and tams >= 0:
                eep_index.append(tams)
            if eep == 'TRGB' and trgb >= 0:
                eep_index.append(trgb)
            if eep == 'ZACHEB' and zacheb >= 0:
                eep_index.append(zacheb)
            if eep == 'TACHEB' and tacheb >= 0:
                eep_index.append(tacheb)
            if eep == 'TPAGB' and tpagb >= 0:
                eep_index.append(tpagb)
            if eep == 'PAGB' and pagb >= 0:
                eep_index.append(pagb)
            if eep == 'WDCS' and wdcs >= 0:
                eep_index.append(wdcs)
            if eep == 'CBURN' and cburn >= 0:
                eep_index.append(cburn)

        # some bookkeeping
        eep_counter = 0
        self.num_primary = len(eep_index)
        self.num_secondary = EEP_INTERVAL*(len(eep_index)-1)
        self.num_eeps = self.num_primary + self.num_secondary
        self.eeps = np.zeros((self.num_eeps), tr.dtype.descr)

        # assign the primary EEPs in the correct location for the EEP track
        for eep in range(self.num_primary):
            self.eeps[eep_counter] = tr[eep_index[eep]]
            eep_counter += EEP_INTERVAL + 1

        # assign the secondary EEPs in the correct location for the EEP track
        for eep in range(1, self.num_primary):
            eep_counter = (eep-1) * (EEP_INTERVAL + 1)
            # fill in secondary
            lo = eep_index[eep-1]
            hi = eep_index[eep]
            dm = (metric[hi] - metric[lo])/(EEP_INTERVAL+1)
            m = metric[lo].item()
            for j in range(1, EEP_INTERVAL+1):
                m += dm
                y = []
                for i, name in enumerate(names):
                    y.append(pchip(x=metric[lo:hi], y=tr[name][lo:hi])(m))
                self.eeps[eep_counter + j] = tuple(y)

    # the following function definitions are for primary EEP determinations
    def _PreMS(self, tr, Dfrac=0.01, guess=0):
        PreMS = -1
        for i in range(len(tr)):
            if tr['center_h2'][i] < Dfrac*tr['center_h2'][0]:
                PreMS = i
                break
        return PreMS

    def _ZAMS(self, tr, dXc=0.001, guess=0):
        ZAMS = -1
        for i in range(len(tr)):
            if abs(tr['center_h1'][i]-tr['center_h1'][0]) > dXc:
                ZAMS = i
                break
        return ZAMS

    def _IAMS(self, tr, Xc=0.1, guess=0):
        IAMS = -1
        for i in range(len(tr)):
            if tr['center_h1'][i] < Xc:
                IAMS = i
                break
        return IAMS

    def _TAMS(self, tr, Xc=0.00001, guess=0):
        TAMS = -1
        for i in range(len(tr)):
            if tr['center_h1'][i] < Xc:
                TAMS = i
                break
        return TAMS

    def _TRGB(self, tr, guess=0):
        Yc_min = tr['center_he4'][guess] - 0.01
        L_He_max = -99.
        Tc_min = 99.
        TRGB = -1
        TRGB1 = 0
        TRGB2 = 0
        for i in range(guess, len(tr)):
            if tr['center_he4'][i] > Yc_min:
                if tr['log_LHe'][i] > L_He_max:
                    L_He_max = tr['log_LHe'][i]
                    TRGB1 = i
                if tr['log_center_T'][i] < Tc_min:
                    Tc_min = tr['log_center_T'][i]
                    TRGB2 = i
        return max(TRGB, min(TRGB1, TRGB2))

    def _ZACHEB(self, tr, guess=0):
        ZACHEB = -1
        Yc_min = max(0.9, tr['center_he4'][guess] - 0.03)
        L_He_max = -99.
        Tc_min = 99.
        ZACHEB1 = 0
        ZACHEB = 0
        for i in range(guess, len(tr)):
            if tr['center_he4'][i] > Yc_min and tr['log_LHe'][i] > L_He_max:
                L_He_max = tr['log_LHe'][i]
                ZACHEB1 = i
        for i in range(ZACHEB1, len(tr)):
            if tr['center_he4'][i] > Yc_min and tr['log_center_T'][i] < Tc_min:
                Tc_min = tr['log_center_T'][i]
                ZACHEB = i
        return ZACHEB

    def _TACHEB(self, tr, Yc_min=0.001, guess=0):
        TACHEB = -1
        for i in range(guess, len(tr)):
            if tr['center_he4'][i] < Yc_min:
                TACHEB = i
                break
        return TACHEB

    def _TPAGB(self, tr, guess=0):
        TPAGB = -1
        He_shell_min = 0.1
        Yc_min = 1.0e-6
        for i in range(guess, len(tr)):
            He_shell_mass = tr['he_core_mass'][i] - tr['c_core_mass'][i]
            if tr['center_he4'][i] < Yc_min and He_shell_mass < He_shell_min:
                TPAGB = i
                break
        return TPAGB

    def _PAGB(self, tr, guess=0):
        PAGB = -1
        core_mass_frac_min = 0.8
        Tc_now = tr['log_center_T'][guess]
        Tc_end = tr['log_center_T'][-1]

        # check for low-inter / high mass split
        if Tc_now > Tc_end:
            for i in range(guess, len(tr)):
                core_mass_frac = tr['c_core_mass'][i] / tr['star_mass'][i]
                if core_mass_frac > core_mass_frac_min:
                    PAGB = i
                    break

        return PAGB

    def _WDCS(self, tr, gamma=10., guess=0):
        WDCS = -1
        for i in range(guess, len(tr)):
            if tr['center_gamma'][i] > gamma:
                WDCS = i
                break
        return WDCS

    def _CBurn(self, tr, XC12=0.1, guess=0):
        CBURN = -1
        XY_min = 1.0E-6
        for i in range(guess, len(tr)):
            Xc = tr['center_h1'][i]
            Yc = tr['center_he4'][i]
            C12 = tr['center_c12'][i]
            if Xc < XY_min and Yc < XY_min and C12 < XC12:
                CBURN = i
                break
        return CBURN

    # this function computes the distance metric along the evolutionary track
    # it is made up of several terms whose weights can be adjusted.  Currently
    # use the H-R and age information only.
    # other terms can be added, must be "monotonic increasing"
    def _metric_function(self, tr):
        term1 = tr['log_Teff']
        term2 = tr['log_L']
        term3 = np.log10(tr['star_age'])
        term4 = tr['log_center_Rho']
        weight1 = 2.0
        weight2 = 0.125
        weight3 = 1.0
        weight4 = 0.0
        # etc.

        metric = np.zeros(len(tr))

        for i in range(1, len(tr)):
            metric[i] = metric[i-1] + \
                weight1*pow(term1[i]-term1[i-1], 2) + \
                weight2*pow(term2[i]-term2[i-1], 2) + \
                weight3*pow(term3[i]-term3[i-1], 2) + \
                weight4*pow(term4[i]-term4[i-1], 2)

        return metric
