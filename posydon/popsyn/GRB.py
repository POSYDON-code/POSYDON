# This code is currently not being actively used, and is not covered by unit tests.
# To test this file, remove it from the 'omit' blocks in POSYDON/setup.cfg.

"""Functions to compute gamma-ray burst (GRB) properties.

This code is currently not being actively used, and is not covered by unit
tests. To test this file, remove it from the 'omit' blocks in POSYDON/setup.cfg.
"""

__author__ = ['Simone Bavera <Simone.Bavera@unige.ch>']

import numpy as np

from posydon.utils.constants import Msun, clight
from posydon.utils.posydonwarning import Pwarn

GRB_PROPERTIES = ['GRB1', 'S1_eta', 'S1_f_beaming', 'S1_E_GRB', 'S1_E_GRB_iso', 'S1_L_GRB_iso',
                  'GRB2', 'S2_eta', 'S2_f_beaming', 'S2_E_GRB', 'S2_E_GRB_iso', 'S2_L_GRB_iso']

def get_GRB_properties(df, GRB_efficiency, GRB_beaming, E_GRB_iso_min=0.):
    """Compute GRB properties for a given binary population.

    Parameters
    ----------
    df : pandas DataFrame
        The binary population dataframe. It must contain the radiated disk
        mass columns 'S1_m_disk_radiated' and 'S2_m_disk_radiated'.
    GRB_efficiency : float
        The GRB energy conversion efficiency factor in [0, 1].
    GRB_beaming : float or str
        The GRB beaming factor. Either a float value in [0, 1] or a string
        selecting a beaming factor distribution ('Goldstein+15', 'Pescalli+15').
    E_GRB_iso_min : float, optional
        The minimum isotropic equivalent energy in erg to consider a GRB.

    Returns
    -------
    df : pandas DataFrame
        The binary population dataframe with the GRB properties added.
    """

    # check if radiated disk mass is avaialble to compute GRB properties
    if ('S1_m_disk_radiated' not in df or
        'S2_m_disk_radiated' not in df):
        raise ValueError('m_disk_radiated not found in the dataframe!')

    # reminder the user about the threshold cut
    if E_GRB_iso_min != 0.:
        Pwarn("We only consider GRBs with isotropic equivalent energy larger "
              f"than {E_GRB_iso_min} erg!", "ApproximationWarning")

    # define some methods

    def compute_eta(i, sel):
        """Set GRB energy conversion efficiency factor.

        Parameters
        ----------
        i : int
            The star component index (1 or 2).
        sel : pandas Series of bool
            Selection mask of the binaries to consider.
        """
        if not isinstance(GRB_efficiency, float):
            raise ValueError(f'GRB efficiency {GRB_efficiency} '
                             'must be float value in [0,1] range!')
        if GRB_efficiency < 0. or GRB_efficiency > 1.:
            raise ValueError(f'GRB efficiency {GRB_efficiency} '
                             'must be float value in [0,1] range!')
        # preset all efficiency factors to NaN and set them
        # only for the selected binaries emitting a GRB
        df[f'S{i}_eta'] = [np.nan] * df.shape[0]
        df.loc[sel,f'S{i}_eta'] = GRB_efficiency

    def compute_beaming_factor(i, sel):
        """Compute GRB beaming factor.

        Parameters
        ----------
        i : int
            The star component index (1 or 2).
        sel : pandas Series of bool
            Selection mask of the binaries to consider.
        """
        # TODO: check out Github copilot suggestions:
        # 'Fong+15', 'Ghirlanda+16', 'Lamb+19'?

        df[f'S{i}_f_beaming'] = [np.nan] * df.shape[0]
        if isinstance(GRB_beaming, float):
            df.loc[sel, f'S{i}_f_beaming'] = GRB_beaming
        elif GRB_beaming == 'Goldstein+15':
            n = df.loc[sel,f'S{i}_f_beaming'].shape[0]
            # generate a opening angle from a lognormal distribution
            # as in Table 3 of https://arxiv.org/pdf/1512.04464.pdf
            theta = 10**(np.random.normal(0.77,0.37,n))/180.*np.pi # radian
            # TODO: truncate the angle between 0 and pi/2
            # compute the solid angle, factor 2 to accounts for the two emispheres
            solid_angle = 4*np.pi*(1.-np.cos(theta))
            df.loc[sel,f'S{i}_f_beaming'] = solid_angle/(4*np.pi)
        # TODO: this option is not usable at the moment as it requires
        # redshift information which is not available in the dataframe
        # the refshift is compute by rate_calculation.py after GRB properties
        # are computed, fix this.
        # elif GRB_beaming == 'Lloyd-Ronning+19':
        #     if z is None:
        #         raise ValueError('Redshift not specified!')
        #     else:
        #         z = df.loc[sel,'z_formation']
        #     n = df.loc[sel,f'S{i}_f_beaming'].shape[0]
        #     # generate a opening angle from a lognormal distribution
        #     # as in Table 3 of https://arxiv.org/pdf/1512.04464.pdf
        #     theta =10**(np.random.normal(0.77*(1+z)**(-0.75),0.37,n))/180.*np.pi # radian
        #     # TODO: truncate the angle between 0 and pi/2
        #     # compute the solid angle, factor 2 to accounts for the two emispheres
        #     solid_angle = 4*np.pi*(1.-np.cos(theta))
        #     df.loc[sel,f'S{i}_f_beaming'] = solid_angle/(4*np.pi)
        elif GRB_beaming == 'Pescalli+15':
            # TODO: check the reference and the equation
            def f_B(E, eta):
                """Compute the beaming factor from Pescalli+15.

                Parameters
                ----------
                E : array
                    The GRB energy.
                eta : float
                    The GRB efficiency.

                Returns
                -------
                array
                    The beaming factors.
                """
                xi=0.5
                t=25.
                C = (7e46)**xi
                fb = np.ones(len(E))
                fb[E>0.] = (C*(eta*E[E>0.]/t)**(-xi))**(1./(1-xi))
                return fb
            E = df.loc[sel,f'S{i}_E_GRB']
            eta = df.loc[sel,f'S{i}_eta']
            df.loc[sel,f'S{i}_f_beaming'] = f_B(E,eta)
        else:
            raise ValueError(f'GRB beaming factor {GRB_beaming} not supported!')

    def compute_E_GRB(i, sel):
        """Compute the GRB energy.

        Parameters
        ----------
        i : int
            The star component index (1 or 2).
        sel : pandas Series of bool
            Selection mask of the binaries to consider.
        """
        df[f'S{i}_E_GRB'] = np.zeros(df.shape[0])
        m_disk_radiated = df.loc[sel,f'S{i}_m_disk_radiated']
        eta = df.loc[sel,f'S{i}_eta']
        df.loc[sel,f'S{i}_E_GRB'] = eta * m_disk_radiated * Msun * clight**2 # erg

    def compute_E_GRB_iso(i, sel):
        """Compute the GRB isotropic equivalent energy.

        Parameters
        ----------
        i : int
            The star component index (1 or 2).
        sel : pandas Series of bool
            Selection mask of the binaries to consider.
        """
        df[f'S{i}_E_GRB_iso'] = np.zeros(df.shape[0])
        E_GRB = df.loc[sel,f'S{i}_E_GRB']
        f_beaming = df.loc[sel,f'S{i}_f_beaming']
        df.loc[sel,f'S{i}_E_GRB_iso'] = E_GRB/f_beaming # erg

    def compute_L_GRB_iso(i, sel):
        """Compute the GRB isotropic equivalent luminosity.

        Parameters
        ----------
        i : int
            The star component index (1 or 2).
        sel : pandas Series of bool
            Selection mask of the binaries to consider.
        """
        df[f'S{i}_L_GRB_iso'] = np.zeros(df.shape[0])
        # compute the GRB isotropic equivalent luminosity assuming a
        # constant average timescale of 25s from Pescalli+15
        # TODO: improve this assumption
        E_GRB_iso = df.loc[sel,f'S{i}_E_GRB_iso']
        df.loc[sel,f'S{i}_L_GRB_iso'] = E_GRB_iso/25. # erg/s

    def flag_GRB_systems(i, sel):
        """Flag all systems emitting at least one GRB.

        Parameters
        ----------
        i : int
            The star component index (1 or 2).
        sel : pandas Series of bool
            Selection mask of the binaries to consider.
        """
        # flag all systems emitting at least one GRB
        df[f'GRB{i}'] = [False] * df.shape[0]
        df.loc[sel,f'GRB{i}'] = df.loc[sel,f'S{i}_E_GRB_iso'] >= E_GRB_iso_min


    # loop over the two star components and compute GRB properties
    # only for stars that formed a disk
    for i in range(1,3):
        # compute energies only for non-zero disk masses
        sel = df[f'S{i}_m_disk_radiated'] > 0.
        compute_eta(i, sel)
        compute_E_GRB(i, sel)
        compute_beaming_factor(i, sel)
        compute_E_GRB_iso(i, sel)
        compute_L_GRB_iso(i, sel)
        flag_GRB_systems(i, sel)

    return df
