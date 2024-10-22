
"""Determins XRB type and calculate X-ray luminosities."""


__authors__ = [
    "Devina Misra <devina.misra@unige.ch>",
]


import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import newton
from scipy.integrate import quad
from posydon.utils import constants as const
import copy
import warnings
from scipy.interpolate import PchipInterpolator
from posydon.utils.common_functions import eddington_limit


PATH_TO_POSYDON = os.environ.get("PATH_TO_POSYDON")


# def xrb_type(binary, idx=-1):
#     if binary.state in ['detached', 'RLO1', 'RLO2']:
#         if binary.star_1.state in ['NS', 'BH']:
#             accretor = binary.star_1
#             donor = binary.star_2
#             don_rel_RL = binary.rl_relative_overflow_2
#         elif binary.star_2.state in ['NS', 'BH']:
#             accretor = binary.star_2
#             donor = binary.star_1
#             don_rel_RL = binary.rl_relative_overflow_1
#         else:
#             xrb_type = 0
#             return xrb_type
#     else:
#         xrb_type = 0
#         return xrb_type

#     if don_rel_RL is None:
#         xrb_type = 0
#         return xrb_type

#     don_RL = 10**donor.log_R / (don_rel_RL + 1)
#     if ( (donor.state == 'H-rich_Core_H_burning') & (donor.mass >= 3.0) &
#         (binary.state == 'detached') & (don_RL <= (100.0 * 10 ** donor.log_R) ) & 
#         (donor.surf_avg_omega_div_omega_crit >= 0.7) & (binary.orbital_period >= 10.0) & (binary.orbital_period <= 300.0) ) :
#             xrb_type = 3
#     elif ((binary.state =='RLO1') & (binary.state =='RLO2')):
#             xrb_type = 1
#             print(binary.state , donor.state, accretor.state)
#     else:
#             xrb_type = 2

#     return xrb_type


def x_ray_luminosity(binary, idx=-1):
    """ Calculate the geometrical beaming of a super-Eddington accreting source.
    Compute the super-Eddington isotropic-equivalent accretion rate and the
    beaming factor of a star. This does not change the intrinsic accretion onto
    the accretor and is an observational effect due to the inflated structure
    of the accretion disc that beams the outgoing X-ray emission. This is
    important for observing super-Eddington X-ray sources
    (e.g. ultraluminous X-ray sources). In the case of a BH we are assuming that it
    has zero spin which is not a good approximation for high accretion rates.
    Parameters
    ----------
    binary : BinaryStar
        The binary object.
    Returns
    -------
    list
        The super-Eddington isotropic-equivalent accretion rate and beaming
        factor respcetively in solar units.
    References
    ----------
    .. [1] Shakura, N. I. & Sunyaev, R. A. 1973, A&A, 24, 337
    .. [2] King A. R., 2008, MNRAS, 385, L113
    """
    if binary.star_1.state in ['NS', 'BH']:
        accretor = binary.star_1
        donor = binary.star_2
        don_rel_RL = binary.rl_relative_overflow_2
    elif binary.star_2.state in ['NS', 'BH']:
        accretor = binary.star_2
        donor = binary.star_1
        don_rel_RL = binary.rl_relative_overflow_1
    else:
        return 0.0, 1.0

    if ((binary.lg_mtransfer_rate is None) and (accretor.lg_mdot is None)):
        return 0.0, 1.0
    elif (binary.lg_mtransfer_rate is not None):
        mdot = 10**binary.lg_mtransfer_rate
    else:
        mdot = 10**accretor.lg_mdot

    if don_rel_RL is None:
        return 0.0, 1.0

    mdot_edd = eddington_limit(binary, idx=-1)[0]
    eta = eddington_limit(binary, idx=-1)[1]
    don_RL = 10**donor.log_R / (don_rel_RL + 1)

    if ( (donor.state == 'H-rich_Core_H_burning') & (donor.mass >= 3.0) &
        (binary.state == 'detached') & (don_RL <= (100.0 * 10 ** donor.log_R) ) & 
        (donor.surf_avg_omega_div_omega_crit >= 0.7) & (binary.orbital_period >= 10.0) & (binary.orbital_period <= 300.0) ): # Be-XRBs
        Lx = 10.0**35 * 10.0**(4.53 - (1.50 * np.log10(binary.orbital_period)) )
        b = 1.0
    else: # RLO and wind
        if mdot >= mdot_edd:
            if (mdot > 150.0 * mdot_edd):
                b = 3.2e-3
            elif (mdot > 8.5 * mdot_edd):
                b = 73. / (mdot / mdot_edd)**2
            else:
                b = 1.
            Lx = mdot_edd*eta * const.clight**2*(1 + np.log((10**binary.lg_mtransfer_rate) / mdot_edd))/b
        else:
            b = 1.
            Lx = 10**accretor.lg_mdot*eta*const.clight**2
    return Lx, b


# def beaming(binary):
#     """Calculate the geometrical beaming of a super-Eddington accreting source.

#     Compute the super-Eddington isotropic-equivalent accretion rate and the
#     beaming factor of a star. This does not change the intrinsic accretion onto
#     the accretor and is an observational effect due to the inflated structure
#     of the accretion disc that beams the outgoing X-ray emission. This is
#     important for observing super-Eddington X-ray sources
#     (e.g. ultraluminous X-ray sources). In case of a BH we are assuming that it
#     has zero spin which is not a good approximation for high accretion rates.

#     Parameters
#     ----------
#     binary : BinaryStar
#         The binary object.

#     Returns
#     -------
#     list
#         The super-Eddington isotropic-equivalent accretion rate and beaming
#         factor respcetively in solar units.

#     References
#     ----------
#     .. [1] Shakura, N. I. & Sunyaev, R. A. 1973, A&A, 24, 337
#     .. [2] King A. R., 2008, MNRAS, 385, L113

#     """
#     mdot_edd = eddington_limit(binary, idx=-1)[0]

#     rlo_mdot = 10**binary.lg_mtransfer_rate

#     if rlo_mdot >= mdot_edd:
#         if rlo_mdot > 8.5 * mdot_edd:
#             # eq. 8 in King A. R., 2009, MNRAS, 393, L41-L44
#             b = 73 / (rlo_mdot / mdot_edd)**2
#         else:
#             b = 1
#         # Shakura, N. I. & Sunyaev, R. A. 1973, A&A, 24, 337
#         mdot_beam = mdot_edd * (1 + np.log(rlo_mdot / mdot_edd)) / b
#     else:
#         b = 1
#         mdot_beam = 10**binary.lg_mtransfer_rate

#     return mdot_beam, b