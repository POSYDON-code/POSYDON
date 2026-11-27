"""Natal kicks of stars after their supernova.

The kick velocity, its direction and the position of the star in its
orbit at the moment of the supernova. The orbit after the kick is
computed in orbit.py.
"""

from abc import ABC, abstractmethod

import numpy as np
import scipy as sp

from posydon.utils.common_functions import is_None_or_nan
from posydon.utils.posydonwarning import Pwarn


def get_kick_prescription(**kwargs):
    """Return the Natal_Kick object for the given SN step settings."""
    return Natal_Kick(**kwargs)


# kick_normalisation values that are now kick_prescription options
DEPRECATED_KICK_NORMALISATIONS = ("asym_ej", "linear")


def replace_deprecated_kick_normalisation(config):
    """Return config with a deprecated kick_normalisation replaced.

    A kick_normalisation of 'asym_ej' or 'linear' becomes the
    kick_prescription, and the normalisation is set to 'one'.
    """
    option = config['kick_normalisation']
    if option not in DEPRECATED_KICK_NORMALISATIONS:
        return config
    Pwarn("kick_normalisation 'asym_ej' and 'linear' are deprecated, use "
          "kick_prescription instead. Setting kick normalization to unity.",
          "DeprecationWarning")
    new_config = config.copy()
    new_config['kick_prescription'] = option
    new_config['kick_normalisation'] = 'one'
    return new_config


class Kick_Velocity_Base(ABC):
    """Kick velocity before normalisation (option kick_prescription)."""
    def __init__(self,**kwargs):
        self.RNG = kwargs['RNG']
        self.verbose = kwargs['verbose']

    def __call__(self,star,sigma,mean):
        """Return the kick velocity of the star in km/s."""
        velocity = self._draw(star,sigma,mean)
        if self.verbose:
            print(f"{type(self).__name__}: kick velocity before "
                  f"normalisation {velocity:.2f} km/s")
        return velocity

    @abstractmethod
    def _draw(self,star,sigma,mean):
        """Draw the kick velocity of the star in km/s."""
        raise NotImplementedError

    def __repr__(self):
        return f"{type(self).__name__}()"


class Kick_Velocity_Maxwellian(Kick_Velocity_Base):
    """Maxwellian distribution with dispersion sigma.

    Hobbs, G., Lorimer, D. R., Lyne, A. G., & Kramer, M. (2005), MNRAS, 360,
    974, with sigma = 265 km/s.
    """
    def _draw(self,star,sigma,mean):
        # sigma is None only when no kick is drawn, so this is a fallback
        if sigma is None:
            sigma = 265.0
        return sp.stats.maxwell.rvs(loc=0., scale=sigma, size=1,
                                    random_state=self.RNG)[0]


class Kick_Velocity_LogNormal(Kick_Velocity_Base):
    """Log-normal distribution with shape sigma and median mean.

    Disberg, P. & Mandel, I. (2025), arXiv:2505.22102.
    """
    def _draw(self,star,sigma,mean):
        # sigma is None only when no kick is drawn, so this is a fallback
        if sigma is None:
            sigma = 0.68
        if mean is None:
            mean = np.exp(5.60)
        return sp.stats.lognorm.rvs(s=sigma, scale=mean, size=1,
                                    random_state=self.RNG)[0]


class Kick_Velocity_AsymEj(Kick_Velocity_Base):
    """Kick from asymmetric ejecta (gravitational tug-boat mechanism).

    Janka, H.-T. (2017), ApJ, 837, 84.
    """
    def _draw(self,star,sigma,mean):
        # fraction of the SN explosion energy that is kinetic energy of the gas
        f_kin = 0.1
        beta = 0.1      # fraction of the ejecta mass that is neutrino heated
        epsilon = 1
        alpha_ej = 0.01
        M_NS = star.mass                                # neutron star mass
        M_rembar = (((3*M_NS/20 + 1)**2) - 1)/0.3       # remnant baryonic mass
        M_ej = abs(star.mass_history[-1] - M_rembar)    # ejecta mass
        return (211*(f_kin*beta*epsilon)**(1/2)*(alpha_ej/0.1)*(M_ej/0.1)
                * (M_NS/1.5)**(-1))


class Kick_Velocity_Linear(Kick_Velocity_Base):
    """Kick linear in the ejecta-to-remnant mass ratio.

    Richards, S. M., Eldridge, J. J., Briel, M. M., Stevance, H. F., &
    Willcox, R. (2022), arXiv:2208.02407.
    """
    def _draw(self,star,sigma,mean):
        M_ej = star.co_core_mass_history[-1] - star.mass    # ejecta mass
        M_rem = star.mass                                   # remnant mass
        # best-fit parameters
        alpha = 115
        beta = 15
        return alpha * (M_ej/M_rem) + beta


KICK_VELOCITY_PRESCRIPTIONS = {
    "maxwellian": Kick_Velocity_Maxwellian,
    "log_normal": Kick_Velocity_LogNormal,
    "asym_ej": Kick_Velocity_AsymEj,
    "linear": Kick_Velocity_Linear,
}


# Factor multiplying the kick velocity (option kick_normalisation).
# 1 - f_fb is Eq. 21 of Fryer et al. (2012), ApJ, 749, 91.
def norm_one_minus_fallback(star):
    """1 - f_fb."""
    return 1.0 - star.f_fb

def norm_one_over_mass(star):
    """1.4 Msun / M for BHs, 1 for NSs."""
    if star.state == 'BH':
        return 1.4 / star.mass
    return 1.0

def norm_NS_one_minus_fallback_BH_one(star):
    """1 - f_fb for NSs, 1 for BHs."""
    if star.state == 'BH':
        return 1.0
    return 1.0 - star.f_fb

def norm_one(star):
    """1: no normalisation."""
    return 1.0

def norm_zero(star):
    """0: no natal kicks."""
    return 0.0

KICK_NORMALISATIONS = {
    "one_minus_fallback": norm_one_minus_fallback,
    "one_over_mass": norm_one_over_mass,
    "NS_one_minus_fallback_BH_one": norm_NS_one_minus_fallback_BH_one,
    "one": norm_one,
    "zero": norm_zero,
}


class Natal_Kick:
    """Natal kick of a star after its supernova.

    The kick velocity is drawn from the distribution set by
    kick_prescription, with sigma and mean depending on the SN type and the
    remnant, and multiplied by the factor set by kick_normalisation. Its
    direction (polar and azimuthal angles) and the position in the orbit
    (mean anomaly) are drawn uniformly. Values already set on the star
    (natal_kick_velocity, natal_kick_azimuthal_angle, natal_kick_polar_angle,
    natal_kick_mean_anomaly) are used instead of drawing them, and drawn
    values are stored on the star. Values that are None or NaN count as not
    set.

    Kicks are only drawn for ECSN, CCSN, PPISN and PISN; WDs get no kick and
    are handled by the SN step. With kick=False the kick velocity is 0 (also
    if one was set on the star), but the angles and the mean anomaly are
    still used: the instantaneous SN mass loss changes the orbit anyway.
    """
    def __init__(self,**kwargs):
        self.RNG = kwargs['RNG']
        self.kick_prescription = kwargs['kick_prescription']
        self.kick_normalisation = kwargs['kick_normalisation']
        if self.kick_prescription not in KICK_VELOCITY_PRESCRIPTIONS:
            raise ValueError(
                f"Invalid kick_prescription {self.kick_prescription!r}. "
                "Valid options are: "
                f"{', '.join(KICK_VELOCITY_PRESCRIPTIONS)}.")
        if self.kick_normalisation not in KICK_NORMALISATIONS:
            raise ValueError(
                f"Invalid kick_normalisation {self.kick_normalisation!r}. "
                f"Valid options are: {', '.join(KICK_NORMALISATIONS)}.")
        velocity = KICK_VELOCITY_PRESCRIPTIONS[self.kick_prescription]
        self._velocity = velocity(**kwargs)
        self._normalisation = KICK_NORMALISATIONS[self.kick_normalisation]
        self.kick = kwargs['kick']
        self.no_kick_direct_collapse = kwargs['no_kick_direct_collapse']
        self.verbose = kwargs['verbose']
        # (sigma, mean) of the kick distribution by SN type / remnant
        self._sigma_mean = {
            "ECSN": (kwargs['sigma_kick_ECSN'], kwargs['mean_kick_ECSN']),
            "NS": (kwargs['sigma_kick_CCSN_NS'], kwargs['mean_kick_CCSN_NS']),
            "BH": (kwargs['sigma_kick_CCSN_BH'], kwargs['mean_kick_CCSN_BH']),
        }

    def __call__(self,star):
        """Return (v_kick, cos_theta, phi, mean_anomaly) of the star's kick.

        v_kick is in km/s, cos_theta is the cosine of the polar angle, and
        phi and mean_anomaly are in radians (see orbit.py for the geometry).
        Random values are always drawn in the same order: velocity,
        azimuthal angle, polar angle, mean anomaly.
        """
        if not self.kick:
            star.natal_kick_velocity = 0.0      # no natal kicks at all
            source = "kick=False"
        elif is_None_or_nan(star.natal_kick_velocity):
            star.natal_kick_velocity = self._draw_velocity(star)
            source = "drawn"
        else:
            source = "set on the star"
        v_kick = star.natal_kick_velocity

        if is_None_or_nan(star.natal_kick_azimuthal_angle):
            star.natal_kick_azimuthal_angle = self.RNG.uniform(0, 2 * np.pi)
        phi = star.natal_kick_azimuthal_angle

        if is_None_or_nan(star.natal_kick_polar_angle):
            cos_theta = self.RNG.uniform(-1, 1)
            star.natal_kick_polar_angle = np.arccos(cos_theta)
        else:
            cos_theta = np.cos(star.natal_kick_polar_angle)

        if is_None_or_nan(star.natal_kick_mean_anomaly):
            star.natal_kick_mean_anomaly = self.RNG.uniform(0, 2 * np.pi)
        elif not isinstance(star.natal_kick_mean_anomaly, float):
            raise ValueError(
                "mean_anomaly must be a single float value."
                f"\n mean_anomaly = {star.natal_kick_mean_anomaly}")
        mean_anomaly = star.natal_kick_mean_anomaly

        if self.verbose:
            print(f"Natal_Kick: v_kick={v_kick:.2f} km/s ({source}), polar "
                  f"angle={star.natal_kick_polar_angle:.3f}, azimuthal "
                  f"angle={phi:.3f}, mean anomaly={mean_anomaly:.3f}")
        return v_kick, cos_theta, phi, mean_anomaly

    def _draw_velocity(self,star):
        """Draw the kick velocity of the star in km/s."""
        if star.SN_type == "ECSN":
            sigma, mean = self._sigma_mean["ECSN"]
        elif star.SN_type in ("CCSN", "PPISN", "PISN"):
            if star.state == "massless_remnant":
                self._report("no kick on a massless remnant")
                return 0.0
            if star.state not in ("NS", "BH"):
                raise ValueError("CCSN/PPISN/PISN only for NS/BH.")
            if (self.no_kick_direct_collapse and star.state == "BH"
                    and star.f_fb == 1.0):
                self._report("no kick for a direct-collapse BH")
                return 0.0
            sigma, mean = self._sigma_mean[star.state]
        else:
            raise ValueError(f"Natal kicks are only drawn for ECSN, CCSN, "
                             f"PPISN and PISN, not for SN type "
                             f"{star.SN_type!r}.")
        if sigma is None:
            self._report(f"no kick: sigma is None for this {star.SN_type}")
            return 0.0
        factor = self._normalisation(star)
        self._report(f"{star.SN_type} {star.state}: sigma={sigma}, "
                     f"mean={mean}, normalisation factor {factor:.3f}")
        return factor * self._velocity(star, sigma, mean)

    def _report(self,message):
        if self.verbose:
            print(f"Natal_Kick: {message}")

    def __repr__(self):
        return (f"Natal_Kick(kick_prescription={self.kick_prescription!r}, "
                f"kick_normalisation={self.kick_normalisation!r})")
