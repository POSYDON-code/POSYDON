__author__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Max Briel <max.briel@unige.ch>",
]

from posydon.utils.constants import Zsun

from astropy.cosmology import Planck15 as cosmology
from astropy import constants as const
import numpy as np
import scipy as sp
from astropy.cosmology import z_at_value
from scipy.interpolate import interp1d
from astropy import units as u


DEFAULT_MODEL = {
    "delta_t": 100,  # Myr
    "SFR": "IllustrisTNG",
    "sigma_SFR": None,
    "Z_max": 5.0,
    "select_one_met": False,
    "dlogZ": None,  # e.g, [np.log10(0.0142/2),np.log10(0.0142*2)]
    "Zsun": Zsun,
    "compute_GRB_properties": False,
    "GRB_beaming": 1.0,  # e.g., 0.5, 'Goldstein+15'
    "GRB_efficiency": 0.0,  # e.g., 0.01
    "E_GRB_iso_min": 0.0,  # e.g., 1e51 erg
}


def get_shell_comoving_volume(z_hor_i, z_hor_f, sensitivity="infinite"):
    """Compute comoving volume corresponding to a redshift shell.

    Parameters
    ----------
    z_hor_i : double
        Cosmological redshift. Lower bound of the integration.
    z_hor_f : double
        Cosmological redshift. Upper bound of the integration.
    sensitivity : string
        hoose which GW detector sensitivity you want to use. At the moment
        only 'infinite' is available, i.e. p_det = 1.

    Returns
    -------
    double
        Retruns the comoving volume between the two shells z_hor_i
        and z_hor_f in Gpc^3.

    """
    c = const.c.to("Gpc/yr").value  # Gpc/yr
    H_0 = cosmology.H(0).to("1/yr").value  # km/Gpc*s

    def E(z):
        Omega_m = cosmology.Om0
        Omega_L = 1 - cosmology.Om0
        return np.sqrt(Omega_m * (1.0 + z) ** 3 + Omega_L)

    def f(z, sensitivity):
        if sensitivity == "infinite":
            return (
                1.0
                / (1.0 + z)
                * 4
                * np.pi
                * c
                / H_0
                * (get_comoving_distance_from_redshift(z) * 10 ** (-3.0)) ** 2.0
                / E(z)
            )
        else:
            # TODO: peanut-shaped antenna patter comoving volume calculation
            raise ValueError("Sensitivity not supported!")

    return sp.integrate.quad(f, z_hor_i, z_hor_f, args=(sensitivity))[0]  # Gpc^3


def get_comoving_distance_from_redshift(z):
    """Compute the comoving distance from redshift.

    Parameters
    ----------
    z : double
        Cosmological redshift.

    Returns
    -------
    double
        Comoving distance in Mpc corresponding to the redhisft z.

    """
    return cosmology.comoving_distance(z).value  # Mpc


def get_cosmic_time_from_redshift(z):
    """Compute the cosmic time from redshift.

    Parameters
    ----------
    z : double
        Cosmological redshift.

    Returns
    -------
    double
        Return age of the cosmic time in Gyr given the redshift z.

    """
    return cosmology.age(z).value  # Gyr


def redshift_from_cosmic_time_interpolator():
    """Interpolator to compute the cosmological redshift given the cosmic time.

    Returns
    -------
    object
        Returns the trained interpolator object.

    """
    # astropy z_at_value method is too slow to compute z_mergers efficinty
    # we must implement interpolation
    t = np.linspace(1e-2, cosmology.age(1e-08).value * 0.9999999, 1000)
    z = np.zeros(1000)
    for i in range(1000):
        z[i] = z_at_value(cosmology.age, t[i] * u.Gyr)
    f_z_m = interp1d(t, z, kind="cubic")
    return f_z_m


def get_redshift_from_cosmic_time(t_cosm):
    """Compute the cosmological redshift given the cosmic time..

    Parameters
    ----------
    t_cosm : array doubles
        Cosmic time to which you want to know the redhisft.

    Returns
    -------
    array doubles
        Cosmolgocial redshift corresponding to the cosmic time.

    Note
    ----
    - The function uses the interpolator redshift_from_cosmic_time_interpolator,
    which is created each time the function is called. `z_at_value` from astropy
    can be used for single values, but it is too slow for arrays.
    """
    return redshift_from_cosmic_time_interpolator(t_cosm)


def get_redshift_bin_edges(delta_t):
    """Compute the redshift bin edges.

    Parameters
    ----------
    delta_t : double
        Time interval in Myr.

    Returns
    -------
    array doubles
        Redshift bin edges.

    """
    n_redshift_bin_centers = int(cosmology.age(0).to("Myr").value / delta_t)
    # generate t_birth at the middle of each self.delta_t bin
    t_birth_bin = [cosmology.age(0.0).value]
    for i in range(n_redshift_bin_centers + 1):
        t_birth_bin.append(t_birth_bin[i] - delta_t * 1e-3)  # Gyr
    # compute the redshift
    z_birth_bin = []
    for i in range(n_redshift_bin_centers):
        # do not count first edge, we add z=0. later
        z_birth_bin.append(z_at_value(cosmology.age, t_birth_bin[i + 1] * u.Gyr))
    # add the first and last bin edge at z=0. and z=inf.=100
    z_birth_bin = np.array([0.0] + z_birth_bin + [100.0])

    return z_birth_bin


def get_redshift_bin_centers(delta_t):
    """Compute the redshift bin centers.

    Parameters
    ----------
    delta_t : double
        Time interval in Myr.

    Returns
    -------
    array doubles
        Redshift bin centers.
    """
    # generate t_birth at the middle of each self.delta_t bin
    t_birth_bin = [cosmology.age(0.0).value]
    t_birth = []
    n_redshift_bin_centers = int(cosmology.age(0).to("Myr").value / delta_t)
    for i in range(n_redshift_bin_centers + 1):
        t_birth.append(t_birth_bin[i] - delta_t * 1e-3 / 2.0)  # Gyr
        t_birth_bin.append(t_birth_bin[i] - delta_t * 1e-3)  # Gyr
    t_birth = np.array(t_birth)
    # compute the redshift
    z_birth = []
    for i in range(n_redshift_bin_centers + 1):
        # z_at_value is from astopy.cosmology
        z_birth.append(z_at_value(cosmology.age, t_birth[i] * u.Gyr))
    z_birth = np.array(z_birth)

    return z_birth
