"""Classical X-ray-binary accretion and luminosity prescriptions.

The functions in this module operate on physical properties rather than on
``BinaryStar`` objects.  Inputs may be scalars or broadcast-compatible NumPy
arrays.  Unless noted otherwise, masses, radii, luminosities, periods, mass
rates, and X-ray luminosities are expressed in ``Msun``, ``Rsun``, ``Lsun``,
days, ``Msun yr^-1``, and ``erg s^-1``, respectively.

This deliberately small module implements classical Roche-lobe-overflow,
wind-fed, and Be-XRB prescriptions.  It does not select X-ray binaries from a
population, apply an instrumental bandpass, or implement advanced  accretion models
(e.g. GRRMHD, magnetic neutron-star etc).

Invalid numerical elements are returned as ``numpy.nan`` so that one bad
population row does not abort a vectorized calculation. Invalid named
prescriptions and scalar configuration values raise :class:`ValueError`.
"""

import numpy as np

from posydon.utils import constants as const

__authors__ = [
    "Anastasios Fragos <Anastasios.Fragkos@unige.ch>",
    "Arnaud Aguet <Arnaud.Aguet@etu.unige.ch>",
    "Devina Misra <devina.misra@ntnu.no>"
]

__all__ = [
    "accretion_luminosity",
    "be_xray_luminosity",
    "black_hole_radiative_efficiency",
    "bondi_hoyle_accretion_rate",
    "eddington_limit_from_properties",
    "neutron_star_radiative_efficiency",
    "wind_velocity",
]


# TODO(POSYDON#886): add a worked population-synthesis tutorial for this API.


def _return_scalar_if_scalar(value, scalar):
    """Return a Python scalar when every input to a calculation was scalar."""
    array = np.asarray(value)
    if scalar:
        return array.item()
    return array


def black_hole_radiative_efficiency(spin):
    r"""Calculate the Novikov--Thorne radiative efficiency of a Kerr BH.

    Parameters
    ----------
    spin : float or array-like
        Dimensionless Kerr spin in the closed interval ``[-1, 1]``. Invalid
        elements return ``nan``.

    Returns
    -------
    float or numpy.ndarray
        Dimensionless radiative efficiency.

    Notes
    -----
    Positive and negative spins represent prograde and retrograde discs,
    respectively. The calculation uses the specific energy at the innermost
    stable circular orbit (ISCO),

    .. math::

        \eta_{\rm BH} = 1 - E_{\rm ISCO}.

    The ISCO radius follows the analytic Kerr expression of Bardeen, Press &
    Teukolsky (1972). The thin-disc efficiency is the Novikov--Thorne value;
    photon capture and other departures from a radiatively efficient thin disc
    are not included.

    References
    ----------
    .. [1] Bardeen, J. M., Press, W. H., & Teukolsky, S. A. 1972, ApJ, 178,
       347. https://doi.org/10.1086/151796
    .. [2] Novikov, I. D., & Thorne, K. S. 1973, in *Black Holes*, eds. C.
       DeWitt & B. DeWitt, 343.
    .. [3] Bambi, C., Malafarina, D., & Tsukamoto, N. 2014, Phys. Rev. D,
       89, 127302. https://doi.org/10.1103/PhysRevD.89.127302
    """
    scalar = np.ndim(spin) == 0
    spin = np.asarray(spin, dtype=float)
    valid = np.isfinite(spin) & (spin >= -1.0) & (spin <= 1.0)
    safe_spin = np.where(valid, spin, 0.0)

    z1 = (1.0 + np.cbrt(1.0 - safe_spin**2)
          * (np.cbrt(1.0 + safe_spin) + np.cbrt(1.0 - safe_spin)))
    z2 = np.sqrt(3.0 * safe_spin**2 + z1**2)
    direction = np.where(safe_spin >= 0.0, 1.0, -1.0)
    r_isco = (3.0 + z2
              - direction * np.sqrt((3.0 - z1)
                                    * (3.0 + z1 + 2.0 * z2)))
    energy_isco = np.sqrt(1.0 - 2.0 / (3.0 * r_isco))
    efficiency = np.where(valid, 1.0 - energy_isco, np.nan)
    return _return_scalar_if_scalar(efficiency, scalar)


def neutron_star_radiative_efficiency(mass, radius):
    r"""Calculate a neutron star's Newtonian accretion efficiency.

    Parameters
    ----------
    mass : float or array-like
        Neutron-star mass in ``Msun``.
    radius : float or array-like
        Neutron-star radius in ``Rsun``.

    Returns
    -------
    float or numpy.ndarray
        Dimensionless efficiency ``G M / (R c^2)``. Invalid elements return
        ``nan``.

    Notes
    -----
    This is the Newtonian surface-accretion estimate

    .. math::

        \eta_{\rm NS} = \frac{G M_{\rm NS}}{R_{\rm NS}c^2}.

    It excludes relativistic surface corrections, magnetic truncation,
    propeller effects, and emission anisotropy. Those extensions require
    physical inputs that are not generally present in POSYDON populations.

    References
    ----------
    .. [1] Frank, J., King, A., & Raine, D. 2002, *Accretion Power in
       Astrophysics*, 3rd ed., Cambridge University Press.
       https://doi.org/10.1017/CBO9781139164245
    """
    scalar = np.ndim(mass) == 0 and np.ndim(radius) == 0
    mass, radius = np.broadcast_arrays(np.asarray(mass, dtype=float),
                                       np.asarray(radius, dtype=float))
    valid = (np.isfinite(mass) & np.isfinite(radius)
             & (mass > 0.0) & (radius > 0.0))
    efficiency = np.full(mass.shape, np.nan, dtype=float)
    np.divide(const.standard_cgrav * mass * const.Msun,
              radius * const.Rsun * const.clight**2,
              out=efficiency, where=valid)
    return _return_scalar_if_scalar(efficiency, scalar)


def eddington_limit_from_properties(accretor_mass, donor_surface_h1,
                                    efficiency):
    r"""Calculate the Eddington mass-accretion rate and luminosity.

    Parameters
    ----------
    accretor_mass : float or array-like
        Accretor mass in ``Msun``.
    donor_surface_h1 : float or array-like
        Donor surface hydrogen mass fraction.
    efficiency : float or array-like
        Dimensionless radiative efficiency of the accretor.

    Returns
    -------
    mdot_edd : float or numpy.ndarray
        Eddington accretion rate in ``Msun yr^-1``.
    luminosity_edd : float or numpy.ndarray
        Eddington luminosity in ``erg s^-1``.

    Notes
    -----
    For fully ionized material, the electron-scattering opacity is
    approximated by :math:`\kappa = 0.2(1+X)` in cgs units. The function uses

    .. math::

        L_{\rm Edd} = \frac{4\pi G M c}{0.2(1+X)}, \qquad
        \dot M_{\rm Edd} = \frac{L_{\rm Edd}}{\eta c^2}.

    Here :math:`X` is the donor surface hydrogen mass fraction and
    :math:`\eta` is supplied by the caller.

    References
    ----------
    .. [1] Rybicki, G. B., & Lightman, A. P. 1979, *Radiative Processes in
       Astrophysics*, Wiley.
    .. [2] Misra, D., et al. 2023, A&A, 672, A99.
       https://doi.org/10.1051/0004-6361/202244929
    """
    scalar = all(np.ndim(x) == 0 for x in
                 (accretor_mass, donor_surface_h1, efficiency))
    mass, surface_h1, efficiency = np.broadcast_arrays(
        np.asarray(accretor_mass, dtype=float),
        np.asarray(donor_surface_h1, dtype=float),
        np.asarray(efficiency, dtype=float),
    )
    valid = (np.isfinite(mass) & np.isfinite(surface_h1)
             & np.isfinite(efficiency) & (mass > 0.0)
             & (surface_h1 >= 0.0) & (surface_h1 <= 1.0)
             & (efficiency > 0.0))

    luminosity_edd = np.full(mass.shape, np.nan, dtype=float)
    luminosity_edd[valid] = (
        4.0 * const.pi * const.standard_cgrav * mass[valid] * const.Msun
        * const.clight / (0.2 * (1.0 + surface_h1[valid]))
    )
    mdot_edd = (luminosity_edd / (efficiency * const.clight**2)
                * const.secyer / const.Msun)
    return (_return_scalar_if_scalar(mdot_edd, scalar),
            _return_scalar_if_scalar(luminosity_edd, scalar))


# TODO(POSYDON#887): consolidate the wind prescription duplicated inside
# ``common_functions.bondi_hoyle`` after its behavior is reviewed in #888.
def wind_velocity(donor_mass, donor_radius, donor_luminosity,
                  donor_wind_mass_loss, donor_surface_h1,
                  donor_he_core_mass, scheme="Kudritzki+2000"):
    r"""Estimate the terminal wind velocity of a detached donor.

    Parameters
    ----------
    donor_mass, donor_radius, donor_luminosity : float or array-like
        Donor properties in ``Msun``, ``Rsun``, and ``Lsun``.
    donor_wind_mass_loss : float or array-like
        Positive wind mass-loss magnitude in ``Msun yr^-1``.
    donor_surface_h1 : float or array-like
        Surface hydrogen mass fraction.
    donor_he_core_mass : float or array-like
        Helium-core mass in ``Msun``.
    scheme : {'Hurley+2002', 'Kudritzki+2000'}
        Terminal-to-escape-speed prescription.

    Returns
    -------
    float or numpy.ndarray
        Wind velocity in ``cm s^-1``. Invalid elements return ``nan``.

    Notes
    -----
    Both named schemes scale the terminal wind speed to the surface escape
    speed. ``Hurley+2002`` uses the mass- and evolutionary-state factors of
    Hurley, Tout & Pols (2002). In the absence of a discrete BSE stellar-type
    flag, a positive helium-core mass together with
    :math:`R>900\,R_\odot` is used as a proxy for an extended H-rich giant;
    this approximation is tracked for review in POSYDON issue #888.

    ``Kudritzki+2000`` adopts temperature-dependent ratios of terminal to
    escape velocity. Hot stripped donors with surface hydrogen fraction below
    0.4 use the empirical mass-loss-rate relation of Sander & Vink (2020).

    References
    ----------
    .. [1] Hurley, J. R., Tout, C. A., & Pols, O. R. 2002, MNRAS, 329, 897.
       https://doi.org/10.1046/j.1365-8711.2002.05038.x
    .. [2] Kudritzki, R.-P., & Puls, J. 2000, ARA&A, 38, 613.
       https://doi.org/10.1146/annurev.astro.38.1.613
    .. [3] Sander, A. A. C., & Vink, J. S. 2020, MNRAS, 499, 873.
       https://doi.org/10.1093/mnras/staa2712
    """
    if scheme not in {"Hurley+2002", "Kudritzki+2000"}:
        raise ValueError("scheme must be 'Hurley+2002' or 'Kudritzki+2000'")

    inputs = (donor_mass, donor_radius, donor_luminosity,
              donor_wind_mass_loss, donor_surface_h1, donor_he_core_mass)
    scalar = all(np.ndim(value) == 0 for value in inputs)
    mass, radius, luminosity, mdot_wind, surface_h1, he_core_mass = (
        np.broadcast_arrays(*[np.asarray(value, dtype=float)
                              for value in inputs])
    )
    valid = (np.isfinite(mass) & np.isfinite(radius)
             & np.isfinite(luminosity) & np.isfinite(mdot_wind)
             & np.isfinite(surface_h1) & np.isfinite(he_core_mass)
             & (mass > 0.0) & (radius > 0.0) & (luminosity > 0.0)
             & (mdot_wind > 0.0) & (surface_h1 >= 0.0)
             & (surface_h1 <= 1.0) & (he_core_mass >= 0.0))

    effective_temperature = np.full(mass.shape, np.nan, dtype=float)
    effective_temperature[valid] = (
        luminosity[valid] * const.Lsun
        / (4.0 * const.pi * (radius[valid] * const.Rsun)**2
           * const.boltz_sigma)
    )**0.25

    if scheme == "Kudritzki+2000":
        factor = np.where(effective_temperature >= 21000.0, 2.65,
                          np.where(effective_temperature <= 10000.0,
                                   1.0, 1.4))
    else:
        beta_h = np.where(
            (he_core_mass > 0.0) & (radius > 900.0), 0.125,
            np.where(mass > 120.0, 7.0,
                     np.where(mass < 1.4, 0.5,
                              0.5 + (mass - 1.4) * 6.5 / 118.6)),
        )
        beta_he = np.where(
            mass > 120.0, 7.0,
            np.where(mass < 10.0, 0.125,
                     0.125 + (mass - 10.0) * 6.875 / 110.0),
        )
        factor = np.sqrt(np.where(surface_h1 <= 0.01, beta_he, beta_h))

    velocity = np.full(mass.shape, np.nan, dtype=float)
    velocity[valid] = (factor[valid]
                       * np.sqrt(2.0 * const.standard_cgrav
                                 * mass[valid] * const.Msun
                                 / (radius[valid] * const.Rsun)))

    stripped_hot = valid & (surface_h1 < 0.4) & (effective_temperature > 1e4)
    log_mdot = np.full(mass.shape, np.nan, dtype=float)
    log_mdot[valid] = np.log10(mdot_wind[valid])
    high_rate_slope = (3.7 - 3.25) / (-2.5 + 5.25)
    low_rate_slope = (3.25 - 3.75) / (-5.25 + 7.25)
    slope = np.where(log_mdot >= -5.25, high_rate_slope, low_rate_slope)
    velocity[stripped_hot] = (
        10.0**(slope[stripped_hot] * log_mdot[stripped_hot]
               + 3.25 + 5.25 * slope[stripped_hot]) * const.km2cm
    )
    return _return_scalar_if_scalar(velocity, scalar)


# TODO(POSYDON#887): consolidate this calculation with the legacy
# ``common_functions.bondi_hoyle`` adapter after resolving #829 and #888.
def bondi_hoyle_accretion_rate(accretor_mass, donor_mass,
                               donor_wind_mass_loss, orbital_separation,
                               eccentricity, wind_speed, alpha=1.5):
    r"""Calculate an orbit-averaged Bondi--Hoyle--Lyttleton accretion rate.

    The analytic Boffin--Jorissen/Hurley orbit average assumes a constant wind
    speed and an orbital period short compared with changes in the stellar
    properties. It replaces random orbital-phase sampling and is deterministic.

    Parameters
    ----------
    accretor_mass : float or array-like
        Accretor mass in ``Msun``.
    donor_mass : float or array-like
        Donor mass in ``Msun``.
    donor_wind_mass_loss : float or array-like
        Positive wind mass-loss magnitude in ``Msun yr^-1``.
    orbital_separation : float or array-like
        Semi-major axis in ``Rsun``.
    eccentricity : float or array-like
        Orbital eccentricity in the half-open interval ``[0, 1)``.
    wind_speed : float or array-like
        Wind velocity in ``cm s^-1``.
    alpha : float, optional
        Dimensionless BHL efficiency. Must be positive.

    Returns
    -------
    float or numpy.ndarray
        Orbit-averaged accretion rate in ``Msun yr^-1``. Invalid physical
        input elements return ``nan``.

    Notes
    -----
    The implemented prescription is

    .. math::

        \left\langle\dot M_{\rm BHL}\right\rangle =
        \frac{\alpha}{2\sqrt{1-e^2}}
        \left(\frac{G M_{\rm acc}}{a v_{\rm wind}^2}\right)^2
        \left[1 + \frac{G(M_{\rm acc}+M_{\rm don})}
        {a v_{\rm wind}^2}\right]^{-3/2}\dot M_{\rm wind}.

    It assumes an isotropic wind with constant speed and stellar properties
    over an orbit. The orbital-velocity correction is evaluated at the
    semi-major axis; this is not a direct hydrodynamical phase integration.
    The relation to the legacy random-phase implementation remains a physics
    question under discussion in POSYDON issue #829.

    References
    ----------
    .. [1] Bondi, H., & Hoyle, F. 1944, MNRAS, 104, 273.
       https://doi.org/10.1093/mnras/104.5.273
    .. [2] Boffin, H. M. J., & Jorissen, A. 1988, A&A, 205, 155.
    .. [3] Hurley, J. R., Tout, C. A., & Pols, O. R. 2002, MNRAS, 329, 897.
       https://doi.org/10.1046/j.1365-8711.2002.05038.x
    """
    if not np.isscalar(alpha) or not np.isfinite(alpha) or alpha <= 0.0:
        raise ValueError("alpha must be a positive finite scalar")

    inputs = (accretor_mass, donor_mass, donor_wind_mass_loss,
              orbital_separation, eccentricity, wind_speed)
    scalar = all(np.ndim(value) == 0 for value in inputs)
    mass, donor_mass, mdot_wind, separation, eccentricity, speed = (
        np.broadcast_arrays(*[np.asarray(value, dtype=float)
                              for value in inputs])
    )
    valid = (np.isfinite(mass) & np.isfinite(donor_mass)
             & np.isfinite(mdot_wind)
             & np.isfinite(separation) & np.isfinite(eccentricity)
             & np.isfinite(speed) & (mass > 0.0) & (donor_mass > 0.0)
             & (mdot_wind >= 0.0)
             & (separation > 0.0) & (eccentricity >= 0.0)
             & (eccentricity < 1.0) & (speed > 0.0))

    rate = np.full(mass.shape, np.nan, dtype=float)
    mass_cgs = mass * const.Msun
    separation_cgs = separation * const.Rsun
    orbital_speed_squared = (const.standard_cgrav
                             * (mass + donor_mass) * const.Msun
                             / separation_cgs)
    velocity_correction = (1.0 + orbital_speed_squared / speed**2)**1.5
    denominator = (2.0 * separation_cgs**2 * speed**4
                   * np.sqrt(1.0 - eccentricity**2)
                   * velocity_correction)
    fraction = np.full(mass.shape, np.nan, dtype=float)
    np.divide(alpha * (const.standard_cgrav * mass_cgs)**2,
              denominator, out=fraction, where=valid)
    rate[valid] = fraction[valid] * mdot_wind[valid]
    return _return_scalar_if_scalar(rate, scalar)


def accretion_luminosity(accretion_rate, accretor_mass, donor_surface_h1,
                         efficiency, beaming_floor=3.2e-3):
    r"""Calculate classical sub- and super-Eddington luminosities.

    The intrinsic luminosity is ``eta * mdot * c^2`` below the Eddington
    limit. Above it, the Shakura--Sunyaev logarithmic luminosity enhancement
    and King beaming prescription are used. The returned luminosity is the
    isotropic-equivalent bolometric luminosity; no detector bandpass or
    bolometric correction is applied.

    Parameters
    ----------
    accretion_rate : float or array-like
        Accretion rate in ``Msun yr^-1``.
    accretor_mass : float or array-like
        Accretor mass in ``Msun``.
    donor_surface_h1 : float or array-like
        Donor surface hydrogen mass fraction.
    efficiency : float or array-like
        Dimensionless accretor radiative efficiency.
    beaming_floor : float, optional
        Minimum King beaming factor. Must lie in ``(0, 1]``.

    Returns
    -------
    luminosity, beaming_factor, eddington_ratio : tuple
        Scalars or arrays containing ``erg s^-1``, a dimensionless beaming
        factor, and ``mdot / mdot_edd``.

    Notes
    -----
    For :math:`\dot m=\dot M/\dot M_{\rm Edd}\leq 1`, the luminosity is
    :math:`L=\eta\dot M c^2`. Above the Eddington rate, the classical
    Shakura--Sunyaev logarithmic enhancement is used:

    .. math::

        L_{\rm iso} = \frac{L_{\rm Edd}}{b}
        \left(1+\ln\dot m\right).

    The King beaming factor is unity through :math:`\dot m=8.5` and
    :math:`b=73/\dot m^2` above that transition, bounded below by
    ``beaming_floor``. The result is an isotropic-equivalent *bolometric*
    luminosity. Detector bandpasses and bolometric corrections belong in the
    calling analysis.

    References
    ----------
    .. [1] Shakura, N. I., & Sunyaev, R. A. 1973, A&A, 24, 337.
    .. [2] King, A. R. 2008, MNRAS, 385, L113.
       https://doi.org/10.1111/j.1745-3933.2008.00444.x
    .. [3] Misra, D., et al. 2023, A&A, 672, A99.
       https://doi.org/10.1051/0004-6361/202244929
    """
    if (not np.isscalar(beaming_floor) or not np.isfinite(beaming_floor)
            or not 0.0 < beaming_floor <= 1.0):
        raise ValueError("beaming_floor must be a finite scalar in (0, 1]")

    inputs = (accretion_rate, accretor_mass, donor_surface_h1, efficiency)
    scalar = all(np.ndim(value) == 0 for value in inputs)
    mdot, mass, surface_h1, efficiency = np.broadcast_arrays(
        *[np.asarray(value, dtype=float) for value in inputs]
    )
    mdot_edd, luminosity_edd = eddington_limit_from_properties(
        mass, surface_h1, efficiency)
    mdot_edd = np.asarray(mdot_edd, dtype=float)
    luminosity_edd = np.asarray(luminosity_edd, dtype=float)
    valid = (np.isfinite(mdot) & (mdot >= 0.0)
             & np.isfinite(mdot_edd) & (mdot_edd > 0.0))

    Eddington_ratio = np.full(mdot.shape, np.nan, dtype=float)
    Eddington_ratio[valid] = mdot[valid] / mdot_edd[valid]
    beaming = np.full(mdot.shape, np.nan, dtype=float)
    beaming[valid] = 1.0
    strongly_super_eddington = valid & (Eddington_ratio > 8.5)
    beaming[strongly_super_eddington] = np.clip(
        73.0 / Eddington_ratio[strongly_super_eddington]**2,
        beaming_floor, 1.0)

    luminosity = np.full(mdot.shape, np.nan, dtype=float)
    sub_eddington = valid & (Eddington_ratio <= 1.0)
    mdot_cgs = mdot * const.Msun / const.secyer
    luminosity[sub_eddington] = (
        efficiency[sub_eddington] * mdot_cgs[sub_eddington]
        * const.clight**2
    )
    super_eddington = valid & (Eddington_ratio > 1.0)
    luminosity[super_eddington] = (
        luminosity_edd[super_eddington]
        * (1.0 + np.log(Eddington_ratio[super_eddington]))
        / beaming[super_eddington]
    )
    return tuple(_return_scalar_if_scalar(value, scalar)
                 for value in (luminosity, beaming, Eddington_ratio))


def be_xray_luminosity(orbital_period):
    r"""Calculate the empirical Be-XRB period--luminosity relation.

    Parameters
    ----------
    orbital_period : float or array-like
        Orbital period in days. The caller is responsible for identifying Be
        systems; non-positive or non-finite periods return ``nan``.

    Returns
    -------
    float or numpy.ndarray
        X-ray luminosity in ``erg s^-1``.

    Notes
    -----
    The empirical relation used in POSYDON XLF work is

    .. math::

        L_X = 10^{35}\,10^{4.53-1.5\log_{10}(P_{\rm orb}/{\rm day})}
        \;{\rm erg\,s^{-1}}.

    The function applies only the period--luminosity relation. Identifying Be
    donors, imposing an applicable period interval, and applying a duty cycle
    remain responsibilities of the calling population analysis.

    References
    ----------
    .. [1] Dai, H.-L., Liu, X.-W., & Li, X.-D. 2006, ApJ, 653, 1410.
       https://doi.org/10.1086/508735
    .. [2] Misra, D., et al. 2023, A&A, 672, A99.
       https://doi.org/10.1051/0004-6361/202244929
    """
    scalar = np.ndim(orbital_period) == 0
    period = np.asarray(orbital_period, dtype=float)
    valid = np.isfinite(period) & (period > 0.0)
    luminosity = np.full(period.shape, np.nan, dtype=float)
    luminosity[valid] = (1.0e35
                         * 10.0**(4.53 - 1.5 * np.log10(period[valid])))
    return _return_scalar_if_scalar(luminosity, scalar)
