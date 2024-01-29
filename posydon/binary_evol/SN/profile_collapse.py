"""Collapse the profile of a star object into a BH.

This script is based on the physics explained in Appendix D of Bavera+2020.

"""

import numpy as np
from scipy import integrate
import posydon.utils.constants as const

from posydon.utils.gridutils import find_index_nearest_neighbour

__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Scott Coughlin <scottcoughlin2014@u.northwestern.edu>",
    "Devina Misra <devina.misra@unige.ch>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]


__credits__ = [
    'Aldo Batta <aldobatta@gmail.com>',
]


def get_initial_BH_properties(star, mass_collapsing, mass_central_BH,
                              neutrino_mass_loss, max_neutrino_mass_loss,
                              verbose):
    """Collapse directly the center of the star and return useful quantities.

    Parameters
    ----------
    star : object
        Star object of a collapsing star containing the MESA profile.
    mass_collapsing : float
        Remnant barionic mass in M_sun collapsing to form the BH. This is the
        mass left to collapse after applying a supernova prescriptions,
        see e.g., rapid and delayed mechanisms of Fryer et al. (2012).
    mass_central_BH : float
        Mass of the central stellar layers (in M_sun) collasping directly to
        form a proto BH.
    neutrino_mass_loss : float
        Mass (in M_sun) lost through neutrinos in the formation of the
        central BH.
    max_neutrino_mass_loss : float
        Maximum mass (in M_sun) lost thorugh neutrinos.
    verbose : bool
        If `True`, it prints some informations.

    Returns
    -------
    mass_initial_BH : float
        Mass of the initial BH in units of g.
    a_initial_BH : float
        Dimensionless spin of the initial BH.
    J_initial_BH : float
        Angular momentum of the initial BH in g*cm^2/s.
    angular_frequency_i : array floats
        Shell's angular frequencies in s^-1 collapsing onto the
        initially-formed BH.
    enclosed_mass_i : array floats
        Shell's enclosed masses in g collapsing onto the initially formed BH.
    radius_i : array floats
        Shell's radii in cm collapsing onto the initially formed BH.
    density_i : array floats
        Shell's densities in g/cm^3 collapsing onto the initially formed BH.
    dm_i : array floats
        Shell's masses in g collapsing onto the initially formed BH.
    dm_i : array floats
        Shell's width in cm collapsing onto the initially formed BH.

    """
    if neutrino_mass_loss < 0.:
        raise ValueError(
            'Something went wrong, neutrino_mass_loss must be positive!')

    if neutrino_mass_loss > max_neutrino_mass_loss:
        raise ValueError(
            'Something went wrong, ',
            'max_neutrino_mass_loss = {:2.2f} Msun'.format(
                max_neutrino_mass_loss),
            'while neutrino_mass_loss = {:2.2f} Msun'.format(
                neutrino_mass_loss), 'was passed!')

    # load units and convert to CGS units
    G = const.standard_cgrav
    c = const.clight
    Mo = const.Msun
    Ro = const.Rsun
    mass_central_BH *= Mo
    neutrino_mass_loss *= Mo

    # read star quantities
    enclosed_mass_all = star.profile['mass'][::-1]*Mo  # cell outer total mass
    radius_all = star.profile['radius'][::-1] * Ro  # cell outer radius
    density_all = 10**(star.profile['logRho'][::-1])  # cell density
    angular_frequency_all = star.profile['omega'][::-1]  # cell angular freq.
    if 'he4' in star.profile.dtype.names:
        # he3_all = star.profile['he3'][::-1]  # he4 mass fraction
        he4_all = star.profile['he4'][::-1]  # he4 mass fraction
        he3_all = np.zeros(len(enclosed_mass_all))  # ENHANCEMENT support he3
    else:
        he4_all = np.zeros(len(enclosed_mass_all))
        he3_all = np.zeros(len(enclosed_mass_all))

    # cut the ejected layers of the profile due to the SN event:
    # index of the layers up to m_rembar from Fryer prescription
    # (+1 for the range)
    if enclosed_mass_all[-1] / Mo <= mass_collapsing:
        # This catches the case that all the star's profile is callapsed.
        # Note that the 'mass' of the MESA profile is the enclosed mass of that
        # shell; the mass of the whole star is then
        #     star.profile['mass'][::-1][-1] + dm,
        # where dm is the mass of the last shell.
        i_rem = len(enclosed_mass_all)
    else:
        i_rem = np.argmax(enclosed_mass_all/Mo > mass_collapsing) + 1

    enclosed_mass = enclosed_mass_all[:i_rem]
    radius = radius_all[:i_rem]
    density = density_all[:i_rem]
    angular_frequency = angular_frequency_all[:i_rem]
    he3 = he3_all[:i_rem]
    he4 = he4_all[:i_rem]

    # max he mass ejected in the SN
    dm_SN = enclosed_mass_all[i_rem+1:] - enclosed_mass_all[i_rem:-1]
    max_he_mass_ejected_SN = sum((he3_all[i_rem:-1] + he4_all[i_rem:-1])*dm_SN)

    # find index containing the mass MBH_0 (in CGS)
    index_initial_BH = find_index_nearest_neighbour(enclosed_mass,
                                                    mass_central_BH)
    # mass of the initial BH collapsing directy assuming that
    # neutrino_mass_loss is lost thorugh neutrinos in the formation of
    # the central BH, note: this neutrinos are carring away angular
    # momentum proportional to the neutrino_mass_loss/mass_central_BH
    mass_initial_BH = enclosed_mass[index_initial_BH] - neutrino_mass_loss

    # Cut input arrays to consider only the shells AFTER the formation of the
    # initial BH note add 1 to get the range after the value index_initial_BH
    angular_frequency_i = angular_frequency[index_initial_BH + 1:]
    enclosed_mass_i = enclosed_mass[index_initial_BH + 1:]
    radius_i = radius[index_initial_BH + 1:]
    density_i = density[index_initial_BH + 1:]
    he3_i = he3[index_initial_BH + 1:]
    he4_i = he4[index_initial_BH + 1:]

    # shell's mass
    dm = enclosed_mass[1:] - enclosed_mass[:-1]
    # shell's width
    dr = radius[1:] - radius[:-1]

    # dm has len(enclosed_mass_initial_BH) = len(dM_initial_BH)-1
    dm_i = dm[index_initial_BH:]
    dr_i = dr[index_initial_BH:]

    # Compute the angular momentum of the initial BH by solving eq. 1 of
    # Batta & Ramirez 2019:
    # J_i_BH = int Omega(r) r^2 sin^2(t) dm dr
    #        = iint 2pi Omega(r) r^4 sin^3(t) rho(r) dt dr
    #        = 2pi int_0^pi sin^3(t) dt int_0^M_core Omega(r) r^4 rho(r) dr
    #        =: 2pi * temp1 * temp2
    def f_temp1(t):
        return np.sin(t)**3

    temp1 = integrate.quad(f_temp1, 0, np.pi)[0]

    # f_nu_AM is a rescaling the J_initial_BH. This account for the
    # angular momentum lost thorugh neutrinos, which under the assumption
    # of efficient AM transport is really low.
    f_nu_AM = mass_initial_BH / enclosed_mass[index_initial_BH]

    f_temp2 = (f_nu_AM*density[:index_initial_BH+1]
               * angular_frequency[:index_initial_BH+1]
               * radius[:index_initial_BH+1]**4)
    temp2 = integrate.simps(f_temp2, x=radius[:index_initial_BH + 1])

    J_initial_BH = 2 * np.pi * temp1 * temp2

    # dimensless spinn of the initial BH
    a_initial_BH = J_initial_BH * c / (G * mass_initial_BH**2)

    if verbose:
        print('')
        print('Initializing the BH properties for 1D collapse.')
        print('')
        print('The BH formed from the ',
              round(mass_central_BH / Mo,
                    2), 'innermost solar masses has mass',
              round(mass_initial_BH / Mo,
                    2), 'M_sun (', neutrino_mass_loss / Mo,
              'M_sun were lost thorugh neutrinos) and a = Jc/GM^2 =',
              round(a_initial_BH, 4))

    return [
        mass_initial_BH, a_initial_BH, J_initial_BH, angular_frequency_i,
        enclosed_mass_i, radius_i, density_i, dm_i, dr_i,
        he3_i, he4_i, max_he_mass_ejected_SN
    ]


def compute_isco_properties(a, m_BH):
    """Compute the BH innermost stable circular orbit (ISCO) parameters.

    Parameters
    ----------
    a : float
        Dimnesionless BH spin.
    m_BH : float
        Mass of the BH in g.

    Returns
    -------
    r_isco : float
        Radius of the ISCO in CGS units (cm).
    j_isco : float
        Specific angular momentum at ISCO in CGS units (cm^2/s).
    efficiency : float
        Orbital energy efficiency at ISCO.

    """
    # load units
    G = const.standard_cgrav
    c = const.clight

    # eq. 3/4 in Batta & Ramirez-Ruiz (2019)
    z1 = 1 + (((1 - a**2)**(1. / 3.))
              * ((1 + a)**(1. / 3.) + (1 - a)**(1. / 3.)))
    z2 = (3 * a**2 + z1**2)**(1. / 2.)
    r_isco = (3 + z2 - np.sqrt((3 - z1) * (3 + z1 + 2 * z2)))

    # note that eq. 5 is the equivalent of what is below expet for the fact
    # that what is below is not defined for a=1
    j_isco = (
        (G * m_BH / c) * (r_isco**2 - 2 * a * r_isco**(0.5) + a**2)
        / (r_isco**(3./4.) * (r_isco**(3./2.) - 3 * r_isco**0.5 + 2 * a)**0.5)
    )
    # fraction of disk's mass accreted by BH
    # (1-efficiency) of the rest mass is "radiated" away
    efficiency = np.sqrt(1. - 2. / (3 * r_isco))

    # assign CGS units
    r_isco = (G * m_BH / c**2) * r_isco

    return r_isco, j_isco, efficiency


def do_core_collapse_BH(star,
                        mass_collapsing,
                        mass_central_BH=2.51,
                        neutrino_mass_loss=None,
                        max_neutrino_mass_loss=0.5,
                        verbose=False):
    """Do the core collapse of a star object with MESA profile provided.

    Parameters
    ----------
    star : object
        Star object of a collapsing star containing the MESA profile.
    mass_collapsing : float
        Remnant barionic mass in M_sun collapsing to form the BH. This is the
        mass left to collapse after applying a supernova prescriptions,
        see e.g. rapid and delayed mechanisms of Fryer et al. (2012).
    mass_central_BH : float
        Mass of the central stellar layers (in M_sun) collasping directly to
        form a proto BH.
    neutrino_mass_loss : float
        Mass (in M_sun) lost thorugh neutrinos in the formation of the central
        BH.
    max_neutrino_mass_loss : float
        Maximum mass (in M_sun) lost thorugh neutrinos.
    verbose : bool
        If `True`, it prints some informations.

    Returns
    -------
    M_BH_total : float
        Mass of the final BH in M_sun.
    a_BH_total : float
        Dimensionless spin of the final BH.
    M_BH_array : array floats
        BH mass evelution in g.
    a_BH_array : array floats
        Dimensionless spin evolution.

    J_accreted_array : array floats
        Angular momentum accreted from a given shell by the BH in CGS units.
    J_total_array : array floats
        Total angular momentum in accreted shells plus BH's initial angular
        momentum in CGS units.
    J_disk_shell_array : array floats
        Angular momentum accreted from the shell's part collapsing to form a
        disk in CGS units.
    radiation_eff_array : array floats
        Fraction of accretion disk radiated away, this is one minus accretion
        efficiency.
    r_isco_array : array floats
        Radius of the innermost stable circular orbit in cm.
    j_isco_array : array floats
        Specific angular momentum at the innermost stable circular orbit in
        CGS.
    M_direct_collapse_array : array floats
        Cumulative mass accreted through direct collapse in g.
    M_disk_array : array floats
        Cumulative mass accreted thorugh the disk in g.
    dm_direct_array : array floats
        Shell's mass accreted through direct collapse in g.
    dm_disk_array : array floats
        Shell's mass accreted thorugh the disk in g.
    j_shell_array : array floats
        Shell's specific angular momentum in CGS.
    M_total_array : array floats
        Cumulative mass of shells and initial BH in g.
    a_star_array : array floats
        Dimensionless spin parameter of the star.

    """
    # convert to CGS units
    G = const.standard_cgrav
    c = const.clight
    Mo = const.Msun

    # ===========================================================
    # Assumes a BH is already formed with a certain spin and mass
    # its spin parameter is a = Jc/GM^2 and could be > 1
    # ===========================================================

    # Extract results from Get_InitialBHProperties()
    (M_BH, a_BH, J_BH, Omega_shells, enclosed_mass, radius_shells, rho_shells,
     dm, dr, he3, he4, max_he_mass_ejected) = get_initial_BH_properties(
         star, mass_collapsing, mass_central_BH, neutrino_mass_loss,
         max_neutrino_mass_loss, verbose)

    # check that there is matter falling onto the BH
    if len(enclosed_mass) == 0:
        arr = np.array([np.nan])
        return [
            M_BH / Mo, a_BH, np.nan, np.nan
            #M_BH / Mo, a_BH, arr, arr, arr, arr, arr, arr,
            #arr, arr, arr, np.array([0.]), arr, arr, arr, arr,
            #arr, arr
        ]

    # shell's specific angular momentum at equator
    j_shells = Omega_shells * radius_shells**2

    if a_BH > 0.99:
        a_BH = 0.99
        J_BH = a_BH * G * M_BH**2 / c

    # get the initial r_isco, j_isco and orbital energy convertion efficiency
    r_isco, j_isco, eff = compute_isco_properties(a_BH, M_BH)

    # Initialize integrated quantities
    M_direct_collapse = M_BH
    M_disk = 0.
    M_total = M_BH
    dm_disk = 0.
    dm_direct = 0.
    J_total = J_BH

    # =======================================================
    # Initialize lists that will contain the BH's properties
    # as a function of the collapsed shells and information
    # of mass fraction of material with low j and high j

    J_accreted_array = [J_BH]  # angular momentum accreted by BH
    J_total_array = [
        J_total
    ]  # total angular momentum in accreted shells + BH's initial J
    M_BH_array = [M_BH]  # BH mass
    a_BH_array = [a_BH]  # BH's spin
    radiation_eff_array = [1 - eff]  # fraction of accretion disk radiated away
    r_isco_array = [r_isco]  # radius of ISCO
    j_isco_array = [j_isco
                    ]  # specific angular momentum at ISCO (prograde orbits)
    M_direct_collapse_array = [
        M_direct_collapse
    ]  # integrated mass accreted through direct collapse
    M_disk_array = [M_disk]  # integrated mass accreted through disk
    dm_disk_array = [dm_disk]  # mass in shell with j > j_isco (forms a disk)
    dm_direct_array = [dm_direct
                       ]  # mass in shell with j < j_isco (direct collapse)
    M_total_array = [M_total]  # integrated mass (shells + initial BH)
    j_shell_array = [j_shells[0]]  # shell's specific angular momentum
    a_star_array = [a_BH]  # star's spin parameter
    J_disk_shell_array = [0.]  # angular momentum of the shell's disk

    # compute BH properties as each shell collapses
    for i, value in enumerate(dm):

        # get r_isco, j_isco, and orbital energy at isco
        r_isco, j_isco, eff = compute_isco_properties(a_BH, M_BH)

        # shell properties
        j_shell = j_shells[i]
        Omega_shell = Omega_shells[i]
        r_shell = radius_shells[i]
        dr_shell = dr[i]
        dm_shell = dm[i]
        rho_shell = rho_shells[i]
        he3_shell = he3[i]
        he4_shell = he4[i]

        # Determine if the specific angular momentum of the shell can form
        # a disk or not and update BH's properties

        # All mass collapses directly to the BH
        if j_shell < j_isco:

            # eq. 9 of Batta & Ramirez-Ruiz (2019) for theta<theta_disk
            # J_shell = 2 int Omega(r) sin^3(t) rho r^4 dr dt dphi
            #   = 4 pi Omega_r rho_r int_0^pi/2 sin^3(t) dt int_r-dr^r r^4 dr
            #   = 4 pi Omega_r rho_r temp1 temp2

            def f_temp1(x):
                return np.sin(x)**3

            temp1 = integrate.quad(f_temp1, 0, np.pi / 2)[0]

            def f_temp2(x):
                return x**4

            temp2 = integrate.quad(f_temp2, r_shell - dr_shell, r_shell)[0]

            # angular momentum of entire shell: J_shell=J_direct
            J_direct = 4 * np.pi * (Omega_shell * rho_shell) * temp1 * temp2

            # Update BH's angular momentum content and mass
            J_BH = J_BH + J_direct
            M_BH = M_BH + dm_shell

            # mass accreted through direct collapse
            dm_direct = dm_shell
            # mass forming a disk
            dm_disk = 0.0

            # Integrated mass from direct collapse
            M_direct_collapse = M_direct_collapse + dm_direct

            # angular momentum accreted from the disk
            J_disk = 0.  # There is no disk

        # Portions of the shells with theta>theta_disk are forming a disk
        else:

            # theta_disk is the angle within which
            # j_shell = Omega_shell * r^2 * sin(t))^2 < j_isco
            # contains all material that will collapse directly to the BH
            # eq. 7 of  Batta & Ramirez-Ruiz (2019)
            theta_disk = np.arcsin(np.sqrt(j_isco / (Omega_shell*r_shell**2)))

            # eq. 9 of Batta & Ramirez-Ruiz (2019) for theta<theta_disk
            # J_shell = J_direct + J_disk

            def f_temp1(x):
                return np.sin(x)**3

            temp1 = integrate.quad(f_temp1, 0, theta_disk)[0]

            def f_temp2(x):
                return x**4

            temp2 = integrate.quad(f_temp2, r_shell - dr_shell, r_shell)[0]

            J_direct = 4 * np.pi * (Omega_shell * rho_shell) * temp1 * temp2

            # shell's mass fraction that will collapse directly
            # eq. 8 of Batta & Ramirez-Ruiz (2019)
            # (typo corrected, see Bavera et al. (2020))
            dm_direct = (1.0 - np.cos(theta_disk)) * dm_shell
            # shell's mass fraction that will form an accretion disk
            # note that we assume that only mass collapsing directy loses
            # neutrinos
            dm_disk = np.cos(theta_disk) * dm_shell

            # angular momentum accreted thorugh the disk
            # eq. 8 of Batta & Ramirez-Ruiz (2019)
            J_disk = j_isco * dm_disk

            # Update BH's angular momentum content and mass
            J_BH = J_BH + J_disk + J_direct
            M_BH = M_BH + dm_disk * eff + dm_direct

            # Integrated mass from direct collapse
            M_direct_collapse = M_direct_collapse + dm_direct
            # Integrated mass from accretion disk
            # eq. 8 of Batta & Ramirez-Ruiz (2019): only fraction eff is
            # accreted, see Thorne (1974)
            M_disk = M_disk + dm_disk * eff

            # max He mass tht can be ejected during the disk formation
            max_he_mass_ejected += dm_shell * (he3_shell + he4_shell)

        # Update BH's spin parameter
        a_BH = J_BH * c / (G * M_BH**2.)

        # Shell's angular momentum (same as J_direct if j < j_isco ) assuming
        J_shell = 8 * np.pi / 3. * temp2 * Omega_shell * rho_shell
        M_total = M_total + dm_direct + dm_disk
        J_total = J_total + J_shell
        a_star = J_total * c / (G * M_total**2)

        # Check if a > 1: this should not happen!
        if a_BH > 1:
            raise ValueError(
                "We got a={:.5g} from shell {} containing {:.5g} M_sun".format(
                    a_BH, i, dm_shell / Mo))

        # Append all quantities to the arrays
        J_accreted_array.append(J_BH)
        M_BH_array.append(M_BH)
        a_BH_array.append(a_BH)
        radiation_eff_array.append(1. - eff)
        r_isco_array.append(r_isco)
        j_isco_array.append(j_isco)
        M_direct_collapse_array.append(M_direct_collapse)
        M_disk_array.append(M_disk)
        dm_disk_array.append(dm_disk)
        dm_direct_array.append(dm_direct)
        j_shell_array.append(j_shell)
        M_total_array.append(M_total)
        J_total_array.append(J_total)
        a_star_array.append(a_star)
        if j_shell < 1.00 * j_isco:
            J_disk_shell_array.append(0.)
        else:
            J_disk_shell_array.append((J_shell - J_direct) / dm_disk)

    # BH mass from the collapse of the entire star in Msun
    M_BH_total = M_BH_array[-1] / Mo
    # BH spin from the collapse of the entire star
    a_BH_total = a_BH_array[-1]
    # BH disk mass accreted in Msun
    m_disk_accreted = M_disk_array[-1] / Mo
    # BH disk mass radiate in Msun
    m_disk_radiated = sum(np.array(dm_disk_array)*np.array(radiation_eff_array))/Mo
    # max He mass tht can be ejected during the disk formation
    max_he_mass_ejected /= Mo

    if verbose:
        print('')
        print('Evolving BH properties in 1D collapse')
        print('')
        print('The BH formed from the collapse of the entire star with mass',
              round(enclosed_mass[-1] / Mo, 2), 'M_sun has mass',
              round(M_BH_total, 2), 'M_sun', 'and a = ', round(a_BH_total, 3),
              '.')
        print('A disk of total mass', round(M_disk_array[-1] / Mo, 2), 'M_sun',
              'was formed around the BH.')
        print('')

    return [
        M_BH_total,
        a_BH_total,
        m_disk_accreted, 
        m_disk_radiated,
        # np.array(M_BH_array),
        # np.array(a_BH_array),
        # np.array(J_accreted_array),
        # np.array(J_total_array),
        # np.array(J_disk_shell_array),
        # np.array(radiation_eff_array),
        # np.array(r_isco_array),
        # np.array(j_isco_array),
        # np.array(M_direct_collapse_array),
        # np.array(M_disk_array),
        # np.array(dm_direct_array),
        # np.array(dm_disk_array),
        # np.array(j_shell_array),
        # np.array(M_total_array),
        # np.array(a_star_array),
        # max_he_mass_ejected,
    ]
