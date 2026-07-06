__version__ = "1.4.0"

__author__ = "Arnaud Aguet <arnaud.aguet@etu.unige.ch><arnaud552000@gmail.com>"


import os

from scipy.optimize import newton                                                            # type : ignore
import matplotlib.pyplot as plt                                                              # type : ignore
import numpy as np
import pandas as pd                                                                          # type : ignore
import warnings                                                                              # type : ignore

import posydon                                                                               # type : ignore
import posydon.utils.constants as const                                                      # type : ignore 
from posydon.popsyn.synthetic_population import Population
from posydon.utils.common_functions import orbital_separation_from_period, roche_lobe_radius # type : ignore
from posydon.utils.posydonwarning import Pwarn, POSYDONWarning                               # type : ignore

# Convert POSYDON warnings to exceptions to stop execution on warnings
warnings.filterwarnings('error', category=POSYDONWarning)

#===== Functions for popualtion selection and extraction =====
def Systems_selection_for_XRB(path_file_h5, filename_out, verbose=True):
    """
    Selects compact-object binaries with non-compact companions and exports the new population.

    Keeps systems that have one compact object (NS/BH) and a companion that is
    neither a compact object nor a massless remnant, and that are not in excluded system states.

    The output file(s) can be later used in create_transient_population() with X-ray binary selection function. 

    Parameters
    ----------
    path_file_h5 : str or list of str
        Path(s) to the input .h5 file(s) containing binary systems.

    filename_out : str or list of str
        Name(s) for the file(s) that will be exported (without the extension).

    verbose : bool, optional
        If True, print detailed information about the process. Default is True.

    Returns
    -------
    None
        The function exports selected systems to .h5 file(s) at the same location as the input file(s),
        named as `filename_out.h5`.

    Examples
    --------
    For a single file:

    >>> Systems_selection_for_XRB('path/1e+00_Zsun_population.h5', 'BHNS_selected')
    ======================================================= 
    📁 Input  : 1e+00_Zsun_population.h5 
    💾 Output : BHNS_selected.h5 
    ─────────────────────────────────────────────────────── 
    Total systems (input):         10000000 
    Selected systems (output):     17500 
    Removed systems:               9982500 
    ─────────────────────────────────────────────────────── 
    Execution time:                8min 36sec 
    ======================================================= 

    For a list of files:

    >>> files_list = ['path/1e+00_Zsun_population.h5', 'path/1e-01_Zsun_population.h5', 'path/1e-04_Zsun_population.h5']
    >>> filename_out_list = ['NewName1', 'NewName2', 'NewName3']
    >>> Systems_selection_for_XRB(files_list, filename_out_list)
    

    """
    
    #===== Input file(s) validation =====
    if isinstance(path_file_h5, str):
        if path_file_h5.strip() == "":
            Pwarn("file path must be a non-empty string", category="InappropriateValueWarning")
            return None
        files = [path_file_h5]
    elif isinstance(path_file_h5, list):
        # List of file paths
        if not all(isinstance(f, str) and f.strip() != "" for f in path_file_h5):
            Pwarn("All elements in file list must be non-empty strings", 
                  category="InappropriateValueWarning")
            return None
        files = path_file_h5
    else:
        Pwarn("file must be a string or a list of strings (file paths)", category="InappropriateValueWarning")
        return None
    
    # Check that all files are .h5 files
    if not all(f.endswith('.h5') for f in files):
        Pwarn("All input files must be .h5 files", category="InappropriateValueWarning")
        return None
    
    #===== Output filename(s) validation =====
    if isinstance(filename_out, str):
        if filename_out.strip() == "":
            Pwarn("filename_out must be a non-empty string", category="InappropriateValueWarning")
            return None
        filenames = [filename_out]
    elif isinstance(filename_out, list):
        # List of filenames
        if not all(isinstance(fn, str) and fn.strip() != "" for fn in filename_out):
            Pwarn("All elements in filename_out list must be non-empty strings", 
                  category="InappropriateValueWarning")
            return None
        filenames = filename_out
    else:
        Pwarn("filename_out must be a string or a list of strings", category="InappropriateValueWarning")
        return None
    
    #===== Length consistency check =====
    if len(files) != len(filenames):
        Pwarn(f"Number of input files ({len(files)}) must match number of output filenames ({len(filenames)})", 
              category="InappropriateValueWarning")
        return None
    
    #===== Binaries extraction and export =====
    compact_types = {'NS', 'BH'}
    objects_to_avoid = {'NS', 'BH', 'massless_remnant', 'WD'}  # Avoid double CO, WD, massless remnants companions
    bad_system_states = {'disrupted', 'ERR', 'initial_RLOF'}

    out_fn = []

    import time

    for fh5, fn in zip(files, filenames):
        t_start = time.time()

        if verbose:
            print(f"\nProcessing file {fh5} in Population() ...")

        tmp_pop = Population(fh5)
        tmp_length = len(tmp_pop)
        tmp_oneline = tmp_pop.oneline[['state_f', 'S1_state_f', 'S2_state_f', 'event_f']]

        state1, state2, state_sys = tmp_oneline['S1_state_f'], tmp_oneline['S2_state_f'], tmp_oneline['state_f']
        
        mask = (
            (state1.isin(compact_types) & ~state2.isin(objects_to_avoid) & ~state_sys.isin(bad_system_states))
            | (state2.isin(compact_types) & ~state1.isin(objects_to_avoid) & ~state_sys.isin(bad_system_states))
        )
        idx = tmp_oneline.index[mask]

        del tmp_oneline

        # Export the selection to the same directory as the input file
        output_path = os.path.join(os.path.dirname(fh5) or '.', f"{fn}.h5")

        if verbose:
            print(f"Exporting selected systems to {output_path} ...")

        tmp_pop.export_selection(list(idx), output_path, overwrite=True)
        out_fn.append(output_path)

        del tmp_pop
        
        t_end = time.time()
        elapsed_time = t_end - t_start

        if elapsed_time < 60:
            time_str = f"{elapsed_time:.2f} seconds"
        else:
            time_str = f"{int(elapsed_time//60)} min {elapsed_time % 60:.1f} seconds"
        if verbose:
            pct = 100 * len(idx) / tmp_length if tmp_length > 0 else 0
            removed = tmp_length - len(idx)

            w = len(f"{tmp_length:,}")
            print(
                f"\n{'='*55}\n"
                f"  📁 Input  : {os.path.basename(fh5)}\n"
                f"  💾 Output : {os.path.basename(output_path)}\n"
                f"{'─'*55}\n"
                f"  {'Total systems (input):':<30} {tmp_length:>{w},}\n"
                f"  {'Selected systems (output):':<30} {len(idx):>{w},}  ({pct:.1f}%)\n"
                f"  {'Removed systems:':<30} {removed:>{w},}  ({100-pct:.1f}%)\n"
                f"{'─'*55}\n"
                f"  {'Execution time:':<30} {time_str:>{w}}\n"
                f"{'='*55}"
                )


#===== Functions for luminosity and XLF calculations =====

def Eddington_limit(accretor_mass, Donor_surface_H1, eta):
    """
    Compute the Eddington mass accretion rate from the accretor mass, its accretion efficiency and the content of hydrogen at the surface of the donor star.
    Also compute the Eddington luminosity.

    Parameters
    ----------
    accretor_mass : float
        Mass of the accretor in [M⊙].
    Donor_surface_H1 : float
        Hydrogen mass fraction of the donor star's surface (surface_h1_f).
    eta : float
        Radiative efficiency of the accretor (NS or BH), dimensionless.

    Returns
    -------
    Mdot_Edd_gps : float
        Eddington accretion rate in [g/s].
    L_edd : float 
        Eddington luminosity in [erg/s].
    
    """
    L_edd = 4.0 * np.pi * const.standard_cgrav * accretor_mass * const.Msun * const.clight /(0.2 * (1.0 + Donor_surface_H1)) # [erg/s]

    Mdot_Edd_gps = L_edd / (eta * const.clight**2) # [g/s]

    return Mdot_Edd_gps, L_edd


# Efficiencies 
def BH_efficiency(Spin): 
    """
    Compute the radiative efficiency of a Kerr black hole (or Schwarzschild BH if spin = 0)
    as a function of the spin parameter *a* and the mass according to Novikov-Thorne radiative efficiency [[1]_] .
    
    Parameters
    ---------- 
    Spin : float
        Spin of the black hole. Must be given as dimensionless spin "*a**" [c*J/(GM^2)].
        Its value can be given either positively (prograde accretion) or negatively (retrograde accretion) 
        (Spin  ∈ [-0.998, 0.998] with ± 0.998 being the Thorne limit).

    Returns
    -------
    tuple
        A tuple with two floats: (eta, r_ISCO).
        - eta (float): Radiative efficiency of the Kerr black hole, dimensionless.
        - r_ISCO (float): Radius of the ISCO in units of GM / c^2.

    References
    ----------
    .. [1] Bambi, C., Malafarina, D., & Tsukamoto, N. (2014).
            "Note on the radiation emitted by a thin accretion disk around a rotating compact object."
            *Phys. Rev. D*, 89(12), 127302. https://arxiv.org/abs/1406.2181
    """
    spin = np.asarray(Spin, dtype=float)

    # To control that the spin is in the range that gives the maximal efficiency
    a_star = np.clip(spin, -0.998, 0.998)

    #=== Compute efficiency [1]=== 
    
    # Use sign = +1 for spin ≥ 0 (prograde) and -1 otherwise (retrograde)
    sign = np.where(a_star >= 0.0, 1.0, -1.0)

    Z1 = 1 + (1 - a_star**2)**(1/3) * ((1 + a_star)**(1/3) + (1 - a_star)**(1/3))
    Z2 = np.sqrt(3 * a_star**2 + Z1**2)

    # Radius at the innermost stable circular orbit (ISCO)
    r_ISCO = 3 + Z2 - sign * np.sqrt((3 - Z1)*(3 + Z1 + 2 * Z2)) # Dimensionless

    # Energy at ISCO in units of mc^2 (i.e E_ISCO < 1, E = 1 means particle energy at infinity)
    E_ISCO = (r_ISCO**(3/2) - 2 * r_ISCO**(1/2) + a_star)/(r_ISCO**(3/4) * np.sqrt(r_ISCO**(3/2) - 3 * r_ISCO**(1/2) + 2 * a_star))
    
    # Radiative efficiency
    eta = 1 - E_ISCO

    return eta, r_ISCO

def NS_efficiency(Mass,log_radius, B_surf=None, accretor_logMdot=None):
    """
    Compute the radiative efficiency of a neutron star as a function of its mass and radius. 
    
    If there is a magnetic field, it will compare the magnetospheric radius with the radius of the neutron star to compute the accretion radius accordingly to [[1]_]. 
    The strength of the magntic field affects the magnetospheric radius of the NS and therefore, if the magnetospheric radius is larger than the NS radius, 
    the efficiency is lowered because the accretion disk is truncated and accretion happens along the magnetic lines up to the poles on smaller surface area.
    This results as a beaming effect of the X-ray emission due to smaller accretion area.

    Parameters
    ----------
    Mass : float
        Mass of the neutron star in [M⊙]
    log_radius : float
        Log10(Radius of the neutron star) in [R⊙]
    B_surf : float, optional
        Surface magnetic field strength of the neutron star in [G]. If provided along with accretor_logMdot, the magnetosphere radius will be computed.
    accretor_logMdot : float, optional
        Log10(Mass accretion rate onto the neutron star) in [M⊙/yr]. Required if B_surf is provided.
    
    Returns
    -------
    eta : float or ndarray
            Radiative efficiency of the neutron star, dimensionless.
    R : float or ndarray
            Accretion radius used for the efficiency computation in [cm]. It is either the NS radius or the magnetosphere radius if B_surf is provided and R_mag > R_NS.

    References
    ----------
    .. [1] Misra et al. (2024), Astronomy & Astrophysics, 682, A69,
           **Exploring the nature of Ultra-luminous X-ray sources across stellar population ages using detailed binary evolution calculations**
           https://doi.org/10.1051/0004-6361/202347880

    """
    
    #===== Mass =====
    mass_arr = np.asarray(Mass, dtype=float)
    
    mass = np.where(~np.isfinite(mass_arr) | (mass_arr <= 0.0), 1.45, mass_arr) 
        # Default 1.45 [M⊙] for Neutron stars if not given

    M = mass * const.Msun # [gr]

    #===== Radius =====
    logR_arr = np.asarray(log_radius, dtype=float)
    
    R_Ns_cm = np.where(np.isfinite(logR_arr), (10.0**(logR_arr)) * const.Rsun, 1.25e6) # Default radius [cm] for Neutron stars
        

    # If there is a magnetic field, compute the magnetosphere radius R_mag (POSYDON do not provide B field for now)
    R_mag = np.zeros_like(R_Ns_cm)
    if B_surf is not None and accretor_logMdot is not None:
        B_arr = np.asarray(B_surf, dtype=float)
        acc_log_mdot_arr = np.asarray(accretor_logMdot, dtype=float)
        mu = B_arr * R_Ns_cm**3  # Magnetic dipole moment [G cm^3]
        Mdot_NS = 10**acc_log_mdot_arr * const.Msun / const.secyer  # [g/s]
        R_mag = (mu**4 / (2 * const.standard_cgrav * M * Mdot_NS**2))**(1/7)  # [cm]
    
    # Choose the larger radius between R_NS and R_mag for the accretion radius.
    # If R_mag > R_Ns_cm, then magnetic field truncates the accretion disk, otherwise magnetosphere radius is within the NS radius
    R = np.where(R_mag > R_Ns_cm, R_mag, R_Ns_cm)


    #===== Efficiency computation =====
    eta = const.standard_cgrav * M /(R * const.clight**2) # dimensionless ~0.1-0.3
    
    return eta, R


# For winds (Bondi-Hoyle Littleton accretion)
def Stefan_Boltzmann_law(Luminosity, Radius):
    """
    Compute the effective temperature from the luminosity and radius of 
    a star using the Stefan-Boltzmann law.
    
    Parameters
    ----------
    Luminosity : float
        Luminosity of the star in [L⊙]
    Radius : float
        Radius of the star in [R⊙]

    Returns
    -------
    Teff : float
        Effective temperature of the star in [K].

    """

    #===== Check inputs =====
    if Luminosity <= 0.0 or Radius <= 0.0:
        Pwarn(f"Warning: Luminosity and Radius must be positive in L⊙ and R⊙. Returning NaN for Teff.",
              category="InappropriateValueWarning")
        return np.nan
    
    #===== Calculation =====
    Teff = ((Luminosity * const.Lsun) / (4.0 * np.pi * (Radius * const.Rsun)**2.0 * const.boltz_sigma))**(1.0/4.0)
    
    return Teff

def Wind_velocity_BHL(donor_mass, donor_logR, donor_logL, donor_log_wind_mdot, 
                      donor_surface_H1, donor_He_core_mass, 
                      Wind_velo_scheme ='Kudritzki+2000'):
    """
    Compute the wind velocity of the donor star in a detached binary system using either the prescription from :cite:`Hurley et al. 2002` or :cite:`Kudritzki & Puls 2000`. 
    
    This wind velocity is relevant for the calculation of the BHL accretion rate in **BHL_mdot_acc()** for detached systems, 
    where the accretion happens through winds. 

    The states of the binary systems must already be filtered as 'detached'. The filtering happens in the function L_bolo(), where the wind velocity is only computed for systems in detached state. 

    Parameters
    ----------
    donor_mass : array of float
        Mass of the donor star in [M⊙].
    donor_logR : array of float
        Log10(Radius of the donor star) in [R⊙].
    donor_logL : array of float
        Log10(Luminosity of the donor star) in [L⊙].
    donor_log_wind_mdot : array of float
        Log10(Mass loss rate of the donor star through winds) in [M⊙/yr].
    donor_surface_H1 : array of float
        Hydrogen mass fraction at the surface of the donor star.
    donor_He_core_mass : array of float
        Mass of the helium core of the donor star in [M⊙].
    Wind_velo_scheme : str, Default: 'Kudritzki+2000'
        Wind velocity calculation scheme to use. It can be either 'Hurley+2002' or 'Kudritzki+2000'.
        - 'Hurley+2002': It uses the prescription from Hurley et al. 2002, which depends on the mass, He core mass of the donor star, its radius as well as the abundance of hydrogen at surface. It distinguishes between H-rich and He-rich stars.
        - 'Kudritzki+2000' : It uses the prescription from Kudritzki & Puls 2000, which depends on the effective temperature of the donor star.

    Returns
    -------
    V_wind_out : array of float
        Wind velocity [cm/s].

    Raises
    ------
    ValueError
        If the wind velocity calculation scheme is not recognized.

    References
    ----------
    .. [1] Hurley, J. R., Tout, C. A., & Pols, O. R. 2002, MNRAS, 329, 897
    .. [2] Kudritzki, R.-P., & Puls, J. 2000, ARA&A, 38, 613
    .. [3] Sander A. A. C., Vink J. S., 2020, MNRAS, 499, 873

    """

    if Wind_velo_scheme not in ['Hurley+2002', 'Kudritzki+2000']:
        Pwarn(f"Wind_velo_scheme '{Wind_velo_scheme}' not recognized. Choose 'Hurley+2002' or 'Kudritzki+2000'.", 
              category="InappropriateValueWarning")
        return np.full_like(donor_mass, np.nan, dtype=float)

    #===== Loading arrays =====
    donor_mass_arr = np.asarray(donor_mass, dtype=float)
    donor_logR_arr = np.asarray(donor_logR, dtype=float)
    donor_logL_arr = np.asarray(donor_logL, dtype=float)
    donor_log_wind_mdot_arr = np.asarray(donor_log_wind_mdot, dtype=float)
    donor_surface_H1_arr = np.asarray(donor_surface_H1, dtype=float)
    donor_He_core_mass_arr = np.asarray(donor_He_core_mass, dtype=float)

    arrays = [donor_mass_arr, donor_logR_arr, donor_logL_arr, donor_log_wind_mdot_arr,
              donor_surface_H1_arr, donor_He_core_mass_arr]
    shapes_non_scalar = {arr.shape for arr in arrays if arr.shape != ()}
    any_scalar = any(arr.shape == () for arr in arrays)

    # If shapes mismatch or scalars are mixed with arrays, warn and bail out
    if (any_scalar and shapes_non_scalar) or len(shapes_non_scalar) > 1:
        Pwarn(f"Wind_velocity_BHL: shape mismatch in inputs: {[arr.shape for arr in arrays]}. Returning NaN.",
              category="InappropriateValueWarning")
        if shapes_non_scalar:
            shape = next(iter(shapes_non_scalar))
            return np.full(shape, np.nan, dtype=float)
        return np.nan

    #===== Calculations =====
    V_wind_out = np.zeros_like(donor_mass_arr, dtype=float)

    for idx in np.ndindex(donor_mass_arr.shape):
        donor_mass_i = donor_mass_arr[idx]
        donor_logR_i = donor_logR_arr[idx]
        donor_logL_i = donor_logL_arr[idx]
        donor_log_wind_mdot_i = donor_log_wind_mdot_arr[idx]
        donor_surface_H1_i = donor_surface_H1_arr[idx]
        donor_He_core_mass_i = donor_He_core_mass_arr[idx]

        Teff = Stefan_Boltzmann_law(10.0**donor_logL_i, 10.0**donor_logR_i) # [K]

        #===== Wind velocity calculation schemes =====
        # [1] Hurley et al. 2002
        if Wind_velo_scheme == 'Hurley+2002':
            # For H-rich stars
            if donor_surface_H1_i > 0.01:
                if donor_He_core_mass_i > 0.0 and 10.0**donor_logR_i > 900.0:
                    beta = 0.125
                elif donor_mass_i > 120.0:
                    beta = 7.0
                elif donor_mass_i < 1.4:
                    beta = 0.5
                else:
                    beta = 0.5 + (donor_mass_i - 1.4) / (120.0 - 1.4) * (6.5)
            else:  # For He-rich stars
                if donor_mass_i > 120.0:
                    beta = 7.0
                elif donor_mass_i < 10.0:
                    beta = 0.125
                else:
                    beta = 0.125 + (donor_mass_i - 10.0) / (120.0 - 10.0) * (6.875)
            f_m = float(np.sqrt(beta))

        # [2] Kudritzki & Puls 2000
        elif Wind_velo_scheme == 'Kudritzki+2000':
            if Teff >= 21000:
                f_m = 2.65
            elif Teff <= 10000:
                f_m = 1.0
            else:
                f_m = 1.4
        else:
            raise ValueError(f"Invalid scheme '{Wind_velo_scheme}'. Choose 'Hurley+2002' or 'Kudritzki+2000'.")

        V_esc = np.sqrt(2.0 * const.standard_cgrav * donor_mass_i * const.Msun / (10**donor_logR_i * const.Rsun)) # [cm/s]
        V_wind_val = f_m * V_esc  # [cm/s]

        # [3] Sander A. A. C., Vink J. S., 2020, MNRAS, 499, 873 
        if donor_surface_H1_i < 0.4 and Teff >= 1e4:
            if donor_log_wind_mdot_i >= -5.25:
                slope = (3.7 - 3.25) / (-2.5 + 5.25)
            else:
                slope = (3.25 - 3.75) / (-5.25 + 7.25)
            V_wind_val = 10 ** (slope * donor_log_wind_mdot_i + 3.25 + 5.25 * slope) * 1e5  # Convert from km/s to cm/s
        else:
            pass

        V_wind_out[idx] = V_wind_val

    return V_wind_out

def BHL_mdot_acc(Macc_Msun, donor_mass_Msun, accretor_type, donor_log_radius,
                 donor_log_wind_mdot, a_Rsun, ecc, V_wind, spin, Donor_RL_radius, **kwargs):
    """
    Compute the `Bondi-Hoyle-Lyttleton (BHL)`accretion rate for a detached binary system where the accretion happens through winds.
    These winds are assumed to be isotropic and the accretion is assumed to be spherically symmetric.

    Parameters
    ----------
    Macc_Msun : array of float
        Mass of the accretor in [M⊙].
    donor_mass_Msun : array of float
        Mass of the donor star in [M⊙].
    accretor_type : array of str
        Type of the accretor. It should be 'NS' or 'BH'.
    donor_log_radius : array of float
        Log10(Radius of the donor star) in [R⊙].
    donor_log_wind_mdot : array of float
        Log10(Mass loss rate of the donor star through winds) in [M⊙/yr].
    a_cm : array of float
        Orbital separation in [cm].
    ecc : array of float
        Eccentricity of the orbit (dimensionless).
    V_wind : array of float
        Wind velocity of the donor star in [cm/s]. Comes from the function **Wind_velocity_BHL()**.
    spin : array of float
        Spin of the accretor. Only required for BH. It is the dimensionless spin parameter a* = cJ/(GM^2).
    Donor_RL_radius : array of float
        Roche lobe radius of the donor star in [R⊙]. It is used to check if the donor star is close to filling its Roche lobe, which can increase the accretion efficiency due to the formation of a disk around the accretor in wind-fed systems.
    **kwargs : dict, optional
        - alpha_BHL (float, *Default = 1.5*) : Dimensionless efficiency parameter for the BHL accretion rate and it should be in the range [1.0, 2.0]. 
            It accounts for the efficiency of the accretion process and can be used to calibrate the BHL accretion rate with observations or more detailed simulations.
        - Wind_Disk_BH_formation (str, *Default = None*) : If 'Senk2021' or 'HiraiMandel21', it check if some criteria are met for the infalling material to form a disk around the BH before being accreted, which can increase the accretion efficiency [[3]_] [[5]_]. `Senk2021` checks if the circularization radius of the infalling material is larger than the ISCO radius of the BH, while `HiraiMandel21` checks if the donor star fills at least 80% of its Roche lobe. 

    Returns
    -------
    Mdot_BHL_gps : array of float
        BHL accretion rate in [g/s].

    References
    ----------
    .. [1] Bondi, H., & Hoyle, F. 1944, MNRAS, 104, 273
    .. [2] Frank, J., King, A., & Raine, D. (2002). "Accretion Power in Astrophysics: Third Edition" Chapter 4.9 pages 73-74.
           *Cambridge University Press*. https://doi.org/10.1017/CBO9781139164417
    .. [3] Sen, K. ,Xu, X. -T., Langer, N., El Mellah, I. , Schurmann, C., Quast, M., 2021, 
           *X-ray emission from BH+O star in binaries expected to descend from the observed galactic WR+O binaries*,
           A&A 652, A138 (2021). https://doi.org/10.1051/0004-6361/202141214
           Equation 10
    .. [4] Bambi, C., Malafarina, D., & Tsukamoto, N. (2014).
            "Note on the radiation emitted by a thin accretion disk around a rotating compact object."
            *Phys. Rev. D*, 89(12), 127302. https://arxiv.org/abs/1406.2181
    .. [5] Hirai, R. & Mandel, I. (2021).  
           *Conditions for accretion disc formation and observability of wind-accreting X-ray binaries.*. 
           Publications of the Astronomical Society of Australia, 38, e056.  
           https://doi.org/10.1017/pasa.2021.53
    """

    #===== Optional parameters =====
    allowed_keys = {'alpha_BHL', 'Wind_Disk_BH_formation'}
    unknown = sorted(set(kwargs.keys()) - allowed_keys)
    
    alpha_BHL = kwargs.get('alpha_BHL', 1.5)
    if not (1.0 <= alpha_BHL <= 2.0):
        Pwarn("alpha_BHL must be in the range [1.0, 2.0]", category="InappropriateValueWarning")

    wind_disk_BH = kwargs.get('Wind_Disk_BH_formation', None)
    if not isinstance(wind_disk_BH, str) or wind_disk_BH not in [None, 'Senk2021', 'HiraiMandel21']:
        Pwarn("Wind_Disk_BH_formation must be None, 'Senk2021' or 'HiraiMandel21'.", category="InappropriateValueWarning")
        wind_disk_BH = None

    if unknown:
        Pwarn(f"BHL_mdot_acc: unknown parameter(s): {unknown} (PID={os.getpid()})", category="InappropriateValueWarning")

    #===== Array conversion =====
    acc_mass_arr = np.asarray(Macc_Msun, dtype=float)
    don_mass_arr = np.asarray(donor_mass_Msun, dtype=float)
    acc_type_arr = np.asarray(accretor_type, dtype=object)
    logR_donor_arr = np.asarray(donor_log_radius, dtype=float)
    logMdot_w_arr = np.asarray(donor_log_wind_mdot, dtype=float)
    a_arr_cm = np.asarray(a_Rsun, dtype=float) * const.Rsun # Convert from [R⊙] to [cm]
    ecc_arr = np.asarray(ecc, dtype=float)
    v_wind_arr = np.asarray(V_wind, dtype=float)
    spin_arr = np.asarray(spin, dtype=float)
    donor_RL_arr = np.asarray(Donor_RL_radius, dtype=float) * const.Rsun # [cm]

    #===== Calculation Mdot BHL =====
    Mdot_BHL_gps = np.zeros_like(acc_mass_arr, dtype=float)

    for idx in np.ndindex(acc_mass_arr.shape):

        #=== Needed parameters ===
        Macc_i = float(acc_mass_arr[idx])
        Mdon_i = float(don_mass_arr[idx])
        acc_type_i = acc_type_arr[idx]
        Rdon_i = float(10**logR_donor_arr[idx] * const.Rsun) # Convert from [R⊙] to [cm]
        logMdot_w_i = float(logMdot_w_arr[idx])
        mdot_w_i = float(10**logMdot_w_i)
        a_i = float(a_arr_cm[idx])
        ecc_i = float(np.clip(ecc_arr[idx], 0.0, 0.999)) if np.isfinite(ecc_arr[idx]) else 0.0
        v_wind_i = float(v_wind_arr[idx])
        spin_i = float(spin_arr[idx])
        donor_RL_i = float(donor_RL_arr[idx])

        if not np.isfinite(a_i) or a_i <= 0.0 or not np.isfinite(v_wind_i) or v_wind_i < 0.0:
            Mdot_BHL_gps[idx] = 0.0
            continue

        # Masses to cgs
        Macc_cgs = Macc_i * const.Msun
        Mdon_cgs = Mdon_i * const.Msun

        #=== Orbital motion parameters ===

        # Orbital frequency [rad/s]
        Omega = np.sqrt(const.standard_cgrav * (Mdon_cgs + Macc_cgs) / Rdon_i**3)
        # Random orbital phase [s]
        t0 = np.random.rand() * (2.0 * np.pi / Omega)
        # Eccentric anomaly E at time t0 from Kepler's equation
        E = newton(lambda x: x - ecc_i * np.sin(x) - Omega * t0, np.pi/2, maxiter=100)
        b = a_i * np.sqrt(1.0 - ecc_i**2)
        r_vec = np.array([a_i * (np.cos(E) - ecc_i), b * np.sin(E)]) # Vectorial form radius
        v_dir = np.array([-a_i * np.sin(E), b * np.cos(E)]) # Vectorial form velocity (derivative of r_vec)
        R = np.linalg.norm(r_vec, axis=0) # Separation between stars [cm]
        v_dir_norm = np.linalg.norm(v_dir, axis=0)
        if R > 0.0 and v_dir_norm > 0.0:
            cos_angle = np.dot(v_dir, r_vec) / (v_dir_norm * R)
            cos_angle = np.clip(cos_angle, -1.0, 1.0)
        else:
            cos_angle = 0.0

        V_orb = np.sqrt(const.standard_cgrav * (Mdon_cgs + Macc_cgs) * (2.0 / R - 1.0 / a_i)) # Orbital velocity [cm/s]

        # Total relative velocity [cm] of winds [2]
        V_rel = np.sqrt(V_orb**2 + v_wind_i**2 + 2.0 * V_orb * v_wind_i * cos_angle) 

        #=== Mdot BHL calculation ===

        Mdot_w_cgs = mdot_w_i * const.Msun / const.secyer

        # From [1], calculation of BHL accretion rate
        if v_wind_i == 0.0:
            frac = 1.0
            Pwarn(
                f"Warning: Wind velocity is zero, setting accretion fraction to 1.0.",
                category="InappropriateValueWarning",
            )
        elif V_rel == 0.0:
            frac = 0.0
            Pwarn(
                f"Warning: Orbital velocity is zero, setting BHL accretion rate to 0.0 g/s.",
                category="InappropriateValueWarning",
            )
        elif R == 0.0:
            frac = 0.0
            Pwarn(
                f"Warning: Orbital separation is zero, setting BHL accretion rate to 0.0 g/s.",
                category="InappropriateValueWarning",
            )
        else:
            frac = alpha_BHL * (const.standard_cgrav * Macc_cgs)**2 / (2 * V_rel**3 * v_wind_i * R**2)

        Mdot_BHL = frac * Mdot_w_cgs

        # From [3], check if an accretion disk can form around the BH accretor
        if (wind_disk_BH == 'SenK2021')and (acc_type_i == 'BH'):
            eta_bh, r_ISCO_bh = BH_efficiency(spin_i)
            gamma = r_ISCO_bh / 6.0 # From [4]
            q = Mdon_cgs / Macc_cgs
            Rdisk_div_RISCO = (
                (2/3) * (eta_bh / (1 + q))**2 * (V_orb / const.clight)**(-2)
                * (1 + (v_wind_i / V_orb)**2)**(-4) / gamma
            )
            if Rdisk_div_RISCO < 1.0:
                Mdot_BHL = 1e-99 # Avoid division by 0

        # [5] Hirai & Mandel 2021
        if (wind_disk_BH == 'HiraiMandel21'):
            mask_hirai = (Rdon_i / donor_RL_i < 0.8) & (acc_type_i == 'BH')
            Lx = np.where(mask_hirai, 10**(-54), Lx)
    

        Mdot_BHL_gps[idx] = float(Mdot_BHL)

    return Mdot_BHL_gps if Mdot_BHL_gps.shape else float(Mdot_BHL_gps)


# Luminosities calculations
def Be_XRB_Lx(period, filtered = True):
    """
    Compute the X-ray luminosity of a system containing a Be star as companion. The donor star should have [[1]_]:
    - Surface velocity higher than 70% of its critical velocity
    - Orbital period between 10 and 300 days
    - The system state should be `detached` and not undergoing RLO
    - Donor mass of at least 6 M⊙
    - A decretion disc 100 times bigger than the star itself and it must go beyond the Roche-Lobe radius of the donor star

    The input parameter **filtered** is to check whether the population has already been filtered to contain only Be systems according to the description. 
    If the population is not already filtered, returns Lx = 0.0. 
    The `XRB_selection_function()`can be used to filter the population to contain only Be systems based on the criteria described above.

    Be stars have fast rotation surface which, when the surface velocity is close to the critical velocity,
    can form a decretion disk around them. The compact object (NS or BH) can then accrete material from this disk
    if it passes through the Roche lobe of the Be star and thus leading to X-ray emission when accreted.

    Parameters
    ----------
    orbital_period_days : float
        Orbital period of the binary system in [days].
    filtered :  bool (**Default** = True)
        If the population has already been filtered to contain only Be systems according to the description. 
        If the population is not already filtered, returns Lx = 0.0. 

    Returns
    -------
    Lx : float
        X-ray luminosity of the system [erg/s].
    
    References
    ----------
    .. [1] Misra, D., Kovlakas, K., Fragos, T., Lazzarini, M., Bavera, S. S.,  
           Lehmer, B. D., Zezas, A., Zapartas, E., Xing, Z., Andrews, J. J., Dotter, A.,  
           Rocha, K. A., Srivastava, P. M., Sun, M. (2023).  
           *X-ray luminosity function of high-mass X-ray binaries: Studying the signatures of different physical processes using detailed binary evolution calculations*. 
           https://arxiv.org/abs/2209.05505v2 
    .. [2] Dai, Hai-Lang ;  Liu, Xi-Wei ;  Li, Xiang-Dong (2006), Equation (11). 
           *Exploration of the Ps-Porb Relation for Wind-fed X-Ray Pulsars*.
           https://doi.org/10.1086/508735
    """

    if not filtered: 
        Pwarn(f"Be_XRB_Lx: Population not filtered for Be systems. Returning Lx = 0.0 for all systems.", 
              category="InappropriateValueWarning", stacklevel=2)
        Lx = np.zeros_like(period, dtype=float)
        return Lx

    if np.any(period <= 0.0):
        Pwarn(f"Be_XRB_Lx: Orbital period should be positive. Returning NaN for non-positive periods.", 
              category="InappropriateValueWarning", stacklevel=2)
        # [2] Equation 11
        Lx = np.where(period > 0.0, (10.0**35) * 10.0**(4.53 - (1.50 * np.log10(period))), np.nan)

    return Lx

def Luminosity_GRRMHD_systems(Accretor_LogMdot, spin, Mdot_Edd, L_edd, **kwargs):
    """
    Calculate the luminosity of Black holes accreting material through Roche lobe overflow.

    The population must have been simulated by using General Relativistic Radiative Magneto-HydroDynamic (GRRMHD) grids (`step_BH_HMS_RLO` and `step_NS_HMS_RLO`).

    Parameters
    ----------
    Accretor_LogMdot : array-like
        Logarithm (base 10) of the mass accretion rate in units of solar masses per year (Msun/yr).
    spin : array-like
        Dimensionless spin parameter of the black hole.
    Mdot_Edd : array-like
        Eddington mass accretion rate in grams per second (g/s).
    L_edd : array-like
        Eddington luminosity in erg/s.
    **kwargs : dict, optional
        - `spin_split` (float, *Default* = 0.75): A threshold value for the spin parameter to determine which interpolation curve to use for the beaming factor.
    
    Returns
    -------
    L_obs : pandas Series
        Observed luminosity in erg/s, accounting for GRRMHD effects and random viewing angles.
    b : pandas Series
        Beaming factor due to GRRMHD effects, as a function of the random viewing angle.
    Weight_times_beaming : bool
        A flag indicating whether the weights of the systems should be multiplied by the beaming factor when plotting the cumulative luminosity function (CCDF). 
        For GRRMHD systems, this is set to False because the luminosity already includes the beaming factor, so we do not need to multiply by it again in `plot_ccdf()`.
    
    References
    ----------
    .. [1] Tom M. Kwan, Dai Lixin, Zepei Xing, Tao Ji, Tassos Fragos, M. Middleton
           *Strongly magnetized super-Eddington acretion flow  around black holes I: Effects of the black hole mass, spin and mass accretion rate*
           In preparation, 2026

    Notes
    -----
    - The luminosity is calculated by interpolating the beaming factor from GRRMHD simulations as a function of the mass accretion rate and the spin of the black hole. 
    For the viewing angle, a random value of the cosine of the angle is drawn from a uniform distribution between 1 and 0.6428 (corresponding to an angle of 50 degrees - cos(50*np.pi/180)), 
    which corresponds to an isotropic distribution of the viewing angles. Then the drawn value is passed in the arccosine function to get the angle in radians, 
    giving  more edge-on views (higher angles) than face-on views (lower angles). This is what we expect for a population of sources with random orientations in the sky. 
    The beaming factor is then interpolated from the GRRMHD simulations as a function of the mass accretion rate and the spin of the black hole, 
    and it is used to compute the observed luminosity as L_obs = b * L_edd * (Mdot / Mdot_Edd).
    - `spin_split` is used to determine which interpolation curve to use for the beaming factor. All spin values below `spin_split` will use 
    the interpolation curve for spin = 0.0, while all spin values above `spin_split` will use the interpolation curve for spin = 0.9.
    The relation mass accretion rate-spin of black hole being not linear, to spin up a black hole, e.g. from 0.4 to 0.6, it does not require the same amount of mass
    accretion as from 0.6 to 0.8, and therefore the beaming factor does not change linearly with the spin of the black hole. 
    Giving a high value for `spin_split` would better split the interpolation curves in the region of the spin where the beaming factor changes 
    more rapidly, which is at high spin, and therefore it would give a better interpolation 
    of the beaming factor as a function of the spin of the black hole.
    - We also interpolate the beaming factor according to the mass accretion rate. Two interpolation curves are used for each spin value,
    one for low mass accretion rates (mdot ~10), and one for high mass accretion rates (mdot ~200). If the mass accretion rate of a system is between these two values,
    we interpolate the beaming factor between the two curves according to the mass accretion rate of the system. 
    If the mass accretion rate is below 10, we use the low mass accretion rate curve, and if it is above 200, we use the high mass accretion rate curve.
    - Systems with mass accretion rates above 2000 are not correctly interpolated, so we set an upper limit to 2000 to avoid incorrect values.
    
    
    """
    
    spin_split = kwargs.get('spin_split', 0.75)
    if not (0.0 <= spin_split <= 0.9):
        Pwarn(f"spin_split should be between 0.0 and 0.9. Using default value of 0.8.", category="InappropriateValueWarning")
        spin_split = 0.75
    
    Spin = np.asarray(spin, dtype=float)
    Mdot_Edd = np.asarray(Mdot_Edd, dtype=float)
    L_edd = np.asarray(L_edd, dtype=float)

    #===== Convert LogMdot from Msun/yr to g/s =====
    Acc_Log_Mdot_MsunYr = np.asarray(Accretor_LogMdot) 
    Acc_Mdot_gps = 10**Acc_Log_Mdot_MsunYr * const.Msun / const.secyer

    # File path to this csv file containig data points of b VS theta for different spins and mass accretion rates. 
    pathfile = 'Data_points_b_VS_theta.csv'
    # Column names in Data_points_b_VS_theta_Kwan.csv file are :
    # 'x_a0_low', 'y_a0_low', 'x_a0_high', 'y_a0_high' for spineless BH
    # 'x_a09_low', 'y_a09_low', 'x_a09_high', 'y_a09_high' for 0.9 spin BH
    # Low refers to mdot ~ 10 and high refers to mdot ~200

    #=============================
    # Interpolation beaming factor 
    # ============================

    from scipy.interpolate import CubicSpline
    df = pd.read_csv(pathfile)
    
    rng = np.random.default_rng(seed=33)
    # Random viewing angle[rad] from an isotropic distribution of the cosine 
    theta = np.arccos(rng.uniform(np.cos(50*np.pi/180), 1, size=len(Accretor_LogMdot)))

    from scipy.interpolate import CubicSpline

    cs0_low = CubicSpline(df['x_a0_low'].dropna(), df['y_a0_low'].dropna())
    cs0_high = CubicSpline(df['x_a0_high'].dropna(), df['y_a0_high'].dropna())
    cs09_low = CubicSpline(df['x_a09_low'].dropna(), df['y_a09_low'].dropna())
    cs09_high = CubicSpline(df['x_a09_high'].dropna(), df['y_a09_high'].dropna())

    # --- Evaluate splines ---
    b0_low  = cs0_low(theta)
    b0_high = cs0_high(theta)
    b9_low  = cs09_low(theta)
    b9_high = cs09_high(theta)

    #===== Min and max mdot for interpolation from [1], table 2 =====
    # All in Eddington units
    mdot_min_a0 = 8
    mdot_min_a09 = 7
    mdot_max = 192 

    Mdot_acc = np.divide(Acc_Mdot_gps, Mdot_Edd, 
                         out=np.full_like(Acc_Mdot_gps, np.nan, dtype=float), 
                         where=Mdot_Edd > 0)
    # As mentionned in [1], above mdot ~ 2000, value must not be correctly interpolated, 
    # so we set an upper limit to 2000 to avoid incorrect values.
    Mdot_acc = np.clip(Mdot_acc, 1e-10, 2000)
    log_mdot = np.log10(Mdot_acc)

    #===== Radiative luminosity =====
    # Spinless L_rad for spin <= spin_split, Spin dependent L_rad for spin > spin_split
    L_rad = np.where(Spin <= spin_split,  0.19 * Mdot_acc**(0.95), 0.73 * Mdot_acc**(1.13)) * L_edd

    # --- Weights ---
    w0 = (np.log10(mdot_max) - log_mdot) / (np.log10(mdot_max) - np.log10(mdot_min_a0))
    w9 = (np.log10(mdot_max) - log_mdot) / (np.log10(mdot_max) - np.log10(mdot_min_a09))

    w0 = np.clip(w0, 0, 1)
    w9 = np.clip(w9, 0, 1)

    # --- Interpolation ---
    b0 = 10**(w0 * np.log10(b0_low) + (1 - w0) * np.log10(b0_high))
    b9 = 10**(w9 * np.log10(b9_low) + (1 - w9) * np.log10(b9_high))

    #===== Final beaming factor =====
    b = np.where(Spin <= spin_split, b0, b9)
    
    Weight_times_beaming = False 

    L_obs = L_rad * b # erg/s

    bad = (~np.isfinite(L_obs)) | (L_obs <= 0.0)
    if np.any(bad): # To avoid a error in log10
        L_obs = np.where(bad <= 0.0, 1e-54, L_obs)
        
    return L_obs, b, Weight_times_beaming

def L_bolo(used_grid, bin_state, lg_mtransfer_rate, acc_log_mdot, donor_log_wind_mdot, acc_state, donor_log_L, 
           Macc, Mdonor, Acc_log_radius, Donor_log_radius, Donor_surface_H1, donor_He_core_mass, 
           Donor_RL_radius, a_Rsun, ecc, spin, B_surface, **kwargs):
    """ 
    This function calculates the bolometric luminosity (L_bolo) of a binary system based on its final 
    state and various physical parameters of the accretor and donor stars. 
    
    It calculates the geometrical beaming factor of a super-Eddington accreting system following 
    the prescription of King et al. (2001) and King (2009) or General Relativistic, Radiation Magneto-Hydrodynamics (GRRMHD) model. 
    This is an observational effect due to the inflated structure of the accretion disk that beams 
    the outgoing X-ray emission. The function also checks for the presence of high wind rates based 
    on the donor star filling at least 80% of its Roche lobe following the prescription of 
    *Hirai & Mandel (2021)* to compute the luminosity in wind-fed systems.

    Requires that the input arrays have been through the `XRB_selection_function()` to ensure that the 
    systems are X-ray binaries containing a NS or BH accretor and a non-degenerate donor star. 
    The binary state is either 'RLO1', 'RLO2' for Roche lobe overflow systems or 'detached' wind-fed systems.

    Parameters
    ----------
    used_grid : str, default='classical'
        Whether simulations have been performed using classical grids or with GRRMHD grids for RLO-BH systems. It can be either 'classical' or 'GRRMHD'.
    bin_state : array of str
        The final state of the binary system (e.g. 'RLO1', 'RLO2', 'detached', etc.).
    lg_mtransfer_rate : array of float
        Log10 of the mass transfer rate through Roche Lobe Overflow (RLO) in [M⊙/yr].
    acc_log_mdot : array of float
        Log10 of the mass changing rate of the accretor in [M⊙/yr].
    donor_log_wind_mdot : array of float
        Log10 of the absolute wind mass loss rate of the donor star in [M⊙/yr].
    acc_state : array of str
        The type of the accretor (e.g. 'NS' or 'BH').
    donor_log_L : array of float
        Log10 of the luminosity of the donor star in [L⊙].
    Macc : array of float
        Mass of the accretor in [M⊙].
    Mdonor : array of float
        Mass of the donor star in [M⊙].
    Acc_log_radius : array of float
        Log10 of the radius of the accretor in [R⊙].
    Donor_log_radius : array of float
        Log10 of the radius of the donor star in [R⊙].
    Donor_surface_H1 : array of float
        Hydrogen mass fraction at the surface of the donor star (surface_h1_f).
    donor_He_core_mass : array of float
        Mass of the helium core of the donor star in [M⊙].
    Donor_RL_radius : array of float
        Roche lobe radius of the donor star in [R⊙].
    a_Rsun : array of float
        Orbital separation in [R⊙].
    ecc : array of float
        Orbital eccentricity (dimensionless).
    spin : array of float
        Spin of the accretor (dimensionless).
    B_surface : array of float
        Surface magnetic field of the accretor in [G]. Only required for NS accretors.
    **kwargs : dict, optional
        Keywords arguments are passted to the downstream functions.

    Returns
    -------
    Tuple of arrays
        Lx, Edd_state, beaming, Weight_times_beaming
        - Lx: X-ray luminosity of the system in [erg/s].
        - Edd_state: 'Sub-Eddington' or 'Super-Eddington' depending on whether the accretion rate is below or above the Eddington limit.
        - beaming: Beaming factor (b) applied to the luminosity in the super-Eddington regime, dimensionless.
        - Weight_times_beaming: A flag indicating whether the weights of the systems should be multiplied by the beaming factor when plotting the cumulative luminosity function (CCDF). 
        For GRRMHD systems, it is set to False because the luminosity already includes the beaming factor, so we do not need to multiply by it again in `plot_ccdf()`. 
        For classical systems, it is set to True because the luminosity does not include the beaming factor, so we need to multiply by it in `plot_ccdf()`.

    References
    ----------
    .. [1] Misra, D., Kovlakas, K., Fragos, T., Lazzarini, M., Bavera, S. S.,  
           Lehmer, B. D., Zezas, A., Zapartas, E., Xing, Z., Andrews, J. J., Dotter, A.,  
           Rocha, K. A., Srivastava, P. M., Sun, M. (2023).  
           *X-ray luminosity function of high-mass X-ray binaries: Studying the signatures of different physical processes using detailed binary evolution calculations*.  
           https://arxiv.org/abs/2209.05505v2
    .. [2] King, A. R., Davies, M. B., Ward, M. J., Fabbiano, G., & Elvis, M. (2001).  
           *Ultraluminous X-Ray Sources in External Galaxies.*. 
           The Astrophysical Journal, 552(2), L109. 
           https://doi.org/10.1086/320343
    .. [3] King, A. R. (2008).  
           *Accretion rates and beaming in ultraluminous X-ray sources*. 
           Monthly Notices of the Royal Astronomical Society, 385(1), L113–L115. 
           https://doi.org/10.1111/j.1745-3933.2008.00444.x
    .. [4] Shakura, N. I. & Sunyaev, R. A. (1973).  
           *Black holes in binary systems. Observational appearance.*. 
           Astronomy and Astrophysics, 24, 337.  
           https://ui.adsabs.harvard.edu/abs/1973A%26A....24..337S
    
    """

    #===== Optional parameters =====
    allowed_keys = {'alpha_BHL', 'Wind_Disk_BH_formation', 'Wind_velo_scheme', 'spin_split'}
    unknown = sorted(set(kwargs.keys()) - allowed_keys)
    if unknown:
            Pwarn(f"L_bolo(): unknown parameter(s): {unknown}", category="InappropriateValueWarning")

    # Passed to BHL_mdot_acc function
    BHL_kwargs = {key: kwargs[key] for key in ['alpha_BHL', 'Wind_Disk_BH_formation'] if key in kwargs}

    GRRMHD_kwargs = {key: kwargs[key] for key in ['spin_split'] if key in kwargs}

    # Passed to Wind_velocity_BHL function
    Wind_velo_scheme = kwargs.get('Wind_velo_scheme', 'Kudritzki+2000')

    # Grid type 
    if used_grid not in ['GRRMHD', 'classical']:
        Pwarn(f"used_grid '{used_grid}' not recognized. Choose 'GRRMHD' or 'classical'.", category="InappropriateValueWarning")
        used_grid = 'classical'

    #===== Array conversion =====
    bin_state = np.asarray(bin_state, dtype=object)
    acc_state = np.asarray(acc_state, dtype=object)
    donor_log_L = np.asarray(donor_log_L, dtype=float)
    Macc = np.asarray(Macc, dtype=float)
    Mdonor = np.asarray(Mdonor, dtype=float)
    Acc_log_radius = np.asarray(Acc_log_radius, dtype=float)
    Donor_log_radius = np.asarray(Donor_log_radius, dtype=float)
    Donor_surface_H1 = np.asarray(Donor_surface_H1, dtype=float)
    donor_He_core_mass = np.asarray(donor_He_core_mass, dtype=float)
    Donor_RL_radius = np.asarray(Donor_RL_radius, dtype=float)
    a_Rsun = np.asarray(a_Rsun, dtype=float)
    ecc = np.asarray(ecc, dtype=float)
    spin = np.asarray(spin, dtype=float)
    lg_mtransfer_rate = np.asarray(lg_mtransfer_rate, dtype=float)
    acc_log_mdot = np.asarray(acc_log_mdot, dtype=float)
    donor_lgmdot_wind = np.asarray(donor_log_wind_mdot, dtype=float)
    B_surf = np.asarray(B_surface, dtype=float)

    #===== Masks =====
    mask_RLO = (bin_state == 'RLO1') | (bin_state == 'RLO2')
    mask_detached = (bin_state == 'detached')
    is_BH = (acc_state == 'BH')

    # Whether to multiply the weights of the systems by the beaming factor when plotting the cumulative luminosity function (CCDF).
    Weight_times_beaming = np.full_like(bin_state, True, dtype=bool)

    #===== Eddington parameters =====
    eta_ns, _ = NS_efficiency(Macc, Acc_log_radius, B_surf=B_surf, accretor_logMdot=acc_log_mdot)
    eta_bh, _ = BH_efficiency(spin)
    eta = np.where(is_BH, eta_bh, eta_ns)

    Medd, L_edd = Eddington_limit(Macc, Donor_surface_H1, eta)
    mask_valid = np.isfinite(Medd) & (Medd > 0.0)

    Mdot = np.full_like(eta, np.nan, dtype=float)
    b = np.ones_like(Mdot, dtype=float)

    ratio = np.zeros_like(Mdot, dtype=float)

    #===== Wind calculation for detached systems =====
    Vwind = Wind_velocity_BHL(Mdonor, Donor_log_radius, donor_log_L, donor_lgmdot_wind, 
                              Donor_surface_H1, donor_He_core_mass, Wind_velo_scheme=Wind_velo_scheme)
    Mdot_BHL = BHL_mdot_acc(Macc, Mdonor, acc_state, Donor_log_radius, donor_lgmdot_wind, a_Rsun, ecc, Vwind, spin, Donor_RL_radius, **BHL_kwargs)

    #Winds-NS/BH
    Mdot = np.where(mask_detached, Mdot_BHL, Mdot)

    #===== Grid splitting =====
    if used_grid == 'classical':
        
        #===== Mass transfer rate through RLO =====
        Mdot_RLO_Msunyr = np.where(np.isfinite(lg_mtransfer_rate), 10.0**lg_mtransfer_rate, 0.0) # [M⊙/yr]
        Mdot_RLO_gps = Mdot_RLO_Msunyr * const.Msun / const.secyer  # [g/s]

        #===== Mdot calculation =====
        # RLO
        Mdot = np.where(mask_RLO, Mdot_RLO_gps, Mdot)

        ratio[mask_valid] = Mdot[mask_valid] / Medd[mask_valid]
    
        #===== Luminosity calculation =====
        # Follows King et al. (2008) [2] & [3]
        Lx_sub = Mdot * eta * const.clight**2 # Sub-Eddington case

        # Beaming
        mask_beam = ratio > 8.5
        if np.any(mask_beam):
            b[mask_beam] = np.maximum(73.0 / (ratio[mask_beam]**2), 3.2e-3)

        Lx_super = (L_edd / b) * (1.0 + np.log(np.maximum(ratio, 1.0)))

        Lx = np.where(ratio <= 1.0, Lx_sub, Lx_super)
        Lx = np.where(Lx==0.0, 1e-54, Lx) # To avoid zero luminosity and division by zero
        
        Edd_state = np.where(ratio <= 1.0, 'Sub-Eddington', 'Super-Eddington')
        beaming = b

    else: # GRRMHD grid for RLO systems

        # RLO (NS + BH fallback, BH overwritten below by GRRMHD prescription)
        Mdot = np.where(mask_RLO, 10**acc_log_mdot * const.Msun / const.secyer, Mdot)

        ratio[mask_valid] = Mdot[mask_valid] / Medd[mask_valid]

        is_Super_Edd = ratio >1.0
        is_RLO_BH_Sup_Edd = mask_RLO & is_BH & is_Super_Edd # For GRRMHD systems RLO-BH systems only

        # For RLO-NS systems, Wind-BH/NS systems
        Lx = eta * Mdot * const.clight**2

        if np.any(is_RLO_BH_Sup_Edd):
            Lx_RLO_BH, b_RLO_BH, W_x_b = Luminosity_GRRMHD_systems(acc_log_mdot[is_RLO_BH_Sup_Edd],
                                                                   spin[is_RLO_BH_Sup_Edd],
                                                                   Medd[is_RLO_BH_Sup_Edd], 
                                                                   L_edd[is_RLO_BH_Sup_Edd],
                                                                   **GRRMHD_kwargs)            
            Weight_times_beaming[is_RLO_BH_Sup_Edd] = W_x_b          
            b[is_RLO_BH_Sup_Edd] = np.asarray(b_RLO_BH)            
            Lx[is_RLO_BH_Sup_Edd] = np.asarray(Lx_RLO_BH)                
            
        Lx = np.where(Lx==0.0, 1e-54, Lx) # To avoid zero luminosity and division by zero
        
        Edd_state = np.where(Mdot / np.where(Medd > 0, Medd, np.inf) <= 1.0, 
                             'Sub-Eddington', 'Super-Eddington')
        beaming = b
        

    return Lx, Edd_state, beaming, Weight_times_beaming

# Selection function to be passed to create_transient_population method
def XRB_selection_function(history_chunk, oneline_chunk, formation_channels_chunk, **kwargs):
    """ 
    A X-ray binary selection function to create a population of XRBs, where we store some important information (Only valid for **Population class** object).

    This function identifies X-ray binaries as an accretor and a donor star and returns their properties.
    Requires to have been through the System_selection_for_XRB() function to only have BH-star or NS-star systems.
    The function will also compute the X-ray luminosity of each selected binary system using the **L_bolo()** function 
    and **Be_XRB_Lx()** function for Be systems. L_bolo() function further splits the calculation according to the grid used for the simulations (classical or GRRMHD) 
    and the binary state (RLO or detached).

    This function has to be passed to **Population_name.create_transient_population(XRB_selection_function,'Transient_population_name')** method 
    to create a population of XRBs from a binary population synthesis simulation.

    
    Parameters
    ----------
    history_chunk : pd.DataFrame, optional
        A chunk of the history DataFrame from population class binaries. 
        Specific lines need to be provided. 
    oneline_chunk : pd.DataFrame
        A chunk of the oneline DataFrame from population class binaries.
    formation_channels_chunk : pd.DataFrame, optional
        A chunk of the formation channels DataFrame of the binaries. Default is None.
    **kwargs : dict, optional
        Additional arguments passed :
        - alpha_BHL (float, default = 1.5) : efficiency parameter for Bondi-Hoyle-Lyttleton accretion.
        - Wind_velo_scheme (str, default = 'Hurley+2002') : Wind velocity scheme to use in Wind_velocity_BHL function. 
        Can be either 'Hurley+2002' _[[4]] or 'Kudritzki+2000' _[[5]].
        - Wind_Disk_BH_formation (str, default = None) : Whether to check for disk formation criteria around BH accretor _[[2]] _[[3]].
        - used_grid (str, default = 'classical') : Whether simulations have been performed using GRRMHD grids for RLO-BH systems or with classical grids. 
        It can be either 'GRRMHD' or 'classical'.
        - N_GRRMHD_stack (int, default = 1) : Number of GRRMHD systems to stack in order to have a smoother interpolation of the beaming factor as a function of the spin of the black hole.
        This is only used if `used_grid` is set to 'GRRMHD'. Only black holes systems accreting via Roche lobe overflow are considered. 
        This stack is useful for these systems because we select a random viewing angle to compute the beaming factor, 
        and therefore by stacking several systems we have a smoother distribution of the beaming factor as a function of the spin of the black hole.

    Returns
    -------
    DF_XRBs : pd.DataFrame
        A DataFrame containing the selected XRBs with accretor/donor data and X-ray luminosities associated to each binary.
    
    References
    ----------
    .. [1] Misra, D., Kovlakas, K., Fragos, T., Lazzarini, M., Bavera, S. S.,  
           Lehmer, B. D., Zezas, A., Zapartas, E., Xing, Z., Andrews, J. J., Dotter, A.,  
           Rocha, K. A., Srivastava, P. M., Sun, M. (2023).  
           *X-ray luminosity function of high-mass X-ray binaries: Studying the signatures of different physical processes using detailed binary evolution calculations*.  
           https://arxiv.org/abs/2209.05505v2
    .. [2] Sen, K. ,Xu, X. -T., Langer, N., El Mellah, I. , Schurmann, C., Quast, M., A&A 652, A138 (2021)
           https://doi.org/10.1051/0004-6361/202141214
           Section 2.3.2, Equation 10
    .. [3] Hirai, R. & Mandel, I. (2021).  
           *Conditions for accretion disc formation and observability of wind-accreting X-ray binaries.*  
           Publications of the Astronomical Society of Australia, 38, e056.   
           https://doi.org/10.1017/pasa.2021.53
    .. [4] Hurley, J. R., Tout, C. A., & Pols, O. R. 2002, MNRAS, 329, 897
    .. [5] Kudritzki, R.-P., & Puls, J. 2000, ARA&A, 38, 613

    Notes
    -----
    - The function is designed to be used with the Population class in posydon, and it assumes to contains the necessary 
    columns in the oneline DataFrame to identify XRBs and compute their properties. The required columns include, but are 
    not limited to : 
    - state_f, S1_state_f. S2_state_f, 
    - S1_mass_f, S2_mass_f, S1_log_R_f, S2_log_R_f, S1_log_L_f, S2_log_L_f, 
    - S1_surface_h1_f, S2_surface_h1_f, S1_He_core_mass_f, S2_He_core_mass_f,
    - S1_surf_avg_omega_div_omega_crit_f, S2_surf_avg_omega_div_omega_crit_f,
    - S1_spin_f, S2_spin_f, S1_lg_mdot_f, S2_lg_mdot_f, S1_log_wind_mdot_f, S2_log_wind_mdot_f,
    - lg_mtransfer_rate_f, orbital_period_f, time_f, metallicity


    Examples
    --------
    >>> from posydon.popsyn.synthetic_population import Population
    
    >>> Systems_selection_XRB(path_to_h5_file/File.h5, 'FileName')

    >>> selec_pop = Population('FileName.h5')
    >>> Transient_pop = selec_pop.create_transient_population(XRB_selection_function, 'Transient_Name')

    # Or, if you need to use kwargs in the selection function : 
    >>> Transient_pop = selec_pop.create_transient_population(lambda h, o, f: XRB_selection_function(h, o, f, **kwargs), 'Transient_Name')

    """
    #===== Control and copy of the dataframes =====
    if not isinstance(oneline_chunk, pd.DataFrame):
        oneline_chunk = oneline_chunk[:]
    else:
        oneline_chunk = oneline_chunk.copy()

    if history_chunk is not None:
        if not isinstance(history_chunk, pd.DataFrame):
            history_chunk = history_chunk[:]
        else:        
            history_chunk = history_chunk.copy()

    if formation_channels_chunk is not None:
        if not isinstance(formation_channels_chunk, pd.DataFrame):
            formation_channels_chunk = formation_channels_chunk[:]
        else:
            formation_channels_chunk = formation_channels_chunk.copy()

    #===== Kwargs control =====
    allowed_keys = {'alpha_BHL', 'Wind_velo_scheme', 'Wind_Disk_BH_formation', 'used_grid', 'spin_split', 'N_GRRMHD_stack'}
    unknown = sorted(set(kwargs.keys()) - allowed_keys)
    if unknown:
        Pwarn(f"XRB_selection_function: unknown parameter(s): {unknown}."
              f"These will be passed to downstream functions.", 
              category="InappropriateValueWarning")
        
    N_GRRMHD_stack = kwargs.get('N_GRRMHD_stack', 1)
    used_grid = kwargs.get('used_grid', 'classical')
    kwargs = {k: v for k, v in kwargs.items() if k not in ['used_grid', 'N_GRRMHD_stack']}

    #===== Filtering Be-systems =====
    orbital_separation = orbital_separation_from_period(oneline_chunk['orbital_period_f'], oneline_chunk['S1_mass_f'], oneline_chunk['S2_mass_f']) # [R⊙]

    ecc = oneline_chunk['eccentricity_f'].to_numpy()
    Donor_RL_radius_periastron = roche_lobe_radius(oneline_chunk['S2_mass_f'].to_numpy(), oneline_chunk['S1_mass_f'].to_numpy(), orbital_separation * (1 - ecc)) # [R⊙]
    Donor_RL_periastron_series = pd.Series(Donor_RL_radius_periastron, index=oneline_chunk.index)

    mask_be = ((oneline_chunk['S2_state_f'] == 'H-rich_Core_H_burning')
                & (oneline_chunk['S2_mass_f'] >= 6.0)
                & (oneline_chunk['state_f'] == 'detached')
                & (oneline_chunk['S2_surf_avg_omega_div_omega_crit_f'] >= 0.7)
                & (oneline_chunk['orbital_period_f'] >= 10.0)
                & (oneline_chunk['orbital_period_f'] <= 300.0)
                & (Donor_RL_radius_periastron <= (100.0 * 10**(oneline_chunk['S2_log_R_f'])))
                )

    Df_Be = oneline_chunk[mask_be].copy()
    Df_Be['Binary_state'] = 'Be system'
    Df_Be['System_type'] = 'Be'
    Df_Be['metallicity'] = oneline_chunk.loc[mask_be, 'metallicity']
    Df_Be['time_f'] = oneline_chunk.loc[mask_be, 'time_f']
    Df_Be['Lx'] = Be_XRB_Lx(oneline_chunk.loc[mask_be, 'orbital_period_f'].to_numpy())

    if formation_channels_chunk is not None:
        Df_Be['formation_channel'] = formation_channels_chunk.loc[Df_Be.index, 'channel']

    Df_Be['Donor_Roche_Lobe_radius'] = Donor_RL_periastron_series.loc[Df_Be.index].to_numpy()
    
    Final_state = Df_Be['S1_state_f']
    is_S1_compact = Final_state.isin(['NS', 'BH'])

    Df_Be['Accretor_state'] = np.where(is_S1_compact, Df_Be['S1_state_f'], Df_Be['S2_state_f'])
    Df_Be['Donor_state'] = np.where(is_S1_compact, Df_Be['S2_state_f'], Df_Be['S1_state_f'])
    Df_Be['Weight_x_beaming'] = True
    Df_Be['Beaming'] = 0.1 # Assumed as a duty cycle of 10% for Be systems
    Df_Be['N_GRRMHD_stack'] = 1 # Be systems are not affected by the GRRMHD grid, so we set it to 1 for all of them.

    #===== RLO / Winds systems =====
    Rest_of_systems = oneline_chunk.loc[~mask_be].copy()
    RLO_mask = Rest_of_systems['state_f'].isin(['RLO1', 'RLO2'])
    RLO_systems = Rest_of_systems.loc[RLO_mask].copy()
    Wind_mask = Rest_of_systems['state_f'] == 'detached'
    Wind_systems = Rest_of_systems.loc[Wind_mask].copy()


    #==========================================================
    # Monte-Carlo stacking of GRRMHD RLO-BH systems
    #==========================================================

    if used_grid == 'GRRMHD' and N_GRRMHD_stack > 1:

        # Identify BH systems among RLO systems
        is_BH_RLO = ((RLO_systems['S1_state_f'] == 'BH') | (RLO_systems['S2_state_f'] == 'BH'))

        # Systems to stack
        RLO_GRRMHD = RLO_systems.loc[is_BH_RLO].copy()

        # NS-RLO systems
        RLO_normal = RLO_systems.loc[~is_BH_RLO].copy()

        # Duplicate GRRMHD systems
        RLO_GRRMHD_stack = pd.concat(
            [RLO_GRRMHD.copy() for _ in range(N_GRRMHD_stack)],
            ignore_index=False)

        RLO_GRRMHD_stack['stack_factor'] = N_GRRMHD_stack

        # Normal systems, NS-RLO
        RLO_normal['stack_factor'] = 1

        # Final RLO dataframe
        RLO_systems = pd.concat(
            [RLO_normal, RLO_GRRMHD_stack],
            ignore_index=False)

    else:

        RLO_systems['stack_factor'] = 1

    Wind_systems['stack_factor'] = 1

    Donor_RL_radius = roche_lobe_radius(oneline_chunk['S2_mass_f'].to_numpy(), oneline_chunk['S1_mass_f'].to_numpy(), orbital_separation) # [R⊙]
    Donor_RL_series = pd.Series(Donor_RL_radius, index=oneline_chunk.index)
    Rest_Donor_RL_radius = Donor_RL_series.loc[~mask_be]

    rows = []

    for name, df in {'RLO' : RLO_systems, 'Wind' : Wind_systems}.items():

        final_state = df['S1_state_f']
        is_S1_compact = final_state.isin(['NS', 'BH'])

        Macc = np.where(is_S1_compact, df['S1_mass_f'], df['S2_mass_f'])
        Mdonor = np.where(is_S1_compact, df['S2_mass_f'], df['S1_mass_f'])
        Acc_log_radius = np.where(is_S1_compact, df['S1_log_R_f'], df['S2_log_R_f'])
        Donor_log_radius = np.where(is_S1_compact, df['S2_log_R_f'], df['S1_log_R_f'])
        acc_state = np.where(is_S1_compact, df['S1_state_f'], df['S2_state_f'])
        donor_state = np.where(is_S1_compact, df['S2_state_f'], df['S1_state_f'])
        acc_log_mdot = np.where(is_S1_compact, df['S1_lg_mdot_f'], df['S2_lg_mdot_f'])
        donor_log_wind_mdot = np.where(is_S1_compact, df['S2_lg_wind_mdot_f'], df['S1_lg_wind_mdot_f'])
        Donor_surface_H1 = np.where(is_S1_compact, df['S2_surface_h1_f'], df['S1_surface_h1_f'])
        donor_He_core_mass = np.where(is_S1_compact, df['S2_he_core_mass_f'], df['S1_he_core_mass_f'])
        spin_acc = np.where(is_S1_compact, df['S1_spin_f'], df['S2_spin_f'])
        Donor_log_lum = np.where(is_S1_compact, df['S2_log_L_f'], df['S1_log_L_f'])

        bin_state = df['state_f']
        lg_mtransfer_rate = df['lg_mtransfer_rate_f']
        a = orbital_separation_from_period(df['orbital_period_f'], df['S1_mass_f'], df['S2_mass_f']) # [R⊙]
        Donor_RL_radius = Rest_Donor_RL_radius.loc[df.index].to_numpy()
        ecc = df['eccentricity_f']

        # If the magnteic field of the NS is available, we can compute the magnetosphere radius and thus the accretion efficiency more accurately. 
        # If not, we assume that the accretion radius is the NS radius. For the moment, magnetic field is not available in simulations, so we set it to None.
        # But if it becomes available, it can be added to the input dataframe and this part of the code can be updated.
        B_surface = None 

        Lx, Edd_state, beaming, Weight_times_beaming = L_bolo(
            used_grid, bin_state, lg_mtransfer_rate, acc_log_mdot, donor_log_wind_mdot, acc_state,
            Donor_log_lum, Macc, Mdonor, Acc_log_radius, Donor_log_radius,
            Donor_surface_H1, donor_He_core_mass, Donor_RL_radius, a, ecc, spin_acc, B_surface,
            **kwargs
        )


        out = pd.DataFrame(index=df.index)
        out['metallicity'] = df['metallicity']
        out['time'] = df['time_f']
        out['Binary_state'] = bin_state
        out['Accretor_state'] = acc_state
        out['Accretor_spin'] = spin_acc
        out['Accretor_mass'] = Macc
        out['Accretor_LogR'] = Acc_log_radius
        out['Donor_state'] = donor_state
        out['Donor_mass'] = Mdonor
        out['Donor_LogR'] = Donor_log_radius
        out['Donor_Roche_lobe_radius'] = Donor_RL_radius
        out['Donor_LogL'] = Donor_log_lum
        out['Donor_H1_at_surface'] = Donor_surface_H1
        out['Donor_He_Core_mass'] = donor_He_core_mass
        out['Accretor_LogMdot'] = acc_log_mdot
        out['Donor_LogMdot_wind'] = donor_log_wind_mdot
        out['lg_mtransfer_rate'] = lg_mtransfer_rate
        out['Orbital_separation'] = a
        out['Eccentricity'] = ecc
        out['System_type'] = name
        out['Eddington_state'] = Edd_state
        out['Beaming'] = beaming
        out['Weight_x_beaming'] = Weight_times_beaming
        out['stack_factor'] = df['stack_factor']
        out['Lx'] = Lx
        if formation_channels_chunk is not None:
            out['formation_channel'] = formation_channels_chunk.loc[df.index, 'channel']

        rows.append(out)

    final_df = pd.concat(rows + [Df_Be], axis=0, ignore_index=False)

    keep_cols = ['metallicity', 'time', 'lg_mtransfer_rate',
    'Binary_state','Accretor_state','Accretor_spin','Accretor_mass','Accretor_LogR',
    'Donor_state','Donor_mass','Donor_LogR','Donor_Roche_lobe_radius','Donor_LogL',
    'Donor_H1_at_surface','Donor_He_Core_mass','Accretor_LogMdot','Donor_LogMdot_wind',
    'Orbital_separation','Eccentricity','System_type','Beaming', 'Weight_x_beaming', 'stack_factor', 
    'Eddington_state','Lx', 'formation_channel']

    final_df = final_df[keep_cols].copy()

    return final_df 


#===== Plotting functions -Cumulative XLF ===== 

def Stat_XLF(transient_population, Lmin=1e30, SFR = 1.0, stat_kwargs=None):
    """
    Compute statistics for a given transient population of XRBs, including the weigthed cumulative number of XRBs above a given luminosity.

    Parameters
    ----------
    transient_population : TransientPopulation object
        A transient population of XRBs.
    Lmin : float, optional
        Minimum luminosity threshold for XRBs [erg/s]. Default is 1e30.
    SFR : float, optional
        Star formation rate [M⊙/yr]. Default is 1.0. (E.g. Milky-Way SFR ~ 1.5 M⊙/yr. Starburst galaxies SFR ~ 10 M⊙/yr)
    stat_kwargs : dict, optional
        A dictionary of population parameters to override the default ones used in the transient population in order to weight the population accordingly.
        See Notes for allowed keys. Default is *None*.

    Returns
    -------
    Weights : pd.Series
        Weights associated with each XRB in the population.

    References
    ----------
    .. [1] Misra, D., Kovlakas, K., Fragos, T., Lazzarini, M., Bavera, S. S.,  
           Lehmer, B. D., Zezas, A., Zapartas, E., Xing, Z., Andrews, J. J., Dotter, A.,  
           Rocha, K. A., Srivastava, P. M., Sun, M. (2023).  
           *X-ray luminosity function of high-mass X-ray binaries: Studying the signatures of different physical processes using detailed binary evolution calculations*.  
           https://arxiv.org/abs/2209.05505v2

    Notes
    -----
    The function uses the calculate_model_weights method of the TransientPopulation class to compute weights based on the specified population parameters.
    For that, you can provide a dictionary of population parameters to override the default ones used in the transient population in order to weight the population accordingly.
    Therefore, the allowed keys in the pop_params dictionary are:
    - 'number_of_binaries'
    - 'binary_fraction_scheme', 'binary_fraction_const'
    - 'star_formation'
    - 'max_simulation_time'
    - 'primary_mass_scheme', 'primary_mass_min', 'primary_mass_max'
    - 'secondary_mass_scheme', 'secondary_mass_min', 'secondary_mass_max'
    - 'orbital_scheme', 'orbital_period_scheme', 'orbital_period_min', 'orbital_period_max'
    - 'eccentricity_scheme'
    - 'q_min', 'q_max'
    
    Example
    -------
    >>> stat_kwargs = {'primary_mass_min': 0.01, 'primary_mass_max': 200.0, 'q_min': 0.0, 'q_max': 1.0,}
    >>> Weight, IMF_correction = Stat_XLF(transient_population, stat_kwargs=stat_kwargs)

    """

    #===== Statistics Kwargs =====

    # Keys allowed in calculate_model_weights()
    from posydon.popsyn.binarypopulation import saved_ini_parameters

    _allowed_keys = set(saved_ini_parameters) | {"q_min", "q_max"}

    Pop_params = transient_population.ini_params.copy()

    if stat_kwargs is not None:
        if not isinstance(stat_kwargs, dict):
            raise TypeError("stat_kwargs must be a dict or None.")
        stat_keys = set(stat_kwargs)
        if not stat_keys.issubset(_allowed_keys):
            bad = sorted(stat_keys - _allowed_keys)
            raise ValueError(
                f"Stat_XLF: Invalid population parameter(s): {bad}. Allowed keys: {sorted(_allowed_keys)}"
            )
        for key, value in stat_kwargs.items():
            Pop_params[key] = value


    T_window_yr = Pop_params['max_simulation_time'] # Years
    Total_underlying_mass = SFR * T_window_yr  # Total underlying mass in M⊙

    if Total_underlying_mass <= 0.0:
        raise ValueError("plot_ccdf: Total_underlying_mass must be > 0.0")

    base_w = transient_population.calculate_model_weights(model_weights_identifier="IMF",
                                                          population_parameters=Pop_params)
    
    if isinstance(base_w, pd.DataFrame):
        # prefer the named column if present, else take first column
        if "IMF" in base_w.columns:
            weights_series = base_w["IMF"].to_numpy(float)
        else:
            weights_series = base_w.iloc[:, 0].to_numpy(float)
    else:
        weights_series = np.asarray(base_w, float)  # assume Series-like 

    Weights = weights_series * Total_underlying_mass / transient_population.population['stack_factor'].to_numpy(dtype=float)

    return Weights


def plot_ccdf(transient_population, ax, weights, mask, label, SFR = 1.0, **ccdf_kwargs):
    """
    Compute the cumulative X-ray luminosity function (CCDF) of a given transient populations. 
    This function returns a *time-averaged* CCDF, i.e. the expected number of X-ray binaries brighter than L at a random time, 
    assuming constant star formation over a time window T_window_Myr.\n
    Masks split the population into different sub-populations to be plotted separately (e.g. RLO, Wind-fed, Be-XRBs, etc.).

    Parameters
    ----------
    transient_population : TransientPopulation object
        The transient population.
    ax : matplotlib.axes.Axes
        The axes on which to plot the CCDF (_fig , ax = plt.subplots(*args, **kwargs)_).
    weights : array-like
        The weights associated with each value.
    mask : array-like
        A boolean mask to filter the values.
    label : str
        The label for the CCDF plot.
    SFR : float, optional
        Star formation rate [M⊙/yr]. Default is 1.0. (E.g. Milky-Way SFR ~ 1.5 M⊙/yr. Starburst galaxies SFR ~ 10 M⊙/yr)
    linestyle : str
        The line style for the CCDF plot.
    **kwargs : dict, optional
        Additional keyword arguments for customization.
        - norm_factor (float, default=0.5) : Normalization factor L_x/L_bolo ~ 0.5 from Misra et al.[1] page 4
        - extend_to_edges (bool, default=True) : Whether to extend the CCDF to the edges of the plot.
    
        Any other keyword arguments valid for `matplotlib.pyplot` can be passed.

    Returns
    -------
    None

    References
    ----------
    .. [1] Misra, D., Kovlakas, K., Fragos, T., Lazzarini, M., Bavera, S. S.,  
           Lehmer, B. D., Zezas, A., Zapartas, E., Xing, Z., Andrews, J. J., Dotter, A.,  
           Rocha, K. A., Srivastava, P. M., Sun, M. (2023).  
           *X-ray luminosity function of high-mass X-ray binaries: Studying the signatures of different physical processes using detailed binary evolution calculations*.  
           https://arxiv.org/abs/2209.05505v2
    Notes
    -----
    - Weights are assumed to already include IMF, binary fraction, and population 
    synthesis corrections

    Example
    -------
    >>> fig, ax = plt.subplots()
    >>> plot_ccdf(transient_population, ax, weights, mask, label='Label name', color='blue', linestyle='-')
    """

    #===== kwargs =====
    controlled_keys = {'norm_factor', 'extend_to_edges'}

    plot_keys = dict(label=label)
    extra_plot_keys = {key: var for key, var in ccdf_kwargs.items() if key not in controlled_keys}
    plot_keys.update(extra_plot_keys)
 
    extend_to_edges = ccdf_kwargs.get('extend_to_edges', True) # Whether to extend the CCDF to the edges of the plot
    norm_factor = ccdf_kwargs.get('norm_factor', 0.5) # Normalization factor L_x/L_bolo ~ 0.5 from Misra et al.[1] page 4

    if norm_factor <= 0.0:
        raise ValueError("plot_ccdf: norm_factor must be > 0.0")

    #===== Compute the CCDF =====
    L = np.log10(transient_population.population['Lx'].to_numpy(dtype = float) * norm_factor)

    W_x_b = transient_population.population['Weight_times_beaming'].to_numpy(dtype=bool)
    b = transient_population.population['Beaming'].to_numpy(dtype = float)
    b = np.where(np.isfinite(b) & (b >0), b, 1.0)

    selection = mask & np.isfinite(L) & (L >= 35)

    if not np.any(selection):
        return

    L_selec = L[selection]

    # If we need to multiply by the beaming factor, we do it here, Otherwise, 
    # we assume that the weight already includes beaming (GRRMHD BH-RLO systems only).
    Weight_select = weights[selection] * np.where(W_x_b[selection], b[selection], 1.0)

    order = np.argsort(L_selec)
    L_sorted = L_selec[order]
    Weights_sorted = Weight_select[order]

    N = np.cumsum(Weights_sorted[::-1])[::-1]  # Cumulative sum from the end to the beginning

    #===== Plotting the CCDF =====
    # Extension of the CCDF to the edges of the plot, with a horizontal line at y = N(Lmin) for L < Lmin, 
    # and a vertical line at y = 0 for L > Lmax 
    # (if the xlim of the plot is larger than the range of L values in the population).
    if extend_to_edges:
        cur_xlim = ax.get_xlim()
        if cur_xlim == (0.0, 1.0):
            xmin, xmax = float(L_sorted[0]), float(L_sorted[-1])
        else:
            xmin, xmax = cur_xlim
        ymin = np.nextafter(0, 1)
        
    # Total vlaue above the smallest L
        y0 = float(np.sum(Weights_sorted))

        x_ext = np.r_[xmin, L_sorted, xmax]
        y_ext = np.r_[y0, N, ymin]

        ax.plot(x_ext, y_ext, **plot_keys)
    else:
        ax.plot(L_sorted, N, **plot_keys)


def Plot_cumulative_XLF(transient_population, weights, SFR = 1.0, title=None, ax=None, figsize=None, show=False, **plot_kwargs):
    '''
    Plot the cumulative X-ray luminosity function (CCDF) of a given transient population of XRBs.

    Parameters
    ----------
    transient_population : TransientPopulation object
        The transient population.
    weights : array-like
        The weights associated with each value.
    SFR : float, optional
        Star formation rate [M⊙/yr]. Default is 1.0. (E.g. Milky-Way SFR ~ 1.5 M⊙/yr. Starburst galaxies SFR ~ 10 M⊙/yr)
    title : str, optional
        The title of the plot. If None, no title will be set. Default is None.
    ax : matplotlib.axes.Axes, optional
        The axes on which to plot the CCDF. If None, a new figure and axes will be created. Default is None. show must be false if ax is not None. 
    figsize : tuple, optional
        The size of the figure to create if ax is None. Default is (6.3, 6.3 * 0.618), which corresponds to a golden ratio aspect.
    show : bool, optional
        Whether to display the plot. Default is False. ax must be None for show=True, otherwise a warning will be raised and the plot will not be shown.
    **plot_kwargs : dict, optional
        Additional keyword arguments for customization of the plot.
        See Notes for allowed keys.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the plot.
    ax : matplotlib.axes.Axes
        The axes object containing the plot.

    Notes
    -----
    Allowed keys in plot_kwargs:
    - norm_factor (float, default=0.5) : Normalization factor L_x/L_bolo ~ 0.5 from Misra et al.[1] page 4
    - extend_to_edges (bool, default=True) : Whether to extend the CCDF to the edges of the plot.
    
    Any other keyword arguments valid for `matplotlib.pyplot` can be passed and will be forwarded to the plot_ccdf function.
    '''

    figsize = (6.3, 6.3 * 0.618) if figsize is None else figsize

    if ax is None:
        fig , ax = plt.subplots(figsize=figsize, dpi=300)
    else: 
        fig = ax.figure
    
    met = transient_population.population['metallicity'][0] if 'metallicity' in transient_population.population.columns else None
    title = f"{title} (Z={met:.2e} $Z_{{\odot}}$)" if title and met is not None else title


    Lumino_range = (np.log10(transient_population.population['Lx'] * 0.5) >= 35) & (np.log10(transient_population.population['Lx'] * 0.5) <= 43) #1e35) #1e43)
    RLO = transient_population.population['System_type'] == 'RLO'
    Be_XRB = transient_population.population['Binary_state'] == 'Be system' 
    Wind = transient_population.population['System_type'] == 'Wind'

    BH = transient_population.population['Accretor_state'] == 'BH'
    NS = transient_population.population['Accretor_state'] == 'NS'

    H1_rich = transient_population.population['Donor_state'].str.contains('H-rich|accreted_He', regex=True, na=False)
    He_rich = transient_population.population['Donor_state'].str.contains('stripped_He|accreted_He_Core_He_burning', regex=True, na=False)

    # Default colors to directly compare with Figure 1 of Misra et al. (A&A 672, A99 (2023)) 
    plot_ccdf(transient_population, ax, weights, Lumino_range, 'Total XRB population', SFR = SFR, linestyle='-', color='#1F77B4', **plot_kwargs)
    plot_ccdf(transient_population, ax, weights, Lumino_range & Be_XRB, 'Be-XRBs', SFR = SFR, linestyle=':', color='#808080',**plot_kwargs)
    plot_ccdf(transient_population, ax, weights, Lumino_range & RLO & BH & He_rich, 'RLO BH-He-rich', SFR = SFR, linestyle='-', color='#B760C4', **plot_kwargs)
    plot_ccdf(transient_population, ax, weights, Lumino_range & RLO & BH & H1_rich, 'RLO BH-H-rich', SFR = SFR, linestyle='-', color='#540B0E', **plot_kwargs)
    plot_ccdf(transient_population, ax, weights, Lumino_range & RLO & NS & He_rich, 'RLO NS-He-rich', SFR = SFR, linestyle='-', color='#F9E979', **plot_kwargs)
    plot_ccdf(transient_population, ax, weights, Lumino_range & RLO & NS & H1_rich, 'RLO NS-H-rich', SFR = SFR, linestyle='-', color='#DA1E37', **plot_kwargs)
    plot_ccdf(transient_population, ax, weights, Lumino_range & Wind & BH & He_rich, 'Wind BH-He-rich', SFR = SFR, linestyle='--', color='#B760C4', **plot_kwargs)
    plot_ccdf(transient_population, ax, weights, Lumino_range & Wind & BH & H1_rich, 'Wind BH-H-rich', SFR = SFR, linestyle='--', color='#540B0E', **plot_kwargs)
    plot_ccdf(transient_population, ax, weights, Lumino_range & Wind & NS & He_rich, 'Wind NS-He-rich', SFR = SFR, linestyle='--', color='#F9E979', **plot_kwargs)
    plot_ccdf(transient_population, ax, weights, Lumino_range & Wind & NS & H1_rich , 'Wind NS-H-rich', SFR = SFR, linestyle='--', color='#DA1E37', **plot_kwargs)


    ax.set_title(title if title else r'XLF of XRBs population', pad = 100, fontsize=11)

    ax.set_xlabel(r'$\log_{10}\left(L_{\mathrm{X},\,0.5-8\,\mathrm{keV}} \,/\, \mathrm{erg\,s^{-1}}\right)$', fontsize=11)
    ax.set_xlim(35, 43)
    ax.set_xticks(np.arange(35, 43, 1))

    ax.set_ylabel(r'$N(> L_{\mathrm{X},\,0.5-8\,\mathrm{keV}})\,/\,(\mathrm{M}_\odot\,\mathrm{yr}^{-1})$', fontsize=11)
    ax.set_yscale('log')
    ax.set_ylim(2e-4, 4e2)


    ax.tick_params(axis='both', which='both', labelsize=9, top=False, right=False)

    ax.legend(loc='lower center', 
              bbox_to_anchor=(0.5, 1.02), 
              bbox_transform=ax.transAxes, 
              ncol=3, 
              frameon=True, 
              fancybox=True, 
              edgecolor='black', 
              handlelength=1.5, 
              borderpad=0.4, 
              fontsize=9
              )
    
    fig.subplots_adjust(top=0.8)  # ajuster l'espace supérieur pour la légende

    if show:
        plt.show()
    
    return fig, ax

def plot_Lehmer19_XLF(ax, path_fits, SFR_total=45.4, color='#FFB2B2', label='Lehmer+19', **kwargs):
    """
    Superpose the observed XLF from Lehmer+19 (table7) on an existing axis. This function uses the total SFR of the observed population of 38 galaxies from `Lehmer+19` to normalize the XLF.
    The total SFR is 45.4 M⊙/yr, which is the sum of the SFRs of the 38 galaxies in the sample (see Lehmer+19, section 3.2, Table 1).
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis on which to plot.
    path_fits : str
        Path to the J_ApJS_243_3_table7.dat.fits file.
    SFR_total : float
        Total SFR of the galaxies from Lehmer+19 in [M⊙/yr].
        Typical value from the paper : ~45.4 M⊙/yr (sum of SFRs of the 38 galaxies).
    color : str
        Color of the curve. Default : '#FFB2B2' (light pink).
    label : str
        Label of the curve. Default : 'Lehmer+19'.
    **kwargs : dict
        Additional keyword arguments for customization of the plot.
    
    Returns
    -------
    None

    References
    ----------
    .. [1] Bret D. Lehmer et al 2019 ApJS 243 3
           *X-Ray Binary Luminosity Function Scaling Relations for Local Galaxies Based on Subgalactic Modeling*
           DOI 10.3847/1538-4365/ab22a8

    Notes
    -----
    The function filters the sources classified as "Loc == 1" (see Table 7 and Appendix "X-Ray Point-source Catalog" of `Lehmer+19` for more details on the classification),
    and plots the cumulative XLF N(>L) normalized by the total SFR of the sample (Table 1).

    Examples
    --------
    >>> fig, ax = plt.subplots()
    >>> plot_Lehmer19_XLF(ax, 'path/to/J_ApJS_243_3_table7.dat.fits', SFR_total=45.4, label='Lehmer+19')

    Combining it with Plot_cumulative_XLF: 
    >>> fig, ax = Plot_cumulative_XLF(transient_population, weights, title='XLF of XRBs population')
    >>> plot_Lehmer19_XLF(ax, 'path/to/J_ApJS_243_3_table7.dat.fits', SFR_total=45.4, label='Lehmer+19')
    >>> ax.set_title('XLF of XRBs population with Lehmer+19 data')
    >>> plt.show()
    """
    from astropy.io import fits
    from astropy.stats import poisson_conf_interval

    # Fits loading
    hdul = fits.open(path_fits)
    data = hdul[1].data

    # Filter only loc == 1, this only contains systems used to construct 
    # the XLF in Lehmer+19 and contains only systems with ≥ 1e35 erg/s.
    mask = data['Loc'] == 1

    sources = data[mask]
    logL = sources['logLFB']  # log10(L [erg/s])
    valid = np.isfinite(logL) & (logL >= 35)

    logL_valid = np.sort(logL[valid])

    counts = np.arange(len(logL_valid), 0, -1)

    low, high = poisson_conf_interval(counts, interval='frequentist-confidence', sigma=1)


    # XLF cumulative N(>L) normalized by the total SFR of the sample
    N_low = low / SFR_total
    N_high = high / SFR_total

    ymin = ax.get_ylim()[0]  
    xmin = ax.get_xlim()[0]

    # To extend the curve to the edges of the plot.
    x_ext = np.r_[xmin, logL_valid, logL_valid[-1]]
    y_low_ext = np.r_[N_low[0], N_low, ymin]
    y_high_ext = np.r_[N_high[0], N_high, ymin]

    fill_lw = kwargs.pop('linewidth', 0.8)
    ax.fill_between(x_ext, y_low_ext, y_high_ext, color=color, label=label, linewidth=fill_lw, **kwargs)

    ax.fill_betweenx([ymin, N_high[-1]], 
                     logL_valid[-1] - 0.15,
                     logL_valid[-1] + 0.02,
                     color=color)
    
    hdul.close()


#===== Other functions used for my work =====

def find_slope_XLF(Transient_population, verbose=False):
    """
    This function calculates the slope of the X-ray luminosity function of black holes accreting via Roche lobe overflow at super-Eddington rates.
    It follows a power-law in log scale.

    Parameters
    ----------
    Transient_population : DataFrame
            The population of transient sources, containing at least the columns 'System_type', 'Accretor_state', 'Eddington_state', and 'Lx'.

    Returns
    -------
    slope : float
            The slope of the cumulative X-ray luminosity function (N(>L) ~ L^slope).
    intercept : float
            The intercept of the linear fit in log-log space.
    alpha_diff : float
            The differential slope of the X-ray luminosity function (dN/dL ~ L^-alpha_diff), where alpha_diff = 1 - slope.
    """

    df = Transient_population.population.copy()

    # Mask RLO-BH super-Eddington
    mask = ((df['System_type'] == 'RLO') & 
           (df['Accretor_state'] == 'BH') &
           (df['Eddington_state'] == 'Super-Eddington')
           )
    
    Lx = df.loc[mask, 'Lx'].to_numpy(dtype=float)
    Lx = Lx[Lx > 0]

    if len(Lx) < 2: 
        return None, None, None
    # CCDF
    L_sorted = np.sort(Lx)
    N_cum = np.arange(len(L_sorted), 0, -1)

    # Fit range
    sel = (L_sorted >= 1e39) & (L_sorted <= 1e44)

    if sel.sum() < 2:
        print('Not enough points in the range 1e39-1e44 to perform a fit.')
        return None, None, None
    
    logL = np.log10(L_sorted[sel])
    logN = np.log10(N_cum[sel])

    # Linear fit in log-log
    slope, intercept = np.polyfit(logL, logN, 1)

    #  dN/dL ~ L^{-alpha}
    alpha_diff = 1.0 - slope

    if verbose:
        print(f"Linear fit (N(>L) vs L): slope = {slope:.4f}, intercept = {intercept:.4f}")
        print(f"Differential index (alpha) of the power law (dN/dL ~ L^-alpha) = {alpha_diff:.4f}")

    return slope, intercept, alpha_diff


#===== Functions improved from common_functions.py =====

# For functions that had a "todo" comment in common_functions.py.

def Orbital_separation_from_period(Mass1, Mass2, period_days):
    """
    Calculate the orbital separation "a" given the period of the system from the 3rd Kepler's law. 

    Parameters 
    ---------
    Mass1 : float 
        Mass of the central body/accretor where Mass1 > Mass2 in [M⊙]
    Mass2 : float
        Mass of the orbiting object/donor where Mass2 < Mass1 in [M⊙]
    period_days : float 
        Period of the orbiting companion in [days]
    
    Returns
    -------
    a_cm : float
        Orbital separation of the binary system in [Rsun]
    """
    if Mass1 is None or Mass2 is None or period_days is None:
        Pwarn("Mass1, Mass2, and period_days must not be None.", 'InappropriateValueWarning')
        return np.nan, np.nan

    #===== Control parameters =====
    mass1 = np.asarray(Mass1, dtype=float)
    mass2 = np.asarray(Mass2, dtype=float)
    period_days = np.asarray(period_days, dtype=float)

    if np.any(mass1 < 0.0) or np.any(mass2 < 0.0):
        Pwarn("Mass1 and Mass2 must be non-negative values.",'InappropriateValueWarning')
        return np.full_like(period_days, np.nan, dtype=float), np.full_like(period_days, np.nan, dtype=float)

    if np.any(period_days <= 0.0) or np.any(np.isnan(period_days)):
        Pwarn("Orbital period must be a positive value.",'InappropriateValueWarning')
        return np.full_like(period_days, np.nan), np.full_like(period_days, np.nan)

    #===== Calculation =====
    period_seconds = period_days * const.day2sec
    M1_gr = mass1 * const.Msun
    M2_gr = mass2 * const.Msun
    a_cm = (const.standard_cgrav * (M1_gr + M2_gr) * period_seconds**2 /(4 * np.pi**2))**(1/3) # [cm]
    a_Rsun = a_cm / const.Rsun  # [Rsun]

    return a_Rsun


def Orbital_velocity_from_separation_and_period(separation, period_days):
    """
    Calculate the orbital velocity of a two body system given the orbital separation and period.

    Parameters
    ----------
    separation : float or array of float
        Orbital separation of a two body system in [R⊙]. 
    period_days : float or array of float
        Orbital period of the two body system in [days].

    Returns
    -------
    V_orb : float or array of float
        Orbital velocity of the two body system in [km/s].
    """
    
    separation = np.asarray(separation, dtype=float)
    period_days = np.asarray(period_days, dtype=float)

    if separation is None or np.any(separation <= 0.0) or np.any(np.isnan(separation)):
        Pwarn("separation must not be None and must be a positive value.", 
              'InappropriateValueWarning')
        return np.nan
    
    if period_days is None or np.any(period_days <= 0.0) or np.any(np.isnan(period_days)):
        Pwarn("period_days must not be None and must be a positive value.", 
              'InappropriateValueWarning')
        return np.nan

    V_orb = 2 * np.pi * (separation * const.Rsun) / (period_days * const.day2sec)  / 1e5  # [km/s]

    return V_orb


def Orbital_period_from_separation(separation, Mass1, Mass2):
    """

    Apply the Third Kepler law.

    Parameters
    ----------
    separation : float
        Orbital separation (semi-major axis) in Rsun.
    Mass1 : float
        Mass of one of the stars in solar units.
    Mass2 : float
        Mass of the other star in solar units.

    Returns
    -------
    P_days : float
        The orbital period in days.
    """

    if Mass1 is None or Mass2 is None or separation is None:
        Pwarn("Mass1, Mass2, and separation must not be None.", 'InappropriateValueWarning')
        return np.nan, np.nan

    mass1 = np.asarray(Mass1, dtype=float)
    mass2 = np.asarray(Mass2, dtype=float)
    sep = np.asarray(separation, dtype=float)

    if np.any(mass1 < 0.0) or np.any(mass2 < 0.0):
        Pwarn("Mass1 and Mass2 must be non-negative values.",'InappropriateValueWarning')
        return np.full_like(sep, np.nan, dtype=float)

    if np.any(separation < 0.0) or np.any(np.isnan(separation)):
        Pwarn("Orbital separation must be a non-negative value.",'InappropriateValueWarning')
        return np.full_like(sep, np.nan, dtype=float)
    
    P = 2 * np.pi * np.sqrt((sep * const.Rsun)**3 / (const.standard_cgrav * (mass1 + mass2) * const.Msun)) # [seconds]
    P_days = P / const.day2sec  # [days]

    return P_days
