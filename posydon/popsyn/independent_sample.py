"""Generate the initial parameters for a binary population."""


__authors__ = [
    "Devina Misra <devina.misra@unige.ch>",
    "Jeffrey Andrews <jeffrey.andrews@ufl.edu>",
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Scott Coughlin <scottcoughlin2014@u.northwestern.edu>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]


import numpy as np
from scipy.stats import truncnorm

from posydon.popsyn import IMFs, distributions
from posydon.popsyn.Moes_distributions import Moe_17_PsandQs
from posydon.utils.common_functions import rejection_sampler

_gen_Moe_17_PsandQs = None


def generate_independent_samples(orbital_scheme='period', **kwargs):
    """Randomly generate a population of binaries at ZAMS.

    Parameters
    ----------
    orbital_scheme : str (default: 'period')
        The scheme to use to get either orbital periods or separations
    **kwargs : dictionary
        kwargs from BinaryPopulation class

    Returns
    -------
    orbital_scheme_set : ndarray of floats
        Randomly drawn orbital separations/periods depending on the scheme
    eccentricity_set : ndarray of floats
        Randomly drawn eccentricities
    m1_set : ndarray of floats
        Randomly drawn primary masses
    m2_set : ndarray of floats
        Randomly drawn secondary masses

    """
    global _gen_Moe_17_PsandQs

    primary_mass_min = kwargs.get("primary_mass_min", 7)
    primary_mass_max = kwargs.get("primary_mass_max", 120)
    secondary_mass_min = kwargs.get("secondary_mass_min", 0.35)
    secondary_mass_max = kwargs.get("secondary_mass_max", 120)
    if primary_mass_max < primary_mass_min:
        raise ValueError("primary_mass_max must be larger than primary_mass_min.")
    if secondary_mass_max < secondary_mass_min:
        raise ValueError("secondary_mass_max must be larger than secondary_mass_min.")

    # Generate primary masses
    m1_set = generate_primary_masses(**kwargs)

    if use_Moe_17_PsandQs(orbital_scheme=orbital_scheme, **kwargs):
        # initialize generator for Moe+17-PsandQs
        if _gen_Moe_17_PsandQs is None:
            _gen_Moe_17_PsandQs = Moe_17_PsandQs(**kwargs)
        # use same defaults as generate_primary_masses
        M1_min = kwargs.get("primary_mass_min", 7)
        M1_max = kwargs.get("primary_mass_max", 120)
        # generate samples
        m2_set, orbital_scheme_set, eccentricity_set, metallicity_set\
         = _gen_Moe_17_PsandQs(m1_set, M_min=M1_min, M_max=M1_max,
                               all_binaries=False)

        ecc_scheme = kwargs.get('eccentricity_scheme', 'zero')
        if ecc_scheme != 'Moe+17-PsandQs':
            # ensure that kwargs hold the current retrieved or
            # set value of eccentricity_scheme
            kwargs['eccentricity_scheme'] = ecc_scheme
            eccentricity_set = generate_eccentricities(**kwargs)
    else:

        # Generate secondary masses
        m2_set = generate_secondary_masses(m1_set, **kwargs)

        if orbital_scheme == 'separation':
            # Generate orbital separations
            orbital_scheme_set = generate_orbital_separations(**kwargs)
        elif orbital_scheme == 'period':
            # Generate orbital periods
            orbital_scheme_set = generate_orbital_periods(m1_set, **kwargs)
        else:
            raise ValueError("Allowed orbital schemes are separation or"
                             " period.")

        # Generate eccentricities
        eccentricity_set = generate_eccentricities(**kwargs)

        # Calculate the binary fraction
        RNG = kwargs.get('RNG', np.random.default_rng())
        binary_fraction = generate_binary_fraction(m1=m1_set, **kwargs)

        # Set values to nan for single stars
        idx = np.where(RNG.uniform(size = len(m1_set)) > binary_fraction)[0]
        orbital_scheme_set[idx] = np.nan
        eccentricity_set[idx] = np.nan
        m2_set[idx] = np.nan

    return orbital_scheme_set, eccentricity_set, m1_set, m2_set


def use_Moe_17_PsandQs(secondary_mass_scheme='', orbital_scheme='',
                       orbital_period_scheme='', eccentricity_scheme='',
                       **kwargs):
    """Check whether Moe & Di Stefano (2017) [1]_ should be used for the
    initial sampling.

    References
    ----------
    .. [1] Moe, M. and Di Stefano, R., “Mind Your Ps and Qs: The Interrelation
        between Period (P) and Mass-ratio (Q) Distributions of Binary Stars”,
        <i>The Astrophysical Journal Supplement Series</i>, vol. 230, no. 2,
        Art. no. 15, IOP, 2017. doi:10.3847/1538-4365/aa6fb6.
    """
    return ((secondary_mass_scheme=='Moe+17-PsandQs')
            or ((orbital_scheme=='period')
                and (orbital_period_scheme=='Moe+17-PsandQs'))
            or (eccentricity_scheme=='Moe+17-PsandQs'))


def generate_orbital_periods(primary_masses,
                             number_of_binaries=1,
                             orbital_period_min=0.35,
                             orbital_period_max=10**3.5,
                             orbital_period_scheme='Sana+12_period_extended',
                             **kwargs):
    """Randomaly generate orbital periods for a sample of binaries."""
    RNG = kwargs.get('RNG', np.random.default_rng())

    # Check inputs

    # Sana H., et al., 2012, Science, 337, 444
    if orbital_period_scheme == 'Sana+12_period_extended':
        period_dist = distributions.Sana12Period(
            p_min=orbital_period_min,
            p_max=orbital_period_max
        )
        orbital_periods = period_dist.rvs(size=number_of_binaries, m1=primary_masses, rng=RNG)
    else: # pragma: no cover
        raise ValueError("You must provide an allowed orbital period scheme.")

    return orbital_periods


def generate_orbital_separations(number_of_binaries=1,
                                 orbital_separation_min=5,
                                 orbital_separation_max=1e5,
                                 log_orbital_separation_mean=None,
                                 log_orbital_separation_sigma=None,
                                 orbital_separation_scheme='log_uniform',
                                 **kwargs):
    """Generate random orbital separations.

    Use the scheme defined in this particular instance of BinaryPopulation.

    Parameters
    ----------
    number_of_binaries : int
        Number of binaries that require randomly sampled orbital separations
    orbital_separation_min : float
        Minimum orbital separation in solar radii
    orbital_separation_max : float
        Maximum orbital separation in solar radii
    log_orbital_separation_mean : float
        Mean of the lognormal distribution.
    log_orbital_separation_sigma : float
        Standard deviation of the lorgnormal distribution.
    orbital_separation_scheme : string
        Distribution from which the orbital separations are randomly drawn

    Returns
    -------
    orbital_separations : ndarray of floats
        Randomly drawn orbital separations
    """
    RNG = kwargs.get('RNG', np.random.default_rng())

    orbital_separation_scheme_options = ['log_uniform', 'log_normal']

    # Check inputs
    if orbital_separation_scheme not in orbital_separation_scheme_options:
        raise ValueError("You must provide an allowed "
                         "orbital separation scheme.")

    if orbital_separation_scheme == 'log_uniform':
        sep_dist = distributions.LogUniform(
            min=orbital_separation_min,
            max=orbital_separation_max
        )
        orbital_separations = sep_dist.rvs(size=number_of_binaries, rng=RNG)

    if orbital_separation_max < orbital_separation_min:
        raise ValueError("`orbital_separation_max` must be "
                         "larger than the orbital_separation_min.")

    elif orbital_separation_scheme == 'log_normal':
        if (log_orbital_separation_mean is None
                or log_orbital_separation_sigma is None):
            raise ValueError(
                "For the `log_normal separation` scheme you must give "
                "`log_orbital_separation_mean`, "
                "`log_orbital_separation_sigma`.")

        sep_dist = distributions.LogNormalSeparation(
            mean=log_orbital_separation_mean,
            sigma=log_orbital_separation_sigma,
            min=orbital_separation_min,
            max=orbital_separation_max
        )
        orbital_separations = sep_dist.rvs(size=number_of_binaries, rng=RNG)

    else:  # pragma: no cover
        pass

    return orbital_separations


def generate_eccentricities(number_of_binaries=1,
                            eccentricity_scheme='thermal',
                            **kwargs):
    """Generate random eccentricities.

    Use the scheme defined in this particular instance of BinaryPopulation.

    Parameters
    ----------
    number_of_binaries : int
        Number of binaries that require randomly sampled orbital separations
    eccentricity_scheme : string
        Distribution from which eccentricities are randomly drawn
    **kwargs : dictionary
        kwargs from BinaryPopulation class

    Returns
    -------
    eccentricities : ndarray of floats
        Randomly drawn eccentricities

    """
    RNG = kwargs.get('RNG', np.random.default_rng())

    eccentricity_scheme_options = ['thermal', 'uniform', 'zero']

    if eccentricity_scheme not in eccentricity_scheme_options:
        raise ValueError("You must provide an allowed eccentricity scheme.")

    if eccentricity_scheme == 'thermal':
        ecc_dist = distributions.ThermalEccentricity()
        eccentricities = ecc_dist.rvs(size=number_of_binaries, rng=RNG)
    elif eccentricity_scheme == 'uniform':
        ecc_dist = distributions.UniformEccentricity()
        eccentricities = ecc_dist.rvs(size=number_of_binaries, rng=RNG)
    elif eccentricity_scheme == 'zero':
        ecc_dist = distributions.ZeroEccentricity()
        eccentricities = ecc_dist.rvs(size=number_of_binaries, rng=RNG)
    else: # pragma: no cover
        # This should never be reached
        pass

    return eccentricities


def generate_primary_masses(number_of_binaries=1,
                            primary_mass_min=7,
                            primary_mass_max=120,
                            primary_mass_scheme='Salpeter',
                            **kwargs):
    """Generate random primary masses.

    Use the scheme defined in this particular instance of BinaryPopulation.

    Parameters
    ----------
    number_of_binaries : int
        Number of binaries that require randomly sampled orbital separations
    primary_mass_min : float
        Minimum primary mass
    primary_mass_max : float
        Maximum primary mass
    primary_mass_scheme : string
        Distribution from which the primary masses are randomly drawn

    Returns
    -------
    primary_masses : ndarray of floats
        Randomly drawn primary masses

    """
    RNG = kwargs.get('RNG', np.random.default_rng())

    primary_mass_scheme_options = ['Salpeter', 'Kroupa1993', 'Kroupa2001']

    if primary_mass_scheme not in primary_mass_scheme_options:
        raise ValueError("You must provide an allowed primary mass scheme.")

    # Salpeter E. E., 1955, ApJ, 121, 161
    if primary_mass_scheme == 'Salpeter':
        imf = IMFs.Salpeter(alpha=2.35, m_min=primary_mass_min, m_max=primary_mass_max)
        primary_masses = imf.rvs(size=number_of_binaries, rng=RNG)

    # Kroupa P., Tout C. A., Gilmore G., 1993, MNRAS, 262, 545
    elif primary_mass_scheme == 'Kroupa1993':
        imf = IMFs.Kroupa1993(alpha=2.7, m_min=primary_mass_min, m_max=primary_mass_max)
        primary_masses = imf.rvs(size=number_of_binaries, rng=RNG)

    # Kroupa P., 2001, MNRAS, 322, 231
    elif primary_mass_scheme == 'Kroupa2001':
        imf = IMFs.Kroupa2001(m_min=primary_mass_min, m_max=primary_mass_max)
        primary_masses = imf.rvs(size=number_of_binaries, rng=RNG)
    else: # pragma: no cover
        pass

    return primary_masses


def generate_secondary_masses(primary_masses,
                              number_of_binaries=1,
                              secondary_mass_min=0.35,
                              secondary_mass_max=120,
                              secondary_mass_scheme='flat_mass_ratio',
                              **kwargs):
    """Generate random secondary masses.

    Use the scheme defined in this particular instance of BinaryPopulation.

    Parameters
    ----------
    primary_masses : ndarray of floats
        Previously drawn primary masses
    number_of_binaries : int
        Number of binaries that require randomly sampled orbital separations
    secondary_mass_min : float
        Minimum secondary mass
    secondary_mass_max : float
        Maximum secondary mass
    secondary_mass_scheme : string
        Distribution from which the secondary masses are randomly drawn

    Returns
    -------
    secondary_masses : ndarray of floats
        Randomly drawn secondary masses

    """
    RNG = kwargs.get('RNG', np.random.default_rng())

    secondary_mass_scheme_options = ['flat_mass_ratio', 'q=1']

    # Input parameter checks
    if secondary_mass_scheme not in secondary_mass_scheme_options:
        raise ValueError("You must provide an allowed secondary mass scheme.")

    if np.min(primary_masses) < secondary_mass_min:
        raise ValueError("`secondary_mass_min` is "
                         "larger than some primary masses")

    # Generate secondary masses
    if secondary_mass_scheme == 'flat_mass_ratio':
        # Calculate mass ratio bounds for each primary mass
        q_min = np.maximum(secondary_mass_min / primary_masses, 0.05)
        q_max = np.minimum(secondary_mass_max / primary_masses, 1.0)

        # Sample mass ratios using the distribution class
        # For mass-dependent bounds, we need to sample individually
        mass_ratios = np.zeros(number_of_binaries)
        for i in range(number_of_binaries):
            q_dist = distributions.FlatMassRatio(q_min=q_min[i], q_max=q_max[i])
            mass_ratios[i] = q_dist.rvs(size=1, rng=RNG)[0]

        secondary_masses = primary_masses * mass_ratios

    if secondary_mass_scheme == 'q=1':
        secondary_masses = primary_masses

    return secondary_masses

def generate_binary_fraction(m1=None, binary_fraction_const=1,
                             binary_fraction_scheme='const', **kwargs):
    """
    Getting the binary fraction depending on the scheme. The two possible
    option are a constant binary fraction and a binary fraction based on the
    values given in Moe and Di Stefano (2017).

    Parameters:
    --------------------
    binary scheme: string
        Determines if the value of the binary fraction will be constant or not
    binary fraction const: int
        Gives the value the constant value of the binary if the constant scheme
        is choosen.

    Returns
    ------------------
    binary fraction: int

    """
    binary_fraction_scheme_options = ['const','Moe+17-massdependent']

    if m1 is None:
        raise ValueError("There was not a primary mass provided in the inputs. Unable to return a binary fraction")
    elif not isinstance(m1,np.ndarray):
        m1 = np.asarray(m1)
    binary_fraction = np.zeros_like(m1, dtype=float)
    # Input parameter checks
    if binary_fraction_scheme not in binary_fraction_scheme_options:
        raise ValueError("You must provide an allowed binary fraction scheme.")

    if binary_fraction_scheme == 'const':
        binary_fraction = binary_fraction_const

    elif binary_fraction_scheme == 'Moe+17-massdependent':
        binary_fraction[(m1 > 16)] = 0.94
        binary_fraction[(m1 <= 16) & (m1 > 9)] = 0.84
        binary_fraction[(m1 <= 9) & (m1 > 5)] = 0.76
        binary_fraction[(m1 <= 5) & (m1 > 2)] = 0.59
        binary_fraction[(m1 <= 2)] = 0.4

    else:  # pragma: no cover
        pass

    return binary_fraction
