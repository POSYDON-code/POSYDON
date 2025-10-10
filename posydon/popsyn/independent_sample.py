"""Generate the initial parameters for a binary population."""


__authors__ = [
    "Devina Misra <devina.misra@unige.ch>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Scott Coughlin <scottcoughlin2014@u.northwestern.edu>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]


import numpy as np
from scipy.stats import truncnorm
from posydon.utils.common_functions import rejection_sampler
from posydon.popsyn.Moes_distributions import Moe_17_PsandQs

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

    # Generate primary masses
    m1_set = generate_primary_masses(**kwargs)

    if use_Moe_17_PsandQs(orbital_scheme=orbital_scheme, **kwargs):
        # initialize generator for Moe+17-PsandQs
        if _gen_Moe_17_PsandQs is None:
            _gen_Moe_17_PsandQs = Moe_17_PsandQs()
        # use same defaults as generate_primary_masses
        M1_min = kwargs.get("primary_mass_min", 7)
        M1_max = kwargs.get("primary_mass_max", 120)
        # generate samples
        m2_set, orbital_scheme_set, eccentricity_set, metallicity_set\
         = _gen_Moe_17_PsandQs(m1_set, M_min=M1_min, M_max=M1_max,
                               all_binaries=False)
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
        # compute periods as if all M1 <= 15Msun (where pi = 0.0)
        orbital_periods_M_lt_15 = 10**RNG.uniform(
            low=np.log10(orbital_period_min),
            high=np.log10(orbital_period_max),
            size=number_of_binaries)

        # compute periods as if all M1 > 15Msun
        def pdf(logp):
            pi = 0.55
            beta = 1 - pi
            A = np.log10(10**0.15)**(-pi)*(np.log10(10**0.15)
                                           - np.log10(orbital_period_min))
            B = 1./beta*(np.log10(orbital_period_max)**beta
                         - np.log10(10**0.15)**beta)
            C = 1./(A + B)
            pdf = np.zeros(len(logp))

            for j, logp_j in enumerate(logp):
                # for logP<=0.15 days, the pdf is uniform
                if np.log10(orbital_period_min) <= logp_j and logp_j < 0.15:
                    pdf[j] = C*0.15**(-pi)

                # original Sana H., et al., 2012, Science, 337, 444
                elif 0.15 <= logp_j and logp_j < np.log10(orbital_period_max):
                    pdf[j] = C*logp_j**(-pi)

                else:
                    pdf[j] = 0.

            return pdf

        orbital_periods_M_gt_15 = 10**(rejection_sampler(
            size=number_of_binaries,
            x_lim=[np.log10(orbital_period_min), np.log10(orbital_period_max)],
            pdf=pdf))

        orbital_periods = np.where(primary_masses <= 15.0,
                                   orbital_periods_M_lt_15,
                                   orbital_periods_M_gt_15)
        #TODO orbital for primary mass< current limit set orbital period = 1e8
        orbital_periods = np.where(primary_masses <= 4.00,
                                   1e8,
                                   orbital_periods)
        
    else:
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
        orbital_separations = 10**RNG.uniform(
            low=np.log10(orbital_separation_min),
            high=np.log10(orbital_separation_max),
            size=number_of_binaries)

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

        # Set limits for truncated normal distribution
        a_low = (np.log10(orbital_separation_min)
                 - log_orbital_separation_mean) / log_orbital_separation_sigma
        a_high = (np.log10(orbital_separation_max)
                  - log_orbital_separation_mean) / log_orbital_separation_sigma

        # generate orbital separations from a truncted normal distribution
        log_orbital_separations = truncnorm.rvs(
            a_low, a_high,
            loc=log_orbital_separation_mean,
            scale=log_orbital_separation_sigma,
            size=number_of_binaries,
            random_state=RNG)
        orbital_separations = 10**log_orbital_separations

    else:
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
        eccentricities = np.sqrt(RNG.uniform(size=number_of_binaries))
    elif eccentricity_scheme == 'uniform':
        eccentricities = RNG.uniform(size=number_of_binaries)
    elif eccentricity_scheme == 'zero':
        eccentricities = np.zeros(number_of_binaries)
    else:
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

    primary_mass_scheme_options = ['Salpeter', 'Kroupa1993', 'Kroupa2001','GPL_IMF']

    if primary_mass_scheme not in primary_mass_scheme_options:
        raise ValueError("You must provide an allowed primary mass scheme.")

    # Salpeter E. E., 1955, ApJ, 121, 161
    if primary_mass_scheme == 'Salpeter':
        alpha = 2.35
        normalization_constant = (1.0-alpha) / (primary_mass_max**(1-alpha)
                                                - primary_mass_min**(1-alpha))
        random_variable = RNG.uniform(size=number_of_binaries)
        primary_masses = (random_variable*(1.0-alpha)/normalization_constant
                          + primary_mass_min**(1.0-alpha))**(1.0/(1.0-alpha))

    # Kroupa P., Tout C. A., Gilmore G., 1993, MNRAS, 262, 545
    elif primary_mass_scheme == 'Kroupa1993':
        alpha = 2.7
        normalization_constant = (1.0-alpha) / (primary_mass_max**(1-alpha)
                                                - primary_mass_min**(1-alpha))
        random_variable = RNG.uniform(size=number_of_binaries)
        primary_masses = (random_variable*(1.0-alpha)/normalization_constant
                          + primary_mass_min**(1.0-alpha))**(1.0/(1.0-alpha))

    # Kroupa P., 2001, MNRAS, 322, 231
    elif primary_mass_scheme == 'Kroupa2001':
        alpha = 2.3
        normalization_constant = (1.0-alpha) / (primary_mass_max**(1-alpha)
                                                - primary_mass_min**(1-alpha))
        random_variable = RNG.uniform(size=number_of_binaries)
        primary_masses = (random_variable*(1.0-alpha)/normalization_constant
                          + primary_mass_min**(1.0-alpha))**(1.0/(1.0-alpha))
    
    #General power-law IMF (GPL_IMF)
    elif primary_mass_scheme == 'GPL_IMF':
        alpha1 = kwargs.get('alpha1',2.35)
        alpha2 = kwargs.get('alpha2',2.35)
        m1 = kwargs.get('m1',primary_mass_min)

        normalization_constant = (1.0/(alpha1 + 1.0) *(m1**(alpha1 + 1.0 ) - primary_mass_min**(alpha1 + 1.0 ) ) 
                            + ((m1**(alpha1-alpha2))/(alpha2 +1.0 )) * (primary_mass_max**(alpha2+1.0 ) - m1**(alpha2 + 1.0 )))**(-1.0)

        random_variable = RNG.uniform(size=number_of_binaries)
        #The lower part of the imf 
        def f1(x):
            return ((alpha1 + 1.0)/normalization_constant * x + primary_mass_min**(alpha1 + 1.0 ) )**(1.0/(1.0+alpha1))
        #The upper part of the imf
        def f2(x):
            return (((alpha2 + 1)/( m1**(alpha1-alpha2)))*((x/normalization_constant) - (m1**(alpha1 +1 ) 
                                - primary_mass_min**(alpha1 +1 ))/(alpha1+ 1)) + m1**(alpha2+1))**(1/(alpha2 +1.0 ))

        x1 = (normalization_constant/(alpha1 +1)* ( m1**(alpha1 +1 ) - primary_mass_min**(alpha1 +1 )))
        
        primary_masses = np.where(random_variable < x1, f1(random_variable), f2(random_variable))

    else:
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
        mass_ratio_min = np.max([secondary_mass_min / primary_masses,
                                 np.ones(len(primary_masses))*0.05], axis=0)
        mass_ratio_max = np.min([secondary_mass_max / primary_masses,
                                 np.ones(len(primary_masses))], axis=0)
        secondary_masses = (
            (mass_ratio_max - mass_ratio_min) * RNG.uniform(
                size=number_of_binaries) + mass_ratio_min) * primary_masses

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
    binary_fraction_scheme_options = ['const','Moe+17-massdependent','Gotberg']

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
    elif binary_fraction_scheme == 'Gotberg':
        binary_fraction = 0.09 + 0.63*np.log10(m1)
    else: 
        pass
    return binary_fraction

    
