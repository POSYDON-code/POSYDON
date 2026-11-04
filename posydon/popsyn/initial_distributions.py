"""Build the initial-condition distributions from population parameters.

These factories are the single place where the population parameters
(e.g. `primary_mass_scheme`, `orbital_period_min`) are turned into
distribution objects. They are used both by the population sampler
(`posydon.popsyn.independent_sample`) and by the normalisation and
reweighting of populations (`posydon.popsyn.norm_pop`), so that a
population is reweighted with the same distributions it was sampled from.
"""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import inspect

import posydon.popsyn.IMFs as IMFs
from posydon.popsyn.distributions import (
    ConditionalFlatMassRatio,
    ConstantBinaryFraction,
    FlatMassRatio,
    LogNormalSeparation,
    LogPowerLawPeriod,
    LogUniform,
    PeriodFromSeparation,
    PowerLawMassRatio,
    PowerLawPeriod,
    Sana12Period,
)


def get_IMF_distribution(kwargs):
    '''get the IMF distribution

    Supported schemes are the (non-abstract) IMF classes in the IMFs module,
    e.g. `Salpeter`, `Kroupa1993`, `Kroupa2001`, `Chabrier2003`.
    Additional parameters can be passed to the scheme
    by adding them to the kwargs dictionary with the scheme name as the key
    Example:

    ```
    kwargs = {
        'primary_mass_scheme': 'Salpeter',
        'primary_mass_min': 0.1,
        'primary_mass_max': 100,
        'Salpeter': {'alpha': 2.35}
    }
    ```

    Parameters
    ----------
    kwargs : dict
        Dictionary containing the parameters.
        `primary_mass_scheme`, `primary_mass_min`, `primary_mass_max` are required.
        Additional parameters are required depending on the scheme.

    Returns
    -------
    imf : IMFs.IMFBase
        IMF distribution with `pdf(m1)`.

    Raises
    ------
    ValueError
        If `primary_mass_scheme` is not an IMF class of the IMFs module.
    '''

    primary_mass_scheme = kwargs.get('primary_mass_scheme', '')
    # dynamically retrieve the IMF class from the IMFs module
    imf_class = (getattr(IMFs, primary_mass_scheme, None)
                 if isinstance(primary_mass_scheme, str) else None)
    if not _is_IMF_class(imf_class):
        supported = sorted(name for name, obj in vars(IMFs).items()
                           if _is_IMF_class(obj))
        raise ValueError(f"The primary_mass_scheme '{primary_mass_scheme}' "
                         f"is not recognized. Supported schemes: {supported}")

    scheme_kwargs = kwargs.get(primary_mass_scheme, {})
    return imf_class(m_min=kwargs['primary_mass_min'],
                     m_max=kwargs['primary_mass_max'],
                     **scheme_kwargs)

def _is_IMF_class(obj):
    '''Check whether obj is a concrete IMF class from the IMFs module.'''
    return (inspect.isclass(obj)
            and issubclass(obj, IMFs.IMFBase)
            and not inspect.isabstract(obj))

def get_mass_ratio_distribution(kwargs):
    """Build the mass ratio distribution from the simulation parameters

    The same distribution classes are used by the population sampler, so
    the PDF, its support and the sampled mass ratios agree.

    Supported schemes:
    - `flat_mass_ratio` for `secondary_mass_scheme`
        Requires either of the following parameter sets:
        - `secondary_mass_min`, `secondary_mass_max` and `primary_mass_min`
          (q bounds depend on m1, as in the sampler)
        - `q_min` and `q_max` (fixed q bounds)
    - `power_law_mass_ratio` for `secondary_mass_scheme`
        Requires the following parameters:
        - `mass_ratio_slope`: exponent alpha in q^alpha
        - `q_min` (optional, default 0.05)
        - `q_max` (optional, default 1.0)

    Parameters
    ----------
    kwargs : dict
        Dictionary containing the simulation parameters

    Returns
    -------
    q_dist : MassRatioDistribution
        Distribution with `pdf(q, m1)` and `support(m1)`.

    Raises
    ------
    KeyError
        If only one of `q_min` and `q_max` is given for `flat_mass_ratio`.
    ValueError
        If the `secondary_mass_scheme` is not supported (such as the
        sampler's `q=1`), if `primary_mass_min` is below
        `secondary_mass_min` (as in the sampler), or if the bounds are
        invalid.
    """
    scheme = kwargs.get('secondary_mass_scheme')
    if scheme == 'flat_mass_ratio':
        if 'q_min' in kwargs or 'q_max' in kwargs:
            # fixed bounds; a KeyError is raised if only one is given
            return FlatMassRatio(q_min=kwargs['q_min'], q_max=kwargs['q_max'])
        # bounds depend on m1, as in the sampler
        if kwargs['primary_mass_min'] < kwargs['secondary_mass_min']:
            raise ValueError("`secondary_mass_min` is larger than "
                             "`primary_mass_min`.")
        return ConditionalFlatMassRatio(m2_min=kwargs['secondary_mass_min'],
                                        m2_max=kwargs['secondary_mass_max'])
    elif scheme == 'power_law_mass_ratio':
        return PowerLawMassRatio(alpha=kwargs['mass_ratio_slope'],
                                 q_min=kwargs.get('q_min', 0.05),
                                 q_max=kwargs.get('q_max', 1.0))
    else:
        raise ValueError(f"The secondary_mass_scheme '{scheme}' is not "
                         "supported for normalisation. Supported schemes: "
                         "['flat_mass_ratio', 'power_law_mass_ratio']")

def get_binary_fraction_distribution(kwargs):
    '''get the binary fraction distribution

    Supported schemes:
    - `const` with `binary_fraction_const`
        requires the following parameters to be present:
        - `binary_fraction_const`

    Parameters
    ----------
    kwargs : dict
        Dictionary containing the parameters

    Returns
    -------
    f_b_dist : ConstantBinaryFraction
        Binary fraction distribution with `pdf(binary)` and `f_b`.

    Raises
    ------
    ValueError
        If the `binary_fraction_scheme` is not supported (such as the
        sampler's `Moe+17-massdependent`).
    '''
    scheme = kwargs.get('binary_fraction_scheme')
    if scheme == 'const':
        return ConstantBinaryFraction(f_b=kwargs['binary_fraction_const'])
    else:
        raise ValueError(f"Binary fraction scheme not recognized: '{scheme}'. "
                         "Supported schemes: ['const']")

def get_orbital_distribution(kwargs):
    '''get the orbital period distribution

    All returned distributions have `pdf_per_dex(P, m1, q)`, the probability
    density per dex, i.e. normalised to 1 in log10(P) space, such that they
    can be compared between schemes.

    Supported schemes:
    - `period` for `orbital_scheme` with `orbital_period_scheme`:
        - `Sana+12_period_extended`
        - `power_law`: dN/dlog10(P) ∝ P^`power_law_slope`
        - `log_power_law`: dN/dlog10(P) ∝ (log10(P))^`power_law_slope`
    - `separation` for `orbital_scheme` with `orbital_separation_scheme`,
      converted to a period distribution with `PeriodFromSeparation`:
        - `log_uniform`
        - `log_normal` with `log_orbital_separation_mean` and
          `log_orbital_separation_sigma`

    Parameters
    ----------
    kwargs : dict
        Dictionary containing the simulation parameters

    Returns
    -------
    orbital_dist : object
        Period distribution with `pdf_per_dex(P, m1, q)`.
    '''
    if (kwargs['orbital_scheme'] == 'period'):
        if kwargs['orbital_period_scheme'] == 'Sana+12_period_extended':
            return Sana12Period(
                p_min=kwargs['orbital_period_min'],
                p_max=kwargs['orbital_period_max'],
            )
        elif kwargs['orbital_period_scheme'] == 'power_law':
            return PowerLawPeriod(
                p_min=kwargs['orbital_period_min'],
                p_max=kwargs['orbital_period_max'],
                slope=kwargs['power_law_slope'],
            )
        elif kwargs['orbital_period_scheme'] == 'log_power_law':
            return LogPowerLawPeriod(
                p_min=kwargs['orbital_period_min'],
                p_max=kwargs['orbital_period_max'],
                slope=kwargs['power_law_slope'],
            )
        else:
            raise ValueError("Orbital period scheme not recognized")
    elif (kwargs['orbital_scheme'] == 'separation'):
        if kwargs['orbital_separation_scheme'] == 'log_uniform':
            separation = LogUniform(
                min=kwargs['orbital_separation_min'],
                max=kwargs['orbital_separation_max'],
            )
        elif kwargs['orbital_separation_scheme'] == 'log_normal':
            separation = LogNormalSeparation(
                mean=kwargs['log_orbital_separation_mean'],
                sigma=kwargs['log_orbital_separation_sigma'],
                min=kwargs['orbital_separation_min'],
                max=kwargs['orbital_separation_max'],
            )
        else:
            raise ValueError("Orbital separation scheme not recognized")
        return PeriodFromSeparation(separation)
    else:
        raise ValueError("Orbital scheme not recognized. "
                         "Please specify either 'period' or 'separation' "
                         "and a valid scheme for that orbital scheme.")
