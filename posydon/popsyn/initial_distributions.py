"""Build the initial-condition distributions from the population parameters.

The population sampler (`independent_sample`) and the reweighting
(`norm_pop`) both use these functions, so a population is reweighted with
the distributions it was sampled from.
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
    '''Get the IMF from `primary_mass_scheme`.

    Any concrete IMF class in `posydon.popsyn.IMFs` is supported. Its
    parameters can be given under the scheme name, e.g.
    `kwargs['Salpeter'] = {'alpha': 2.35}`.

    Parameters
    ----------
    kwargs : dict
        Population parameters, including `primary_mass_scheme`,
        `primary_mass_min` and `primary_mass_max`.

    Returns
    -------
    IMFs.IMFBase
        The IMF.

    Raises
    ------
    ValueError
        If `primary_mass_scheme` is not an IMF class.
    '''

    primary_mass_scheme = kwargs.get('primary_mass_scheme', '')
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
    """Get the mass ratio distribution from `secondary_mass_scheme`.

    Supported schemes:
    - `flat_mass_ratio`: with `secondary_mass_min` and `secondary_mass_max`
      (bounds depend on m1, as in the sampler), or with `q_min` and `q_max`.
    - `power_law_mass_ratio`: with `mass_ratio_slope` and optional `q_min`
      (default 0.05) and `q_max` (default 1.0).

    Parameters
    ----------
    kwargs : dict
        Population parameters.

    Returns
    -------
    MassRatioDistribution
        The mass ratio distribution.

    Raises
    ------
    KeyError
        If only one of `q_min` and `q_max` is given for `flat_mass_ratio`.
    ValueError
        If the scheme is not supported, or if `primary_mass_min` is below
        `secondary_mass_min`.
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
    '''Get the binary fraction from `binary_fraction_scheme`.

    Only `const`, with `binary_fraction_const`, is supported.

    Parameters
    ----------
    kwargs : dict
        Population parameters.

    Returns
    -------
    ConstantBinaryFraction
        The binary fraction.

    Raises
    ------
    ValueError
        If the scheme is not supported.
    '''
    scheme = kwargs.get('binary_fraction_scheme')
    if scheme == 'const':
        return ConstantBinaryFraction(f_b=kwargs['binary_fraction_const'])
    else:
        raise ValueError(f"Binary fraction scheme not recognized: '{scheme}'. "
                         "Supported schemes: ['const']")

def get_orbital_distribution(kwargs):
    '''Get the orbital period distribution from `orbital_scheme`.

    All returned distributions provide `pdf_per_dex(P, m1, q)`.

    Supported schemes:
    - `period`, with `orbital_period_min`, `orbital_period_max` and
      `orbital_period_scheme`: `Sana+12_period_extended`, `power_law` or
      `log_power_law` (both with `power_law_slope`).
    - `separation`, with `orbital_separation_min`, `orbital_separation_max`
      and `orbital_separation_scheme`: `log_uniform` or `log_normal` (with
      `log_orbital_separation_mean` and `log_orbital_separation_sigma`),
      converted to periods with `PeriodFromSeparation`.

    Parameters
    ----------
    kwargs : dict
        Population parameters.

    Returns
    -------
    object
        The orbital period distribution.

    Raises
    ------
    ValueError
        If a scheme is not recognized.
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
                mean=kwargs.get('log_orbital_separation_mean'),
                sigma=kwargs.get('log_orbital_separation_sigma'),
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
