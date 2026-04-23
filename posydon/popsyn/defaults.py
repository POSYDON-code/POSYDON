"""Default binary population parameters."""


__authors__ = [
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
]


from posydon.utils.constants import age_of_universe

default_kwargs = {
    # Random Number Generation
    'entropy': None,   # `None` uses system entropy (recommended)

    # Size of the population
    'number_of_binaries': 100,
    'metallicities' : [1.0], # in Zsolar

    # Star Formation History
    'star_formation': 'constant',
    # burst_time: 0,
    # max_time: 0,
    # min_time: age_of_universe,
    'max_simulation_time': age_of_universe,

    # Orbital Period
    'orbital_scheme': "period",    # alternatively: 'period'

    # if orbital_scheme == separation
    'orbital_separation_scheme': 'log_uniform',
    'orbital_separation_min': 5.0,  # in Rsun
    'orbital_separation_max': 1.0e5,  # in Rsun
    'log_orbital_seperation_mean': None,
    'log_orbital_seperation_sigma': None,

    # if orbital_scheme == period
    'orbital_period_scheme': 'Sana+12_period_extended',
    'orbital_period_min': 0.75,  # in days
    'orbital_period_max': 6000,  # in days
    # 'log_orbital_seperation_mean': None,
    # 'log_orbital_seperation_sigma': None,

    # Eccentricity
    'eccentricity_scheme': 'zero',

    # Primary mass
    'primary_mass_scheme': 'Kroupa2001',
    'primary_mass_min': 7.0,
    'primary_mass_max': 120.0,

    # Secondary mass
    'secondary_mass_scheme': 'flat_mass_ratio',
    'secondary_mass_min': 0.35,
    'secondary_mass_max': 120.0,

    #Binary fraction
    'binary_fraction_const' : 1,
    'binary_fraction_scheme' : 'const',
}
