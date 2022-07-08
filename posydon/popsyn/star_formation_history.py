"""Implements the selection of different star-formation history scenarios."""


__authors__ = [
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Devina Misra <devina.misra@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
]


import numpy as np
from posydon.utils.constants import age_of_universe
from posydon.utils import (rejection_sampler, histogram_sampler,
                           read_histogram_from_file)


SFH_SCENARIOS = ["burst", "constant", "custom_linear", "custom_log10",
                 "custom_linear_histogram", "custom_log10_histogram"]


def get_formation_times(N_binaries, star_formation='constant', **kwargs):
    """Get formation times of binaries in a population based on a SFH scenario.

    Parameters
    ----------
    N_binaries : int
        Number of formation ages to produce.
    star_formation : str, {constant, burst}
        Constant - random formation times from a uniform distribution.
        Burst - all stars are born at the same time.
    burst_time : float, 0 (years)
        Sets birth time in years.
    min_time : float, 0 (years)
        If constant SF, sets minimum of random sampling.
    max_time : float, age_of_universe (years)
        If constant SF, sets maximum of random sampling.
    RNG : <class, np.random.Generator>
        Random generator instance.

    Returns
    -------
    array
        The formation times array.
    """
    RNG = kwargs.get('RNG', np.random.default_rng())

    scenario = star_formation.lower()

    if scenario == 'burst':
        burst_time = kwargs.get('burst_time', 0.0)
        return np.ones(N_binaries)*burst_time

    max_time_default = kwargs.get('max_simulation_time', age_of_universe)
    max_time = kwargs.get('max_time', max_time_default)

    if scenario == 'constant':
        min_time = kwargs.get('min_time', 0.0)
        return RNG.uniform(size=N_binaries, low=min_time, high=max_time)

    if scenario in ["custom_linear", "custom_log10"]:
        custom_ages_file = kwargs.get('custom_ages_file')
        x, y = np.loadtxt(custom_ages_file, unpack=True)
        current_binary_ages = rejection_sampler(x, y, N_binaries)
        if "log10" in scenario:
            current_binary_ages = 10.0 ** current_binary_ages
        return max_time - current_binary_ages

    if scenario in ["custom_linear_histogram", "custom_log10_histogram"]:
        custom_ages_file = kwargs.get('custom_ages_file')
        x, y = read_histogram_from_file(custom_ages_file)
        current_binary_ages = histogram_sampler(x, y)
        if "log10" in scenario:
            current_binary_ages = 10.0 ** current_binary_ages
        return max_time - current_binary_ages

    raise ValueError(
        "Unknown star formation scenario '{}' given. Valid options: {}".
        format(star_formation, ",".join(SFH_SCENARIOS)))
