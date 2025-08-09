"""Gravitational wave radiation."""


import posydon.utils.constants as const


def calculate_gravitational_radiation(p):
    """Calculate the effect of gravitational wave radiation on binaries.

    Calculate orbital evolution (change in the orbital separation and
    eccentricity) due to gravitational wave radiation. Equations used are taken
    from eq. 35 and 36 of Junker et al. 1992, 254, 146.

    Parameters
    ----------
    p : tuple
        A tuple containing M_pri (the primary star's mass in Msun), M_sec (the
        secondary star's mass in Msun), a (the orbital separation in Rsun), and
        e (the binary's eccentricity).

    Returns
    -------
    da_gr, de_gr : floats
        The derivatives of the orbital separation in Rsun/s and eccentricity
        in 1/s.

    """
    M_pri, M_sec, a, e = p

    v = (M_pri * M_sec / (M_pri + M_sec) ** 2)

    da_gr = (
        (-2 * const.clight / 15) * (v / ((1 - e ** 2) ** (9 / 2)))
        * (const.standard_cgrav * (M_pri + M_sec) * const.msol
           / (a * const.rsol * const.clight ** 2)) ** 3
        * ((96 + 292 * e ** 2 + 37 * e ** 4) * (1 - e ** 2)
           - (1 / 28 * const.standard_cgrav * (M_pri + M_sec) * const.msol
           / (a * const.rsol * const.clight ** 2))
           * ((14008 + 4707 * v)+(80124 + 21560 * v) * e ** 2
           + (17325 + 10458) * e ** 4 - 0.5 * (5501 - 1036 * v) * e ** 6))
        ) * const.secyer / const.rsol

    de_gr = (
        (-1 / 15) * ((v * const.clight ** 3) / (
            const.standard_cgrav * (M_pri + M_sec) * const.msol))
        * (const.standard_cgrav * (M_pri + M_sec) * const.msol / (
            a * const.rsol * const.clight ** 2)) ** 4
        * (e / (1 - e ** 2) ** (7 / 2))
        * ((304 + 121 * e ** 2) * (1 - e ** 2)
           - (1 / 56 * const.standard_cgrav * (M_pri + M_sec) * const.msol
           / (a * const.rsol * const.clight ** 2))
           * (8 * (16705 + 4676 * v) + 12 * (9082 + 2807 * v) * e ** 2
           - (25211 - 3388 * v) * e ** 4))
        ) * const.secyer

    return da_gr, de_gr
