import posydon.utils.constants as const

def default_gravrad(a, e, primary, secondary, verbose = False):
    """
        Calculates the change in orbital separation and eccentricity 
    due to gravitational wave radiation, according to: 
    
        Junker, W., & Schafer, G. 1992, MNRAS, 254, 146

    Parameters
    ----------
        a : float
            The current orbital separation. [Rsolar]

        e : float
            The current orbital eccentricity.

        primary : SingleStar object
            A single star object, representing the primary (more evolved) star 
            in the binary and containing its properties.
        
        secondary : SingleStar object
            A single star object, representing the secondary (less evolved) star 
            in the binary and containing its properties.

        verbose : bool
            True if we want to print stuff.
    
    Returns
    -------
        da : float
            The change in orbital separation for a time step in 
        step_detached's solve_ivp(). [Rsolar/yr]

        de : float
            The change in orbital eccentricity for a time step in 
        step_detached's solve_ivp().

    """

    m1 = primary.latest["mass"]
    m2 = secondary.latest["mass"]

    v = (m1 * m2 / (m1 + m2) ** 2)
    
    da_gr = (
        (-2 * const.clight / 15) * (v / ((1 - e ** 2) ** (9 / 2)))
        * (const.standard_cgrav * (m1 + m2) * const.msol
        / (a * const.rsol * const.clight ** 2)) ** 3
        * ((96 + 292 * e ** 2 + 37 * e ** 4) * (1 - e ** 2)
        - (1 / 28 * const.standard_cgrav * (m1 + m2) * const.msol
        / (a * const.rsol * const.clight ** 2))
        * ((14008 + 4707 * v)+(80124 + 21560 * v) * e ** 2
        + (17325 + 10458) * e ** 4 - 0.5 * (5501 - 1036 * v) * e ** 6))
        ) * const.secyer / const.rsol
    # eq. (35) of Junker et al. 1992, 254, 146

    de_gr = (
        (-1 / 15) * ((v * const.clight ** 3) / (
            const.standard_cgrav * (m1 + m2) * const.msol))
        * (const.standard_cgrav * (m1 + m2) * const.msol / (
            a * const.rsol * const.clight ** 2)) ** 4
        * (e / (1 - e ** 2) ** (7 / 2))
        * ((304 + 121 * e ** 2) * (1 - e ** 2)
        - (1 / 56 * const.standard_cgrav * (m1 + m2) * const.msol
        / (a * const.rsol * const.clight ** 2))
        * (8 * (16705 + 4676 * v) + 12 * (9082 + 2807 * v) * e ** 2
        - (25211 - 3388 * v) * e ** 4))
        ) * const.secyer
    # eq. (36) of Junker et al. 1992, 254, 146

    if verbose:
        print("da, de_gr = ", da_gr, de_gr)

    da = da_gr
    de = de_gr

    return da, de