"""Wind loss."""

__authors__ = [
    "Zepei Xing <Zepei.Xing@unige.ch>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
]


def calculate_wind_loss(M_1, M_2, Mdot_1, a):
    """Calculate the effect on the orbit due to wind mass loss.

    Parameters
    ----------
    M1 : float
        Mass of the current star
    M2 : float
        Mass of the companion star
    Mdot_1 : float
        Mass loss rate of the current star

    Returns
    -------
    da_mt : double
        Change in orbital separation due to mass loss

    """
    q = M_1 / M_2
    k11 = (1 / (1 + q)) * (Mdot_1 / M_1)
    k21 = Mdot_1 / M_1
    k31 = Mdot_1 / (M_1 + M_2)

    # This is simplified to da_mt = -a * Mdot/(M+Macc), for only (negative)
    # wind Mdot from star M.
    da_mt = a * (2 * k11 - 2 * k21 + k31)

    return da_mt
