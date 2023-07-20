"""Various utility functions used for the manipulating grid data."""

import os
import gzip
import warnings

import numpy as np
import pandas as pd

from posydon.utils.constants import clight, Msun, Rsun
from posydon.utils.constants import standard_cgrav as cgrav
from posydon.utils.constants import secyer as secyear


__authors__ = [
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Scott Coughlin <scottcoughlin2014@u.northwestern.edu>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]


def join_lists(A, B):
    """Make a joint list of A and B without repetitions and keeping the order.

    Parameters
    ----------
    A : iterable
        The baseline list to be extended.
    B : iterable
        The second list from which additional elements will be appened.

    Returns
    -------
    list
        The joint list by merging A and B.

    """
    result = []
    for element in A:
        result.append(element)
    A_set = set(A)
    for element in B:
        if element not in A_set:
            result.append(element)
    return result


def read_MESA_data_file(path, columns):
    """Read specific columns from a MESA output file to an array.

    Note that this function also works with `.gz` files.

    Parameters
    ----------
    path : str
        The path to the file.
    columns : list of strings
        The list of names of the columns to be returned.

    Returns
    -------
    array or None
        The array containing the MESA data, or `None` if path was `None`.

    """
    if path is None:
        return None
    elif os.path.exists(path):
        try:
            return np.atleast_1d(np.genfromtxt(path, skip_header=5, names=True,
                                               usecols=columns,
                                               invalid_raise=False))
        except:
            warnings.warn("Problems with reading file "+path)
            return None
    else:
        return None


def read_EEP_data_file(path, columns):
    """Read an EEP file (can be `.gz`) - similar to `read_MESA_data_file()`."""
    try:
        return np.atleast_1d(np.genfromtxt(path, skip_header=11, names=True,
                                           usecols=columns, invalid_raise=False))
    except:
        warnings.warn("Problems with reading file "+path)
        return None


def fix_He_core(history):
    """Make He core mass/radius at least equal to CO core mass/radius."""
    if history is not None:
        columns = history.dtype.names
        if "he_core_mass" in columns and "co_core_mass" in columns:
            history["he_core_mass"] = np.maximum(history["he_core_mass"],
                                                 history["co_core_mass"])
        if "he_core_radius" in columns and "co_core_radius" in columns:
            history["he_core_radius"] = np.maximum(history["he_core_radius"],
                                                   history["co_core_radius"])
    return history


def add_field(a, descr):
    """Return a new array that is like `a`, but has additional fields.

    The contents of `a` are copied over to the appropriate fields in the new
    array, whereas the new fields are uninitialized. The arguments are not
    modified.

    Parameters
    ----------
    a : structured numpy array
        The base array.
    descr : str
        A numpy type description of the new fields.

    Returns
    -------
    type
        Description of returned object.

    Example
    -------
    >>> sa = numpy.array([(1, 'Foo'), (2, 'Bar')],
                         dtype=[('id', int), ('name', 'S3')])
    >>> sa.dtype.descr == numpy.dtype([('id', int), ('name', 'S3')])
    True
    >>> sb = add_field(sa, [('score', float)])
    >>> sb.dtype.descr == numpy.dtype([('id', int), ('name', 'S3'),
                                       ('score', float)])
    True
    >>> numpy.all(sa['id'] == sb['id'])
    True
    >>> numpy.all(sa['name'] == sb['name'])
    True

    """
    if a.dtype.fields is None:
        raise ValueError("`A' must be a structured numpy array")
    b = np.empty(a.shape, dtype=a.dtype.descr + descr)
    for name in a.dtype.names:
        b[name] = a[name]
    return b


def get_cell_edges(grid_x, grid_y):
    """Return the edges of grid cells to be used for plotting with pcolormesh.

    Parameters
    ----------
    grid_x : array
        Center values of grid cells for X.
    grid_y : array
        Center values of grid cells for Y.

    Returns
    -------
    array
        X grid values.
    array
        Y grid values.

    """
    q, p = grid_x, grid_y

    shift = (q[1]-q[0])/2
    new_q = q-shift

    last = q[-1]+shift
    new_q = list(new_q)
    new_q.append(last)

    shift = (p[1]-p[0])/2
    new_p = p-shift
    last = p[-1]+shift
    new_p = list(new_p)
    new_p.append(last)

    # make mesh
    QQ, PP = np.meshgrid(new_q, new_p)
    return QQ, PP


def find_nearest(val, array):
    """Find the element of `array` closest to the value `val`."""
    nearest_idx = (abs(val-array)).argmin()
    return array[nearest_idx]


def find_index_nearest_neighbour(array, value):
    """Find the index of `array` closest to the value."""
    idex = np.argmin(np.abs(array - value))
    return idex


def get_final_proposed_points(proposed_x, grid_x, proposed_y, grid_y):
    """Map proposed simulation points to the center of a correspoding cell.

    Only one mapped value is returned if more than one proposed point falls
    within one cell.

    Parameters
    ----------
    proposed_x : array
        Proposed points for X.
    grid_x : array
        The values of X at the center of the grid cells.
    proposed_y : array
        Proposed points for Y.
    grid_y : type
        The values of Y at the center of the grid cells.

    Returns
    -------
    array
        Mapped X values.
    array
        Mapped Y values.

    """
    mapped_x = []
    mapped_y = []

    for i, x in enumerate(proposed_x):
        mapped_x.append(find_nearest(x, grid_x))

    for i, y in enumerate(proposed_y):
        mapped_y.append(find_nearest(y, grid_y))

    # filtering: I only want one unique mapped value in each cell
    coords = np.array([mapped_x, mapped_y])
    unique_coord = np.unique(coords.T, axis=0)
    mapped_x, mapped_y = unique_coord.T[0], unique_coord.T[1]
    print('Number of unique mapped coordinates:', len(mapped_x))

    return mapped_x, mapped_y


def T_merger_P(P, m1, m2):
    """Merger time given initial orbital period and masses of binary.

    Parameters
    ----------
    P : float
        Orbital period (days).
    m1 : float
        Mass of first star (Msun).
    m2 : float
        Mass of second star (Msun).

    Returns
    -------
    float
        Merger time (Gyr)

    """
    return T_merger_a(kepler3_a(P, m1, m2), m1, m2)


def beta_gw(m1, m2):
    """Evaluate the "beta" (equation 5.9) from Peters (1964).

    Parameters
    ----------
    m1 : float
        Mass of the first star.
    m2 : type
        Mass of the second star.

    Returns
    -------
    float
        Peters' beta in cgs units.

    """
    return 64. / 5. * cgrav**3 / clight**5 * m1 * m2 * (m1 + m2)


def kepler3_a(P, m1, m2):
    """Calculate the semimajor axis of a binary from its period and masses.

    Parameters
    ----------
    P : float
        Orbital period (days).
    m1 : float
        Mass of first star.
    m2 : type
        Mass of second star.

    Returns
    -------
    float
        Semi-major axis (Rsun) using Kepler's third law.

    """
    return ((P * 24.0 * 3600.0)**2.0
            * cgrav * (m1 + m2) * Msun / (4.0 * np.pi**2))**(1.0 / 3.0) / Rsun


def T_merger_a(a, m1, m2):
    """Merger time given initial orbital separation and masses of binary.

    Parameters
    ----------
    a : float
        Orbital separation (Rsun).
    m1 : float
        Mass of first star (Msun).
    m2 : float
        Mass of second star (Msun).

    Returns
    -------
    float
        Merger time (Gyr) following Peters (1964), eq. (5.10).

    """
    return (a * Rsun)**4 / (4. * beta_gw(m1, m2) * Msun**3) / (secyear * 1.0e9)


def convert_output_to_table(
    output_file, binary_history_file=None,
    star1_history_file=None, star2_history_file=None, column_names=[
        "result", "CE_flag", "M_1f(Msun)", "M_2f(Msun)",
        "Porb_f(d)", "tmerge(Gyr)",
        "log_L_1", "log_T_1", "He_core_1(Msun)", "C_core_1(Msun)",
        "log_L_2", "log_T_2", "He_core_2(Msun)", "C_core_2(Msun)"]):
    """Convert output of a run, to a pandas dataframe.

    Note that this function also works with `.gz` files.

    Parameters
    ----------
    output_file : str
        Path to MESA output file.
    binary_history_file : str
        Path to binary history data.
    star1_history_file : str
        Path to history data for star1.
    star2_history_file : str
        Path to history data for star1.
    column_names : list of str
        Which columns to return.

    Returns
    -------
    pandas dataframe
        The table containing the requested columns from the run's histories.

    """
    values = {}
    # read out.txt file (simulation terminal output)
    if not os.path.isfile(output_file):
        raise ValueError("File {0} does not exist".format(output_file))

    elif os.stat(output_file).st_size == 0:
        raise ValueError('The output file does not have any data.')

    if binary_history_file is not None:
        binary_history = pd.DataFrame(np.genfromtxt(binary_history_file,
                                                    skip_header=5, names=True))
    else:
        warnings.warn(
            "You have not supplied a binary history file to parse. "
            "This will cause all binary history columns to be dashes.")

    if star1_history_file is not None:
        star1_history = pd.DataFrame(np.genfromtxt(star1_history_file,
                                                   skip_header=5, names=True))
    else:
        warnings.warn(
            "You have not supplied a star1 history file to parse. "
            "This will cause all star1 history columns to be dashes.")

    if star2_history_file is not None:
        star2_history = pd.DataFrame(np.genfromtxt(star2_history_file,
                                                   skip_header=5, names=True))
    else:
        warnings.warn("You have not supplied a star2 history file to parse. "
                      "This will cause all star2 history columns to be dashes")

    # TODO: is this necessary to be done with `popen`?
    # outdata = os.popen('cat ' + output_file).read()
    if output_file.endswith(".gz"):
        with gzip.open(output_file, "rt") as f:
            outdata = f.read()
    else:
        with open(output_file, "r") as f:
            outdata = f.read()

    # Checking for individual terminating conditions

    # carbon depletion
    depleted_c = outdata.find('Terminate due to primary depleting carbon')
    # evolved into contact
    contact = outdata.find("Terminate because accretor "
                           "(r-rl)/rl > accretor_overflow_terminate")
    # CE happened
    CE_flag = outdata.find("TURNING ON CE")
    # initial RLOF overflow, too compact orbit
    overflow = outdata.find("model is overflowing at ZAMS")
    # CE merger
    CE_merge = outdata.find("below CE_xa_diff_to_terminate")
    # CE merger
    CE_merge2 = outdata.find("CE_terminate_when_core_overflows")

    if overflow != -1:
        values["result"] = "ZAMS_RLOF"
    else:
        if CE_flag != -1:
            values["CE_flag"] = 1
        else:
            values["CE_flag"] = 0

        if (contact != -1) and (CE_flag == -1):
            values["result"] = "contact"
        elif ((CE_flag != -1) and (CE_merge != -1 or CE_merge2 != -1
                                   or contact != -1)):
            values["result"] = "CE_merger"

        elif depleted_c != -1:
            if star1_history_file is not None:
                hist = star1_history
                values["log_L_1"] = hist["log_L"].iloc[-1]
                values["log_T_1"] = hist["log_Teff"].iloc[-1]
                values["He_core_1(Msun)"] = hist["he_core_mass"].iloc[-1]
                values["C_core_1(Msun)"] = hist["c_core_mass"].iloc[-1]

            if star2_history_file is not None:
                hist = star2_history
                values["log_L_2"] = hist["log_L"].iloc[-1]
                values["log_T_2"] = hist["log_Teff"].iloc[-1]
                values["He_core_2(Msun)"] = hist["he_core_mass"].iloc[-1]
                values["C_core_2(Msun)"] = hist["c_core_mass"].iloc[-1]

            if binary_history_file is not None:
                tmerge = T_merger_P(
                    binary_history["period_days"].iloc[-1],
                    binary_history["star_1_mass"].iloc[-1],
                    binary_history["star_2_mass"].iloc[-1])
                max_lg_mtransfer_rate = binary_history[
                    "lg_mtransfer_rate"].max()

                values["M_1f(Msun)"] = binary_history["star_1_mass"].iloc[-1]
                values["M_2f(Msun)"] = binary_history["star_2_mass"].iloc[-1]
                values["Porb_f(d)"] = binary_history["period_days"].iloc[-1]
                values["tmerge(Gyr)"] = tmerge

            if max_lg_mtransfer_rate < -5:
                values["result"] = "no_interaction"
            elif CE_flag != -1:
                if tmerge > 13.8:
                    values["result"] = "CE_ejection"
                else:
                    values["result"] = "CE_ejection_m"
            else:
                if tmerge > 13.8:
                    values["result"] = "stable_MT_to_wide_binary"
                else:
                    values["result"] = "stable_MT_to_merger"
        else:
            values["result"] = "error"

    # If there is any reason that a column cannot be reported based on
    # the MESA simulation we fill this value with a "-"
    for column in column_names:
        if column not in values.keys():
            values[column] = "-"

    return pd.DataFrame.from_dict(values, orient='index').transpose().apply(
        pd.to_numeric, errors='ignore')


def clean_inlist_file(inlist, **kwargs):
    """Get dictionary of possible parameters/values from the default inlist.

    Parameters
    ----------
    inlist : str
        Path to inlist file.
    **kwargs : dict
        TODO

    Returns
    -------
    dict
        Dictionary of parameters and values from inlist.

    """
    section = kwargs.pop('section', None)

    with open(inlist, 'r') as f:
        # clean inlist into nice list
        param_value_list = []
        for i in f.readlines():
            # strip away all the comments and whitespace, etc
            param_and_value = i.strip('\n').strip().split('!')[0].strip()
            # check that this line actually has a parameter and value par
            if ((param_and_value) and ('=' in param_and_value
                                       or '&' in param_and_value)):
                param_value_list.append(param_and_value)

    # does this inlist have multiple sections?
    dict_of_parameters = {k: {} for k in param_value_list if '&' in k}
    if not dict_of_parameters:
        # mesa default files do not have sections,
        # because controls and jobs are in separate files
        dict_of_parameters[section] = {
            i.split('=', 1)[0].strip(): i.split('=', 1)[1].strip()
            for i in param_value_list
            if (i.split('=', 1)[1].strip() not in ["''", "'.'"])}
    else:
        # probably other people are sending the same inlist with
        # both job and controls sections
        for i in param_value_list:
            if '&' in i:
                current_key = i
            elif (i.split('=', 1)[1].strip() not in ["''", "'.'"]):
                dict_of_parameters[current_key][i.split('=', 1)[0].strip()] = \
                    i.split('=', 1)[1].strip()

    return dict_of_parameters

def get_new_grid_name(path, compression, create_missing_directories=False):
    """Get the name of a new grid slice based on the path and the compression.

    Parameters
    ----------
    path : str
        Path to grid slice data.
    compression : str
        Compression value. (Directory to put the new grid slice in.)
    create_missing_directories : bool
        Flag to create missing directories.

    Returns
    -------
    grid_output
        File name for the new grid slice.

    """
    grid_path, grid_name = os.path.split(os.path.normpath(path))
    output_path = os.path.join(grid_path, compression)
    grid_output = os.path.join(output_path, grid_name+'.h5')
    if create_missing_directories:
        # check that LITE/ or ORIGINAL/ directory exists
        if not os.path.isdir(output_path):
            os.makedirs(output_path)
    return grid_output


