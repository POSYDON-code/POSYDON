"""Module for post-processing POSYDON grids."""

from posydon.binary_evol.binarystar import BinaryStar
from posydon.binary_evol.singlestar import SingleStar
from posydon.binary_evol.SN.step_SN import StepSN
from posydon.binary_evol.flow_chart import STAR_STATES_CC
from posydon.utils.common_functions import (
    calculate_Patton20_values_at_He_depl,
    CEE_parameters_from_core_abundance_thresholds,
    check_state_of_star)
from posydon.grids.SN_MODELS import get_SN_MODEL, SN_MODELS
from posydon.visualization.combine_TF import combine_TF12, TF1_POOL_STABLE
from posydon.visualization.plot_defaults import DEFAULT_MARKERS_COLORS_LEGENDS
import numpy as np
from tqdm import tqdm
import copy
from posydon.utils.posydonwarning import (Pwarn, Catch_POSYDON_Warnings)


__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]


__credits__ = [
    'Konstantinos Kovlakas <kkovlakas@physics.uoc.gr>'
]


CC_quantities = ['state', 'SN_type', 'f_fb', 'mass', 'spin',
                 'm_disk_accreted', 'm_disk_radiated',
                 'CO_interpolation_class', 'M4', 'mu4',
                 'h1_mass_ej', 'he4_mass_ej']

def assign_core_collapse_quantities_none(EXTRA_COLUMNS, star_i,
                                         SN_MODEL_NAME=None):
    """Assign None values to all core collapse properties.

        Parameters
        ----------
        EXTRA_COLUMNS : dict of lists
            Each new column should be a list in this dictionary. This
            dictionary will get modified by this function.
        star_i : int (1 or 2)
            Index of the star.
        SN_MODEL_NAME : str or list of str, or None (default: None)
            SN_MODEL key like the ones given in SN_MODELS.py.

    """
    if (not isinstance(star_i, int) or (star_i<1) or (star_i>2)):
        raise ValueError("'star_i' should be 1 or 2.")
    if SN_MODEL_NAME is None:
        for NAME in SN_MODELS.keys():
            for quantity in CC_quantities:
                EXTRA_COLUMNS[f'S{star_i}_{NAME}_{quantity}'].append(None)
    elif isinstance(SN_MODEL_NAME, str):
        for quantity in CC_quantities:
            EXTRA_COLUMNS[f'S{star_i}_{SN_MODEL_NAME}_{quantity}'].append(None)
    elif isinstance(SN_MODEL_NAME, list):
        for NAME in SN_MODEL_NAME:
            for quantity in CC_quantities:
                EXTRA_COLUMNS[f'S{star_i}_{NAME}_{quantity}'].append(None)
    else:
        raise TypeError("'SN_MODEL_NAME' should be a string, a list of "
                        "strings, or None.")

def print_CC_quantities(star, SN_MODEL_NAME=None):
    """Print quantities at core collapse.

        Parameters
        ----------
        star : SingleStar object
            The star whose quantities to print.
        SN_MODEL_NAME : str or None (default: None)
            Name of the supernova model. If None, print pre SN values.

    """
    format_string = "{:<50} {:<33} {:12} {:10} {:15} {:10} {:25} {:25} {:25} "\
                    + "{:25} {:25} {:25}"
    format_val_preSN = "{:<50} {:<33} {:12} {:10} {:7.2f} {:12.2f} {:25} "\
                       + "{:25} {:25} {:25} {:25} {:25}"
    format_val = "{:<50} {:<33} {:12} {:1.2f} {:13.2f} {:12.2f} {:20.2f} "\
                 + "{:20.2f} {:20.2f} {:20.2f} {:20.2f} {:20.2f}"
    if SN_MODEL_NAME is None:
        print('')
        print(format_string.format("mechanism", "state", "SN type", "f_fb",
                                   "mass [Msun]", "spin",
                                   "m_disk_accreted [Msun]",
                                   "m_disk_radiated [Msun]", "M4 [m/Msun]",
                                   "mu4 [(dm/Msun)/(dr/1000km)]",
                                   "h1_mass_ej [Msun]", "he4_mass_ej [Msun]"))
        print('')
        try:
            print(format_val_preSN.format('PRE SN STAR', star.state, '', '',
                                          star.mass, star.spin, '', '', '', '',
                                          '', ''))
        except Exception as e:
            Pwarn("Failed to print star values!\nWarning in preSN: "
                  "{}".format(e), "InappropriateValueWarning")
        print('')
    else:
        try:
            checked_quantities_for_None = {}
            for quantity in ["spin", "M4", "mu4", "h1_mass_ej", "he4_mass_ej"]:
                if getattr(star, quantity) is None:
                    checked_quantities_for_None[quantity] = np.nan
                else:
                    checked_quantities_for_None[quantity] = getattr(star,
                                                                    quantity)
            print(format_val.format(SN_MODEL_NAME, star.state, star.SN_type,
                                    star.f_fb, star.mass,
                                    checked_quantities_for_None["spin"],
                                    star.m_disk_accreted, star.m_disk_radiated,
                                    checked_quantities_for_None["M4"],
                                    checked_quantities_for_None["mu4"],
                                    checked_quantities_for_None["h1_mass_ej"],
                                    checked_quantities_for_None["he4_mass_ej"])
                  )
        except Exception as e:
            Pwarn("Failed to print star values!\nWarning in "
                  "{}: {}".format(SN_MODEL_NAME, e),
                  "InappropriateValueWarning")

def post_process_grid(grid, index=None, star_2_CO=True, SN_MODELS=SN_MODELS,
                      single_star=False, verbose=False):
    """Compute post processed quantity of any grid.

    This function post processes any supported grid and computes:
    - Core collapse quantities for 5 prescriptions given the fiducial POSYDON
    assumption given in SN_MODEL:
      for each prescrition we store the compact object state (WD/NS/BH),
      SN type (WD, ECSN, CCSN, PISN, PPISN), fallback mass fraction f_gb,
      compact object mass and spin.
    - Core masses at he-depletion (used in Patton core collapse)
    - Mass envelopes for common envelope step.

    Parameters
    ----------
    grid : PSyGrid
        MESA grid in PSyGrid format.
    index : None, int, or list with two int (default: None)
        If None, loop over all indicies otherwise provide a range, e.g. [10,20]
        or a single index, e.g. 42.
    star_2_CO : bool (default: True)
        If 'False' star 2 is not a compact object.
    SN_MODELS : list of dict (default are the models defined in SN_MODELS.py)
        List of supported core collapse model assumptions.
    single_star : bool (default: False)
        If `True` the PSyGrid contains single stars.
    verbose : bool (default: False)
        If `True` print the results of each core collapse on screen.

    Returns
    -------
    MESA_dirs: list
        List containing the path to each run corresponding to the post
        processed values. This is used to ensure one to one mapping when
        appending the extra columns back to a grid.
    EXTRA_COLUMNS: dict
        Dictionary containing all post processed quantities.

    """
    EXTRA_COLUMNS = {}
    if single_star: # only star 1 columns in case of single star grid
        stars = [1]
    else:
        stars = [1, 2]
        EXTRA_COLUMNS['mt_history'] = []
    for star in stars:
        # core masses at He depletion. stellar states and composition
        for quantity in ['avg_c_in_c_core_at_He_depletion',
                         'co_core_mass_at_He_depletion',
                         'state', 'surface_other', 'center_other']:
            EXTRA_COLUMNS[f'S{star}_{quantity}'] = []
        # common envelope quantities
        for quantity in ['lambda_CE', 'm_core_CE', 'r_core_CE']:
            for val in [1, 10, 30, 'pure_He_star_10']:
                EXTRA_COLUMNS[f'S{star}_{quantity}_{val}cent'] = []
        # Core collapse qunatities: [state, SN_type, f_fb, mass, spin]
        for SN_MODEL_NAME in SN_MODELS.keys():
            for quantity in CC_quantities:
                EXTRA_COLUMNS[f'S{star}_{SN_MODEL_NAME}_{quantity}'] = []

    # loop over all gird, index or range
    if index is None:
        indicies = range(len(grid.MESA_dirs))
        MESA_dirs = grid.MESA_dirs
    elif isinstance(index, int):
        indicies = range(index, index+1)
        MESA_dirs = grid.MESA_dirs[index:index+1]
    elif isinstance(index, list):
        if len(index) != 2:
            raise ValueError('Index range should have dim=2!')
        indicies = range(index[0], index[1])
        MESA_dirs = grid.MESA_dirs[index[0]:index[1]]
    else:
        raise TypeError('The argument `index` should be None, an integer, or '
                        'a list of two integers.')

    for i in tqdm(indicies):

        if not single_star:
            # initilise binary
            binary = BinaryStar.from_run(grid[i], history=True, profiles=True)
            stars = [binary.star_1, binary.star_2]
            stars_CO = [False, star_2_CO]
            interpolation_class = grid.final_values['interpolation_class'][i]
            IC = grid.final_values['interpolation_class'][i]
            TF2 = grid.final_values['termination_flag_2'][i]
        else:
            star = SingleStar.from_run(grid[i], history=True, profile=True)
            stars = [star]
            stars_CO = [False]
            IC = 'no_MT'
            TF2 = 'no_RLOF'
        TF1 = grid.final_values['termination_flag_1'][i]

        # compute properties
        for j, star in enumerate(stars):
            if not stars_CO[j] and IC in ['no_MT', 'stable_MT', 'unstable_MT',
                                          'stable_reverse_MT']:
                # stellar states
                EXTRA_COLUMNS['S%s_state' % (j+1)].append(check_state_of_star(
                    star, star_CO=False))
                # core masses at he depletion
                with Catch_POSYDON_Warnings(record=True) as cpw:
                    if ((star.avg_c_in_c_core_at_He_depletion is None) or
                        (star.co_core_mass_at_He_depletion is None)):
                        calculate_Patton20_values_at_He_depl(star)
                EXTRA_COLUMNS[f'S{j+1}_avg_c_in_c_core_at_He_depletion'].append(
                    star.avg_c_in_c_core_at_He_depletion)
                EXTRA_COLUMNS[f'S{j+1}_co_core_mass_at_He_depletion'].append(
                    star.co_core_mass_at_He_depletion)
                # CE quantities
                with Catch_POSYDON_Warnings(record=True) as cpw:
                    try:
                        CEE_parameters_from_core_abundance_thresholds(star)
                    except Exception as ex:
                        print(ex)
                        print(f'The exception was raised by {grid.MESA_dirs[i]} '
                              f'in CEE_parameters_from_core_abundance_thresholds(star_{j+1}).')
                for quantity in ['lambda_CE', 'm_core_CE', 'r_core_CE']:
                    for val in [1, 10, 30, 'pure_He_star_10']:
                        EXTRA_COLUMNS[f'S{j+1}_{quantity}_{val}cent'].append(
                                            getattr(star, f'{quantity}_{val}cent'))
                # aboundances
                try:
                    s_o = (1. - star.surface_h1 - star.surface_he4 - star.surface_c12
                           - star.surface_n14 - star.surface_o16)
                    c_o = (1. - star.center_h1 - star.center_he4 - star.center_c12
                           - star.center_n14 - star.center_o16)
                except TypeError as ex:
                    s_o = 0.
                    c_o = 0.
                    print(ex)
                    print(f'The error was raised by {grid.MESA_dirs[i]} '
                           f'while accessing aboundances in star_{j+1}.')
                EXTRA_COLUMNS['S%s_surface_other' % (j+1)].append(s_o)
                EXTRA_COLUMNS['S%s_center_other' % (j+1)].append(c_o)
            else:
                # fill everything with Nones
                if IC == 'initial_MT' or IC == 'not_converged':
                    EXTRA_COLUMNS['S%s_state' % (j+1)].append(None)
                else:
                    # CO states are classified and used in mesa step
                    try:
                        EXTRA_COLUMNS['S%s_state' % (j+1)].append(check_state_of_star(star, star_CO=True))
                    except TypeError as ex:
                        EXTRA_COLUMNS['S%s_state' % (j+1)].append(None)
                        print(ex)
                        print(f'The error was raised by {grid.MESA_dirs[i]} '
                               f'in check_state_of_star(star_{j+1}) with IC={IC}.')
                for quantity in ['avg_c_in_c_core_at_He_depletion', 'co_core_mass_at_He_depletion',
                                 'surface_other', 'center_other', 'lambda_CE', 'm_core_CE', 'r_core_CE']:
                    if 'CE' in quantity:
                        for val in [1, 10, 30, 'pure_He_star_10']:
                            EXTRA_COLUMNS[f'S{j+1}_{quantity}_{val}cent'].append(None)
                    else:
                        EXTRA_COLUMNS[f'S{j+1}_{quantity}'].append(None)

        # core collapse quantities
        if not single_star:
            if interpolation_class in ['no_MT', 'stable_MT',
                                       'stable_reverse_MT']:
                if (star_2_CO or (TF1 in TF1_POOL_STABLE and
                    ('primary' in TF1 or 'Primary' in TF1))):
                    star = binary.star_1
                    star_i = 1
                    assign_core_collapse_quantities_none(EXTRA_COLUMNS, 2,
                     list(SN_MODELS.keys()))
                elif (TF1 in TF1_POOL_STABLE and
                    ('secondary' in TF1 or 'Secondary' in TF1)):
                    star = binary.star_2
                    star_i = 2
                    assign_core_collapse_quantities_none(EXTRA_COLUMNS, 1,
                     list(SN_MODELS.keys()))
                elif TF1 == 'gamma_center_limit':
                    if (binary.star_1.center_gamma is not None and
                        binary.star_1.center_gamma >= 10.):
                        star = binary.star_1
                        star_i = 1
                        assign_core_collapse_quantities_none(EXTRA_COLUMNS, 2,
                         list(SN_MODELS.keys()))
                    elif (binary.star_2.center_gamma is not None and
                        binary.star_2.center_gamma >= 10.):
                        star = binary.star_2
                        star_i = 2
                        assign_core_collapse_quantities_none(EXTRA_COLUMNS, 1,
                         list(SN_MODELS.keys()))
                    else:
                        assign_core_collapse_quantities_none(EXTRA_COLUMNS, 1,
                         list(SN_MODELS.keys()))
                        assign_core_collapse_quantities_none(EXTRA_COLUMNS, 2,
                         list(SN_MODELS.keys()))
                        Pwarn(f'{grid.MESA_dirs[i]} ended with '
                              'TF1=gamma_center_limit however the star has '
                              'center_gamma < 10. This star cannot go through '
                              'step_SN appending NONE compact object '
                              'properties!', "InappropriateValueWarning")
                        star_i = 0
                else:
                    assign_core_collapse_quantities_none(EXTRA_COLUMNS, 1,
                     list(SN_MODELS.keys()))
                    assign_core_collapse_quantities_none(EXTRA_COLUMNS, 2,
                     list(SN_MODELS.keys()))
                    Pwarn(f'{grid.MESA_dirs[i]} ended with TF={TF1} and '
                          f'IC={interpolation_class}. This star cannot go '
                          'through step_SN appending NONE compact object '
                          'properties!', "InappropriateValueWarning")
                    star_i = 0

                if star_i == 0:
                    pass
                elif star.state in STAR_STATES_CC:
                    if verbose:
                        print_CC_quantities(star)

                    for SN_MODEL_NAME in SN_MODELS.keys():
                        SN_MODEL = get_SN_MODEL(SN_MODEL_NAME)
                        mechanism = SN_MODEL['mechanism']+SN_MODEL['engine']
                        SN_MODEL['allow_spin_None'] = True
                        SN = StepSN(**SN_MODEL)
                        star_copy = copy.copy(star)
                        try:
                            flush = False
                            SN.collapse_star(star_copy)
                            for quantity in CC_quantities:
                                if quantity in ['state', 'SN_type']:
                                    if not isinstance(getattr(star_copy, quantity), str):
                                        flush = True
                                        Pwarn(f'{SN_MODEL_NAME} {mechanism} {quantity} is not a string!', "InappropriateValueWarning")
                                elif quantity != 'CO_interpolation_class':
                                    if quantity in ['spin', 'M4', 'mu4', "h1_mass_ej", "he4_mass_ej"]:
                                        if ((not isinstance(getattr(star_copy, quantity), float))
                                            and (getattr(star_copy, quantity) != None)):
                                            flush = True
                                            Pwarn(f'{SN_MODEL_NAME} {mechanism} {quantity} is not a float nor None!', "InappropriateValueWarning")
                                    elif not isinstance(getattr(star_copy, quantity), float):
                                        flush = True
                                        Pwarn(f'{SN_MODEL_NAME} {mechanism} {quantity} is not a float!', "InappropriateValueWarning")
                        except Exception as e:
                            flush = True
                            if verbose:
                                print('')
                                print(f'Error during {SN_MODEL_NAME} {mechanism} core collapse prescrition!')
                                print(e)
                                print('TF1:', TF1)
                                print('interpolation class:',  interpolation_class)
                                print('run directory:', grid.MESA_dirs[i])
                                print('')
                        if flush:
                            assign_core_collapse_quantities_none(EXTRA_COLUMNS,
                                                                 star_i,
                                                                 SN_MODEL_NAME)
                        else:
                            for quantity in CC_quantities:
                                if quantity != 'CO_interpolation_class':
                                    EXTRA_COLUMNS[f'S{star_i}_{SN_MODEL_NAME}_{quantity}'].append(
                                    getattr(star_copy, quantity))
                                else:
                                    if getattr(star_copy, 'state') == 'BH' and 'case' in TF2 and '1' in TF2 and '2' in TF2:
                                        EXTRA_COLUMNS[f'S{star_i}_{SN_MODEL_NAME}_{quantity}'].append(
                                        getattr(star_copy, 'state')+'_reverse_MT')
                                    else:
                                        EXTRA_COLUMNS[f'S{star_i}_{SN_MODEL_NAME}_{quantity}'].append(
                                        getattr(star_copy, 'state'))
                            if verbose:
                                print_CC_quantities(star_copy, f'{SN_MODEL_NAME}_{mechanism}')
                else:
                    # star not explodable
                    assign_core_collapse_quantities_none(EXTRA_COLUMNS, star_i,
                     list(SN_MODELS.keys()))

            else:
                # inital_RLOF, unstable_MT not_converged
                assign_core_collapse_quantities_none(EXTRA_COLUMNS, 1,
                                                     list(SN_MODELS.keys()))
                assign_core_collapse_quantities_none(EXTRA_COLUMNS, 2,
                                                     list(SN_MODELS.keys()))

        else:
            if star.state in STAR_STATES_CC:
                if verbose:
                    print_CC_quantities(star)

                for SN_MODEL_NAME in SN_MODELS.keys():
                    SN_MODEL = get_SN_MODEL(SN_MODEL_NAME)
                    mechanism = SN_MODEL['mechanism']+SN_MODEL['engine']
                    SN_MODEL['allow_spin_None'] = True
                    SN = StepSN(**SN_MODEL)
                    star_copy = copy.copy(star)
                    try:
                        flush = False
                        SN.collapse_star(star_copy)
                        for quantity in CC_quantities:
                            if quantity in ['state', 'SN_type']:
                                if not isinstance(getattr(star_copy, quantity), str):
                                    flush = True
                                    Pwarn(f'{SN_MODEL_NAME} {mechanism} {quantity} is not a string!', "InappropriateValueWarning")
                            elif quantity != 'CO_interpolation_class':
                                if quantity in ['spin', 'M4', 'mu4', "h1_mass_ej", "he4_mass_ej"]:
                                    if ((not isinstance(getattr(star_copy, quantity), float))
                                        and (getattr(star_copy, quantity) != None)):
                                        flush = True
                                        Pwarn(f'{SN_MODEL_NAME} {mechanism} {quantity} is not a float nor None!', "InappropriateValueWarning")
                                elif not isinstance(getattr(star_copy, quantity), float):
                                    flush = True
                                    Pwarn(f'{SN_MODEL_NAME} {mechanism} {quantity} is not a float!', "InappropriateValueWarning")
                    except Exception as e:
                        flush = True
                        if verbose:
                            print('')
                            print(f'Error during {SN_MODEL_NAME} {mechanism} core collapse prescrition!')
                            print(e)
                            print('TF1:', TF1)
                            print('run directory:', grid.MESA_dirs[i])
                            print('')
                    if flush:
                        assign_core_collapse_quantities_none(EXTRA_COLUMNS, 1,
                                                             SN_MODEL_NAME)
                    else:
                        for quantity in CC_quantities:
                            if quantity != 'CO_interpolation_class':
                                EXTRA_COLUMNS[f'S1_{SN_MODEL_NAME}_{quantity}'].append(
                                 getattr(star_copy, quantity))
                            else:
                                EXTRA_COLUMNS[f'S1_{SN_MODEL_NAME}_{quantity}'].append(
                                 getattr(star_copy, 'state'))
                        if verbose:
                            print_CC_quantities(star_copy, f'{SN_MODEL_NAME}_{mechanism}')
            else:
                assign_core_collapse_quantities_none(EXTRA_COLUMNS, 1,
                                                     list(SN_MODELS.keys()))

        # add MT history column by combining TF1 and TF2
        if not single_star:
            combined_TF12 = combine_TF12([IC], [TF2])
            mt_history = DEFAULT_MARKERS_COLORS_LEGENDS['combined_TF12'][combined_TF12[0]][3]
            EXTRA_COLUMNS['mt_history'].append(mt_history)

        # check dataset completeness
        n_control = len(EXTRA_COLUMNS['S1_state'])
        for key in EXTRA_COLUMNS.keys():
            if n_control != len(EXTRA_COLUMNS[key]): # pragma: no cover
                raise ValueError(
                    '%s has not the correct dimension! Error occoured after '
                    'collapsing binary index=%s' % (key, i))

    # to avoid confusion rename core-collaspe compact object state
    # "SN_MODEL_NAME_state" to "SN_MODEL_NAME_CO_type"
    for SN_MODEL_NAME in SN_MODELS.keys():
        EXTRA_COLUMNS[f'S1_{SN_MODEL_NAME}_CO_type'] = EXTRA_COLUMNS.pop(
                f'S1_{SN_MODEL_NAME}_state')
        if f'S2_{SN_MODEL_NAME}_state' in EXTRA_COLUMNS:
            EXTRA_COLUMNS[f'S2_{SN_MODEL_NAME}_CO_type'] = EXTRA_COLUMNS.pop(
                f'S2_{SN_MODEL_NAME}_state')

    return MESA_dirs, EXTRA_COLUMNS

def add_post_processed_quantities(grid, MESA_dirs_EXTRA_COLUMNS,
                                  EXTRA_COLUMNS):
    """Append post processed quantity to a grid.

    This function appends the quantities computed in post_process_grid to any
    grid. Note that this function ensures you can append the quantities only if
    the grid follows the order of MESA_dirs_EXTRA_COLUMNS.

    Parameters
    ----------
    grid : PSyGrid
        MESA grid in PSyGrid format.
    MESA_dirs_EXTRA_COLUMNS: list
        List containing the path to each run corresponding to the post
        processed values. This is used to ensure one to one mapping when
        appending the extra columns back to a grid.
    EXTRA_COLUMNS: dict
        Dictionary containing all post processed quantities.

    """
    # check correspondance of EXTRA_COLUMNS with grid
    if not np.all(np.array(grid.MESA_dirs)
                  == np.array(MESA_dirs_EXTRA_COLUMNS)):
        raise ValueError(
            'EXTRA_COLUMNS do not follow the correct order of grid!')

    for column in EXTRA_COLUMNS.keys():
        if (("state" in column) or ("type" in column) or ("class" in column)
            or (column == 'mt_history')):
            values = np.asarray(EXTRA_COLUMNS[column], str)
        else:
            values = np.asarray(EXTRA_COLUMNS[column], float)
        grid.add_column(column, values, overwrite=True)
