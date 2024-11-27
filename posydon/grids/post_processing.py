"""Module for post-processing POSYDON grids."""

from posydon.binary_evol.binarystar import BinaryStar
from posydon.binary_evol.singlestar import SingleStar
from posydon.binary_evol.SN.step_SN import StepSN
from posydon.binary_evol.flow_chart import STAR_STATES_CC
from posydon.utils.common_functions import (
    calculate_Patton20_values_at_He_depl,
    CEE_parameters_from_core_abundance_thresholds,
    check_state_of_star)
import numpy as np
from tqdm import tqdm
import copy


__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
]


__credits__ = [
    'Konstantinos Kovlakas <kkovlakas@physics.uoc.gr>'
]


MODEL = {
    "mechanism": None,
    "engine": None,
    "PISN": "Marchant+19",
    "ECSN": "Podsiadlowksi+04",
    "max_neutrino_mass_loss": 0.5,
    "kick": True,
    "sigma_kick_CCSN_NS": 265.0,
    "sigma_kick_CCSN_BH": 265.0,
    "sigma_kick_ECSN": 20.0,
    "max_NS_mass": 2.5,
    "use_interp_values": False,
    "use_profiles": True,
    "use_core_masses": False,
    "approx_at_he_depletion": False,
    "verbose": False,
}


def post_process_grid(grid, index=None, star_2_CO=True, MODEL=MODEL,
                      single_star=False, verbose=False):
    """Compute post processed quantitiy of any grid.

    This function post processes any supported grid and computes:
    - Core collapse quantities for 5 prescriptions given the fiducial POSYDON
    assumption given in MODEL plus:
      A: direct collapse
      B: Fryer+12-rapid
      C: Fryer+12-delayed
      D: Sukhbold+16-engine, N20
      E: Patton&Sukhbold20-engine, N20
      for each prescrition we store the compact object state (WD/NS/BH),
      SN type (WD, ECSN, CCSN, PISN, PPISN), fallback mass fraction f_gb,
      compact object mass and spin.
    - Core masses at he-depletion (used in Patton core collapse)
    - Mass envelopes for common envelope step.

    Parameters
    ----------
    grid : PSyGrid
        MESA grid in PSyGrid format.
    index : None, touple or int
        If None, loop over all indicies otherwise provide a range, e.g. [10,20]
        or a index, e.g. 42.
    star_2_CO : bool
        If 'False' star 2 is not a compact object.
    MODEL : dict
        Core collapse model assumptions.
    single_star : bool
        If `True` the PSyGrid contains single stars.
    verbose : bool
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
    CORE_COLLAPSES = [['direct', ''], ['Fryer+12-rapid', ''],
                      ['Fryer+12-delayed', ''], ['Sukhbold+16-engine', 'N20'],
                      ['Patton&Sukhbold20-engine', 'N20']]

    EXTRA_COLUMNS = {
        # Core collapse qunatities: [state, SN_type, f_fb, mass, spin]
        'S1_direct': [],
        'S1_Fryer+12-rapid': [],
        'S1_Fryer+12-delayed': [],
        'S1_Sukhbold+16-engineN20': [],
        'S1_Patton&Sukhbold20-engineN20': [],
        'S2_direct': [],
        'S2_Fryer+12-rapid': [],
        'S2_Fryer+12-delayed': [],
        'S2_Sukhbold+16-engineN20': [],
        'S2_Patton&Sukhbold20-engineN20': [],
        # core masses at He depletion
        'S1_avg_c_in_c_core_at_He_depletion': [],
        'S2_avg_c_in_c_core_at_He_depletion': [],
        'S1_co_core_mass_at_He_depletion': [],
        'S2_co_core_mass_at_He_depletion': [],
        # common envelope quantities
        'S1_lambda_CE_1cent': [],
        'S1_lambda_CE_10cent': [],
        'S1_lambda_CE_30cent': [],
        'S1_lambda_CE_pure_He_star_10cent': [],
        'S1_m_core_CE_1cent': [],
        'S1_m_core_CE_10cent': [],
        'S1_m_core_CE_30cent': [],
        'S1_m_core_CE_pure_He_star_10cent': [],
        'S1_r_core_CE_1cent': [],
        'S1_r_core_CE_10cent': [],
        'S1_r_core_CE_30cent': [],
        'S1_r_core_CE_pure_He_star_10cent': [],
        'S2_lambda_CE_1cent': [],
        'S2_lambda_CE_10cent': [],
        'S2_lambda_CE_30cent': [],
        'S2_lambda_CE_pure_He_star_10cent': [],
        'S2_m_core_CE_1cent': [],
        'S2_m_core_CE_10cent': [],
        'S2_m_core_CE_30cent': [],
        'S2_m_core_CE_pure_He_star_10cent': [],
        'S2_r_core_CE_1cent': [],
        'S2_r_core_CE_10cent': [],
        'S2_r_core_CE_30cent': [],
        'S2_r_core_CE_pure_He_star_10cent': [],
        # stellar states
        'S1_state': [],
        'S2_state': [],
        # composition
        'S1_surface_other': [],
        'S1_center_other': [],
        'S2_surface_other': [],
        'S2_center_other': [],
    }

    # remove star 2 columns in case of single star grid
    if single_star:
        for key in list(EXTRA_COLUMNS.keys()):
            if 'S2' in key:
                del EXTRA_COLUMNS[key]

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

    for i in tqdm(indicies):

        if not single_star:
            # initilise binary
            binary = BinaryStar.from_run(grid[i], history=True, profiles=True)
            stars = [binary.star_1, binary.star_2]
            stars_CO = [False, star_2_CO]
            interpolation_class = grid.final_values['interpolation_class'][i]
            IC = grid.final_values['interpolation_class'][i]
        else:
            star = SingleStar.from_run(grid[i], history=True, profile=True)
            stars = [star]
            stars_CO = [False]
            IC = 'no_MT'
        TF1 = grid.final_values['termination_flag_1'][i]

        # He depeltion and CE quantities
        for j, star in enumerate(stars):
            if not stars_CO[j] and IC in ['no_MT', 'stable_MT', 'unstable_MT']:
                EXTRA_COLUMNS['S%s_state' % (j+1)].append(check_state_of_star(
                    star, star_CO=False))
                calculate_Patton20_values_at_He_depl(star)
                EXTRA_COLUMNS['S%s_avg_c_in_c_core_at_He_depletion' % (j+1)].append(star.avg_c_in_c_core_at_He_depletion)
                EXTRA_COLUMNS['S%s_co_core_mass_at_He_depletion' % (j+1)].append(star.co_core_mass_at_He_depletion)
                CEE_parameters_from_core_abundance_thresholds(star)
                EXTRA_COLUMNS['S%s_m_core_CE_1cent' % (j+1)].append(star.m_core_CE_1cent)
                EXTRA_COLUMNS['S%s_m_core_CE_10cent' % (j+1)].append(star.m_core_CE_10cent)
                EXTRA_COLUMNS['S%s_m_core_CE_30cent' % (j+1)].append(star.m_core_CE_30cent)
                EXTRA_COLUMNS['S%s_m_core_CE_pure_He_star_10cent' % (j+1)].append(star.m_core_CE_pure_He_star_10cent)
                EXTRA_COLUMNS['S%s_r_core_CE_1cent' % (j+1)].append(star.r_core_CE_1cent)
                EXTRA_COLUMNS['S%s_r_core_CE_10cent' % (j+1)].append(star.r_core_CE_10cent)
                EXTRA_COLUMNS['S%s_r_core_CE_30cent' % (j+1)].append(star.r_core_CE_30cent)
                EXTRA_COLUMNS['S%s_r_core_CE_pure_He_star_10cent' % (j+1)].append(star.r_core_CE_pure_He_star_10cent)
                EXTRA_COLUMNS['S%s_lambda_CE_1cent' % (j+1)].append(star.lambda_CE_1cent)
                EXTRA_COLUMNS['S%s_lambda_CE_10cent' % (j+1)].append(star.lambda_CE_10cent)
                EXTRA_COLUMNS['S%s_lambda_CE_30cent' % (j+1)].append(star.lambda_CE_30cent)
                EXTRA_COLUMNS['S%s_lambda_CE_pure_He_star_10cent' % (j+1)].append(star.lambda_CE_pure_He_star_10cent)
                s_o = 1. - star.surface_h1 - star.surface_he4 - star.surface_c12 - star.surface_n14 - star.surface_o16
                c_o = 1. - star.center_h1 - star.center_he4 - star.center_c12 - star.center_n14 - star.center_o16
                EXTRA_COLUMNS['S%s_surface_other' % (j+1)].append(s_o)
                EXTRA_COLUMNS['S%s_center_other' % (j+1)].append(c_o)
            else:
                if IC == 'initial_MT':
                    EXTRA_COLUMNS['S%s_state' % (j+1)].append(None)
                else:
                    # CO states are classified and used in mesa step
                    EXTRA_COLUMNS['S%s_state' % (j+1)].append(check_state_of_star(star, star_CO=True))
                EXTRA_COLUMNS['S%s_avg_c_in_c_core_at_He_depletion' % (j+1)].append(None)
                EXTRA_COLUMNS['S%s_co_core_mass_at_He_depletion' % (j+1)].append(None)
                EXTRA_COLUMNS['S%s_m_core_CE_1cent' % (j+1)].append(None)
                EXTRA_COLUMNS['S%s_m_core_CE_10cent' % (j+1)].append(None)
                EXTRA_COLUMNS['S%s_m_core_CE_30cent' % (j+1)].append(None)
                EXTRA_COLUMNS['S%s_m_core_CE_pure_He_star_10cent' % (j+1)].append(None)
                EXTRA_COLUMNS['S%s_r_core_CE_1cent' % (j+1)].append(None)
                EXTRA_COLUMNS['S%s_r_core_CE_10cent' % (j+1)].append(None)
                EXTRA_COLUMNS['S%s_r_core_CE_30cent' % (j+1)].append(None)
                EXTRA_COLUMNS['S%s_r_core_CE_pure_He_star_10cent' % (j+1)].append(None)
                EXTRA_COLUMNS['S%s_lambda_CE_1cent' % (j+1)].append(None)
                EXTRA_COLUMNS['S%s_lambda_CE_10cent' % (j+1)].append(None)
                EXTRA_COLUMNS['S%s_lambda_CE_30cent' % (j+1)].append(None)
                EXTRA_COLUMNS['S%s_lambda_CE_pure_He_star_10cent' % (j+1)].append(None)
                EXTRA_COLUMNS['S%s_surface_other' % (j+1)].append(None)
                EXTRA_COLUMNS['S%s_center_other' % (j+1)].append(None)

        if not single_star:
            # core collapse quantities
            if interpolation_class in ['no_MT', 'stable_MT']:
                if star_2_CO or TF1 == 'Primary has depleted central carbon':
                    star = binary.star_1
                    s = 'S1_'
                    for m in CORE_COLLAPSES:
                        EXTRA_COLUMNS['S2_'+m[0]+m[1]].append([None]*5)
                elif TF1 == 'Secondary has depleted central carbon':
                    star = binary.star_2
                    s = 'S2_'
                    for m in CORE_COLLAPSES:
                        EXTRA_COLUMNS['S1_'+m[0]+m[1]].append([None]*5)
                elif TF1 == 'gamma_center_limit':
                    if binary.star_1.center_gamma >= 10.:
                        star = binary.star_1
                        s = 'S1_'
                        for m in CORE_COLLAPSES:
                            EXTRA_COLUMNS['S2_'+m[0]+m[1]].append([None]*5)
                    elif binary.star_2.center_gamma >= 10.:
                        star = binary.star_2
                        s = 'S2_'
                        for m in CORE_COLLAPSES:
                            EXTRA_COLUMNS['S1_'+m[0]+m[1]].append([None]*5)
                    else:
                        raise ValueError('No star has center_gamma < 10')
                else:
                    raise ValueError('TF1 = %s not supported' % TF1)

                if verbose:
                    print("{:<30} {:<33} {:12} {:10} {:15} {:10}".format(
                        "mechanism", "state", "SN type", "f_fb",
                        "mass [Msun]", "spin"))
                    print('')
                    print("{:<30} {:<33} {:10} {:10} {:7.2f} {:14.2f}".format(
                        'PRE SN STAR', star.state, '',
                        '', star.mass, star.spin))
                    print('')

                for m in CORE_COLLAPSES:
                    MODEL["mechanism"] = m[0]
                    MODEL["engine"] = m[1]
                    SN = StepSN(**MODEL)
                    try:
                        star_copy = copy.copy(star)
                        SN.collapse_star(star_copy)
                        EXTRA_COLUMNS[s+m[0]+m[1]].append(
                            [star_copy.state, star_copy.SN_type,
                             star_copy.f_fb, star_copy.mass, star_copy.spin])
                        if verbose:
                            print(
                                "{:<30} {:<33} {:10} {:6.2f} {:11.2f} {:14.2f}"
                                .format(MODEL["mechanism"], star_copy.state,
                                        star_copy.SN_type, star_copy.f_fb,
                                        star_copy.mass, star_copy.spin))
                    except Exception as e:
                        EXTRA_COLUMNS[s+m[0]+m[1]].append([None]*5)

                        if verbose:
                            print('')
                            print('Error during %s core collapse prescrition!'
                                  % m[0])
                            print(e)
                            print('TF1', TF1)
                            print('interpolation class',  interpolation_class)

            else:    # inital_RLOF, unstable_MT not_convergedd
                if TF1 == 'Primary has depleted central carbon':
                    raise ValueError(
                        'Primary reached carbon depletion but was not '
                        'collapsed! This should not happen!')
                if TF1 == 'Secondary has depleted central carbon':
                    raise ValueError(
                        'Secondary reached carbon depletion but was not '
                        'collapsed! This should not happen!')

                for m in CORE_COLLAPSES:
                    EXTRA_COLUMNS['S1_'+m[0]+m[1]].append([None]*5)
                    EXTRA_COLUMNS['S2_'+m[0]+m[1]].append([None]*5)

        else:
            if star.state in STAR_STATES_CC:
                if verbose:
                    print("{:<30} {:<33} {:12} {:10} {:15} {:10}".format(
                        "mechanism", "state", "SN type", "f_fb",
                        "mass [Msun]", "spin"))
                    print('')
                    print("{:<30} {:<33} {:10} {:10} {:7.2f} {:14.2f}".format(
                        'PRE SN STAR', star.state,
                        '', '', star.mass, star.spin))
                    print('')

                for m in CORE_COLLAPSES:
                    MODEL["mechanism"] = m[0]
                    MODEL["engine"] = m[1]
                    SN = StepSN(**MODEL)
                    try:
                        star_copy = copy.copy(star)
                        SN.collapse_star(star_copy)
                        EXTRA_COLUMNS['S1_'+m[0]+m[1]].append([
                            star_copy.state, star_copy.SN_type, star_copy.f_fb,
                            star_copy.mass, star_copy.spin])
                        if verbose:
                            print(
                                "{:<30} {:<33} {:10} {:6.2f} {:11.2f} {:14.2f}"
                                .format(MODEL["mechanism"], star_copy.state,
                                        star_copy.SN_type, star_copy.f_fb,
                                        star_copy.mass, star_copy.spin))
                    except Exception as e:
                        EXTRA_COLUMNS[s+m[0]+m[1]].append([None]*5)

                        if verbose:
                            print('')
                            print('Error during %s core collapse prescrition!'
                                  % m[0])
                            print(e)
                            print('TF1', TF1)
                            print('interpolation class',  interpolation_class)
            else:
                for m in CORE_COLLAPSES:
                    EXTRA_COLUMNS['S1_'+m[0]+m[1]].append([None]*5)

        # check dataset completeness
        n_control = len(EXTRA_COLUMNS['S1_direct'])
        for key in EXTRA_COLUMNS.keys():
            if n_control != len(EXTRA_COLUMNS[key]):
                raise ValueError(
                    '%s has not the correct dimension! Error occoured after '
                    'collapsing binary index=%s' % (key, i))

    return MESA_dirs, EXTRA_COLUMNS


def add_post_processed_quantities(grid, MESA_dirs_EXTRA_COLUMNS, EXTRA_COLUMNS,
                                  verbose=False):
    """Append post processed quantity to a grid.

    This function appends the quantities computed in post_process_grid to any
    grid. Note that this function ensures you can append the quantities only if
    the grid follows the order of MESA_dirs_EXTRA_COLUMNS.

    Parameters
    ----------
    grid : PSyGrid
        MESA grid in PSyGrid format.
    MESA_dirs: list
        List containing the path to each run corresponding to the post
        processed values. This is used to ensure one to one mapping when
        appending the extra columns back to a grid.
    EXTRA_COLUMNS: dict
        Dictionary containing all post processed quantities.
    verbose : bool
        If `True` print the results of each core collapse on screen.

    """
    # check correspondance of EXTRA_COLUMNS with grid
    if not np.all(np.array(grid.MESA_dirs)
                  == np.array(MESA_dirs_EXTRA_COLUMNS)):
        raise ValueError(
            'EXTRA_COLUMNS do not follow the correct order of grid!')

    CC_keys = ['state', 'SN_type', 'f_fb', 'mass', 'spin']

    for key in EXTRA_COLUMNS.keys():
        if ('direct' in key or 'Fryer+12-rapid' in key
                or 'Fryer+12-delayed' in key or 'Sukhbold+16-engineN20' in key
                or 'Patton&Sukhbold20-engineN20' in key):
            for j, key_CC in enumerate(CC_keys):
                column = '%s_%s' % (key, key_CC)
                values = []
                for i in range(len(EXTRA_COLUMNS[key])):
                    values.append(EXTRA_COLUMNS[key][i][j])
                if "state" in column or "SN_type" in column:
                    values = np.asarray(values, str)
                else:
                    values = np.asarray(values, float)
                grid.add_column(column, values, overwrite=True)
        else:
            column = key
            values = []
            for i in range(len(EXTRA_COLUMNS[key])):
                values.append(EXTRA_COLUMNS[key][i])
            if "state" in column:
                values = np.asarray(values, str)
            else:
                values = np.asarray(values, float)
            grid.add_column(column, values, overwrite=True)
