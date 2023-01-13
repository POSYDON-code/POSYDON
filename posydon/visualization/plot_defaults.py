"""Default values for plotting routines.

This file contains the default parameters of each plotting routines.
PLOT_PROPERTIES:
    Each default plotting variable is contained in this dictionary.
DEFAULT_MARKERS_COLORS_LEGENDS:
    Each termination condition is associated with a marker shape, size, color
    and label.
"""


__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Devina Misra <devina.misra@unige.ch>",
    "Scott Coughlin <scottcoughlin2014@u.northwestern.edu>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
]


from matplotlib import checkdep_usetex


PLOT_PROPERTIES = {
    'show_fig': False,
    'close_fig': True,
    'path_to_file': './',
    'fname': None,
    'figsize': (3.38, 3.38),
    'bbox_inches': 'tight',
    'dpi': 300,
    'log10_x': False,
    'log10_y': False,
    'log10_z': False,
    'xmin': None,
    'xmax': None,
    'ymin': None,
    'ymax': None,
    'zmin': None,
    'zmax': None,
    'title': None,
    'rcParams': {"text.usetex": checkdep_usetex(True),
                 "font.family": "serif",
                 "font.sans-serif": ["Computer Modern Roman"]},
    'title_font_dict': {'fontsize': 10},
    'title_loc': 'center',
    'xlabel': None,
    'xlabel_kwargs': {
        'fontsize': 10
    },
    'ylabel': None,
    'ylabel_kwargs': {
        'fontsize': 10
    },
    'marker_size': 4,
    'hspace': None,
    'wspace': None,
    'const_R_lines': False,
    'colorbar': {
        'label': None,
        'label_size': 10,
        'orientation': 'horizontal',
        'fraction': 0.15,
        'pad': 0.2,
        'shrink': 1,
        'aspect': 20,
        'anchor': (0.0, 0.5),
        'panchor': (1.0, 0.5),
        'extend': 'neither'
    },
    'legend1D': {
        'title': None,
        'lines_legend': None,
        'title_font_size': 8,
        'loc': 'upper right',
        'ncol': 1,
        'borderaxespad': None,
        'handletextpad': None,
        'columnspacing': 0.4,
        'prop': {
            'size': 8
        },
        'shrink_box': 1,
        'bbox_to_anchor': None
    },
    'legend2D': {
        'title': 'Termination flags',
        'title_font_size': 8,
        'loc': 'center left',
        'ncol': 1,
        'borderaxespad': None,
        'handletextpad': None,
        'columnspacing': 0.4,
        'prop': {
            'size': 8
        },
        'shrink_box': 0.85,
        'bbox_to_anchor': (1, 0.5)
    }
}

list_of_colors = ['#a6611a',
                  '#dfc27d',
                  [(31/255, 119/255, 180/255)],
                  [(255/255, 127/255, 14/255)]]

TF1_label_stable = 'Reached end life'
TF1_label_initial = 'Initial RLOF'
TF1_label_unstable = 'Unstable RLOF'
color_unstable = 'black'

DEFAULT_MARKERS_COLORS_LEGENDS = {
    'termination_flag_1': {
        'terminate due to primary depleting carbon (inverse sn?)':
            ['s', 2, None, TF1_label_stable],
        'Primary has depleted central carbon':
            ['s', 2, None, TF1_label_stable],
        'Secondary has depleted central carbon':
            ['o', 2, None, TF1_label_stable],
        'offcenter neon ignition for primary':
            ['s', 2, None, TF1_label_stable],
        'offcenter neon ignition for secondary':
            ['o', 2, None, TF1_label_stable],
        'forced_initial_RLO':
            ['.', 1, color_unstable, TF1_label_initial],
        'overflow from L1 at ZAMS':
            ['.', 1, color_unstable, TF1_label_initial],
        'Terminate because of overflowing initial model':
            ['.', 1, color_unstable, TF1_label_initial],
        'overflow from L2 point for q>1 at ZAMS':
            ['.', 1, color_unstable, TF1_label_initial],
        'overflow from L2 surface for q<1 at ZAMS':
            ['.', 1, color_unstable, TF1_label_initial],
        'overflow from L2 surface for q>1 at ZAMS':
            ['.', 1, color_unstable, TF1_label_initial],
        'overflow from L2 surface for q<1':
            ['D', 1, color_unstable, TF1_label_unstable],
        r'overflow from L2 (D_L2) distance for q(=Macc/Mdon)>1, '
        'donor is star 1':
            ['D', 1, color_unstable, TF1_label_unstable],
        r'overflow from L2 (R_L2) surface for q(=Macc/Mdon)<1, '
        'donor is star 1':
            ['D', 1, color_unstable, TF1_label_unstable],
        r'overflow from L2 (R_L2) surface for q(=Macc/Mdon)>1, '
        'donor is star 1':
            ['D', 1, color_unstable, TF1_label_unstable],
        r'overflow from L2 (D_L2) distance for q(=Macc/Mdon)<1, '
        'donor is star 1':
            ['D', 1, color_unstable, TF1_label_unstable],
        'overflow from L2 (R_L2) surface for q(=Macc/Mdon)<1, donor is star 2':
            ['D', 1, color_unstable, TF1_label_unstable],
        'reached maximum mass transfer rate: 10.0d0':
            ['D', 1, color_unstable, TF1_label_unstable],
        'Reached maximum mass transfer rate: 1d-1':
            ['D', 1, color_unstable, TF1_label_unstable],
        'Reached the critical mt rate':
            ['D', 1, color_unstable, TF1_label_unstable],
        'Reached TPAGB':
            ['s', 2, None, TF1_label_initial],
        'Both stars fill their Roche Lobe and at least one of them is off MS':
            ['D', 1, color_unstable, TF1_label_unstable],
        'Terminate due to L2 overflow during case A':
            ['D', 1, color_unstable, TF1_label_unstable],
        'Reached maximum mass transfer rate: Exceeded photon trapping radius':
            ['D', 1, color_unstable, TF1_label_unstable],
        'Terminate because accretor (r-rl)/rl > accretor_overflow_terminate':
            ['D', 1, color_unstable, TF1_label_unstable],
        'logQ_limit':
            ['x', 1, 'red', 'logQ_limit'],
        'logQ_min_limit':
            ['x', 1, 'red', 'logQ_limit'],
        'min_timestep_limit':
            ['x', 1, 'red', 'Not converged'],
        'reach cluster timelimit':
            ['x', 1, 'red', 'Not converged'],
        'no termination code':
            ['x', 1, 'red', 'no termination code'],
        'envelope_mass_limit':
            ['s', 2, None, TF1_label_stable],
        'gamma_center_limit':
            ['s', 2, None, TF1_label_stable],
        'max_age':
            ['s', 2, None, TF1_label_stable],
        'Initial RLOF':
            ['.', 1, 'black', TF1_label_initial],
        'Not converged':
            ['x', 1, 'red', 'Not converged'],
        'ignored_no_BH':
            ['.', 1, color_unstable, TF1_label_initial],
        'ignored_no_RLO':
            ['.', 1, color_unstable, TF1_label_initial],
        'unknown':
            ['+', 1, 'red', 'unknown'],

    },

    'termination_flag_2': {

        'initial_RLOF':
            ['+', 1, 'black', 'initial RLOF'],
        'forced_initial_RLO':
            ['+', 1, 'black', 'initial RLOF'],
        'ignored_no_BH':
            ['+', 1, 'black', 'initial RLOF'],
        'ignored_no_RLO':
            ['+', 1, 'black', 'initial RLOF'],
        'no_RLOF':
            ['.', 1, 'black', 'no RLOF'],
        'contact_during_MS':
            ['v', 1, 'black', 'contact during MS'],
        'L2_RLOF':
            ['^', 1, 'black', 'L2 RLOF'],

        'caseA_from_star1':
            ['s', 2, 'tab:blue', 'case A from star1'],
        'caseA/B_from_star1':
            ['s', 2, 'tab:green', 'case A/B from star1'],
        'caseA/B/A_from_star1':
            ['s', 2, 'lightgrey', 'case A/B/A from star1'],
        'caseA/B/C_from_star1':
            ['s', 2, 'yellow', 'case A/B/C from star1'],
        'caseA/C_from_star1':
            ['s', 2, 'yellow', 'case A/C from star1'],
        'caseA/B/BB_from_star1':
            ['s', 2, 'tab:red', 'case A/B/BB from star1'],
        'caseA/B/C/BB_from_star1':
            ['s', 2, 'tab:cyan', 'case A/B/C/BB from star1'],
        'caseB_from_star1':
            ['s', 2, 'tab:purple', 'case B from star1'],
        'caseB/A_from_star1':
            ['s', 2, 'tab:olive', 'case B/A from star1'],
        'caseB/BB_from_star1':
            ['s', 2, 'tab:pink', 'case B/BB from star1'],
        'caseB/C/BB_from_star1':
            ['s', 2, 'tab:gray', 'case B/C/BB from star1'],
        'caseB/C_from_star1':
            ['s', 2, 'tab:orange', 'case B/C from star1'],
        'caseC_from_star1':
            ['s', 2, 'black', 'case C from star1'],
        'caseC/BB_from_star1':
            ['s', 2, 'brown', 'case C/BB from star1'],

        'caseA_from_star2':
            ['o', 2, 'tab:blue', 'case A from star2'],
        'caseA/B_from_star2':
            ['o', 2, 'tab:green', 'case A/B from star2'],
        'caseA/B/A_from_star2':
            ['s', 2, 'lightgrey', 'case A/B/A from star2'],
        'caseA/B/C_from_star2':
            ['s', 2, 'yellow', 'case A/B/C from star2'],
        'caseA/B/BB_from_star2':
            ['o', 2, 'tab:red', 'case A/B/BB from star2'],
        'caseA/B/C/BB_from_star2':
            ['o', 2, 'tab:cyan', 'case A/B/C/BB from star2'],
        'caseB_from_star2':
            ['o', 2, 'tab:purple', 'case B from star2'],
        'caseB/A_from_star2':
            ['o', 2, 'tab:olive', 'case B/A from star2'],
        'caseB/BB_from_star2':
            ['o', 2, 'tab:pink', 'case B/BB from star2'],
        'caseB/C/BB_from_star2':
            ['o', 2, 'tab:gray', 'case B/C/BB from star2'],
        'caseB/C_from_star2':
            ['o', 2, 'tab:orange', 'case B/C from star2'],
        'caseC_from_star2':
            ['o', 2, 'black', 'case C from star2'],
        'caseC/BB_from_star2':
            ['s', 2, 'brown', 'case C/BB from star2'],

        'caseBA_from_star1':
            ['s', 2, 'tab:blue', 'case BA from star1'],
        'caseBB_from_star1':
            ['s', 2, 'tab:green', 'case BB from star1'],
        'caseBA/BB_from_star1':
            ['s', 2, 'tab:red', 'case BA/BB from star1'],

        '_from_star1':
            ['x', 1, 'black', 'unknown'],
        '_from_star2':
            ['x', 1, 'black', 'unknown'],

    },

    'termination_flag_3': {
        'H-rich_non_burning':
            ['v', 2, 'tab:orange', 'H-rich H non-buring'],
        'H-rich_Core_H_burning':
            ['s', 2, 'tab:olive', 'H-rich core H buring'],
        'H-rich_Shell_H_burning':
            ['s', 2, 'tab:red', 'H-rich shell H buring'],
        'H-rich_Core_C_burning':
            ['s', 2, 'tab:pink', 'H-rich core C burning'],
        'H-rich_Central_C_depletion':
            ['s', 2, 'tab:brown', 'H-rich central C depletion'],
        'H-rich_Core_He_burning':
            ['s', 2, 'tab:blue', 'H-rich core He burning'],
        'H-rich_Central_He_depleted':
            ['s', 2, 'tab:green', 'H-rich shell He burning'],
        'stripped_He_Core_He_burning':
            ['o', 2, 'tab:blue', 'stripped He-star core He burning'],
        'stripped_He_Central_He_depleted':
            ['o', 2, 'tab:green', 'stripped He-star shell He burning'],
        'stripped_He_Core_C_burning':
            ['o', 2, 'tab:pink', 'stripped He-star core C burning'],
        'stripped_He_Central_C_depletion':
            ['o', 2, 'tab:brown', 'stripped He-star C depletion'],
        'stripped_He_non_burning':
            ['o', 2, 'gray', 'stripped He-star non burning'],
        'stripped_He_Core_H_burning':
            ['x', 1, 'black', 'unknown'],
        'stripped_He_Shell_H_burning':
            ['x', 1, 'black', 'unknown'],
        'undetermined_evolutionary_state':
            ['x', 1, 'black', 'unknown'],
        'BH':
            ['*', 1, 'black', 'BH'],
        'NS':
            ['*', 1, 'tab:gray', 'NS'],
        'ignored_no_BH':
            ['s', 2, 'tab:olive', 'H-rich core H buring'],
        'ignored_no_RLO':
            ['s', 2, 'tab:olive', 'H-rich core H buring'],
    },
    'termination_flag_4': {
        'H-rich_non_burning':
            ['v', 2, 'tab:orange', 'H-rich H non-buring'],
        'H-rich_Core_H_burning':
            ['s', 2, 'tab:olive', 'H-rich core H buring'],
        'H-rich_Shell_H_burning':
            ['s', 2, 'tab:red', 'H-rich shell H buring'],
        'H-rich_Core_C_burning':
            ['s', 2, 'tab:pink', 'H-rich core C burning'],
        'H-rich_Central_C_depletion':
            ['s', 2, 'tab:brown', 'H-rich central C depletion'],
        'H-rich_Core_He_burning':
            ['s', 2, 'tab:blue', 'H-rich core He burning'],
        'H-rich_Near_Central_C_depletion':
            ['s', 2, 'tab:purple', 'H-rich near central C depletion'],
        'H-rich_Central_He_depleted':
            ['s', 2, 'tab:green', 'H-rich shell He burning'],
        'stripped_He_Core_He_burning':
            ['o', 2, 'tab:blue', 'stripped He-star core He burning'],
        'stripped_He_Central_He_depleted':
            ['o', 2, 'tab:green', 'stripped He-star shell He burning'],
        'stripped_He_Core_C_burning':
            ['o', 2, 'tab:pink', 'stripped He-star core C burning'],
        'stripped_He_Central_C_depletion':
            ['o', 2, 'tab:brown', 'stripped He-star C depletion'],
        'stripped_He_non_burning':
            ['o', 2, 'gray', 'stripped He-star non burning'],
        'stripped_He_Core_H_burning':
            ['x', 1, 'black', 'unknown'],
        'stripped_He_Shell_H_burning':
            ['x', 1, 'black', 'unknown'],
        'undetermined_evolutionary_state':
            ['x', 1, 'black', 'unknown'],
        'BH':
            ['*', 1, 'black', 'BH'],
        'NS':
            ['*', 1, 'tab:gray', 'NS'],
        'ignored_no_BH':
            ['s', 2, 'tab:olive', 'H-rich core H buring'],
        'ignored_no_RLO':
            ['s', 2, 'tab:olive', 'H-rich core H buring'],
    },

    'combined_TF12': {
        'Stable contact':
            ['s', 2, list_of_colors[3], 'Stable contact phase'],
        'Stable case A':
            ['s', 2, list_of_colors[2], 'Stable RLOF during MS'],
        'Stable case AB':
            ['s', 2, list_of_colors[1], 'Stable RLOF during postMS'],
        'Stable case ABB':
            ['s', 2, list_of_colors[0], 'Stable RLOF during stripped He star'],
        'Stable case B':
            ['s', 2, list_of_colors[1], 'Stable RLOF during postMS'],
        'Stable case C':
            ['s', 2, list_of_colors[1], 'Stable RLOF during postMS'],
        'Stable case BA':
            ['s', 2, list_of_colors[0], 'Stable RLOF during stripped He star'],
        'Stable case BB':
            ['s', 2, list_of_colors[0], 'Stable RLOF during stripped He star'],
        'Stable case BC':
            ['s', 2, list_of_colors[1], 'Stable RLOF during postMS'],
        'Unstable contact':
            ['D', 1, list_of_colors[3], 'Unstable contact phase'],
        'Unstable case A':
            ['D', 1, list_of_colors[2], 'Unstable RLOF during MS'],
        'Unstable case AB':
            ['D', 1, list_of_colors[1], 'Unstable RLOF during postMS'],
        'Unstable case ABB':
            ['D', 1, list_of_colors[0],
             'Unstable RLOF during stripped He star'],
        'Unstable case B':
            ['D', 1, list_of_colors[1], 'Unstable RLOF during postMS'],
        'Unstable case C':
            ['D', 1, list_of_colors[1], 'Unstable RLOF during postMS'],
        'Unstable case BA':
            ['D', 1, list_of_colors[0],
             'Unstable RLOF during stripped He star'],
        'Unstable case BB':
            ['D', 1, list_of_colors[0],
             'Unstable RLOF during stripped He star'],
        'Unstable case BC':
            ['D', 1, list_of_colors[1], 'Unstable RLOF during postMS'],
        'Unstable L2 RLOF':
            ['D', 1, list_of_colors[0],
             'Unstable RLOF during stripped He star'],
        'Initial RLOF':
            ['.', 1, 'black', 'Initial RLOF'],
        'no_RLOF':
            ['s', 2, 'lightgrey', 'no RLOF'],
        'Not converged':
            ['x', 1, 'red', 'Not converged'],
        'unknown':
            ['+', 1, 'green', 'unknown'],
        },
    'debug': {
        'terminate due to primary depleting carbon (inverse sn?)':
            ['s', 2, None, TF1_label_stable],
        'Primary has depleted central carbon':
            ['s', 2, None, TF1_label_stable],
        'Secondary has depleted central carbon':
            ['o', 2, None, TF1_label_stable],
        'offcenter neon ignition for primary':
            ['s', 2, None, TF1_label_stable],
        'offcenter neon ignition for secondary':
            ['o', 2, None, TF1_label_stable],
        'overflow from L1 at ZAMS':
            ['.', 1.5, None, TF1_label_initial],
        'Terminate because of overflowing initial model':
            ['.', 1.5, None, TF1_label_initial],
        'overflow from L2 point for q>1 at ZAMS':
            ['.', 1.5, None, TF1_label_initial],
        'overflow from L2 surface for q<1 at ZAMS':
            ['.', 1.5, None, TF1_label_initial],
        'overflow from L2 surface for q>1 at ZAMS':
            ['.', 1.5, None, TF1_label_initial],
        'overflow from L2 surface for q<1':
            ['D', 1, None, TF1_label_unstable],
        r'overflow from L2 (D_L2) distance for q(=Macc/Mdon)>1, '
        'donor is star 1':
            ['D', 1, None, TF1_label_unstable],
        r'overflow from L2 (R_L2) surface for q(=Macc/Mdon)<1, '
        'donor is star 1':
            ['D', 1, None, TF1_label_unstable],
        r'overflow from L2 (R_L2) surface for q(=Macc/Mdon)>1, '
        'donor is star 1':
            ['D', 1, None, TF1_label_unstable],
        r'overflow from L2 (D_L2) distance for q(=Macc/Mdon)<1, '
        'donor is star 1':
            ['D', 1, None, TF1_label_unstable],
        'overflow from L2 (R_L2) surface for q(=Macc/Mdon)<1, '
        'donor is star 2':
            ['D', 1, None, TF1_label_unstable],
        'reached maximum mass transfer rate: 10.0d0':
            ['D', 1, None, TF1_label_unstable],
        'Reached maximum mass transfer rate: 1d-1':
            ['D', 1, None, TF1_label_unstable],
        'Reached the critical mt rate':
            ['D', 1, None, TF1_label_unstable],
        'Reached TPAGB':
            ['s', 2, None, TF1_label_initial],
        'Both stars fill their Roche Lobe and at least one of them is off MS':
            ['D', 1, None, TF1_label_unstable],
        'Terminate due to L2 overflow during case A':
            ['D', 1, None, TF1_label_unstable],
        'Reached maximum mass transfer rate: Exceeded photon trapping radius':
            ['D', 1, None, TF1_label_unstable],
        'logQ_limit':
            ['x', 1, None, 'logQ_limit'],
        'logQ_min_limit':
            ['x', 1, None, 'logQ_limit'],
        'min_timestep_limit':
            ['x', 1, None, 'Not converged'],
        'reach cluster timelimit':
            ['x', 1, None, 'Not converged'],
        'no termination code':
            ['x', 1, None, 'no termination code'],
        'envelope_mass_limit':
            ['s', 2, None, TF1_label_stable],
        'gamma_center_limit':
            ['s', 2, None, TF1_label_stable],
        'max_age':
            ['s', 2, None, TF1_label_stable],
        'Initial RLOF':
            ['.', 1.5, 'black', TF1_label_initial],
        'Not converged':
            ['x', 1, None, 'Not converged'],
        'ignored_no_BH':
            ['.', 1.5, color_unstable, TF1_label_initial],
        'ignored_no_RLO':
            ['.', 1.5, color_unstable, TF1_label_initial],
        'unknown':
            ['+', 1, 'green', 'unknown'],
    },
    'interpolation_class': {
        'initial_MT':
            ['o', 2, 'tab:blue', 'initial_MT'],
        'no_MT':
            ['o', 2, 'tab:pink', 'no_MT'],
        'not_converged':
            ['o', 2, 'tab:red', 'not_converged'],
        'stable_MT':
            ['o', 2, 'tab:orange', 'stable_MT'],
        'unstable_MT':
            ['o', 2, 'tab:purple', 'unstable_MT']
    },
    'SN_type': {
        'CCSN':
            ['o', 2, 'tab:blue', 'CCSN'],
        'ECSN':
            ['o', 2, 'tab:orange', 'ECSN'],
        'PPISN':
            ['o', 2, 'tab:red', 'PPISN'],
        'WD':
            ['o', 2, 'tab:purple', 'WD'],
        'None':
            ['o', 2, 'black', 'intial MT / unstable MT / not converged'],
    },
    'state': {
        'BH':
            ['o', 2, 'tab:blue', 'BH'],
        'NS':
            ['o', 2, 'tab:orange', 'NS'],
        'WD':
            ['o', 2, 'tab:purple', 'WD'],
        'None':
            ['o', 2, 'black', 'intial MT / unstable MT / not converged'],
    }
}


DEFAULT_LABELS = {
    # extra
    'mass_ratio':
        [r'$q$', r'$\log_{10}(q)$'],

    # history1/history2
    'star_age':
        [r'$t \, [\mathrm{yr}]$', r'$\log_{10}(t / \mathrm{yr})$'],
    'star_mass':
        [r'$M_\mathrm{\star} \, [M_\odot]$',
         r'$\log_{10}(M_\mathrm{\star} / M_\odot)$'],
    'he_core_mass':
        [r'$M_\mathrm{He-core} \, [M_\odot]$',
         r'$\log_{10}(M_\mathrm{He-core} / M_\odot)$'],
    'c_core_mass':
        [r'$M_\mathrm{C-core} \, [M_\odot]$',
         r'$\log_{10}(M_\mathrm{C-core} / M_\odot)$'],
    'o_core_mass':
        [r'$M_\mathrm{O-core} \, [M_\odot]$',
         r'$\log_{10}(M_\mathrm{O-core} / M_\odot)$'],
    'co_core_mass':
        [r'$M_\mathrm{CO-core} \, [M_\odot]$',
         r'$\log_{10}(M_\mathrm{CO-core} / M_\odot)$'],
    'he_core_radius':
        [r'$R_\mathrm{He-core} \, R_\odot$',
         r'$\log_{10}(R_\mathrm{He-core} / R_\odot)$'],
    'c_core_radius':
        [r'$R_\mathrm{C-core} \, [M_\odot]$',
         r'$\log_{10}(R_\mathrm{C-core} / R_\odot)$'],
    'o_core_radius':
        [r'$R_\mathrm{O-core} \, [M_\odot]$',
         r'$\log_{10}(R_\mathrm{O-core} / R_\odot)$'],
    'co_core_radius':
        [r'$R_\mathrm{CO-core} \, [M_\odot]$',
         r'$\log_{10}(R_\mathrm{CO-core} / R_\odot)$'],
    'center_h1':
        [r'${}^1H_\mathrm{center}$', r'$\log_{10}({}^1H_\mathrm{center})$'],
    'center_he4':
        [r'${}^4He_\mathrm{center}$', r'$\log_{10}({}^4He_\mathrm{center})$'],
    'center_c12':
        [r'${}^{12}C_\mathrm{center}$',
         r'$\log_{10}({}^{12}C_\mathrm{center})$'],
    'center_n14':
        [r'${}^{14}N_\mathrm{center}$',
         r'$\log_{10}({}^{14}N_\mathrm{center})$'],
    'center_o16':
        [r'${}^{16}O_\mathrm{center}$',
         r'$\log_{10}({}^{16}O_\mathrm{center})$'],
    'surface_h1':
        [r'${}^1H_\mathrm{surface}$', r'$\log_{10}({}^1H_\mathrm{surface})$'],
    'surface_he4':
        [r'${}^4He_\mathrm{surface}$',
         r'$\log_{10}({}^4He_\mathrm{surface})$'],
    'surface_c12':
        [r'${}^{12}C_\mathrm{surface}$',
         r'$\log_{10}({}^{12}C_\mathrm{surface})$'],
    'surface_n14':
        [r'${}^{14}N_\mathrm{surface}$',
         r'$\log_{10}({}^{14}N_\mathrm{surface})$'],
    'surface_o16':
        [r'${}^{16}O_\mathrm{surface}$',
         r'$\log_{10}({}^{16}O_\mathrm{surface})$'],
    'c12_c12':
        [r'$c12_c12\,[L_\odot]$', r'$\log_{10}(c12_c12/L_\odot)$'],
    'log_LH':
        [r'$\log_{10}(L_\mathrm{H}/L_\odot)$', 'log_log_LH'],
    'log_LHe':
        [r'$\log_{10}(L_\mathrm{He}/L_\odot)$',
         'log_log_LHe'],
    'log_LZ':
        [r'$\log_{10}(L_\mathrm{Z}/L_\odot)$',
         'log_log_LZ'],
    'log_Lnuc':
        [r'$\log_{10}(L_\mathrm{nuc}/L_\odot)$',
         'log_log_Lnuc'],
    'log_Teff':
        [r'$\log_{10}(T_\mathrm{eff}/\mathrm{K})$', 'log_log_Teff'],
    'log_L':
        [r'$\log_{10}(L_\mathrm{surf}/L_\odot)$',
         'log_log_Lsurf'],
    'log_R':
        [r'$\log_{10}(R/R_\odot)$', r'$\log_{10}(\log_{10}(R/R_\odot))$'],
    'center_gamma':
        ['center_gamma', 'log_center_gamma'],
    'avg_c_in_c_core':
        ['avg_c_in_c_core', 'log_avg_c_in_c_core'],
    'surf_avg_omega':
        [r'$\omega_\mathrm{s}\,[\mathrm{yr}^{-1}]$',
         r'$\log_{10}(\omega_\mathrm{s}/\mathrm{yr}^{-1})'],
    'surf_avg_omega_div_omega_crit':
        [r'$(\omega_\mathrm{s}/\omega_\mathrm{s,crit})$',
         r'$\log_{10}(\omega_\mathrm{s}/\omega_\mathrm{s,crit})$'],
    'total_moment_of_inertia':
        [r'$I_\mathrm{tot}\,[g\,\mathrm{cm}^{2}]$',
         r'$\log_{10}(I_\mathrm{tot}/(g\,\mathrm{cm}^{2}))$'],
    'spin_parameter':
        [r'$a_\star$', r'$\log_{10}(a_\star)$'],
    'log_total_angular_momentum':
        [r'$\log_{10}(J_\mathrm{tot}/(g\,\mathrm{cm}^{2}\mathrm{s}^{-1}]))$',
         'log_log_total_angular_momentum'],
    'conv_env_top_mass':
        [r'$M_\mathrm{top-conv-env}\,[M_\odot]$',
         r'$\log_{10}(M_\mathrm{top-conv-env}/M_\odot)$'],
    'conv_env_bot_mass':
        [r'$M_\mathrm{bot-conv-env}\,[M_\odot]$',
         r'$\log_{10}(M_\mathrm{bot-conv-env}/M_\odot)$'],
    'conv_env_top_radius':
        [r'$R_\mathrm{top-conv-env}\,[R_\odot]$',
         r'$\log_{10}(R_\mathrm{top-conv-env}/M_\odot)$'],
    'conv_env_bot_radius':
        [r'$R_\mathrm{bot-conv-env}\,[R_\odot]$',
         r'$\log_{10}(R_\mathrm{bot-conv-env}/M_\odot)$'],
    'conv_env_turnover_time_g':
        [r'$t^\mathrm{turnover-g}_{conv-env}\,[\mathrm{yr}]$',
         r'$\log_{10}(t^\mathrm{turnover-g}_{conv-env}/\mathrm{yr})$'],
    'conv_env_turnover_time_l_b':
        [r'$t^\mathrm{turnover-l-b}_{conv-env}\,[\mathrm{yr}]$',
         r'$\log_{10}(t^\mathrm{turnover-l-b}_{conv-env}/\mathrm{yr})$'],
    'conv_env_turnover_time_l_t':
        [r'$t^\mathrm{turnover-l-t}_{conv-env}\,[\mathrm{yr}]$',
         r'$\log_{10}(t^\mathrm{turnover-l-t}_{conv-env}/\mathrm{yr})$'],
    'envelope_binding_energy':
        [r'$E_\mathrm{env}\,[\mathrm{erg}]$',
         r'$\log_{10}(E_\mathrm{env}/(\mathrm{erg})$'],
    'mass_conv_reg_fortides':
        [r'$M_\mathrm{conv-reg}\,[M_\odot]$',
         r'$\log_{10}(M_\mathrm{conv-reg}/M_\odot)$'],
    'thickness_conv_reg_fortides':
        [r'$dR_\mathrm{conv-reg}\,[R_\odot]$',
         r'$\log_{10}(dR_\mathrm{conv-reg}/R_\odot)$'],
    'radius_conv_reg_fortides':
        [r'$R_\mathrm{conv-reg}\,[R_\odot]$',
         r'$\log_{10}(R_\mathrm{conv-reg}/R_\odot)$'],
    'lambda_CE_1cent':
        [r'$\lambda_\mathrm{CE,1\%}$',
         r'$\log_{10}(\lambda_\mathrm{CE,1\%})$'],
    'lambda_CE_10cent':
        [r'$\lambda_\mathrm{CE,10\%}$',
         r'$\log_{10}(\lambda_\mathrm{CE,10\%})$'],
    'lambda_CE_30cent':
        [r'$\lambda_\mathrm{CE,30\%}$',
         r'$\log_{10}(\lambda_\mathrm{CE,30\%})$'],
    'lambda_CE_pure_He_star_10cent':
        [r'$\lambda_\mathrm{CE-He-star,10\%}$',
         r'$\log_{10}(\lambda_\mathrm{CE-He-star,30\%})$'],

    # binary_history
    'model_number':
        ['model_number', 'log_model_number'],
    'age':
        [r'$t \, [\mathrm{yr}]$', r'$\log_{10}(t / \mathrm{yr})$'],
    'star_1_mass':
        [r'$M_\mathrm{1} \, [M_\odot]$',
         r'$\log_{10}(M_\mathrm{1} / M_\odot)$'],
    'star_2_mass':
        [r'$M_\mathrm{2} \, [M_\odot]$',
         r'$\log_{10}(M_\mathrm{2} / M_\odot)$'],
    'period_days':
        [r'$P_\mathrm{orb} \, [\mathrm{days}]$',
         r'$\log_{10}(P_\mathrm{orb} / \mathrm{days})$'],
    'binary_separation':
        [r'$A \, [R_\odot]$',
         r'$\log_{10}(A / R_\odot)$'],
    'rl_relative_overflow_1':
        ['relative RL overflow 1', 'log10 relative RL overflow 1'],
    'rl_relative_overflow_2':
        ['relative RL overflow 2', 'log10 relative RL overflow 2'],
    'lg_mtransfer_rate':
        [r'$\log_{10}(\dot{M}/(M_\odot\,\mathrm{yr}^{-1}))$',
         'log_log_mtransfer_rate'],
    'lg_system_mdot_1':
        [r'$\log_{10}(\dot{M}_1/(M_\odot\,\mathrm{yr}^{-1}))$',
         'log_log_system_mdot_1'],
    'lg_system_mdot_2':
        [r'$\log_{10}(\dot{M}_2/(M_\odot\,\mathrm{yr}^{-1}))$',
         'log_log_system_mdot_2'],
    'lg_wind_mdot_1':
        [r'$\log_{10}(\dot{M}_\mathrm{1,wind}/(M_\odot\,\mathrm{yr}^{-1}))$',
         'log_log_wind_mdot_1'],
    'lg_wind_mdot_2':
        [r'$\log_{10}(\dot{M}_\mathrm{2,wind}/(M_\odot\,\mathrm{yr}^{-1}))$',
         'log_log_wind_mdot_2'],
    'lg_mstar_dot_1':
        [r'$\log_{10}(\dot{M}_\mathrm{\star,1}/(M_\odot\,\mathrm{yr}^{-1}))$',
         'log_log_star_dot_1'],
    'lg_mstar_dot_2':
        [r'$\log_{10}(\dot{M}_\mathrm{\star,2}/(M_\odot\,\mathrm{yr}^{-1}))$',
         'log_log_star_dot_1'],
    'xfer_fraction':
        ['xfer_fraction', 'log_xfer_fraction'],
    'trap_radius':
        [r'$R_\mathrm{trap}\,[R_\odot]$',
         r'$\log_{10}(R_\mathrm{trap}/R_\odot$)'],
    'acc_radius':
        [r'$R_\mathrm{acc}\,[\mathrm{cm}]$',
         r'$\log_{10}(R_\mathrm{acc}/\mathrm{cm})$'],
    't_sync_rad_1':
        [r'$t^1_\mathrm{sync-rad}\,[\mathrm{s}]$',
         r'$\log_{10}(t^1_\mathrm{sync-rad}/\mathrm{s})$'],
    't_sync_conv_1':
        [r'$t^1_\mathrm{conv}\,[\mathrm{s}]$',
         r'$\log_{10}(t^1_\mathrm{conv}/\mathrm{s})$'],
    't_sync_rad_2':
        [r'$t^2_\mathrm{sync-rad}\,[\mathrm{s}]$',
         r'$\log_{10}(t^2_\mathrm{sync-rad}/\mathrm{s})$'],
    't_sync_conv_2':
        [r'$t^2_\mathrm{conv}\,[\mathrm{s}]$',
         r'$\log_{10}(t^2_\mathrm{conv}/\mathrm{s})$'],
}
