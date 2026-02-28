"""Module for combining termination flags for visualization purposes."""


__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]


import re

import numpy as np

# TODO: rename the varaible, these are termination flags that indicates
# the star has reached the end of the evolution
TF1_POOL_STABLE = [
    'Primary has depleted central carbon',
    'Secondary has depleted central carbon',
    'Primary got stopped before central carbon depletion',
    'Secondary got stopped before central carbon depletion',
    'Primary enters pair-instability regime',
    'Secondary enters pair-instability regime',
    'Primary enters pulsational pair-instability regime',
    'Secondary enters pulsational pair-instability regime',
    # 'offcenter neon ignition for primary',
    # 'offcenter neon ignition for secondary',
    'envelope_mass_limit',
    'gamma_center_limit',
    # 'Reached TPAGB',
    'max_age'
    ]

TF1_POOL_UNSTABLE = [
    'overflow from L2 (D_L2) distance for q(=Macc/Mdon)>1, donor is star 1',
    'overflow from L2 (R_L2) surface for q(=Macc/Mdon)<1, donor is star 1',
    'overflow from L2 (R_L2) surface for q(=Macc/Mdon)>1, donor is star 1',
    'overflow from L2 (D_L2) distance for q(=Macc/Mdon)<1, donor is star 1',
    'overflow from L2 (R_L2) surface for q(=Macc/Mdon)<1, donor is star 2',
    'overflow from L2 (R_L2) surface for q(=Macc/Mdon)>1, donor is star 2',
    'Reached maximum mass transfer rate: 1d-1',
    'Reached maximum mass transfer rate: Exceeded photon trapping radius',
    'Both stars fill their Roche Lobe and at least one of them is off MS',
    'Terminate due to L2 overflow during case A',
    'Both stars fill their Roche Lobe and t_kh > t_acc',
    'overflow from L2, t_kh > t_acc and w > w_crit_lim, donor is star 1',
    'overflow from L2, t_kh > t_acc and w > w_crit_lim, donor is star 2'
    ]

TF1_POOL_INITIAL_RLO = [
    'forced_initial_RLO',
    'Terminate because of overflowing initial model',
    ]

TF1_POOL_ERROR = [
    'logQ_min_limit',
    'min_timestep_limit',
    'reach cluster timelimit',
    ]

TF2_POOL_NO_RLO = [
    'no_RLOF',
    ]

TF2_POOL_INITIAL_RLO = [
    'initial_RLOF',
    ]

TF2_POOL_UNSTABLE = [
    "L2_RLOF"
    ]

TF2_POOL_UNKNOWN = [
    'None'
    ]

# Legacy pool kept for backward compatibility with old grid data
TF2_POOL_CONT = [
    'contact_during_MS'
    ]

TF2_POOL_A = [
    'case_A1',
    'case_A2'
    ]

TF2_POOL_B = [
    'case_B1',
    'case_B2'
    ]

TF2_POOL_C = [
    'case_C1',
    'case_C2',
    ]

TF2_POOL_BA = [
    'case_BA1',
    'case_BA2'
    ]

TF2_POOL_BB = [
#    'case_A1/B1/BB1',
#    'case_A2/B2/BB2',
#    'case_B1/BB1',
#    'case_B2/BB2',
#    'case_B1/C1/BB1',
#    'case_B2/C2/BB2',
#    'case_BA1/BB1',
#    'case_BA2/BB2',
#    'case_C1/BB1',
#    'case_C2/BB2',
    'case_BB1',
    'case_BB2'
    ]

#TF2_POOL_AB = [
#    'case_A1/B1',
#    'case_A2/B2',
#    'case_A1/B1/C1',
#    'case_A2/B2/C2',
#    'case_A1/B1/C1/BB1',
#    'case_A2/B2/C2/BB2',
#    'case_A1/B1/A1',
#    'case_A2/B2/A2'
#    ]

#TF2_POOL_AC = [
#    'case_A1/C1',
#    'case_A2/C2'
#    ]

#TF2_POOL_ABB = [
#    'case_A1/BB1',
#    'case_A2/BB2'
#    ]

#TF2_POOL_BC = [
#    'case_B1/C1',
#    'case_B2/C2'
#    ]


def _strip_contact_markers(tf2):
    """Strip contact 'c' markers from a TF2 string.

    E.g., 'case_Ac1' -> 'case_A1', 'case_Ac1/Bc2' -> 'case_A1/B2'
    """
    return re.sub(r'c([12])', r'\1', tf2)


def _has_contact(tf2):
    """Check if any MT case in the TF2 string has a contact marker."""
    return bool(re.search(r'c[12]', tf2))


def combine_TF12(IC, TF2, verbose=False):
    """Get the combination of interpolation classion and termination flag 2."""
    N = len(IC)
    TF12 = np.array(['unknown']*N, dtype='U30')

    # mask unknown cases during history TF2
    TF2 = [(tf[1:] if tf.startswith("?case") else tf) for tf in TF2]

    for i in range(N):
        if IC[i] == 'initial_MT':
            TF12[i] = 'Initial RLOF'
        elif IC[i] == 'not_converged':
            TF12[i] = 'Not converged'
        elif IC[i] == 'no_MT':
            TF12[i] = 'no_RLOF'
        elif IC[i] == 'stable_reverse_MT':
            TF12[i] = 'Reverse stable MT'
        elif IC[i] == 'stable_MT':
            # Check for contact markers in the TF2 string
            is_contact = _has_contact(TF2[i])
            contact_prefix = 'contact ' if is_contact else ''
            # Strip contact markers for base classification
            tf2_base = _strip_contact_markers(TF2[i])
            if '1' in tf2_base and '2' in tf2_base:
                TF12[i] = 'Reverse stable MT'
            elif tf2_base in TF2_POOL_CONT:
                # Legacy: old grid data with 'contact_during_MS'
                TF12[i] = 'Stable contact'
            elif tf2_base in TF2_POOL_A:
                TF12[i] = 'Stable {}case A'.format(contact_prefix)
            elif tf2_base in TF2_POOL_B:
                TF12[i] = 'Stable {}case B'.format(contact_prefix)
            elif tf2_base in TF2_POOL_C:
                TF12[i] = 'Stable {}case C'.format(contact_prefix)
            elif tf2_base in TF2_POOL_BA:
                TF12[i] = 'Stable {}case BA'.format(contact_prefix)
            elif tf2_base in TF2_POOL_BB:
                TF12[i] = 'Stable {}case BB'.format(contact_prefix)
#            elif tf2_base in TF2_POOL_AB:
#                TF12[i] = 'Stable {}case AB'.format(contact_prefix)
#            elif tf2_base in TF2_POOL_AC:
#                TF12[i] = 'Stable {}case AC'.format(contact_prefix)
#            elif tf2_base in TF2_POOL_ABB:
#                TF12[i] = 'Stable {}case ABB'.format(contact_prefix)
#            elif tf2_base in TF2_POOL_BC:
#                TF12[i] = 'Stable {}case BC'.format(contact_prefix)
            elif 'case_nonburning' in tf2_base:
                TF12[i] = 'Stable {}case n'.format(contact_prefix)
            elif '/' in tf2_base:
                # multiple cases, split them up
                MTcases = tf2_base.split('/')
                # record first case
                TF12[i] = 'Stable {}'\
                           .format(contact_prefix)+MTcases[0]\
                           .replace('_',' ')[:-1]
                # record last case
                if 'nonburning' in MTcases[-1]:
                    TF12[i] += 'n'
                else:
                    TF12[i] += MTcases[-1][:-1]
        elif IC[i] == 'unstable_MT':
            # Check for contact markers in the TF2 string
            is_contact = _has_contact(TF2[i])
            contact_prefix = 'contact ' if is_contact else ''
            # Strip contact markers for base classification
            tf2_base = _strip_contact_markers(TF2[i])
            if '1' in tf2_base and '2' in tf2_base:
                TF12[i] = 'Reverse unstable MT'
            elif tf2_base in TF2_POOL_CONT:
                # Legacy: old grid data with 'contact_during_MS'
                TF12[i] = 'Unstable contact'
            elif tf2_base in TF2_POOL_A:
                TF12[i] = 'Unstable {}case A'.format(contact_prefix)
            elif tf2_base in TF2_POOL_B:
                TF12[i] = 'Unstable {}case B'.format(contact_prefix)
            elif tf2_base in TF2_POOL_C:
                TF12[i] = 'Unstable {}case C'.format(contact_prefix)
            elif tf2_base in TF2_POOL_BA:
                TF12[i] = 'Unstable {}case BA'.format(contact_prefix)
            elif tf2_base in TF2_POOL_BB:
                TF12[i] = 'Unstable {}case BB'.format(contact_prefix)
#            elif tf2_base in TF2_POOL_AB:
#                TF12[i] = 'Unstable {}case AB'.format(contact_prefix)
#            elif tf2_base in TF2_POOL_AC:
#                TF12[i] = 'Unstable {}case AC'.format(contact_prefix)
#            elif tf2_base in TF2_POOL_ABB:
#                TF12[i] = 'Unstable {}case ABB'.format(contact_prefix)
#            elif tf2_base in TF2_POOL_BC:
#                TF12[i] = 'Unstable {}case BC'.format(contact_prefix)
            elif tf2_base in TF2_POOL_UNSTABLE:
                TF12[i] = "Unstable L2 RLOF"
            elif 'case_nonburning' in tf2_base:
                TF12[i] = 'Unstable {}case n'.format(contact_prefix)
            elif '/' in tf2_base:
                # multiple cases, split them up
                MTcases = tf2_base.split('/')
                # record first case
                TF12[i] = 'Unstable {}'\
                           .format(contact_prefix)+MTcases[0]\
                           .replace('_',' ')[:-1]
                # record last case
                if 'nonburning' in MTcases[-1]:
                    TF12[i] += 'n'
                else:
                    TF12[i] += MTcases[-1][:-1]

        # catch if something is missing from the logic
        if TF12[i] == 'unknown':
            if verbose:
                print(i, IC[i], TF2[i])

    return TF12
