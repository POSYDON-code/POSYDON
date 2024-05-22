"""Module for combining termination flags for visualization purposes."""


__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]


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

def combine_TF12(IC, TF2, verbose=False):
    """Get the combination of interpolation classion and termination flag 2."""
    N = len(IC)
    TF12 = np.array(['unknown']*N, dtype='U25')

    # mask unknown cases during history TF2
    TF2 = [(tf[1:] if tf.startswith("?case") else tf) for tf in TF2]

    for i in range(N):
        if IC[i] == 'initial_MT':
            TF12[i] = 'Initial RLOF'
        elif IC[i] == 'not_converged':
            TF12[i] = 'Not converged'
        elif IC[i] == 'no_MT':
            TF12[i] = 'no_RLOF'
        elif IC[i] == 'stable_MT':
            if '1' in TF2[i] and '2' in TF2[i]:
                TF12[i] = 'Reverse stable MT'
            elif TF2[i] in TF2_POOL_CONT:
                TF12[i] = 'Stable contact'
            elif TF2[i] in TF2_POOL_A:
                TF12[i] = 'Stable case A'
            elif TF2[i] in TF2_POOL_B:
                TF12[i] = 'Stable case B'
            elif TF2[i] in TF2_POOL_C:
                TF12[i] = 'Stable case C'
            elif TF2[i] in TF2_POOL_BA:
                TF12[i] = 'Stable case BA'
            elif TF2[i] in TF2_POOL_BB:
                TF12[i] = 'Stable case BB'
#            elif TF2[i] in TF2_POOL_AB:
#                TF12[i] = 'Stable case AB'
#            elif TF2[i] in TF2_POOL_AC:
#                TF12[i] = 'Stable case AC'
#            elif TF2[i] in TF2_POOL_ABB:
#                TF12[i] = 'Stable case ABB'
#            elif TF2[i] in TF2_POOL_BC:
#                TF12[i] = 'Stable case BC'
            elif 'case_nonburning' in TF2[i]:
                TF12[i] = 'Stable case n'
            elif '/' in TF2[i]:
                # multiple cases, split them up
                MTcases = TF2[i].split('/')
                # record first case
                TF12[i] = 'Stable '+MTcases[0].replace('_',' ')[:-1]
                # record last case
                if 'nonburning' in MTcases[-1]:
                    TF12[i] += 'n'
                else:
                    TF12[i] += MTcases[-1][:-1]
        elif IC[i] == 'unstable_MT':
            if '1' in TF2[i] and '2' in TF2[i]:
                TF12[i] = 'Reverse unstable MT'
            elif TF2[i] in TF2_POOL_CONT:
                TF12[i] = 'Unstable contact'
            elif TF2[i] in TF2_POOL_A:
                TF12[i] = 'Unstable case A'
            elif TF2[i] in TF2_POOL_B:
                TF12[i] = 'Unstable case B'
            elif TF2[i] in TF2_POOL_C:
                TF12[i] = 'Unstable case C'
            elif TF2[i] in TF2_POOL_BA:
                TF12[i] = 'Unstable case BA'
            elif TF2[i] in TF2_POOL_BB:
                TF12[i] = 'Unstable case BB'
#            elif TF2[i] in TF2_POOL_AB:
#                TF12[i] = 'Unstable case AB'
#            elif TF2[i] in TF2_POOL_AC:
#                TF12[i] = 'Unstable case AC'
#            elif TF2[i] in TF2_POOL_ABB:
#                TF12[i] = 'Unstable case ABB'
#            elif TF2[i] in TF2_POOL_BC:
#                TF12[i] = 'Unstable case BC'
            elif TF2[i] in TF2_POOL_UNSTABLE:
                TF12[i] = "Unstable L2 RLOF"
            elif 'case_nonburning' in TF2[i]:
                TF12[i] = 'Unstable case n'
            elif '/' in TF2[i]:
                # multiple cases, split them up
                MTcases = TF2[i].split('/')
                # record first case
                TF12[i] = 'Unstable '+MTcases[0].replace('_',' ')[:-1]
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
