"""Module for combining termination flags for visualization purposes."""


__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]


import numpy as np

TF1_POOL_STABLE = ['Primary has depleted central carbon',
                   'Secondary has depleted central carbon',
                   'Primary enters pair-instability regime',
                   'Secondary enters pair-instability regime',
                   'Primary enters pulsational pair-instability regime',
                   'Secondary enters pulsational pair-instability regime',
                   'offcenter neon ignition for primary',
                   'offcenter neon ignition for secondary',
                   'envelope_mass_limit',
                   'gamma_center_limit',
                   'Reached TPAGB',
                   'max_age']

TF1_POOL_UNSTABLE = [
    'Unstable',
    'overflow from L2 (D_L2) distance for q(=Macc/Mdon)>1, donor is star 1',
    'overflow from L2 (R_L2) surface for q(=Macc/Mdon)<1, donor is star 1',
    'overflow from L2 (R_L2) surface for q(=Macc/Mdon)>1, donor is star 1',
    'overflow from L2 (D_L2) distance for q(=Macc/Mdon)<1, donor is star 1',
    'overflow from L2 (R_L2) surface for q(=Macc/Mdon)<1, donor is star 2',
    'overflow from L2 (R_L2) surface for q(=Macc/Mdon)>1, donor is star 2',
    'reached maximum mass transfer rate: 10.0d0',
    'Reached maximum mass transfer rate: 1d-1',
    'Reached the critical mt rate',
    'Reached maximum mass transfer rate: Exceeded photon trapping radius',
    'Both stars fill their Roche Lobe and at least one of them is off MS',
    'Terminate due to L2 overflow during case A']

TF1_POOL_INITIAL_RLO = ['Initial RLOF',
                        'forced_initial_RLO',
                        'overflow from L1 at ZAMS',
                        'Terminate because of overflowing initial model',
                        'overflow from L2 point for q>1 at ZAMS',
                        'overflow from L2 surface for q<1 at ZAMS',
                        'overflow from L2 surface for q>1 at ZAMS',
                        'ignored_no_BH',
                        'ignored_scrubbed']

TF1_POOL_ERROR = ['logQ_limit',
                  'min timestep',
                  'min_timestep_limit',
                  'cluster timelimit',
                  'reach cluster timelimit',
                  'no termination code',
                  'fe_core_infall_limit']

TF2_POOL_NO_RLO = [
    'no_RLOF',
    'undetermined_flag_mass_transfer_from_evolutionary_states_from_star1'
]

TF2_POOL_INITIAL_RLO = [
    'initial_RLOF',
    'ignored_scrubbed',
    'ignored_no_BH'
]

TF2_UNSTABLE_L2 = [
    "L2_RLOF"
]


def combine_TF12(IC, TF2, verbose=False):
    """Get the combination of interpolation classion and termination flag 2."""
    N = len(IC)
    TF12 = np.array(['unknown']*N, dtype='U25')

    # INITIAL RLO
    TF2_pool_SCONT = ['contact_during_MS']

    TF2_pool_SA = ['caseA_from_star1',
                   'caseA_from_star2']

    TF2_pool_SAB = ['caseA/B_from_star1',
                    'caseA/B_from_star2',
                    'caseA/B/C_from_star1',
                    'caseA/B/C/BB_from_star1',
                    'caseA/B/A_from_star1']

    TF2_pool_SABB = ['caseA/BB_from_star1',
                     'caseA/BB_from_star2']

    TF2_pool_SB = ['caseB_from_star1',
                   'caseB_from_star2']

    TF2_pool_SC = ['caseC_from_star1',
                   'caseC_from_star2']

    TF2_pool_SBC = ['caseB/C_from_star1',
                    'caseB/C_from_star2']

    TF2_pool_SBB = ['caseB/C/BB_from_star1',
                    'caseB/C/BB_from_star2',
                    'caseB/BB_from_star1',
                    'caseB/BB_from_star2',
                    'caseA/B/BB_from_star1',
                    'caseA/B/BB_from_star2',
                    'caseBA/BB_from_star1',
                    'caseB/A_from_star1',
                    'caseC/BB_from_star1',
                    'caseBB_from_star1',
                    'caseBB_from_star2']

    TF2_pool_SBA = ['caseBA_from_star1',
                    'caseBA_from_star2']

    TF2_unknown = ['?_from_star1', '?_from_star2']

    TF2 = [(tf[1:] if tf.startswith("?case") else tf) for tf in TF2]

    for i in range(N):
        if IC[i] == 'initial_MT':
            TF12[i] = 'Initial RLOF'
        elif IC[i] == 'not_converged':
            TF12[i] = 'Not converged'
        elif IC[i] == 'no_MT':
            TF12[i] = 'no_RLOF'
        elif IC[i] == 'stable_MT':
            if TF2[i] in TF2_pool_SCONT:
                TF12[i] = 'Stable contact'
            elif TF2[i] in TF2_pool_SA:
                TF12[i] = 'Stable case A'
            elif TF2[i] in TF2_pool_SAB:
                TF12[i] = 'Stable case AB'
            elif TF2[i] in TF2_pool_SABB:
                TF12[i] = 'Stable case ABB'
            elif TF2[i] in TF2_pool_SB:
                TF12[i] = 'Stable case B'
            elif TF2[i] in TF2_pool_SC:
                TF12[i] = 'Stable case C'
            elif TF2[i] in TF2_pool_SBC:
                TF12[i] = 'Stable case BC'
            elif TF2[i] in TF2_pool_SBB:
                TF12[i] = 'Stable case BB'
            elif TF2[i] in TF2_pool_SBA:
                TF12[i] = 'Stable case BA'
        elif IC[i] == 'unstable_MT':
            if TF2[i] in TF2_pool_SCONT:
                TF12[i] = 'Unstable contact'
            elif TF2[i] in TF2_pool_SA:
                TF12[i] = 'Unstable case A'
            elif TF2[i] in TF2_pool_SAB:
                TF12[i] = 'Unstable case AB'
            elif TF2[i] in TF2_pool_SABB:
                TF12[i] = 'Unstable case ABB'
            elif TF2[i] in TF2_pool_SB:
                TF12[i] = 'Unstable case B'
            elif TF2[i] in TF2_pool_SC:
                TF12[i] = 'Unstable case C'
            elif TF2[i] in TF2_pool_SBB:
                TF12[i] = 'Unstable case BB'
            elif TF2[i] in TF2_pool_SBC:
                TF12[i] = 'Unstable case BC'
            elif TF2[i] in TF2_UNSTABLE_L2:
                TF12[i] = "Unstable L2 RLOF"
            elif TF2[i] in TF2_pool_SBA:
                TF12[i] = 'Unstable case BA'

        # catch if something is missing from the logic
        if TF12[i] == 'unknown' or TF2[i] in TF2_unknown:
            if verbose:
                print(i, IC[i], TF2[i])

    return TF12
