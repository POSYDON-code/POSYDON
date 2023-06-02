"""
Module setting the default evolution flow and allowing alterations.

This file contains the offical POSYDON flow chart binary evolution. To support
future development and deal with complexity it is build dynamically.
"""


__authors__ = [
    "Devina Misra <devina.misra@unige.ch>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Tassos Fragos <Anastasios.Fragkos@unige.ch>",
    "Zepei Xing <Zepei.Xing@unige.ch>",
]


STAR_STATES_ALL = [
    'WD',
    'NS',
    'BH',
    'H-rich_Core_H_burning',
    'H-rich_Core_He_burning',
    'H-rich_Shell_H_burning',
    'H-rich_Central_He_depleted',
    'H-rich_Shell_He_burning',
    'H-rich_Core_C_burning',
    'H-rich_Central_C_depletion',
    'H-rich_non_burning',
    'stripped_He_Core_He_burning',
    'stripped_He_Central_He_depleted',
    'stripped_He_Central_C_depletion',
    'stripped_He_non_burning',
    'massless_remnant'
]

STAR_STATES_CO = ['BH', 'NS', 'WD','massless_remnant']

STAR_STATES_NOT_CO = STAR_STATES_ALL.copy()
[STAR_STATES_NOT_CO.remove(x) for x in STAR_STATES_CO]

STAR_STATES_H_RICH = STAR_STATES_NOT_CO.copy()
[STAR_STATES_H_RICH.remove(x) for x in ['stripped_He_Core_He_burning',
                                        'stripped_He_Central_He_depleted',
                                        'stripped_He_Central_C_depletion',
                                        'stripped_He_non_burning']]

STAR_STATES_HE_RICH = STAR_STATES_NOT_CO.copy()
[STAR_STATES_HE_RICH.remove(x) for x in ['H-rich_Core_H_burning',
                                         'H-rich_Core_He_burning',
                                         'H-rich_Shell_H_burning',
                                         'H-rich_Central_He_depleted',
                                         'H-rich_Shell_He_burning',
                                         'H-rich_Core_C_burning',
                                         'H-rich_Central_C_depletion']]

STAR_STATES_C_DEPLETION = [st for st in STAR_STATES_ALL if "C_depletion" in st]

# these states can be evolved through MESA grids
STAR_STATES_H_RICH_EVOLVABLE = list(set(STAR_STATES_H_RICH)
                                    - set(STAR_STATES_C_DEPLETION))

STAR_STATES_HE_RICH_EVOLVABLE = list(set(STAR_STATES_HE_RICH)
                                     - set(STAR_STATES_C_DEPLETION))
# CE ejcetion happens istantanously, the star does not readjust before
# we infer the state, if core_definition_H_fraction=0.1 then surface_h1=0.1
# and the state is H-rich_non_burning which we stil want to evolve thorugh
# the step_CO_HeMS
STAR_STATES_HE_RICH_EVOLVABLE.append('H-rich_non_burning')

BINARY_STATES_ALL = [
    'initially_single_star',
    'detached',
    'RLO1',
    'RLO2',
    'contact',
    'disrupted',
    'merged',
    'initial_RLOF'
]

BINARY_EVENTS_ALL = [
    None,
    'CC1',
    'CC2',
    'ZAMS',
    'oRLO1',
    'oRLO2',
    'oCE1',
    'oCE2',
    'oDoubleCE1',
    'oDoubleCE2',
    'CO_contact',
    'redirect',
    'MaxTime_exceeded',
    'maxtime',
    'oMerging1',
    'oMerging2'
]


# dynamically construct the flow chart
POSYDON_FLOW_CHART = {}

# mesa grid ZAMS
STAR_STATES_ZAMS = ['H-rich_Core_H_burning']
BINARY_STATES_ZAMS = ['detached']

for b in BINARY_STATES_ZAMS:
    for s1 in STAR_STATES_ZAMS:
        for s2 in STAR_STATES_ZAMS:
            POSYDON_FLOW_CHART[(s1, s2, b, 'ZAMS')] = 'step_HMS_HMS'

# ZAMS very wide binary that falls outside the grid
# and has been returned by step_HMS_HMS
for b in BINARY_STATES_ZAMS:
    for s1 in STAR_STATES_ZAMS:
        for s2 in STAR_STATES_ZAMS:
            POSYDON_FLOW_CHART[(s1, s2, b, 'redirect')] = 'step_detached'


# stripped_He star on a detached binary another H- or stripped_He star
# This will be the outcome of a CE.
for s1 in STAR_STATES_NOT_CO:
    for s2 in STAR_STATES_NOT_CO:
        POSYDON_FLOW_CHART[(s1, s2, 'detached', None)] = 'step_detached'
        POSYDON_FLOW_CHART[(s2, s1, 'detached', None)] = 'step_detached'


# H-rich star on a detached binary with a compact object
for s1 in STAR_STATES_H_RICH:
    for s2 in STAR_STATES_CO:
        POSYDON_FLOW_CHART[(s1, s2, 'detached', None)] = 'step_detached'
        POSYDON_FLOW_CHART[(s2, s1, 'detached', None)] = 'step_detached'


# H-rich star roche-lobe overflow onto a compact object
for s1 in STAR_STATES_H_RICH_EVOLVABLE:
    for s2 in STAR_STATES_CO:
        POSYDON_FLOW_CHART[(s1, s2, 'RLO1', 'oRLO1')] = 'step_CO_HMS_RLO'
        POSYDON_FLOW_CHART[(s2, s1, 'RLO2', 'oRLO2')] = 'step_CO_HMS_RLO'

# H-rich star on a detached binary with a compact object
# that fall outside the grid and has been returned by step_CO_HMS_RLO
for s1 in STAR_STATES_H_RICH:
    for s2 in STAR_STATES_CO:
        POSYDON_FLOW_CHART[(s1, s2, 'detached', "redirect")] = 'step_detached'
        POSYDON_FLOW_CHART[(s2, s1, 'detached', "redirect")] = 'step_detached'

# stripped_He star on a detached binary with a compact object
for s1 in STAR_STATES_HE_RICH_EVOLVABLE:
    for s2 in STAR_STATES_CO:
        POSYDON_FLOW_CHART[(s1, s2, 'detached', None)] = 'step_CO_HeMS'
        POSYDON_FLOW_CHART[(s2, s1, 'detached', None)] = 'step_CO_HeMS'


# stripped_He star on a detached binary with a compact object
# that fall outside the grid and has been returned by step_CO_HeMS
for s1 in STAR_STATES_HE_RICH:
    for s2 in STAR_STATES_CO:
        POSYDON_FLOW_CHART[(s1, s2, 'detached', "redirect")] = 'step_detached'
        POSYDON_FLOW_CHART[(s2, s1, 'detached', "redirect")] = 'step_detached'


# Binaries that go to common envelope

for s1 in STAR_STATES_NOT_CO:
    for s2 in STAR_STATES_ALL:
        POSYDON_FLOW_CHART[(s1, s2, 'RLO1', 'oCE1')] = 'step_CE'
        POSYDON_FLOW_CHART[(s1, s2, 'RLO1', 'oDoubleCE1')] = 'step_CE'
        POSYDON_FLOW_CHART[(s2, s1, 'RLO2', 'oCE2')] = 'step_CE'
        POSYDON_FLOW_CHART[(s2, s1, 'RLO2', 'oDoubleCE2')] = 'step_CE'
        POSYDON_FLOW_CHART[(s2, s1, 'contact', 'oCE1')] = 'step_CE'
        POSYDON_FLOW_CHART[(s2, s1, 'contact', 'oCE2')] = 'step_CE'
        POSYDON_FLOW_CHART[(s2, s1, 'contact', 'oDoubleCE1')] = 'step_CE'
        POSYDON_FLOW_CHART[(s2, s1, 'contact', 'oDoubleCE2')] = 'step_CE'


# core collapse
STAR_STATES_CC = [
    'H-rich_Central_C_depletion',
    'H-rich_Central_He_depleted',
    'stripped_He_Central_He_depleted',
    'stripped_He_Central_C_depletion',
    # catch runs with gamma center limit which map to WD
    'stripped_He_non_burning',
    'H-rich_non_burning',
    'H-rich_Shell_H_burning',
    ]


BINARY_STATES_CC = BINARY_STATES_ALL.copy()
#BINARY_STATES_CC = BINARY_STATES_ALL.copy()
#BINARY_STATES_CC.remove('disrupted')

for b in BINARY_STATES_CC:
    for s1 in STAR_STATES_CC:
        for s2 in STAR_STATES_ALL:
            POSYDON_FLOW_CHART[(s1, s2, b, 'CC1')] = 'step_SN'
            POSYDON_FLOW_CHART[(s2, s1, b, 'CC2')] = 'step_SN'

# Double compact objects. These can either be send to the orbital evolution
# due to GR step, or end the evolution by setting the 'step_dco' to end
for s1 in STAR_STATES_CO:
    for s2 in STAR_STATES_CO:
        POSYDON_FLOW_CHART[(s1, s2, 'detached', None)] = 'step_dco'
        POSYDON_FLOW_CHART[(s2, s1, 'detached', None)] = 'step_dco'


# catch states to be ended
for b in ['initial_RLOF',
          'detached (Integration failure)',
          'detached (GridMatchingFailed)', 'RLO2 (OutsideGrid)']:
    for s1 in STAR_STATES_ALL:
        for s2 in STAR_STATES_ALL:
            for e in BINARY_EVENTS_ALL:
                POSYDON_FLOW_CHART[(s1, s2, b, e)] = 'step_end'
                POSYDON_FLOW_CHART[(s1, s2, b, e)] = 'step_end'

for b in ['initially_single_star']:
    for s1 in STAR_STATES_ALL:
        for s2 in STAR_STATES_ALL:
            for e in BINARY_EVENTS_ALL:
                if e == 'ZAMS':
                    POSYDON_FLOW_CHART[(s1, s2, b, e)] = 'step_initially_single'
                    POSYDON_FLOW_CHART[(s2, s1, b, e)] = 'step_initially_single'


for s1 in STAR_STATES_CO:
    for s2 in ['massless_remnant']:
        for e in BINARY_EVENTS_ALL:
            POSYDON_FLOW_CHART[(s1, s2, b, e)] = 'step_end'
            POSYDON_FLOW_CHART[(s2, s1, b, e)] = 'step_end'

BINARY_EVENTS_OF_SN_OR_AFTER_DETACHED = BINARY_EVENTS_ALL.copy()
[BINARY_EVENTS_OF_SN_OR_AFTER_DETACHED.remove(x) for x in ['CC1','CC2','MaxTime_exceeded','maxtime']]

for b in ['disrupted']:
    for s1 in STAR_STATES_ALL:
        for s2 in STAR_STATES_ALL:
            for e in BINARY_EVENTS_OF_SN_OR_AFTER_DETACHED:
                POSYDON_FLOW_CHART[(s1, s2, b, e)] = 'step_disrupted'
                POSYDON_FLOW_CHART[(s2, s1, b, e)] = 'step_disrupted'
# if we have two compcat objects in a disrupted binary, we stop the evolution.
for b in ['disrupted']:
    for s1 in STAR_STATES_CO:
        for s2 in STAR_STATES_CO:
            for e in BINARY_EVENTS_ALL:
                POSYDON_FLOW_CHART[(s1, s2, b, e)] = 'step_end'
                POSYDON_FLOW_CHART[(s2, s1, b, e)] = 'step_end'


for b in ['merged']:
    for s1 in STAR_STATES_ALL:
        for s2 in STAR_STATES_ALL:
            for e in ['oMerging1', 'oMerging2']:
                POSYDON_FLOW_CHART[(s1, s2, b, e)] = 'step_end' #'step_merged'
                POSYDON_FLOW_CHART[(s2, s1, b, e)] = 'step_end' #'step_merged'

# catch initial_RLO states
for s1 in STAR_STATES_ALL:
    for s2 in STAR_STATES_ALL:
        POSYDON_FLOW_CHART[(s1, s2, 'initial_RLOF', 'ZAMS')] = 'step_end'
        POSYDON_FLOW_CHART[(s2, s1, 'initial_RLOF', 'ZAMS')] = 'step_end'

# catch all maxtime
for b in BINARY_STATES_ALL:
    for s1 in STAR_STATES_ALL:
        for s2 in STAR_STATES_ALL:
            for e in ['maxtime', 'MaxTime_exceeded', 'CO_contact']:
                POSYDON_FLOW_CHART[(s1, s2, b, e)] = 'step_end'


def flow_chart(FLOW_CHART=POSYDON_FLOW_CHART, CHANGE_FLOW_CHART=None):
    """Generate the flow chart.

    Default step nomenclature:
    - step_SN :  StepSN in posydon/bianry_evol/SN/step_SN.py
    - step_HMS_HMS : MS_MS_step in posydon/binary_evol/mesa_step.py
    - step_CO_HMS_RLO : CO_HMS_RLO_step in posydon/binary_evol/mesa_step.py
    - step_CO_HeMS : CO_HeMS_step in posydon/binary_evol/mesa_step.py
    - step_detached : detached_step in posydon/binary_evol/detached_step.py
    - step_CE : StepCEE in posydon/binary_evol/CE/step_CEE.py
    - step_dco : DoubleCO in posydon/binary_evol/double_CO.py
    - step_end


    Parameters
    ----------
    FLOW_CHART : dict
        Flow chart that maps the tuple (star_1_state, star_1_state,
        binary_state, binary_event) to a step.
    CHANGE_FLOW_CHART : dict
        Flow chart to change the default FLOW_CHART. E.g.
        CHANGE_FLOW_CHART = {('NS', 'NS', 'detached', None) : 'step_end',
                             ('NS', 'HeMS', 'RLO1', 'pRLO1') : 'step_my_RLO1'}

    Returns
    -------
    dict
        Flow chart.

    """
    if CHANGE_FLOW_CHART is not None:
        for key in CHANGE_FLOW_CHART.keys():
            if key in FLOW_CHART.keys():
                # this assume CHANGE_FLOW_CHART[key] = 'step_XXX'
                FLOW_CHART[key] = CHANGE_FLOW_CHART[key]

    return FLOW_CHART
