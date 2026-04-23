"""An example of a modified function of the flow chart."""


__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]

import posydon.binary_evol.flow_chart as fc
from posydon.binary_evol.flow_chart import STAR_STATES_NORMALSTAR, flow_chart


def modified_flow_chart(FLOW_CHART=None, CHANGE_FLOW_CHART=None):
    """Derive a flow chart from POSYDON's standard one and modify it.

    Parameters
    ----------
    FLOW_CHART : dict or None
    CHANGE_FLOW_CHART : dict or None

    Returns
    -------
    dict
        Modified flow chart.
    """

    # get POSYDON's default flow chart
    if FLOW_CHART is None:
        if CHANGE_FLOW_CHART is None:
            MY_FLOW_CHART = fc.flow_chart()
        else:
            MY_FLOW_CHART = fc.flow_chart(CHANGE_FLOW_CHART=CHANGE_FLOW_CHART)
    else:
        if CHANGE_FLOW_CHART is None:
            MY_FLOW_CHART = fc.flow_chart(FLOW_CHART=FLOW_CHART)
        else:
            MY_FLOW_CHART = fc.flow_chart(FLOW_CHART=FLOW_CHART,
                                          CHANGE_FLOW_CHART=CHANGE_FLOW_CHART)
    # modify the flow chart
    for key in MY_FLOW_CHART.keys():
        # here changing to go always from ZAMS to step_detached
        if ((len(key)>3) and (key[3]=='ZAMS') and
            (key[2] in fc.BINARY_STATES_ZAMS)): # check for event
            MY_FLOW_CHART[key] = 'step_detached'
        # here changing all starting RLO to end instead
        if ((len(key)>2) and ('RLO' in key[2])): # check for state
            MY_FLOW_CHART[key] = 'step_end'
    # adding new entries to the flow to send RLO to end which would otherwise
    # be jumped over in the HMS-HMS step
    for s1 in STAR_STATES_NORMALSTAR:
        for s2 in STAR_STATES_NORMALSTAR:
            MY_FLOW_CHART[(s1, s2, 'RLO1', 'oRLO1')] = 'step_end'
            MY_FLOW_CHART[(s2, s1, 'RLO1', 'oRLO1')] = 'step_end'
            MY_FLOW_CHART[(s1, s2, 'RLO2', 'oRLO2')] = 'step_end'
            MY_FLOW_CHART[(s2, s1, 'RLO2', 'oRLO2')] = 'step_end'

    return MY_FLOW_CHART
