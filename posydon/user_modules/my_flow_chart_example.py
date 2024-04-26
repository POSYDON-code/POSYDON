import posydon.binary_evol.flow_chart as fc
from posydon.binary_evol.flow_chart import flow_chart

def modified_flow_chart(FLOW_CHART=None, CHANGE_FLOW_CHART=None):
    """Derive a flow chart from POSYDON's standard one

    Parameters
    ----------
    FLOW_CHART : dict
    CHANGE_FLOW_CHART : dict or None
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

    return MY_FLOW_CHART
