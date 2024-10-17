from posydon.binary_evol.flow_chart import (flow_chart, STAR_STATES_NORMALSTAR)

def end_flow_chart(FLOW_CHART=None):
    FLOW_CHART = {}
    for s1 in STAR_STATES_NORMALSTAR:
        for s2 in STAR_STATES_NORMALSTAR:
            FLOW_CHART[(s1, s2, 'detached', 'ZAMS')] = 'step_end'
    return FLOW_CHART



class my_CE_step(object):
    """Compute a fake CE event."""

    def __init__(self, verbose=False):
        self.verbose = verbose

    def __call__(self, binary):

        if self.verbose:
            print('The orbital separation post CE is half the pre CE orbital separation!')

        # Determine which star is the donor and which is the companion
        if binary.event in ["oCE1", "oDoubleCE1"]:
            donor_star = binary.star_1
            comp_star = binary.star_2
        elif binary.event in ["oCE2", "oDoubleCE2"]:
            donor_star = binary.star_2
            comp_star = binary.star_1
        else:
            raise ValueError("CEE does not apply if `event` is not "
                            "`oCE1`, 'oDoubleCE1' or `oCE2`, 'oDoubleCE1'")

        binary.orbital_period /= 2.
        donor_star.mass = donor_star.he_core_mass # lose envelope
        donor_star.state = donor_star.state.replace('H-rich', 'stripped_He')
        binary.state = 'detached'
        binary.event = None