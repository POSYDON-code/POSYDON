"""Merging and isolated evolution step."""


__authors__ = [
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>"
]


from posydon.binary_evol.DT.step_isolated import IsolatedStep
from posydon.binary_evol.flow_chart import STAR_STATES_H_RICH, STAR_STATES_HE_RICH
from posydon.config import PATH_TO_POSYDON_DATA

LIST_ACCEPTABLE_STATES_FOR_HMS = ["H-rich_Core_H_burning"]
LIST_ACCEPTABLE_STATES_FOR_HeMS = ["stripped_He_Core_He_burning"]

LIST_ACCEPTABLE_STATES_FOR_POSTMS = STAR_STATES_H_RICH.copy()
[LIST_ACCEPTABLE_STATES_FOR_POSTMS.remove(x) for x in LIST_ACCEPTABLE_STATES_FOR_HMS]

LIST_ACCEPTABLE_STATES_FOR_POSTHeMS = STAR_STATES_HE_RICH.copy()
[LIST_ACCEPTABLE_STATES_FOR_POSTHeMS.remove(x) for x in LIST_ACCEPTABLE_STATES_FOR_HeMS]

class InitiallySingleStep(IsolatedStep):
    """
    Prepare a runaway star to do an an isolated_step)
    """

    def __init__(self,
        grid_name_Hrich=None,
        grid_name_strippedHe=None,
        path=PATH_TO_POSYDON_DATA,
        *args, **kwargs):

        super().__init__(
        grid_name_Hrich=grid_name_Hrich,
        grid_name_strippedHe=grid_name_strippedHe,
        *args,
        **kwargs)


    def __call__(self,binary):
        if binary.state != "initially_single_star":
            raise AttributeError("sent to InitiallySingleStep without the binary.state being initially_single_star")

        binary.event == None
        super().__call__(binary)
