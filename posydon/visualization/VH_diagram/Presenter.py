"""Data visualization in VH diagram."""


__authors__ = [
    "Maxime Rambosson <Maxime.Rambosson@etu.unige.ch>",
]


from .SimulationModel import SimulationModel
from .GraphVisualizer import CaseInfos, PointInfos, StateInfos, columnTYPE

from .MainWindow import MainWindow
from .PresenterMode import PresenterMode

from posydon.utils.common_functions import orbital_separation_from_period
from posydon.utils.common_functions import PATH_TO_POSYDON

from datetime import datetime
import os
import numpy as np


def to_megayears(nb):
    """Convert yr to Myr."""
    return nb / 1000000


def combine_simplified_data(*, data_to, data_from):
    """Combine 2 dictionaries of SimplifiedInfo.

    data_to[].state_after take the value of data_from[].state_before.

    Parameters
    ----------
    data_to : dictionary of SimplifiedInfo
        Informations who will be udpated.
    data_from : dictionary of SimplifiedInfo
        Informations who will be combined.

    """
    for key in data_to.keys():
        data_to[key].state_after = data_from[key].state_before


def to_simplified_data(row):
    """Convert dictionnary of SimplifiedInfo with values as `.state_before`.

    Parameters
    ----------
    row : dictionnary
        Dictionnary of value to convert.

    Returns
    -------
    dictionnary of SimplifiedInfo
        Dictionnary of SimplifiedInfo containing the value as
        SimplifiedInfo.state_before.

    """
    data = {}

    for key in row.keys():
        data[key] = SimplifiedInfo(row[key])

    # Side effect
    data["state"] = SimplifiedInfo(
        row["event"] if (str(row["event"]) != "nan") else row["state"]
    )

    return data


def file_exist(filename):
    """Check if file exists."""
    return os.path.exists(filename)


def get_star_state_filename(state, *, suffix=""):
    """Get the name of the file illustring the state of one star.

    Parameters
    ----------
    state : str
        State name of the star.

    Returns
    -------
    str
        Name of the file illustrating this state.

    """
    if state.count("_") < 2:
        return state + suffix

    parts = state.split("_")

    return parts[0] + "_" + suffix


def get_event_state_filename(left_star_state, event_state, right_star_state, *,
                             suffix=""):
    """Get the name of the file illustrating the event with the 2 given star.

    Parameters
    ----------
    left_star_state : str
        State's name of the left star.
    event_state : str
        State's name of the event.
    right_star_state : type
        State's name of the right star.

    Returns
    -------
    str
        Name of the file illustrating this complete state.

    """
    left_star_filename = get_star_state_filename(left_star_state)
    if left_star_filename[-1] != "_":
        left_star_filename += "_"

    return (
        left_star_filename
        + event_state
        + "_"
        + get_star_state_filename(right_star_state)
        + suffix
    )


def equal_with_epsilon(nb_1, nb_2, epsilon):
    """Check if 2 numbers are equal up to epsilon.

    Parameters
    ----------
    nb_1 : int, float, double
        Value 1 to compare.
    nb_2 : int, float, double
        Value 2 compare.
    epsilon : int, float, double
        Difference allowed.

    Returns
    -------
    bool
        True if numbers are equal up to epsilon, False if not.

    """
    return abs(nb_1 - nb_2) < epsilon


def get_max_distance(data):
    """Get max distance from a dictionary of SimplifiedInfo.

    Parameters
    ----------
    data : dictionnary of SimplifiedInfo
        Dictionnary, full of SimplifiedInfo, containing all infos

    Returns
    -------
    float
        Max distance of this row, accroding to before/after state

    """
    if data["separation"].state_after is not None:
        if np.isnan(data["separation"].state_after) or (
            data["separation"].state_before > data["separation"].state_after
        ):
            return (
                data["separation"].state_before
                if (not np.isnan(data["separation"].state_before))
                else 0
            )
        else:
            return data["separation"].state_after
    else:
        if np.isnan(data["separation"].state_before):
            return 0
        else:
            return data["separation"].state_before


def calculate_separation(data):
    """For each row in dataframe, calculate the separation.

    Parameters
    ----------
    data : pandas.DataFrame
        Complete dataframe without 'separation'.

    Returns
    -------
    Array of float
        Array of calculated separations, len(array) == nb rows.

    """
    if (("orbital_period" not in data) or ("S1_mass" not in data)
            or ("S2_mass" not in data)):
        return None

    separations = []

    for index, row in data.iterrows():
        separations.append(
            orbital_separation_from_period(row["orbital_period"],
                                           row["S1_mass"], row["S2_mass"]))
    return separations


class SimplifiedInfo:
    """Simple class to store an before/after info and simplfy dataframe.

    Parameters
    ----------
    state_before : any value
        Value for the older state.
    state_after : any value
        Value for the younger state.

    Attributes
    ----------
    state_before
    state_after
    """

    def __init__(self, state_before=None, state_after=None):
        """Initialize a SimplifiedInfo instance."""
        self.state_before = state_before
        self.state_after = state_after


class Presenter:
    """Charged to setup window and format data before calling display funcs."""

    # PATH_TO_DRAWS = os.path.dirname(__file__) + "\\draws\\%s.png"
    PATH_TO_DRAWS = os.path.join(
        PATH_TO_POSYDON, "posydon", "visualization", "VH_diagram", "draws"
    )
    PATH_TO_SCREENS = os.path.join(os.getcwd(), "screens")

    _formating_dict = {
        "binary_index": "",
        "state": "",
        "event": "",
        "time": r"$%.2f\, Mys$",
        "separation": "",
        "orbital_period": r"$p = %.2f\, d$",
        "eccentricity": r"$e = %.2f $",
        "V_sys": "",
        "rl_relative_overflow_1": "",
        "rl_relative_overflow_2": "",
        "lg_mtransfer_rate": "",
        "mass_transfer_case": "",
        "nearest_neighbour_distance": "",
        "S1_state": "",
        "S1_metallicity": "",
        "S1_mass": r"$m=%3.2f\, \mathrm{M}_\odot$",
        "S1_log_R": "",
        "S1_log_L": "",
        "S1_lg_mdot": "",
        "S1_lg_system_mdot": "",
        "S1_lg_wind_mdot": "",
        "S1_he_core_mass": "",
        "S1_he_core_radius": "",
        "S1_c_core_mass": "",
        "S1_c_core_radius": "",
        "S1_o_core_mass": "",
        "S1_o_core_radius": "",
        "S1_center_h1": "",
        "S1_center_he4": "",
        "S1_center_c12": "",
        "S1_center_n14": "",
        "S1_center_o16": "",
        "S1_surface_h1": "",
        "S1_surface_he4": "",
        "S1_surface_c12": "",
        "S1_surface_n14": "",
        "S1_surface_o16": "",
        "S1_log_LH": "",
        "S1_log_LHe": "",
        "S1_log_LZ": "",
        "S1_log_Lnuc": "",
        "S1_c12_c12": "",
        "S1_avg_c_in_c_core": "",
        "S1_surf_avg_omega": "",
        "S1_surf_avg_omega_div_omega_crit": "",
        "S1_total_moment_of_inertia": "",
        "S1_log_total_angular_momentum": "",
        "S1_spin": "",
        "S1_profile": "",
        "S2_state": "",
        "S2_metallicity": "",
        "S2_mass": r"$m=%3.2f\, \mathrm{M}_\odot$",
        "S2_log_R": "",
        "S2_log_L": "",
        "S2_lg_mdot": "",
        "S2_lg_system_mdot": "",
        "S2_lg_wind_mdot": "",
        "S2_he_core_mass": "",
        "S2_he_core_radius": "",
        "S2_c_core_mass": "",
        "S2_c_core_radius": "",
        "S2_o_core_mass": "",
        "S2_o_core_radius": "",
        "S2_center_h1": "",
        "S2_center_he4": "",
        "S2_center_c12": "",
        "S2_center_n14": "",
        "S2_center_o16": "",
        "S2_surface_h1": "",
        "S2_surface_he4": "",
        "S2_surface_c12": "",
        "S2_surface_n14": "",
        "S2_surface_o16": "",
        "S2_log_LH": "",
        "S2_log_LHe": "",
        "S2_log_LZ": "",
        "S2_log_Lnuc": "",
        "S2_c12_c12": "",
        "S2_avg_c_in_c_core": "",
        "S2_surf_avg_omega": "",
        "S2_surf_avg_omega_div_omega_crit": "",
        "S2_total_moment_of_inertia": "",
        "S2_log_total_angular_momentum": "",
        "S2_spin": "",
        "S2_profile": "",
    }

    def __init__(self, filename, path="./"):
        """Initialize a Presenter instance."""
        self._model = SimulationModel(filename=filename, path=path)
        self._model.load_csv()

        self._main_window = MainWindow()

        self._current_index = None
        self._current_data = None

        self._visualizer = self._main_window.start_visualization()
        self._visualizer.detail_required().connect(
            self._set_visualisation_detailed)
        self._visualizer.reduce_required().connect(
            self._set_visualisation_reduced)
        self._visualizer.simplify_required().connect(
            self._set_visualisation_simplified)
        self._visualizer.diagram_required().connect(
            self._set_visualisation_diagram)
        self._visualizer.distance_representation_required().connect(
            self._set_distance_representation)
        breakpoint()
        self._visualizer.save_required().connect(self.screen)

        self._present_mode = PresenterMode.DETAILED
        self._visualizer.options().set_showed_mode(self._present_mode)

        # maybe change it by flags if more option
        self._distance_representation = True
        self._visualizer.options().set_distance_representation(
            self._distance_representation
        )

    def screen(self):
        """Take a screen of the displayed window."""
        os.makedirs(self.PATH_TO_SCREENS, exist_ok=True)

        filepath = os.path.join(
            self.PATH_TO_SCREENS,
            datetime.today().strftime("%Y-%m-%d-%H-%M-%S")
            + "_" + str(self._current_index) + ".png")

        self._visualizer().update()
        self._visualizer().repaint()
        self._visualizer().saveAsPicture(filepath)

        return filepath

    def _prepare_corresponding_data(self, index):
        """Get binary data and calculate missing columns.

        Parameters
        ----------
        index : int
            Binary index to load.

        """
        if self._current_index == index:
            return

        self._current_data = self._model.get_by_binary_index(index)

        if "separation" not in self._current_data:
            separations = calculate_separation(self._current_data)
            if separations is not None:
                self._current_data["separation"] = separations

        self._current_index = index

    def present(self, index, mode=PresenterMode.DETAILED):
        """Preset the binary."""
        self._present_mode = mode
        self._visualizer.options().set_showed_mode(self._present_mode)

        self._prepare_corresponding_data(index)
        self._update_visualisation()

        self._main_window.show()

    def close(self):
        """Close the window."""
        self._main_window.close()

    def _update_visualisation(self):
        """Update the display according to current_index and present_mode."""
        if self._current_index is None:
            return

        self._visualizer().reset()

        if self._present_mode == PresenterMode.DETAILED:
            self._prepare_basic_columns()
            self._detailed_presentation(self._current_data)
        elif self._present_mode == PresenterMode.REDUCED:
            self._prepare_basic_columns()
            self._reduced_presentation(self._current_data)
        elif self._present_mode == PresenterMode.SIMPLIFIED:
            self._prepare_basic_columns()
            self._simplified_presentation(self._current_data)
        elif self._present_mode == PresenterMode.DIAGRAM:
            self._prepare_diagram_columns()
            self._digram_presentation(self._current_data)

    def _prepare_basic_columns(self):
        """Add columns for Detailled/Reduced/Simplified View."""
        self._time_id = self._visualizer().add_column(columnTYPE.TIMELINE)
        self._visualizer().get_column(self._time_id).set_title("TIME")

        self._S1_id = self._visualizer().add_column(columnTYPE.CONNECTED)
        self._visualizer().get_column(self._S1_id).set_title("S1")

        self._event_id = self._visualizer().add_column(columnTYPE.CONNECTED)
        self._visualizer().get_column(self._event_id).set_title("EVENT/STATE")

        self._S2_id = self._visualizer().add_column(columnTYPE.CONNECTED)
        self._visualizer().get_column(self._S2_id).set_title("S2")

    def _prepare_diagram_columns(self):
        """Add columns for Diagram View."""
        self._time_id = self._visualizer().add_column(columnTYPE.TIMELINE)
        self._visualizer().get_column(self._time_id).set_title("TIME")

        self._state_id = self._visualizer().add_column(columnTYPE.CONNECTED)
        self._visualizer().get_column(self._state_id).set_title("STATE")

    def _set_visualisation_detailed(self):
        self._present_mode = PresenterMode.DETAILED
        self._update_visualisation()

    def _set_visualisation_reduced(self):
        self._present_mode = PresenterMode.REDUCED
        self._update_visualisation()

    def _set_visualisation_simplified(self):
        self._present_mode = PresenterMode.SIMPLIFIED
        self._update_visualisation()

    def _set_visualisation_diagram(self):
        self._present_mode = PresenterMode.DIAGRAM
        self._update_visualisation()

    def _set_distance_representation(self, state):
        self._distance_representation = state
        self._update_visualisation()

    def _prepare_basic_data(self, row):
        """Prepare widget to display as basic (non-simplified) data.

        Parameters
        ----------
        row : dictionary
            Dataframe's row needing to be formated.

        Returns
        -------
        Array of Infos
            Array containing derivated struct of Infos to create correspondong
            widget.

        """
        time_data = PointInfos(
            self._time_id,
            self._formating_dict["time"] % to_megayears(row["time"]))

        s1_data = CaseInfos(
            self._S1_id,
            centered_txt=row["S1_state"],
            bot_left_txt=self._formating_dict["S1_mass"] % row["S1_mass"],
        )

        s1_data.border_width = 2
        event_data = CaseInfos(self._event_id)
        event_data.border_width = 2
        event_data.centered_text = (str(row["event"])
                                    if (str(row["event"]) != "nan")
                                    else str(row["state"]))

        if "eccentricity" in row:
            event_data.bot_right_text = (
                self._formating_dict["eccentricity"] % row["eccentricity"]
            )

        if "orbital_period" in row:
            event_data.bot_left_text = (
                self._formating_dict["orbital_period"] % row["orbital_period"]
            )

        s2_data = CaseInfos(
            self._S2_id,
            centered_txt=row["S2_state"],
            bot_right_txt=self._formating_dict["S2_mass"] % row["S2_mass"],
        )

        s2_data.border_width = 2

        return [time_data, s1_data, event_data, s2_data]

    def _set_simplified_infos(self, dataset, name, infos, side=False):
        """Format the information 'data_name' for simplfied view.

        This is done according to whether it exists or not, and if there is a
        before/after.

        Parameters
        ----------
        dataset : dictionary of SimplifiedInfo
            Dictionary containing all infos to format.
        name : str
            Name of the info to format.
        infos : type
            Struct to create coresponding widget.
        side : bool
            Indicate on which side display the information, False = Left

        """
        if name not in dataset:
            return

        if dataset[name].state_after is not None:
            if not side:
                infos.bot_left_text = (self._formating_dict[name]
                                       % dataset[name].state_after)
                infos.top_left_text = (self._formating_dict[name]
                                       % dataset[name].state_before)
            else:
                infos.bot_right_text = (self._formating_dict[name]
                                        % dataset[name].state_after)
                infos.top_right_text = (self._formating_dict[name]
                                        % dataset[name].state_before)
        else:
            if not side:
                infos.bot_left_text = (self._formating_dict[name]
                                       % dataset[name].state_before)
            else:
                infos.bot_right_text = (self._formating_dict[name]
                                        % dataset[name].state_before)

    def _prepare_simplified_line(self, data):
        """Prepare widget informations to display as simplified data.

        Parameters
        ----------
        data : dictionary of SimplifiedInfo
            Infos formated to be displayed in widgets.

        Returns
        -------
        Array of Infos
            Array containing derivated struct of Infos to create corresponding
            widgets.

        """
        time_data = PointInfos(
            self._time_id,
            self._formating_dict["time"]
            % to_megayears(data["time"].state_before))

        s1_data = CaseInfos(self._S1_id)
        s1_data.border_width = 2

        if data["S1_state"].state_after is not None:
            s1_data.centered_text = data["S1_state"].state_after
        else:
            s1_data.centered_text = data["S1_state"].state_before

        self._set_simplified_infos(data, "S1_mass", s1_data, False)

        event_data = CaseInfos(self._event_id)
        event_data.border_width = 2
        if data["state"].state_after is not None:
            event_data.centered_text = data["state"].state_after
        else:
            event_data.centered_text = data["state"].state_before

        self._set_simplified_infos(data, "eccentricity", event_data, False)
        self._set_simplified_infos(data, "orbital_period", event_data, True)

        s2_data = CaseInfos(self._S2_id)
        s2_data.border_width = 2
        if data["S2_state"].state_after is not None:
            s2_data.centered_text = data["S2_state"].state_after
        else:
            s2_data.centered_text = data["S2_state"].state_before

        self._set_simplified_infos(data, "S2_mass", s2_data, True)

        return [time_data, s1_data, event_data, s2_data]

    def _get_distance_representation(self, simplified_distance, max_distance):
        """Normalize 1 of the 2 distances in simplified_distance.

        Use the max_distance, according to whether it exists or not,
        if it's NaN or 0.

        Parameters
        ----------
        simplified_distance : SimplifiedInfo
            Struct containing distances of the current state.
        max_distance : float
            Description of parameter `max_distance`.

        Returns
        -------
        float
            Normalized current star distance.

        """
        if max_distance == 0:
            return 0

        if simplified_distance.state_after is not None and not np.isnan(
            simplified_distance.state_after
        ):
            if simplified_distance.state_after < 1:
                simplified_distance.state_after = 1
            if max_distance < 1:
                return np.log(simplified_distance.state_after)

            return (np.log(simplified_distance.state_after)
                    / np.log(max_distance))

        else:
            if np.isnan(simplified_distance.state_before):
                return 0
            elif simplified_distance.state_before < 1:
                simplified_distance.state_before = 1

            if max_distance < 1:
                return np.log(simplified_distance.state_before)

            return (np.log(simplified_distance.state_before)
                    / np.log(max_distance))

    def _set_diagram_infos(self, dataset, data_name, diagram_infos):
        """Format the information 'data_name' for diagram view.

        This is done according to whether it exists or not, and if there is
        a before/after.

        Parameters
        ----------
        dataset : dictionary of SimplifiedInfo
            Dictionary containing all infos to format.
        data_name : str
            Name of the info to format.
        diagram_infos : StateInfos
            Struct to create coresponding widget.

        """
        if data_name not in dataset:
            diagram_infos.top_texts.append("")
            diagram_infos.bot_texts.append("")
            return

        diagram_infos.top_texts.append(self._formating_dict[data_name]
                                       % dataset[data_name].state_before)

        if dataset[data_name].state_after is not None:
            diagram_infos.bot_texts.append(self._formating_dict[data_name]
                                           % dataset[data_name].state_after)
        else:
            diagram_infos.bot_texts.append("")

    def _prepare_diagram_line(self, data, max_distance):
        """Prepare widget information to display basic (non-simplified) data.

        Parameters
        ----------
        data : dictionary of SimplifiedInfo
            Infos formated to be displayed in widgets.
        max_distance : float
            Max distance between 2 star, used to set normalized distance
            for this state.

        Returns
        -------
        Array of Infos
            Array containing derivated struct of Infos to create corresponding
            widgets.

        """
        time_info = PointInfos(
            self._time_id,
            self._formating_dict["time"] % to_megayears(
                data["time"].state_before),
        )

        state_info = StateInfos(self._state_id)
        if self._distance_representation and "separation" in data:
            state_info.distance = self._get_distance_representation(
                data["separation"], max_distance
            )
        else:
            state_info.distance = 1

        state_before_filename = os.path.join(
            self.PATH_TO_DRAWS,
            get_event_state_filename(
                data["S1_state"].state_before,
                data["state"].state_before,
                data["S2_state"].state_before,
                suffix=".png",
            ),
        )
        state_after_filename = os.path.join(
            self.PATH_TO_DRAWS,
            get_event_state_filename(
                data["S1_state"].state_before,
                data["state"].state_before,
                data["S2_state"].state_before,
                suffix=".png",
            ),
        )
        if (data["state"].state_after != "detached"
                and file_exist(state_after_filename)):
            state_info.event_filename = state_after_filename
        elif data["state"].state_before != "detached" and file_exist(
            state_before_filename
        ):
            state_info.event_filename = state_before_filename
        else:
            if data["S1_state"].state_after is not None:
                state_info.S1_filename = os.path.join(
                    self.PATH_TO_DRAWS,
                    get_star_state_filename(
                        data["S1_state"].state_after, suffix=".png"
                    ),
                )
            else:
                state_info.S1_filename = os.path.join(
                    self.PATH_TO_DRAWS,
                    get_star_state_filename(
                        data["S1_state"].state_before, suffix=".png"
                    ),
                )

            if data["S2_state"].state_after is not None:
                state_info.S2_filename = os.path.join(
                    self.PATH_TO_DRAWS,
                    get_star_state_filename(
                        data["S2_state"].state_after, suffix=".png"
                    ),
                )
            else:
                state_info.S2_filename = os.path.join(
                    self.PATH_TO_DRAWS,
                    get_star_state_filename(
                        data["S2_state"].state_before, suffix=".png"
                    ),
                )

        self._set_diagram_infos(data, "S1_mass", state_info)

        self._set_diagram_infos(data, "eccentricity", state_info)
        self._set_diagram_infos(data, "orbital_period", state_info)

        self._set_diagram_infos(data, "S2_mass", state_info)

        return [time_info, state_info]

    def _detailed_presentation(self, data):
        """Format and display data in the Detailled View.

        Parameters
        ----------
        data : pandas.DataFrame
            Data to display in detailled view.

        """
        for index, row in data.iterrows():
            if row["event"] == "END":
                continue

            formatted_data = self._prepare_basic_data(row)
            self._visualizer().add_line(formatted_data)

    def _reduce_data(self, data):
        """Reduce the dataframe given by folowing hardcoded rules.

        If the row has same time (within 10,000 years) as the previous one,
        each info is combined with the previous info in a SimplifiedInfo struct
        giving a dictionary with same row's key, full of SimplifiedInfo, which
        are the simplified row.

        Parameters
        ----------
        data : pandas.DataFrame
            Data to simplify.

        Returns
        -------
        Array of dictionary of SimplifiedInfo
            Each dictionary hold same keys of 1 data's row, but with
            SimplifiedInfo as value.
            Array length =/= row's number in data because some infos are
            combined in SimplifiedInfo.
        """
        reduced_data = []

        for index, row in data.iterrows():

            if row["event"] == "END" or row["event"] == "redirect":
                continue

            new_data = to_simplified_data(row)

            if not reduced_data:
                reduced_data.append(new_data)
            else:
                if (
                    equal_with_epsilon(
                        reduced_data[-1]["time"].state_before,
                        new_data["time"].state_before,
                        10000,
                    )
                    and (reduced_data[-1]["state"].state_before != "detached")
                    and (new_data["state"].state_before != "detached")
                ):  # Need to combine
                    combine_simplified_data(
                        data_to=reduced_data[-1], data_from=new_data
                    )
                else:
                    reduced_data.append(new_data)

        return reduced_data

    def _reduced_presentation(self, data):
        """Format and display data in the Reduced View.

        Parameters
        ----------
        data : pandas.DataFrame
            Data to display in reduced view.

        """
        reduced_data = self._reduce_data(data)

        for line_data in reduced_data:
            self._visualizer().add_line(
                self._prepare_simplified_line(line_data))

    def _simplify_data(self, data):
        """Simplify the dataframe given by folowing hardcoded rules.

        If the row isn't interesting, each info is combined with the previous
        info in a SimplifiedInfo struct, giving a dictionary with same row's
        key, full of SimplifiedInfo, which are the simplified row.

        Parameters
        ----------
        data : pandas.DataFrame
            Data to simplify.

        Returns
        -------
        Array of dictionary of SimplifiedInfo
            Each dictionary hold same keys of 1 data's row, but with
            SimplifiedInfo as value.
            Array length =/= row's number in data because some infos are
            combined in SimplifiedInfo.

        """
        simplified_data = []

        for index, row in data.iterrows():
            if row["event"] == "END" or row["event"] == "redirect":
                continue

            new_data = to_simplified_data(row)

            if not simplified_data:
                simplified_data.append(new_data)
            else:
                if (
                    equal_with_epsilon(
                        simplified_data[-1]["time"].state_before,
                        new_data["time"].state_before,
                        10000,
                    )
                    and simplified_data[-1]["state"].state_before != "detached"
                    and new_data["state"].state_before != "detached"
                ) or (
                    simplified_data[-1]["state"].state_before == "detached"
                    and new_data["state"].state_before == "detached"
                ):  # Need to combine
                    combine_simplified_data(
                        data_to=simplified_data[-1], data_from=new_data
                    )
                else:
                    simplified_data.append(new_data)

        return simplified_data

    def _simplified_presentation(self, data):
        """Format and display data in the Simplified View.

        Parameters
        ----------
        data : pandas.DataFrame
            Data to display in simplified view.

        """
        simplified_data = self._simplify_data(data)

        for data in simplified_data:
            self._visualizer().add_line(self._prepare_simplified_line(data))

    def _digram_presentation(self, data):
        """Format and display data in the Diagram View.

        Parameters
        ----------
        data : pandas.DataFrame
            Data to display in diagram view.

        """
        simplified_data = self._simplify_data(data)

        max_distance_data = max(simplified_data, key=get_max_distance)
        max_distance = get_max_distance(max_distance_data)

        for line_data in simplified_data:
            self._visualizer().add_line(
                self._prepare_diagram_line(line_data, max_distance)
            )

        if (
            simplified_data[-1]["state"].state_after == "disrupted"
            or simplified_data[-1]["state"].state_after == "merged"
        ):
            aditional_info = CaseInfos(self._state_id)
            aditional_info.border_width = 2
            aditional_info.centered_text = simplified_data[-1][
                "state"].state_after
            aditional_info.connected = False
            self._visualizer().add_line([aditional_info])