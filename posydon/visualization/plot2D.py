"""Plotting class for 2D (MESA) psygrids.

The 2D visualization plotting class allows to plot 2D slices of PsyGrid
objects. The PsyGrid object is composed of 2D/3D/4D MESA grid run with POSYDON
and post processed with the psygrid object into an h5 file.
"""


__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
]


import numpy as np
import matplotlib.pyplot as plt
from posydon.utils.gridutils import add_field
from posydon.visualization.plot_defaults import DEFAULT_MARKERS_COLORS_LEGENDS
from posydon.visualization.plot_defaults import PLOT_PROPERTIES
from posydon.visualization.plot_defaults import DEFAULT_LABELS
from posydon.visualization.combine_TF import combine_TF12
import copy


class plot2D(object):
    """Plotting class for 2D (MESA) grids."""

    def __init__(
        self,
        psygrid,
        x_var_str,
        y_var_str,
        z_var_str=None,
        selected_star_history_for_z_var=1,
        termination_flag="termination_flag_1",
        grid_3D=False,
        slice_3D_var_str=None,
        slice_3D_var_range=None,
        grid_4D=False,
        slice_4D_var_str=None,
        slice_4D_var_range=None,
        extra_grid=None,
        slice_at_RLO=False,
        MARKERS_COLORS_LEGENDS=None,
        verbose=False,
        **kwargs
    ):
        """Read a PsyGrid object and plot a 2D slice of x vs y.

        Parameters
        ----------
        psygrid : object
            PsyGrid object containing a 2D/3D/4D MESA grid.
        x_var_str : str
            String of the initial value to plot on the x axis. Allowed strings
            are `psygrid.initial_values.dtype.names`.
        y_var_str : str
            String of the initial value to plot on the y axis. Allowed strings
            are `psygrid.initial_values.dtype.names`.
        z_var_str : str
            String of the initial value to plot on the z axis (displayed as
            a color). Allowed strings are
            `psygrid.final_values.dtype.names`, `psygrid.history1.dtype.names`
            or `psygrid.history2.dtype.names` depending on
            "selected_star_history_for_z_var" value, and
            `psygrid.binary_history.dtype.names`.
        selected_star_history_for_z_var: int
            Accepted valuess: 1 or 2. In case z_var_str is an attribute of
            history1 or history2, then selected_star_history_for_z_var
            determines which of the two to select.
        termination_flag : str
            Termination flag to display, allowed values are:
            "termination_flag_1", "termination_flag_2", "termination_flag_3",
            "termination_flag_4", "all".
        grid_3D : bool
            If `True`, the psygrid object is a 3D grid and needs to be sliced.
        slice_3D_var_str : str
            Variable along which the 3D space will be sliced. Allowed values
            are `psygrid.initial_values.dtype.names`.
        slice_3D_var_range : tuple
            Range between which you want to slice the variable slice_3D_var_str
            e.g., `(2.5,3.)`.
        grid_4D : bool
            If `True`, the psygrid object is a 4D grid and needs to be sliced.
        slice_4D_var_str : str
            Variable along which the 4D space will be sliced. Allowed values
            are `psygrid.initial_values.dtype.names`.
        slice_4D_var_range : tuple
            Range between which you want to slice the variable slice_4D_var_str
            e.g., `(2.5,3.)`.
        extra_grid : object or array of objects
            If subset of the grid was rerun a or an extention was added, one
            can overlay the new psygrid by passing it here.
        slice_at_RLO : bool
            If `True`, the object plots the tracks until onset of Roche Lobe
            overflow.
        MARKERS_COLORS_LEGENDS : dict
            Each termination flag is associated with a marker shape, size,
            color and label (cf. `MARKERS_COLORS_LEGENDS` in
            `plot_defaults.py`).
        DEFAULT_LABELS : dict
            Each varaible is associated to an axis label. (cf. `DEFAULT_LABELS`
            in `plot_defaults.py`).
        verbose : bool
            If `True`, the object reports by printing to standard output.
        **kwargs : dict
            Dictionary containing extra visualisation options (cf.
            `PLOT_PROPERTIES` in `plot_defaults.py`.

        """
        self.psygrid = psygrid

        # info 4D/3D parameter space
        self.grid_3D = grid_3D
        self.slice_3D_var_str = slice_3D_var_str
        self.slice_3D_var_range = slice_3D_var_range
        self.grid_4D = grid_4D
        self.slice_4D_var_str = slice_4D_var_str
        self.slice_4D_var_range = slice_4D_var_range
        self.slice_at_RLO = slice_at_RLO
        self.verbose = verbose

        # store the extra psygrid
        if extra_grid is None:
            self.extra_grid = extra_grid
        if isinstance(extra_grid, list):
            raise ValueError(
                "We support only one extra psygrid at the moment!")
        else:
            self.extra_grid = extra_grid

        # read kwargs
        for key in kwargs:
            if key not in PLOT_PROPERTIES:
                raise ValueError(key + " is not a valid parameter name!")

        for varname in PLOT_PROPERTIES:
            default_value = PLOT_PROPERTIES[varname]
            if (
                varname not in ["colorbar", "legend1D", "legend2D"]
                or varname not in kwargs.keys()
            ):
                setattr(self, varname, kwargs.get(varname, default_value))
            else:
                temp_var = {}
                for sub_varname in PLOT_PROPERTIES[varname]:
                    default_value = PLOT_PROPERTIES[varname][sub_varname]
                    temp_var[sub_varname] = kwargs[varname].get(
                        sub_varname, default_value
                    )

                setattr(self, varname, temp_var)

        # plotting fonts
        plt.rcParams.update(self.rcParams)

        # store the initial/final values
        self.initial_values = self.psygrid.initial_values
        # add extra properties to initial_values
        self.add_properties_to_initial_values()
        self.initial_values_str = self.initial_values.dtype.names
        self.final_values = self.psygrid.final_values
        # add extra properties to final_values
        if termination_flag in ["combined_TF12", "debug"]:
            self.add_properties_to_final_values(termination_flag)
        if termination_flag == 'termination_flag_2':
            # remmove ? from strings
            TF2 = self.final_values['termination_flag_2']
            TF2_clean = [TF.replace('?', '') for TF in TF2]
            self.final_values['termination_flag_2'] = TF2_clean
        self.final_values_str = self.final_values.dtype.names

        # x, y and z variables must exist
        if x_var_str not in self.initial_values_str and not self.slice_at_RLO:
            raise ValueError(
                "x_var_str = {} is not available in psygrid.initial_values".
                format(x_var_str))
        elif (
            x_var_str not in self.psygrid[0].binary_history.dtype.names
            and self.slice_at_RLO
        ):
            raise ValueError("x_var_str = {} is not available in "
                             "psygrid.binary_history".format(x_var_str))
        else:
            self.x_var_str = x_var_str
        if y_var_str not in self.initial_values_str and not self.slice_at_RLO:
            raise ValueError("y_var_str = {} is not available in "
                             "psygrid.initial_values".format(y_var_str))
        elif (
            y_var_str not in self.psygrid[0].binary_history.dtype.names
            and self.slice_at_RLO
        ):
            raise ValueError("y_var_str = {} is not available in "
                             "psygrid.binary_history".format(y_var_str))
        else:
            self.y_var_str = y_var_str

        if selected_star_history_for_z_var in [1, 2]:
            self.selected_star_history_for_z_var = (
                selected_star_history_for_z_var)
        else:
            raise ValueError(
                "selected_star_history_for_z_var should be either 1 or 2"
            )
        if z_var_str is not None:
            if not self.slice_at_RLO:
                if isinstance(z_var_str, np.ndarray):
                    self.z_var_str = None
                    self.z_var = z_var_str
                    self.history = False
                    self.binary_history = False
                elif z_var_str in self.final_values_str:
                    self.z_var_str = z_var_str
                    self.history = False
                    self.binary_history = False
                elif 'relative_change' in z_var_str:
                    self.z_var_str = z_var_str
                    self.history = False
                    self.binary_history = False
                    self.add_properties_to_final_values(None)
                elif (self.selected_star_history_for_z_var == 1
                      and z_var_str in self.psygrid[0].history1.dtype.names):
                    self.z_var_str = z_var_str
                    self.history = True
                    self.binary_history = False
                elif (self.selected_star_history_for_z_var == 2
                      and z_var_str in self.psygrid[0].history2.dtype.names):
                    self.z_var_str = z_var_str
                    self.history = True
                    self.binary_history = False
                elif z_var_str in self.psygrid[0].binary_history.dtype.names:
                    self.z_var_str = z_var_str
                    self.history = False
                    self.binary_history = True
                else:
                    raise ValueError(
                        "z_var_str = {} is not available in "
                        "psygrid.final_values or psygrid.history1/2 or "
                        "psygrid.binary_history".format(z_var_str)
                    )

            else:
                if self.selected_star_history_for_z_var == 1 and \
                  z_var_str in self.psygrid[0].history1.dtype.names:
                    self.z_var_str = z_var_str
                    self.history = True
                    self.binary_history = False
                elif (self.selected_star_history_for_z_var == 2
                      and z_var_str in self.psygrid[0].history2.dtype.names):
                    self.z_var_str = z_var_str
                    self.history = True
                    self.binary_history = False
                elif z_var_str in self.psygrid[0].binary_history.dtype.names:
                    self.z_var_str = z_var_str
                    self.history = False
                    self.binary_history = True
                else:
                    raise ValueError(
                        "z_var_str = {} is not available in psygrid.history1/2"
                        " or psygrid.binary_history".format(z_var_str)
                    )

        else:
            self.z_var_str = None
            self.history = None
            self.binary_history = None

        # get values to plot
        if termination_flag in [
            "termination_flag_1",
            "termination_flag_2",
            "termination_flag_3",
            "termination_flag_4",
            "combined_TF12",
            "debug",
            "interpolation_class"
        ] or ('SN_type' in termination_flag or 
              'CO_type' in termination_flag or 
              'state' in termination_flag):
            self.all_termination_flags = False
            if 'SN_type' in termination_flag:
                self.update_markers_colors_legends('SN_type',
                                                   MARKERS_COLORS_LEGENDS)
            elif 'state' in termination_flag or 'CO_type' in termination_flag:
                self.update_markers_colors_legends('state',
                                                   MARKERS_COLORS_LEGENDS)
            else:
                self.update_markers_colors_legends(termination_flag,
                                                   MARKERS_COLORS_LEGENDS)
            self.update_values_to_plot(termination_flag)
            self.extra_grid_termination_flag = termination_flag
        elif termination_flag == "all":
            self.all_termination_flags = True
        else:
            raise ValueError('termination_flag can only be 1,2,3,4 or "all"!')

    def __call__(self):
        """Generate the plot when the class is called."""
        fig = plt.figure(figsize=self.figsize)

        if self.all_termination_flags:

            ax1 = plt.subplot(2, 2, 1)
            self.update_markers_colors_legends("termination_flag_1")
            self.update_values_to_plot("termination_flag_1")
            self.plot_panel(ax1)

            ax2 = plt.subplot(2, 2, 2)
            self.update_markers_colors_legends("termination_flag_2")
            self.update_values_to_plot("termination_flag_2")
            self.plot_panel(ax2)

            ax3 = plt.subplot(2, 2, 3)
            self.update_markers_colors_legends("termination_flag_3")
            self.update_values_to_plot("termination_flag_3")
            self.plot_panel(ax3)

            ax4 = plt.subplot(2, 2, 4)
            self.update_markers_colors_legends("termination_flag_4")
            self.update_values_to_plot("termination_flag_4")
            self.plot_panel(ax4)

            # adjust spacing
            plt.subplots_adjust(wspace=self.wspace, hspace=self.hspace)

            # save figure
            if self.fname is not None:
                fig.savefig(self.path_to_file + self.fname,
                            dpi=self.dpi, bbox_inches=self.bbox_inches)

            # show figure
            if self.show_fig:
                plt.show()

            # close figure
            if self.close_fig:
                plt.close(fig)
            else:
                return fig

        else:
            ax = plt.subplot(111)
            self.plot_panel(ax)

            # add extra layer of grid on top of the plot
            if self.extra_grid is not None:
                # switch values to extra psygrid and update values to plot
                self.initial_values = self.extra_grid.initial_values
                self.final_values = self.extra_grid.final_values
                self.add_properties_to_initial_values()
                self.update_values_to_plot(self.extra_grid_termination_flag)

                self.plot_panel(ax, extra_grid_call=True)

        # add title
        self.set_title(fig)

        # save figure
        if self.fname is not None:
            fig.savefig(self.path_to_file + self.fname,
                        dpi=self.dpi, bbox_inches=self.bbox_inches)

        # show figure
        if self.show_fig:
            plt.show()

        # close figure
        if self.close_fig:
            plt.close(fig)
        else:
            return fig

    def plot_panel(self, ax, extra_grid_call=False):
        """Plot the 2D pannel.

        Parameters
        ----------
        ax : object
            matplotlib figure axes.
        extra_grid_call : bool
            If `True`, one ore more extra grids are passed.

        """
        scatters = []
        scatters_legend = []

        # plot figure by looping over termination_flag
        sc_last = None
        for flag in self.termination_flag_str:
            selection = self.termination_flag == flag
            if self.MARKERS_COLORS_LEGENDS[flag][2] is not None:
                if self.slice_at_RLO:
                    for i in range(len(self.x_var[selection])):
                        if not isinstance(self.x_var_oRLO[selection][i],
                                          float):
                            if (not any(np.isnan(
                                self.x_var_oRLO[selection][i]))
                                    and not any(np.isnan(
                                        self.y_var_oRLO[selection][i]))):
                                plt.plot(
                                    self.x_var[selection][i],
                                    self.y_var[selection][i],
                                    marker=".",
                                    color="black",
                                )
                                plt.plot(
                                    self.x_var_oRLO[selection][i],
                                    self.y_var_oRLO[selection][i],
                                    color="black",
                                )
                            sc = ax.scatter(
                                self.x_var_oRLO[selection][i][-1],
                                self.y_var_oRLO[selection][i][-1],
                                marker=self.MARKERS_COLORS_LEGENDS[flag][0],
                                linewidths=self.MARKERS_COLORS_LEGENDS[flag][
                                    1],
                                c=self.MARKERS_COLORS_LEGENDS[flag][2],
                                s=self.marker_size,
                            )
                        else:
                            plt.plot(
                                self.x_var[selection][i],
                                self.y_var[selection][i],
                                marker=".",
                                color="black",
                            )
                            sc = ax.scatter(
                                self.x_var[selection][i],
                                self.y_var[selection][i],
                                marker=self.MARKERS_COLORS_LEGENDS[flag][0],
                                linewidths=self.MARKERS_COLORS_LEGENDS[
                                    flag][1],
                                c=self.MARKERS_COLORS_LEGENDS[flag][2],
                                s=self.marker_size,
                            )
                else:
                    sc = ax.scatter(
                        self.x_var[selection],
                        self.y_var[selection],
                        marker=self.MARKERS_COLORS_LEGENDS[flag][0],
                        linewidths=self.MARKERS_COLORS_LEGENDS[flag][1],
                        c=self.MARKERS_COLORS_LEGENDS[flag][2],
                        s=self.marker_size,
                    )
            else:
                if self.z_var is not None:
                    if self.slice_at_RLO:
                        for i in range(len(self.x_var[selection])):
                            if not isinstance(self.x_var_oRLO[selection][i],
                                              float):
                                if not any(
                                    np.isnan(self.x_var_oRLO[selection][i])
                                ) and not any(np.isnan(
                                        self.y_var_oRLO[selection][i])):
                                    plt.plot(
                                        self.x_var[selection][i],
                                        self.y_var[selection][i],
                                        marker=".",
                                        color="black",
                                    )
                                    plt.plot(
                                        self.x_var_oRLO[selection][i],
                                        self.y_var_oRLO[selection][i],
                                        color="black",
                                    )
                                sc = ax.scatter(
                                    self.x_var_oRLO[selection][i][-1],
                                    self.y_var_oRLO[selection][i][-1],
                                    marker=self.MARKERS_COLORS_LEGENDS[
                                        flag][0],
                                    linewidths=self.MARKERS_COLORS_LEGENDS[
                                        flag][1],
                                    c=self.z_var[selection][i],
                                    s=self.marker_size,
                                    alpha=0.5,
                                    vmin=self.zmin,
                                    vmax=self.zmax,
                                )
                            else:
                                plt.plot(
                                    self.x_var[selection][i],
                                    self.y_var[selection][i],
                                    marker=".",
                                    color="black",
                                )
                                sc = ax.scatter(
                                    self.x_var[selection][i],
                                    self.y_var[selection][i],
                                    marker=self.MARKERS_COLORS_LEGENDS[
                                        flag][0],
                                    linewidths=self.MARKERS_COLORS_LEGENDS[
                                        flag][1],
                                    c=self.z_var[selection][i],
                                    s=self.marker_size,
                                    alpha=0.5,
                                    vmin=self.zmin,
                                    vmax=self.zmax,
                                )

                    else:
                        sc = ax.scatter(
                            self.x_var[selection],
                            self.y_var[selection],
                            marker=self.MARKERS_COLORS_LEGENDS[flag][0],
                            linewidths=self.MARKERS_COLORS_LEGENDS[flag][1],
                            c=self.z_var[selection],
                            s=self.marker_size,
                            vmin=self.zmin,
                            vmax=self.zmax,
                        )
                    sc_last = sc
            # collect scatters for legend
            if self.MARKERS_COLORS_LEGENDS[flag][3] not in scatters_legend:
                scatters.append(sc)
                scatters_legend.append(self.MARKERS_COLORS_LEGENDS[flag][3])
        if sc_last is not None and not extra_grid_call:
            self.set_color_bar(sc_last)

        # add labels and legend
        self.set_xlabel()
        self.set_ylabel()
        self.set_xlim()
        self.set_ylim()
        self.set_legend(ax, scatters, scatters_legend)

    def add_properties_to_initial_values(self):
        """Add extra initial values."""
        # add the column mass_ratio
        old_initial_values = copy.copy(self.initial_values)
        mass_ratio = (old_initial_values["star_2_mass"]
                      / old_initial_values["star_1_mass"])
        new_initial_values = add_field(old_initial_values,
                                       [("mass_ratio", "<f8")])
        new_initial_values["mass_ratio"] = mass_ratio
        self.initial_values = new_initial_values

    def add_properties_to_final_values(self, termination_flag=None):
        """Add extra initial values."""
        old_initial_values = copy.copy(self.initial_values)
        old_final_values = copy.copy(self.final_values)
        if termination_flag == "combined_TF12":
            combined_TF12 = combine_TF12(
                old_final_values['interpolation_class'],
                old_final_values['termination_flag_2'],
                self.verbose)
            new_final_values = add_field(old_final_values,
                                         [("combined_TF12", "<U70")])
            new_final_values["combined_TF12"] = combined_TF12
        elif termination_flag == "debug":
            new_final_values = add_field(old_final_values, [("debug", "<U70")])
            new_final_values["debug"] = old_final_values['termination_flag_1']
        elif 'relative_change' in self.z_var_str:
            key = self.z_var_str.split('relative_change_')[1]
            relative_change_key = (
                (old_final_values[key] - old_initial_values[key])
                / old_initial_values[key])
            new_final_values = add_field(old_final_values,
                                         [(self.z_var_str, "<f8")])
            new_final_values[self.z_var_str] = relative_change_key
        self.final_values = new_final_values

    def update_markers_colors_legends(
        self, termination_flag, MARKERS_COLORS_LEGENDS=None
    ):
        """Udpdate markers, colors and legend.

        Parameters
        ----------
        termination_flag : string
            Termination flag to display, allowed values are:
            "termination_flag_1", "termination_flag_2", "termination_flag_3",
            "termination_flag_4".
        MARKERS_COLORS_LEGENDS : dict
            Each termination flag is associated with a marker shape, size,
            color and label.
        """
        if MARKERS_COLORS_LEGENDS is None:
            self.MARKERS_COLORS_LEGENDS = DEFAULT_MARKERS_COLORS_LEGENDS[
                termination_flag]
        else:
            self.MARKERS_COLORS_LEGENDS = MARKERS_COLORS_LEGENDS[
                termination_flag]

    def update_values_to_plot(self, termination_flag):
        """Update all values to plot.

        Parameters
        ----------
        termination_flag : string
            Termination flag to display, allowed values are:
            "termination_flag_1", "termination_flag_2", "termination_flag_3",
            "termination_flag_4".
        """
        # get termination flags
        self.termination_flag = self.final_values[termination_flag]

        # save values to plot
        self.get_x_var()
        self.get_y_var()
        if self.z_var_str is not None:
            self.get_z_var()
        elif not hasattr(self, 'z_var'):
            self.z_var = None

        # if 4D space: slice it to 2D
        if self.grid_4D and self.grid_3D:
            slice = np.logical_and(
                np.logical_and(
                    self.initial_values[self.slice_4D_var_str]
                    >= self.slice_4D_var_range[0],
                    self.initial_values[self.slice_4D_var_str]
                    <= self.slice_4D_var_range[1],
                ),
                np.logical_and(
                    self.initial_values[self.slice_3D_var_str]
                    >= self.slice_3D_var_range[0],
                    self.initial_values[self.slice_3D_var_str]
                    <= self.slice_3D_var_range[1],
                ),
            )
            self.x_var = self.x_var[slice]
            self.y_var = self.y_var[slice]
            if self.z_var is not None:
                self.z_var = self.z_var[slice]
            if self.slice_at_RLO:
                self.x_var_oRLO = self.x_var_oRLO[slice]
                self.y_var_oRLO = self.y_var_oRLO[slice]

            self.termination_flag = self.termination_flag[slice]
            if self.verbose:
                print("The 4D space was sliced along {} = {} and {} = {}.".
                      format(self.slice_4D_var_str, self.slice_4D_var_range,
                             self.slice_3D_var_str, self.slice_3D_var_range))
                print("")
                print("Total values to plot {}".
                      format(len(self.termination_flag)))

        # if 3D space: slice it to 2D
        if (not self.grid_4D) and self.grid_3D:
            slice = np.logical_and(
                self.initial_values[self.slice_3D_var_str]
                >= self.slice_3D_var_range[0],
                self.initial_values[self.slice_3D_var_str]
                <= self.slice_3D_var_range[1],
            )
            self.x_var = self.x_var[slice]
            self.y_var = self.y_var[slice]
            if self.z_var is not None:
                self.z_var = self.z_var[slice]
            if self.slice_at_RLO:
                self.x_var_oRLO = self.x_var_oRLO[slice]
                self.y_var_oRLO = self.y_var_oRLO[slice]
            self.termination_flag = self.termination_flag[slice]
            if self.verbose:
                print(
                    "The 3D space was sliced along {} = {}.".format(
                        self.slice_3D_var_str, self.slice_3D_var_range
                    )
                )
                print("")
                print("Total values to plot {}.".
                      format(len(self.termination_flag)))

        # find all different termination flags
        self.termination_flag_str = np.unique(self.termination_flag)

        # fix max min color bar
        if self.z_var is not None:
            if self.zmin is None:
                not_nan = np.invert(np.isnan(self.z_var))
                self.zmin = min(self.z_var[not_nan])

            if self.zmax is None:
                not_nan = np.invert(np.isnan(self.z_var))
                self.zmax = max(self.z_var[not_nan])

    def get_x_var(self):
        """Get x value to plot."""
        if self.log10_x:
            self.x_var = np.log10(self.initial_values[self.x_var_str])
        else:
            self.x_var = self.initial_values[self.x_var_str]

        if self.slice_at_RLO:
            values = []
            for run in self.psygrid:
                # failed runs are stored as signle values or None, not arrays
                if run.binary_history is None:
                    values.append(np.nan)
                elif isinstance(run.binary_history[self.x_var_str],
                                np.ndarray):
                    # index of onset of RLO
                    indicies_RLO = np.where(
                        (run.binary_history["rl_relative_overflow_1"] >= -0.05)
                        & (run.binary_history["lg_mtransfer_rate"] >= -12)
                    )[0]
                    if len(indicies_RLO) > 0:
                        index_oRLO = indicies_RLO[0]
                        if self.log10_x:
                            values.append(np.log10(run.binary_history[
                                self.x_var_str][: index_oRLO + 1]))
                        else:
                            values.append(run.binary_history[
                                self.x_var_str][: index_oRLO + 1])
                    else:
                        # the array is empty, no RLO
                        values.append(np.nan)
                else:
                    values.append(np.nan)
            self.x_var_oRLO = np.array(values, dtype=object)

    def get_y_var(self):
        """Get y value to plot."""
        if self.log10_y:
            self.y_var = np.log10(self.initial_values[self.y_var_str])
        else:
            self.y_var = self.initial_values[self.y_var_str]

        if self.slice_at_RLO:
            values = []
            for run in self.psygrid:
                # failed runs are stored as signle values or None, not arrays
                if run.binary_history is None:
                    values.append(np.nan)
                elif isinstance(run.binary_history[self.y_var_str],
                                np.ndarray):
                    # index of onset of RLO
                    indicies_RLO = np.where(
                        (run.binary_history["rl_relative_overflow_1"] >= -0.05)
                        & (run.binary_history["lg_mtransfer_rate"] >= -12)
                    )[0]
                    if len(indicies_RLO) > 0:
                        index_oRLO = indicies_RLO[0]
                        if self.log10_y:
                            values.append(np.log10(run.binary_history[
                                self.y_var_str][: index_oRLO + 1]))
                        else:
                            values.append(run.binary_history[
                                self.y_var_str][: index_oRLO + 1])
                    else:
                        # the array is empty, no RLO
                        values.append(np.nan)
                else:
                    values.append(np.nan)
            self.y_var_oRLO = np.array(values, dtype=object)

    def get_z_var(self):
        """Get z value to plot."""
        if self.history is None:
            raise ValueError("Something went wrong!")

        # read the final values from history1
        if not self.slice_at_RLO:
            if self.history:
                final_values = []
                for run in self.psygrid:
                    if self.selected_star_history_for_z_var == 1:
                        history = run.history1
                    elif self.selected_star_history_for_z_var == 2:
                        history = run.history2
                    else:
                        raise ValueError(
                            "wrong selected_star_history_for_z_var")
                    # failed runs are stored as signle values or None
                    # and not arrays
                    if history is None:
                        final_values.append(np.nan)
                    elif isinstance(history[self.z_var_str], np.ndarray):
                        final_values.append(history[self.z_var_str][-1])
                    else:
                        final_values.append(np.nan)
                final_values = np.array(final_values)
                if self.log10_z:
                    self.z_var = np.log10(final_values)
                else:
                    self.z_var = final_values

            elif self.binary_history:
                final_values = []
                for run in self.psygrid:
                    # failed runs are stored as signle values or None
                    # and not arrays
                    if run.binary_history is None:
                        final_values.append(np.nan)
                    elif isinstance(run.binary_history[self.z_var_str],
                                    np.ndarray):
                        final_values.append(
                            run.binary_history[self.z_var_str][-1])
                    else:
                        final_values.append(np.nan)
                final_values = np.array(final_values)
                if self.log10_z:
                    self.z_var = np.log10(final_values)
                else:
                    self.z_var = final_values

            # read final values from final_values
            else:
                if self.log10_z:
                    self.z_var = np.log10(self.final_values[self.z_var_str])
                else:
                    self.z_var = self.final_values[self.z_var_str]

        # take the z_var at oRLO
        else:
            if self.history is None and self.binary_history is None:
                raise ValueError("Something went wrong!")

            if self.history:
                values = []
                for run in self.psygrid:
                    if self.selected_star_history_for_z_var == 1:
                        history = run.history1
                    elif self.selected_star_history_for_z_var == 2:
                        history = run.history2
                    else:
                        raise ValueError(
                            "wrong selected_star_history_for_z_var")
                    # failed runs are stored as signle values or None
                    # and not arrays
                    if history is None:
                        values.append(np.nan)
                    elif isinstance(history[self.z_var_str], np.ndarray):
                        # index of onset of RLO
                        indicies_RLO = np.where(
                            (run.binary_history["rl_relative_overflow_1"]
                             >= -0.05)
                            & (run.binary_history["lg_mtransfer_rate"] >= -12)
                        )[0]
                        if len(indicies_RLO) > 0:
                            index_oRLO = indicies_RLO[0]
                            values.append(history[self.z_var_str][index_oRLO])
                        else:
                            # the array is empty, no RLO
                            values.append(np.nan)
                    else:
                        values.append(np.nan)
                values = np.array(values)
                if self.log10_z:
                    self.z_var = np.log10(values)
                else:
                    self.z_var = values

            elif self.binary_history:
                values = []
                for run in self.psygrid:
                    # failed runs are stored as signle values or None
                    # and not arrays
                    if run.binary_history is None:
                        values.append(np.nan)
                    elif isinstance(run.binary_history[self.z_var_str],
                                    np.ndarray):
                        # index of onset of RLO
                        indicies_RLO = np.where(
                            (run.binary_history["rl_relative_overflow_1"]
                             >= -0.05)
                            & (run.binary_history["lg_mtransfer_rate"] >= -12)
                        )[0]
                        if len(indicies_RLO) > 0:
                            index_oRLO = indicies_RLO[0]
                            values.append(
                                run.binary_history[self.z_var_str][index_oRLO]
                            )
                        else:
                            # the array is empty, no RLO
                            values.append(np.nan)
                    else:
                        values.append(np.nan)
                values = np.array(values)
                if self.log10_z:
                    self.z_var = np.log10(values)
                else:
                    self.z_var = values

    def set_title(self, fig):
        """Add title.

        Parameters
        ----------
        fig : object
            matplotlib figure object.
        """
        if self.title is not None and not self.all_termination_flags:
            plt.title(self.title,
                      fontdict=self.title_font_dict, loc=self.title_loc)
        elif self.title is not None and self.all_termination_flags:
            fig.suptitle(self.title, fontdict=self.title_font_dict)

    def set_xlabel(self):
        """Add x label."""
        if self.xlabel is not None:
            plt.xlabel(self.xlabel, **self.xlabel_kwargs)
        else:
            if self.log10_x:
                plt.xlabel(DEFAULT_LABELS[self.x_var_str][1],
                           **self.xlabel_kwargs)
            else:
                plt.xlabel(DEFAULT_LABELS[self.x_var_str][0],
                           **self.xlabel_kwargs)

    def set_ylabel(self):
        """Add y label."""
        if self.ylabel is not None:
            plt.ylabel(self.ylabel, **self.ylabel_kwargs)
        else:
            if self.log10_y:
                plt.ylabel(DEFAULT_LABELS[self.y_var_str][1],
                           **self.ylabel_kwargs)
            else:
                plt.ylabel(DEFAULT_LABELS[self.y_var_str][0],
                           **self.ylabel_kwargs)

    def set_xlim(self):
        """Set x axes limits."""
        if self.xmin is not None and self.xmax is not None:
            plt.xlim(self.xmin, self.xmax)

    def set_ylim(self):
        """Set y axes limits."""
        if self.ymin is not None and self.ymax is not None:
            plt.ylim(self.ymin, self.ymax)

    def set_legend(self, ax, scatters, scatters_legend):
        """Add legend.

        Parameters
        ----------
        ax : object
            matplotlib figure axes.
        scatters : object
            matplotlib scatter object.
        scatters_legend : list of str
            List of strings which will be used as labels.

        """
        if self.legend2D["title"] is not None:

            # defailt: shrink current axis by 20% and put tje legend to
            # the right of the current axis
            box = ax.get_position()
            ax.set_position(
                [box.x0, box.y0, box.width * self.legend2D["shrink_box"],
                 box.height])

            ax.legend(
                scatters,
                scatters_legend,
                borderaxespad=self.legend2D["borderaxespad"],
                handletextpad=self.legend2D["handletextpad"],
                columnspacing=self.legend2D["columnspacing"],
                title=self.legend2D["title"],
                title_fontsize=self.legend2D["title_font_size"],
                prop=self.legend2D["prop"],
                loc=self.legend2D["loc"],
                ncol=self.legend2D["ncol"],
                bbox_to_anchor=self.legend2D["bbox_to_anchor"],
            )

    def set_color_bar(self, scatter):
        """Add colorbar.

        Parameters
        ----------
        scatters : object
            matplotlib scatter object.
        """
        if self.colorbar["label"] is not None:
            label = self.colorbar["label"]
        elif isinstance(self.z_var_str, str):
            z_var_str = self.z_var_str.replace('S1_', '').replace('S2_', '')
            if z_var_str in DEFAULT_LABELS.keys():
                if self.log10_z:
                    label = DEFAULT_LABELS[z_var_str][1]
                else:
                    label = DEFAULT_LABELS[z_var_str][0]
            else:
                label = None
        else:
            label = None

        cbar = plt.colorbar(
            mappable=scatter,
            orientation=self.colorbar["orientation"],
            fraction=self.colorbar["fraction"],
            pad=self.colorbar["pad"],
            shrink=self.colorbar["shrink"],
            aspect=self.colorbar["aspect"],
            anchor=self.colorbar["anchor"],
            panchor=self.colorbar["panchor"],
            extend=self.colorbar["extend"],
        )
        
        cbar.set_label(label=label, size=self.colorbar["label_size"])
        cbar.mappable.set_clim(self.colorbar["vmin"], vmax = self.colorbar["vmax"])
        
