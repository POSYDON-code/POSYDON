"""Plotting class for 1D (MESA) psygrids.

The 1D visualization plotting class allows to plot 1D tracks of PsyGrid
objects. The PsyGrid object is composed of nD MESA grid run with POSYDON
and post processed with the psygrid object into an h5 file.
"""


__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]


import matplotlib.pyplot as plt
import numpy as np
from labellines import labelLine
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

import posydon.utils.constants as const
from posydon.visualization.plot_defaults import (
    DEFAULT_LABELS,
    DEFAULT_MARKERS_COLORS_LEGENDS,
    PLOT_PROPERTIES,
)


class plot1D(object):
    """Plotting class for 1D (MESA) grids."""

    def __init__(
        self,
        run,
        x_var_str,
        y_var_str,
        z_var_str=None,
        history="binary_history",
        star_states=None,
        HR=False,
        verbose=False,
        **kwargs
    ):
        """Read a PsyGrid object and plot a 1D track of x vs y.

        Parameters
        ----------
        run : int or list of int
            Index or list of indeces of the PsyGrid object you would like to
            plot.
        x_var_str : str
            String of values to plot on the x axis. Allowed strings are the
            one in `psygrid.history.dtype.names` where "history" needs to be
            chosen accordingly.
        y_var_str : str or list of str
            String or list of stringvalues to plot on the y axis. Allowed
            strings are the one in `psygrid.history.dtype.names` where
            "history" needs to be chosen accordingly.
        z_var_str : str
            String of values to plot on the z axis (displayed with a color).
            Allowed strings are the one in `psygrid.history.dtype.names` where
            "history" needs to be chosen accordingly.
        history : str
            The x, y, z variables are read from either: "binary_history",
            "history1", "history2".
        HR : bool
            If `True`, an HR diagram will be plotted.
        verbose : bool
            If `True`, the object reports by printing to standard output.
        **kwargs : dict
            Dictionary containing extra visualisation options (cf.
            `PLOT_PROPERTIES` in `plot_defaults.py`.

        """
        self.verbose = verbose
        self.HR = HR
        self.run = run
        self.star_states = star_states

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

        # history
        if history not in ["binary_history", "history1", "history2"]:
            raise ValueError(
                "history is must be either binary_history or history1/2!")

        # store the history values
        if isinstance(self.run, list):
            self.history = []
            for run in self.run:
                self.history.append(run[history])
        else:
            self.history = [self.run[history]]
        self.n_runs = len(self.history)
        self.history_str = ()
        for i in range(self.n_runs):
            if self.history[i] is not None:
                self.history_str = self.history[i].dtype.names
                break

        # x, y and z variables must exist
        if isinstance(x_var_str, str):
            if x_var_str not in self.history_str:
                raise ValueError(
                    "x_var_str = {} is not available in run.{}".format(
                        x_var_str, self.history
                    )
                )
            else:
                self.x_var_str = x_var_str
        elif x_var_str is not None:
            raise ValueError(
                "x_var_str = {} is not available a str type!".format(x_var_str)
            )
        if isinstance(y_var_str, str):
            if y_var_str not in self.history_str:
                raise ValueError(
                    "y_var_str = {} is not available in run.{}".format(
                        y_var_str, self.history
                    )
                )
            else:
                self.y_var_str = y_var_str
                self.number_of_plots = 1
        elif isinstance(y_var_str, list):
            for var_str in y_var_str:
                if var_str not in self.history_str:
                    raise ValueError(
                        "y_var_str = {} is not available in run.{}".
                        format(var_str, self.history))
            self.y_var_str = y_var_str
            self.number_of_plots = len(y_var_str)
        elif x_var_str is not None:
            raise ValueError(
                "y_var_str = {} is not available a str type!".format(y_var_str)
            )
        else:
            self.number_of_plots = 1
        if z_var_str is None:
            self.z_var_str = None
        elif isinstance(z_var_str, str):
            if z_var_str not in self.history_str:
                raise ValueError(
                    "z_var_str = {} is not available in run.{}".format(
                        z_var_str, self.history
                    )
                )
            else:
                self.z_var_str = z_var_str
                self.number_of_plots = 1
        else:
            raise ValueError(
                "z_var_str = {} is not str type!".format(z_var_str))

    def __call__(self):
        """Generate the plot when the class is called."""
        if not self.HR:
            fig = plt.figure(figsize=self.figsize)

            if self.number_of_plots > 1:
                i = 1
                for var_str in self.y_var_str:
                    ax = plt.subplot(self.number_of_plots, 1, i)
                    self.update_values_to_plot(var_str)
                    self.plot_panel(ax, i, var_str)
                    i += 1

                # adjust spacing
                plt.subplots_adjust(wspace=self.wspace, hspace=self.hspace)
            elif self.number_of_plots == 1:
                if self.z_var_str is None:
                    ax = plt.subplot(111)
                    self.update_values_to_plot(self.y_var_str)
                    self.plot_panel(ax, 1, self.y_var_str)
                else:
                    ax = plt.subplot(111)
                    self.update_values_to_plot(self.y_var_str)
                    self.plot_panel(ax, 1, self.y_var_str)
            else:
                raise ValueError(
                    "number_of_plots={} must be a positive integer!".format(
                        self.number_of_plots
                    )
                )

            # add title
            self.set_title(fig)

            # save figure
            if self.PdfPages is not None:
                self.PdfPages.savefig(figure=fig, dpi=self.dpi,
                                      bbox_inches=self.bbox_inches)
            elif self.fname is not None:
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
            self.HR_diagram()

    def plot_panel(self, ax, i, y_var_str):
        """Plot the 1D pannel.

        Parameters
        ----------
        ax : object
            matplotlib figure axes.
        i : int
            Index of the run to plot.
        y_var_str : str
            String or list of stringvalues to plot on the y axis.

        """
        if self.z_var_str is None:

            lines = []
            for j in range(self.n_runs):
                (line,) = ax.plot(self.x_var[j], self.y_var[j])
                lines.append(line)

            # add labels and legend
            if i == self.number_of_plots:
                self.set_xlabel()
                self.set_legend(ax, lines)
            else:
                ax.set_xticklabels([])

            self.set_ylabel(i)
            self.set_xlim()
            self.set_ylim()

        else:

            lines = []
            for j in range(self.n_runs):
                (line,) = ax.plot(self.x_var[j], self.y_var[j], zorder=1)
                sc = ax.scatter(
                    self.x_var[j],
                    self.y_var[j],
                    c=self.z_var[j],
                    s=self.marker_size,
                    vmin=self.zmin,
                    vmax=self.zmax,
                    zorder=2,
                )
                lines.append(line)

            # add labels and legend
            if i == self.number_of_plots:
                self.set_xlabel()
                self.set_legend(ax, lines)
                self.set_color_bar(sc)
            else:
                ax.set_xticklabels([])

            self.set_ylabel(i)
            self.set_xlim()
            self.set_ylim()

    def update_values_to_plot(self, y_var_str):
        """Update all values to plot.

        Parameters
        ----------
        y_var_str : str
            String or list of stringvalues to plot on the y axis.

        """
        # save values to plot
        if self.log10_x:
            self.x_var = []
            for history in self.history:
                self.x_var.append(np.log10(history[self.x_var_str]))
        else:
            self.x_var = []
            for history in self.history:
                self.x_var.append(history[self.x_var_str])
        if self.log10_y:
            self.y_var = []
            for history in self.history:
                self.y_var.append(np.log10(history[y_var_str]))
        else:
            self.y_var = []
            for history in self.history:
                self.y_var.append(history[y_var_str])

        if self.z_var_str is not None:
            if self.log10_z:
                self.z_var = []
                for history in self.history:
                    self.z_var.append(np.log10(history[self.z_var_str]))
            else:
                self.z_var = []
                for history in self.history:
                    self.z_var.append(history[self.z_var_str])

            # ensure to have the minimal value for many tracks scatter plots
            if self.zmin is None:
                self.zmin = min(self.z_var[0])
                for z_var in self.z_var:
                    if self.zmin > min(z_var):
                        self.zmin = min(z_var)
            if self.zmax is None:
                self.zmax = max(self.z_var[0])
                for z_var in self.z_var:
                    if self.zmax < max(z_var):
                        self.zmax = max(z_var)

    def set_title(self, fig):
        """Add title.

        Parameters
        ----------
        fig : object
            matplotlib figure object.
        """
        if self.title is not None and self.number_of_plots == 1:
            plt.title(self.title, fontdict=self.title_font_dict,
                      loc=self.title_loc)
        elif self.title is not None and self.number_of_plots > 1:
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

    def set_ylabel(self, i):
        """Add y label."""
        if self.ylabel is not None:
            if isinstance(self.ylabel, str):
                plt.ylabel(self.ylabel, **self.ylabel_kwargs)
            elif isinstance(self.ylabel, list):
                plt.ylabel(self.ylabel[i - 1], **self.ylabel_kwargs)
        else:
            if self.log10_y:
                if isinstance(self.y_var_str, str):
                    plt.ylabel(DEFAULT_LABELS[self.y_var_str][1],
                               **self.ylabel_kwargs)
                elif isinstance(self.y_var_str, list):
                    plt.ylabel(
                        DEFAULT_LABELS[self.y_var_str[i - 1]][1],
                        **self.ylabel_kwargs
                    )
            else:
                if isinstance(self.y_var_str, str):
                    plt.ylabel(DEFAULT_LABELS[self.y_var_str][0],
                               **self.ylabel_kwargs)
                elif isinstance(self.y_var_str, list):
                    plt.ylabel(
                        DEFAULT_LABELS[self.y_var_str[i - 1]][0],
                        **self.ylabel_kwargs
                    )

    def set_xlim(self):
        """Set x axes limits."""
        if self.xmin is not None and self.xmax is not None:
            plt.xlim(self.xmin, self.xmax)

    def set_ylim(self):
        """Set y axes limits."""
        if self.ymin is not None and self.ymax is not None:
            plt.ylim(self.ymin, self.ymax)

    def set_legend(self, ax, lines):
        """Add legend.

        Parameters
        ----------
        ax : object
            matplotlib figure axes.
        lines : object
            matplotlib lines object.
        """
        if self.legend1D["lines_legend"] is not None:

            # defailt: shrink current axis by 20% and put tje legend to
            # the right of the current axis
            box = ax.get_position()
            ax.set_position(
                [box.x0, box.y0, box.width * self.legend1D["shrink_box"],
                 box.height])

            ax.legend(
                lines,
                self.legend1D["lines_legend"],
                borderaxespad=self.legend1D["borderaxespad"],
                handletextpad=self.legend1D["handletextpad"],
                columnspacing=self.legend1D["columnspacing"],
                title=self.legend1D["title"],
                title_fontsize=self.legend1D["title_font_size"],
                prop=self.legend1D["prop"],
                loc=self.legend1D["loc"],
                ncol=self.legend1D["ncol"],
                bbox_to_anchor=self.legend1D["bbox_to_anchor"],
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

        plt.colorbar(
            mappable=scatter,
            orientation=self.colorbar["orientation"],
            fraction=self.colorbar["fraction"],
            pad=self.colorbar["pad"],
            shrink=self.colorbar["shrink"],
            aspect=self.colorbar["aspect"],
            anchor=self.colorbar["anchor"],
            panchor=self.colorbar["panchor"],
            extend=self.colorbar["extend"],
        ).set_label(label=label, size=self.colorbar["label_size"])

    def lines_constant_radius(self):
        """Constant radius lines for the HR diagram."""

        def luminosity(Teff, R):
            return (
                4
                * np.pi
                * const.boltz_sigma
                * Teff ** 4
                * (R * const.Rsun) ** 2
                / const.Lsun
            )

        # constan radius line
        Teff = np.logspace(3., 6., 100)
        return (
            [r"$0.001\,R_\odot$", r"$0.01\,R_\odot$", r"$0.1\,R_\odot$",
             r"$1\,R_\odot$", r"$10\,R_\odot$", r"$100\,R_\odot$",
             r"$1000\,R_\odot$"],
            Teff,
            [
                luminosity(Teff, R=0.001),
                luminosity(Teff, R=0.01),
                luminosity(Teff, R=0.1),
                luminosity(Teff, R=1),
                luminosity(Teff, R=10),
                luminosity(Teff, R=100),
                luminosity(Teff, R=1000),
            ],
        )

    def HR_diagram(self):
        """Plot and HR diagram."""
        if "log_Teff" not in self.history_str:
            raise ValueError(
                "You cannot plot an HR diagram without log_Teff in history!"
            )

        if "log_L" not in self.history_str:
            raise ValueError(
                "You cannot plot an HR diagram without log_L in history!")

        fig = plt.figure(figsize=self.figsize)
        ax = plt.subplot(111)

        # flip axes
        plt.gca().invert_xaxis()

        if self.const_R_lines:
            handels, Teff, luminosities = self.lines_constant_radius()
            for i, L in enumerate(luminosities):
                slice_T = np.logical_and(
                    np.log10(Teff) > self.xmin * 0.98,
                    np.log10(Teff) < self.xmax * 1.02)
                slice_L = np.logical_and(
                    np.log10(L) > self.ymin * 0.98,
                    np.log10(L) < self.ymax * 1.02)
                slice = np.logical_and(slice_T, slice_L)
                if len(Teff[slice]) > 0:
                    (line,) = plt.plot(
                        np.log10(Teff[slice]),
                        np.log10(L[slice]),
                        "-.",
                        color="gray",
                        linewidth=0.5,
                        zorder=1.,
                    )

                    if self.xmin < min(np.log10(Teff[slice])):
                        x = np.log10(Teff[slice])[0] * 1.05
                    else:
                        x = self.xmin * 1.05
                    labelLine(line, x, label=handels[i],
                              align=True, fontsize=5, zorder=1.5)

        lines = []
        for j in range(self.n_runs):
            if self.star_states is None:
                (line,) = ax.plot(
                    self.history[j]["log_Teff"],
                    self.history[j]["log_L"],
                )
                lines.append(line)
            else:
                points = np.array(
                    [self.history[j]["log_Teff"],
                     self.history[j]["log_L"]]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                convention = DEFAULT_MARKERS_COLORS_LEGENDS[
                    'termination_flag_4']
                states_values = [convention[key][2]
                                 for key in self.star_states[j]]
                lc = LineCollection(segments,
                                    colors=states_values, linewidth=1)
                ax.add_collection(lc)
                # plot with a marker the endpoint of the evolion
                # the code does not work plus the majority of stars do not have
                # the points in the EEPs
                # end = np.logical_or(self.star_states[j] ==
                #                     'H-rich_Core_C_depleted',
                #                     self.star_states[j]
                #                     == 'stripped_He_Core_C_depleted')
                # end_x = self.history[j]["log_Teff"][end]
                # end_y = self.history[j]["log_L"][end]
                # ax.plot(end_x, end_y, marker='o', markersize=10,
                #         color=convention['H-rich_Core_C_depleted'][2])
                if j == 0:
                    custom_lines = []
                    custom_legend = []
                    key_skip = ['undetermined_evolutionary_state',
                                'BH', 'NS',
                                'ignored_no_binary_history', 'ignored_no_RLO',
                                'H-rich_non_burning',
                                'stripped_He_non_burning',
                                'accreted_He_Core_H_burning']
                    for key in convention.keys():
                        if key in key_skip:
                            continue
                        custom_lines.append(Line2D([0], [0],
                                                   color=convention[key][2]))
                        edited_key = key.replace('_', ' ')
                        edited_key = edited_key.replace('Core', 'core')
                        edited_key = edited_key.replace('Shell', 'shell')
                        edited_key = edited_key.replace('Central', 'central')
                        custom_legend.append(edited_key)
                    ax.legend(custom_lines, custom_legend,
                              borderaxespad=self.legend1D["borderaxespad"],
                              handletextpad=self.legend1D["handletextpad"],
                              columnspacing=self.legend1D["columnspacing"],
                              title=self.legend1D["title"],
                              title_fontsize=self.legend1D["title_font_size"],
                              prop=self.legend1D["prop"],
                              loc=self.legend1D["loc"],
                              ncol=self.legend1D["ncol"],
                              bbox_to_anchor=self.legend1D["bbox_to_anchor"])

            if (self.history[j] is not None and
                "star_mass" in self.history[j].dtype.names):
                plt.text(
                    self.history[j]["log_Teff"][0] * 1.1,
                    self.history[j]["log_L"][0],
                    r'$%3.1f \, M_\odot$' % self.history[j]["star_mass"][0],
                    fontsize=5
                )

            if self.xmin is None:
                self.xmin = min(self.history[j]["log_Teff"]) * 0.98
            elif min(self.history[j]["log_Teff"]) < self.xmin:
                self.xmin = min(self.history[j]["log_Teff"]) * 0.98
            if self.xmax is None:
                self.xmax = max(self.history[j]["log_Teff"]) * 1.02
            elif max(self.history[j]["log_Teff"]) > self.xmax:
                self.xmax = max(self.history[j]["log_Teff"]) * 1.02

            if self.ymin is None:
                self.ymin = min(self.history[j]["log_L"]) * 0.98
            elif min(self.history[j]["log_L"]) < self.ymin:
                self.ymin = min(self.history[j]["log_L"]) * 0.98
            if self.ymax is None:
                self.ymax = max(self.history[j]["log_L"]) * 1.02
            elif max(self.history[j]["log_L"]) > self.ymax:
                self.ymax = max(self.history[j]["log_L"]) * 1.02

        self.xlabel = r"$\log_{10}(T_\mathrm{eff}/K)$"
        self.ylabel = r"$\log_{10}(L/L_\odot)$"

        self.set_title(fig)
        self.set_xlabel()
        self.set_ylabel(1)
        self.set_xlim()
        plt.gca().invert_xaxis()
        self.set_ylim()
        self.set_legend(ax, lines)
        # save figure
        if self.PdfPages is not None:
            self.PdfPages.savefig(figure=fig, dpi=self.dpi,
                                  bbox_inches=self.bbox_inches)
        elif self.fname is not None:
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
