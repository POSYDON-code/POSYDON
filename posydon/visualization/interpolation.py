""" Module to evaluate IFInterpolator class """

__authors__ = [
    "Philipp Moura Srivastava <philipp.msrivastava@gmail.com>"
]


from posydon.visualization.plot2D import plot2D
import matplotlib.pyplot as plt
import numpy as np

class EvaluateIFInterpolator:
    """ Class that is helpful for evaluating interpolation performance
    """
    def __init__(self, interpolator, test_grid):
        """ Initialize the EvaluateIFInterpolator class

        Parameters
        ----------
        interpolator : IFInterpolator
            Interpolator that user wants to test
        test_grid : PSyGrid
            Grid object containing testing tracks

        """

        self.interpolator = interpolator

        self.test_grid = test_grid

        # assuming that in_keys are the same for all interpolators
        self.in_keys = interpolator.interpolators[0].in_keys
        self.out_keys = []

        for interp in interpolator.interpolators: # getting out keys
            self.out_keys.extend(interp.out_keys)

        self.__compute_errs()

    def __compute_errs(self):
        """ Method that computes both interpolation and classification errors """

        iv = np.array(self.test_grid.initial_values[self.in_keys].tolist()) # initial values
        fv = np.array(self.test_grid.final_values[self.out_keys].tolist()) # final values
        ic = self.test_grid.final_values["interpolation_class"] # final values

        ivalid_inds = np.where(
            (self.test_grid.final_values["interpolation_class"] != "not_converged") &
            (self.test_grid.final_values["interpolation_class"] != "initial_MT")
        )

        i = self.interpolator.test_interpolator(iv) # interpolated
        # ifv = fv[ivalid_inds]

        self.errs = {}

        self.errs["relative"] = np.abs((fv - i) / fv)
        self.errs["absolute"] = np.abs(fv - i)
        self.errs["valid_inds"] = ivalid_inds
        

        cvalid_inds = np.where(
            self.test_grid.final_values["interpolation_class"] != "not_converged"
        )[0]

        self.cfv = self.test_grid.final_values[cvalid_inds]

        c = self.interpolator.test_classifiers(iv[cvalid_inds]) # classifying

        # computing confusion matrices
        self.matrices = {}

        for key, value in c.items():

            labels = self.cfv[key]

            classes = self.__find_labels(key)

            matrix = {}

            for _class in classes:
                class_inds = np.where(labels == _class)[0]
                pred_classes, counts = np.unique(value[class_inds], return_counts = True)

                row = {c: 0 for c in classes}

                for pred_class, count in zip(pred_classes, counts):
                    row[pred_class] = count / len(class_inds)
            
                matrix[_class] = row

            self.matrices[key] = matrix

        # saving classes
        self.c = c
        self.ic = ic
                

    def __format(self, s, title = False):
        """ Method that formats keys for plots 

        Parameters
        ----------
        s : str
            string to be formatted

        """

        return s.replace("_", " ").title()

    def __find_labels(self, key):
        """ Method that finds labels in classifier

        Parameters
        ----------
        key : str
            name of the classifier
            
        Returns
        -------
        list of class labels

        """

        labels = None

        for interp in self.interpolator.interpolators:
            
            if key in interp.classifiers.keys():
                labels = interp.classifiers[key].labels # finding labels and saving

        return labels

    def __clean_errs(self, errs):

        errs = errs[~np.isnan(errs).any(axis = 1)] # dropping nans
        errs = errs[~np.isinf(errs).any(axis = 1)] # dropping infs

        return errs


    def violin_plots(self, err_type = "relative", keys = None, save_path = None):
        """ Method that plots distribution of specified error for given keys and
        optionally saves it.

        Parameters
        ----------

        err_type : str
            Either relative or absolute, default is relative
        keys : list
            A list of keys for which the errors will be shown, by default is all of them
        save_path: str
            The path where the figure should be saved to

        """

        if keys is None:
            keys = self.out_keys

        k_inds = [self.out_keys.index(key) for key in keys]
        
        errs = self.__clean_errs(self.errs[err_type].T[k_inds].T)

        n_tracks = self.test_grid.final_values["star_1_mass"].shape[0]

        print(f"Using {errs.shape[0]} out of {n_tracks} tracks in violin plots")

        stable_inds = np.where(self.ic == "stable_MT")
        no_inds = np.where(self.ic == "no_MT")
        unstable_inds = np.where(self.ic == "unstable_MT")

        stable_errs = self.__clean_errs(self.errs[err_type].T[k_inds].T[stable_inds])
        no_errs = self.__clean_errs(self.errs[err_type].T[k_inds].T[no_inds])
        unstable_errs = self.__clean_errs(self.errs[err_type].T[k_inds].T[unstable_inds])


        plt.rcParams.update({"font.size": 32, "font.family": "stixgeneral"})

        fig, axs = plt.subplots(1, 1,
                                figsize = (24, 10),
                                tight_layout = True)
        
        rel_plot = axs.violinplot(np.log10(errs + 1.0e-8), showmedians = True, points = 1000)
        stable_plot = axs.violinplot(np.log10(stable_errs + 1.0e-8), showmedians = True, points = 1000)
        no_plot = axs.violinplot(np.log10(no_errs + 1.0e-8), showmedians = True, points = 1000)
        unstable_plot = axs.violinplot(np.log10(unstable_errs + 1.0e-8), showmedians = True, points = 1000)

        axs.set_title(f"Distribution of {err_type.capitalize()} Errors")
        axs.set_xticks(np.arange(1, len(keys) + 1), 
            labels = [
                f"{self.__format(ec)} ({(med * 100):.2f}%)" for ec, med in zip(keys, np.nanmedian(errs, axis = 0))
            ], rotation = 20)

        axs.set_ylabel("Errors in Log 10 Scale")
        axs.grid(axis = "y")

        def halve_paths(field, color, right = True):
            
            for i, path in enumerate(field.get_paths()):

                # getting mean
                m = np.mean(path.vertices[:, 0])

                first = m if right == True else -np.inf
                second = np.inf if right == True else m

                # modify the paths to not go further left than the center
                field.get_paths()[i].vertices[:, 0] = np.clip(path.vertices[:, 0], first, second)
                field.set_edgecolor(color)

        def customize_violinplot(plot, color, outlined = False, right = True):

            halve_paths(plot["cmins"], color, right = right)
            halve_paths(plot["cmaxes"], color, right = right)
            halve_paths(plot["cmedians"], color, right = right)

            for pc in plot["bodies"]:

                halve_paths(pc, color, right = right)


                pc.set_facecolor(color if not outlined else "None")
                pc.set_edgecolor(color)
                pc.set_linewidth(4)
                pc.set_alpha(0.75)
        
        customize_violinplot(rel_plot, "coral")
        customize_violinplot(stable_plot, "#1e90ff", True, False)
        customize_violinplot(no_plot, "crimson", True, False)
        customize_violinplot(unstable_plot, "olive", True, False)

        axs.legend(
            [rel_plot["bodies"][0], stable_plot["bodies"][0], no_plot["bodies"][0], unstable_plot["bodies"][0]], 
            ["Relative Error", "Stable MT Error", "No MT Error", "Unstable MT Error"],
            bbox_to_anchor = (0, 1.02, 1, 0.2),
            loc = "lower left",
            mode = "expand",
            ncol = 4
        )

        plt.show()

        if save_path is not None:
            fig.save(save_path)


    def confusion_matrix(self, key, params = {}, save_path = None):
        """ Method that plots confusion matrices to evaluate classification

        Parameters
        ----------
        key : str
            The key for the classifier of interest
        params : dict
            Extra params to pass to matplolib, x_labels (list), y_labels (list), title (str)
        save_path : str
            The path where the figure should be saved to

        """

        if key not in self.matrices.keys():
            raise Exception("Key not in List of Matrices")

        arr_mat = []

        for k, value in self.matrices[key].items():
            arr_mat.append(list(value.values()))

        figsize = params["figsize"] if "figsize" in params.keys() else (4, 8)

        fig, ax = plt.subplots(1, 1, figsize = figsize, constrained_layout = True)

        im = ax.imshow(arr_mat)

        x_axis = [self.__format(x) for x in self.matrices[key].keys()] if "x_axis" not in params.keys() else params["x_axis"]
        y_axis = [self.__format(y) for y in self.matrices[key][list(self.matrices[key].keys())[0]].keys()] if "y_axis" not in params.keys() else params["y_axis"]
        title = f"Confusion Matrix for {self.__format(key)}" if "title" not in params.keys() else params["title"]

        # Show all ticks and label them with the respective list entries
        ax.set_xticks(np.arange(len(x_axis)), labels = x_axis)
        ax.set_yticks(np.arange(len(y_axis)), labels = y_axis)

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation = 45, ha = "right",
                rotation_mode = "anchor")

        # Loop over data dimensions and create text annotations.
        for i in range(len(arr_mat)):
            for j in range(len(arr_mat[i])):
                text = ax.text(j, i, f"{100 * arr_mat[i][j]:.2f}",
                            ha = "center", va = "center", color = "w" if arr_mat[i][j] < 0.9 else "black")

        ax.set_xlabel("Predicted")
        ax.set_ylabel("Actual")
        ax.set_title(title)

        cax = ax.inset_axes([1.1, 0.1, 0.05, 0.8])

        fig.colorbar(im, ax = ax, cax = cax, pad = 1)

        fig.tight_layout()

        if save_path is not None: # saving
            fig.save(save_path)

    def classifiers(self):
        """ Method that lists classifiers available """
        classes = self.interpolator.test_classifiers(
            np.array([
            self.test_grid.initial_values[self.in_keys].tolist()
            ])[0]
        )

        return list(classes.keys())

    def keys(self):
        """ Method that lists out keys available """

        out_keys = []

        for interpolator in self.interpolator.interpolators:
            out_keys.extend(interpolator.out_keys)

        return out_keys

    def decision_boundaries(self):

        pass

    def plot2D(self, q):

        PLOT_PROPERTIES_TF12 = {
            'figsize': (4.5, 4.),
            'show_fig' : True,
            'close_fig' : True,
            'path_to_file': '/projects/b1119/prs5019/grid_slices',
            'fname': 'q_%1.1f_TF12.png'%q,
            'title' : f"Grid Slice q = {q}",
            'log10_x' : True,
            'log10_y' : True,
        }

        clean_errs = self.__clean_errs(self.errs["relative"].T[[0]].T)

        PLOT_PROPERTIES_TF1 = PLOT_PROPERTIES_TF12
        PLOT_PROPERTIES_TF1['fname'] = 'q_%1.1f_TF1_max.png'%q
        PLOT_PROPERTIES_TF1['figsize'] = (4.5,6.)
        PLOT_PROPERTIES_TF1['colorbar'] = {'label' : "Relative Error", 'vmin' : clean_errs.T[0].min(), 'vmax' : 0.00001}

        slice_3D_var_range = (q-0.05,q+0.05)
        mass_ratio = self.test_grid.initial_values["star_2_mass"] / self.test_grid.initial_values["star_1_mass"]

        slice = (mass_ratio >= slice_3D_var_range[0]) & (mass_ratio <= slice_3D_var_range[1])

        max_omega = [max(self.test_grid[i].history2['surf_avg_omega_div_omega_crit'])
                    if self.test_grid[i].history2 is not None else np.nan
                    for i in range(len(self.test_grid.MESA_dirs))]

        slice_errs = np.nanmean(self.errs["relative"], axis = 1)
        slice_errs = [slice_errs[i] if in_slice and i in self.errs["valid_inds"][0] else np.nan for i, in_slice in enumerate(slice)]

        fig = self.test_grid.plot2D('star_1_mass', 'period_days', np.array(slice_errs),
                    termination_flag='termination_flag_1',
                    grid_3D=True, slice_3D_var_str='mass_ratio',
                    slice_3D_var_range=slice_3D_var_range,
                    verbose=False, **PLOT_PROPERTIES_TF1)


        




        


    
