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

        i = self.interpolator.test_interpolator(iv) # interpolated

        self.errs = {}

        self.errs["relative"] = np.abs((fv - i) / fv)
        self.errs["absolute"] = np.abs(fv - i)

        valid_inds = np.where(
            self.test_grid.final_values["interpolation_class"] != "not_converged"
        )[0]
        cfv = self.test_grid.final_values[valid_inds]

        c = self.interpolator.test_classifiers(iv[valid_inds]) # classifying

        # computing confusion matrices
        self.matrices = {}

        for key, value in c.items():

            labels = cfv[key]

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
        
        errs = self.errs[err_type].T[k_inds].T

        errs = errs[~np.isnan(errs).any(axis = 1)] # dropping nans
        print(k_inds, np.nanmedian(errs, axis = 1))

        plt.rcParams.update({"font.size": 18, "font.family": "Times New Roman"})

        fig, axs = plt.subplots(1, 1,
                                figsize = (24, 10),
                                tight_layout = True)
        
        parts = axs.violinplot(np.log10(errs + 1.0e-8), showmedians = True, points = 1000)
        axs.set_title(f"Distribution of {err_type.capitalize()} Errors")
        axs.set_xticks(np.arange(1, len(keys) + 1), 
            labels = [
                f"{self.__format(ec)} ({(med * 100):.2f})" for ec, med in zip(keys, np.nanmedian(errs, axis = 1))
            ], rotation = 20)

        axs.set_ylabel("Errors in Log 10 Scale")
        axs.grid(axis = "y")

        for pc in parts["bodies"]:
            pc.set_facecolor("#D43F3A")
            pc.set_edgecolor("black")
            pc.set_alpha(0.85)

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

    def plot2D(self, x_var_str, y_var_str, z_var_str = None,
               termination_flag = "termination_flag_1",
               grid_3D = True, slice_3D_var_str = "star_2_mass", slice_3D_var_range = None,
               grid_4D = False, slice_4D_var_str = None, slice_4D_var_range = None,
               extra_grid = None, slice_at_RLO = False,
               MARKERS_COLORS_LEGENDS = None,
               verbose = False, **kwargs):
        """Plot a 2D slice of x_var_str vs y_var_str of one or more runs.

        Parameters
        ----------
        x_var_str : str
            String of the initial value to plot on the x axis. Allowed strings
            are `psygrid.initial_values.dtype.names`.
        y_var_str : str
            String of the initial value to plot on the y axis. Allowed strings
            are `psygrid.initial_values.dtype.names`.
        z_var_str : str
            String of the initial value to plot on the z axis (displayed as
            a color). Allowed strings are
            `psygrid.final_values.dtype.names`, `psygrid.history1.dtype.names`,
            `psygrid.binary_history.dtype.names`.
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
        slice_4D_var_range : toople
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
        verbose : bool
            If `True`, the object reports by printing to standard output.
        **kwargs : dict
            Dictionary containing extra visualisation options (cf.
            `PLOT_PROPERTIES` in `plot_defaults.py`.

        """
        plot = plot2D(psygrid = self.test_grid,
                      x_var_str = x_var_str,
                      y_var_str = y_var_str,
                      z_var_str = z_var_str,
                      termination_flag = termination_flag,
                      grid_3D = grid_3D,
                      slice_3D_var_str = slice_3D_var_str,
                      slice_3D_var_range = slice_3D_var_range,
                      grid_4D = grid_4D,
                      slice_4D_var_str = slice_4D_var_str,
                      slice_4D_var_range = slice_4D_var_range,
                      extra_grid = extra_grid,
                      slice_at_RLO = slice_at_RLO,
                      MARKERS_COLORS_LEGENDS = MARKERS_COLORS_LEGENDS,
                      verbose = verbose,
                      **kwargs)

        plot()


        




        


    
