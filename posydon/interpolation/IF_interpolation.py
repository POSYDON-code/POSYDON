"""Module for performing initial-final inteprolation.

We showcase the initial-final interpolator which plays a critical role in the
evolving binary populations. To use the initial-final interpolator we first
import the IFInterpolator class from the POSYDON library.


  # importing interpolator
  from posydon.interpolation.IF_interpolation import IFInterpolator

I. Loading a Pretrained Interpolator
------------------------------------

To load a pretrained interpolator we need to
pass in the filename argument into the IFInterpolator
constructor which specifies the path to a .pkl file where
the pretrained interpolator can be loaded from. POSYDON provides
various pretrained models whose corresponding .pkl files
can be found in the data directory of the POSYDON repository.


  model = IFInterpolator()

  model.load("path/to/file.pkl")


II. Training the Interpolator
-------------------------

The interpolator can be trained using an instance of the PSyGrid
class which can be constructed by running ones own simulations or
by loading a simulation from an h5 file stored in the data directory
of the POSYDON repository. For more details on the PSyGrid class
please visit the PSyGrid documentation. The IFInterpolator
class relies on the BaseIFInterpolator class to perform the interpolation
so parameters to construct instances of the BaseIFInterpolator classes are
required to construct the IFInterpolator. These parameters are:

1. interp_method: the interpolation method to be used can be either
linear, 1NN, or a list specifying which interpolation method to
be used for each type of track. If interp_classes is specified and this
parameter is not a list then the interpolator will use the specified method
for all classes.

2. interp_classes: a list of classes that the simulation tracks can
fall into. Usually specified as the mass transfer type. This only needs
be specified if interp_method is a list.

3. class_method: the classification method to be used, either kNN or
1NN.

4. in_keys: the keys to be used as the input to the interpolator, by default
these are star_1_mass, star_2_mass, and period_days.

5. out_keys: the keys for which the interpolator is supposed to provide values,
by default all keys are used.

6. in_scaling: The scalings for the input keys, by default these scalings are
optimized through Monte Carlo Cross Validation.

7. out_scaling: The scalings for the output keys, by default these scalings
are optimized through Monte Carlo Cross Validation.

8. c_keys: A list of strings specifying which classifiers are to be trained

9. c_key: A string specifying by which class the interpolator should
interpolate binaries. Only to be specified in the MCInterpolator case.

For most applications specifying only the first four parameters is recommended.


    from posydon.grids.psygrid import PSyGrid
    from posydon.interpolation.IF_interpolation import IFInterpolator

    grid = PSyGrid("path/to/h5/file.h5") # loading grid from h5 file

    interp = IFInterpolator(grid = grid, interpolators = [
        {
            "interp_method": ["linear", "linear", "linear"],
            "interp_classes": ["no_MT", "stable_MT", "unstable_MT"],
            "out_keys": first,
            "class_method": "kNN",
            "c_keys": ["interpolation_class"],
            "c_key": "interpolation_class"
        },
        {
            "interp_method": ["linear", "linear", "linear", "linear"],
            "interp_classes": ["None", "BH", "WD", "NS"],
            "out_keys": second,
            "class_method": "kNN",
            "c_keys": ['S1_direct_state'],
            "c_key": 'S1_direct_state'
        },
        {
            "interp_method": ["linear", "linear", "linear", "linear"],
            "interp_classes": ["None", "BH", "WD", "NS"],
            "out_keys": third,
            "class_method": "kNN",
            "c_keys": ['S1_Fryer+12-rapid_state'],
            "c_key": 'S1_Fryer+12-rapid_state'
        },
        {
            "interp_method": ["linear", "linear", "linear", "linear"],
            "interp_classes": ["None", "BH", "WD", "NS"],
            "out_keys": fourth,
            "class_method": "kNN",
            "c_keys": ['S1_Fryer+12-delayed_state'],
            "c_key": 'S1_Fryer+12-delayed_state'
        },
        {
            "interp_method": ["linear", "linear", "linear", "linear"],
            "interp_classes": ["None", "BH", "WD", "NS"],
            "out_keys": fifth,
            "class_method": "kNN",
            "c_keys": ['S1_Sukhbold+16-engineN20_state'],
            "c_key": 'S1_Sukhbold+16-engineN20_state'
        },
        {
            "interp_method": ["linear", "linear", "linear", "linear"],
            "interp_classes": ["None", "BH", "WD", "NS"],
            "out_keys": sixth,
            "class_method": "kNN",
            "c_keys": ['S1_Patton&Sukhbold20-engineN20_state'],
            "c_key": 'S1_Patton&Sukhbold20-engineN20_state'
        }
    ]) # constructing IFInterpolator

    interp.train() # training interpolator


III. Using the Interpolator
---------------------------

Once the interpolator has been trained or loaded from a .pkl file it can be
used to accomplish various tasks which most commonly are to classify a track
into its class given an input vector and or to approximate a final vector given
an input vector.


    from posydon.binary_evol.binarystar import BinaryStar
    from posydon.binary_evol.singlestar import SingleStar


    # creating binary, refer to BinaryStar documentation
    binary = BinaryStar(**binary_params,
                        star_1=SingleStar(**star1_params),
                        star_2=SingleStar(**star2_params))

    # evaluating returns a tuple of dictionaries
    interpolation, classification = interp.evaluate(binary)

Finally a trained interpolator can be easily saved by specifying a path to a
.pkl file where the interpolator will be saved to.


   model.save("path/to/file.pkl") # saving interpolator



"""


__authors__ = [
    "Juanga Serra Perez <jgserra@northwestern.edu>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Philipp Moura Srivastava <philipp.msrivastava@northwestern.edu>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
]


import os
import pickle
import warnings
from datetime import date
# POSYDON
from posydon.grids.psygrid import PSyGrid
from posydon.interpolation.data_scaling import DataScaler
# Maths
import numpy as np
# Machine Learning
from scipy.interpolate import LinearNDInterpolator
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import confusion_matrix
# Plotting
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap
from posydon.visualization.plot_defaults import DEFAULT_LABELS
from posydon.interpolation.constraints import (
    find_constraints_to_apply, sanitize_interpolated_quantities)


# INITIAL-FINAL INTERPOLATOR
class IFInterpolator:
    """Class handling initial-final interpolation.

    Class handling initial-final interpolation with support to interpolate
    certain keys by differing classes.
    """

    def __init__(self, grid=None, interpolators=None):
        """Initialize the IFInterpolator class.

        Parameters
        ----------
        grid : PSyGrid
            The training grid
        interpolators : list of dictionaries
            Contain parameters for the BaseIFIneterpolators to be constructed

        """
        self.interpolators = []
        self.interpolator_parameters = interpolators
        self.grid = grid

        if self.interpolator_parameters is not None:
            out_keys = []

            for params in self.interpolator_parameters:

                if ("out_keys" not in params
                        and len(self.interpolator_parameters) > 1):
                    raise Exception("Overlapping out keys between different "
                                    "interpolators are not permited!")
                elif "out_keys" not in params:
                    continue

                for key in params["out_keys"]:
                    if(key in out_keys):
                        raise Exception(
                            f"Overlapping out keys between different "
                            f"interpolators are not permited! ({key} in more "
                            f"than one set of out_keys)")
                    else:
                        out_keys.append(key)
        

    def train(self):
        """Train the interpolator(s) on the PSyGrid used for construction."""
        for interpolator in self.interpolator_parameters:
            self.interpolators.append(BaseIFInterpolator(grid=self.grid,
                                                         **interpolator))

        self.interp_in_q = self.interpolators[0].interp_in_q

    def evaluate(self, binary, sanitization_verbose=False):
        """Get the output vector approximation.

        Get the output vector approximation as well as all of the
        classifications given a Binary Star

        Parameters
        ----------
        binary : BinaryStar
            A class which is an input vector which are to be classified and for
            which an output vector will be approximated.

        Returns
        -------
        Output space approximation as a tuple containing two dictionary

        """
        ynums = {}
        ycats = {}

        for interpolator in self.interpolators:
            ynum, ycat = interpolator.evaluate(binary, sanitization_verbose)

            ynums = {**ynums, **ynum}
            ycats = {**ycats, **ycat}

        return ynums, ycats


    def test_interpolator(self, initial_values, sanitization_verbose = False):
        """ Method that can take in a 2-D numpy array for more efficient use of the interpolator

        Parameters
        ----------
        initial_values : numpy array
            A numpy array containing the in-key values of the binaries to be evolved 

        Return
        ------
        Interpolated values
        """
        new_values = np.array(initial_values, copy = True)

        if self.interp_in_q:
            new_values.T[1] = new_values.T[1] / new_values.T[0]

        final_values = []

        for interpolator in self.interpolators:
            final_values.append(interpolator.test_interpolator(new_values).T)

        return np.concatenate(final_values).T

    def test_classifiers(self, initial_values):
        """ Method that can take in a 2-D numpy array for more efficient use of classifiers

        Parameters
        ----------
        initial_values : numpy array
            A numpy array containing the in-key values of the binaries to be classified 

        Return
        ------
        Classified values
        """
        new_values = np.array(initial_values, copy = True)

        if self.interp_in_q:
            new_values.T[1] = new_values.T[1] / new_values.T[0]

        classes = {}

        for interpolator in self.interpolators:

            classes = {**classes, **interpolator.test_classifiers(new_values)}
            
        return classes



    def load(self, filename):
        """Load a saved IFInterpolator from a .pkl file.

        Parameters
        ----------
        filename : string
            Path to the .pkl file

        """
        with open(filename, 'rb') as f:
            self.interpolators = pickle.load(f)

        self.interp_in_q = self.interpolators[0].interp_in_q

    def save(self, filename):
        """Save the interpolator to a .pkl file.

        Parameters
        ----------
        filename : string
            Path where .pkl file will be saved

        """
        with open(filename, "wb") as f:
            pickle.dump(self.interpolators, f)


class BaseIFInterpolator:
    """Class handling the initial-final interpolation for POSYDON."""

    def __init__(self, grid=None, in_keys=None, out_keys=None, in_scaling=None,
                 out_scaling=None, filename=None, interp_method="linear",
                 interp_classes=None, class_method="kNN",
                 c_keys=None, c_key=None):
        """Initialize the BaseIFInterpolation.

        Initialize the BaseIFInterpolation and if 'filename' is provided load
        a pretrained interpolator. If no filename is provided, the
        interpolator(s) and classifier(s) are trained upon initialization.

        Parameters
        ----------
        grid : PSyGrid
            An instance of the PSyGrid class.
        in_keys : list of strings
            The keys to be used as the input to the interpolator, by default
            these are star_1_mass, star_2_mass, and period_days.
        out_keys : list of strings
            The keys for which the interpolator is supposed to provide values,
            by default all keys are used.
        in_scaling : list of strings
            The scalings for the input keys, by default these scalings are
            optimized through Monte Carlo Cross Validation.
        out_scaling : list of strings
            The scalings for the output keys, by default these scalings
            are optimized through Monte Carlo Cross Validation.
        filename : string
            A path to a .pkl file containing a saved interpolator
        interp_method : string or list of strings
            The interpolation method to be used can be either
            linear, 1NN, or a list specifying which interpolation method to be
            used for each type of track. If interp_classes is specified and
            this parameter is not a list then the interpolator will use the
            specified method for all classes.
        interp_classes : list of strings
            A list of classes that the simulation tracks can
            fall into. Usually specified as the mass transfer type. This only
            needs be specified if interp_method is a list.
        class_method : string or list of strings
            The classification method to be used, either kNN or
            1NN.
        c_keys: list of strings
            The keys for which a classifier should be trained. At a minimum
            should contain the keys needed for multi-class interpolation if it
            is requested

        """
        self.in_keys, self.out_keys = in_keys, out_keys
        self.out_nan_keys = []
        self.in_scaling, self.out_scaling = in_scaling, out_scaling
        self.n_in, self.n_out = 0, 0
        self.in_scalers, self.out_scalers = [], []
        self.interp_in_q = False
        self.N = []
        self.valid = None
        self.c_key = c_key if c_key is not None else "interpolation_class"

        self.interp_method = interp_method
        self.interpolator = None
        if interp_classes is not None:
            if not isinstance(interp_classes, list):
                raise Exception("interp_classes must be a list of "
                                "valid interpolation methods.")
            if isinstance(self.interp_method, list):
                if len(self.interp_method) != len(interp_classes):
                    raise Exception("No. of interpolation methods must "
                                    "match no. of interpolation classes.")
            else:
                self.interp_method = [self.interp_method] * len(interp_classes)
            self.interp_classes = interp_classes
        self.class_method = class_method
        self.classifiers = None

        if filename is not None:
            self.load(filename)

        else:
            if grid is None:
                raise Exception(
                    "grid must be specified to create an IF interpolator."
                    " Conversely, load an existing model specifying filename.")

            if self.in_keys is None:
                self.in_keys = ['star_1_mass', 'star_2_mass', 'period_days']
            if self.out_keys is None:
                self.out_keys = [
                    key for key in grid.final_values.dtype.names
                    if key != "model_number"
                    and (type(grid.final_values[key][0]) != np.str_)
                    and any(~np.isnan(grid.final_values[key]))
                ]
                self.out_nan_keys = [
                    key for key in grid.final_values.dtype.names
                    if type(grid.final_values[key][0]) != np.str_
                    and all(np.isnan(grid.final_values[key]))
                ]

            self.constraints = find_constraints_to_apply(self.out_keys)

            if c_keys is None:
                c_keys = [key for key in grid.final_values.dtype.names
                          if type(grid.final_values[key][0]) == np.str_]

            self.classifiers = {key: None for key in c_keys}
            self.n_in, self.n_out = len(self.in_keys), len(self.out_keys)
            self.interp_in_q = self._interpIn_q(grid)
            self.XT, self.YT = self._grid2array(grid)
            self.valid = self._setValid(
                grid.final_values[self.c_key], self.XT)
            self.N = np.sum(self.valid >= 0)

            analize_nans(self.out_keys, self.YT, self.valid)
            analize_nones(self.classifiers, self.valid, grid)

            if (self.interp_method == 'linear'
                    or isinstance(self.interp_method, list)):
                print("\nFilling missing values (nans) with 1NN")
                self._fillNans()
            elif self.interp_method == '1NN':
                self.YT[np.isnan(self.YT)] = -100

            if (self.in_scaling is None) or (self.out_scaling is None):
                if self.interp_method == '1NN':
                    self.in_scaling, self.out_scaling = (self._bestInScaling(),
                                                         ['none']*self.n_out)
                else:
                    self.in_scaling, self.out_scaling = self._bestScaling()

            self.X_scaler = MatrixScaler(self.in_scaling,
                                         self.XT[self.valid >= 0, :])
            self.Y_scaler = MatrixScaler(self.out_scaling,
                                         self.YT[self.valid > 0, :])

            if self.class_method == "kNN":
                options = {'nfolds': 3, 'p_test': 0.05, 'nmax': 10}
            else:
                options = {}
            self.train_classifiers(grid, method=self.class_method, **options)

            self.train_interpolator(
                ic=grid.final_values[self.c_key])

    def save(self, filename):
        """Save complete interpolation model.

        Parameters
        ----------
        filename : str
            path/name of '.pkl' file where the model will be saved.

        """
        SAVE_ATTRS = self.__dict__
        DONT_SAVE = []
        myattrs = {key: SAVE_ATTRS[key]
                   for key in SAVE_ATTRS if key not in DONT_SAVE}

        with open(filename, 'wb') as f:
            pickle.dump(myattrs, f)

    def load(self, filename):
        """Load interpolation model, which can only be used for predictions.

        Parameters
        ----------
        filename : str
            path/name of pickle file to be loaded.

        """
        with open(filename, 'rb') as f:
            myattrs = pickle.load(f)
            for key in myattrs:
                setattr(self, key, myattrs[key])

    def _grid2array(self, grid):
        """Convert a PSyGrid to a numpy array.

        Parameters
        ----------
        grid : PSyGrid
            a loaded PSyGrid object

        Returns
        -------
        Two numpy arrays with the first being the input space and
        the second being the output space

        """
        self.n_in, self.n_out = len(self.in_keys), len(self.out_keys)
        self.N = len(grid.initial_values[self.in_keys[0]])
        X = np.empty((self.N, self.n_in))
        Y = np.empty((self.N, self.n_out))
        for i in range(self.n_in):
            X[:, i] = grid.initial_values[self.in_keys[i]]
        for i in range(self.n_out):
            Y[:, i] = grid.final_values[self.out_keys[i]]

        if self.interp_in_q:
            i_m1 = self.in_keys.index('star_1_mass')
            i_m2 = self.in_keys.index('star_2_mass')
            X[:, i_m2] = X[:, i_m2] / X[:, i_m1]

        return X, Y

    def _binary2array(self, binary):
        """Convert a binary instance to a numpy array.

        Parameters
        ----------
        binary : BinaryStar
            An initialized instance of the Binary star class

        Returns
        -------
        A binary as a numpy array

        """
        var2 = binary.star_2.mass
        if self.interp_in_q:
            var2 = binary.star_2.mass / binary.star_1.mass
        Xt = np.array([[binary.star_1.mass, var2, binary.orbital_period]])

        return Xt

    def _setValid(self, ic, X):
        """Set binary tracks as (in)valid depending on termination flags."""
        valid = np.zeros(len(ic), dtype=int)

        for flag in ['not_converged', 'ignored_no_BH', 'ignored_no_RLO']:
            which = (ic == flag)
            valid[which] = -1
            print(f"Discarded {np.sum(which)} binaries with "
                  f"interpolation_class = [{flag}]")

        which = np.isnan(np.sum(X, axis=1))
        valid[which] = -1
        print(f"Discarded {np.sum(which)} binaries with nans in input values.")

        for i, flag in enumerate(self.interp_classes):
            valid[ic == flag] = i + 1

        if(self.interp_in_q):   # if HMS-HMS grid, take out q = 1
            for i, iv in enumerate(X):
                if(iv[1] == 1):
                    valid[i] = -1

        return valid

    def train_interpolator(self, ic=None):
        """Train the interpolator.

        Parameters
        ----------
        ic : numpy array
            Contains the interpolation class for each track in the grid used
            for training the interpolator.

        """
        XTn = self.X_scaler.normalize(self.XT[self.valid > 0, :])
        YTn = self.Y_scaler.normalize(self.YT[self.valid > 0, :])

        if self.interp_method == "linear":
            self.interpolator = LinInterpolator()
            self.interpolator.train(XTn, YTn)
        elif self.interp_method == "1NN":
            self.interpolator = NNInterpolator()
            self.interpolator.train(XTn, YTn)
        elif isinstance(self.interp_method, list):
            self.interpolator = MC_Interpolator(
                self.classifiers[self.c_key],
                self.interp_classes, self.interp_method)
            self.interpolator.train(XTn, YTn, ic[self.valid > 0])

    def test_interpolator(self, Xt):
        """Use the interpolator to approximate output vector.

        Parameters
        ----------
        Xt : numpy array
            A list of input vectors for which the output vectors are
            to be approximated.

        Returns
        -------
        Output space approximation as numpy array

        """
        Xtn = self.X_scaler.normalize(Xt)
        Ypredn = self.interpolator.predict(Xtn)

        Ypredn = np.array([
            list(sanitize_interpolated_quantities(
                dict(zip(self.out_keys, track)),
                self.constraints, verbose=False).values())
            for track in self.Y_scaler.denormalize(Ypredn)
        ])

        return Ypredn

    def train_classifiers(self, grid, method='kNN', **options):
        """Train the classifiers.

        Paramaters
        ----------

        grid : PSyGrid
            The training grid
        method : string
            Specifies the classification method, default is kNN

        """
        for key in self.classifiers:
            self.train_classifier(grid, key, method, **options)

    def test_classifiers(self, Xt):
        """Use the classifiers.

        Parameters
        ----------
        Xt : numpy array
            A list of input vectors which are to be classified.

        """
        return {key: self.test_classifier(key, Xt)
                if self.classifiers[key] is not None else None
                for key in self.classifiers}

    def train_classifier(self, grid, key, method='kNN', **options):
        """Train a specific classifier.

        Paramaters
        ----------
        grid : PSyGrid
            The training grid
        method : string
            Specifies the classification method, default is kNN

        """
        v = self.valid >= 0
        wnn = grid.final_values[key] != 'None'
        if any(wnn[v]):
            if (method == 'kNN') or (method == '1NN'):
                self.classifiers[key] = KNNClassifier()
            else:
                raise ValueError("Wrong method name.")
            # If we want to exclude Nones as a class:
            # XTn = self.X_scaler.normalize(self.XT[v & wnn, :])
            # y = grid.final_values[key][v & wnn]
            XTn = self.X_scaler.normalize(self.XT[v, :])
            y = grid.final_values[key][v]
            uy = np.unique(y)
            print(f"\nTraining {method} classifier for {key} "
                  f"[nclasses = {len(uy)}]...")
            for c in uy:
                print(f" - [{np.sum(y == c)}] {c}")
            if method == '1NN':
                self.classifiers[key].train(XTn, y, K=1)
            elif method == 'kNN':
                self.classifiers[key].xtrain(XTn, y)

    def test_classifier(self, key, Xt):
        """Use just one specific classifier.

        Parameters
        ----------
        key : string
            Name of the classifier
        Xt : numpy array
            A list of input vectors which are to be classified.

        Returns
        -------
        Output space approximation as numpy array

        """
        Xn = self.X_scaler.normalize(Xt)

        return self.classifiers[key].predict(Xn)

    def test_classifier_prob(self, key, Xt):
        """Test the classifier.

        Use just one specific classifier to get the probabilities of the input
        being in a class

        Parameters
        ----------
        key : string
            Name of the classifier
        Xt : numpy array
            A list of input vectors which are to be classified.

        Returns
        -------
        Output space approximation as numpy array

        """
        Xn = self.X_scaler.normalize(Xt)
        return self.classifiers[key].predict_prob(Xn)

    def evaluate_mat(self, Xt):
        """Get the output vector approximation.

        Get the output vector approximation as well as all of the
        classifications given an input vector.

        Parameters
        ----------
        Xt : numpy array
            A list of input vectors which are to be classified and for which
            output vectors are to be approximated.

        Returns
        -------
        Output space approximations as numpy arrays

        """
        if len(Xt.shape) == 1:
            Xt = Xt.reshape((1, -1))
        if Xt.shape[1] != self.n_in:
            raise Exception("Wrong dimensions. Xt should have as many "
                            "columns as it was trained with.")
        # if binary classified as 'initial_MT', set numerical quantities to nan
        ynum, ycat = self.test_interpolator(Xt), self.test_classifiers(Xt)
        if self.class_method != '1NN':
            ynum[ycat[self.c_key] == 'initial_MT', :] = np.nan
        if self.interp_method == '1NN':
            ynum[ynum == - 100] = np.nan

        return ynum, ycat

    def evaluate(self, binary, sanitization_verbose=False):
        """Get the output vector approximation.

        Get the output vector approximation as well as all of the
        classifications given a Binary Star.

        Parameters
        ----------
        binary : BinaryStar
            A class which is an input vector which are to be classified and for
            which an output vector will be approximated.

        Returns
        -------
        Output space approximation as a tuple containing two dictionary

        """
        ynum, ycat = self.evaluate_mat(self._binary2array(binary))
        for key in ycat:
            if ycat[key] is not None:
                ycat[key] = ycat[key][0]

        yt = dict(zip(self.out_keys, ynum.flatten()))  # for a single binary
        yt.update(dict(zip(self.out_nan_keys,
                           [np.nan] * len(self.out_nan_keys))))

        # yt = sanitize_interpolated_quantities(yt,
        #                                       verbose=sanitization_verbose)

        # yt = np.append(ynum, np.array([np.nan] * len(self.out_nan_keys)))

        return yt, ycat

    def _interpIn_q(self, grid):
        """Find if the (regular) grid was uniformly sampled in M_2 or q.

        Parameters
        ----------
        grid : posydon.grids.psygrid.PSyGrid
            Grid of binaries to convert to data matrix.

        Returns
        -------
        bool
            True if the grid was sampled in q.

        """
        m1 = grid.initial_values['star_1_mass'].copy()
        m2 = grid.initial_values['star_2_mass'].copy()
        # Make sure nans do not affect
        wnan = np.isnan(m1) | np.isnan(m2)
        m1, m2 = m1[~wnan], m2[~wnan]

        tol = 1

        return (m2[m1 > 0.95 * m1.max()].min()
                / m2[m1 < 1.05 * m1.min()].min() > 1 + tol)

    def _bestInScaling(self):
        """Find the best scaling for the input space."""
        in_scaling = []  # I assume inputs are positive-valued
        for i in range(self.n_in):
            if (np.abs(np.mean(self.XT[:, i]) - np.median(self.XT[:, i]))
                    / np.std(self.XT[:, i])
                    < np.abs(np.mean(np.log10(self.XT[:, i]))
                             - np.median(np.log10(self.XT[:, i])))
                    / np.std(np.log10(self.XT[:, i]))):
                in_scaling.append('min_max')
            else:
                in_scaling.append('log_min_max')

        return in_scaling

    def _bestScaling(self, unique_in=True):
        """Find the best scaling for both input and output space."""
        # Decide best scaling linear/log with cross validation
        nfolds = 5
        p_test = 0.15

        # if False, the scaling for the inputs will be output-dependent
        unique_in = True

        if unique_in:
            r = 1
        else:
            r = 2 ** self.n_in

        err = np.nan * np.ones((r * 2, self.n_out, nfolds))
        which_abs = np.abs(self.YT[self.valid > 0, :]).min(axis=0) == 0

        out_scalings = [['min_max'] * self.n_out, ['log_min_max'] * self.n_out]

        for i, key in enumerate(self.out_keys):
            y = self.YT[self.valid > 0, i]
            if np.nanmin(y) == np.nanmax(y):
                out_scalings[0][i] = 'none'
                out_scalings[1][i] = 'none'
            elif np.nanmax(y) < 0:
                out_scalings[1][i] = 'neg_log_min_max'
            elif np.nanmin(y) <= 0:
                out_scalings[1][i] = 'min_max'

        in_scaling = self._bestInScaling()

        XT = self.XT[self.valid > 0, :]
        YT = self.YT[self.valid > 0, :]
        for j in range(2):          # normal and log when possible
            for i in range(r):      # valid also with metallicity included
                if r > 1:
                    in_scaling = [
                        'min_max'
                        if i % (2 ** (j + 1)) < 2 ** j else 'log_min_max'
                        for j in range(self.n_in)
                    ]

                xs, ys = (MatrixScaler(in_scaling, XT),
                          MatrixScaler(out_scalings[j], YT))
                XTn, YTn = xs.normalize(XT), ys.normalize(YT)

                np.random.seed(0)

                for fold in range(nfolds):
                    iTrain, itest = xval_indices(np.sum(self.valid > 0),
                                                 percent_test=p_test)

                    X_T, Y_T = XTn[iTrain, :], YTn[iTrain, :]
                    X_t, Y_t = XT[itest, :], YT[itest, :]

                    interp = LinearNDInterpolator(X_T, Y_T)
                    ypred_i = ys.denormalize(interp(xs.normalize(X_t)))

                    err[i + r * j, ~which_abs, fold] = np.nanpercentile(
                        np.abs((ypred_i[:, ~which_abs] - Y_t[:, ~which_abs])
                               / Y_t[:, ~which_abs]), 90, axis=0)
                    err[i + r * j, which_abs, fold] = np.nanpercentile(
                        np.abs(ypred_i[:, which_abs] - Y_t[:, which_abs]), 90,
                        axis=0)

        where_min = np.nanargmin(np.nanmean(err, axis=2), axis=0)

        out_scaling = []
        for i in range(self.n_out):
            if where_min[i] < r:
                out_scaling.append(out_scalings[0][i])
            else:
                out_scaling.append(out_scalings[1][i])
                where_min[i] -= r

        return in_scaling, out_scaling

    def _fillNans(self):
        """Fill nan values i numerical magnitudes with 1NN."""
        for i in range(self.n_out):
            wnan = np.isnan(self.YT[:, i]) | np.isinf(self.YT[:, i])
            if any(np.isinf(self.YT[:, i])):
                print(f"inftys: {np.sum(np.isinf(self.YT[:, i]))}, "
                      f"{self.out_keys[i]}")
            if any(wnan[self.valid > 0]):
                k1r = KNeighborsRegressor(n_neighbors=1)
                wT = (~wnan) & (self.valid > 0)
                xs = MatrixScaler(self._bestInScaling(),
                                  self.XT[self.valid >= 0, :])
                k1r.fit(xs.normalize(self.XT[wT, :]), self.YT[wT, i])
                wt = wnan & (self.valid > 0)
                self.YT[wt, i] = k1r.predict(xs.normalize(self.XT[wt, :]))


# BASE INTERPOLATORS


class Interpolator:
    """Class used as a shell parent class for all Interpolators."""

    def __init__(self):
        """Initialize the Interpolator."""
        self.interpolator = None
        self.D = None

    def train(self, XT, YT):
        """Train the interpolator.

        Parameters
        ----------
        XT : numpy array
            List of input vectors
        YT : numpy array
            List of corresponding output vectors

        """
        self.D = XT.shape[1]

    def predict(self, Xt):
        """Interpolate and approximate output vectors given input vectors.

        Parameters
        ----------
        XT : numpy array
            List of input vectors

        """
        if self.interpolator is None:
            raise Exception("Train Interpolator first.")
        if Xt.shape[1] != self.D:
            raise Exception("Wrong input dimension.")

    def train_error(self, XT, YT):
        """Calculate approximation error given testing data.

        Parameters
        ----------
        XT : numpy array
            List of input vectors
        YT : numpy array
            List of corresponding output vectors
        """
        error = np.abs(YT - self.predict(XT))

        print("\nMean Abs. Train Error:")
        print(error.mean(axis=0))


class NNInterpolator(Interpolator):
    """Class implementing Nearest Neighbor interpolation."""

    def train(self, XT, YT):
        """Train the interpolator.

        Parameters
        ----------
        XT : numpy array
            List of input vectors
        YT : numpy array
            List of corresponding output vectors

        """
        super().train(XT, YT)
        self.interpolator = KNeighborsRegressor(n_neighbors=1)
        self.interpolator.fit(XT, YT)
        # self.train_error(XT, YT)

    def predict(self, Xt):
        """Interpolate and approximate output vectors given input vectors.

        Parameters
        ----------
        XT : numpy array
            List of input vectors

        Returns
        -------
        Output space approximation as numpy array

        """
        super().predict(Xt)
        return self.interpolator.predict(Xt)


class LinInterpolator(Interpolator):
    """Class implementing Linear Interpolation."""

    def train(self, XT, YT):
        """Train the interpolator.

        Parameters
        ----------
        XT : numpy array
            List of input vectors
        YT : numpy array
            List of corresponding output vectors

        """
        super().train(XT, YT)
        self.interpolator = [LinearNDInterpolator(XT, YT)]
        OoH = KNeighborsRegressor(n_neighbors=1)
        OoH.fit(XT, YT)
        self.interpolator.append(OoH)
        # self.train_error(XT, YT)

    def predict(self, Xt):
        """Interpolate and approximate output vectors given input vectors.

        Parameters
        ----------
        XT : numpy array
            List of input vectors

        Returns
        -------
        Output space approximation as numpy array

        """
        super().predict(Xt)
        Ypred = self.interpolator[0](Xt)
        wnan = np.isnan(Ypred[:, 0])
        # 1NN interpolation for binaries out of hull
        if np.any(wnan):
            Ypred[wnan, :] = self.interpolator[1].predict(Xt[wnan, :])
        if np.any(wnan):
            warnings.warn(f"1NN interpolation used for {np.sum(wnan)} "
                          "binaries out of hull.")
        return Ypred


class MC_Interpolator:
    """Class implementing class-wise interpolation."""

    def __init__(self, classifier, classes, methods):
        """Initialize the class-wise interpolation.

        classifier : KNNClassifier
            The classifier that is used to classify input vectors
        classes : List of strings
            The classes in question
        methods : List of strings
            The methods to be used to interpolate between each class of tracks
        """
        self.interpolators = [None] * len(classes)
        self.D = classifier.D
        self.classifier = classifier
        self.classes = []
        for c in classes:
            if isinstance(c, list):
                self.classes.append(c)
            else:
                self.classes.append([c])
        self.M = None
        if not isinstance(methods, list):
            self.methods = [methods] * len(self.classes)
        else:
            if len(methods) != len(self.classes):
                raise Exception("No. of interpolators must match "
                                "no. of classes.")
            self.methods = methods
        for i, method in enumerate(self.methods):
            if method == "linear":
                self.interpolators[i] = LinInterpolator()
            elif method == "1NN":
                self.interpolators[i] = NNInterpolator()
            else:
                raise AttributeError(
                    "Valid interpolation methods are 'linear' and '1NN'.")

    def train(self, XT, YT, z):
        """Train the interpolator.

        Parameters
        ----------
        XT : numpy array
            List of input vectors
        YT : numpy array
            List of corresponding output vectors

        """
        self.M = YT.shape[1]
        for i in range(len(self.classes)):
            which = np.zeros_like(z, dtype=bool)
            for j in range(len(self.classes[i])):
                which += z == self.classes[i][j]
            self.interpolators[i].train(XT[which, :], YT[which, :])

    def predict(self, Xt):
        """Interpolate and approximate output vectors given input vectors.

        Parameters
        ----------
        XT : numpy array
            List of input vectors

        Returns
        -------
        Output space approximation as numpy array

        """
        zpred = self.classifier.predict(Xt)
        Ypred = np.ones((Xt.shape[0], self.M)) * np.nan
        for i in range(len(self.classes)):
            which = np.zeros_like(zpred, dtype=bool)
            for j in range(len(self.classes[i])):
                which += zpred == self.classes[i][j]
            Ypred[which, :] = self.interpolators[i].predict(Xt[which, :])
        return Ypred


# BASE CLASSIFIERS


class Classifier:
    """Class used as a shell parent class for all classifiers."""

    def __init__(self):
        """Initialize the Classifier."""
        self.classifier = None
        self.D = None
        self.labels = None
        self.method = None

    def train(self, XT, yT):
        """Train the classifier.

        Parameters
        ----------
        XT : numpy array
            List of input vectors
        YT : numpy array
            List of corresponding classes

        """
        assert XT.shape[0] == len(yT)
        self.D = XT.shape[1]

    def predict(self, Xt):
        """Classify and approximate classes given input vectors.

        Parameters
        ----------
        XT : numpy array
            List of input vectors

        Returns
        -------
        Output space approximation as numpy array

        """
        if self.classifier is None:
            raise Exception("Train Classifier first.")
        if Xt.shape[1] != self.D:
            raise Exception("Wrong input dimension.")

    def predict_prob(self, Xt):
        """Classify and get probability of input vector belonging to any class.

        Parameters
        ----------
        XT : numpy array
            List of input vectors

        Returns
        -------
        Output space approximation as numpy array

        """
        if self.classifier is None:
            raise Exception("Train Classifier first.")
        if Xt.shape[1] != self.D:
            raise Exception("Wrong input dimension.")

    def train_error(self, XT, yT):
        """Calculate approximation error given testing data.

        Parameters
        ----------
        XT : numpy array
            List of input vectors
        YT : numpy array
            List of corresponding classes
        """
        error = np.sum(yT == self.predict(XT))

        print("\nOverall Train Accuracy:")
        print(error / XT.shape[0])


class KNNClassifier(Classifier):
    """Class implementing K - Nearest Neighbors classifier."""

    def train(self, XT, yT, K=5):
        """Train the classifier.

        Parameters
        ----------
        XT : numpy array
            List of input vectors
        YT : numpy array
            List of corresponding classes

        """
        super().train(XT, yT)
        self.classifier = KNeighborsClassifier(n_neighbors=K,
                                               weights='distance')
        self.classifier.fit(XT, yT)
        self.labels = self.classifier.classes_
        self.method = str(K)+'-NN'
        # self.train_error(XT, yT)

    def predict(self, Xt):
        """Classify and approximate classes given input vectors.

        Parameters
        ----------
        XT : numpy array
            List of input vectors

        Returns
        -------
        Output space approximation as numpy array

        """
        super().predict(Xt)
        return self.classifier.predict(Xt)

    def predict_prob(self, Xt):
        """Classify and get probability of input vector belonging to any class.

        Parameters
        ----------
        XT : numpy array
            List of input vectors

        Returns
        -------
        Output space approximation as numpy array

        """
        super().predict_prob(Xt)
        return self.classifier.predict_proba(Xt)

    def xtrain(self, XT, yT, **opts):
        """Perform cross validation to find optimal k to use in classification.

        Parameters
        ----------
        XT : numpy array
            List of input vectors
        YT : numpy array
            List of corresponding classes

        """
        nfolds = opts.get('nfolds', 10)
        p_test = opts.get('p_test', 0.1)
        nmax = opts.get('nmax', 20)

        if len(np.unique(yT)) > 5:
            print("Too many classes for xval training. Set K=2.")
            n_opt = 2
        else:
            acc = np.zeros(nmax)
            for n in range(1, nmax + 1):
                # print(f"  - xval for k = {n}")
                for i in range(nfolds):
                    iTrain, itest = xval_indices(
                        XT.shape[0], percent_test=p_test, labels=yT)
                    Xi, yi = XT[iTrain, :], yT[iTrain]
                    ki = KNeighborsClassifier(n_neighbors=n,
                                              weights='distance')
                    ki.fit(Xi, yi)
                    acc[n - 1] += balanced_accuracy_score(
                        yT[itest], ki.predict(XT[itest, :])) / nfolds

            n_opt = np.argmax(acc) + 1

        self.train(XT, yT, K=n_opt)


# MATRIX SCALING


class MatrixScaler:
    """Class used to scale input and output space."""

    def __init__(self, norms, XT):
        """Initialize the Scaler with desired scalings."""
        if len(norms) != XT.shape[1]:
            raise Exception("The no. of columns in XT must be equal "
                            "to the length of norms.")
        self.N = XT.shape[1]
        self.scalers = []
        for i in range(self.N):
            self.scalers.append(DataScaler())
            which = ~np.isnan(XT[:, i])
            self.scalers[i].fit(XT[which, i], method=norms[i])

    def normalize(self, X):
        """Scale input X."""
        assert X.shape[1] == self.N
        Xn = np.empty_like(X)
        for i in range(self.N):
            Xn[:, i] = self.scalers[i].transform(X[:, i])

        return Xn

    def denormalize(self, Xn):
        """Unscale input X."""
        assert Xn.shape[1] == self.N
        X = np.empty_like(Xn)
        for i in range(self.N):
            X[:, i] = self.scalers[i].inv_transform(Xn[:, i])
        return X


# MISCELLANEOUS FUNCTIONS


def xval_indices(N, percent_test=0.15, labels=None):
    """Perform Monte Carlo Cross Validation; stratified if labels are provided.

    Parameters
    ----------
    N : int
        number of samples in the data set.
    percent_test : float
        Percentage of samples to be in the test sets. (0 < percent_test < 0.5)

    Returns
    -------
    (np.ndarray, np.ndarray)
        1D int vectors containing the indices for train and test splits,
        respectively.

    """
    if labels is None:
        indices = np.random.permutation(N)
        NT = int(np.round(N * (1 - percent_test)))
        iT, it = indices[:NT], indices[NT:]
    else:
        assert len(labels) == N
        indices = np.arange(N)
        iT, it = np.empty(0, dtype=int), np.empty(0, dtype=int)
        for label in np.unique(labels):
            ind_l = indices[labels == label]
            NT_l = int(np.round(len(ind_l) * (1 - percent_test)))
            iT = np.hstack((iT, ind_l[:NT_l]))
            it = np.hstack((it, ind_l[NT_l:]))

        np.random.shuffle(iT), np.random.shuffle(it)

    return iT, it


# PERFORMANCE ASSESSMENT


def train_and_assess(grid_train, grid_test, ux2, path, classes):
    """Train interpolators and classes and create assesment plots.

    Parameters
    ----------
    grid_train : string
        Path to training grid
    grid_test : string
        Path to testing grid
    ux2 : list of floats
        List containing the slices at which plots are to be generated. Slices
        are the 3rd variable in the input space, that is the mass of the second
        star for CO-HeMS and CO-HMS_RLO grids and the mass ratio for the
        HMS-HMS grid
    path : string
        The path where assesment plots will be saved
    classes : list of strings
        List containing the classes to be used in the MC-Interpolator

    """
    # check folders exist
    path2obj = os.path.join(path, 'interpolation_objects')
    if not os.path.isdir(path2obj):
        os.mkdir(path2obj)
    path2figs = os.path.join(path, "plots")
    if not os.path.isdir(path2figs):
        os.mkdir(path2figs)
        os.mkdir(os.path.join(path2figs, 'classif'))
        os.mkdir(os.path.join(path2figs, 'interp'))
    else:
        if not os.path.isdir(os.path.join(path2figs, 'classif')):
            os.mkdir(os.path.join(path2figs, 'classif'))
        if not os.path.isdir(os.path.join(path2figs, 'interp')):
            os.mkdir(os.path.join(path2figs, 'interp'))

    grid_T = PSyGrid()
    grid_T.load(grid_train)

    grid_t = PSyGrid()
    grid_t.load(grid_test)

    dat = date.today().strftime("%y%m%d")

    interp_method, class_method = "linear", "kNN"
    mlin = IFInterpolator(grid=grid_T, interp_method=interp_method,
                          class_method=class_method)
    name = "IF_" + dat + '_' + interp_method + "_" + class_method + ".pkl"
    mlin.save(os.path.join(path2obj, name))

    interp_method, class_method = "1NN", "1NN"
    m1NN = IFInterpolator(grid=grid_T, interp_method=interp_method,
                          class_method=class_method)
    name = "IF_" + dat + '_' + interp_method + "_" + class_method + ".pkl"
    m1NN.save(os.path.join(path2obj, name))

    interp_method, class_method = "linear", "kNN"
    # classes = ['no_MT', 'stable_MT', 'unstable_MT']
    mMC = IFInterpolator(grid=grid_T, interp_method=interp_method,
                         class_method=class_method, interp_classes=classes)
    interp_method = "linear3c"
    name = "IF_" + dat + "_" + interp_method + "_" + class_method + ".pkl"
    mMC.save(os.path.join(path2obj, name))

    assess_models([m1NN, mlin, mMC], grid_T, grid_t, ux2, path=path2figs)

    return


def assess_models(models, grid_T, grid_t, ux2, path='./'):
    """Assess models using a grid.

    Parameters
    ----------
    models : psi.Interpolator or list of psi.Interpolator
        Models to be evaluated.
    grid_T : psg.PSyGrid
        Train grid
    grid_t : psg.PSyGrid
        Test grid

    """
    if not isinstance(models, list):
        models = [models]

    path2prmaps = os.path.join(path, 'classif')
    path2interp = os.path.join(path, 'interp')

    # Classes on which the interpolator was trained
    uic = np.unique(
        grid_T.final_values['interpolation_class'][models[0].valid >= 0])

    # Test binaries
    Xt, Yt = models[0]._grid2array(grid_t)
    validt = models[0]._setValid(grid_t.final_values['interpolation_class'],
                                 Xt)
    Ytpred, Ztpred = [], []
    interp_methods, class_methods = [], []
    for m in models:
        ynum, ycat = m.evaluate_mat(Xt[validt >= 0, :])
        Ytpred.append(ynum)
        Ztpred.append(ycat)
        class_methods.append(m.class_method)
        if isinstance(m.interp_method, str):
            interp_methods.append(m.interp_method)
        else:
            interp_methods.append('-'.join(m.interp_method))

    # Assess classification results
    print("\n\n#### Classification ######################################")

    for key in models[0].classifiers:
        print(f"\n   Classification for [{key}]:")
        if models[0].classifiers[key] is not None:
            for c in uic:
                f = grid_T.final_values['interpolation_class'] == c
                print(
                    f"\t{c}: [{np.sum(grid_T.final_values[key][f] == 'None')}"
                    f"/{np.sum(f)}] NaNs")

            # Ground Truth
            zt_gt = grid_t.final_values[key][validt >= 0]
            # Confusion Matrices
            for i, m in enumerate(models):
                print(f"\n\t  {m.classifiers[key].method}")
                bas = balanced_accuracy_score(zt_gt, Ztpred[i][key])
                print(f"\t  Balanced Accuracy Score: {100 * bas:.3f}")
                oa = np.sum(Ztpred[i][key] == zt_gt) / len(zt_gt)
                print(f"\t  Overall Accuracy: {100 * oa:.3f}")
                CM = confusion_matrix(zt_gt, Ztpred[i][key],
                                      m.classifiers[key].labels)
                savename = os.path.join(
                    path2prmaps, key + '_' + class_methods[i] + '.png')
                fig, ax = plot_conf_matrix(
                    100 * CM, m.classifiers[key].method,
                    key, m.classifiers[key].labels, savename=savename)
                plt.close(fig)

                ZT = grid_T.final_values[key][m.valid >= 0]
                zt_gt = grid_t.final_values[key][validt >= 0]
                fig, ax = plot_mc_classifier(
                    m, key, m.XT[m.valid >= 0, :], ZT, ux2, Xt=Xt[validt >= 0],
                    zt=zt_gt, zt_pred=Ztpred[i][key], path=path2prmaps)
        else:
            print("\tNaN")

    # Assess interpolation results
    print("\n\n#### Interpolation ######################################")
    for m in models:
        plot_interpolation(m, ['star_1_mass', 'period_days'], ux2[2:4], ux2,
                           scales=['log', 'log'], path=path2interp)

    Ygt = Yt[validt >= 0, :]
    err_r, err_a = [], []
    for i in range(len(models)):
        err_a.append(np.abs(Ytpred[i] - Ygt))
        err_r.append(err_a[i] / (np.abs(Ygt)))

    # Analize regression errors
    # Nt = np.sum(validt >= 0)
    ic_gt = grid_t.final_values['interpolation_class'][validt >= 0]
    w_interp = ic_gt != 'initial_MT'

    YT = models[0].YT[models[0].valid > 0, :]
    keys = models[0].out_keys

    for ind in range(len(keys)):
        print("\n==================================================")
        print(f"{ind}: {keys[ind]}")
        print("==================================================\n")

        print("                             e_r                    "
              "                             e_a")
        print("        ============================================"
              "        ============================================")
        w_sMT, w_uMT, w_nMT = [], [], []
        for i, m in enumerate(models):
            normi = m.out_scaling[ind]
            ic_t = Ztpred[i]['interpolation_class']
            correct = w_interp == (ic_t != 'initial_MT')
            w_sMT.append(((ic_gt == 'stable_MT') & (ic_t == 'stable_MT')))
            w_uMT.append(((ic_gt == 'unstable_MT') & (ic_t == 'unstable_MT')))
            w_nMT.append(((ic_gt == 'no_MT') & (ic_t == 'no_MT')))

            w = [w_sMT, w_uMT, w_nMT, correct]
            percent = [50, 90]
            tab = np.zeros((len(percent), 2 * len(w)))
            for k in range(4):
                for l in range(len(percent)):
                    tab[l, k] = np.nanpercentile(err_r[i][w[k][i], ind],
                                                 percent[l])
                    tab[l, k + 4] = np.nanpercentile(err_a[i][w[k][i], ind],
                                                     percent[l])

            print(f"classif.: [{class_methods[i]}] --- "
                  f"interp: [{interp_methods[i]}] --- scaling: [{normi}]")
            print("        stable MT  unstable MT    no MT        total"
                  "        stable MT  unstable MT    no MT        total")
            for l in range(len(percent)):
                print(f"  p{str(percent[l])}\t{tab[l, 0]:.2e}\t{tab[l, 1]:.2e}"
                      f"\t{tab[l, 2]:.2e}\t{tab[l, 3]:.2e}"
                      f"\t\t{tab[l, 4]:.2e}\t{tab[l, 5]:.2e}\t{tab[l, 6]:.2e}"
                      f"\t{tab[l, 7]:.2e}")

        print(f"range (training samples): [{YT[:, ind].min()}, "
              f"{YT[:, ind].max()}]")

    plot_interp_error(Ygt, Ytpred, w[0:3], interp_methods, keys,
                      path=path2interp)

    return


def analize_nans(out_keys, YT, valid):
    """Find nans in input numerical variables."""
    print("\nNans in numerical variables:")
    for i, key in enumerate(out_keys):
        n = []
        for v in range(4):
            n.append(np.sum(np.isnan(YT[valid == v, i])))
        if np.sum(np.array(n)) > 0:
            print(f"[i_MT] = {n[0]:5d}, [no_MT] = {n[1]:5d}, "
                  f"[st_MT] = {n[2]:5d}, [u_MT] = {n[3]:5d} : {i:3d} {key}")

    n = []
    for v in range(4):
        n.append(np.sum(valid == v))

    print(f"\n[i_MT] = {n[0]:5d}, [no_MT] = {n[1]:5d}, [st_MT] = {n[2]:5d}, "
          f"[u_MT] = {n[3]:5d}")


def analize_nones(classifiers, valid, grid):
    """Find None values in input categorical variables."""
    print("\nNones in categorical variables:")
    for key in classifiers:
        n = []
        y = grid.final_values[key]
        for v in range(4):
            n.append(np.sum(y[valid == v] == 'None'))
        if np.sum(np.array(n)) > 0:
            print(f"[i_MT] = {n[0]:5d}, [no_MT] = {n[1]:5d}, "
                  f"[st_MT] = {n[2]:5d}, [u_MT] = {n[3]:5d} : {key}")

    n = []
    for v in range(4):
        n.append(np.sum(valid == v))

    print(f"\n[i_MT] = {n[0]:5d}, [no_MT] = {n[1]:5d}, "
          f"[st_MT] = {n[2]:5d}, [u_MT] = {n[3]:5d}")


# PLOTTING FUNCTIONS


def heatmap(CM, ax, **heat_kwargs):
    """Plot confusion matrix as heatmap."""
    lw = heat_kwargs.pop('linewidth', 1.5)
    ax.matshow(CM, **heat_kwargs)
    mi, ma = CM.flatten().min(), CM.flatten().max()
    for i in range(CM.shape[0]):
        for j in range(CM.shape[1]):
            col = 'w'
            if CM[i, j] > (0.7 * (ma - mi) + mi):
                col = 'k'
            ax.text(x=j, y=i, s=str(round(CM[i, j], 2)), va='center',
                    ha='center', size='large', c=col)

    for i in range(CM.shape[0] - 1):
        for j in range(CM.shape[1] - 1):
            ax.plot([i + 0.5] * 2, [-0.5, CM.shape[1] - 0.5], 'w',
                    linewidth=lw)
            ax.plot([-0.5, CM.shape[0] - 0.5], [i + 0.5] * 2, 'w',
                    linewidth=lw)


CONFMAT_KWARGS = {
    'cmap': 'cividis',
    'linewidth': 1.5,
}


def plot_conf_matrix(CM, method, varname, labels, savename=None, **kwargs):
    """Plot confusion matrix for classification assessment."""
    heat_kwargs = {}
    for arg in CONFMAT_KWARGS:
        heat_kwargs[arg] = kwargs.get(arg, CONFMAT_KWARGS[arg])

    fig, ax = plt.subplots(1, 1, figsize=(len(labels),) * 2)
    s = np.sum(CM, axis=1)[:, np.newaxis]
    s[s == 0] = 1

    heatmap(100 * CM / s, ax, **heat_kwargs)

    lmax = 0
    for lab in labels:
        llab = len(lab)
        if llab > lmax:
            lmax = llab

    labels_short = [lab if len(lab) < 12 else lab[:11] for lab in labels]
    alpha = 0
    if lmax > 9:
        alpha = 20

    ax.xaxis.tick_top()
    plt.xticks(np.arange(len(labels)), labels_short, rotation=alpha)
    plt.yticks(np.arange(len(labels)), labels_short, rotation=90 - alpha,
               va='center')

    ax.set(title='predicted class', xlabel=method + ' -- ' + varname,
           ylabel='actual class')

    if lmax > 11:
        handles = [Line2D([0], [0], marker=r'$\mathcal{C}_' + str(i) + '$',
                   color='w', label=str(s[i]) + ' ' + labels[i],
                   markerfacecolor='k', markeredgecolor='None', markersize=11)
                   for i in range(len(labels))]
        plt.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc='upper left',
                   borderaxespad=0.)

    if savename is not None:
        plt.savefig(savename, bbox_inches='tight')

    return fig, ax


PLT_CLASS_KWARGS = {
        'N': 150,
        'figsize': (7, 5),
        'tight_layout': True,
        'dpi': 150,
        's': 10,  # marker size
        'linewidths': 0.05
}


def plot_mc_classifier(m, key, XT, zT, ux2, Xt=None, zt=None, zt_pred=None,
                       path=None, **pltargs):
    """Plot slices illustrating decision boundaries and classification prob."""
    N = pltargs.pop('N', PLT_CLASS_KWARGS['N'])
    marker_size = pltargs.pop('s', PLT_CLASS_KWARGS['s'])
    fig_kwargs = {}
    sct_kwargs = {}
    for arg in PLT_CLASS_KWARGS:
        if arg in ['figsize', 'dpi', 'tight_layout']:
            fig_kwargs[arg] = pltargs.get(arg, PLT_CLASS_KWARGS[arg])
        elif arg in ['s', 'linewidths']:
            sct_kwargs[arg] = pltargs.get(arg, PLT_CLASS_KWARGS[arg])

    classes = m.classifiers[key].labels
    mycolors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
                'tab:brown', 'tab:pink', 'tab:olive', 'tab:cyan', 'tab:gray',
                'lime', 'm', 'aqua', 'bisque', 'salmon',
                'gold', 'palegreen', 'mistyrose', 'yellow', 'darkviolet']
    MAP = ListedColormap(mycolors[:len(classes)])

    for i in range(len(ux2)):

        fig, ax = plt.subplots(**fig_kwargs)
        # Which training and test binaries belong to the ith slice
        whichT = is_in_ball(XT[:, 1], ux2[i],
                            ux2, m.in_scaling[1].startswith('log'))
        xx, yy, X = test_mesh(XT[whichT, :], N)
        proba = (m.test_classifier_prob(key, X).max(axis=1)).reshape(xx.shape)
        zpred = m.test_classifier(key, X).reshape(xx.shape)

        for k, c in enumerate(classes):
            zpred[zpred == c] = k
        _ = ax.contourf(xx, yy, zpred.astype(float), alpha=.6, cmap=MAP)
        contf = ax.contourf(xx, yy, proba, alpha=0.3, cmap='gist_gray')
        fig.colorbar(contf)

        for k, c in enumerate(classes):
            color = [MAP(k)]
            wk = zT[whichT] == c
            ax.scatter(XT[whichT, 0][wk], XT[whichT, 2][wk], c=color,
                       marker='.', **sct_kwargs)

        if Xt is not None:
            whicht = is_in_ball(Xt[:, 1], ux2[i], ux2,
                                m.in_scaling[1].startswith('log'))
            for k, c in enumerate(classes):
                color = [MAP(k)]
                wk = zt[whicht] == c
                ax.scatter(Xt[whicht, 0][wk], Xt[whicht, 2][wk], c=color,
                           marker='*', s=1.5 * marker_size, linewidths=0.4)
            which_miss = whicht & (zt != zt_pred)
            ax.scatter(Xt[which_miss, 0], Xt[which_miss, 2], marker='*',
                       s=1.5 * marker_size, label='missclassified',
                       facecolor='None', edgecolor='k', linewidths=0.4)

        ax.set_xscale('log')
        ax.set_yscale('log')

        if m.interp_in_q:
            var2 = DEFAULT_LABELS['mass_ratio'][0]
        else:
            var2 = DEFAULT_LABELS['star_2_mass'][0]

        ax.set(xlabel=DEFAULT_LABELS['star_1_mass'][1],
               ylabel=DEFAULT_LABELS['period_days'][1])
        ax.set(title=key + ' (' + var2 + ' = ' + str(round(ux2[i], 2)) + ')')

        legend_elements = [
            Line2D([0], [0], marker='.', color='None', label='train binary',
                   markerfacecolor='k', markeredgecolor='k',
                   markersize=0.5 * marker_size, mew=0.4),
            Line2D([0], [0], marker='*', color='None', label='test binary',
                   markerfacecolor='gray', markeredgecolor='None',
                   markersize=0.75 * marker_size, mew=0.4),
            Line2D([0], [0], marker='*', color='None', label='misclassified',
                   markerfacecolor='gray', markeredgecolor='k',
                   markersize=0.75 * marker_size, mew=0.4)]
        leg1 = plt.legend(handles=legend_elements, loc='lower right', ncol=3,
                          fontsize='small')
        legend_el_class = [
            Patch(facecolor=mycolors[j], edgecolor='k', label=classes[j])
            for j in range(len(classes))
        ]
        plt.legend(handles=legend_el_class, loc='upper right',
                   fontsize='small')
        ax.add_artist(leg1)
        if path is not None:
            savename = os.path.join(path,
                                    key + '_' + m.classifiers[key].method
                                    + '_' + str(round(ux2[i], 2)) + '.png')
            fig.savefig(savename)
        plt.close(fig)

    return fig, ax


def is_in_ball(x, x0, ux, log):
    """Helper function for plot_mc_classifier."""
    if log:
        d = np.diff(np.log10(ux))[0] / 2
        w = (np.log10(x) > np.log10(x0) - d) & (np.log10(x) < np.log10(x0) + d)
    else:
        d = np.diff(ux)[0] / 2
        w = (x > x - d) & (x < x + d)
    return w


def test_mesh(XT, N):
    """Helper function for plot_mc_classifier creating mesh for test data."""
    a, b = np.log10(XT[:, 0].min()), np.log10(XT[:, 0].max())
    delta = (b - a) / 50
    x_min, x_max = 10 ** (a - delta), 10 ** (b + delta)
    a, b = np.log10(XT[:, 2].min()), np.log10(XT[:, 2].max())
    delta = (b - a) / 50  # 50 is around 2 times the # of points in the axis
    y_min, y_max = 10 ** (a - delta), 10 ** (b + delta)

    xx, yy = np.meshgrid(np.logspace(np.log10(x_min), np.log10(x_max), N),
                         np.logspace(np.log10(y_min), np.log10(y_max), N))
    X = np.vstack((xx.ravel(), XT[0, 1] * np.ones(N ** 2), yy.ravel())).T

    return xx, yy, X


# Default arguments for matplotlib
PLT_INTERP_KWARGS = {
        'N': 300,
        'figsize': (7, 5),
        'tight_layout': True,
        'dpi': 150,
        's': 10,  # marker size
        'linewidths': 0.05
}


def plot_interpolation(m, keys, v2, ux2, scales=None, path=None, **pltargs):
    """Plot the interpolation errors by key.

    Parameters
    ----------
    m : IFInterpolator
        A trained instance of the IFInterpolator
    keys : List of strings
        Variables for which interpolation errors will be plotted

    """
    N = pltargs.pop('N', PLT_INTERP_KWARGS['N'])
    # marker_size = pltargs.pop('s', PLT_INTERP_KWARGS['s'])
    fig_kwargs = {}
    sct_kwargs = {}
    for arg in PLT_INTERP_KWARGS:
        if arg in ['figsize', 'dpi', 'tight_layout']:
            fig_kwargs[arg] = pltargs.get(arg, PLT_INTERP_KWARGS[arg])
        elif arg in ['s', 'linewidths']:
            sct_kwargs[arg] = pltargs.get(arg, PLT_INTERP_KWARGS[arg])
    keys = list(keys)
    XT = m.XT[m.valid >= 0, :]
    for v in v2:
        # Which training and test binaries belong to the ith slice
        whichT = is_in_ball(XT[:, 1], v, ux2,
                            m.in_scaling[1].startswith('log'))
        xx, yy, X = test_mesh(XT[whichT, :], N)
        ynum, ycat = m.evaluate_mat(X)

        for i, key in enumerate(keys):
            fig, ax = plt.subplots(**fig_kwargs)
            Ct = ynum[:, m.out_keys.index(key)].reshape(xx.shape)
            clab = DEFAULT_LABELS[key][0]
            if scales is not None:
                if scales[i] == 'log':
                    Ct = np.log10(Ct)
                    clab = DEFAULT_LABELS[key][1]

            c = ax.pcolormesh(xx, yy, Ct, cmap='viridis')
            cbar = fig.colorbar(c, ax=ax)
            cbar.set_label(clab)

            ax.set_xscale('log')
            ax.set_yscale('log')

            if m.interp_in_q:
                var2 = DEFAULT_LABELS['mass_ratio'][0]
            else:
                var2 = DEFAULT_LABELS['star_2_mass'][0]

            ax.set(xlabel=DEFAULT_LABELS['star_1_mass'][1],
                   ylabel=DEFAULT_LABELS['period_days'][1])
            ax.set(title=var2 + ' = ' + str(round(v, 2)))

            if path is not None:
                if not isinstance(m.interp_method, list):
                    method = m.interp_method
                else:
                    method = '-'.join(m.interp_method)
                name = (key + '_' + m.classifiers['interpolation_class'].method
                        + '_' + method + '_' + str(round(v, 2)) + '.png')
                name = os.path.join(path, name)
                print(name)
                fig.savefig(name)
            plt.close(fig)


def plot_interp_error(Ygt, Ys, ind_MT, methods, keys, path='.'):
    """Make violin plots for the error distributions of each of the variables.

    Parameters
    ----------
    Ygt : numpy array
        Ground truth of the final values of testing tracks
    Ys : numpy array
        Model's approximation of the final values of testing tracks
    methods : List of strings
        List that indicates the interpolation methods
    keys : List of strings
        List indicating for which variables error distributions will be plotted

    """
    labels = ['stable MT', 'unstable MT', 'no MT']
    col_lab = ['tab:purple', 'tab:olive', 'tab:pink']
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
    errors = ['$e_d (e_a/range)$', '$e_r$', '$e_a$']

    for ind in range(len(keys)):
        for er in errors:
            fig, axs = plt.subplots(1, 2, figsize=(12, 4))

            for j, m in enumerate(methods):
                absdif = np.abs(Ygt[:, ind] - Ys[j][:, ind])
                if er == '$e_a/range$':
                    e_d = absdif / (np.nanmax(Ygt[:, ind])
                                    - np.nanmin(Ygt[:, ind]))
                elif er == '$e_r$':
                    e_d = absdif / (np.abs(Ygt[:, ind]) + 1e-16)
                else:
                    e_d = absdif

                col = colors[j]
                for k in range(len(ind_MT)):
                    if np.sum(ind_MT[k][j]) > 0:
                        parts = axs[0].violinplot(
                            np.log10(e_d[ind_MT[k][j]]),
                            positions=[k], showextrema=True)
                        for pc in parts['bodies']:
                            pc.set_facecolor(col)
                            # pc.set_alpha(0.5)
                        parts['cmins'].set_edgecolor(col)
                        parts['cmaxes'].set_edgecolor(col)
                        parts['cbars'].set_edgecolor(col)
                        axs[0].scatter(k, np.log10(
                            np.nanmedian(e_d[ind_MT[k][j]])),
                                       color=col, marker='o')
                        axs[0].scatter(k, np.log10(
                            np.nanpercentile(e_d[ind_MT[k][j]], 90)),
                                       color=col, marker='x')
                        if j == 2:
                            xx = Ygt[:, ind][ind_MT[k][j]]
                            yy = e_d[ind_MT[k][j]]
                            axs[1].scatter(xx, yy, c=col_lab[k],
                                           label=labels[k], alpha=0.4)
                            axs[1].plot([np.nanmin(xx), np.nanmax(xx)],
                                        [np.nanpercentile(yy, 90)] * 2,
                                        c=col_lab[k])
                            axs[1].set_title(keys[ind])
                            axs[1].set_xlabel('ground truth')
                            axs[1].set_ylabel('error')

            axs[0].set_xticks(np.arange(0, len(labels)))
            axs[0].set_xticklabels(labels)
            axs[0].set_xlim(-0.75, len(labels) - 0.25)
            axs[0].set_title(keys[ind])
            axs[0].set_ylabel(er)

            legend_m = [
                Patch(facecolor=colors[i], edgecolor='k', label=methods[i])
                for i in range(len(methods))]
            axs[0].legend(handles=legend_m, loc='upper right',
                          fontsize='small')
            axs[1].legend(loc='upper right', fontsize='small')
            plt.savefig(os.path.join(
                path, str(ind) + '_' + keys[ind] + '_' + er[1:4] + '.png'))
            plt.close()

    absdif = np.abs(Ygt - Ys[j])
    e_r = absdif / (np.abs(Ygt) + 1e-16)
    order = np.argsort(np.nanpercentile(e_r[ind_MT[0][j], :], 90, axis=0))
    ncol = 10
    nrows = len(order) // ncol
    fig, axs = plt.subplots(nrows, 1, figsize=(12, 60))
    axs[0].set_title('purple: st. MT, green: unst. MT, pink no MT')
    for k in range(len(ind_MT)):
        if np.sum(ind_MT[k][j]) > 0:
            for i in range(nrows):
                indices = order[i * ncol:np.min([(i + 1) * ncol, len(order)])]
                parts = axs[i].violinplot(
                    np.log10(e_r[ind_MT[k][j], :][:, indices]),
                    showextrema=True, showmedians=True)
                for pc in parts['bodies']:
                    pc.set_facecolor(col_lab[k])
                parts['cmins'].set_edgecolor(col_lab[k])
                parts['cmaxes'].set_edgecolor(col_lab[k])
                parts['cbars'].set_edgecolor(col_lab[k])
                parts['cmedians'].set_edgecolor(col_lab[k])
                plt.sca(axs[i])
                plt.xticks(np.arange(1, 11), [str(i) + ' ' + keys[i]
                                              for i in indices], rotation=20)
                axs[i].set_ylim(-7, 1)
                axs[i].set_ylabel(errors[1])
        for i in range(nrows):
            axs[i].grid(axis='y')
    fig.tight_layout()
    plt.savefig(os.path.join(path, '0_' + 'errors.png'))
