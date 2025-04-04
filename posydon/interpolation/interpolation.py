"""Definition of the psyTackInterp class."""

__authors__ = [
    "Juanga Serra Perez <jgserra@northwestern.edu>",
    "Philipp Moura Srivastava <philipp.msrivastava@northwestern.edu>"
]


import pickle
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from .data_scaling import DataScaler
from posydon.grids.psygrid import PSyGrid
from posydon.grids.MODELS import MODELS
from posydon.utils.interpolators import interp1d


class psyTrackInterp:
    """Perform track interpolation for POSYDON."""

    def __init__(self, path, in_keys=None, interp_in_q=False, verbose=False):
        """Initialize the psyTrackInterp class.

        Parameters
        ----------
        path : string
            The path to the training grid for the interpolator
        in_keys : List of strings
            A list indicating the variables used for defining the
            input space, default is None
        interp_in_q : boolean
            Indicates whether or not mass ratio is to be used, default is false
        verbose : boolean
            Indicates whether functionality should be run in verbose mode,
            default is false

        """
        self.grid = PSyGrid()

        # Load grid data
        self.grid.load(path)
        self.grid_path = path

        if in_keys is None:
            self.in_keys = ['star_1_mass', 'star_2_mass', 'period_days']
        else:
            self.in_keys = in_keys
        self.n_in = len(self.in_keys)

        self.interp_in_q = interp_in_q
        self.method = None
        self.verbose = verbose

        self.in_scaling = ['log']*3
        if self.interp_in_q:
            self.in_scaling[1] = 'none'
        self.scalers = None

    def save(self, filename):
        """Save complete interpolation model.

        Parameters
        ----------
        filename : str
            path/name of '.pkl' file where the model will be saved.

        """
        SAVE_ATTRS = self.__dict__
        DONT_SAVE = ['grid']
        myattrs = {key: SAVE_ATTRS[key]
                   for key in SAVE_ATTRS if key not in DONT_SAVE}

        with open(filename, 'wb') as f:
            pickle.dump(myattrs, f)

    def load(self, filename):
        """Load interpolation model to be used for predictions.

        filename : str
            path/name of pickle file to be loaded.

        """
        with open(filename, 'rb') as f:
            myattrs = pickle.load(f)
            for key in myattrs:
                setattr(self, key, myattrs[key])

    def close(self):
        """Close any loaded psygrids."""
        if hasattr(self, 'grid'):
            self.grid.close()

    def train(self, method='NearestNeighbor'):
        """Training method."""
        self.method = method

        XT = self._grid2array(self.grid)
        XTn = self._normalize_fit(XT, self.in_scaling)

        # mask out binaries with any nan value
        # mask = np.logical_and(self.grid.final_values['interpolation_class']
        #                       != 'not_converged',
        # np.invert([pd.isna(XTn)[x,:].any() for x in range(XTn.shape[0])]))
        mask = (self.grid.final_values['interpolation_class']
                != 'not_converged') & pd.notna(np.sum(XT, axis=1))
        self.valid_ind = np.arange(XT.shape[0])[mask]

        if self.method == 'NearestNeighbor':
            self.model = NearestNeighbors(
                n_neighbors=1, algorithm='kd_tree').fit(XTn[mask])
            # I think it should pick KD_tree as inputs are low dimensional
        else:
            self.method = None
            raise ValueError("Method %s not supported." % self.method)

    def evaluate(self, binary, print_dist=False):
        """Evaluate given a binary object."""
        Xt = self._binary2array(binary)
        Xtn = self._normalize(Xt)

        if self.method == 'NearestNeighbor':
            distance, index = self.model.kneighbors(Xtn)
            closest_run = self.grid[self.valid_ind[index[0, 0]]]
            tflags = [closest_run.final_values['interpolation_class'],
                      closest_run.final_values['termination_flag_2'],
                      closest_run.final_values['mt_history'],
                      ]
            dist = np.array([binary.star_1.mass
                             - closest_run.initial_values['star_1_mass'],
                             binary.star_2.mass
                             - closest_run.initial_values['star_2_mass'],
                             binary.orbital_period
                             - closest_run.initial_values['period_days']])
            if print_dist:
                print("Nearest neighbor:")
                print(f"Distance in transformed domain -- d={distance[0,0]}")
                print(f"star_1_mass diff = {dist[0]}")
                print(f"star_2_mass diff = {dist[1]}")
                print(f"orb_period diff = {dist[2]}")

            if self.verbose:
                print('------------------------------------------------------')
                index_in_grid = self.valid_ind[index[0, 0]]
                print('Index to closest MESA binary track in the grid: \n',
                      index_in_grid)
                print('Path to closest MESA binary track: \n',
                      self.grid.MESA_dirs[index_in_grid])
                print('------------------------------------------------------')
        else:
            raise AttributeError("No track interpolation method was trained.")
        return closest_run, dist, tflags

    def _grid2array(self, grid):
        """Convert grid object into data matrix.

        For this, the object will extract the columns indicated by the fields
        given in in_keys and out_keys, stacking them into a the matrix
        `self.XYT`. This matrix will have size N x (self.n_in).

        Parameters
        ----------
        grid : posydon.grids.psygrid.PSyGrid
            Grid of binaries to convert to data matrix.

        Returns
        -------
        numpy.ndarray
            data matrix with self.n_in columns that correspond to the inputs.

        """
        Xn = np.empty((len(grid.initial_values[self.in_keys[0]]), self.n_in))
        for i in range(self.n_in):
            Xn[:, i] = grid.initial_values[self.in_keys[i]]

        if self.interp_in_q:
            i_m1 = self.in_keys.index('star_1_mass')
            i_m2 = self.in_keys.index('star_2_mass')
            Xn[:, i_m2] = Xn[:, i_m2] / Xn[:, i_m1]

        return Xn

    def _binary2array(self, binary):

        var2 = binary.star_2.mass
        if self.interp_in_q:
            var2 = binary.star_2.mass / binary.star_1.mass
        Xt = np.array([[binary.star_1.mass, var2, binary.orbital_period]])

        return Xt

    def _normalize_fit(self, X, norms):
        """Normalize matrix X given the scaling options given by scaling.

        The normalized matrix will be stored as a class attribute.
        The function also stores DataScaler objects that allow for later
        denormalization of variables.

        """
        if not (X.shape[1] == self.n_in) or not (len(norms) == self.n_in):
            raise ValueError("The number of columns in X must match the length of norms.")

        self.scalers = []
        Xn = np.empty_like(X)
        for i in range(self.n_in):
            self.scalers.append(DataScaler())
            Xn[:, i] = self.scalers[i].fit_and_transform(X[:, i],
                                                         method=norms[i])

        return Xn

    def _normalize(self, X):
        """Normalize using learned scalers (_normalize) and apply to matrix X.

        X must have self.n_in columns and be a 2D nparray.

        Parameters
        ----------
        X : numpy.ndarray
            The input matrix to be normalized and of size [N x self.n_in].

        Returns
        -------
        numpy.ndarray
            The normalized test matrix of the same size as X where each column
            has been scaled independently.

        """
        assert (len(X.shape) == 2) and (X.shape[1] == self.n_in)
        Xn = np.empty_like(X)
        for i in range(self.n_in):
            Xn[:, i] = self.scalers[i].transform(X[:, i])
        return Xn


def scaler(x, xmin, xmax):
    """Perform min-max scaling.

    Parameters
    ----------
    x : array_like
        The array that requires min-max rescaling.
    xmin : float
        The lower limit of the range to scale.
    xmax : float
        The upper limit of the range to scale.

    Returns
    -------
    ndarray
        The rescaled array of `x`.

    """
    x = np.asanyarray(x)
    return (x - xmin) / (xmax - xmin)


def inv_scaler(x, xmin, xmax):
    """Restore the original array of min-max rescaled array.

    Parameters
    ----------
    x : array_like
        The array which have been min-max rescaled.
    xmin : float
        The lower limit of the min-max range.
    xmax : float
        The upper limit of the min-max range.

    Returns
    -------
    ndarray
        The original array of `x` that have been min-max rescaled.

    """
    x = np.asanyarray(x)
    return x * (xmax - xmin) + xmin


def fscale(u, lim=None, inv=False):
    """Successively perform min-max rescaling on the last dimension of `u`.

    Parameters
    ----------
    u : array_like
        The array which last dimension requires min-max rescaling.
    lim : sequence of tuple
        The tuples have length 2 and contains the limits on
        the min-max rescaling.
        There should be one tuple for each array in the last dimension.
    inv : bool
        If `False` then the standard min-max rescaling will be performed.
        If `True` then the inverse min-max rescaling
        to restore and array will be performed.

    Returns
    -------
    ndarray
        The rescaled/original array with the same shape as `u`.

    """
    scale = inv_scaler if inv else scaler
    return (np.transpose(
        [scale(u.take(i, -1), *lim[i])
         for i in range(u.shape[-1])]) if lim is not None else u)


class GRIDInterpolator():
    """Class to interpolate between single star MESA tracks.

    Attributes
    ----------
    path : str
        The path to the directory that contains the h5 grid.

    keys : tuple of str
        Contains valid keys for accessing the data.
    """

    def __init__(self, path, verbose=False):
        """Initialize the GRIDInterpolator."""
        self.path = path
        self.verbose = verbose

        grid = PSyGrid()
        self.grid = grid
        grid.load(path)

        grid_mass = []
        for i in range(len(grid)):
            grid_mass.append(grid[i].history1['star_mass'][0])
        self.grid_mass = np.array(grid_mass)

        self.grid_data = dict()
        self.grid_final_values = dict()
        self.grid_profile = dict()

        self.keys = ('age',
                     'mass',
                     'mdot',
                     'conv_mx1_top_r',
                     'conv_mx1_bot_r',
                     'mass_conv_reg_fortides',
                     'thickness_conv_reg_fortides',
                     'radius_conv_reg_fortides',
                     'surface_he3',
                     'he_core_mass',
                     'c_core_mass',
                     'o_core_mass',
                     'he_core_radius',
                     'c_core_radius',
                     'o_core_radius',
                     'center_h1',
                     'center_he4',
                     'center_c12',
                     'center_n14',
                     'center_o16',
                     'surface_h1',
                     'surface_he4',
                     'surface_c12',
                     'surface_n14',
                     'surface_o16',
                     'c12_c12',
                     'center_gamma',
                     'avg_c_in_c_core',
                     'surf_avg_omega',
                     'surf_avg_omega_div_omega_crit',
                     'log_LH',
                     'log_LHe',
                     'log_LZ',
                     'log_Lnuc',
                     'log_Teff',
                     'log_L',
                     'log_R',
                     'inertia',
                     'spin_parameter',
                     'log_total_angular_momentum',
                     'conv_env_top_mass',
                     'conv_env_bot_mass',
                     'conv_env_top_radius',
                     'conv_env_bot_radius',
                     'conv_env_turnover_time_g',
                     'conv_env_turnover_time_l_b',
                     'conv_env_turnover_time_l_t',
                     'envelope_binding_energy',
                     'mass_conv_reg_fortides',
                     'thickness_conv_reg_fortides',
                     'radius_conv_reg_fortides',
                     'lambda_CE_1cent',
                     'lambda_CE_10cent',
                     'lambda_CE_30cent',
                     'co_core_mass',
                     'co_core_radius',
                     'lambda_CE_pure_He_star_10cent',
                     'total_mass_h1',
                     'total_mass_he4')

        self.translate = {
            'age': 'star_age',
            'mass': 'star_mass',
            'he_core_mass': 'he_core_mass',
            'c_core_mass': 'c_core_mass',
            'o_core_mass': 'o_core_mass',
            'he_core_radius': 'he_core_radius',
            'c_core_radius': 'c_core_radius',
            'o_core_radius': 'o_core_radius',
            'center_h1': 'center_h1',
            'center_he4': 'center_he4',
            'center_c12': 'center_c12',
            'center_n14': 'center_n14',
            'center_o16': 'center_o16',
            'surface_h1': 'surface_h1',
            'surface_c12': 'surface_c12',
            'surface_n14': 'surface_n14',
            'surface_o16': 'surface_o16',
            'surface_he4': 'surface_he4',
            'c12_c12': 'c12_c12',
            'avg_c_in_c_core': 'avg_c_in_c_core',
            'center_gamma': 'center_gamma',
            'surf_avg_omega': 'surf_avg_omega',
            'surf_avg_omega_div_omega_crit': 'surf_avg_omega_div_omega_crit',
            'log_LH': 'log_LH',
            'log_LHe': 'log_LHe',
            'log_LZ': 'log_LZ',
            'log_Lnuc': 'log_Lnuc',
            'log_Teff': 'log_Teff',
            'log_L': 'log_L',
            'log_R': 'log_R',
            'inertia': 'total_moment_of_inertia',
            'spin_parameter': 'spin_parameter',
            'log_total_angular_momentum': 'log_total_angular_momentum',
            'conv_env_top_mass': 'conv_env_top_mass',
            'conv_env_bot_mass': 'conv_env_bot_mass',
            'conv_env_top_radius': 'conv_env_top_radius',
            'conv_env_bot_radius': 'conv_env_bot_radius',
            'conv_env_turnover_time_g': 'conv_env_turnover_time_g',
            'conv_env_turnover_time_l_b': 'conv_env_turnover_time_l_b',
            'conv_env_turnover_time_l_t': 'conv_env_turnover_time_l_t',
            'envelope_binding_energy': 'envelope_binding_energy',
            'mass_conv_reg_fortides': 'mass_conv_reg_fortides',
            'thickness_conv_reg_fortides': 'thickness_conv_reg_fortides',
            'radius_conv_reg_fortides': 'radius_conv_reg_fortides',
            'lambda_CE_1cent': 'lambda_CE_1cent',
            'lambda_CE_10cent': 'lambda_CE_10cent',
            'lambda_CE_30cent': 'lambda_CE_30cent',
            'co_core_mass': 'co_core_mass',
            'co_core_radius': 'co_core_radius',
            'lambda_CE_pure_He_star_10cent': 'lambda_CE_pure_He_star_10cent',
            'mdot': 'star_mdot',
            'conv_mx1_top_r': 'conv_mx1_top_r',
            'conv_mx1_bot_r': 'conv_mx1_bot_r',
            'mass_conv_reg_fortides': 'mass_conv_reg_fortides',
            'thickness_conv_reg_fortides': 'thickness_conv_reg_fortides',
            'radius_conv_reg_fortides': 'radius_conv_reg_fortides',
            'surface_he3': 'surface_he3',
            'total_mass_h1': 'total_mass_h1',
            'total_mass_he4': 'total_mass_he4',
            }

        # processed keys
        self.final_keys = (
            'S1_avg_c_in_c_core_at_He_depletion',
            'S1_co_core_mass_at_He_depletion',
            'S1_m_core_CE_1cent',
            'S1_m_core_CE_10cent',
            'S1_m_core_CE_30cent',
            'S1_m_core_CE_pure_He_star_10cent',
            'S1_r_core_CE_1cent',
            'S1_r_core_CE_10cent',
            'S1_r_core_CE_30cent',
            'S1_r_core_CE_pure_He_star_10cent'
        )

        # core collapse keys
        keys = []
        for MODEL_NAME in MODELS.keys():
            for key in ['CO_type', 'SN_type', 'f_fb', 'mass', 'spin',
                        'm_disk_accreted', 'm_disk_radiated','M4', 'mu4',
                        'h1_mass_ej', 'he4_mass_ej']:
                keys.append('S1_' + MODEL_NAME + '_' + key )
        self.final_keys += tuple(keys)

        self.profile_keys = (
            'radius',
            'mass',
            'logRho',
            'omega',
            'energy',
            'x_mass_fraction_H',
            'y_mass_fraction_He',
            'z_mass_fraction_metals',
            'neutral_fraction_H',
            'neutral_fraction_He',
            'avg_charge_He'
        )

    def load_grid(self, *args):
        """Load the requested data to `grid_data`.

        Parameters
        ----------
        *args
            Associated initial masses
            which corresponding data should be loaded.

        """
        grid = self.grid
        for m in args:
            idx = np.argmax(m == self.grid_mass)
            data = grid[idx].history1
            final_value = grid[idx].final_values
            profile = grid[idx].final_profile1
            self.grid_data[m] = dict()
            self.grid_final_values[m] = dict()
            self.grid_profile[m] = dict()
            for key in self.keys:
                self.grid_data[m][key] = data[self.translate[key]]
            for key in self.final_keys:
                self.grid_final_values[m][key] = final_value[key]
            for key in self.profile_keys:
                self.grid_profile[m][key] = profile[key]

    def close(self):
        """Close any loaded psygrids."""
        if hasattr(self, 'grid'):
            self.grid.close()

    def get(self, key, M_new):
        """Perform linear interpolation between specific time-series.

        Parameters
        ----------
        key : str
             The specific time-series described by the key.
             Valid keys are in the `keys` attribute.
        M_new : float
            The associated initial mass
            which time-series requires interpolation.

        Returns
        -------
        ndarray
            The interpolated ZAMS time-series specified by `key`
            associated with the initial mass `M_new`.

        """
        if M_new in self.grid_mass:
            try:
                kvalue = self.grid_data[M_new][key]
            except KeyError:
                self.load_grid(M_new)
                kvalue = self.grid_data[M_new][key]

            log_L = self.grid_data[M_new]['log_L']
            log_LH = self.grid_data[M_new]['log_LH']
            i_zams = np.argmax(np.log10(0.99) + log_L <= log_LH)
        else:
            mass_low = np.max(self.grid_mass[M_new > self.grid_mass])
            mass_high = np.min(self.grid_mass[M_new < self.grid_mass])
            weight = (M_new - mass_low) / (mass_high - mass_low)

            try:
                kvalue_low = self.grid_data[mass_low][key]
            except KeyError:
                self.load_grid(mass_low)
                kvalue_low = self.grid_data[mass_low][key]
            try:
                kvalue_high = self.grid_data[mass_high][key]
            except KeyError:
                self.load_grid(mass_high)
                kvalue_high = self.grid_data[mass_high][key]

            i_cut = min(len(kvalue_low), len(kvalue_high))
            kvalue = kvalue_low[:i_cut] + weight * (kvalue_high[:i_cut]
                                                    - kvalue_low[:i_cut])
            log_L_low = self.grid_data[mass_low]['log_L']
            log_L_high = self.grid_data[mass_high]['log_L']

            log_LH_low = self.grid_data[mass_low]['log_LH']
            log_LH_high = self.grid_data[mass_high]['log_LH']

            log_L = log_L_low[:i_cut] + weight * (log_L_high[:i_cut]
                                                  - log_L_low[:i_cut])
            log_LH = log_LH_low[:i_cut] + weight * (log_LH_high[:i_cut]
                                                    - log_LH_low[:i_cut])
            i_zams = np.argmax(np.log10(0.99) + log_L <= log_LH)

        if key == 'age':
            kvalue = kvalue - kvalue[i_zams]

        return kvalue[i_zams:]

    def get_final_values(self, key, M_new):
        if M_new in self.grid_mass:
            try:
                kvalue = self.grid_final_values[M_new][key]
            except KeyError:
                self.load_grid(M_new)
                kvalue = self.grid_final_values[M_new][key]
        else:
            mass_low = np.max(self.grid_mass[M_new > self.grid_mass])
            mass_high = np.min(self.grid_mass[M_new < self.grid_mass])

            try:
                kvalue_low = self.grid_final_values[mass_low][key]
            except KeyError:
                self.load_grid(mass_low)
                kvalue_low = self.grid_final_values[mass_low][key]

            while pd.isna(kvalue_low):
                # escape if no lower mass is available
                if np.sum(mass_low > self.grid_mass) == 0:
                    break
                mass_low = np.max(self.grid_mass[mass_low > self.grid_mass])
                try:
                    kvalue_low = self.grid_final_values[mass_low][key]
                except KeyError:
                    self.load_grid(mass_low)
                    kvalue_low = self.grid_final_values[mass_low][key]

            try:
                kvalue_high = self.grid_final_values[mass_high][key]
            except KeyError:
                self.load_grid(mass_high)
                kvalue_high = self.grid_final_values[mass_high][key]

            while pd.isna(kvalue_high):
                # escape if no higher mass is available
                if np.sum(mass_high < self.grid_mass) == 0:
                    break
                mass_high = np.min(self.grid_mass[mass_high < self.grid_mass])
                try:
                    kvalue_high = self.grid_final_values[mass_high][key]
                except KeyError:
                    self.load_grid(mass_high)
                    kvalue_high = self.grid_final_values[mass_high][key]

            weight = (M_new - mass_low) / (mass_high - mass_low)

            if ((weight == 0) | ((kvalue_high - kvalue_low) == 0)
                    | np.isinf(weight) | np.isinf((kvalue_high - kvalue_low))):
                kvalue = kvalue_low
            else:
                kvalue = kvalue_low + weight * (kvalue_high - kvalue_low)

        return kvalue

    def get_final_state(self, key, M_new):
        if M_new in self.grid_mass:
            try:
                kstate = self.grid_final_values[M_new][key]
            except KeyError:
                self.load_grid(M_new)
                kstate = self.grid_final_values[M_new][key]
        else:
            mass_low = np.max(self.grid_mass[M_new > self.grid_mass])
            mass_high = np.min(self.grid_mass[M_new < self.grid_mass])

            try:
                kstate_low = self.grid_final_values[mass_low][key]
            except KeyError:
                self.load_grid(mass_low)
                kstate_low = self.grid_final_values[mass_low][key]

            try:
                kstate_high = self.grid_final_values[mass_high][key]
            except KeyError:
                self.load_grid(mass_high)
                kstate_high = self.grid_final_values[mass_high][key]

            mass_center = mass_low + (mass_high - mass_low)/2

            if M_new < mass_center:
                kstate = kstate_low
            else:
                kstate = kstate_high

        return kstate

    def get_profile(self, key, M_new):
        grid = self.grid
        if M_new in self.grid_mass:
            try:
                kvalue = self.grid_profile[M_new][key]
                idx = np.argmax(M_new == self.grid_mass)
                profile_old = grid[idx].final_profile1
            except KeyError:
                self.load_grid(M_new)
                kvalue = self.grid_profile[M_new][key]
                idx = np.argmax(M_new == self.grid_mass)
                profile_old = grid[idx].final_profile1
        else:
            mass_low = np.max(self.grid_mass[M_new > self.grid_mass])
            mass_high = np.min(self.grid_mass[M_new < self.grid_mass])

            try:
                kvalue_low = self.grid_profile[mass_low][key]
            except KeyError:
                self.load_grid(mass_low)
                kvalue_low = self.grid_profile[mass_low][key]
            try:
                kvalue_high = self.grid_profile[mass_high][key]
            except KeyError:
                self.load_grid(mass_high)
                kvalue_high = self.grid_profile[mass_high][key]

            m_cor_low = (self.grid_profile[mass_low]['mass']
                         / self.grid_profile[mass_low]['mass'][0])
            m_cor_high = (self.grid_profile[mass_high]['mass']
                          / self.grid_profile[mass_high]['mass'][0])
            if m_cor_low[-1] < m_cor_high[-1]:
                m_cor = m_cor_high
                idx = np.argmax(mass_high == self.grid_mass)
                profile_old = grid[idx].final_profile1
                f = interp1d(m_cor_low, kvalue_low)
                kvalue_low = f(m_cor)
            else:
                m_cor = m_cor_low
                idx = np.argmax(mass_low == self.grid_mass)
                profile_old = grid[idx].final_profile1
                f = interp1d(m_cor_high, kvalue_high)
                kvalue_high = f(m_cor)

            weight = (M_new - mass_low) / (mass_high - mass_low)
            kvalue = kvalue_low + weight * (kvalue_high - kvalue_low)

        return kvalue, profile_old

    def get_masses_gridfiles(self):
        """Return the masses of the grid files.

        Returns
        -------
        list of floats
        """
        return self.grid_mass
