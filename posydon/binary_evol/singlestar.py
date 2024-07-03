"""The star object describe the star current and past state.

The star object contains the current and past states of the star. Only
parameters in the STARPROPERTIES list are stored in the history. The
current parameter value of the star object is accessed as, e.g. `star.mass`,
while his past history with `star.mass_history`.
"""


__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Nam Tran <tranhn03@gmail.com>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Devina Misra <devina.misra@unige.ch>",
]


import numpy as np
import pandas as pd
from posydon.utils.common_functions import check_state_of_star
from posydon.grids.MODELS import MODELS



STARPROPERTIES = [
    'state',            # the evolutionary state of the star. For more info see
                        # `posydon.utils.common_functions.check_state_of_star`
    'metallicity',      # initial mass fraction of metals
    'mass',             # mass (solar units)
    'log_R',            # log10 of radius (solar units)
    'log_L',            # log10 luminosity (solar units)
    'lg_mdot',          # log10 of the absolulte mass change of the star due to
                        # winds and RLO coming from binary history (Msun/yr)
    'lg_system_mdot',   # log10 of the system's absolute mass loss from around
                        # the star due to inneficient mass transfer (Msun/yr)
    'lg_wind_mdot',     # log10 of the absolute wind mass loss of the star from
                        # `lg_wind_mdot_1/2` in binary_history (Msun/yr)
    'he_core_mass',     # in solar units
    'he_core_radius',   # in solar units
    'c_core_mass',      # in solar units
    'c_core_radius',    # in solar units
    'o_core_mass',      # in solar units
    'o_core_radius',    # in solar units
    'co_core_mass',     # in solar units
    'co_core_radius',   # in solar units
    # Central mass fractions
    'center_h1',
    'center_he4',
    'center_c12',
    'center_n14',
    'center_o16',
    'surface_h1',
    # Mass fractions at the surface
    'surface_he4',
    'surface_c12',
    'surface_n14',
    'surface_o16',
    # Energy production from nuclear burning (solar units)
    'log_LH',           # H burning
    'log_LHe',          # He burning
    'log_LZ',           # non-H and non-He burning
    'log_Lnuc',         # total nuclear burning energy production
    'c12_c12',          # Carbon burning through c12-c12
    'center_gamma',
    'avg_c_in_c_core',  # average carbon mass fraction at the carbon core
    'surf_avg_omega',   # surface average rotational angular velocity (rad/sec)
    'surf_avg_omega_div_omega_crit',    # ratio of `surf_avg_omega` to the
                                        # surface critical rotation limit,
                                        # taking into account centrifugal force
                                        # & radiation pressure (Eddington lim.)
    'total_moment_of_inertia',      # total moment of inertia (gr*cm^2)
    'log_total_angular_momentum',   # log10 of total ang. momentum (gr*cm^2/s)
    'spin',                         # the dimesionless spin of the star, if it
                                    # is a compact object, which is equal to
                                    # c*J/(GM^2).
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
    'lambda_CE_pure_He_star_10cent',
    'profile',  # the profile of the star, including extended information of
                # its internal structure, for a specific timestep, usually for
                # the end of the previous step including MESA psygrid.
]

# attributes read from single-star grid runs
STAR_ATTRIBUTES_FROM_STAR_HISTORY_SINGLE = {
    'state': None,                      # to be computed after loading
    'metallicity': None,                # from initial values
    'mass': "star_mass",
    'spin': 'spin_parameter',
    'profile': None
}

def properties_massless_remnant():
    PROPERTIES_MASSLESS = {}
    for key in STARPROPERTIES:
        PROPERTIES_MASSLESS[key] =  np.nan
    PROPERTIES_MASSLESS["state"] = "massless_remnant"
    PROPERTIES_MASSLESS["mass"] = 0.0
    return PROPERTIES_MASSLESS

def convert_star_to_massless_remnant(star):
    for key in STARPROPERTIES:
        setattr(star, key, properties_massless_remnant()[key])
    return star

class SingleStar:
    """Class describing a single star."""

    def __init__(self, **kwargs):
        """Initialize the star.

        Arguments
        ---------
        **kwargs : dict
          List of initialization parameters for a star.
        """
        # Set the initial star properties
        for item in STARPROPERTIES:
            setattr(self, item, kwargs.pop(item, None))
            setattr(self, item + '_history', [getattr(self, item)])
        for item, val in kwargs.items():
            setattr(self, item, val)

        # store extra values in the star object without a history

        # these quantities are updated in step_SN.py
        if not hasattr(self, 'natal_kick_array'):
            self.natal_kick_array = [None] * 4
        if not hasattr(self, 'spin_orbit_tilt_first_SN'):
            self.spin_orbit_tilt_first_SN = None
        if not hasattr(self, 'spin_orbit_tilt_second_SN'):
            self.spin_orbit_tilt_second_SN = None
        if not hasattr(self, 'f_fb'):
            self.f_fb = None
        if not hasattr(self, 'SN_type'):
            self.SN_type = None
        if not hasattr(self, 'm_disk_accreted'):
            self.m_disk_accreted = None
        if not hasattr(self, 'm_disk_radiated'):
            self.m_disk_radiated = None
        if not hasattr(self, 'diff_mass'):
            self.diff_mass = None
        if not hasattr(self, 'diff_log_R'):
            self.diff_log_R = None
        if not hasattr(self, 'diff_center_he4'):
            self.diff_center_he4 = None
        if not hasattr(self, 'diff_surface_he4'):
            self.diff_surface_he4 = None
        if not hasattr(self, 'diff_surface_h1'):
            self.diff_surface_h1 = None
        if not hasattr(self, 'diff_he_core_mass'):
            self.diff_he_core_mass = None
        if not hasattr(self, 'diff_center_c12'):
            self.diff_center_c12 = None
        if not hasattr(self, 'star_state_for_diff_matching'):
            self.star_state_for_diff_matching = None
        if not hasattr(self, 'binary_state_for_diff_matching'):
            self.binary_state_for_diff_matching = None




        # the following quantities are updated in mesa_step.py

        # common envelope quantities
        for quantity in ['m_core_CE', 'r_core_CE']:
            for val in [1, 10, 30, 'pure_He_star_10']:
                if not hasattr(self, f'{quantity}_{val}cent'):
                    setattr(self, f'{quantity}_{val}cent', None)

        # core masses at He depletion
        for quantity in ['avg_c_in_c_core_at_He_depletion',
                         'co_core_mass_at_He_depletion']:
            setattr(self, quantity, None)

        # core collapse quantities
        for MODEL_NAME in MODELS.keys():
            if not hasattr(self, MODEL_NAME):
                setattr(self, MODEL_NAME, None)

    def append_state(self):
        """Append the new version of the star to the end of the star state."""
        for item in STARPROPERTIES:
            getattr(self, item + '_history').append(getattr(self, item))

    def restore(self, i=0, hooks=None):
        """Restore the SingleStar() object to its i-th state, keeping the star history before the i-th state.

        Parameters
        ----------
        i : int
            Index of the star object history to reset the star to. By default
            i == 0, i.e. the star will be restored to its initial state.
        hooks : list
            List of extra hooks associated with the SimulationProperties() of the BinaryStar()
            object containing this SingleStar(), if applicable. This parameter is
            automatically set when restoring a BinaryStar() object.
        """
        if hooks is None:
            hooks = []

        # Move current star properties to the ith step, using its history
        for p in STARPROPERTIES:
            setattr(self, p, getattr(self, '{}_history'.format(p))[i])

            ## delete the star history after the i-th index
            setattr(self, p + '_history', getattr(self, p + '_history')[0:i+1])

        ## if running with extra hooks, restore any extra hook columns
        for hook in hooks:

            if hasattr(hook, 'extra_star_col_names'):
                extra_columns = getattr(hook, 'extra_star_col_names')

                for col in extra_columns:
                    setattr(self, col, getattr(self, col)[0:i+1])


    def to_df(self, **kwargs):
        """Return history parameters from the star in a DataFrame.

        By default all parameters in STARPROPERTIES are included.

        Parameters
        ----------
        extra_columns : dict( 'name':dtype, .... )
            Extra star history parameters to return in DataFrame that are not
            included in STARPROPERTIES. All columns must have an
            associated pandas data type.
            Can be used in combination with `only_select_columns`.
            Assumes names have no suffix.
        ignore_columns : list
            Names of STARPROPERTIES parameters to ignore.
            Assumes names have `_history` suffix.
        only_select_columns : list
            Names of the only columns to include.
            Can be used in combination with `extra_columns`.
            Assumes names have `_history` suffix.
        include_profile : bool
            Include the star's profile in the dataframe (NOT RECOMMENDED)
        null_value : float, optional
            Replace all None values with something else (for saving).
            Default is np.NAN.
        prefix : str, optional
            Prefix to all column names. (e.g. 'star_1', 'S1')
            Default has no prefix.

        Returns
        -------
        pandas DataFrame
        """
        extra_cols_dict = kwargs.get('extra_columns', {})
        extra_columns = list(extra_cols_dict.keys())
        extra_columns_dtypes_user = list(extra_cols_dict.values())

        all_keys = ([key+'_history' for key in STARPROPERTIES]
                    + extra_columns)

        ignore_cols = list(kwargs.get('ignore_columns', []))

        # Do we want to include stellar profiles in output df
        include_profile = kwargs.get('include_profile', False)
        if not include_profile:
            ignore_cols.append('profile')

        keys_to_save = [i for i in all_keys if not
                        ((i.split('_history')[0] in ignore_cols)
                         or (i in ignore_cols))]

        if bool(kwargs.get('only_select_columns')):
            user_keys_to_save = list(kwargs.get('only_select_columns'))
            keys_to_save = ([key+'_history' for key in user_keys_to_save]
                            + extra_columns)
        try:
            # shape of data_to_save (history columns , time steps)
            data_to_save = [getattr(self, key) for key in keys_to_save]


            col_lengths = [len(x) for x in data_to_save]
            max_col_length = np.max(col_lengths)

            # If a singlestar fails, usually history cols have diff lengths.
            # This should append NAN to create even columns.
            all_equal_length_cols = len(set(col_lengths)) == 1
            if not all_equal_length_cols:
                for col in data_to_save:
                    col.extend([np.nan] * abs(max_col_length - len(col)))


            where_none = np.array(
                [[True if var is None else False for var in column]
                 for column in data_to_save], dtype=bool)


        except AttributeError as err:
            raise AttributeError(
                str(err) + "\n\nAvailable attributes in SingleStar: \n{}".
                format(self.__dict__.keys()))

        # casting into object array keeps things general
        star_data = np.array(data_to_save, dtype=object)
        # Convert None to np.NAN by default
        star_data[where_none] = kwargs.get('null_value', np.NAN)
        # sets rows as time steps, columns as history output
        star_data = np.transpose(star_data)

        # add star prefix and remove the '_history' at the end of all column names
        prefix = kwargs.get('prefix', "")
        column_names = [prefix + name.split('_history')[0]
                        for name in keys_to_save]

        star_df = pd.DataFrame(star_data, columns=column_names)

        # try to coerce data types automatically
        star_df = star_df.infer_objects()

        return star_df

    def to_oneline_df(self, history=True, prefix='', **kwargs):
        """Convert SingleStar into a single row DataFrame.

        By default, initial final values of history are used with the `to_df`
        method. Any scalar values can also be added.

        Parameters
        ----------
        scalar_names : list of str
            Names of any values to be added to the oneline DataFrame.
        history : bool
            Include the history initial-final values from to_df method.
        prefix : str
            Any prefix to go at the beginning of all columns names.
        **kwargs
            All options for the to_df method.

        Returns
        -------
        oneline_df
        """
        scalar_names = kwargs.get('scalar_names', [])
        if history:
            output_df = self.to_df(**kwargs)
            # first last row, all cols
            initial_final_data = output_df.values[[0, -1], :]

            oneline_names = (
                [prefix + s + '_i' for s in list(output_df.columns)]
                + [prefix + s + '_f' for s in list(output_df.columns)])
            oneline_data = [d for d in initial_final_data.flatten()]
            oneline_df = pd.DataFrame(data=[oneline_data],
                                      columns=oneline_names)
        else:
            oneline_df = pd.DataFrame()

        for name in scalar_names:
            if hasattr(self, name):
                if name == 'natal_kick_array':
                    natal_kick_array = getattr(self, name)
                    for i in range(4):
                        col_name = prefix+name+'_{}'.format(int(i))
                        oneline_df[col_name] = [natal_kick_array[i]]
                else:
                    oneline_df[prefix+name] = [getattr(self, name)]
        return oneline_df

    def __repr__(self):
        """Return the object representation when print is called."""
        s = '<{}.{} at {}>\n'.format(self.__class__.__module__,
                                     self.__class__.__name__, hex(id(self)))
        for p in STARPROPERTIES:
            s += '{}: {}\n'.format(p, getattr(self, p))

        return s[:-1]

    @staticmethod
    def from_run(run, history=False, profile=False, which_star=None):
        """Create a SingleStar object from a single-star grid run."""
        star = SingleStar()
        star_history = run.history1
        prefix = "S1"
        if star_history is None:
            return star

        n_steps = len(star_history["star_age"]) if history else 1
        h_colnames = star_history.dtype.names

        for attr in STARPROPERTIES:
            # get the corresponding column name (default=same name)
            colname = STAR_ATTRIBUTES_FROM_STAR_HISTORY_SINGLE.get(attr, attr)
            if colname is not None and colname in h_colnames:
                final_value = run.final_values[prefix + "_" + colname]
                if history:
                    col_history = list(star_history[colname])
                else:
                    col_history = [final_value]
            else:
                final_value = None
                col_history = [None] * n_steps
            assert n_steps == len(col_history)
            setattr(star, attr + "_history", col_history)
            setattr(star, attr, final_value)

        try:
            star.metallicity = run.initial_values["Z"]
        except AttributeError:
            star.metallicity = None

        # add values at He depletion
        for colname in run.final_values.dtype.names:
            if "at_He_depletion" in colname:
                if colname[0:3]=="S1_":
                    attr = colname[3:]
                    final_value = run.final_values[colname]
                    setattr(star, attr, final_value)

        star.state_history = [check_state_of_star(star, i=i, star_CO=False)
                              for i in range(n_steps)]
        star.state = star.state_history[-1]

        if profile:
            star.profile_history = [None]*(n_steps-1)+[run.final_profile1]
            star.profile = star.profile_history[-1]

        return star
