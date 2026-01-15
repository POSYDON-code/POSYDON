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
    "Seth Gossage <seth.gossage@northwestern.edu>"
]


import numpy as np
import pandas as pd

from posydon.grids.SN_MODELS import SN_MODELS
from posydon.popsyn.io import (
    EXTRA_STAR_COLUMNS_DTYPES,
    SCALAR_NAMES_DTYPES,
    STARPROPERTIES_DTYPES,
)
from posydon.utils.common_functions import (
    CO_radius,
    check_state_of_star,
    infer_star_state,
)
from posydon.utils.constants import Zsun, zams_table
from posydon.utils.limits_thresholds import THRESHOLD_CENTRAL_ABUNDANCE
from posydon.utils.posydonwarning import Pwarn

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
    # Mass fractions at the surface
    'surface_h1',
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
    'total_mass_h1',   # total mass of Hydrogen throughout the star
    'total_mass_he4',  # total mass of Helium throughout the star
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

        # Initialize composition and radius at least.
        # This is needed to match to a star and for
        # evolution to be possible. Ideally the matcher
        # would check things like burning luminosity and more
        # elements beyong He. So, in the future we might
        # consider making more robust matching criteria.

        # (This only comes into play if a user doesn't set these)
        setattr(self, 'metallicity', kwargs.pop('metallicity', 1.0))
        state = kwargs.get('state', 'H-rich_Core_H_burning')
        CO_states = ['massless_remnant', 'WD', 'NS', 'BH']

        if "natal_kick_array" in kwargs:
            tmp = kwargs['natal_kick_array']
            kwargs['natal_kick_velocity'] = tmp[0]
            kwargs['natal_kick_azimuthal_angle'] = tmp[1]
            kwargs['natal_kick_polar_angle'] = tmp[2]
            kwargs['natal_kick_mean_anomaly'] = tmp[3]
            del kwargs['natal_kick_array']

        if state in CO_states:
            Z_div_Zsun = self.metallicity
            Y = np.nan
            Z = np.nan
            X = np.nan
            LOW_ABUNDANCE = np.nan
            # a low/high value to guess with, seems to work well
            LOW_LOGR_GUESS = np.nan
            HIGH_LOGR_GUESS = np.nan

            # start by guessing a smallish radius and no He core
            default_log_R = np.nan
            default_He_core_mass = np.nan
        else:
            Z_div_Zsun = self.metallicity
            if Z_div_Zsun in zams_table.keys():
                Y = zams_table[Z_div_Zsun]
            else:
                raise KeyError(f"{Z_div_Zsun} is a not defined metallicity")
            Z = Z_div_Zsun*Zsun
            X = 1.0 - Y - Z
            LOW_ABUNDANCE = 1e-6
            # a low/high value to guess with, seems to work well
            LOW_LOGR_GUESS = 0.0
            HIGH_LOGR_GUESS = 4.0

            # start by guessing a smallish radius and no He core
            default_log_R = LOW_LOGR_GUESS
            default_He_core_mass = 0.0

        # MAIN SEQUENCE
        if "Core_H_burning" in state:
            # default HMS ZAMS
            default_core_X = X
            default_core_Y = Y
        # POST-MS
        elif "burning" in state and ("_Shell_" in state or "_non_" in state):
            # default TAMS, either accreted He or H-rich, low core X, high Y
            default_core_X = THRESHOLD_CENTRAL_ABUNDANCE
            default_core_Y = 1.0 - default_core_X - Z
            if "stripped_He" in state:
                default_core_X = LOW_ABUNDANCE
                default_core_Y = LOW_ABUNDANCE
                default_He_core_mass = kwargs.get('mass')

        # ADVANCED BURNING
        elif "Core_He_burning" in state:
            # default to start of He burn
            default_core_X = LOW_ABUNDANCE
            default_core_Y = 1.0 - default_core_X - Z
            # make radius big to encourage match to giant
            default_log_R = HIGH_LOGR_GUESS
            if 'stripped_He' in state:
                # unless it is a stripped He star
                default_log_R = LOW_LOGR_GUESS
            # large He core mass to encourage maximal He core size
            default_He_core_mass = kwargs.get('mass')
        elif "Core_C_burning" in state:
            # default to start of He burn
            default_core_X = LOW_ABUNDANCE
            default_core_Y = LOW_ABUNDANCE
            # make radius big to encourage match to giant
            default_log_R = HIGH_LOGR_GUESS
            # This stays big after He core forms
            default_He_core_mass = kwargs.get('mass')
        elif ("Core_" in state) and ("_depleted" in state):
            # core He or heavier depleted
            # default to end of He burn. Matching does not check
            # other elements havier than He, so this is best
            # we can do.
            default_core_X = LOW_ABUNDANCE
            default_core_Y = LOW_ABUNDANCE
            # This stays big after He core forms
            default_He_core_mass = kwargs.get('mass')
        # COMPACT OBJECT
        elif state in ['WD', 'NS', 'BH']:
            default_core_X = np.nan
            default_core_Y = np.nan
            default_log_R = np.log10(CO_radius(kwargs.get('mass'), state))
            default_He_core_mass = np.nan
            # If a user gives a mass that does not comply with our
            # CO star state logic, you can get weird stuff like a
            # BH or NS turning into a WD.
            inferred_state = infer_star_state(star_mass=kwargs.get('mass'),
                             surface_h1=kwargs.get('surface_h1', np.nan),
                             center_h1=kwargs.get('center_h1', np.nan),
                             center_he4=kwargs.get('center_he4', np.nan),
                             center_c12=kwargs.get('center_c12', np.nan),
                             star_CO=True)
            if state != inferred_state:
                kwargs['state'] = inferred_state
        elif state == 'massless_remnant':
            default_core_X = np.nan
            default_core_Y = np.nan
            default_log_R = np.nan
            default_He_core_mass = np.nan
        else:
            # some state not caught above, default HMS ZAMS
            Pwarn(f"The initial state {state} was not caught in "
                  "SingleStar.__init__() so it was initialized " \
                  "as an H-rich_Core_H_burning star.",
                  "InitializationWarning")
            default_core_X = X
            default_core_Y = Y


        # Done setting initial values, assign everything to SingleStar next

        for item in STARPROPERTIES:

            # set default values when a kwarg is absent
            # matching at least needs these initiailized to non-null values
            if item == 'log_R':
                # initiate log radius
                setattr(self, item, kwargs.pop(item, default_log_R))
            elif item == 'center_h1':
                # initialize surface and center h1
                setattr(self, item, kwargs.pop(item, default_core_X))
            elif item == 'center_he4':
                # initialize surface and center he4
                setattr(self, item, kwargs.pop(item, default_core_Y))
            elif 'core_mass' in item:
                # intiailize all core_mass values to 0
                # they are used in matching, but we will rely on
                # abundances set above to find a good match
                if item == 'he_core_mass':
                    default_core_mass = default_He_core_mass
                else:
                    if state in CO_states:
                        default_core_mass = 0.0
                    else:
                        default_core_mass = np.nan

                setattr(self, item, kwargs.pop(item, default_core_mass))
            elif item == 'metallicity':
                # already set
                pass

            # everything else, set it according to stored dtype:
            else:
                dtype = STARPROPERTIES_DTYPES.get(item, '')
                if dtype == 'float64' or dtype == 'int':
                    default = np.nan
                elif dtype == 'string':
                    default = ''
                else:
                    # if we haven't defined a dtype for this prop
                    default = None
                setattr(self, item, kwargs.pop(item, default))

            setattr(self, item + '_history', [getattr(self, item)])

        for key, val in kwargs.items():
            setattr(self, key, val)
            # do not create history for scalar values
            if key in SCALAR_NAMES_DTYPES.keys() or key in EXTRA_STAR_COLUMNS_DTYPES.keys():
                continue
            setattr(self, key + '_history', [val])

        # store extra values in the star object without a history

        # these quantities are updated in step_SN.py
        if not hasattr(self, 'natal_kick_velocity'):
            self.natal_kick_velocity = None
        if not hasattr(self, 'natal_kick_azimuthal_angle'):
            self.natal_kick_azimuthal_angle = None
        if not hasattr(self, 'natal_kick_polar_angle'):
            self.natal_kick_polar_angle = None
        if not hasattr(self, 'natal_kick_mean_anomaly'):
            self.natal_kick_mean_anomaly = None
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
        if not hasattr(self, 'h1_mass_ej'):
            self.h1_mass_ej = None
        if not hasattr(self, 'he4_mass_ej'):
            self.he4_mass_ej = None
        if not hasattr(self, 'M4'):
            self.M4 = None
        if not hasattr(self, 'mu4'):
            self.mu4 = None
        if not hasattr(self, 'Xi'):
            self.Xi = None
        if not hasattr(self, 'sc'):
            self.sc = None
        if not hasattr(self, 'interp1d'):
            self.interp1d = None

        # the following quantities are updated in mesa_step.py

        # common envelope quantities
        for quantity in ['m_core_CE', 'r_core_CE']:
            for val in [1, 10, 30, 'pure_He_star_10']:
                if not hasattr(self, f'{quantity}_{val}cent'):
                    setattr(self, f'{quantity}_{val}cent', None)

        # core masses at He depletion
        for quantity in ['avg_c_in_c_core_at_He_depletion',
                         'co_core_mass_at_He_depletion']:
            if not hasattr(self, quantity):
                setattr(self, quantity, None)

        # core collapse quantities
        for SN_MODEL_NAME in SN_MODELS.keys():
            if not hasattr(self, SN_MODEL_NAME):
                setattr(self, SN_MODEL_NAME, None)

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
            Default is np.nan.
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
        # Convert None to np.nan by default
        star_data[where_none] = kwargs.get('null_value', np.nan)
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
                # Handle legacy natal_kick_array for backward compatibility
                if name == 'natal_kick_array':
                    Pwarn("The 'natal_kick_array' attribute will be deprecated. "
                            "Please use 'natal_kick_velocity', "
                            "'natal_kick_azimuthal_angle', "
                            "'natal_kick_polar_angle', and "
                            "'natal_kick_mean_anomaly' instead. "
                            "Adding both properties to the DataFrame.",
                            "DeprecationWarning"
                    )
                    # Create array from individual properties
                    natal_kick_array = [
                        getattr(self, 'natal_kick_velocity', None),
                        getattr(self, 'natal_kick_azimuthal_angle', None),
                        getattr(self, 'natal_kick_polar_angle', None),
                        getattr(self, 'natal_kick_mean_anomaly', None)
                    ]
                    for i in range(4):
                        col_name = prefix+name+'_{}'.format(int(i))
                        oneline_df[col_name] = natal_kick_array[i]

                    # also output better named columns
                    oneline_df[prefix+'natal_kick_velocity'] = natal_kick_array[0]
                    oneline_df[prefix+'natal_kick_azimuthal_angle'] = natal_kick_array[1]
                    oneline_df[prefix+'natal_kick_polar_angle'] = natal_kick_array[2]
                    oneline_df[prefix+'natal_kick_mean_anomaly'] = natal_kick_array[3]
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
