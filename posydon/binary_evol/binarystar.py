"""The binary object describes current and past state of the binary.

The binary object is composed of two star objects and contains the current and
past states of the binary. Only parameters in the BINARYPROPERTIES list are
stored in the history.

The current parameter value of the star object is accessed as, e.g.
`binary.orbital_period` while its past history with
`binary.orbital_period_history`.

The two stars are accessed as, e.g. `binary.star_1.mass`
while their past history with `binary.star_1.mass_history`.

"""


__authors__ = [
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
    "Nam Tran <tranhn03@gmail.com>",
    "Philipp Moura Srivastava <philipp.msrivastava@gmail.com>",
    "Devina Misra <devina.misra@unige.ch>",
    "Scott Coughlin <scottcoughlin2014@u.northwestern.edu>",
    "Seth Gossage <seth.gossage@northwestern.edu>"
]


import signal
import copy
import numpy as np
import pandas as pd
from posydon.binary_evol.simulationproperties import SimulationProperties
from posydon.binary_evol.singlestar import SingleStar, STARPROPERTIES
from posydon.utils.common_functions import (
    check_state_of_star, orbital_period_from_separation,
    orbital_separation_from_period, get_binary_state_and_event_and_mt_case)
from posydon.popsyn.io import (clean_binary_history_df, clean_binary_oneline_df, 
                               BINARYPROPERTIES_DTYPES)
from posydon.binary_evol.flow_chart import UNDEFINED_STATES
from posydon.utils.posydonerror import FlowError 


# star property: column names in binary history for star 1 and star 2
STAR_ATTRIBUTES_FROM_BINARY_HISTORY = {
    "mass": ["star_1_mass", "star_2_mass"],
    "lg_mdot": ["lg_mstar_dot_1", "lg_mstar_dot_2"],
    "lg_system_mdot": ["lg_system_mdot_1", "lg_system_mdot_2"],
    "lg_wind_mdot": ["lg_wind_mdot_1", "lg_wind_mdot_2"],
}

# only mention those with names different from the column names in history data
STAR_ATTRIBUTES_FROM_STAR_HISTORY = {
    'state': None,                      # to be computed after loading
    'metallicity': None,                # from initial values
    'mass': None,                       # from binary history
    'lg_mdot': None,                    # from binary history
    'lg_system_mdot': None,             # from binary history
    'lg_wind_mdot': None,               # from binary history
    'spin': 'spin_parameter',
    'profile': None
}


BINARY_ATTRIBUTES_FROM_HISTORY = {
    'state': None,
    'event': None,
    'time': 'age',
    'separation': 'binary_separation',
    'orbital_period': 'period_days',
    'eccentricity': None,
    'V_sys': None,
    'mass_transfer_case': None,
    'nearest_neighbour_distance': None
}


BINARYPROPERTIES = [
    # The state and event of the system. For more information, see
    # `posydon.utils.common_functions.get_binary_state_and_event_and_mt_case()
    'state',                    #
    'event',
    'time',                     # age of the system (yr)
    'separation',               # binary orbital separation (solar radii)
    'orbital_period',           # binary orbital period (days)
    'eccentricity',             # binary eccentricity
    'V_sys',                    # list of the 3 systemic velocity coordinates
    # (R_{star} - R_{Roche_lobe}) / R_{Roche_lobe}...
    'rl_relative_overflow_1',   # ...for star 1
    'rl_relative_overflow_2',   # ...for star 2
    'lg_mtransfer_rate',        # log10 of mass lost from the donor (Msun/yr)
    'mass_transfer_case',       # current mass transfer case of the system.
                                # See `get_binary_state_and_event_and_mt_case`
                                # in `posydon.utils.common_functions`.
    'trap_radius',
    'acc_radius',
    't_sync_rad_1',
    't_sync_conv_1',
    't_sync_rad_2',
    't_sync_conv_2',
    'nearest_neighbour_distance',   # the distance of system from its nearest
                                    # neighbour of MESA binary system  in case
                                    # of interpolation during the the end of
                                    # the previous step including MESA psygrid.
                                    # The distance is normalized in the
                                    # parameter space and limits at which it
                                    # was calculated. See `mesa_step` for more.
]


MAXIMUM_STEP_TIME = 120


def signal_handler(signum, frame):
    """React to a maximum time signal."""
    raise RuntimeError("Binary Step Exceeded Alloted Time: {}".
                    format(MAXIMUM_STEP_TIME))


signal.signal(signal.SIGALRM, signal_handler)


class BinaryStar:
    """A class containing the state and history of a stellar binary."""

    def __init__(self, star_1=None, star_2=None, index=None, properties=None,
                 **binary_kwargs):
        """Initialize a binary star.

        Arguments
        ---------
        properties : SimulationProperties
            Instance of the SimulationProperties class (default: None)
        star_1 : SingleStar
            The first star of the binary.
        star_2 : Star
            The second star of the binary.
        **binary_kwargs : dictionary
            List of initialization parameters for a binary
        """
        # Binary Index
        self.index = index

        # Create the two stars
        self.star_1 = star_1 if star_1 is not None else SingleStar()
        self.star_2 = star_2 if star_2 is not None else SingleStar()

        # Set the initial binary properties
        for item in BINARYPROPERTIES:
            if item == 'V_sys':
                setattr(self, item, binary_kwargs.pop(item, np.array([0.0,
                                                                      0.0,
                                                                      0.0])))
            elif item == 'mass_transfer_case':
                setattr(self, item, binary_kwargs.pop(item, 'None'))
            elif item == 'nearest_neighbour_distance':
                setattr(self, item, binary_kwargs.pop(item, ['None',
                                                             'None',
                                                             'None']))
            else:
                setattr(self, item, binary_kwargs.pop(item, None))
            setattr(self, item + '_history', [getattr(self, item)])

        for key, val in binary_kwargs.items():
            setattr(self, key, val)

        if getattr(self.star_1, "mass") is not None and getattr(self.star_2, "mass") is not None:
            if getattr(self, "separation") is None and getattr(self, "orbital_period") is not None:
                setattr(self, "separation", 
                        orbital_separation_from_period(self.orbital_period, self.star_1.mass, self.star_2.mass))
            elif getattr(self, "orbital_period") is None and getattr(self, "separation") is not None:
                setattr(self, "orbital_period", 
                        orbital_period_from_separation(self.separation, self.star_1.mass, self.star_2.mass))

        if not hasattr(self, 'inspiral_time'):
            self.inspiral_time = None
        if not hasattr(self, 'mass_transfer_case'):
            self.mass_transfer_case = 'None'

        if not hasattr(self, 'true_anomaly_first_SN'):
            self.true_anomaly_SN1 = None
        if not hasattr(self, 'true_anomaly_second_SN'):
            self.true_anomaly_SN2 = None
        if not hasattr(self, 'first_SN_already_occurred'):
            self.first_SN_already_occurred = False

        if not hasattr(self, 'history_verbose'):
            self.history_verbose = False

        # store interpolation_class and mt_history for each step_MESA
        for grid_type in ['HMS_HMS','CO_HMS_RLO','CO_HeMS','CO_HeMS_RLO']:
            if not hasattr(self, f'interp_class_{grid_type}'):
                setattr(self, f'interp_class_{grid_type}', None)
            if not hasattr(self, f'mt_history_{grid_type}'):
                setattr(self, f'mt_history_{grid_type}', None)
            if not hasattr(self, f'culmulative_mt_case_{grid_type}'):
                setattr(self, f'culmulative_mt_case_{grid_type}', None)
        
        # SimulationProperties object - parameters & parameterizations
        if isinstance(properties, SimulationProperties):
            self.properties = properties
        else:
            self.properties = SimulationProperties()

    def evolve(self):
        """Evolve a binary from start to finish."""

        self.properties.pre_evolve(self)

        # Code to make sure start time is less than max_simulation_time
        if self.time > self.properties.max_simulation_time:
            raise ValueError(
                "The binary's birth time ({0}) is greater than "
                "`max_simulation_time` ({1}).".format(
                    self.time, self.properties.max_simulation_time))

        max_n_steps = self.properties.max_n_steps_per_binary
        n_steps = 0
    
        while (self.event != 'END' and self.event != 'FAILED'
                and self.event not in self.properties.end_events
                and self.state not in self.properties.end_states):
            
            signal.alarm(MAXIMUM_STEP_TIME)
            self.run_step()

            n_steps += 1
            if max_n_steps is not None:
                if n_steps > max_n_steps:
                    raise RuntimeError("Exceeded maximum number of steps ({})".format(max_n_steps))
        
        signal.alarm(0)     # turning off alarm
        self.properties.post_evolve(self)

    def run_step(self):
        """Evolve the binary through one evolutionary step."""
        
        total_state = (self.star_1.state, self.star_2.state, self.state, self.event)
        if total_state in UNDEFINED_STATES:
            raise FlowError(f"Binary failed with a known undefined state in the flow:\n{total_state}")

        next_step_name = self.properties.flow.get(total_state)
        if next_step_name is None:
            raise ValueError("Undefined next step given stars/binary states {}.".format(total_state))

        next_step = getattr(self.properties, next_step_name, None)
        if next_step is None:
            raise ValueError("Next step name '{}' does not correspond to a function in "
                             "SimulationProperties.".format(next_step_name))

        self.properties.pre_step(self, next_step_name)
        next_step(self)
        
        self.append_state()
        self.properties.post_step(self, next_step_name)

    def append_state(self):
        """Update the history of the binaries' properties."""

        ## do not append redirect steps to the binary history if history_verbose=False
        if not self.history_verbose and getattr(self, "event") is not None:
            if "redirect" in getattr(self, "event"):
                return
        
        # Append to the binary history lists
        for item in BINARYPROPERTIES:
            getattr(self, item + '_history').append(getattr(self, item))

        # Append to the individual star history lists
        self.star_1.append_state()
        self.star_2.append_state()

    def switch_star(self):
        """Switch stars."""
        self.star_1, self.star_2 = self.star_2, self.star_1

    def restore(self, i=0):
        """Restore the BinaryStar() object to its i-th state, keeping the binary history before the i-th state.

        Parameters
        ----------
        i : int
            The index of the binary object history to reset the binary to.
            By default 0, i.e. the star will be restored to its initial state.

        """
        # Move current binary properties to the ith step, using its history
        for p in BINARYPROPERTIES:            
            setattr(self, p, getattr(self, '{}_history'.format(p))[i])

            ## delete the binary history after the i-th index
            setattr(self, p + '_history', getattr(self, p + '_history')[0:i+1])
                       
        ## if running with extra hooks, restore any extra hook columns
        for hook in self.properties.all_hooks_classes:
            
            if hasattr(hook, 'extra_binary_col_names'):
                extra_columns = getattr(hook, 'extra_binary_col_names')

                for col in extra_columns:
                    setattr(self, col, getattr(self, col)[0:i+1])                    
        
        for star in (self.star_1, self.star_2):
            star.restore(i, hooks=self.properties.all_hooks_classes)
                             

    def reset(self, properties=None):
        """Reset the binary to its ZAMS state.

        Parameters
        ----------
        properties : SimulationProperties
            Instance of the SimulationProperties class (default: None)

        """
        # If provided, update the simulation properties class
        if properties is not None:
            self.properties = SimulationProperties(properties)

        # Use the restore function to move the binary back to its initial state
        self.restore(i=0)

    def update_star_states(self):
        """Update the states of the two stars in the binary."""
        if self.star_1.state != 'massless_remnant':
            self.star_1.state = check_state_of_star(
                self.star_1, star_CO=self.star_1.state in ["WD", "NS", "BH"])
        if self.star_2.state != 'massless_remnant':
            self.star_2.state = check_state_of_star(
                self.star_2, star_CO=self.star_2.state in ["WD", "NS", "BH"])

    def to_df(self, **kwargs):
        """Return history parameters from the binary in a DataFrame.

        Includes star 1 and 2 (S1, S2) data and an extra column 'binary_index'.

        Parameters
        ----------
        extra_columns : dict( 'name':dtype, .... )
            Extra binary parameters to return in DataFrame that are not
            included in BINARYPROPERTIES. All columns must have an
            associated pandas data type.
            Can be used in combination with `only_select_columns`.
            Assumes names have no suffix.
        ignore_columns : list
            Names of binary parameters to ignore.
            Assumes names have `_history` suffix.
        only_select_columns : list
            Names of the only columns to include.
            Assumes names have `_history` suffix.
            Can be used in combination with `extra_columns`.
        null_value : float
            Replace all None values with something else (for saving).
            Default is np.nan.
        include_S1, include_S2 : bool
            Choose to include star 1 or 2 data to the DataFrame.
            The default is to include both.
        S1_kwargs, S2_kwargs : dict
            kwargs to pass to each star's 'to_df' method (extra/ignore columns)

        Returns
        -------
        pandas DataFrame

        """
        extra_binary_cols_dict = kwargs.get('extra_columns', {})
        extra_columns = list(extra_binary_cols_dict.keys())

        # dictionary mapping binary properties (plus extras) to their dtypes
        properties_dtypes = {**BINARYPROPERTIES_DTYPES, **extra_binary_cols_dict}

        all_keys = (["binary_index"]
                    + [key+'_history' for key in BINARYPROPERTIES]
                    + extra_columns)

        ignore_cols = list(kwargs.get('ignore_columns', []))
        keys_to_save = [i for i in all_keys if not (
            (i.split('_history')[0] in ignore_cols) or (i in ignore_cols))]

        if bool(kwargs.get('only_select_columns')):
            user_keys_to_save = list(kwargs.get('only_select_columns'))
            keys_to_save = (["binary_index"]
                            + [key+'_history' for key in user_keys_to_save]
                            + extra_columns)


        try:
            data_to_save = [getattr(self, key) for key in keys_to_save[1:]]
            col_lengths = [len(x) for x in data_to_save]
            max_col_length = np.max(col_lengths)

            # binary_index
            data_to_save.insert(0, [self.index]*max_col_length)

            # If a binary fails, usually history cols have diff lengths.
            # This should append NAN to create even columns.
            all_equal_length_cols = len(set(col_lengths)) == 1

            if not all_equal_length_cols:
                for i, col in enumerate(data_to_save):

                    dkey = all_keys[i]
                    try:
                        dtype = properties_dtypes[dkey]
                    except KeyError:
                        # default back to float64 -> np.nan fillers if can not 
                        # determine what the data type should be. This should
                        # never happen.
                        dtype = 'float64'

                    # check type of data and determine what to fill data column with
                    if dtype == 'string':
                        filler_value = 'None'
                    elif dtype == 'float64':
                        filler_value = np.nan
                    # array-like data
                    elif 'object' in dtype:
                        # get the max length that data elements have
                        ele_lengths = [len(x) for x in col]
                        max_ele_length = np.max(ele_lengths)
                        # array of strings
                        if 'string' in dtype:
                            filler_value = np.array([''] * max_ele_length)
                        # array of float64's
                        elif 'float64' in dtype:
                            filler_value = np.array([np.nan] * max_ele_length)

                    # extend data column with determined filler values
                    col.extend([filler_value] * abs(max_col_length - len(col)))

            where_none = np.array([[True if var is None else False
                                    for var in column]
                                   for column in data_to_save], dtype=bool)

        except AttributeError as err:
            raise AttributeError(
                str(err) + "\n\nAvailable attributes in BinaryStar: \n{}".
                format(self.__dict__.keys()))

        # Convert None to np.nan by default
        bin_data = np.array(data_to_save, dtype=object)
        bin_data[where_none] = kwargs.get('null_value', np.nan)

        bin_data = np.transpose(bin_data)

        # remove the _history at the end of all column names
        column_names = [name.split('_history')[0] for name in keys_to_save]
        bin_df = pd.DataFrame(bin_data, columns=column_names)

        # Add 3 columns for V_sys
        if 'V_sys' in column_names:
            V_sys_x = np.zeros(len(bin_df))
            V_sys_y = np.zeros(len(bin_df))
            V_sys_z = np.zeros(len(bin_df))
            for i in range(len(bin_df)):
                V_sys_x[i] = bin_df.iloc[i]['V_sys'][0]
                V_sys_y[i] = bin_df.iloc[i]['V_sys'][1]
                V_sys_z[i] = bin_df.iloc[i]['V_sys'][2]

            bin_df['V_sys_x'] = copy.deepcopy(V_sys_x)
            bin_df['V_sys_y'] = copy.deepcopy(V_sys_y)
            bin_df['V_sys_z'] = copy.deepcopy(V_sys_z)

            # Lose the V_sys list
            bin_df = bin_df.drop(['V_sys'], axis=1)

        frames = [bin_df]
        if kwargs.get('include_S1', True):
            # we are hard coding the prefix
            frames.append(self.star_1.to_df(
                prefix='S1_', null_value=kwargs.get('null_value', np.nan),
                **kwargs.get('S1_kwargs', {})))
        if kwargs.get('include_S2', True):
            frames.append(self.star_2.to_df(
                prefix='S2_', null_value=kwargs.get('null_value', np.nan),
                **kwargs.get('S2_kwargs', {})))
        binary_df = pd.concat(frames, axis=1)

        binary_df.set_index('binary_index', inplace=True)

        extra_s1_cols_dict = kwargs.get('S1_kwargs', {}).get('extra_columns', {})
        extra_s2_cols_dict = kwargs.get('S2_kwargs', {}).get('extra_columns', {})
        binary_df = clean_binary_history_df(binary_df,
                                            extra_binary_dtypes_user=extra_binary_cols_dict,
                                            extra_S1_dtypes_user=extra_s1_cols_dict,
                                            extra_S2_dtypes_user=extra_s2_cols_dict)
        return binary_df

    @classmethod
    def from_df(cls, dataframe, **kwargs):
        """Convert a binary from a pandas DataFrame to BinaryStar instance.

        Parameters
        ----------
        dataframe : Pandas DataFrame
            data to turn into a BinaryStar instance.
        index : int, optional
            Sets the binary index.
        extra_columns : dict, optional
            Column names to be added directly to binary
            not in BINARYPROPERTIES.

        Returns
        -------
            New instance of BinaryStar

        """
        if isinstance(dataframe, pd.Series):
            dataframe = pd.DataFrame(dataframe.to_dict(), index=[0])

        # split input dataframe into kwargs dicts
        binary_params, star1_params, star2_params = dict(), dict(), dict()
        extra_params = dict()
        extra_columns = kwargs.get('extra_columns', {})
        hist_lengths = []
        for name in list(dataframe.columns):
            if 'S1' in name:
                corr_name = name.split('S1_')[-1] + '_history'
                star1_params[corr_name] = list(dataframe[name])
                hist_lengths.append(len(star1_params[corr_name]))
            elif 'S2' in name:
                corr_name = name.split('S2_')[-1] + '_history'
                star2_params[corr_name] = list(dataframe[name])
                hist_lengths.append(len(star2_params[corr_name]))
            elif name in extra_columns:
                # this assumes extra cols in binary, not star 1 or 2
                extra_params[name] = list(dataframe[name])
                hist_lengths.append(len(extra_params[name]))
            else:
                corr_name = name + '_history'
                binary_params[corr_name] = list(dataframe[name])
                hist_lengths.append(len(binary_params[corr_name]))

        # make sure all history columns have equal length
        assert len(set(hist_lengths)) == 1
        history_length = set(hist_lengths).pop()

        if ('binary_index' in dataframe.index.name
                and not kwargs.get('index', False)):
            binary_index = set(dataframe.index).pop()
        else:
            binary_index = kwargs.get('index', None)

        binary = cls(index=binary_index,
                     star_1=SingleStar(**star1_params),
                     star_2=SingleStar(**star2_params),
                     **binary_params)

        # set extra history columns directly
        for key, val in extra_params.items():
            setattr(binary, key, val)

        # set some orbital parameters that should exist by hand
        bp_keys = binary_params.keys()
        if 'eccentricity_history' not in bp_keys:
            setattr(binary, 'eccentricity_history', [0]*history_length)
            setattr(binary, 'eccentricity', 0)
        if ('separation_history' not in bp_keys
                and 'orbital_period_history' in bp_keys):
            separation = orbital_separation_from_period(
                            np.array(binary.orbital_period_history),
                            np.array(binary.star_1.mass_history),
                            np.array(binary.star_2.mass_history),)
            setattr(binary, 'separation_history', list(separation))
            setattr(binary, 'separation', list(separation)[-1])
        if ('orbital_period_history' not in bp_keys
                and 'seperation_history' in bp_keys):
            period = orbital_period_from_separation(
                            np.array(binary.seperation_history),
                            np.array(binary.star_1.mass_history),
                            np.array(binary.star_2.mass_history),)
            setattr(binary, 'orbital_period_history', list(period))
            setattr(binary, 'orbital_period', list(period)[-1])

        # set the binary, star1, star2 parameters to last history value in df
        for params, pointer in zip(
                [star1_params,  star2_params,  binary_params],
                [binary.star_1, binary.star_2, binary]):
            for key, val in params.items():
                setattr(pointer, key.split('_history')[0], val[-1])

        # make BINARYPROPERTIES, history columns same length if not given
        already_included_cols = [name.split('_history')[0]
                                 for name in binary_params.keys()]
        valid_binaryprop_keys = [p for p in BINARYPROPERTIES
                                 if (p not in already_included_cols)]
        for prop_key in valid_binaryprop_keys:
            last_default_value = getattr(binary, prop_key+'_history')[-1]
            len_default = len(getattr(binary, prop_key+'_history'))
            diff = history_length - len_default
            getattr(binary, prop_key+'_history').extend([last_default_value]
                                                        * diff)

        # make STARPROPERTIES, history columns same length if not given
        for params, star in zip(
                [star1_params, star2_params], [binary.star_1, binary.star_2]):
            already_included_cols = [name.split('_history')[0]
                                     for name in params.keys()]
            valid_starprop_keys = [p for p in STARPROPERTIES
                                   if not(p in already_included_cols)]
            for prop_key in valid_starprop_keys:
                last_default_value = getattr(star, prop_key+'_history')[-1]
                len_default = len(getattr(star, prop_key+'_history'))
                diff = history_length - len_default
                getattr(star, prop_key+'_history').extend([last_default_value]
                                                          * diff)

        return binary

    def to_oneline_df(self, scalar_names=[], history=True, **kwargs):
        """Convert binary into a single row DataFrame."""
        if history:
            bin_kwargs = kwargs.copy()
            bin_kwargs['include_S1'] = False
            bin_kwargs['include_S2'] = False
            output_df = self.to_df(**bin_kwargs)
            initial_final_data = output_df.values[[0, -1], :]  # first/last row
            col_names = list(output_df.columns)

            oneline_names = (['binary_index']
                             + [s + '_i' for s in col_names]
                             + [s + '_f' for s in col_names])

            oneline_data = [self.index] + [d for d in
                                           initial_final_data.flatten()]

            bin_df = pd.DataFrame(data=[oneline_data], columns=oneline_names)
        else:
            bin_df = pd.DataFrame()

        s1_kwargs = kwargs.get('S1_kwargs', {})
        if bool(s1_kwargs):
            s1_df = self.star_1.to_oneline_df(prefix='S1_', **s1_kwargs)
        else:
            s1_df = pd.DataFrame()

        s2_kwargs = kwargs.get('S2_kwargs', {})
        if bool(s2_kwargs):
            s2_df = self.star_2.to_oneline_df(prefix='S2_', **s2_kwargs)
        else:
            s2_df = pd.DataFrame()

        oneline_df = pd.concat([bin_df, s1_df, s2_df], axis=1)

        for name in scalar_names:
            if hasattr(self, name):
                oneline_df[name] = [getattr(self, name)]

        # check for variables set in BinaryPopulation handling safe evolution
        if hasattr(self, 'traceback'):
            oneline_df['FAILED'] = [1]
        else:
            oneline_df['FAILED'] = [0]
        if hasattr(self, 'warnings'):
            oneline_df['WARNING'] = [1]
        else:
            oneline_df['WARNING'] = [0]

        oneline_df.set_index('binary_index', inplace=True)

        # try to coerce data types automatically
        oneline_df = oneline_df.infer_objects()

        # Set data types for all columns explicitly
        # we are assuming you may pass the same kwargs to both to_df and oneline
        extra_binary_cols_dict = kwargs.get('extra_columns', {})
        extra_s1_cols_dict = kwargs.get('S1_kwargs', {}).get('extra_columns', {})
        extra_s2_cols_dict = kwargs.get('S2_kwargs', {}).get('extra_columns', {})
        oneline_df = clean_binary_oneline_df(oneline_df,
                                            extra_binary_dtypes_user=extra_binary_cols_dict,
                                            extra_S1_dtypes_user=extra_s1_cols_dict,
                                            extra_S2_dtypes_user=extra_s2_cols_dict)

        return oneline_df

    @classmethod
    def from_oneline_df(cls, oneline_df, **kwargs):
        """Convert a oneline DataFrame into a BinaryStar.

        The oneline DataFrame is expected to have initial-final
        values from history and any individual values that don't have
        histories.

        Parameters
        ----------
        oneline_df : DataFrame
            A oneline DataFrame describing a binary.
        index : int, None
            Binary index
        extra_columns : dict
            Names of any extra history columns not inlcuded
            in BINARYPROPERTIES

        Returns
        -------
            A new BinaryStar instance.

        """
        if isinstance(oneline_df, pd.Series):
            oneline_df = pd.DataFrame(oneline_df.to_dict(), index=[0])

        binary_params, star1_params, star2_params = dict(), dict(), dict()
        extra_params = dict()
        extra_columns = kwargs.get('extra_columns', {})
        hist_lengths = []
        for name in list(oneline_df.columns):
            if '_f' in name[-2:]:
                continue    # ignore final values
            param_name = name.split('_i')[0]

            special_cases = ['natal_kick_array']
            if any([i in param_name for i in special_cases]):
                continue    # deal with special cases later

            # ignore error and warning
            if name in ['FAILED', 'WARNING']:
                continue

            if 'S1' in name:
                param_name = param_name.split('S1_')[-1]
                ending_str = '_history' if param_name in STARPROPERTIES else ''
                star1_params[param_name + ending_str] = list(oneline_df[name])
                hist_lengths.append(len(list(oneline_df[name])))
            elif 'S2' in name:
                param_name = param_name.split('S2_')[-1]
                ending_str = '_history' if param_name in STARPROPERTIES else ''
                star2_params[param_name + ending_str] = list(oneline_df[name])
                hist_lengths.append(len(list(oneline_df[name])))
            elif param_name in extra_columns:
                # this assumes extra cols in binary, not star 1 or 2
                extra_params[param_name] = list(oneline_df[name])
                hist_lengths.append(len(list(oneline_df[name])))
            else:
                # binary
                ending_str = ('_history'
                              if param_name in BINARYPROPERTIES else '')
                binary_params[param_name + ending_str] = list(oneline_df[name])
                hist_lengths.append(len(list(oneline_df[name])))

        # make sure all history columns have equal length
        assert len(set(hist_lengths)) == 1
        history_length = set(hist_lengths).pop()

        if any(['S1_natal_kick_array' in name for name in oneline_df.columns]):
            natalkick_names = ['S1_natal_kick_array_{}'.format(i)
                               for i in range(4)]
            star1_params['natal_kick_array'] = \
                oneline_df[natalkick_names].values

        if any(['S2_natal_kick_array' in name for name in oneline_df.columns]):
            natalkick_names = ['S2_natal_kick_array_{}'.format(i)
                               for i in range(4)]
            star2_params['natal_kick_array'] = \
                oneline_df[natalkick_names].values

        if ('binary_index' in oneline_df.index.name
                and not kwargs.get('index', None)):
            binary_index = set(oneline_df.index).pop()
        else:
            binary_index = kwargs.get('index', None)

        binary = cls(index=binary_index,
                     star_1=SingleStar(**star1_params),
                     star_2=SingleStar(**star2_params),
                     **binary_params)

        # set extra history columns directly
        for key, val in extra_params.items():
            setattr(binary, key, val)

        # set some orbital parameters that should exist by hand
        bp_keys = binary_params.keys()
        if 'eccentricity_history' not in bp_keys:
            setattr(binary, 'eccentricity_history', [0])
            setattr(binary, 'eccentricity', 0)
        if ('separation_history' not in bp_keys
                and 'orbital_period_history' in bp_keys):
            separation = orbital_separation_from_period(
                            np.array(binary.orbital_period_history),
                            np.array(binary.star_1.mass_history),
                            np.array(binary.star_2.mass_history),)
            setattr(binary, 'separation_history', list(separation))
            setattr(binary, 'separation', list(separation)[-1])
        if ('orbital_period_history' not in bp_keys
                and 'seperation_history' in bp_keys):
            period = orbital_period_from_separation(
                            np.array(binary.seperation_history),
                            np.array(binary.star_1.mass_history),
                            np.array(binary.star_2.mass_history),)
            setattr(binary, 'orbital_period_history', list(period))
            setattr(binary, 'orbital_period', list(period)[-1])

        # set the binary, star1, star2 parameters to last history value in df
        for params, pointer in zip(
                [star1_params,  star2_params,  binary_params],
                [binary.star_1, binary.star_2, binary]):
            for key, val in params.items():
                setattr(pointer, key.split('_history')[0], val[-1])

        # make BINARYPROPERTIES, history columns same length if not given
        already_included_cols = [name.split('_history')[0]
                                 for name in binary_params.keys()]
        valid_binaryprop_keys = [p for p in BINARYPROPERTIES
                                 if not(p in already_included_cols)]
        for prop_key in valid_binaryprop_keys:
            last_default_value = getattr(binary, prop_key+'_history')[-1]
            len_default = len(getattr(binary, prop_key+'_history'))
            diff = history_length - len_default
            getattr(binary, prop_key+'_history').extend([last_default_value]
                                                        * diff)

        # if the binary errored, then set the final event
        if bool(oneline_df['FAILED'].values[-1]):
            setattr(binary, 'event', 'FAILED')
        return binary

    def __repr__(self):
        """Return the object representation when print is called."""
        s = '<{}.{} at {}>\n'.format(self.__class__.__module__,
                                     self.__class__.__name__, hex(id(self)))
        if hasattr(self, "error_message"):
            s += "BINARY FAILED: {}\n".format(self.error_message)
        if hasattr(self, "warning_message"):
            s += "WARNING FOUND: {}\n".format(self.warning_message)
        for p in BINARYPROPERTIES:
            s += '{}: {}\n'.format(p, getattr(self, p))
        for star in (self.star_1, self.star_2):
            s += '\n{}\n'.format(star)
        return s[:-1]

    def __str__(self):
        """Get a printable description of the binary star."""
        s = ''
        gap = ', '
        for name in ['state', 'event']:
            s += str(getattr(self, name)) + gap

        def nan_if_not_int_or_float(value):
            """Return nan if `value` is neither int nor float."""
            if isinstance(value, (float, int)):
                return value
            return np.nan

        orb_p = nan_if_not_int_or_float(self.orbital_period)
        m1 = nan_if_not_int_or_float(self.star_1.mass)
        m2 = nan_if_not_int_or_float(self.star_2.mass)

        s += 'p={0:.2f}'.format(orb_p) + gap
        s += 'S1=({0},M={1:.2f})'.format(self.star_1.state, m1) + gap
        s += 'S2=({0},M={1:.2f})'.format(self.star_2.state, m2)
        return 'BinaryStar(' + s + ')'

    @staticmethod
    def from_run(run, history=False, profiles=False):
        """Create a BinaryStar object from a PSyGrid run."""
        binary = BinaryStar()

        # get the data for the binary
        if run.binary_history is not None:
            n_steps = len(run.binary_history["age"]) if history else 1
            bh_colnames = run.binary_history.dtype.names
            for attr in BINARYPROPERTIES:
                colname = BINARY_ATTRIBUTES_FROM_HISTORY.get(attr, attr)
                if colname is not None and colname in bh_colnames:
                    final_value = run.final_values[colname]
                    if history:
                        col_history = list(run.binary_history[colname])
                    else:
                        col_history = [final_value]
                else:
                    final_value = None
                    col_history = [None] * n_steps
                assert n_steps == len(col_history)
                setattr(binary, attr + "_history", col_history)
                setattr(binary, attr, final_value)
        else:
            return binary    # if no binary history, return defaults

        # get the data for each companion star
        for star_history, star, prefix in zip([run.history1, run.history2],
                                              [binary.star_1, binary.star_2],
                                              ["S1", "S2"]):
            if star_history is not None:
                h_colnames = star_history.dtype.names
                for attr in STARPROPERTIES:
                    # get the corresponding column name (default=same name)
                    colname = STAR_ATTRIBUTES_FROM_STAR_HISTORY.get(attr, attr)
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

        # set metallicities (if defined in the track)...
        try:
            metallicity = run.initial_values["Z"]
        except AttributeError:
            metallicity = None
        # ...and other star parameters taken from the binary history
        for star_index, star in enumerate([binary.star_1, binary.star_2]):
            star.metallicity = metallicity
            star.metallicity_history = [metallicity] * n_steps
            for attr, colnames in STAR_ATTRIBUTES_FROM_BINARY_HISTORY.items():
                colname = colnames[star_index]
                final_value = run.final_values[colname]
                if history:
                    col_history = list(run.binary_history[colname])
                else:
                    col_history = [final_value]
                assert n_steps == len(col_history)
                setattr(star, attr, final_value)
                setattr(star, attr + "_history", col_history)

        # add values at He depletion
        for colname in run.final_values.dtype.names:
            if "at_He_depletion" in colname:
                if colname[0:3]=="S1_":
                    attr = colname[3:]
                    final_value = run.final_values[colname]
                    setattr(binary.star_1, attr, final_value)
                elif colname[0:3]=="S2_":
                    attr = colname[3:]
                    final_value = run.final_values[colname]
                    setattr(binary.star_2, attr, final_value)
                else:
                    attr = colname
                    final_value = run.final_values[colname]
                    setattr(binary, attr, final_value)
        
        # update eccentricity
        binary.eccentricity = 0.0
        binary.eccentricity_history = [0.0] * n_steps

        # update star states
        n_history = len(binary.time_history)
        for star, track in zip([binary.star_1, binary.star_2],
                               [run.history1, run.history2]):
            is_CO = track is None
            state_history = [check_state_of_star(star, i=i, star_CO=is_CO)
                             for i in range(n_history)]
            star.state_history = state_history
            star.state = state_history[-1]

        # update binary state, event and MT case
        binary.state_history = []
        binary.event_history = []
        binary.mass_transfer_case_history = []
        for i in range(n_history):  # step-by-step: previous states matter!
            result = get_binary_state_and_event_and_mt_case(binary, i=i)
            binary.state, binary.event, binary.mass_transfer_case = result
            binary.state_history.append(binary.state)
            binary.event_history.append(binary.event)
            binary.mass_transfer_case_history.append(binary.mass_transfer_case)

        if profiles:
            binary.star_1.profile = run.final_profile1
            binary.star_2.profile = run.final_profile2
            
        return binary
    
    def initial_condition_message(self, ini_params=None):
        """Generate a message with the initial conditions.

        Parameters
        ----------
       
        ini_params : None or iterable of str
            If None take the initial conditions from the binary, otherwise add
            each item of it to the message.

        Returns
        -------
        string
            The message with the binary initial conditions.
        """
    
        if ini_params is None:
            ini_params = ["\nFailed Binary Initial Conditions:\n",
                    f"S1 mass: {self.star_1.mass_history[0]} \n",
                    f"S2 mass: {self.star_2.mass_history[0]} \n",
                    f"S1 state: {self.star_1.state_history[0]} \n",
                    f"S2 state: {self.star_2.state_history[0]}\n",
                    f"orbital period: {self.orbital_period_history[0] } \n",
                    f"eccentricity: {self.eccentricity_history[0]} \n",
                    f"binary state: {self.state_history[0] }\n",
                    f"binary event: {self.state_history[0] }\n",
                    f"S1 natal kick array: {self.star_1.natal_kick_array }\n",
                    f"S2 natal kick array: {self.star_2.natal_kick_array}\n"]
            
        message = ""
        for i in ini_params:
            message += i 

        return message

