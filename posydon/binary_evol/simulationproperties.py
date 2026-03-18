"""Simulation properties for the population class.

This class contains the simulation properties, e.g. flow, steps and `max_time`.
"""


__authors__ = [
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Nam Tran <tranhn03@gmail.com>",
    "Seth Gossage <seth.gossage@northwestern.edu>"
]


import os
import time

from posydon.binary_evol.track_match import TrackMatcher
from posydon.config import PATH_TO_POSYDON_DATA
from posydon.interpolation.interpolation import GRIDInterpolator
from posydon.popsyn.io import simprop_kwargs_from_ini
from posydon.utils.common_functions import convert_metallicity_to_string
from posydon.utils.constants import age_of_universe
from posydon.utils.posydonwarning import Pwarn


class NullStep:
    """An evolution step that does nothing but is used to initialize."""
    pass

class SimulationProperties:
    """Class describing the properties of a population synthesis simulation."""

    # these steps and the flow do not require a metallicity
    ignore_for_met = ["flow", "step_SN", "step_end"]
    steps_with_matching = ["step_detached", "step_CE",
                            "step_isolated", "step_dco",
                            "step_merged", "step_disrupted",
                            "step_initially_single"]

    def __init__(self, flow=({}, {}),
                       step_HMS_HMS = (NullStep(), {}),
                       step_CO_HeMS = (NullStep(), {}),
                       step_CO_HMS_RLO = (NullStep(), {}),
                       step_CO_HeMS_RLO = (NullStep(), {}),
                       step_detached = (NullStep(), {}),
                       step_disrupted = (NullStep(), {}),
                       step_merged = (NullStep(), {}),
                       step_initially_single = (NullStep(), {}),
                       step_dco = (NullStep(), {}),
                       step_SN = (NullStep(), {}),
                       step_CE = (NullStep(), {}),
                       step_end = (NullStep(), {}),
                       extra_hooks = [],
                       **kwargs):

        """Construct the simulation properties object.

        Parameters
        ----------
        flow_chart : dict
            A POSYDON flow_chart dictionary.

        step_HMS_HMS : tuple
            A tuple whose first element is a MesaGridStep class handling
            HMS-HMS evolution, like MS_MS_step. The second element is a
            dictionary of kwargs for that step.

        step_CO_HeMS : tuple
            A tuple whose first element is a MesaGridStep class handling
            CO-HeMS evolution, like CO_HeMS_step. The second element is a
            dictionary of kwargs for that step.

        step_CO_HMS_RLO : tuple
            A tuple whose first element is a MesaGridStep class handling
            CO-HMS-RLO evolution, like CO_HMS_RLO_step. The second element is a
            dictionary of kwargs for that step.

        step_CO_HeMS_RLO : tuple
            A tuple whose first element is a MesaGridStep class handling
            CO-HeMS-RLO evolution, like CO_HeMS_RLO_step. The second element is a
            dictionary of kwargs for that step.

        step_detached : tuple
            A tuple whose first element is a detached_step class handling detached
            evolution. The second element is a dictionary of kwargs for that step.

        step_disrupted : tuple
            A tuple whose first element is a DisruptedStep class handling
            disrupted evolution. The second element is a dictionary of kwargs for
            that step.

        step_merged : tuple
            A tuple whose first element is a MergedStep class handling
            merged evolution. The second element is a dictionary of kwargs for
            that step.

        step_initially_single : tuple
            A tuple whose first element is a InitiallySingleStep class handling
            initially single evolution. The second element is a dictionary of kwargs for
            that step.

        step_dco : tuple
            A tuple whose first element is a DoubleCO class handling
            double CO evolution. The second element is a dictionary of kwargs for
            that step.

        step_SN : tuple
            A tuple whose first element is a StepSN class handling
            supernova evolution. The second element is a dictionary of kwargs for
            that step.

        step_CE : tuple
            A tuple whose first element is a StepCEE class handling
            common envelope evolution. The second element is a dictionary of kwargs for
            that step.

        step_end : tuple
            A tuple whose first element is a step_end class handling
            the end of evolution. The second element is a dictionary of kwargs for
            that step.

        extra_hooks : list of tuples
            Each tuple contains a hooks class and kwargs or the extra step name
            (e.g., 'extra_pre_evolve', 'extra_pre_step', 'extra_post_step',
            'extra_post_evolve') and the corresponding function.
        """

        self.kwargs = kwargs
        verbose = self.kwargs.get('verbose', False)

        # gather kwargs
        self.kwargs["flow"] = flow

        step_kwargs = {"step_HMS_HMS": step_HMS_HMS,
                       "step_CO_HeMS": step_CO_HeMS,
                       "step_CO_HMS_RLO": step_CO_HMS_RLO,
                       "step_CO_HeMS_RLO": step_CO_HeMS_RLO,
                       "step_detached": step_detached,
                       "step_disrupted": step_disrupted,
                       "step_merged": step_merged,
                       "step_initially_single": step_initially_single,
                       "step_dco": step_dco,
                       "step_SN": step_SN,
                       "step_CE": step_CE,
                       "step_end": step_end}

        for key, step_tuple in step_kwargs.items():

            step_class, _ = step_tuple
            if verbose and isinstance(step_class, NullStep):
                Pwarn(f"Step {key} not provided, skipping it.",
                                "StepWarning")

            self.kwargs[key] = step_tuple

        self.kwargs["extra_hooks"] = extra_hooks

        self.default_hooks = EvolveHooks()
        self.all_hooks_classes = [self.default_hooks]
        for item in self.kwargs.get('extra_hooks', []):
            if isinstance(item, tuple):
                if isinstance(item[0], type) and isinstance(item[1], dict):
                    cls, params = item
                    self.all_hooks_classes.append(cls(**params))
                elif isinstance(item[0], str):
                    # setting extra_pre/post_step/evolve methods
                    setattr(self, item[0], item[1])
            else:
                raise ValueError(
                    "`extra_hooks` must be list of tuples with either "
                    "(i) a class deriving from EvolveHooks and a kwargs dict, "
                    "or (ii) the name of the extra function and the callable.")

        # Binary parameters and parameterizations
        self.initial_rotation = 0.0
        self.mass_transfer_efficiency = 1.0

        self.common_envelope_efficiency = 1.0

        # Limits on simulation
        if not hasattr(self, 'max_simulation_time'):
            self.max_simulation_time = age_of_universe
        if not hasattr(self, 'end_events'):
            self.end_events = []
        if not hasattr(self, 'end_states'):
            self.end_states = []

        # for debugging purposes
        if not hasattr(self, 'max_n_steps_per_binary'):
            self.max_n_steps_per_binary = 100

        # Set functions for evolution
        self.all_step_names = []  ## list of strings of all evolutionary steps
        for key, val in self.kwargs.items():
            if "step" not in key:   # skip loading steps
                setattr(self, key, val)
            elif "step" in key:
                self.all_step_names.append(key)
        self.steps_loaded = False

        self.preload_imports()

        # To hold TrackMatcher objects per step, if needed.
        # maybe get rid of this
        self.track_matchers = {}

        # Should possibly be a section in sim props .ini file
        self.grid_path = PATH_TO_POSYDON_DATA
        # These hold GRIDInterpolator objects
        # and associated grid names for ea. metallicity
        # (intended keys are metallicities):
        self.grid_names_Hrich = {}
        self.grids_Hrich = {}
        self.grid_names_strippedHe = {}
        self.grids_strippedHe = {}

    def preload_imports(self):
        """
            Preload the imports of detached_step and MesaGridStep to avoid
        importing them when they are needed when `close()` is called. In
        particular, detached_step imports sklearn, which in turn utilizes
        loky, which invokes its own register.at_exit call. If this happens
        during the `close()` call, which is invoked at shutdown, a
        failure occurs, hence the need for something like this.
        """

        from posydon.binary_evol.CE.step_CEE import StepCEE
        from posydon.binary_evol.DT.step_detached import detached_step
        from posydon.binary_evol.MESA.step_mesa import MesaGridStep

        self._detached_step = detached_step
        self._step_CE = StepCEE
        self._MesaGridStep = MesaGridStep

    @classmethod
    def from_ini(cls, path, metallicity = None, load_steps=False, verbose=False, **override_sim_kwargs):
        """Create a SimulationProperties instance from an inifile.

        Parameters
        ----------
        path : str
            Path to an inifile to load in.

        metallicity : float
            A metallicity (Z) may be provided to automatically assign
            to steps as they are loaded. Should be one of e.g., 2.0, 1.0,
            4.5e-1, 2e-1, 1e-1, 1e-2, 1e-3, 1e-4, corresponding to
            metallicities available in your POSYDON_DATA grids.

        load_steps : bool
            Whether or not evolution steps should be automatically loaded.

        verbose : bool
            Print useful info.

        Returns
        -------
        SimulationProperties
            A new instance of SimulationProperties.
        """

        sim_kwargs = simprop_kwargs_from_ini(path)

        sim_kwargs = {**sim_kwargs, **override_sim_kwargs}

        new_instance = cls(**sim_kwargs)

        if load_steps:
            # Load the steps and required data
            new_instance.load_steps(metallicity=metallicity,
                                    verbose=verbose)

        return new_instance

    def load_steps(self, metallicity=None, verbose=False):
        """Instantiate all step classes and set as instance attributes.

        Parameters
        ----------
        metallicity : float
            A metallicity (Z) may be provided to automatically assign
            to steps as they are loaded. Should be one of e.g., 2.0, 1.0,
            4.5e-1, 2e-1, 1e-1, 1e-2, 1e-3, 1e-4, corresponding to
            metallicities available in your POSYDON_DATA grids.

        verbose : bool
            Print extra information.

        Returns
        -------
        None
        """

        # for every other step, give it a metallicity and load each step
        for name, tup in self.kwargs.items():
            if isinstance(tup, tuple):
                step_kwargs = tup[1]
                metallicity = step_kwargs.get('metallicity', metallicity)
                self.load_a_step(name, tup, metallicity=metallicity, verbose=verbose)

        # track that all steps have been loaded (DEPRECATED?)
        self.steps_loaded = True

    def load_a_step(self, step_name, step_tup=(NullStep, {}), metallicity=None, from_ini='', verbose=False):
        """
        Instantiate and attach a simulation step to this object.

        This method creates an instance of a step class and assigns it as an
        attribute of SimulationProperties using ``step_name`` as the attribute
        name. Step keyword arguments may be provided directly via ``step_tup``
        or loaded from an `.ini` configuration file. Before instantiation,
        step arguments are validated and augmented (e.g., assigning metallicity
        and creating a TrackMatcher if required).

        Parameters
        ----------
        step_name : str
            Name of the evolution step. The created step instance will be
            attached to the object as ``self.<step_name>``. See
            ``SimulationProperties.__init__`` for the standard set of steps.

        step_tup : tuple, optional
            Tuple of the form ``(step_class, kwargs_dict)`` where:

            - ``step_class`` is the class representing the step.
            - ``kwargs_dict`` is a dictionary of keyword arguments used to
            initialize the step.

            Default is ``(NullStep, {})``.

        metallicity : float, optional
            Metallicity (Z) to assign to the step if required and not already
            specified in the step keyword arguments. Default supported values
            are: 2.0, 1.0, 4.5e-1, 2e-1, 1e-1, 1e-2, 1e-3, 1e-4.

        from_ini : str, optional
            Path to an `.ini` file containing step configuration. If provided
            and the file exists, the step class and keyword arguments for
            ``step_name`` are loaded from this file and override ``step_tup``.

        verbose : bool, optional
            If True, print detailed information about step loading and the
            keyword arguments used to instantiate the step.

        Returns
        -------
        None

        Notes
        -----
        - Step keyword arguments are processed by ``self.check_step`` before
        instantiation. This may assign a metallicity and/or attach a
        ``TrackMatcher`` if required for the step.
        - The instantiated step is stored as an attribute of SimulationProperties.
        - After loading, ``self.steps_loaded`` is updated to indicate whether
        all configured steps have been successfully attached.
        """

        if verbose:
            print(f"Loading {step_name}...")

        # grab kwargs from ini file for given step
        if os.path.isfile(from_ini):
            step_tup = simprop_kwargs_from_ini(from_ini, only=step_name)[step_name]

        step_func = step_tup[0]
        step_kwargs = step_tup[1].copy()

        # check to make sure the step has a...
        # 1) metallicity assigned (if needed)
        # 2) TrackMatcher assigned (if needed)
        step_kwargs = self.check_step(metallicity, step_name,
                                      step_kwargs, verbose)

        # Try to load the step
        try:
            setattr(self, step_name, step_func(**step_kwargs))
            if verbose:
                if step_kwargs:
                    print("step_kwargs: ")
                    kw_list = [f"\t{key}: {val}" for key, val in step_kwargs.items()]
                    print("\n".join(kw_list))
                print(f"{step_name} loaded successfully.")
                print()

        except TypeError as e:
            Pwarn(f"Error loading {step_name}: {e}", "StepWarning")
            print(f"Loading {step_name} without arguments.")
            setattr(self, step_name, step_func())

        # check if all steps have been loaded
        self.steps_loaded = all(hasattr(self, name)
                                for name, tup in self.kwargs.items()
                                if isinstance(tup, tuple))

    def check_step(self, metallicity, step_name, step_kwargs, verbose=False):
        """
        Validate and update configuration for an evolution step.

        This method ensures that a valid metallicity is assigned to the step
        (unless the step is excluded from metallicity handling) and that a
        corresponding TrackMatcher exists if the step requires track matching.
        If a TrackMatcher for the `(metallicity, step_name)` combination does
        not yet exist, it is created and stored.

        Parameters
        ----------
        metallicity : float or None
            Default metallicity value to use for the step if not explicitly
            provided in ``step_kwargs``.
        step_name : str
            Name of the pipeline step being checked.
        step_kwargs : dict
            Keyword arguments for the step. This dictionary may be modified
            in-place to include validated metallicity and/or a TrackMatcher
            instance.
        verbose : bool, optional
            If True, print the keyword arguments used to construct the
            TrackMatcher.

        Returns
        -------
        dict
            The updated ``step_kwargs`` dictionary, containing a validated
            ``metallicity`` entry and potentially a ``track_matcher`` object.

        Notes
        -----
        - If metallicity is not provided for a step that requires it, a warning
        is issued and a default value of ``Z = 1.0`` (solar metallicity) is used.
        - TrackMatcher objects are stored in ``self.track_matchers`` and reused
        for repeated `(metallicity, step_name)` combinations.
        """

        # check/assign metallicity for the step
        if step_name not in self.ignore_for_met:
            metallicity = step_kwargs.get('metallicity', metallicity)
            if metallicity is None:
                Pwarn(f"{step_name} not assigned a metallicity. "
                    "Defaulting to Z = Zsun (solar).",
                    "MissingValueWarning")
                metallicity = 1.0
            step_kwargs['metallicity'] = float(metallicity)

        # each metallicity/step combo could require
        # a unique TrackMatcher, so check for that
        matcher_key = (metallicity, step_name)
        if step_name in self.steps_with_matching:
            matcher_needed = matcher_key not in self.track_matchers
            if matcher_needed:
                # create TrackMatcher if needed
                step_kwargs, matcher_kwargs = TrackMatcher.separate_kwargs(step_kwargs)
                self.create_track_matcher(metallicity, step_name, matcher_kwargs)

            if verbose:
                kw_list = [f"\t{key}: {val}" for key, val in matcher_kwargs.items()]
                print(f"matcher_kwargs: \n" + "\n".join(kw_list))
            step_kwargs['track_matcher'] = self.track_matchers[matcher_key]

        return step_kwargs

    def create_track_matcher(self, metallicity, step_name, matcher_kwargs):
        """
        Create and store a TrackMatcher for a given metallicity and step.

        This method ensures that the required stellar evolution grids
        (H-rich and stripped-He) are loaded for the specified metallicity.
        If the corresponding GRIDInterpolator objects do not yet exist,
        they are created and cached. The interpolators are then passed to
        a TrackMatcher instance, which is stored internally.

        Parameters
        ----------
        metallicity : float
            Stellar metallicity used to select the appropriate grid files.
        step_name : str
            Identifier for the evolutionary step associated with this
            TrackMatcher.
        matcher_kwargs : dict
            Keyword arguments used to initialize the TrackMatcher. This
            dictionary will be updated in-place with the following keys:
            'grid_Hrich' and 'grid_strippedHe'.

        Notes
        -----
        - GRIDInterpolator objects are created only once per metallicity
        and reused for subsequent TrackMatcher creations.
        - The created TrackMatcher is stored in ``self.track_matchers``
        using the key ``(metallicity, step_name)``.
        """

        z_str = convert_metallicity_to_string(metallicity)
        # set up GRIDInterpolator objects
        # (only if one hasn't been created already for a given metallicity)
        if metallicity not in self.grids_Hrich:
            self.grid_names_Hrich[metallicity] = os.path.join('single_HMS',
                                                        z_str+'_Zsun.h5')
            grid_path_Hrich = os.path.join(self.grid_path,
                                            self.grid_names_Hrich[metallicity])
            self.grids_Hrich[metallicity] = GRIDInterpolator(grid_path_Hrich)
        if metallicity not in self.grids_strippedHe:
            self.grid_names_strippedHe[metallicity] = os.path.join('single_HeMS',
                                                        z_str+'_Zsun.h5')
            grid_path_strippedHe = os.path.join(self.grid_path,
                                                self.grid_names_strippedHe[metallicity])
            self.grids_strippedHe[metallicity] = GRIDInterpolator(grid_path_strippedHe)

        # Create TrackMatcher object as needed, passing GRIDInterpolator references
        matcher_kwargs['grid_Hrich'] = self.grids_Hrich[metallicity]
        matcher_kwargs['grid_strippedHe'] = self.grids_strippedHe[metallicity]
        self.track_matchers[(metallicity, step_name)] = TrackMatcher(**matcher_kwargs)

    def close(self):
        """Close hdf5 files before exiting."""

        all_step_funcs = [getattr(self, key) for key, val in
                          self.__dict__.items() if 'step_' in key]
        for step_func in all_step_funcs:
            if isinstance(step_func, self._MesaGridStep):
                step_func.close()

        for metallicity in self.grids_Hrich:
            self.grids_Hrich[metallicity].close()
        for metallicity in  self.grids_strippedHe:
            self.grids_strippedHe[metallicity].close()

    def pre_evolve(self, binary):
        """Functions called before a binary evolves.

        Uses all extra hooks classes or extra functions.

        Parameters
        ----------
        binary : instance of <class, BinaryStar>
            The binary before evolution starts.

        Returns
        -------
        binary : instance of <class, BinaryStar>
        """
        for hooks in self.all_hooks_classes:
            hooks.pre_evolve(binary)
        if hasattr(self, 'extra_pre_evolve'):
            self.extra_pre_evolve(binary)
        return binary

    def pre_step(self, binary, step_name):
        """Prepare for step.

        The method is called before every evolution step; uses all extra hooks
        classes or extra functions (except for undefined next step errors).

        Parameters
        ----------
        binary : instance of <class, BinaryStar>
            The binary before evolution starts.
        step_name : str
            The name of the step about to be called (as defined in the flow).

        Returns
        -------
        binary : instance of <class, BinaryStar>

        """
        for hooks in self.all_hooks_classes:
            hooks.pre_step(binary, step_name)
        if hasattr(self, 'extra_pre_step'):
            self.extra_pre_step(binary, step_name)
        return binary

    def post_step(self, binary, step_name):
        """Finalize step.

        The method is called after every evolution step; uses all extra hooks
        classes or extra functions (except for undefined next step errors).

        Parameters
        ----------
        binary : instance of <class, BinaryStar>
            The binary before evolution starts.
        step_name : str
            The name of the step about to be called (as defined in the flow).

        Returns
        -------
        binary : instance of <class, BinaryStar>

        """
        ## do not call extra step hooks if history_verbose=False
        if not binary.history_verbose and binary.event is not None:
            if "redirect" in binary.event:
                return binary

        for hooks in self.all_hooks_classes:
            hooks.post_step(binary, step_name)
        if hasattr(self, 'extra_post_step'):
            self.extra_post_step(binary, step_name)
        return binary

    def post_evolve(self, binary):
        """Finalize the evolution of the binary.

        The method is called after a binary exits the evolution loop.
        Uses all extra hooks classes or extra functions.

        Parameters
        ----------
        binary : instance of <class, BinaryStar>
            The binary after evolution is ended.

        Returns
        -------
        binary : instance of <class, BinaryStar>

        """
        for hooks in self.all_hooks_classes:
            hooks.post_evolve(binary)
        if hasattr(self, 'extra_post_evolve'):
            self.extra_post_evolve(binary)
        return binary


class EvolveHooks:
    """Base class for hooking into binary evolution."""

    def __init__(self):
        """
        Add any new output columns to the hooks constructor.
        Example for extra binary columns:
            self.extra_binary_col_names = ["column_name_1", "column_name_2"]
        Example for extra star columns:
            self.extra_star_col_names = ["column_name_1", "column_name_2"]
        """
        pass

    def pre_evolve(self, binary):
        """Perform actions before a binary evolves."""
        return binary

    def pre_step(self, binary, step_name):
        """Perform actions before every evolution step."""
        return binary

    def post_step(self, binary, step_name):
        """Perform acctions after every evolution step."""
        return binary

    def post_evolve(self, binary):
        """Perform actions after a binary exits the evolution loop."""
        return binary


class TimingHooks(EvolveHooks):
    """Add history column 'step_times' (time taken by step) to each binary.

    Example
    -------
    >>> pop.to_df(extra_columns={'step_times': float})
    """
    def __init__(self):
        self.extra_binary_col_names = ["step_times"]

    def pre_evolve(self, binary):
        """Initialize the step time to match history."""
        if not hasattr(binary, 'step_times'):
            binary.step_times = [0.0]
        return binary

    def pre_step(self, binary, step_name):
        """Record the wall time before taking the step."""
        self.step_start_time = time.time()
        return binary

    def post_step(self, binary, step_name):
        """Record the duration of the step."""

        binary.step_times.append(time.time() - self.step_start_time)

        if len(binary.event_history) > len(binary.step_times):
            diff = len(binary.event_history) - len(binary.step_times)
            binary.step_times += [None] * (diff)
        elif len(binary.event_history) < len(binary.step_times):
            last_items = len(binary.event_history)
            binary.step_times = binary.step_times[-(last_items - 1):]

        return binary

    def post_evolve(self, binary):
        """Add None's to step_times to match history rows."""
        if binary.event == 'END' or binary.event == 'FAILED':
            diff = int(len(binary.event_history) - len(binary.step_times))
            binary.step_times += [None] * diff
        return binary


class StepNamesHooks(EvolveHooks):
    """Add history column 'step_name' to each binary.

    Name of evolutionary step as defined in SimulationProperties.

    >>> pop.to_df(extra_columns={'step_names': str})
    """
    def __init__(self):
        self.extra_binary_col_names = ["step_names"]

    def pre_evolve(self, binary):
        """Initialize the step name to match history."""
        if not hasattr(binary, 'step_names'):
            binary.step_names = ['initial_cond']
        return binary

    def pre_step(self, binary, step_name):
        """Do not do anything before the step."""
        return binary

    def post_step(self, binary, step_name):
        """Record the step name."""

        binary.step_names.append(step_name)
        len_binary_hist = len(binary.event_history)
        len_step_names = len(binary.step_names)
        diff = len_binary_hist - len_step_names

        if len_binary_hist > len_step_names:
            binary.step_names += [None] * (diff)
        elif len_binary_hist < len_step_names:
            binary.step_names = binary.step_names[-(len_binary_hist - 1):]

        return binary

    def post_evolve(self, binary):
        """Ensure None's are append to step_names to match rows in history."""
        if binary.event == 'END' or binary.event == 'FAILED':
            diff = int(len(binary.event_history) - len(binary.step_names))
            binary.step_names += [None]*diff
        return binary


class PrintStepInfoHooks(EvolveHooks):
    """Simple example for adding extra print info."""

    def pre_step(self, binary, step_name):
        """Print the step name for each binary, before taking it."""
        print(binary.index, step_name)
        return binary

    def post_evolve(self, binary):
        """Report at the end of the evolution of each binary."""
        print("End evol for binary {}".format(binary.index), end='\n'*2)
        return binary
