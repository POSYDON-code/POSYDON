"""Simulation properties for the population class.

This class contains the simulation properties, e.g. flow, steps and `max_time`.
"""


__authors__ = [
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Nam Tran <tranhn03@gmail.com>",
]


import time
from posydon.utils.constants import age_of_universe
from posydon.utils.x_ray_luminosities import x_ray_luminosity

class SimulationProperties:
    """Class describing the properties of a population synthesis simulation."""

    def __init__(self, properties=None, **kwargs):
        """Construct the simulation properties object.

        Parameters
        ----------
        properties : object
            Simulation Properties class containing, e.g. flow, steps.
        extra_hooks : list of tuples
            Each tuple contains a hooks class and kwargs or the extra step name
            ('extra_pre_evolve', 'extra_pre_step', 'extra_post_step',
            'extra_post_evolve') and the corresponding function.
        """
        self.kwargs = kwargs.copy()
        # Check if binary_properties is passed

        self.default_hooks = EvolveHooks()
        self.all_hooks_classes = [self.default_hooks]
        for item in kwargs.get('extra_hooks', []):
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
        #if not hasattr(self, "verbose_binary_errors"):
        #    self.verbose_binary_errors = False

        # Set functions for evolution
        for key, val in kwargs.items():
            if "step" not in key:   # skip loading steps
                setattr(self, key, val)
        self.steps_loaded = False

    def load_steps(self, verbose=False):
        """Instantiate step classes and set as instance attributes.

        Parameters
        ----------
        verbose: bool
            Print extra information.

        Returns
        -------
        None
        """
        if verbose:
            print('STEP NAME'.ljust(20) + 'STEP FUNCTION'.ljust(25) + 'KWARGS')
        for name, tup in self.kwargs.items():
            if isinstance(tup, tuple):
                if verbose:
                    print(name, tup, end='\n')
                step_func, kwargs = tup
                setattr(self, name, step_func(**kwargs))
        self.steps_loaded = True

    def close(self):
        """Close hdf5 files before exiting."""
        from posydon.binary_evol.MESA.step_mesa import MesaGridStep
        from posydon.binary_evol.DT.step_detached import detached_step
        all_step_funcs = [getattr(self, key) for key, val in
                          self.__dict__.items() if 'step_' in key]
        for step_func in all_step_funcs:
            if isinstance(step_func, MesaGridStep):
                step_func.close()
            elif isinstance(step_func, detached_step):
                for grid_interpolator in [step_func.grid_Hrich, step_func.grid_strippedHe]:
                    grid_interpolator.close()

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


class XrbLuminositiesHooks(EvolveHooks):
    """Idetifies the X-ray binaries and calculates their X-ray luminosities
    """
    def pre_evolve(self, binary):
        """Perform actions before a binary evolves."""
        binary.xrb_luminosity = []
        binary.xrb_beaming=[]
        return binary

    def pre_step(self, binary, step_name):
        """Perform actions before every evolution step."""
        return binary

    def post_step(self, binary, step_name):
        """Perform acctions after every evolution step."""
        # print("in extra hooks, added XRB Lx and beaming = ", x_ray_luminosity(binary))
        binary.xrb_luminosity.append(x_ray_luminosity(binary)[0])
        binary.xrb_beaming.append(x_ray_luminosity(binary)[1])

        return binary

    def post_evolve(self, binary):
        """Ensure None's are append to match rows in history."""

        if binary.event == 'END' or binary.event == 'FAILED':
             binary.xrb_luminosity.append(x_ray_luminosity(binary)[0])
             binary.xrb_beaming.append(x_ray_luminosity(binary)[1])
#            diff = int(len(binary.event_history) - len(binary.Lx))
#            binary.xrb_luminosities += [0.0]*diff
#            binary.xrb_beaming += [1.0]*diff

        return binary

