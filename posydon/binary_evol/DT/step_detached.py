"""Detached evolution step."""


__authors__ = [
    "Devina Misra <devina.misra@unige.ch>",
    "Zepei Xing <Zepei.Xing@unige.ch>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Nam Tran <tranhn03@gmail.com>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
    "Camille Liotine <cliotine@u.northwestern.edu>",
    "Seth Gossage <seth.gossage@northwestern.edu>"
]

import numpy as np
import pandas as pd
import time
from scipy.integrate import solve_ivp

from posydon.config import PATH_TO_POSYDON_DATA
from posydon.binary_evol.binarystar import BINARYPROPERTIES
from posydon.binary_evol.singlestar import STARPROPERTIES
#from posydon.interpolation.data_scaling import DataScaler
from posydon.utils.common_functions import (bondi_hoyle,
                                            orbital_period_from_separation,
                                            roche_lobe_radius,
                                            check_state_of_star,
                                            set_binary_to_failed,
                                            zero_negative_values)
from posydon.binary_evol.flow_chart import (STAR_STATES_CC, 
                                            STAR_STATES_H_RICH_EVOLVABLE,
                                            STAR_STATES_HE_RICH_EVOLVABLE,
                                            UNDEFINED_STATES)
import posydon.utils.constants as const
from posydon.utils.posydonerror import (NumericalError, POSYDONError, 
                                        FlowError, ClassificationError)
from posydon.utils.posydonwarning import Pwarn

from posydon.binary_evol.DT.track_match import TrackMatcher
from posydon.binary_evol.DT.key_library import (DEFAULT_TRANSLATION,
                                                DEFAULT_TRANSLATED_KEYS)

def event(terminal, direction=0):
    """Return a helper function to set attributes for solve_ivp events."""
    def dec(f):
        f.terminal = True
        f.direction = direction
        return f
    return dec

class detached_step:
    """
    Evolve a detached binary.

    The binary will be evolved until Roche-lobe overflow, core-collapse or
    maximum simulation time, using the standard equations that govern the
    orbital evolution.

    Parameters
    ----------
    path : str
        Path to the directory that contains POSYDON data HDF5 files. Defaults 
        to the PATH_TO_POSYDON_DATA environment variable. Used for track 
        matching.

    metallicity : float
        The metallicity of the grid. This should be one of the eight 
        supported metallicities:

            [2e+00, 1e+00, 4.5e-01, 2e-01, 1e-01, 1e-02, 1e-03, 1e-04]

        and this will be converted to a corresponding string (e.g.,
        1e+00 --> "1e+00_Zsun"). Used for track matching. 

    matching_method : str
        Method to find the best match between a star from a previous step and a
        point in a single star evolution track. Options: 
        
            "root": Tries to find a root of two matching quantities. It is 
                    possible to not find one, causing the evolution to fail.

            "minimize": Minimizes the sum of squares of differences of 
                        various quantities between the previous evolution step and 
                        a stellar evolution track. 
                        
        Used for track matching.

    grid_name_Hrich : str
        Name of the single star H-rich grid h5 file, 
        including its parent directory. This is set to 
        (for example):

            grid_name_Hrich = 'single_HMS/1e+00_Zsun.h5'  

        by default if not specified. Used for track matching.

    grid_name_strippedHe : str
        Name of the single star He-rich grid h5 file. This is 
        set to (for example):

            grid_name_strippedHe = 'single_HeMS/1e+00_Zsun.h5'
        
        by default if not specified. Used for track matching.

    list_for_matching_HMS : list
        A list of mixed type that specifies properties of the matching 
        process for HMS stars. Used for track matching.

    list_for_matching_postMS : list
        A list of mixed type that specifies properties of the matching 
        process for postMS stars. Used for track matching.

    list_for_matching_HeStar : list
        A list of mixed type that specifies properties of the matching 
        process for He stars. Used for track matching.

    record_matching : bool
        Whether properties of the matched star(s) should be recorded in the 
        binary evolution history. Used for track matching.

    Attributes
    ----------
    KEYS : list[str]
        Contains keywords corresponding to MESA data column names 
        which are used to extract quantities from the single star 
        evolution grids.
        
    dt : float
        The timestep size, in years, to be appended to the history of the
        binary. None means only the final step. Note: do not select very 
        small timesteps because it may mess with the solving of the ODE.

    n_o_steps_history : int
        Alternatively, we can define the number of timesteps to be appended to
        the history of the binary. None means only the final step. If both `dt`
        and `n_o_steps_history` are different than None, `dt` has priority.

    do_wind_loss : bool
        If True, take into account change of separation due to mass loss 
        from the star.

    do_tides : bool
        If True, take into account change of separation, eccentricity and 
        star spin due to tidal forces.

    do_gravitational_radiation : bool
        If True, take into account change of separation and eccentricity 
        due to gravitational wave radiation.

    do_magnetic_braking : bool
        If True, take into account change of star spin due to magnetic 
        braking.

    magnetic_braking_mode : str
        A string corresponding to the desired magnetic braking prescription.
            -- RVJ83: Rappaport, Verbunt, & Joss 1983 (Default)
            -- M15: Matt et al. 2015
            -- G18: Garraffo et al. 2018
            -- CARB: Van & Ivanova 2019

    do_stellar_evolution_and_spin_from_winds : bool
        If True, take into account change of star spin due to change of its
        moment of inertia during its evolution and due to spin angular 
        momentum loss due to winds.

    RLO_orbit_at_orbit_with_same_am : bool
        Binaries are circularized instaneously when RLO occurs and this 
        option dictates how that is handled. If False (default), place 
        the binary in an orbit with separation equal to the binary's 
        separation at periastron. If True, circularize the orbit assuming 
        that angular momentum is conserved w.r.t. the previously (possibly) 
        eccentric orbit. In the latter case, the star may no longer 
        fill its Roche lobe after circularization, and may be further 
        evolved until RLO commences once again, but without changing the 
        orbit.

    translate : dict
        Dictionary containing data column name (key) translations between 
        POSYDON h5 file PSyGrid data names (items) and MESA data names (keys).

    track_matcher : TrackMatcher object
        The TrackMatcher object performs functions related to matching 
        binary stellar evolution components to single star evolution models.

    verbose : bool
        True if we want to print stuff.

    """

    def __init__(
            self,
            dt=None,
            n_o_steps_history=None,
            do_wind_loss=True,
            do_tides=True,
            do_gravitational_radiation=True,
            do_magnetic_braking=True,
            magnetic_braking_mode="RVJ83",
            do_stellar_evolution_and_spin_from_winds=True,
            RLO_orbit_at_orbit_with_same_am=False,
            record_matching=False,
            verbose=False,
            grid_name_Hrich=None,
            grid_name_strippedHe=None,
            metallicity=None,
            path=PATH_TO_POSYDON_DATA,
            matching_method="minimize",
            list_for_matching_HMS=None,
            list_for_matching_postMS=None,
            list_for_matching_HeStar=None
    ):
        """Initialize the step. See class documentation for details."""
        self.dt = dt
        self.n_o_steps_history = n_o_steps_history
        self.do_wind_loss = do_wind_loss
        self.do_tides = do_tides
        self.do_gravitational_radiation = do_gravitational_radiation
        self.do_magnetic_braking = do_magnetic_braking
        self.magnetic_braking_mode = magnetic_braking_mode
        self.do_stellar_evolution_and_spin_from_winds = (
            do_stellar_evolution_and_spin_from_winds
        )
        self.RLO_orbit_at_orbit_with_same_am = RLO_orbit_at_orbit_with_same_am
        self.verbose = verbose

        if self.verbose:
            print(
                dt,
                n_o_steps_history,
                matching_method,
                do_wind_loss,
                do_tides,
                do_gravitational_radiation,
                do_magnetic_braking,
                magnetic_braking_mode,
                do_stellar_evolution_and_spin_from_winds)

        self.translate = DEFAULT_TRANSLATION

        # these are the KEYS read from POSYDON h5 grid files (after translating
        # them to the appropriate columns)
        self.KEYS = DEFAULT_TRANSLATED_KEYS

        # creating a track matching object
        self.track_matcher = TrackMatcher(grid_name_Hrich = grid_name_Hrich,
                                          grid_name_strippedHe = grid_name_strippedHe,
                                          path=path, metallicity = metallicity,
                                          matching_method = matching_method,
                                          list_for_matching_HMS = list_for_matching_HMS,
                                          list_for_matching_HeStar = list_for_matching_HeStar,
                                          list_for_matching_postMS = list_for_matching_postMS,
                                          record_matching = record_matching,
                                          verbose = self.verbose,)

        return

    def __repr__(self):
        """Return the name of evolution step."""
        return "Detached Step."

    def __call__(self, binary):
        """
            Evolve the binary through detached evolution until RLO or 
        compact object formation.

        Parameters
        ----------
        binary : BinaryStar object
            A BinaryStar object containing a binary system's properties. 
            This binary will be evolved through detached evolution here.

        Raises
        ------
        ValueError
            If the max time is exceeded by the current time of 
            evolution.

        NumericalError
            If numerical integration fails for the binary during 
            the calculation of its detached evolution. We mark 
            the binary as failed in this case.

        FlowError
            If the evolution of H-rich/He-rich stars in RLO onto 
            H-rich/He-rich stars after HMS-HMS is detected. This 
            evolution pathway is not yet supported by our grids. 
            The binary evolution is marked as failed at this 
            point.
                
        ClassificationError
            If stable RLO between two HMS-HMS stars is determined 
            as a result of detached evolution. We mark these 
            binaries as failed.

        POSYDONError
            If both stars are calculated to be ready for collapse 
            as a result of detached evolution, but the two stars 
            differ in mass.
             
        """

        # Get simulation properties and step names
        binary_sim_prop = getattr(binary, "properties")
        all_step_names = getattr(binary_sim_prop, "all_step_names")

        # get the next step's name to display for match recording in data frame
        # (in the event that the total_state is not in the flow, this will be None,
        #  and the binary will be set to fail in BinaryStar().run_step()).
        next_step_name = binary.get_next_step_name()
        
        # match stars to single star models for detached evolution
        primary, secondary, only_CO = self.track_matcher.do_matching(binary, next_step_name)
        
        if only_CO:
            if self.verbose:
                print("Binary system only contains compact objects."
                      "Exiting step_detached, nothing left to do here.")
            return

        secondary.t_max = secondary.interp1d["t_max"]
        primary.t_max = primary.interp1d["t_max"]
        secondary.t_offset = binary.time - secondary.interp1d["t0"]
        primary.t_offset = binary.time - primary.interp1d["t0"]
        max_time = secondary.interp1d["max_time"]

        if (self.ev_rlo1(binary.time, [binary.separation, binary.eccentricity], primary, secondary) >= 0
            or self.ev_rlo2(binary.time, [binary.separation, binary.eccentricity], primary, secondary) >= 0):
            binary.state = "initial_RLOF"
            return
        else:
            if not (max_time - binary.time > 0.0):
                raise ValueError("max_time is lower than the current time. "
                                "Evolution of the detached binary will go to "
                                "lower times.")

            with np.errstate(all="ignore"):

                t_before_ODEsolution = time.time()
                try:
                    res = solve_ivp(self.diffeq, 
                                    events=[self.ev_rlo1, self.ev_rlo2, 
                                            self.ev_max_time1, self.ev_max_time2],
                                    method="Radau", 
                                    t_span=(binary.time, max_time),
                                    y0=[binary.separation, binary.eccentricity,
                                        secondary.omega0, primary.omega0],
                                    args = (primary, secondary),
                                    dense_output=True)
                except Exception:
                    res = solve_ivp(self.diffeq,
                                    events=[self.ev_rlo1, self.ev_rlo2, 
                                            self.ev_max_time1, self.ev_max_time2],
                                    method="RK45",
                                    t_span=(binary.time, max_time),
                                    y0=[binary.separation, binary.eccentricity,
                                        secondary.omega0, primary.omega0],
                                    args=(primary, secondary),
                                    dense_output=True)

            t_after_ODEsolution = time.time()

            if self.verbose:
                ivp_tspan = t_after_ODEsolution - t_before_ODEsolution
                print(f"\nODE solver duration: {ivp_tspan:.6g} sec")
                print("solution of ODE", res)

            if res.status == -1:
                failed_state = binary.state
                set_binary_to_failed(binary)
                raise NumericalError(f"Integration failed for {failed_state} binary.")
                         
            # update binary/star properties after detached evolution
            t = self.get_time_after_evo(res, binary)
            self.update_after_evo(res, t, binary, primary, secondary)

            # check primary/secondary star states
            secondary.state = check_state_of_star(secondary, star_CO=False)
            for timestep in range(-len(t[:-1]), 0):
                secondary.state_history[timestep] = check_state_of_star(secondary, i=timestep, star_CO=False)

            if primary.state == "massless_remnant":
                pass
            elif primary.co:
                mdot_acc = np.atleast_1d(bondi_hoyle(
                    binary, primary, secondary, slice(-len(t), None),
                    wind_disk_criteria=True, scheme='Kudritzki+2000'))
                primary.lg_mdot = np.log10(mdot_acc.item(-1))
                primary.lg_mdot_history[len(primary.lg_mdot_history) - len(t) + 1:] = np.log10(mdot_acc[:-1])
            else:
                primary.state = check_state_of_star(primary, star_CO=False)
                for timestep in range(-len(t[:-1]), 0):
                    primary.state_history[timestep] = check_state_of_star(primary, i=timestep, star_CO=False)

            ## CHECK IF THE BINARY IS IN RLO
            if res.t_events[0] or res.t_events[1]:
                if self.RLO_orbit_at_orbit_with_same_am:
                    # final circular orbit conserves angular momentum
                    # compared to the eccentric orbit
                    binary.separation *= (1 - res.y[1][-1]**2)
                    binary.orbital_period *= (1 - res.y[1][-1]**2) ** 1.5
                else:
                    # final circular orbit is at periastron of the ecc. orbit
                    binary.separation *= (1 - res.y[1][-1])
                    binary.orbital_period *= (1 - res.y[1][-1]) ** 1.5

                abs_diff_porb = np.abs(binary.orbital_period - orbital_period_from_separation(
                                binary.separation, secondary.mass, primary.mass)) / binary.orbital_period

                assert abs_diff_porb < 10 ** (-2),  \
                f"\nabs_diff_porb = {abs_diff_porb:.4f}" + \
                f"\nbinary.orbital_period = {binary.orbital_period:.4f}" +\
                "\norbital_period_from_separation(binary.separation, secondary.mass, primary.mass) =" + \
                f"{orbital_period_from_separation(binary.separation, secondary.mass, primary.mass):.4f}"

                # instantly circularize at RLO
                binary.eccentricity = 0

                if res.t_events[0]:
                    if secondary == binary.star_1:
                        binary.state = "RLO1"
                        binary.event = "oRLO1"
                    else:
                        binary.state = "RLO2"
                        binary.event = "oRLO2"

                elif res.t_events[1]:
                    if secondary == binary.star_1:
                        binary.state = "RLO2"
                        binary.event = "oRLO2"
                    else:
                        binary.state = "RLO1"
                        binary.event = "oRLO1"

                if ('step_HMS_HMS_RLO' not in all_step_names):
                    if ((binary.star_1.state in STAR_STATES_HE_RICH_EVOLVABLE 
                         and binary.star_2.state in STAR_STATES_H_RICH_EVOLVABLE)
                    or (binary.star_1.state in STAR_STATES_H_RICH_EVOLVABLE
                         and binary.star_2.state in STAR_STATES_HE_RICH_EVOLVABLE)):
                        set_binary_to_failed(binary)
                        raise FlowError("Evolution of H-rich/He-rich stars in RLO onto H-rich/He-rich stars after " 
                                    "HMS-HMS not yet supported.") 

                    elif (binary.star_1.state in STAR_STATES_H_RICH_EVOLVABLE
                         and binary.star_2.state in STAR_STATES_H_RICH_EVOLVABLE):
                        set_binary_to_failed(binary)
                        raise ClassificationError("Binary is in the detached step but has stable RLO with two HMS stars - "
                                              "should it have undergone CE (was its HMS-HMS interpolation class unstable MT?)") 


            ## CHECK IF STARS WILL UNDERGO CC
            elif res.t_events[2]:
                # reached t_max of track. End of life (possible collapse) of secondary
                if secondary == binary.star_1:
                    binary.event = "CC1"
                else:
                    binary.event = "CC2"

                self.track_matcher.get_star_final_values(secondary)
                self.track_matcher.get_star_profile(secondary)

                if not primary.co and primary.state in STAR_STATES_CC:
                    # simultaneous core-collapse of the other star as well
                    primary_time = primary.t_max + primary.t_offset - t[-1]
                    secondary_time = secondary.t_max + secondary.t_offset - t[-1]

                    if primary_time == secondary_time:
                        # we manually check if s.t_events[3] should also be happening simultaneously
                        self.track_matcher.get_star_final_values(primary)
                        self.track_matcher.get_star_profile(primary)

                    if primary.mass != secondary.mass:
                        raise POSYDONError(
                            "Both stars are found to be ready for collapse "
                            "(i.e. end of their life) during the detached "
                            "step, but do not have the same mass")

            elif res.t_events[3]:
                # reached t_max of track. End of life (possible collapse) of primary
                if secondary == binary.star_1:
                    binary.event = "CC2"
                else:
                    binary.event = "CC1"

                self.track_matcher.get_star_final_values(primary)
                self.track_matcher.get_star_profile(primary)

            else:  # Reached max_time asked.
                if binary.properties.max_simulation_time - binary.time < 0.0:
                    binary.event = "MaxTime_exceeded"
                else:
                    binary.event = "maxtime"

    @event(True, 1)
    def ev_rlo1(self, t, y, primary, secondary):
        """
            Difference between radius and Roche lobe at a given time. Used 
        to check if there is RLOF mass transfer during the detached binary 
        evolution interpolation. Calculated for the secondary.

        Parameters
        ----------
        t : float
            Time of the evolution, in years.

        y : tuple(float)
            [separation, eccentricity] at that time. Separation should be
            in solar radii.

        primary : SingleStar object
            A single star object, representing the primary (more evolved) star 
            in the binary and containing its properties.
        
        secondary : SingleStar object
            A single star object, representing the secondary (less evolved) star 
            in the binary and containing its properties.

        Returns
        -------
        RL_diff : float
            Difference between stellar radius and 95% of the Roche lobe 
            radius in solar radii.

        """
        pri_mass = primary.interp1d["mass"](t - primary.t_offset)
        sec_mass = secondary.interp1d["mass"](t - secondary.t_offset)

        sep = y[0]
        ecc = y[1]
        
        RL = roche_lobe_radius(sec_mass, pri_mass, (1 - ecc) * sep)
        
        # 95% filling of the RL is enough to assume beginning of RLO,
        # as we do in CO-HMS_RLO grid
        RL_diff = secondary.interp1d["R"](t - secondary.t_offset) - 0.95*RL
        return RL_diff

    @event(True, 1)
    def ev_rlo2(self, t, y, primary, secondary):
        """
            Difference between radius and Roche lobe at a given time. Used 
        to check if there is RLOF mass transfer during the detached binary 
        evolution interpolation. Calculated for the primary.

        Parameters
        ----------
        t : float
            Time of the evolution, in years

        y : tuple(float)
            [separation, eccentricity] at that time. Separation should be
            in solar radii.

        primary : SingleStar object
            A single star object, representing the primary (more evolved) star 
            in the binary and containing its properties.
        
        secondary : SingleStar object
            A single star object, representing the secondary (less evolved) star 
            in the binary and containing its properties.

        Returns
        -------
        RL_diff : float
            Difference between stellar radius and 95% of the Roche lobe 
            radius in solar radii.

        """
        pri_mass = primary.interp1d["mass"](t - primary.t_offset)
        sec_mass = secondary.interp1d["mass"](t - secondary.t_offset)
        
        sep = y[0]
        ecc = y[1]
        
        RL = roche_lobe_radius(pri_mass, sec_mass, (1 - ecc) * sep)
        RL_diff = primary.interp1d["R"](t - primary.t_offset) - 0.95*RL

        return RL_diff

    @event(True, 1)
    def ev_rel_rlo1(self, t, y, primary, secondary):
        """
            Relative difference between radius and Roche lobe. Used to 
        check if there is RLOF mass transfer during the detached binary 
        evolution interpolation. Calculated for the secondary.

        Parameters
        ----------
        t : float
            Time of the evolution, in years.

        y : tuple(float)
            [separation, eccentricity] at that time. Separation should be
            in solar radii.

        primary : SingleStar object
            A single star object, representing the primary (more evolved) star 
            in the binary and containing its properties.
        
        secondary : SingleStar object
            A single star object, representing the secondary (less evolved) star 
            in the binary and containing its properties.

        Returns
        -------
        RL_rel_diff : float
            Relative difference between stellar radius and Roche lobe
            radius.

        """
        pri_mass = primary.interp1d["mass"](t - primary.t_offset)
        sec_mass = secondary.interp1d["mass"](t - secondary.t_offset)
        
        sep = y[0]
        ecc = y[1]
        
        RL = roche_lobe_radius(sec_mass, pri_mass, (1 - ecc) * sep)
        RL_rel_diff = (secondary.interp1d["R"](t - secondary.t_offset) - RL) / RL
        return RL_rel_diff

    @event(True, 1)
    def ev_rel_rlo2(self, t, y, primary, secondary):
        """
            Relative difference between radius and Roche lobe. Used to 
        check if there is RLOF mass transfer during the detached binary 
        evolution interpolation. Calculated for the primary.

        Parameters
        ----------
        t : float
            Time of the evolution, in years.

        y : tuple(float)
            [separation, eccentricity] at that time. Separation should be
            in solar radii.

        primary : SingleStar object
            A single star object, representing the primary (more evolved) star 
            in the binary and containing its properties.
        
        secondary : SingleStar object
            A single star object, representing the secondary (less evolved) star 
            in the binary and containing its properties.

        Returns
        -------
        RL_rel_diff : float
            Relative difference between stellar radius and Roche lobe
            radius.
        """
        pri_mass = primary.interp1d["mass"](t - primary.t_offset)
        sec_mass = secondary.interp1d["mass"](t - secondary.t_offset)

        sep = y[0]
        ecc = y[1]
        
        RL = roche_lobe_radius(pri_mass, sec_mass, (1 - ecc) * sep)
        RL_rel_diff = (primary.interp1d["R"](t - primary.t_offset) - RL) / RL
        return RL_rel_diff

    @event(True, -1)
    def ev_max_time1(self, t, y, primary, secondary):
        return secondary.t_max + secondary.t_offset - t

    @event(True, -1)
    def ev_max_time2(self, t, y, primary, secondary):
        return primary.t_max + primary.t_offset - t
    
    def get_time_after_evo(self, res, binary):
        """
            After detached evolution, this uses the ODESolver result 
        to determine what the current time is.

        Parameters
        ----------
        res : ODESolver object
            This is the ODESolver object produced by SciPy's 
            solve_ivp function that contains calculated values 
            of the stars evolution through the detached step.

        binary: BinaryStar object
            A binary star object, containing the binary system's properties.
        
        Returns
        -------
        t : float or array[float]
            This is the time elapsed as a result of detached 
            evolution in years. This is a float unless the 
            user specifies a timestep (see `n_o_steps_history` 
            or `dt`) to use via the simulation properties ini 
            file, in which case it is an array.
        
        """
        
        if self.dt is not None and self.dt > 0:
            t = np.arange(binary.time, res.t[-1] + self.dt/2.0, self.dt)[1:]
            if t[-1] < res.t[-1]:
                t = np.hstack([t, res.t[-1]])
        elif (self.n_o_steps_history is not None
                and self.n_o_steps_history > 0):
            t_step = (res.t[-1] - binary.time) / self.n_o_steps_history
            t = np.arange(binary.time, res.t[-1] + t_step / 2.0, t_step)[1:]
            if t[-1] < res.t[-1]:
                t = np.hstack([t, res.t[-1]])
        else:  # self.dt is None and self.n_o_steps_history is None
            t = np.array([res.t[-1]])

        return t

    def update_after_evo(self, res, t, binary, primary, secondary):

        """
            Update star and binary properties and interpolators with 
        ODESolver result from detached evolution. This update gives 
        the binary/stars their appropriate values, according to the 
        interpolation after detached evolution.

        Parameters
        ----------
        res : ODESolver object
            This is the ODESolver object produced by SciPy's 
            solve_ivp function that contains calculated values 
            of the stars evolution through the detached step.

        t : float or array[float]
            This is the time elapsed as a result of detached 
            evolution in years. This is a float unless the 
            user specifies a timestep to use via the simulation 
            properties ini file, in which case it is an array.

        binary : BinaryStar object
            A binary star object, containing the binary system's properties.

        primary : SingleStar object
            A single star object, representing the primary (more evolved) star 
            in the binary and containing its properties.
        
        secondary : SingleStar object
            A single star object, representing the secondary (less evolved) star 
            in the binary and containing its properties.

        Warns
        -----
        InappropriateValueWarning
            If trying to compute log angular momentum for object with no spin.

        """

        sep_interp, ecc_interp, omega_interp_sec, omega_interp_pri = res.sol(t)
        mass_interp_sec = secondary.interp1d[self.translate["mass"]]
        mass_interp_pri = primary.interp1d[self.translate["mass"]]

        secondary.interp1d["sep"] = sep_interp
        secondary.interp1d["ecc"] = ecc_interp    
        secondary.interp1d["omega"] = omega_interp_sec
        primary.interp1d["omega"] = omega_interp_pri

        secondary.interp1d["porb"] = orbital_period_from_separation(
            sep_interp, mass_interp_sec(t - secondary.t_offset),
            mass_interp_pri(t - primary.t_offset))
        primary.interp1d["porb"] = orbital_period_from_separation(
            sep_interp, mass_interp_pri(t - primary.t_offset),
            mass_interp_sec(t - secondary.t_offset))

        secondary.interp1d["time"] = t

        for obj, prop in zip([secondary, primary, binary], 
                                [STARPROPERTIES, STARPROPERTIES, BINARYPROPERTIES]):
            
            interp1d = primary.interp1d if obj == primary else secondary.interp1d
            t_offset = primary.t_offset if obj == primary else secondary.t_offset

            for key in prop:
                if key in ["event",
                            "mass_transfer_case",
                            "nearest_neighbour_distance",
                            "state", "metallicity", "V_sys"]:
                    current = getattr(obj, key)
                    # For star objects, the state is calculated further below
                    history = [current] * len(t[:-1])

                # replace the actual surf_avg_w with the effective omega,
                # which takes into account the whole star
                # key  = 'effective_omega' # in rad/sec
                # current = s.y[2][-1] / 3.1558149984e7
                # history_of_attribute = s.y[2][:-1] / 3.1558149984e7
                elif (key in ["surf_avg_omega_div_omega_crit"] and obj != binary):#primary):
                    if obj.co: #primary.co:
                        current = None
                        history = [current] * len(t[:-1])

                    else:
                        # TODO: change `item()` to 0
                        omega_crit_current = np.sqrt(const.standard_cgrav
                            * interp1d[self.translate["mass"]](t[-1] - t_offset).item() * const.msol
                            / (interp1d[self.translate["R"]](t[-1] - t_offset).item() * const.rsol)**3)

                        omega_crit_hist = np.sqrt(const.standard_cgrav
                            * interp1d[self.translate["mass"]](t[:-1] - t_offset) * const.msol
                            / (interp1d[self.translate["R"]](t[:-1] - t_offset) * const.rsol)**3)

                        current = (interp1d["omega"][-1] / const.secyer / omega_crit_current)
                        history = (interp1d["omega"][:-1] / const.secyer / omega_crit_hist)

                        # ensure positive rotation values
                        current = zero_negative_values([current], key)[0]
                        history = zero_negative_values(history, key)

                elif (key in ["surf_avg_omega"] and obj != binary):
                    if obj.co:
                        current = None
                        history = [current] * len(t[:-1])
                    else:
                        current = interp1d["omega"][-1] / const.secyer
                        history = interp1d["omega"][:-1] / const.secyer

                        current = zero_negative_values([current], key)[0]
                        history = zero_negative_values(history, key)
                        
                elif ("rl_relative_overflow_" in key and obj == binary):
                    s = binary.star_1 if "_1" in key[-2:] else binary.star_2
                    s_alt = binary.star_2 if "_1" in key[-2:] else binary.star_1
                    if s.state in ("BH", "NS", "WD","massless_remnant"):
                        current = None
                        history = [current] * len(t[:-1])

                    elif secondary == s:
                        current = self.ev_rel_rlo1(t[-1], [interp1d["sep"][-1], interp1d["ecc"][-1]], primary, secondary)
                        history = self.ev_rel_rlo1(t[:-1], [interp1d["sep"][:-1], interp1d["ecc"][:-1]], primary, secondary)

                    elif secondary == s_alt:
                        current = self.ev_rel_rlo2(t[-1], [interp1d["sep"][-1], interp1d["ecc"][-1]], primary, secondary)
                        history = self.ev_rel_rlo2(t[:-1], [interp1d["sep"][:-1], interp1d["ecc"][:-1]], primary, secondary)
                        
                elif key in ["separation", "orbital_period", "eccentricity", "time"]:
                    current = interp1d[self.translate[key]][-1].item()
                    history = interp1d[self.translate[key]][:-1]

                    current = zero_negative_values([current], key)[0]
                    history = zero_negative_values(history, key)
                    
                elif (key in ["total_moment_of_inertia"] and obj != binary):
                    if obj.co:
                        current = getattr(obj, key)
                        history = [current] * len(t[:-1])
                    else:
                        current = interp1d[self.translate[key]](
                            t[-1] - t_offset).item() * (const.msol * const.rsol**2)

                        history = interp1d[self.translate[key]](
                            t[:-1] - t_offset) * (const.msol * const.rsol**2)
                        
                        current = zero_negative_values([current], key)[0]
                        history = zero_negative_values(history, key)
                    
                elif (key in ["log_total_angular_momentum"] and obj != binary):
                    if obj.co:
                        current = getattr(obj, key)
                        history = [current] * len(t[:-1])
                    else:
                        tot_j = (interp1d["omega"][-1] / const.secyer) \
                                  * (interp1d[self.translate["total_moment_of_inertia"]]( \
                                    t[-1] - t_offset).item() * (const.msol * const.rsol**2))
                        current = np.log10(tot_j) if tot_j > 0.0 else -99
                    
                        tot_j_hist = (interp1d["omega"][:-1] / const.secyer) \
                                       * (interp1d[self.translate["total_moment_of_inertia"]]( \
                                       t[:-1] - t_offset) * (const.msol * const.rsol**2))
                        history = np.where(tot_j_hist > 0, np.log10(tot_j_hist), -99)

                        current = zero_negative_values([current], key)[0]
                        history = zero_negative_values(history, key)
                    
                elif (key in ["spin"] and obj != binary):
                    if obj.co:
                        current = getattr(obj, key)
                        history = [current] * len(t[:-1])
                    else:
                        current = (const.clight
                            * (interp1d["omega"][-1] / const.secyer)
                            * interp1d[self.translate["total_moment_of_inertia"]](
                                    t[-1] - t_offset).item() * (const.msol * const.rsol**2)
                            / (const.standard_cgrav * (interp1d[self.translate["mass"]](
                                    t[-1] - t_offset).item() * const.msol)**2))
                        
                        history = (const.clight 
                            * (interp1d["omega"][:-1] / const.secyer)
                            * interp1d[self.translate["total_moment_of_inertia"]](
                                    t[:-1] - t_offset) * (const.msol * const.rsol**2)
                            / (const.standard_cgrav * (interp1d[self.translate["mass"]](
                                    t[:-1] - t_offset) * const.msol)**2))

                elif (key in ["lg_mdot", "lg_wind_mdot"] and obj != binary):
                    if obj.co:
                        current = None
                        history = [current] * len(t[:-1])
                    else:
                        if interp1d[self.translate[key]](t[-1] - t_offset) == 0:
                            current = -98.99
                        else:
                            current = np.log10(np.abs(interp1d[self.translate[key]](
                                    t[-1] - t_offset))).item()
                            
                        history = np.ones_like(t[:-1])
                        for i in range(len(t)-1):
                            if (interp1d[self.translate[key]](t[i] - t_offset) == 0):
                                history[i] = -98.99
                            else:
                                history[i] = np.log10(np.abs(interp1d[self.translate[key]](
                                            t[i] - t_offset)))
                    
                elif (self.translate[key] in interp1d and obj != binary):
                    if obj.co:
                        current = getattr(obj, key)
                        history = [current] * len(t[:-1])
                    else:
                        current = interp1d[self.translate[key]](t[-1] - t_offset).item()
                        history = interp1d[self.translate[key]](t[:-1] - t_offset)
                        
                elif key in ["profile"]:
                    current = None
                    history = [current] * len(t[:-1])

                else:
                    current = np.nan
                    history = np.ones_like(t[:-1]) * current

                setattr(obj, key, current)
                getattr(obj, key + "_history").extend(history)

    def diffeq(self, t, y, primary, secondary):
        """
            Diff. equation describing the orbital evolution of a detached binary.

        The equation handles wind mass-loss [1]_, tidal [2]_, gravational [3]_
        effects and magnetic braking [4]_, [5]_, [6]_, [7]_, [8]_. It also handles
        the change of the secondary's stellar spin due to its change of moment of
        intertia and due to mass-loss from its spinning surface. It is assumed that
        the mass loss is fully non-conservative. Magnetic braking is fully applied
        to secondary stars with mass less than 1.3 Msun and fully off for stars
        with mass larger then 1.5 Msun. The effect of magnetic braking falls
        linearly for stars with mass between 1.3 Msun and 1.5 Msun.

        TODO: explain new features (e.g., double COs)

        Parameters
        ----------
        t : float
            The age of the system in years

        y : list[float]
            Contains the separation, eccentricity and angular velocity, in Rsolar,
            dimensionless and rad/year units, respectively.

        primary : SingleStar object
            A single star object, representing the primary (more evolved) star 
            in the binary and containing its properties.
        
        secondary : SingleStar object
            A single star object, representing the secondary (less evolved) star 
            in the binary and containing its properties.

        Warns
        -----
        UnsupportedModelWarning
            If an unsupported model or model is unspecified is determined from 
            `magnetic_braking_mode`. In this case, magnetic braking will not 
            be calculated during the detached step.

        Returns
        -------
        result : list[float]
            Contains the change of the separation, eccentricity and angular 
            velocity, in Rsolar, dimensionless and rad/year units, respectively.

        References
        ----------
        .. [1] Tauris, T. M., & van den Heuvel, E. 2006,
            Compact stellar X-ray sources, 1, 623
        .. [2] Hut, P. 1981, A&A, 99, 126
        .. [3] Junker, W., & Schafer, G. 1992, MNRAS, 254, 146
        .. [4] Rappaport, S., Joss, P. C., & Verbunt, F. 1983, ApJ, 275, 713
        .. [5] Matt et al. 2015, ApJ, 799, L23
        .. [6] Garraffo et al. 2018, ApJ, 862, 90
        .. [7] Van & Ivanova 2019, ApJ, 886, L31
        .. [8] Gossage et al. 2021, ApJ, 912, 65

        """
        # secondary star properties used in evolution
        R_sec = secondary.interp1d["R"](t - secondary.t_offset) # Rsol
        L_sec = secondary.interp1d["L"](t - secondary.t_offset) # Lsol
        M_sec = secondary.interp1d["mass"](t - secondary.t_offset) # Msol
        Mdot_sec = secondary.interp1d["mdot"](t - secondary.t_offset) # Msol/yr
        I_sec = secondary.interp1d["inertia"](t - secondary.t_offset) # Msol*Rsol^2
        conv_mx1_top_r_sec = secondary.interp1d["conv_mx1_top_r"](t - secondary.t_offset) # Rsol
        conv_mx1_bot_r_sec = secondary.interp1d["conv_mx1_bot_r"](t - secondary.t_offset) # Rsol
        surface_h1_sec = secondary.interp1d["surface_h1"](t - secondary.t_offset) # dex
        center_h1_sec = secondary.interp1d["center_h1"](t - secondary.t_offset) # dex
        M_env_sec = secondary.interp1d["mass_conv_reg_fortides"](t - secondary.t_offset) # Msol
        DR_env_sec = secondary.interp1d["thickness_conv_reg_fortides"](t - secondary.t_offset) # Rsol
        Renv_middle_sec = secondary.interp1d["radius_conv_reg_fortides"](t - secondary.t_offset) # Rsol
        Idot_sec = secondary.interp1d["Idot"](t - secondary.t_offset) # Msol*Rsol^2/yr
        tau_conv_sec = secondary.interp1d["conv_env_turnover_time_l_b"](t - secondary.t_offset) # yr

        # primary star properties used in evolution
        R_pri = primary.interp1d["R"](t - primary.t_offset)
        L_pri = primary.interp1d["L"](t - primary.t_offset)
        M_pri = primary.interp1d["mass"](t - primary.t_offset)
        Mdot_pri = primary.interp1d["mdot"](t - primary.t_offset)
        I_pri = primary.interp1d["inertia"](t - primary.t_offset)
        conv_mx1_top_r_pri = primary.interp1d["conv_mx1_top_r"](t - primary.t_offset)
        conv_mx1_bot_r_pri = primary.interp1d["conv_mx1_bot_r"](t - primary.t_offset)
        surface_h1_pri = primary.interp1d["surface_h1"](t - primary.t_offset)
        center_h1_pri = primary.interp1d["center_h1"](t - primary.t_offset)
        M_env_pri = primary.interp1d["mass_conv_reg_fortides"](t - primary.t_offset)
        DR_env_pri = primary.interp1d["thickness_conv_reg_fortides"](t - primary.t_offset)
        Renv_middle_pri = primary.interp1d["radius_conv_reg_fortides"](t - primary.t_offset)
        Idot_pri = primary.interp1d["Idot"](t - primary.t_offset)
        tau_conv_pri = primary.interp1d["conv_env_turnover_time_l_b"](t - primary.t_offset)
    
        # detached evolution options
        do_wind_loss = self.do_wind_loss
        do_tides = self.do_tides
        do_gravitational_radiation = self.do_gravitational_radiation
        do_magnetic_braking = self.do_magnetic_braking
        magnetic_braking_mode = self.magnetic_braking_mode
        do_stellar_evolution_and_spin_from_winds = self.do_stellar_evolution_and_spin_from_winds
        verbose = False
                        

        y[0] = np.max([y[0], 0])  # We limit separation to non-negative values
        a = y[0]
        y[1] = np.max([y[1], 0])  # We limit eccentricity to non-negative values
        e = y[1]
        if e > 0 and e < 10.0 ** (-3):
            # we force a negligible eccentricity to become 0
            # for computational stability
            e = 0.0
            if verbose and verbose != 1:
                print("Negligible eccentricity became 0 for "
                      "computational stability")
        y[2] = np.max([y[2], 0])  # We limit omega spin to non-negative values
        Omega_sec = y[2]  # in rad/yr
        y[3] = np.max([y[3], 0])
        Omega_pri = y[3]

        da = 0.0
        de = 0.0
        dOmega_sec = 0.0
        dOmega_pri = 0.0
        #  Mass Loss
        if do_wind_loss:
            q1 = M_sec / M_pri
            k11 = (1 / (1 + q1)) * (Mdot_sec / M_sec)
            k21 = Mdot_sec / M_sec
            k31 = Mdot_sec / (M_pri + M_sec)
            # This is simplified to da_mt = -a * Mdot/(M+Macc), for only (negative)
            # wind Mdot from star M.
            da_mt_sec = a * (2 * k11 - 2 * k21 + k31)

            q2 = M_pri / M_sec
            k12 = (1 / (1 + q2)) * (Mdot_pri / M_pri)
            k22 = Mdot_pri / M_pri
            k32 = Mdot_pri / (M_pri + M_sec)
            da_mt_pri = a * (
                    2 * k12 - 2 * k22 + k32
            )
            if verbose and verbose != 1:
                print("da_mt = ", da_mt_sec, da_mt_pri)

            da = da + da_mt_sec + da_mt_pri

        #  Tidal forces
        if do_tides:
            q1 = M_pri / M_sec
            q2 = M_sec / M_pri
            # P_orb in years. From 3rd Kepler's law, transforming separation from
            # Rsolar to AU, to avoid using constants
            # TODO: we aleady have this function!
            P_orb = np.sqrt((a / const.aursun) ** 3 / (M_pri + M_sec))
            n = 2.0 * const.pi / P_orb  # mean orbital ang. vel. in rad/year
            f1 = (
                    1
                    + (31 / 2) * e ** 2
                    + (255 / 8) * e ** 4
                    + (185 / 16) * e ** 6
                    + (25 / 64) * e ** 8
            )
            f2 = 1 + (15 / 2) * e ** 2 + (45 / 8) * e ** 4 + (5 / 16) * e ** 6
            f3 = 1 + (15 / 4) * e ** 2 + (15 / 8) * e ** 4 + (5 / 64) * e ** 6
            f4 = 1 + (3 / 2) * e ** 2 + (1 / 8) * e ** 4
            f5 = 1 + 3 * e ** 2 + (3 / 8) * e ** 4

            # equilibrium timecale
            if ((pd.notna(M_env_sec) and M_env_sec != 0.0)
                    and (pd.notna(DR_env_sec) and DR_env_sec != 0.0)
                    and (pd.notna(Renv_middle_sec) and Renv_middle_sec != 0.0)):
                # eq. (31) of Hurley et al. 2002, generalized for convective layers
                # not on surface too
                tau_conv_sec = 0.431 * ((M_env_sec * DR_env_sec * Renv_middle_sec
                                        / (3 * L_sec)) ** (1.0 / 3.0))
            else:
                if verbose and verbose != 1:
                    print("something wrong with M_env/DR_env/Renv_middle",
                        M_env_sec, DR_env_sec, Renv_middle_sec)
                tau_conv_sec = 1.0e99
            if ((pd.notna(M_env_pri) and M_env_pri != 0.0)
                    and (pd.notna(DR_env_pri) and DR_env_pri != 0.0)
                    and (pd.notna(Renv_middle_pri) and Renv_middle_pri != 0.0)):
                # eq. (31) of Hurley et al. 2002, generalized for convective layers
                # not on surface too
                tau_conv_pri = 0.431 * ((M_env_pri * DR_env_pri * Renv_middle_pri
                                        / (3 * L_pri)) ** (1.0/3.0))
            else:
                if verbose and verbose != 1:
                    print("something wrong with M_env/DR_env/Renv_middle",
                        M_env_pri, DR_env_pri, Renv_middle_pri)
                tau_conv_pri = 1.0e99
            P_spin_sec = 2 * np.pi / Omega_sec
            P_spin_pri = 2 * np.pi / Omega_pri
            P_tid_sec = np.abs(1 / (1 / P_orb - 1 / P_spin_sec))
            P_tid_pri = np.abs(1 / (1 / P_orb - 1 / P_spin_pri))
            f_conv_sec = np.min([1, (P_tid_sec / (2 * tau_conv_sec)) ** 2])
            f_conv_pri = np.min([1, (P_tid_pri / (2 * tau_conv_pri)) ** 2])
            F_tid = 1.  # not 50 as before
            kT_conv_sec = (
                    (2. / 21) * (f_conv_sec / tau_conv_sec) * (M_env_sec / M_sec)
            )  # eq. (30) of Hurley et al. 2002
            kT_conv_pri = (
                    (2. / 21) * (f_conv_pri / tau_conv_pri) * (M_env_pri / M_pri)
            )
            if kT_conv_sec is None or not np.isfinite(kT_conv_sec):
                kT_conv_sec = 0.0
            if kT_conv_pri is None or not np.isfinite(kT_conv_pri):
                kT_conv_pri = 0.0

                if verbose:
                    print("kT_conv_sec is", kT_conv_sec, ", set to 0.")
                    print("kT_conv_pri is", kT_conv_pri, ", set to 0.")
            # this is the 1/timescale of all d/dt calculted below in yr^-1

            if verbose and verbose != 1:
                print(
                    "Equilibrium tides in deep convective envelope",
                    M_env_sec,
                    DR_env_sec,
                    Renv_middle_sec,
                    R_sec,
                    M_sec,
                    M_env_pri,
                    DR_env_pri,
                    Renv_middle_pri,
                    R_pri,
                    M_pri
                )
                print("convective tiimescales and efficiencies",
                    tau_conv_sec, P_orb, P_spin_sec, P_tid_sec,
                    f_conv_sec,
                    F_tid,
                    tau_conv_pri, P_orb, P_spin_pri, P_tid_pri,
                    f_conv_pri,
                    F_tid,
                    )

            # dynamical timecale
            F_tid = 1
            # E2 = 1.592e-9*M**(2.84) # eq. (43) of Hurley et al. 2002. Deprecated
            R_conv_sec = conv_mx1_top_r_sec - conv_mx1_bot_r_sec
            R_conv_pri = conv_mx1_top_r_pri - conv_mx1_bot_r_pri  # convective core
            # R_conv = conv_mx1_top_r  # convective core
            if (R_conv_sec > R_sec or R_conv_sec <= 0.0
                    or conv_mx1_bot_r_sec / R_sec > 0.1):
                # R_conv = 0.5*R
                # if verbose:
                #     print(
                #         "R_conv of the convective core is not behaving well or "
                #         "we are not calculating the convective core, we make it "
                #         "equal to half of Rstar",
                #         R_conv,
                #         R,
                #         conv_mx1_top_r,
                #         conv_mx1_bot_r,
                #     )
                # we switch to Zahn+1975 calculation of E2
                E21 = 1.592e-9 * M_sec ** (2.84)
            else:
                if R_sec <= 0:
                    E21 = 0
                elif surface_h1_sec > 0.4:
                    E21 = 10.0 ** (-0.42) * (R_conv_sec / R_sec) ** (
                        7.5
                    )  # eq. (9) of Qin et al. 2018, 616, A28
                elif surface_h1_sec <= 0.4:  # "HeStar":
                    E21 = 10.0 ** (-0.93) * (R_conv_sec / R_sec) ** (
                        6.7
                    )  # eq. (9) of Qin et al. 2018, 616, A28
                else:  # in principle we should not go here
                    E21 = 1.592e-9 * M_sec ** (
                        2.84
                    )  # eq. (43) of Hurley et al. 2002 from Zahn+1975 Depreciated
            # kT = 1.9782e4 * np.sqrt(M * R**2 / a**5) * (1 + q)**(5. / 6) * E2
            # eq. (42) of Hurley et al. 2002. Depreciated

            if (R_conv_pri > R_pri or R_conv_pri <= 0.0
                    or conv_mx1_bot_r_pri / R_pri > 0.1):
                E22 = 1.592e-9 * M_pri ** (2.84)
                if verbose and verbose != 1:
                    print(
                        "R_conv of the convective core is not behaving well or we "
                        "are not calculating the convective core, we switch to "
                        "Zahn+1975 calculation of E2",
                        R_conv_sec,
                        R_sec,
                        conv_mx1_top_r_sec,
                        conv_mx1_bot_r_sec,
                        E21,
                        R_conv_pri,
                        R_pri,
                        conv_mx1_top_r_pri,
                        conv_mx1_bot_r_pri,
                        E22
                    )
            else:
                if R_pri <= 0:
                    E22 = 0
                elif surface_h1_pri > 0.4:
                    E22 = 10.0 ** (-0.42) * (R_conv_pri / R_pri) ** (
                        7.5
                    )  # eq. (9) of Qin et al. 2018, 616, A28
                elif surface_h1_pri <= 0.4:  # "HeStar":
                    E22 = 10.0 ** (-0.93) * (R_conv_pri / R_pri) ** (
                        6.7
                    )  # eq. (9) of Qin et al. 2018, 616, A28
                else:  # in principle we should not go here
                    E22 = 1.592e-9 * M_pri ** (
                        2.84
                    )
            kT_rad_sec = (
                np.sqrt(const.standard_cgrav * (M_sec * const.msol)
                        * (R_sec * const.rsol)**2 / (a * const.rsol)**5)
                * (1 + q1) ** (5.0 / 6)
                * E21
                * const.secyer)
            kT_rad_pri = (
                np.sqrt(const.standard_cgrav * (M_pri * const.msol)
                        * (R_pri * const.rsol)**2 / (a * const.rsol)**5)
                * (1 + q2) ** (5.0 / 6)
                * E22
                * const.secyer)
            # this is the 1/timescale of all d/dt calculted below in yr^-1
            if verbose and verbose != 1:
                print(
                    "Dynamical tides in radiative envelope",
                    conv_mx1_top_r_sec,
                    conv_mx1_bot_r_sec,
                    R_conv_sec,
                    E21,
                    conv_mx1_top_r_pri,
                    conv_mx1_bot_r_pri,
                    R_conv_pri,
                    E22,
                    F_tid
                )
            kT_sec = max(kT_conv_sec, kT_rad_sec)
            kT_pri = max(kT_conv_pri, kT_rad_pri)
            if verbose and verbose != 1:
                print("kT_conv/rad of tides is ", kT_conv_sec, kT_rad_sec,
                    kT_conv_pri, kT_rad_pri, "in 1/yr, and we picked the ",
                    kT_sec, kT_pri)

            da_tides_sec = (
                    -6
                    * F_tid
                    * kT_sec
                    * q1
                    * (1 + q1)
                    * (R_sec / a) ** 8
                    * (a / (1 - e ** 2) ** (15 / 2))
                    * (f1 - (1 - e ** 2) ** (3 / 2) * f2 * Omega_sec / n)
            )  # eq. (9) of Hut 1981, 99, 126

            da_tides_pri = (
                    -6
                    * F_tid
                    * kT_pri
                    * q2
                    * (1 + q2)
                    * (R_pri / a) ** 8
                    * (a / (1 - e ** 2) ** (15 / 2))
                    * (f1 - (1 - e ** 2) ** (3 / 2) * f2 * Omega_pri / n)
            )
            de_tides_sec = (
                    -27
                    * F_tid
                    * kT_sec
                    * q1
                    * (1 + q1)
                    * (R_sec / a) ** 8
                    * (e / (1 - e ** 2) ** (13 / 2))
                    * (f3 - (11 / 18) * (1 - e ** 2) ** (3 / 2) * f4 * Omega_sec/n)
            )  # eq. (10) of Hut 1981, 99, 126

            de_tides_pri = (
                    -27
                    * F_tid
                    * kT_pri
                    * q2
                    * (1 + q2)
                    * (R_pri / a) ** 8
                    * (e / (1 - e ** 2) ** (13 / 2))
                    * (f3 - (11 / 18) * (1 - e ** 2) ** (3 / 2) * f4 * Omega_pri/n)
            )

            dOmega_tides_sec = (
                    (3 * F_tid * kT_sec * q1 ** 2 * M_sec * R_sec ** 2 / I_sec)
                    * (R_sec / a) ** 6
                    * n
                    / (1 - e ** 2) ** 6
                    * (f2 - (1 - e ** 2) ** (3 / 2) * f5 * Omega_sec / n)
            )  # eq. (11) of Hut 1981, 99, 126
            dOmega_tides_pri = (
                    (3 * F_tid * kT_pri * q2 ** 2 * M_pri * R_pri ** 2 / I_pri)
                    * (R_pri / a) ** 6
                    * n
                    / (1 - e ** 2) ** 6
                    * (f2 - (1 - e ** 2) ** (3 / 2) * f5 * Omega_pri / n)
            )
            if verbose:
                print("da,de,dOmega_tides = ",
                    da_tides_sec, de_tides_sec, dOmega_tides_sec,
                    da_tides_pri, de_tides_pri, dOmega_tides_pri)
            da = da + da_tides_sec + da_tides_pri
            de = de + de_tides_sec + de_tides_pri
            dOmega_sec = dOmega_sec + dOmega_tides_sec
            dOmega_pri = dOmega_pri + dOmega_tides_pri

        #  Gravitional radiation
        if do_gravitational_radiation:
            v = (M_pri * M_sec / (M_pri + M_sec) ** 2)
            da_gr = (
                (-2 * const.clight / 15) * (v / ((1 - e ** 2) ** (9 / 2)))
                * (const.standard_cgrav * (M_pri + M_sec) * const.msol
                / (a * const.rsol * const.clight ** 2)) ** 3
                * ((96 + 292 * e ** 2 + 37 * e ** 4) * (1 - e ** 2)
                - (1 / 28 * const.standard_cgrav * (M_pri + M_sec) * const.msol
                / (a * const.rsol * const.clight ** 2))
                * ((14008 + 4707 * v)+(80124 + 21560 * v) * e ** 2
                + (17325 + 10458) * e ** 4 - 0.5 * (5501 - 1036 * v) * e ** 6))
                ) * const.secyer / const.rsol
            # eq. (35) of Junker et al. 1992, 254, 146
            de_gr = (
                (-1 / 15) * ((v * const.clight ** 3) / (
                    const.standard_cgrav * (M_pri + M_sec) * const.msol))
                * (const.standard_cgrav * (M_pri + M_sec) * const.msol / (
                    a * const.rsol * const.clight ** 2)) ** 4
                * (e / (1 - e ** 2) ** (7 / 2))
                * ((304 + 121 * e ** 2) * (1 - e ** 2)
                - (1 / 56 * const.standard_cgrav * (M_pri + M_sec) * const.msol
                / (a * const.rsol * const.clight ** 2))
                * (8 * (16705 + 4676 * v) + 12 * (9082 + 2807 * v) * e ** 2
                - (25211 - 3388 * v) * e ** 4))
                ) * const.secyer
            # eq. (36) of Junker et al. 1992, 254, 146
            if verbose:
                print("da,de_gr = ", da_gr, de_gr)
            da = da + da_gr
            de = de + de_gr

        #  Magnetic braking
        if do_magnetic_braking:
            # domega_mb / dt = torque_mb / I is calculated below.
            # All results are in units of [yr^-2], i.e., the amount of change
            # in Omega over 1 year.

            if magnetic_braking_mode == "RVJ83":
                # Torque from Rappaport, Verbunt, and Joss 1983, ApJ, 275, 713
                # The torque is eq.36 of Rapport+1983, with γ = 4. Torque units
                # converted from cgs units to [Msol], [Rsol], [yr] as all stellar
                # parameters are given in units of [Msol], [Rsol], [yr] and so that
                # dOmega_mb/dt is in units of [yr^-2].
                dOmega_mb_sec = (
                    -3.8e-30 * (const.rsol**2 / const.secyer)
                    * M_sec
                    * R_sec**4
                    * Omega_sec**3
                    / I_sec
                    * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
                )
                dOmega_mb_pri = (
                    -3.8e-30 * (const.rsol**2 / const.secyer)
                    * M_pri
                    * R_pri**4
                    * Omega_pri**3
                    / I_pri
                    * np.clip((1.5 - M_pri) / (1.5 - 1.3), 0, 1)
                )
                # Converting units:
                # The constant 3.8e-30 from Rappaport+1983 has units of [cm^-2 s]
                # which need to be converted...
                #
                # -3.8e-30 [cm^-2 s] * (const.rsol**2/const.secyer) -> [Rsol^-2 yr]
                # * M [Msol]
                # * R ** 4 [Rsol^4]
                # * Omega ** 3 [yr^-3]
                # / I [Msol Rsol^2 ]
                #
                # Thus, dOmega/dt comes out to [yr^-2]

            elif magnetic_braking_mode == "M15":

                # Torque prescription from Matt et al. 2015, ApJ, 799, L23
                # Constants:
                # [erg] or [g cm^2 s^-2] -> [Msol Rsol^2 yr^-2]
                K = 1.4e30 * const.secyer**2 / (const.msol * const.rsol**2)
                # m = 0.22
                # p = 2.6
                # Above constants were calibrated as in
                # Gossage et al. 2021, ApJ, 912, 65

                # TODO: I am not sure which constants are used from each reference

                # Below, constants are otherwise as assumed as in
                # Matt et al. 2015, ApJ, 799, L23
                omega_sol = 2.6e-6 * const.secyer   # [s^-1] -> [yr^-1]
                # solar rossby = 2
                # solar convective turnover time = 12.9 days
                # Rossby number saturation threshold = 0.14
                chi = 2.0 / 0.14
                tau_conv_sol = 12.9 / 365.25        # 12.9 [days] -> [yr]

                Prot_pri = 2 * np.pi / Omega_pri    # [yr]
                Rossby_number_pri = Prot_pri / tau_conv_pri
                Prot_sec = 2 * np.pi / Omega_sec    # [yr]
                Rossby_number_sec = Prot_sec / tau_conv_sec

                # critical rotation rate in rad/yr
                Omega_crit_pri = np.sqrt(
                    const.standard_cgrav * M_pri * const.msol
                    / ((R_pri * const.rsol) ** 3)) * const.secyer
                Omega_crit_sec = np.sqrt(
                    const.standard_cgrav * M_sec * const.msol
                    / ((R_sec * const.rsol) ** 3)) * const.secyer

                # omega/omega_c
                wdivwc_pri = Omega_pri / Omega_crit_pri
                wdivwc_sec = Omega_sec / Omega_crit_sec

                gamma_pri = (1 + (wdivwc_pri / 0.072)**2)**0.5
                T0_pri = K * R_pri**3.1 * M_pri**0.5 * gamma_pri**(-2 * 0.22)
                gamma_sec = (1 + (wdivwc_sec / 0.072)**2)**0.5
                T0_sec = K * R_sec**3.1 * M_sec**0.5 * gamma_sec**(-2 * 0.22)

                if (Rossby_number_sec < 0.14):
                    dOmega_mb_sec = (
                        T0_sec * (chi**2.6) * (Omega_sec / omega_sol) / I_sec
                        * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
                    )
                else:
                    dOmega_mb_sec = (
                        T0_sec * ((tau_conv_sec/tau_conv_sol)**2.6)
                            * ((Omega_sec/omega_sol)**(2.6 + 1)) / I_sec
                        * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
                    )

                if (Rossby_number_pri < 0.14):
                    dOmega_mb_pri = (
                        T0_pri * (chi**2.6) * (Omega_pri / 2.6e-6) / I_pri
                        * np.clip((1.5 - M_pri) / (1.5 - 1.3), 0, 1)
                    )
                else:
                    dOmega_mb_pri = (
                        T0_pri * ((tau_conv_pri/tau_conv_sol)**2.6)
                            * ((Omega_pri/omega_sol)**(2.6 + 1)) / I_pri
                        * np.clip((1.5 - M_pri) / (1.5 - 1.3), 0, 1)
                    )

            elif magnetic_braking_mode == "G18":

                # Torque prescription from Garraffo et al. 2018, ApJ, 862, 90
                # a = 0.03
                # b = 0.5
                # [g cm^2] -> [Msol Rsol^2]
                c = 3e41 / (const.msol * const.rsol**2)
                # Above are as calibrated in Gossage et al. 2021, ApJ, 912, 65

                Prot_pri = 2 * np.pi / Omega_pri            # [yr]
                Rossby_number_pri = Prot_pri / tau_conv_pri
                Prot_sec = 2 * np.pi / Omega_sec            # [yr]
                Rossby_number_sec = Prot_sec / tau_conv_sec

                n_pri = (0.03 / Rossby_number_pri) + 0.5 * Rossby_number_pri + 1.0
                n_sec = (0.03 / Rossby_number_sec) + 0.5 * Rossby_number_sec + 1.0

                Qn_pri = 4.05 * np.exp(-1.4 * n_pri)
                Qn_sec = 4.05 * np.exp(-1.4 * n_sec)

                dOmega_mb_sec = (
                        c * Omega_sec**3 * tau_conv_sec * Qn_sec / I_sec
                        * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
                )

                dOmega_mb_pri = (
                        c * Omega_pri**3 * tau_conv_pri * Qn_pri / I_pri
                        * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
                )

            elif magnetic_braking_mode == "CARB":

                # Torque prescription from Van & Ivanova 2019, ApJ, 886, L31
                # Based on files hosted on Zenodo:
                #         https://zenodo.org/record/3647683#.Y_TfedLMKUk,
                # with units converted from [cm], [g], [s] to [Rsol], [Msol], [yr]

                # Constants as assumed in Van & Ivanova 2019, ApJ, 886, L31
                omega_sol = 3e-6 * const.secyer         # [s^-1] -> [yr^-1]
                tau_conv_sol = 2.8e6 / const.secyer     # [s] -> yr
                K2 = 0.07**2

                tau_ratio_sec = tau_conv_sec / tau_conv_sol
                tau_ratio_pri = tau_conv_pri / tau_conv_sol
                rot_ratio_sec = Omega_sec / omega_sol
                rot_ratio_pri = Omega_pri / omega_sol

                # below in units of [Rsol yr^-1]^2
                v_esc2_sec = ((2 * const.standard_cgrav * M_sec / R_sec)
                            * (const.msol * const.secyer**2 / const.rsol**3))
                v_esc2_pri = ((2 * const.standard_cgrav * M_pri / R_pri)
                            * (const.msol * const.secyer**2 / const.rsol**3))
                v_mod2_sec = v_esc2_sec + (2 * Omega_sec**2 * R_sec**2) / K2
                v_mod2_pri = v_esc2_pri + (2 * Omega_pri**2 * R_pri**2) / K2

                # Van & Ivanova 2019, MNRAS 483, 5595 replace the magnetic field
                # with Omega * tau_conv phenomenology. Thus, the ratios
                # (rot_ratio_* and tau_ratio_*) inherently have units of Gauss
                # [cm^-0.5 g^0.5 s^-1] that needs to be converted to [Rsol],
                # [Msol], [yr]. VI2019 assume the solar magnetic field strength is
                # on average 1 Gauss.
                if (abs(Mdot_sec) > 0):
                    R_alfven_div_R3_sec = (
                        R_sec**4 * rot_ratio_sec**4 * tau_ratio_sec**4
                        / (Mdot_sec**2 * v_mod2_sec)
                        * (const.rsol**2 * const.secyer / const.msol**2))
                else:
                    R_alfven_div_R3_sec = 0.0

                if (abs(Mdot_pri) > 0):
                    R_alfven_div_R3_pri = (
                        R_pri**4 * rot_ratio_pri**4 * tau_ratio_pri**4
                        / (Mdot_pri**2 * v_mod2_pri)
                        * (const.rsol**2 * const.secyer / const.msol**2))
                else:
                    R_alfven_div_R3_pri = 0.0

                # Alfven radius in [Rsol]
                R_alfven_sec = R_sec * R_alfven_div_R3_sec**(1./3.)
                R_alfven_pri = R_pri * R_alfven_div_R3_pri**(1./3.)

                dOmega_mb_sec = (
                    (2./3.) * Omega_sec * Mdot_sec * R_alfven_sec**2 / I_sec
                    * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
                )

                dOmega_mb_pri = (
                    (2./3.) * Omega_pri * Mdot_pri * R_alfven_pri**2 / I_pri
                    * np.clip((1.5 - M_sec) / (1.5 - 1.3), 0, 1)
                )

            else:
                Pwarn("WARNING: Magnetic braking is not being calculated in the "
                    "detached step. The given magnetic_braking_mode string ",
                    f"'{magnetic_braking_mode}' does not match the available "
                    "built-in cases. To enable magnetic braking, please set "
                    "magnetc_braking_mode to one of the following strings: "
                    "'RVJ83' for Rappaport, Verbunt, & Joss 1983"
                    "'G18' for Garraffo et al. 2018"
                    "'M15' for Matt et al. 2015"
                    "'CARB' for Van & Ivanova 2019", "UnsupportedModelWarning")

            if verbose and verbose != 1:
                print("magnetic_braking_mode = ", magnetic_braking_mode)
                print("dOmega_mb = ", dOmega_mb_sec, dOmega_mb_pri)
                
            dOmega_sec = dOmega_sec + dOmega_mb_sec
            dOmega_pri = dOmega_pri + dOmega_mb_pri

        if do_stellar_evolution_and_spin_from_winds:
            # Due to the secondary's own evolution, we have:
            # domega_spin/dt = d(Jspin/I)/dt = dJspin/dt * 1/I + Jspin*d(1/I)/dt.
            # These are the two terms calculated below.
            dOmega_spin_wind_sec = (
                    2.0 / 3.0 * R_sec ** 2 * Omega_sec * Mdot_sec / I_sec
            )
            # jshell*Mdot/I : specific angular momentum of a thin sperical shell
            # * mdot  / moment of inertia
            # Omega is in rad/yr here, and R, M in solar (Mdot solar/yr).
            dOmega_deformation_sec = np.min(
                [-Omega_sec * Idot_sec / I_sec, 100]
            )
            # This is term of Jspin*d(1/I)/dt term of the domega_spin/dt. #
            # We limit its increase due to contraction to 100 [(rad/yr)/yr]
            # increase, otherwise the integrator will fail
            # (usually when we have WD formation).
            dOmega_spin_wind_pri = (
                    2.0 / 3.0 * R_pri ** 2 * Omega_pri * Mdot_pri / I_pri
            )
            # jshell*Mdot/I : specific angular momentum of a thin sperical shell
            # * mdot  / moment of inertia
            # Omega is in rad/yr here, and R, M in solar (Mdot solar/yr).
            dOmega_deformation_pri = np.min(
                [-Omega_pri * Idot_pri / I_pri, 100]
            )
            if verbose:
                print(
                    "dOmega_spin_wind , dOmega_deformation = ",
                    dOmega_spin_wind_sec,
                    dOmega_deformation_sec,
                    dOmega_spin_wind_pri,
                    dOmega_deformation_pri,
                )
            dOmega_sec = dOmega_sec + dOmega_spin_wind_sec
            dOmega_pri = dOmega_pri + dOmega_spin_wind_pri
            dOmega_sec = dOmega_sec + dOmega_deformation_sec
            dOmega_pri = dOmega_pri + dOmega_deformation_pri
        if verbose:
            print("a[Ro],e,Omega[rad/yr] have been =", a, e, Omega_sec, Omega_pri)
            print("da,de,dOmega (all) in 1yr is = ",
                da, de, dOmega_sec, dOmega_pri)

        result = [da, de, dOmega_sec, dOmega_pri]

        return result
