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

import time

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

import posydon.utils.constants as const
from posydon.binary_evol.binarystar import BINARYPROPERTIES
from posydon.binary_evol.DT.gravitational_radiation.default_gravrad import (
    default_gravrad,
)
from posydon.binary_evol.DT.key_library import (
    DEFAULT_TRANSLATED_KEYS,
    DEFAULT_TRANSLATION,
)
from posydon.binary_evol.DT.magnetic_braking.prescriptions import (
    CARB_braking,
    G18_braking,
    M15_braking,
    RVJ83_braking,
)
from posydon.binary_evol.DT.tides.default_tides import default_tides
from posydon.binary_evol.DT.track_match import TrackMatcher
from posydon.binary_evol.DT.winds.default_winds import (
    default_sep_from_winds,
    default_spin_from_winds,
)
from posydon.binary_evol.flow_chart import (
    STAR_STATES_CC,
    STAR_STATES_H_RICH_EVOLVABLE,
    STAR_STATES_HE_RICH_EVOLVABLE,
    UNDEFINED_STATES,
)
from posydon.binary_evol.singlestar import STARPROPERTIES
from posydon.config import PATH_TO_POSYDON_DATA

#from posydon.interpolation.data_scaling import DataScaler
from posydon.utils.common_functions import (
    bondi_hoyle,
    check_state_of_star,
    orbital_period_from_separation,
    roche_lobe_radius,
    set_binary_to_failed,
    zero_negative_values,
)
from posydon.utils.posydonerror import (
    ClassificationError,
    FlowError,
    NumericalError,
    POSYDONError,
)
from posydon.utils.posydonwarning import Pwarn


def event(terminal, direction=0):
    """Return a helper function to set attributes for solve_ivp events."""
    def dec(f):
        f.terminal = terminal
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
            matching_tolerance=1e-2,
            matching_tolerance_hard=1e-1,
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
                                          matching_tolerance=matching_tolerance,
                                          matching_tolerance_hard=matching_tolerance_hard,
                                          list_for_matching_HMS = list_for_matching_HMS,
                                          list_for_matching_HeStar = list_for_matching_HeStar,
                                          list_for_matching_postMS = list_for_matching_postMS,
                                          record_matching = record_matching,
                                          verbose = self.verbose)

        # create evolution handler object
        self.init_evo_kwargs()
        self.evo = detached_evolution(**self.evo_kwargs)

        return

    def init_evo_kwargs(self):
        """Store keyword args required to initialize detached evolution, based on step's kwargs."""
        self.evo_kwargs = {
            "primary": None,
            "secondary": None,
            "do_wind_loss": self.do_wind_loss,
            "do_tides": self.do_tides,
            "do_magnetic_braking": self.do_magnetic_braking,
            "magnetic_braking_mode": self.magnetic_braking_mode,
            "do_stellar_evolution_and_spin_from_winds": self.do_stellar_evolution_and_spin_from_winds,
            "do_gravitational_radiation": self.do_gravitational_radiation,
            "verbose": self.verbose,
        }

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

        secondary.t_max = secondary.interp1d.t_max
        primary.t_max = primary.interp1d.t_max
        # set the age offset on the matched track to be the time span
        # from the start of the track to the current age
        # (for these interp1d, x = time)
        secondary.t_offset = binary.time - secondary.interp1d.t0
        secondary.interp1d.offset = secondary.t_offset

        primary.t_offset = binary.time - primary.interp1d.t0
        primary.interp1d.offset = primary.t_offset

        self.max_time = secondary.interp1d.max_time

        # store memory references of primary/secondary
        # for detached evolution
        self.evo.set_stars(primary, secondary, t0 = binary.time)

        if (self.evo.ev_rlo1(binary.time, [binary.separation, binary.eccentricity]) >= 0
            or self.evo.ev_rlo2(binary.time, [binary.separation, binary.eccentricity]) >= 0):
            binary.state = "initial_RLOF"
            return
        else:
            if not (self.max_time - binary.time > 0.0):
                raise ValueError("max_time is lower than the current time. "
                                "Evolution of the detached binary will go to "
                                "lower times.")

            with np.errstate(all="ignore"):

                # solve ODEs for detached evolution
                t_before_ODEsolution = time.time()
                self.res = self.solve_ODEs(binary, primary, secondary)
                t_after_ODEsolution = time.time()

            # clear dictionaries that held current properties during ODE solution
            if hasattr(primary, "latest"):
                del primary.latest
            if hasattr(secondary, "latest"):
                del secondary.latest

            if self.verbose:
                ivp_tspan = t_after_ODEsolution - t_before_ODEsolution
                print(f"\nODE solver duration: {ivp_tspan:.6g} sec")
                print("solution of ODE", self.res)

            if self.res.status == -1:
                failed_state = binary.state
                set_binary_to_failed(binary)
                raise NumericalError(f"Integration failed for {failed_state} binary.")

            # update binary/star properties after detached evolution
            t = self.get_time_after_evo(binary)
            self.update_after_evo(t, binary, primary, secondary)
            self.update_co_stars(t, primary, secondary)

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
            if self.res.t_events[0] or self.res.t_events[1]:

                if self.RLO_orbit_at_orbit_with_same_am:
                    # final circular orbit conserves angular momentum
                    # compared to the eccentric orbit
                    binary.separation *= (1 - self.res.y[1][-1]**2)
                    binary.orbital_period *= (1 - self.res.y[1][-1]**2) ** 1.5
                else:
                    # final circular orbit is at periastron of the ecc. orbit
                    binary.separation *= (1 - self.res.y[1][-1])
                    binary.orbital_period *= (1 - self.res.y[1][-1]) ** 1.5

                abs_diff_porb = np.abs(binary.orbital_period - orbital_period_from_separation(
                                binary.separation, secondary.mass, primary.mass)) / binary.orbital_period


                abs_diff_porb_str = f"\nabs_diff_porb = {abs_diff_porb:.4f}" + \
                    f"\nbinary.orbital_period = {binary.orbital_period:.4f}" +\
                    "\norbital_period_from_separation(binary.separation, secondary.mass, primary.mass) = " + \
                    f"{orbital_period_from_separation(binary.separation, secondary.mass, primary.mass):.4f}"

                if self.verbose:
                    print(abs_diff_porb_str)

                assert abs_diff_porb < 1e-2, \
                        "detached_step: abs_diff_porb >= 1e-2\n" + \
                         abs_diff_porb_str

                # instantly circularize at RLO
                binary.eccentricity = 0

                if self.res.t_events[0]:
                    if secondary == binary.star_1:
                        binary.state = "RLO1"
                        binary.event = "oRLO1"
                    else:
                        binary.state = "RLO2"
                        binary.event = "oRLO2"

                elif self.res.t_events[1]:
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
            elif self.res.t_events[2]:
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

            elif self.res.t_events[3]:
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

    def solve_ODEs(self, binary, primary, secondary):
        """
            Utilizes SciPy's solve_ivp() method to solve a set of
        differential equations that describe the orbital evolution
        (separation and eccentricity) and stellar rotation rate
        evolution during step_detached.

        Parameters
        ----------
        binary : BinaryStar object
            A binary star object, containing the binary system's properties.

        primary : SingleStar object
            A single star object, representing the primary (more evolved) star
            in the binary and containing its properties.

        secondary : SingleStar object
            A single star object, representing the secondary (less evolved) star
            in the binary and containing its properties.

        Returns
        -------
        res : ODESolver object
            This is the ODESolver object produced by SciPy's
            solve_ivp function that contains calculated values
            of the stars evolution through the detached step.

        """
        try:
            res = solve_ivp(self.evo,
                            events=[self.evo.ev_rlo1, self.evo.ev_rlo2,
                                    self.evo.ev_max_time1, self.evo.ev_max_time2],
                            method="Radau",
                            t_span=(binary.time, self.max_time),
                            y0=[binary.separation, binary.eccentricity,
                                secondary.omega0, primary.omega0],
                            dense_output=True)
        except Exception as e:
            if self.verbose:
                print(f"RK45 failed with error {e}, trying with 'RK45' method.")
            res = solve_ivp(self.evo,
                            events=[self.evo.ev_rlo1, self.evo.ev_rlo2,
                                    self.evo.ev_max_time1, self.evo.ev_max_time2],
                            method="RK45",
                            t_span=(binary.time, self.max_time),
                            y0=[binary.separation, binary.eccentricity,
                                secondary.omega0, primary.omega0],
                            dense_output=True)
        return res

    def get_time_after_evo(self, binary):
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
            t = np.arange(binary.time, self.res.t[-1] + self.dt/2.0, self.dt)[1:]
            if t[-1] < self.res.t[-1]:
                t = np.hstack([t, self.res.t[-1]])
        elif (self.n_o_steps_history is not None
                and self.n_o_steps_history > 0):
            t_step = (self.res.t[-1] - binary.time) / self.n_o_steps_history
            t = np.arange(binary.time, self.res.t[-1] + t_step / 2.0, t_step)[1:]
            if t[-1] < self.res.t[-1]:
                t = np.hstack([t, self.res.t[-1]])
        else:  # self.dt is None and self.n_o_steps_history is None
            t = np.array([self.res.t[-1]])

        return t

    def update_after_evo(self, t, binary, primary, secondary):

        """
            Update star and binary properties and interpolators with
        ODESolver result from detached evolution. This update gives
        the binary/stars their appropriate values, according to the
        interpolation after detached evolution.

        Parameters
        ----------
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
        t = np.array(t)

        sep_interp, ecc_interp, omega_interp_sec, omega_interp_pri = self.res.sol(t)

        result_interp_primary = primary.interp1d(t)
        result_interp_secondary = secondary.interp1d(t)

        mass_interp_sec = result_interp_secondary[self.translate["mass"]]
        mass_interp_pri = result_interp_primary[self.translate["mass"]]

        result_interp_secondary["sep"] = sep_interp
        result_interp_secondary["ecc"] = ecc_interp
        result_interp_secondary["time"] = t
        result_interp_secondary["omega"] = omega_interp_sec

        result_interp_primary["sep"] = sep_interp
        result_interp_primary["ecc"] = ecc_interp
        result_interp_primary["time"] = t
        result_interp_primary["omega"] = omega_interp_pri

        p_orb_interp = orbital_period_from_separation(
            sep_interp, mass_interp_sec,
            mass_interp_pri)

        result_interp_primary["porb"] = p_orb_interp
        result_interp_secondary["porb"] = p_orb_interp

        for obj, prop in zip([secondary, primary, binary],
                                [STARPROPERTIES, STARPROPERTIES, BINARYPROPERTIES]):

            # just update orbit and normal stars, COs later
            if obj != binary:
                if obj.co: continue

            interp1d = result_interp_primary if obj == primary else result_interp_secondary

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
                elif (key in ["surf_avg_omega_div_omega_crit"] and obj != binary):
                    # TODO: change `item()` to 0
                    omega_crit_current = np.sqrt(const.standard_cgrav
                        * interp1d[self.translate["mass"]][-1].item() * const.msol
                        / (interp1d[self.translate["R"]][-1].item() * const.rsol)**3)

                    omega_crit_hist = np.sqrt(const.standard_cgrav
                        * interp1d[self.translate["mass"]][:-1] * const.msol
                        / (interp1d[self.translate["R"]][:-1] * const.rsol)**3)

                    current = (interp1d["omega"][-1] / const.secyer / omega_crit_current)
                    history = (interp1d["omega"][:-1] / const.secyer / omega_crit_hist)

                    # ensure positive rotation values
                    current = zero_negative_values([current], key)[0]
                    history = zero_negative_values(history, key)

                elif (key in ["surf_avg_omega"] and obj != binary):

                    current = interp1d["omega"][-1] / const.secyer
                    history = interp1d["omega"][:-1] / const.secyer

                    current = zero_negative_values([current], key)[0]
                    history = zero_negative_values(history, key)

                elif ("rl_relative_overflow_" in key and obj == binary):
                    s = binary.star_1 if "_1" in key[-2:] else binary.star_2
                    s_alt = binary.star_2 if "_1" in key[-2:] else binary.star_1
                    if secondary == s:
                        current = self.evo.ev_rel_rlo1(t[-1], [interp1d["sep"][-1], interp1d["ecc"][-1]])
                        history = self.evo.ev_rel_rlo1(t[:-1], [interp1d["sep"][:-1], interp1d["ecc"][:-1]])

                    elif secondary == s_alt:
                        current = self.evo.ev_rel_rlo2(t[-1], [interp1d["sep"][-1], interp1d["ecc"][-1]])
                        history = self.evo.ev_rel_rlo2(t[:-1], [interp1d["sep"][:-1], interp1d["ecc"][:-1]])

                elif key in ["separation", "orbital_period", "eccentricity", "time"]:
                    current = interp1d[self.translate[key]][-1].item()
                    history = interp1d[self.translate[key]][:-1]

                    current = zero_negative_values([current], key)[0]
                    history = zero_negative_values(history, key)

                elif (key in ["total_moment_of_inertia"] and obj != binary):
                    current = interp1d[self.translate[key]][-1].item() * (const.msol * const.rsol**2)
                    history = interp1d[self.translate[key]][:-1] * (const.msol * const.rsol**2)

                    current = zero_negative_values([current], key)[0]
                    history = zero_negative_values(history, key)

                elif (key in ["log_total_angular_momentum"] and obj != binary):
                    tot_j = (interp1d["omega"][-1] / const.secyer) \
                              * (interp1d[self.translate["total_moment_of_inertia"]][-1].item() \
                              * (const.msol * const.rsol**2))
                    current = np.log10(tot_j) if tot_j > 0.0 else -99

                    tot_j_hist = (interp1d["omega"][:-1] / const.secyer) \
                                   * (interp1d[self.translate["total_moment_of_inertia"]][:-1] \
                                   * (const.msol * const.rsol**2))
                    history = np.where(tot_j_hist > 0, np.log10(tot_j_hist), -99)

                    current = zero_negative_values([current], key)[0]
                    history = zero_negative_values(history, key)

                elif (key in ["spin"] and obj != binary):
                    current = (const.clight
                        * (interp1d["omega"][-1] / const.secyer)
                        * interp1d[self.translate["total_moment_of_inertia"]][-1].item() \
                        * (const.msol * const.rsol**2)
                        / (const.standard_cgrav * (interp1d[self.translate["mass"]][-1].item() \
                        * const.msol)**2))

                    history = (const.clight
                        * (interp1d["omega"][:-1] / const.secyer) \
                        * interp1d[self.translate["total_moment_of_inertia"]][:-1] \
                        * (const.msol * const.rsol**2) \
                        / (const.standard_cgrav * (interp1d[self.translate["mass"]][:-1] \
                        * const.msol)**2))

                    current = zero_negative_values([current], key)[0]
                    history = zero_negative_values(history, key)

                elif (key in ["lg_mdot", "lg_wind_mdot"] and obj != binary):

                    if interp1d[self.translate[key]][-1] == 0:
                        current = -99.0
                    else:
                        current = np.log10(np.abs(interp1d[self.translate[key]][-1])).item()

                    history = np.ones_like(t[:-1])
                    for i in range(len(t)-1):
                        if (interp1d[self.translate[key]][i] == 0):
                            history[i] = -99.0
                        else:
                            history[i] = np.log10(np.abs(interp1d[self.translate[key]][i]))

                elif (self.translate[key] in interp1d and obj != binary):
                    current = interp1d[self.translate[key]][-1].item()
                    history = interp1d[self.translate[key]][:-1]

                elif key in ["profile"]:
                    current = None
                    history = [current] * len(t[:-1])

                else:
                    current = np.nan
                    history = np.ones_like(t[:-1]) * current

                setattr(obj, key, current)
                getattr(obj, key + "_history").extend(history)

    def update_co_stars(self, t, primary, secondary):

        """
            Update compact object properties after detached
        evolution. The properties are updated here using the
        CO star properties from the last step. Often, these
        values are null.

        Parameters
        ----------
        t : float or array[float]
            This is the time elapsed as a result of detached
            evolution in years. This is a float unless the
            user specifies a timestep to use via the simulation
            properties ini file, in which case it is an array.

        primary : SingleStar object
            A single star object, representing the primary (more evolved) star
            in the binary and containing its properties.

        secondary : SingleStar object
            A single star object, representing the secondary (less evolved) star
            in the binary and containing its properties.
        """

        for obj, prop in zip([secondary, primary],
                           [STARPROPERTIES, STARPROPERTIES]):

            # only update compact objects here
            if not obj.co:
                continue

            for key in prop:

                # simply get the current attribute value and update
                # this step's props with it. Detached evolution does not
                # modify these properties for a CO by default, so they
                # typically remain unchanged from the previous step.
                current = getattr(obj, key)
                history = [current] * len(t[:-1])
                setattr(obj, key, current)
                getattr(obj, key + "_history").extend(history)

class detached_evolution:

    def __init__(self, primary=None, secondary=None,
                    do_wind_loss=True,
                    do_tides=True,
                    do_magnetic_braking=True,
                    magnetic_braking_mode="RVJ83",
                    do_stellar_evolution_and_spin_from_winds=True,
                    do_gravitational_radiation=True,
                    verbose=False):

        self.verbose = verbose

        # initially these are typically NoneType
        # and set by set_stars()
        self.primary = primary
        self.secondary = secondary
        # physical properties tracked by interpolators
        self.phys_keys = []

        # also separation and eccentricity
        self.a = np.nan
        self.e = np.nan

        # detached evolution options
        self.do_wind_loss = do_wind_loss
        self.do_tides = do_tides
        self.do_gravitational_radiation = do_gravitational_radiation
        self.do_magnetic_braking = do_magnetic_braking
        self.magnetic_braking_mode = magnetic_braking_mode
        self.do_stellar_evolution_and_spin_from_winds = do_stellar_evolution_and_spin_from_winds

    # timing events for solve_ivp...
    # detects secondary RLO
    @event(True, 1)
    def ev_rlo1(self, t, y):
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
        result_primary = self.primary.interp1d(t)
        result_secondary = self.secondary.interp1d(t)

        pri_mass = result_primary["mass"]
        sec_mass = result_secondary["mass"]

        sep = y[0]
        ecc = y[1]

        RL = roche_lobe_radius(sec_mass, pri_mass, (1 - ecc) * sep)

        # 95% filling of the RL is enough to assume beginning of RLO,
        # as we do in CO-HMS_RLO grid
        RL_diff = result_secondary["R"] - 0.95*RL
        return RL_diff

    # detects primary RLO
    @event(True, 1)
    def ev_rlo2(self, t, y):
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
        result_primary = self.primary.interp1d(t)
        result_secondary = self.secondary.interp1d(t)

        pri_mass = result_primary["mass"]
        sec_mass = result_secondary["mass"]


        sep = y[0]
        ecc = y[1]

        RL = roche_lobe_radius(pri_mass, sec_mass, (1 - ecc) * sep)
        RL_diff = result_primary["R"] - 0.95*RL
        return RL_diff

    # detects secondary RLO via relative difference btwn. R and R_RL
    @event(True, 1)
    def ev_rel_rlo1(self, t, y):
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
        result_primary = self.primary.interp1d(t)
        result_secondary = self.secondary.interp1d(t)

        pri_mass = result_primary["mass"]
        sec_mass = result_secondary["mass"]

        sep = y[0]
        ecc = y[1]

        RL = roche_lobe_radius(sec_mass, pri_mass, (1 - ecc) * sep)
        RL_rel_diff = (result_secondary["R"] - RL) / RL
        return RL_rel_diff

    # detects primary RLO via relative difference btwn. R and R_RL
    @event(True, 1)
    def ev_rel_rlo2(self, t, y):
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
        result_primary = self.primary.interp1d(t)
        result_secondary = self.secondary.interp1d(t)

        pri_mass = result_primary["mass"]
        sec_mass = result_secondary["mass"]

        sep = y[0]
        ecc = y[1]

        RL = roche_lobe_radius(pri_mass, sec_mass, (1 - ecc) * sep)
        RL_rel_diff = (result_primary["R"] - RL) / RL
        return RL_rel_diff

    # detects if the max age in track of secondary is reached
    @event(True, -1)
    def ev_max_time1(self, t, y):
        return self.secondary.t_max - t + self.secondary.t_offset

    # detects if the max age in track of primary is reached
    @event(True, -1)
    def ev_max_time2(self, t, y):
        return self.primary.t_max - t + self.primary.t_offset

    def set_stars(self, primary, secondary, t0=0.0):
        """Sets memory references for primary and secondary star associated with
        this evolution. It is expected that primary/secondary have interp1d
        objects already, as required for detached evolution.

        Parameters
        ----------
        primary : SingleStar object
            A single star object, representing the primary (more evolved) star
            in the binary and containing its properties.

        secondary : SingleStar object
            A single star object, representing the secondary (less evolved) star
            in the binary and containing its properties.

        t0 : float
            The time at the start of detached evolution. Typically should be the
            binary.time prior to detached evolution.

        """

        self.primary = primary
        self.secondary = secondary
        # update list of tracked physical properties
        self.phys_keys = list(secondary.interp1d.keys)
        # except we don't evolve these:
        #for prop in ["t0", "m0", "t_max", "max_time"]:
        #    self.phys_keys.remove(prop)

        # dictionaries to track current properties during evolution
        self.primary.latest = {}
        self.secondary.latest = {}

        self.t = t0

    def update_props(self, t, y):
        """
            Update properties of stars w/ current age during detached evolution.
        This uses the interp1d objects (PchipInterpolator2) of the primary and
        secondary to interpolate stellar properties using the current time of
        integration, t, along a stellar track. Compact objects have no stellar
        track to interpolate along, and simply either have copies of the surviving
        star, or array-like, zeroed arrays.

            The primary and secondary (SingleStar objects) are each given a new
        attribute called `latest` that is a dictionary containing the stellar
        properties of the star at the latest time of integration, t, for access
        in calculations during step_detached evolution.

        Parameters
        ----------
            t : float
                The current time of integration during solve_ivp().

            y : list
                The current set of solutions for orbital separation [Rsol],
            orbital eccentricity, spin of the secondary star [rad/yr], and
            spin of the primary star [rad/yr].
        """

        # updating star properties with interpolated values at current time
        # TODO: update star current/history progressively here rather than after evo?
        primary_data = self.primary.interp1d(t)
        secondary_data = self.secondary.interp1d(t)
        for key in self.phys_keys:
            self.primary.latest[key] = primary_data[key]
            self.secondary.latest[key] = secondary_data[key]

        # update omega, a, e, based on current diffeq solution
        y[0] = np.max([y[0], 0])  # We limit separation to non-negative values
        self.a = y[0]
        y[1] = np.max([y[1], 0])  # We limit eccentricity to non-negative values
        self.e = y[1]
        if self.e > 0 and self.e < 10.0 ** (-3):
            # we force a negligible eccentricity to become 0
            # for computational stability
            self.e = 0.0
            if self.verbose and self.verbose != 1:
                print("Negligible eccentricity became 0 for "
                    "computational stability")

        # update omega
        y[2] = np.max([y[2], 0])  # We limit omega spin to non-negative values
        self.secondary.latest["omega"] = y[2]  # in rad/yr
        y[3] = np.max([y[3], 0])
        self.primary.latest["omega"] = y[3]

        # store current delta(t)/time
        self.t = t

    def __call__(self, t, y):
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

        # update star/orbital props w/ current time during integration
        self.update_props(t, y)
        # initialize deltas for this timestep
        self.da = 0.0
        self.de = 0.0
        self.dOmega_sec = 0.0
        self.dOmega_pri = 0.0

        #  Tidal forces affecting orbit and stellar spins
        if self.do_tides:
            self.tides()
            if self.verbose:
                print(f"After tides: da={self.da:.6e}, de={self.de:.6e}, dO_sec={self.dOmega_sec:.6e}, dO_pri={self.dOmega_pri:.6e}")

        #  Gravitional radiation affecting the orbit
        if self.do_gravitational_radiation:
            self.gravitational_radiation()
            if self.verbose:
                print(f"After GR: da={self.da:.6e}, de={self.de:.6e}")


        #  Magnetic braking affecting stellar spins
        if self.do_magnetic_braking:
            self.magnetic_braking()
            if self.verbose:
                print(f"After MB: dO_sec={self.dOmega_sec:.6e}, dO_pri={self.dOmega_pri:.6e}")

        #  Mass Loss affecting orbital separation
        if self.do_wind_loss:
            self.sep_from_winds()
            if self.verbose:
                print(f"After winds: da={self.da:.6e}")

        # Mass loss affecting stellar spins
        if self.do_stellar_evolution_and_spin_from_winds:
            self.spin_from_winds()
            if self.verbose:
                print(f"After spin winds: dO_sec={self.dOmega_sec:.6e}, dO_pri={self.dOmega_pri:.6e}")

        result = [self.da, self.de, self.dOmega_sec, self.dOmega_pri]

        if self.verbose:
            print(f"Final derivatives: {result}")
            # Check for non-finite values
            if not all(np.isfinite(result)):
                print(f"ERROR: Non-finite derivative detected!")

        return result

    def spin_from_winds(self):

        dOmega_sec_winds, dOmega_pri_winds = default_spin_from_winds(self.a,
                                                                        self.e,
                                                                        self.primary,
                                                                        self.secondary,
                                                                        self.verbose)
        # update spins
        self.dOmega_sec += dOmega_sec_winds
        self.dOmega_pri += dOmega_pri_winds

    def sep_from_winds(self):

        da_winds = default_sep_from_winds(self.a, self.e,
                                            self.primary, self.secondary,
                                            self.verbose)
        # update separation
        self.da += da_winds

    def tides(self):

        da_tides, de_tides, dOmega_sec_tides, dOmega_pri_tides = default_tides(self.a,
                                                                                self.e,
                                                                                self.primary,
                                                                                self.secondary,
                                                                                self.verbose)
        # update orbital params and spin
        self.da += da_tides
        self.de += de_tides
        self.dOmega_sec += dOmega_sec_tides
        self.dOmega_pri += dOmega_pri_tides

    def gravitational_radiation(self):

        da_gr, de_gr = default_gravrad(self.a, self.e,
                                        self.primary, self.secondary,
                                        self.verbose)

        # update orbital params
        self.da += da_gr
        self.de += de_gr

    def magnetic_braking(self):
        # domega_mb / dt = torque_mb / I is calculated below.
        # All results are in units of [yr^-2], i.e., the amount of change
        # in Omega over 1 year.

        if self.magnetic_braking_mode == "RVJ83":
            dOmega_mb_sec, dOmega_mb_pri = RVJ83_braking(self.primary, self.secondary,
                                                            self.verbose)

        elif self.magnetic_braking_mode == "M15":

            dOmega_mb_sec, dOmega_mb_pri = M15_braking(self.primary, self.secondary,
                                                        self.verbose)

        elif self.magnetic_braking_mode == "G18":

            dOmega_mb_sec, dOmega_mb_pri = G18_braking(self.primary, self.secondary,
                                                        self.verbose)

        elif self.magnetic_braking_mode == "CARB":

            dOmega_mb_sec, dOmega_mb_pri = CARB_braking(self.primary, self.secondary,
                                                        self.verbose)

        else:
            Pwarn("WARNING: Magnetic braking is not being calculated in the "
                "detached step. The given magnetic_braking_mode string ",
                f"'{self.magnetic_braking_mode}' does not match the available "
                "built-in cases. To enable magnetic braking, please set "
                "magnetc_braking_mode to one of the following strings: "
                "'RVJ83' for Rappaport, Verbunt, & Joss 1983"
                "'G18' for Garraffo et al. 2018"
                "'M15' for Matt et al. 2015"
                "'CARB' for Van & Ivanova 2019", "UnsupportedModelWarning")

        #if self.verbose:
        #    print("magnetic_braking_mode = ", self.magnetic_braking_mode)
        #    print("dOmega_mb = ", dOmega_mb_sec, dOmega_mb_pri)

        # Sanitise output
        if not np.isfinite(dOmega_mb_sec):
            dOmega_mb_sec = 0.0
        if not np.isfinite(dOmega_mb_pri):
            dOmega_mb_pri = 0.0

        # update spins
        self.dOmega_sec += dOmega_mb_sec
        self.dOmega_pri += dOmega_mb_pri
