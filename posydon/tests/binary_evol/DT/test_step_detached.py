import os
import unittest

from posydon.utils.common_functions import PATH_TO_POSYDON

from posydon.binary_evol.DT.step_detached import detached_step
from posydon.binary_evol.DT.step_detached import diffeq

from posydon.utils import constants as const
from posydon.utils import common_functions as cf

from posydon.binary_evol.binarystar import BinaryStar
from posydon.binary_evol.singlestar import SingleStar
from posydon.binary_evol.simulationproperties import SimulationProperties

PATH_TO_DATA = os.path.join(
    PATH_TO_POSYDON, "posydon/tests/data/POSYDON-UNIT-TESTS/binary_evol/detached/")
#eep_version = "POSYDON"


class TestDetached_step(unittest.TestCase):
    def test_matching1_root(self):
        method = "root"
        matching = detached_step(#grid='POSYDON',
                              path=PATH_TO_DATA,
                              matching_method=method,
                              #eep_version=eep_version,
                              verbose=False)
        get_mist0 = detached_step.get_mist0
        get_track_val = detached_step.get_track_val
        htrack = True
        PROPERTIES_STAR = {
            "mass": 60.0,
            "log_R": 1.0,
            "mdot": -(10.0**(-5)),
            "state": "H-rich_Core_H_burning",
            "center_he4": 0.48,
            "center_h1": 0.5,
            "total_moment_of_inertia": 10.0**57,
            "log_total_angular_momentum": 52,
            "he_core_mass": 0.0,
            "surface_he4": 0.28,
            "surface_h1": 0.7,
        }
        m0, t = get_mist0(matching, SingleStar(**PROPERTIES_STAR),htrack)

        self.assertAlmostEqual(
            m0,
            64.2410914922183,
            places=1,
            msg=
            "Initial mass in MIST matching not exactly what expected. Should be 64.68538051198551",
        )
        self.assertAlmostEqual(
            get_track_val(matching, "mass",htrack, m0, t),
            60.0000000000419,
            places=3,
            msg=
            "Current mass in matching not exactly what expected. Should be 60.0000000000419",
        )
        self.assertAlmostEqual(
            get_track_val(matching, "log_R",htrack, m0, t),
            1.1282247490900794,
            places=1,
            msg=
            "Current log_R in matching not exactly what expected. Should be 1.1282247490900794",
        )
        self.assertAlmostEqual(
            get_track_val(matching, "center_he4",htrack, m0, t),
            0.4861139708172655,
            places=2,
            msg=
            "Current center_he4 in matching not exactly what expected. Should be 0.4861139708172655",
        )
        self.assertAlmostEqual(
            get_track_val(matching, "he_core_mass",htrack, m0, t),
            0.0,
            places=1,
            msg=
            "Current he_core_mass matching not exactly what expected. Should be 0.0",
        )

    def test_matching1_minimize(self):
        method = "minimize"
        matching = detached_step(#grid='POSYDON',
                              path=PATH_TO_DATA,
                              matching_method=method,
                              #eep_version=eep_version,
                              verbose=False)
        get_mist0 = detached_step.get_mist0
        get_track_val = detached_step.get_track_val
        htrack = True
        PROPERTIES_STAR = {
            "mass": 60.0,
            "log_R": 1.0,
            "mdot": -(10.0**(-5)),
            "state": "H-rich_Core_H_burning",
            "center_he4": 0.48,
            "center_h1": 0.5,
            "total_moment_of_inertia": 10.0**57,
            "log_total_angular_momentum": 52,
            "he_core_mass": 0.0,
            "surface_he4": 0.28,
            "surface_h1": 0.7,
        }
        m0, t = get_mist0(matching, SingleStar(**PROPERTIES_STAR),htrack)

        #self.assertAlmostEqual(
        #    m0,
        #    62.88453923015954,
        #    places=
        #    1,  # less accuracy because we try to fit more alternative parameters than "root" method at the same time
        #    msg=
        #    "Initial mass in MIST matching not exactly what expected. Should be 62.78903050084804",
        #)
        #self.assertAlmostEqual(
        #    get_track_val(matching, "mass",htrack, m0, t),
        #    59.96234634155157,
        #    places=1,
        #    msg=
        #    "Current mass in matching not exactly what expected. Should be 59.96234634155157",
        #)
        #self.assertAlmostEqual(
        #    get_track_val(matching, "log_R",htrack, m0, t),
        #    1.0973066851601672,
        #    places=1,
        #    msg=
        #    "Current log_R in matching not exactly what expected. Should be 1.0973066851601672",
        #)
        #self.assertAlmostEqual(
        #    get_track_val(matching, "center_he4",htrack, m0, t),
        #    0.4252496982220549,
        #    places=1,
        #    msg=
        #    "Current center_he4 in matching not exactly what expected. Should be 0.4252496982220549",
        #)
        self.assertAlmostEqual(
            get_track_val(matching, "he_core_mass",htrack, m0, t),
            0.0,
            places=1,
            msg=
            "Current mass he_core_mass matching not exactly what expected. Should be 0.0",
        )

    def test_only_tides(self):
        method = "minimize"
        matching = detached_step(#grid='POSYDON',
                              path=PATH_TO_DATA,
                              matching_method=method,
                              #eep_version=eep_version,
                              verbose=False)
        step_ODE_minimize_hist = detached_step(
            #grid='POSYDON',
            path=PATH_TO_DATA,
            n_o_steps_history=30,
            #eep_version=eep_version,
            matching_method=method,
            verbose=False,
        )
        step_ODE_minimize_hist_onlytides = detached_step(
            #grid='POSYDON',
            path=PATH_TO_DATA,
            n_o_steps_history=30,
            matching_method=method,
            #eep_version=eep_version,
            verbose=False,
            do_wind_loss=False,
            do_tides=True,
            do_gravitational_radiation=False,
            do_magnetic_braking=False,
            do_stellar_evolution_and_spin_from_winds=False,
        )
        PROPERTIES_STAR1 = {"mass": 10.0, "state": "BH"}
        LOW_MS_PROPERTIES_STAR2_non_rot = {
            "mass": 8.0,
            "log_R": 0.6,
            "mdot": -(10.0**(-7)),
            "state": "H-rich_Core_H_burning",
            "center_he4": 0.28,
            "center_h1": 0.7,
            "total_moment_of_inertia": 10.0**57,
            "log_total_angular_momentum": -10.99,  # non-rotating
            "he_core_mass": 0.0,
            "surface_he4": 0.28,
            "surface_h1": 0.7,
        }
        init_orbital_period = 10
        init_separation = cf.orbital_separation_from_period(
            init_orbital_period,
            PROPERTIES_STAR1["mass"],
            LOW_MS_PROPERTIES_STAR2_non_rot["mass"],
        )
        CLOSE_BINARY = {
            "time": 5 * 10.0**6,
            "orbital_period": init_orbital_period,
            "separation": init_separation,
            "state": "detached",
            "eccentricity": 0.0,
            "event": "None",
        }

        binary = BinaryStar(
            star_1=SingleStar(**PROPERTIES_STAR1),
            star_2=SingleStar(**LOW_MS_PROPERTIES_STAR2_non_rot),
            **CLOSE_BINARY)
        binary.properties.max_simulation_time = 10.0**10
        step_ODE_minimize_hist_onlytides(binary)

        self.assertLessEqual(
            getattr(binary, "separation_history")[-1],
            getattr(binary, "separation_history")[0],
            msg=
            "final sepertation with tides only and a non-rotating donor should decrease.",
        )

    def test_tides_vs_tides_and_winds(self):
        method = "minimize"
        matching = detached_step(#grid='POSYDON',
                              path=PATH_TO_DATA,
                              matching_method=method,
                              #eep_version=eep_version,
                              verbose=False)
        step_ODE_minimize_hist_tides_and_winds = detached_step(
            #grid='POSYDON',
            path=PATH_TO_DATA,
            n_o_steps_history=30,
            matching_method=method,
            #eep_version=eep_version,
            verbose=False,
            do_wind_loss=True,
            do_tides=True,
            do_gravitational_radiation=False,
            do_magnetic_braking=False,
            do_stellar_evolution_and_spin_from_winds=False,
        )
        step_ODE_minimize_hist_onlytides = detached_step(
            #grid='POSYDON',
            path=PATH_TO_DATA,
            n_o_steps_history=30,
            matching_method=method,
            #eep_version=eep_version,
            verbose=False,
            do_wind_loss=False,
            do_tides=True,
            do_gravitational_radiation=False,
            do_magnetic_braking=False,
            do_stellar_evolution_and_spin_from_winds=False,
        )

        PROPERTIES_STAR1 = {"mass": 10.0, "state": "BH"}
        LOW_MS_PROPERTIES_STAR2_non_rot = {
            "mass": 8.0,
            "log_R": 0.6,
            "mdot": -(10.0**(-7)),
            "state": "H-rich_Core_H_burning",
            "center_he4": 0.28,
            "center_h1": 0.7,
            "total_moment_of_inertia": 10.0**57,
            "log_total_angular_momentum": -10.99,  # non-rotating
            "he_core_mass": 0.0,
            "surface_he4": 0.28,
            "surface_h1": 0.7,
        }
        init_orbital_period = 10
        init_separation = cf.orbital_separation_from_period(
            init_orbital_period,
            PROPERTIES_STAR1["mass"],
            LOW_MS_PROPERTIES_STAR2_non_rot["mass"],
        )
        CLOSE_BINARY = {
            "time": 5 * 10.0**6,
            "orbital_period": init_orbital_period,
            "separation": init_separation,
            "state": "detached",
            "eccentricity": 0.0,
            "event": "None",
        }

        binary = BinaryStar(
            star_1=SingleStar(**PROPERTIES_STAR1),
            star_2=SingleStar(**LOW_MS_PROPERTIES_STAR2_non_rot),
            **CLOSE_BINARY)
        binary_test = BinaryStar(
            star_1=SingleStar(**PROPERTIES_STAR1),
            star_2=SingleStar(**LOW_MS_PROPERTIES_STAR2_non_rot),
            **CLOSE_BINARY)

        binary.properties.max_simulation_time = 10.0**10
        binary_test.properties.max_simulation_time = 10.0**10

        step_ODE_minimize_hist_tides_and_winds(binary)
        step_ODE_minimize_hist_onlytides(binary_test)

        self.assertLessEqual(
            getattr(binary_test, "separation_history")[-1],
            getattr(binary, "separation_history")[-1],
            msg=
            "final sepertation with tides only and a non-rotating donor should be lower than including winds too that widen the orbit too.",
        )

    # the following tests are out because they need more EEPS MIST models around their mass. If included they should work.
    """
    def test_matching2_root(self):
        method = "root"
        matching = HMS_detached_step(PATH_TO_EEPS, matching_method=method, verbose=True)
        get_mist0 = HMS_detached_step.get_mist0
        get_track_val = HMS_detached_step.get_track_val
        PROPERTIES_STAR = {
            "mass": 20.0,
            "log_R": 2.5,
            "mdot": -(10.0 ** (-5)),
            "state": "PostMS",
            "center_he4": 0.8,
            "center_h1": 0.0,
            "total_moment_of_inertia": 10.0 ** 57,
            "log_total_angular_momentum": 52,
            "he_core_mass": 7.0,
            "surface_he4": 0.2,
            "surface_h1": 0.7,
        }
        m0, t = get_mist0(matching, SingleStar(**PROPERTIES_STAR))

        self.assertAlmostEqual(
            m0,
            23.092430444226363,
            places=5,
            msg="Initial mass in MIST matching not exactly what expected. Should be 23.092430444226363",
        )
        self.assertAlmostEqual(
            get_track_val(matching, "mass", m0, t),
            20.000000000003112,
            places=5,
            msg="Current mass in matching not exactly what expected. Should be 20.000000000003112",
        )
        self.assertAlmostEqual(
            get_track_val(matching, "log_R", m0, t),
            3.006026700371161,
            places=5,
            msg="Current log_R in matching not exactly what expected. Should be 3.006026700371161",
        )
        self.assertAlmostEqual(
            get_track_val(matching, "center_he4", m0, t),
            0.6426242107646186,
            places=5,
            msg="Current center_he4 in matching not exactly what expected. Should be 0.6426242107646186",
        )
        self.assertAlmostEqual(
            get_track_val(matching, "he_core_mass", m0, t),
            6.99999999999913,
            places=5,
            msg="Current mass he_core_mass matching not exactly what expected. Should be 6.99999999999913",
        )

    def test_matching2_minimize(self):
        method = "minimize"
        matching = HMS_detached_step(PATH_TO_EEPS, matching_method=method, verbose=True)
        get_mist0 = HMS_detached_step.get_mist0
        get_track_val = HMS_detached_step.get_track_val
        PROPERTIES_STAR = {
            "mass": 20.0,
            "log_R": 2.5,
            "mdot": -(10.0 ** (-5)),
            "state": "PostMS",
            "center_he4": 0.8,
            "center_h1": 0.0,
            "total_moment_of_inertia": 10.0 ** 57,
            "log_total_angular_momentum": 52,
            "he_core_mass": 7.0,
            "surface_he4": 0.2,
            "surface_h1": 0.7,
        }
        m0, t = get_mist0(matching, SingleStar(**PROPERTIES_STAR))

        self.assertAlmostEqual(
            m0,
            23.441360274390483,
            places=5,
            msg="Initial mass in MIST matching not exactly what expected. Should be 23.441360274390483",
        )
        self.assertAlmostEqual(
            get_track_val(matching, "mass", m0, t),
            21.931950016322705,
            places=5,
            msg="Current mass in matching not exactly what expected. Should be 21.931950016322705",
        )
        self.assertAlmostEqual(
            get_track_val(matching, "log_R", m0, t),
            2.5019520817980974,
            places=5,
            msg="Current log_R in matching not exactly what expected. Should be 2.5019520817980974",
        )
        self.assertAlmostEqual(
            get_track_val(matching, "center_he4", m0, t),
            0.9095107795447867,
            places=5,
            msg="Current center_he4 in matching not exactly what expected. Should be 0.9095107795447867",
        )
        self.assertAlmostEqual(
            get_track_val(matching, "he_core_mass", m0, t),
            6.86182071726521,
            places=5,
            msg="Current mass he_core_mass matching not exactly what expected. Should be 6.86182071726521",
        )
    """


if __name__ == "__main__":
    unittest.main()
