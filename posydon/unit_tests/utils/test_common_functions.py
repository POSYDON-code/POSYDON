"""Unit tests of posydon/utils/common_functions.py
"""

__authors__ = [
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>"
]

# import the module which will be tested
import posydon.utils.common_functions as totest

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns, approx
from inspect import isroutine, isclass
from scipy.interpolate import interp1d
from posydon.binary_evol.binarystar import BinaryStar as BinaryStar_class
from posydon.binary_evol.singlestar import SingleStar as SingleStar_class
from posydon.utils.posydonwarning import (EvolutionWarning,\
    InappropriateValueWarning, ApproximationWarning, InterpolationWarning,\
    ReplaceValueWarning, ClassificationWarning)

@fixture
def BinaryStar():
    # initialize a BinaryStar instance, which is a required argument
    return BinaryStar_class()

@fixture
def SingleStar():
    # initialize a SingleStar instance, which is a required argument
    return SingleStar_class()

@fixture
def star_profile():
    # generate a profile of a test star
    profile = totest.np.empty((4,), dtype=[(f, 'f8') for f in ['mass', 'dm',\
              'radius', 'log_R', 'energy', 'x_mass_fraction_H',\
              'y_mass_fraction_He', 'z_mass_fraction_metals',\
              'neutral_fraction_H', 'neutral_fraction_He', 'avg_charge_He']])
    profile['mass'] = totest.np.array([1.0, 0.5, 0.1, 0.001])
    # get dm as differences from mass and the last mass entry (like next would
    # be 0)
    profile['dm'][:-1] = profile['mass'][:-1] - profile['mass'][1:]
    profile['dm'][-1] = profile['mass'][-1]
    profile['radius'] = totest.np.array([1.0, 0.5, 0.1, 0.001])
    # get log10 of radius
    profile['log_R'] = totest.np.log10(profile['radius'])
    profile['energy'] = totest.np.array([1.0, 0.5, 0.1, 0.001])*1.0e+16
    profile['x_mass_fraction_H'] = totest.np.array([0.7, 0.5, 0.1, 0.001])
    profile['y_mass_fraction_He'] = totest.np.array([0.2, 0.5, 0.8, 0.2])
    # get Z=1-X-Y
    profile['z_mass_fraction_metals'] = 1.0-profile['x_mass_fraction_H']\
                                           -profile['y_mass_fraction_He']
    profile['neutral_fraction_H'] = totest.np.array([1.0, 0.5, 0.1, 0.0])
    profile['neutral_fraction_He'] = totest.np.array([1.0, 0.5, 0.1, 0.0])
    profile['avg_charge_He'] = totest.np.array([0.0, 0.5, 1.1, 1.8])
    return profile

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['ALL_RLO_CASES', 'ALL_STAR_STATES', 'BURNING_STATES',\
                    'CEE_parameters_from_core_abundance_thresholds',\
                    'COMPACT_OBJECTS', 'CO_radius',\
                    'DEFAULT_CE_OPTION_FOR_LAMBDA', 'He_MS_lifetime',\
                    'LG_MTRANSFER_RATE_THRESHOLD', 'LOG10_BURNING_THRESHOLD',\
                    'MT_CASE_A', 'MT_CASE_B', 'MT_CASE_BA', 'MT_CASE_BB',\
                    'MT_CASE_BC', 'MT_CASE_C', 'MT_CASE_NONBURNING',\
                    'MT_CASE_NO_RLO', 'MT_CASE_TO_STR',\
                    'MT_CASE_UNDETERMINED', 'MT_STR_TO_CASE',\
                    'PATH_TO_POSYDON', 'PchipInterpolator',\
                    'PchipInterpolator2', 'Pwarn',\
                    'REL_LOG10_BURNING_THRESHOLD', 'RICHNESS_STATES',\
                    'RL_RELATIVE_OVERFLOW_THRESHOLD',\
                    'STATE_NS_STARMASS_UPPER_LIMIT', 'STATE_UNDETERMINED',\
                    'Schwarzschild_Radius', 'THRESHOLD_CENTRAL_ABUNDANCE',\
                    'THRESHOLD_HE_NAKED_ABUNDANCE', '__authors__',\
                    '__builtins__', '__cached__', '__doc__', '__file__',\
                    '__loader__', '__name__', '__package__', '__spec__',\
                    'beaming', 'bondi_hoyle',\
                    'calculate_H2recombination_energy',\
                    'calculate_Mejected_for_integrated_binding_energy',\
                    'calculate_Patton20_values_at_He_depl',\
                    'calculate_binding_energy', 'calculate_core_boundary',\
                    'calculate_lambda_from_profile',\
                    'calculate_recombination_energy', 'check_state_of_star',\
                    'check_state_of_star_history_array', 'const',\
                    'convert_metallicity_to_string', 'copy',\
                    'cumulative_mass_transfer_flag',\
                    'cumulative_mass_transfer_numeric',\
                    'cumulative_mass_transfer_string', 'eddington_limit',\
                    'flip_stars', 'get_binary_state_and_event_and_mt_case',\
                    'get_binary_state_and_event_and_mt_case_array',\
                    'get_i_He_depl', 'get_internal_energy_from_profile',\
                    'get_mass_radius_dm_from_profile', 'histogram_sampler',\
                    'infer_mass_transfer_case', 'infer_star_state',\
                    'initialize_empty_array',\
                    'inspiral_timescale_from_orbital_period',\
                    'inspiral_timescale_from_separation', 'interp1d',\
                    'inverse_sampler', 'is_number',\
                    'linear_interpolation_between_two_cells', 'newton', 'np',\
                    'orbital_period_from_separation',\
                    'orbital_separation_from_period', 'os', 'pd',\
                    'period_change_stabe_MT', 'period_evol_wind_loss',\
                    'profile_recomb_energy', 'quad',\
                    'read_histogram_from_file', 'rejection_sampler',\
                    'roche_lobe_radius', 'rotate', 'rzams',\
                    'separation_evol_wind_loss', 'set_binary_to_failed',\
                    'spin_stable_mass_transfer', 'stefan_boltzmann_law']
        assert dir(totest) == elements, "There might be added or removed "+\
               "objects without an update on the unit test."

    def test_instance_PATH_TO_POSYDON(self):
        assert isinstance(totest.PATH_TO_POSYDON, str)

    def test_instance_STATE_UNDETERMINED(self):
        assert isinstance(totest.STATE_UNDETERMINED, str)

    def test_instance_BURNING_STATES(self):
        assert isinstance(totest.BURNING_STATES, list)

    def test_instance_RICHNESS_STATES(self):
        assert isinstance(totest.RICHNESS_STATES, list)

    def test_instance_COMPACT_OBJECTS(self):
        assert isinstance(totest.COMPACT_OBJECTS, list)

    def test_instance_ALL_STAR_STATES(self):
        assert isinstance(totest.ALL_STAR_STATES, list)

    def test_instance_MT_CASE_NO_RLO(self):
        assert isinstance(totest.MT_CASE_NO_RLO, int)

    def test_instance_MT_CASE_A(self):
        assert isinstance(totest.MT_CASE_A, int)

    def test_instance_MT_CASE_B(self):
        assert isinstance(totest.MT_CASE_B, int)

    def test_instance_MT_CASE_C(self):
        assert isinstance(totest.MT_CASE_C, int)

    def test_instance_MT_CASE_BA(self):
        assert isinstance(totest.MT_CASE_BA, int)

    def test_instance_MT_CASE_BB(self):
        assert isinstance(totest.MT_CASE_BB, int)

    def test_instance_MT_CASE_BC(self):
        assert isinstance(totest.MT_CASE_BC, int)

    def test_instance_MT_CASE_NONBURNING(self):
        assert isinstance(totest.MT_CASE_NONBURNING, int)

    def test_instance_MT_CASE_UNDETERMINED(self):
        assert isinstance(totest.MT_CASE_UNDETERMINED, int)

    def test_instance_ALL_RLO_CASES(self):
        assert isinstance(totest.ALL_RLO_CASES, set)

    def test_instance_MT_CASE_TO_STR(self):
        assert isinstance(totest.MT_CASE_TO_STR, dict)

    def test_instance_MT_STR_TO_CASE(self):
        assert isinstance(totest.MT_STR_TO_CASE, dict)

    def test_instance_DEFAULT_CE_OPTION_FOR_LAMBDA(self):
        assert isinstance(totest.DEFAULT_CE_OPTION_FOR_LAMBDA, str)

    def test_instance_is_number(self):
        assert isroutine(totest.is_number)

    def test_instance_stefan_boltzmann_law(self):
        assert isroutine(totest.stefan_boltzmann_law)

    def test_instance_rzams(self):
        assert isroutine(totest.rzams)

    def test_instance_roche_lobe_radius(self):
        assert isroutine(totest.roche_lobe_radius)

    def test_instance_orbital_separation_from_period(self):
        assert isroutine(totest.orbital_separation_from_period)

    def test_instance_orbital_period_from_separation(self):
        assert isroutine(totest.orbital_period_from_separation)

    def test_instance_eddington_limit(self):
        assert isroutine(totest.eddington_limit)

    def test_instance_beaming(self):
        assert isroutine(totest.beaming)

    def test_instance_bondi_hoyle(self):
        assert isroutine(totest.bondi_hoyle)

    def test_instance_rejection_sampler(self):
        assert isroutine(totest.rejection_sampler)

    def test_instance_inverse_sampler(self):
        assert isroutine(totest.inverse_sampler)

    def test_instance_histogram_sampler(self):
        assert isroutine(totest.histogram_sampler)

    def test_instance_read_histogram_from_file(self):
        assert isroutine(totest.read_histogram_from_file)

    def test_instance_inspiral_timescale_from_separation(self):
        assert isroutine(totest.inspiral_timescale_from_separation)

    def test_instance_inspiral_timescale_from_orbital_period(self):
        assert isroutine(totest.inspiral_timescale_from_orbital_period)

    def test_instance_spin_stable_mass_transfer(self):
        assert isroutine(totest.spin_stable_mass_transfer)

    def test_instance_check_state_of_star(self):
        assert isroutine(totest.check_state_of_star)

    def test_instance_check_state_of_star_history_array(self):
        assert isroutine(totest.check_state_of_star_history_array)

    def test_instance_get_binary_state_and_event_and_mt_case(self):
        assert isroutine(totest.get_binary_state_and_event_and_mt_case)

    def test_instance_get_binary_state_and_event_and_mt_case_array(self):
        assert isroutine(totest.get_binary_state_and_event_and_mt_case_array)

    def test_instance_CO_radius(self):
        assert isroutine(totest.CO_radius)

    def test_instance_He_MS_lifetime(self):
        assert isroutine(totest.He_MS_lifetime)

    def test_instance_Schwarzschild_Radius(self):
        assert isroutine(totest.Schwarzschild_Radius)

    def test_instance_flip_stars(self):
        assert isroutine(totest.flip_stars)

    def test_instance_set_binary_to_failed(self):
        assert isroutine(totest.set_binary_to_failed)

    def test_instance_infer_star_state(self):
        assert isroutine(totest.infer_star_state)

    def test_instance_infer_mass_transfer_case(self):
        assert isroutine(totest.infer_mass_transfer_case)

    def test_instance_cumulative_mass_transfer_numeric(self):
        assert isroutine(totest.cumulative_mass_transfer_numeric)

    def test_instance_cumulative_mass_transfer_string(self):
        assert isroutine(totest.cumulative_mass_transfer_string)

    def test_instance_cumulative_mass_transfer_flag(self):
        assert isroutine(totest.cumulative_mass_transfer_flag)

    def test_instance_get_i_He_depl(self):
        assert isroutine(totest.get_i_He_depl)

    def test_instance_calculate_Patton20_values_at_He_depl(self):
        assert isroutine(totest.calculate_Patton20_values_at_He_depl)

    def test_instance_CEE_parameters_from_core_abundance_thresholds(self):
        assert isroutine(totest.CEE_parameters_from_core_abundance_thresholds)

    def test_instance_initialize_empty_array(self):
        assert isroutine(totest.initialize_empty_array)

    def test_instance_calculate_core_boundary(self):
        assert isroutine(totest.calculate_core_boundary)

    def test_instance_period_evol_wind_loss(self):
        assert isroutine(totest.period_evol_wind_loss)

    def test_instance_separation_evol_wind_loss(self):
        assert isroutine(totest.separation_evol_wind_loss)

    def test_instance_period_change_stabe_MT(self):
        assert isroutine(totest.period_change_stabe_MT)

    def test_instance_linear_interpolation_between_two_cells(self):
        assert isroutine(totest.linear_interpolation_between_two_cells)

    def test_instance_calculate_lambda_from_profile(self):
        assert isroutine(totest.calculate_lambda_from_profile)

    def test_instance_get_mass_radius_dm_from_profile(self):
        assert isroutine(totest.get_mass_radius_dm_from_profile)

    def test_instance_get_internal_energy_from_profile(self):
        assert isroutine(totest.get_internal_energy_from_profile)

    def test_instance_calculate_H2recombination_energy(self):
        assert isroutine(totest.calculate_H2recombination_energy)

    def test_instance_calculate_recombination_energy(self):
        assert isroutine(totest.calculate_recombination_energy)

    def test_instance_profile_recomb_energy(self):
        assert isroutine(totest.profile_recomb_energy)

    def test_instance_calculate_binding_energy(self):
        assert isroutine(totest.calculate_binding_energy)

    def test_instance_calculate_Mejected_for_integrated_binding_energy(self):
        assert isroutine(totest.\
                         calculate_Mejected_for_integrated_binding_energy)

    def test_instance_PchipInterpolator2(self):
        assert isclass(totest.PchipInterpolator2)

    def test_instance_convert_metallicity_to_string(self):
        assert isroutine(totest.convert_metallicity_to_string)

    def test_instance_rotate(self):
        assert isroutine(totest.rotate)


class TestValues:
    # check that the values fit
    def test_value_PATH_TO_POSYDON(self):
        assert '/' in totest.PATH_TO_POSYDON

    def test_value_STATE_UNDETERMINED(self):
        assert totest.STATE_UNDETERMINED == "undetermined_evolutionary_state"

    def test_value_BURNING_STATES(self):
        for v in ["Core_H_burning", "Core_He_burning", "Shell_H_burning",\
                  "Central_He_depleted", "Central_C_depletion"]:
            # check required values
            assert v in totest.BURNING_STATES, "missing entry"

    def test_value_RICHNESS_STATES(self):
        for v in ["H-rich", "stripped_He", "accreted_He"]:
            # check required values
            assert v in totest.RICHNESS_STATES, "missing entry"

    def test_value_COMPACT_OBJECTS(self):
        for v in ["WD", "NS", "BH","massless_remnant"]:
            # check required values
            assert v in totest.COMPACT_OBJECTS, "missing entry"

    def test_value_ALL_STAR_STATES(self):
        for v in totest.COMPACT_OBJECTS:
            # check required values of COMPACT_OBJECTS
            assert v in totest.ALL_STAR_STATES, "missing entry"
        # check required STATE_UNDETERMINED
        assert totest.STATE_UNDETERMINED in totest.ALL_STAR_STATES,\
               "missing entry"
        for v1 in totest.RICHNESS_STATES:
            for v2 in totest.BURNING_STATES:
                # check required values of combinatios of RICHNESS_STATES and
                # BURNING_STATES
                v = v1+"_"+v2
                assert v in totest.ALL_STAR_STATES, "missing entry"

    def test_value_MT_CASE_NO_RLO(self):
        assert totest.MT_CASE_NO_RLO == 0

    def test_value_MT_CASE_A(self):
        assert totest.MT_CASE_A == 1

    def test_value_MT_CASE_B(self):
        assert totest.MT_CASE_B == 2

    def test_value_MT_CASE_C(self):
        assert totest.MT_CASE_C == 3

    def test_value_MT_CASE_BA(self):
        assert totest.MT_CASE_BA == 4

    def test_value_MT_CASE_BB(self):
        assert totest.MT_CASE_BB == 5

    def test_value_MT_CASE_BC(self):
        assert totest.MT_CASE_BC == 6

    def test_value_MT_CASE_NONBURNING(self):
        assert totest.MT_CASE_NONBURNING == 8

    def test_value_MT_CASE_UNDETERMINED(self):
        assert totest.MT_CASE_UNDETERMINED == 9

    def test_value_ALL_RLO_CASES(self):
        for v in [totest.MT_CASE_A, totest.MT_CASE_B, totest.MT_CASE_C,\
                  totest.MT_CASE_BA, totest.MT_CASE_BB, totest.MT_CASE_BC,\
                  totest.MT_CASE_NONBURNING]:
            # check required values
            assert v in totest.ALL_RLO_CASES, "missing entry"

    def test_value_MT_CASE_TO_STR(self):
        for v in [totest.MT_CASE_NO_RLO, totest.MT_CASE_A, totest.MT_CASE_B,\
                  totest.MT_CASE_C, totest.MT_CASE_BA, totest.MT_CASE_BB,\
                  totest.MT_CASE_BC, totest.MT_CASE_NONBURNING,\
                  totest.MT_CASE_UNDETERMINED]:
            # check required values
            assert v in totest.MT_CASE_TO_STR.keys(), "missing entry"

    def test_value_MT_STR_TO_CASE(self):
        for k,v in totest.MT_STR_TO_CASE.items():
            # check required items
            assert v in totest.MT_CASE_TO_STR.keys()

    def test_value_DEFAULT_CE_OPTION_FOR_LAMBDA(self):
        assert totest.DEFAULT_CE_OPTION_FOR_LAMBDA == "lambda_from_profile_"+\
               "gravitational_plus_internal_minus_recombination"


class TestFunctions:
    @fixture
    def csv_path(self, tmp_path):
        # a temporary csv path for testing
        return totest.os.path.join(tmp_path, "hist.csv")

    @fixture
    def csv_file(self, csv_path):
        # a temporary csv file for testing
        with open(csv_path, "w") as test_file:
            test_file.write("0.0,1.0\n")
            test_file.write("1.0,1.0\n\n")
            test_file.write("2.0,1.0\n")
        return

    @fixture
    def csv_path2(self, tmp_path):
        # a temporary csv path for testing
        return totest.os.path.join(tmp_path, "hist2.csv")

    @fixture
    def csv_file2(self, csv_path2):
        # a temporary csv file for testing
        with open(csv_path2, "w") as test_file:
            test_file.write("#0.0,1.0\n")
            test_file.write("0.0,1.0,2.0\n")
            test_file.write("1.0,1.0,1.0\n")
        return

    # test functions
    def test_is_number(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional argument"):
            totest.is_number()
        # try a string of a float
        assert totest.is_number("1.2e-3")
        # try a non-float string
        assert totest.is_number("test")==False

    def test_stefan_boltzmann_law(self):
        # missing argument
        with raises(TypeError, match="missing 2 required positional"+\
                    " arguments: 'L' and 'R'"):
            totest.stefan_boltzmann_law()
        # examples
        assert totest.stefan_boltzmann_law(1.0, 1.0) ==\
               approx(5.77603658298e+3, abs=6e-9)
        assert totest.stefan_boltzmann_law(2.0, 1.0) ==\
               approx(6.86890380099e+3, abs=6e-9)
        assert totest.stefan_boltzmann_law(1.0, 2.0) ==\
               approx(4.08427463621e+3, abs=6e-9)

    def test_rzams(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'm'"):
            totest.rzams()
        # examples
        assert totest.rzams(1.0) ==\
               approx(0.88824945030, abs=6e-12)
        assert totest.rzams(2.0) ==\
               approx(1.61021038543, abs=6e-12)
        assert totest.rzams(1.0, z=1.0) ==\
               approx(0.85822941705, abs=6e-12)
        assert totest.rzams(1.0, Zsun=1.0) ==\
               approx(0.84963691291, abs=6e-12)

    def test_roche_lobe_radius(self, capsys):
        # missing argument
        with raises(TypeError, match="missing 2 required positional"+\
                    " arguments: 'm1' and 'm2'"):
            totest.roche_lobe_radius()
        # bad input
        for a in [-1.0, totest.np.array([]), totest.np.array([-1.0])]:
            with warns(EvolutionWarning, match="Trying to compute RL radius"+\
                       " for binary with invalid separation"):
                assert totest.np.isnan(totest.roche_lobe_radius(1.0, 1.0,\
                       a_orb=a))
        for m in [0.0, totest.np.array([]), totest.np.array([0.0])]:
            with warns(EvolutionWarning, match="Trying to compute RL radius"+\
                       " for nonexistent object"):
                assert totest.np.isnan(totest.roche_lobe_radius(m, 1.0))
        for m in [0.0, totest.np.array([]), totest.np.array([0.0])]:
            with warns(EvolutionWarning, match="Trying to compute RL radius"+\
                       " for nonexistent companion"):
                assert totest.np.isnan(totest.roche_lobe_radius(1.0, m))
        # examples
        assert totest.roche_lobe_radius(1.0, 1.0) ==\
               approx(0.37892051838, abs=6e-12)
        assert totest.roche_lobe_radius(2.0, 1.0) ==\
               approx(0.44000423753, abs=6e-12)
        assert totest.roche_lobe_radius(1.0, 2.0) ==\
               approx(0.32078812033, abs=6e-12)
        assert totest.roche_lobe_radius(1.0, 1.0, a_orb=2.0) ==\
               approx(0.75784103676, abs=6e-12)
        assert totest.np.allclose(totest.roche_lobe_radius(totest.np.array([1.0, 2.0, 1.0, 1.0]), totest.np.array([1.0, 1.0, 2.0, 1.0]), a_orb=totest.np.array([1.0, 1.0, 1.0, 2.0])), totest.np.array([0.37892051838, 0.44000423753, 0.32078812033, 0.75784103676]))
        # check that roche lobe sum never exceeds orbital separation
        for m in [1.0e+1, 1.0e+2, 1.0e+3, 1.0e+4, 1.0e+5, 1.0e+6, 1.0e+7]:
            assert totest.roche_lobe_radius(m, 1.0)+\
                   totest.roche_lobe_radius(1.0, m) < 1.0

    def test_orbital_separation_from_period(self):
        # missing argument
        with raises(TypeError, match="missing 3 required positional"+\
                    " arguments: 'period_days', 'm1_solar', and 'm2_solar'"):
            totest.orbital_separation_from_period()
        # examples
        assert totest.orbital_separation_from_period(1.0, 15.0, 30.0) ==\
               approx(14.9643417735, abs=6e-11)
        assert totest.orbital_separation_from_period(2.0, 15.0, 30.0) ==\
               approx(23.7544118733, abs=6e-11)
        assert totest.orbital_separation_from_period(1.0, 30.0, 30.0) ==\
               approx(16.4703892879, abs=6e-11)
        assert totest.orbital_separation_from_period(1.0, 15.0, 60.0) ==\
               approx(17.7421890201, abs=6e-11)

    def test_orbital_period_from_separation(self):
        # missing argument
        with raises(TypeError, match="missing 3 required positional"+\
                    " arguments: 'separation', 'm1', and 'm2'"):
            totest.orbital_period_from_separation()
        # examples
        assert totest.orbital_period_from_separation(1.0, 15.0, 30.0) ==\
               approx(1.72773763511e-2, abs=6e-14)
        assert totest.orbital_period_from_separation(2.0, 15.0, 30.0) ==\
               approx(4.88677999159e-2, abs=6e-14)
        assert totest.orbital_period_from_separation(1.0, 30.0, 30.0) ==\
               approx(1.49626468308e-2, abs=6e-14)
        assert totest.orbital_period_from_separation(1.0, 15.0, 60.0) ==\
               approx(1.33829981748e-2, abs=6e-14)

    def test_eddington_limit(self, BinaryStar):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'binary'"):
            totest.eddington_limit()
        # bad input
        with raises(ValueError, match="Eddington limit is being calculated"+\
                    " for a non-CO"):
            totest.eddington_limit(BinaryStar)
        with raises(IndexError, match="index 2 is out of bounds"):
            BinaryStar.star_1.state = 'WD'
            totest.eddington_limit(BinaryStar, idx=2)
        # examples: 1Msun accretor is star1
        BinaryStar.star_1.mass = 1.0
        for CO, r in zip(['WD', 'NS', 'BH'],\
            [(approx(4.01681147088e-17, abs=6e-29),\
              approx(6.40623627946e+7, abs=6e-5)),\
             (approx(2.17747384546e-8, abs=6e-20),\
              approx(0.11817658993, abs=6e-12)),\
             (approx(4.49942509871e-8, abs=6e-20),\
              approx(5.71909584179e-2, abs=6e-14))]):
            BinaryStar.star_1.state = CO
            assert totest.eddington_limit(BinaryStar) == r
        BinaryStar.star_1.state = None
        # examples: 1.2Msun accretor is star2, while donor has X_surf=0.1
        BinaryStar.star_2.mass = 1.2
        for CO, r in zip(['WD', 'NS', 'BH'],\
            [(approx(5.89502616630e-17, abs=6e-29),\
              approx(8.16917025429e+7, abs=6e-5)),\
             (approx(3.39586943808e-8, abs=6e-20),\
              approx(0.14181190792, abs=6e-12)),\
             (approx(8.42046955292e-8, abs=6e-20),\
              approx(5.71909584179e-2, abs=6e-14))]):
            BinaryStar.star_2.state = CO
            BinaryStar.star_1.surface_h1 = 0.1
            assert totest.eddington_limit(BinaryStar) == r
        # examples: 1.3Msun accretor is star1 with history
        BinaryStar.star_1.state = 'BH'
        BinaryStar.star_1.mass_history = totest.np.array([1.3, 1.3])
        BinaryStar.star_1.state_history = totest.np.array([None, None])
        with raises(ValueError, match='COtype must be "BH", "NS", or "WD"'):
            totest.eddington_limit(BinaryStar, idx=0)
        for CO, r in zip(['WD', 'NS', 'BH'],\
            [(approx(3.68044499210e-17, abs=6e-29),\
              approx(9.08923688740e+7, abs=6e-5)),\
             (approx(2.17747384546e-8, abs=6e-20),\
              approx(0.15362956691, abs=6e-12)),\
             (approx(5.84925262833e-8, abs=6e-20),\
              approx(5.71909584179e-2, abs=6e-14))]):
            BinaryStar.star_1.state_history = totest.np.array([None, CO])
            assert totest.eddington_limit(BinaryStar, idx=0) == r

    def test_beaming(self, BinaryStar):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'binary'"):
            totest.beaming()
        # bad input
        with raises(ValueError, match="Eddington limit is being calculated"+\
                    " for a non-CO"):
            totest.beaming(BinaryStar)
        # examples: 1Msun accretor is star1 with no mass-transfer rate
        BinaryStar.star_1.mass = 1.0
        for CO, r in zip(['WD', 'NS', 'BH'],\
            [(totest.np.nan, 1),\
             (totest.np.nan, 1),\
             (totest.np.nan, 1)]):
            BinaryStar.star_1.state = CO
            assert totest.beaming(BinaryStar) == r
        # examples: 1.2Msun accretor is star1
        BinaryStar.star_1.mass = 1.2
        BinaryStar.lg_mtransfer_rate = -7.5
        for CO, r in zip(['WD', 'NS', 'BH'],\
            [(approx(7.80787388850, abs=6e-12),\
              approx(1.04303350643e-16, abs=6e-28)),\
             (approx(2.98994840792e-8, abs=6e-20), 1),\
             (approx(3.16227766017e-8, abs=6e-20), 1)]):
            BinaryStar.star_1.state = CO
            assert totest.beaming(BinaryStar) == r

    def test_bondi_hoyle(self, BinaryStar, monkeypatch):
        def mock_rand(shape):
            return totest.np.zeros(shape)
        # missing argument
        with raises(TypeError, match="missing 3 required positional"+\
                    " arguments: 'binary', 'accretor', and 'donor'"):
            totest.bondi_hoyle()
        # bad input
        with raises(RuntimeError, match="Failed to converge after 100"+\
                    " iterations"):
            totest.bondi_hoyle(BinaryStar, BinaryStar.star_1,\
                               BinaryStar.star_2)
        # examples:
        BinaryStar.separation = 1.0            #a semi-major axis of 1Rsun
        BinaryStar.eccentricity = 0.1          #a small eccentricity
        BinaryStar.star_1.mass = 1.1           #accretor's mass is 1.1Msun
        BinaryStar.star_1.state = 'WD'         #accretor is WD
        BinaryStar.star_2.mass = 1.2           #donor's mass is 1.2Msun
        BinaryStar.star_2.lg_wind_mdot = -10.0 #donor's wind is 10^{-10}Msun/yr
        BinaryStar.star_2.he_core_radius = 0.2 #donor's he-core rad. is 0.2Rsun
        BinaryStar.star_2.log_R = -0.5         #donor's radius is 10^{-0.5}Rsun
        BinaryStar.star_2.surface_h1 = 0.7     #donor's X_surf=0.7
        BinaryStar.star_2.log_L = 0.3          #donor's lum. is 10^{0.3}Lsun
        with raises(UnboundLocalError, match="cannot access local variable"+\
                    " 'f_m' where it is not associated with a value"):
            # undefined scheme
            totest.bondi_hoyle(BinaryStar, BinaryStar.star_1,\
                               BinaryStar.star_2, scheme='')
        monkeypatch.setattr(totest.np.random, "rand", mock_rand)
        assert totest.bondi_hoyle(BinaryStar, BinaryStar.star_1,\
               BinaryStar.star_2) == approx(3.92668160462e-17, abs=6e-29)
        assert totest.bondi_hoyle(BinaryStar, BinaryStar.star_1,\
               BinaryStar.star_2, scheme='Kudritzki+2000') ==\
               approx(3.92668160462e-17, abs=6e-29)
        BinaryStar.star_2.log_R = 1.5          #donor's radius is 10^{1.5}Rsun
        assert totest.bondi_hoyle(BinaryStar, BinaryStar.star_1,\
               BinaryStar.star_2, scheme='Kudritzki+2000') ==\
               approx(3.92668160462e-17, abs=6e-29)
        BinaryStar.star_2.log_R = -1.5         #donor's radius is 10^{-1.5}Rsun
        assert totest.bondi_hoyle(BinaryStar, BinaryStar.star_1,\
               BinaryStar.star_2, scheme='Kudritzki+2000') == 1e-99
        BinaryStar.star_2.surface_h1 = 0.25    #donor's X_surf=0.25
        assert totest.bondi_hoyle(BinaryStar, BinaryStar.star_1,\
               BinaryStar.star_2) == 1e-99
        BinaryStar.star_2.lg_wind_mdot = -4.0 #donor's wind is 10^{-4}Msun/yr
        assert totest.bondi_hoyle(BinaryStar, BinaryStar.star_1,\
               BinaryStar.star_2) == 1e-99
        assert totest.bondi_hoyle(BinaryStar, BinaryStar.star_1,\
               BinaryStar.star_2, wind_disk_criteria=False) ==\
               approx(5.34028698228e-17, abs=6e-29) # form always a disk

    def test_rejection_sampler(self, monkeypatch):
        def mock_uniform(low=0.0, high=1.0, size=1):
            return totest.np.linspace(low, high, num=size)
        def mock_interp1d(x, y):
            if x[0]>x[-1]:
                raise ValueError
            return interp1d(x, y)
        def mock_pdf(x):
            return 1.0-totest.np.sqrt(x)
        # bad input
        with raises(TypeError, match="'>=' not supported between instances"+\
                    " of 'NoneType' and 'float'"):
            totest.rejection_sampler()
        # examples:
        monkeypatch.setattr(totest.np.random, "uniform", mock_uniform)
        monkeypatch.setattr(totest, "interp1d", mock_interp1d)
        assert totest.np.array_equal(totest.rejection_sampler(\
               x=totest.np.array([0.0, 1.0]), y=totest.np.array([0.4, 0.6]),\
               size=5), totest.np.array([0.0, 0.25, 0.5, 0.75, 1.0]))
        assert totest.np.array_equal(totest.rejection_sampler(\
               x=totest.np.array([1.0, 0.0]), y=totest.np.array([0.2, 0.8]),\
               size=5), totest.np.array([0.0, 0.25, 0.5, 0.0, 0.0]))
        assert totest.np.array_equal(totest.rejection_sampler(\
               x_lim=totest.np.array([0.0, 1.0]), pdf=mock_pdf,\
               size=5), totest.np.array([0.0, 0.25, 0.0, 0.0, 0.0]))

    def test_inverse_sampler(self, monkeypatch):
        def mock_uniform(low=0.0, high=1.0, size=1):
            return totest.np.linspace(low, high, num=size)
        # missing argument
        with raises(TypeError, match="missing 2 required positional"+\
                    " arguments: 'x' and 'y'"):
            totest.inverse_sampler()
        # bad input
        with raises(ValueError, match="diff requires input that is at least"+\
                    " one dimensional"):
            totest.inverse_sampler(x=None, y=None)
        # examples:
        monkeypatch.setattr(totest.np.random, "uniform", mock_uniform)
        assert totest.np.allclose(totest.inverse_sampler(\
               x=totest.np.array([0.0, 1.0]), y=totest.np.array([0.4, 0.6]),\
               size=5), totest.np.array([0.0, 0.29128785, 0.54950976,\
               0.78388218, 1.0]))
        with warns(RuntimeWarning,\
                   match="invalid value encountered in divide"):
            assert totest.np.allclose(totest.inverse_sampler(\
                   x=totest.np.array([0.0, 1.0]),\
                   y=totest.np.array([0.5, 0.5]), size=5),\
                   totest.np.array([0.0, 0.25, 0.5, 0.75, 1.0]))

    def test_histogram_sampler(self, monkeypatch):
        def mock_uniform(low=0.0, high=1.0, size=1):
            return totest.np.linspace(low, high, num=size)
        def mock_choice(a, size=None, replace=True, p=None):
            if isinstance(a, int):
                a=totest.np.arange(a)
            sample = []
            for v,q in zip(a,p):
                sample += round(size*q) * [v]
            if len(sample)<size:
                sample += (size-len(sample)) * [a[-1]]
            return totest.np.array(sample[:size])
        # missing argument
        with raises(TypeError, match="missing 2 required positional"+\
                    " arguments: 'x_edges' and 'y'"):
            totest.histogram_sampler()
        # bad input
        with raises(TypeError, match="'>=' not supported between instances"+\
                    " of 'NoneType' and 'float'"):
            totest.histogram_sampler(x_edges=None, y=None)
        # examples:
        monkeypatch.setattr(totest.np.random, "uniform", mock_uniform)
        monkeypatch.setattr(totest.np.random, "choice", mock_choice)
        assert totest.np.allclose(totest.histogram_sampler(\
               x_edges=totest.np.array([0.0, 0.5, 1.0]),\
               y=totest.np.array([0.2, 0.8]), size=5),\
               totest.np.array([0.0, 0.5, 0.66666667, 0.83333333, 1.0]))

    def test_read_histogram_from_file(self, csv_path, csv_file, csv_path2,\
                                      csv_file2):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'path'"):
            totest.read_histogram_from_file()
        # bad input
        with raises(TypeError, match="expected str, bytes or os.PathLike"+\
                    " object, not NoneType"):
            totest.read_histogram_from_file(path=None)
        # examples:
        with raises(RuntimeError, match="More than two lines found in the"+\
                    " histogram document."):
            totest.read_histogram_from_file(path=csv_path)
        assert totest.np.array_equal(totest.read_histogram_from_file(\
               path=csv_path2)[0],totest.np.array([0.0, 1.0, 2.0]))
        assert totest.np.array_equal(totest.read_histogram_from_file(\
               path=csv_path2)[1],totest.np.array([1.0, 1.0, 1.0]))

    def test_inspiral_timescale_from_separation(self):
        # missing argument
        with raises(TypeError, match="missing 4 required positional"+\
                    " arguments: 'star1_mass', 'star2_mass', 'separation',"+\
                    " and 'eccentricity'"):
            totest.inspiral_timescale_from_separation()
        # bad input
        with raises(TypeError) as error_info:
            totest.inspiral_timescale_from_separation(None, None, None, None)
        assert error_info.value.args[0] == "unsupported operand type(s) for"+\
               " *: 'NoneType' and 'float'"
        with raises(ValueError, match="Mass of star 1 is <= 0, which is not"+\
                    " a physical value."):
            totest.inspiral_timescale_from_separation(0.0, 0.5, 0.5, 0.5)
        with raises(ValueError, match="Mass of star 2 is <= 0, which is not"+\
                    " a physical value."):
            totest.inspiral_timescale_from_separation(0.5, 0.0, 0.5, 0.5)
        with raises(ValueError, match="Separation is <= 0, which is not a"+\
                    " physical value."):
            totest.inspiral_timescale_from_separation(0.5, 0.5, 0.0, 0.5)
        with raises(ValueError, match="Eccentricity is < 0, which is not a"+\
                    " physical value."):
            totest.inspiral_timescale_from_separation(0.5, 0.5, 0.5, -0.5)
        with raises(ValueError, match="Eccentricity is >= 1, which is not a"+\
                    " physical value."):
            totest.inspiral_timescale_from_separation(0.5, 0.5, 0.5, 1.5)
        # examples:
        assert totest.inspiral_timescale_from_separation(0.5, 0.5, 0.5, 0.5)\
               == approx(13.44121924697, abs=6e-12)
        assert totest.inspiral_timescale_from_separation(1.5, 0.5, 0.5, 0.5)\
               == approx(2.24020320783, abs=6e-12)
        assert totest.inspiral_timescale_from_separation(0.5, 1.5, 0.5, 0.5)\
               == approx(2.24020320783, abs=6e-12)
        assert totest.inspiral_timescale_from_separation(0.5, 0.5, 1.5, 0.5)\
               == approx(1.08873875900e+3, abs=6e-9)
        assert totest.inspiral_timescale_from_separation(0.5, 0.5, 0.5, 0.0)\
               == approx(37.56647511040, abs=6e-12)
        assert totest.inspiral_timescale_from_separation(0.5, 0.5, 0.5,\
               1.0e-5) == approx(37.56647487803, abs=6e-12)
        assert totest.inspiral_timescale_from_separation(0.5, 0.5, 0.5,\
               1.0-1.0e-5) == approx(2.41297178064e-15, abs=6e-27)

    def test_inspiral_timescale_from_orbital_period(self):
        # missing argument
        with raises(TypeError, match="missing 4 required positional"+\
                    " arguments: 'star1_mass', 'star2_mass',"+\
                    " 'orbital_period', and 'eccentricity'"):
            totest.inspiral_timescale_from_orbital_period()
        # bad input
        with raises(TypeError) as error_info:
            totest.inspiral_timescale_from_orbital_period(None, None, None, None)
        assert error_info.value.args[0] == "unsupported operand type(s) for"+\
               " *: 'NoneType' and 'float'"
        with raises(ValueError, match="Mass of star 1 is <= 0, which is not"+\
                    " a physical value."):
            totest.inspiral_timescale_from_orbital_period(0.0, 0.5, 0.5, 0.5)
        with raises(ValueError, match="Mass of star 2 is <= 0, which is not"+\
                    " a physical value."):
            totest.inspiral_timescale_from_orbital_period(0.5, 0.0, 0.5, 0.5)
        with raises(ValueError, match="Separation is <= 0, which is not a"+\
                    " physical value."):
            totest.inspiral_timescale_from_orbital_period(0.5, 0.5, 0.0, 0.5)
        with raises(ValueError, match="Eccentricity is < 0, which is not a"+\
                    " physical value."):
            totest.inspiral_timescale_from_orbital_period(0.5, 0.5, 0.5, -0.5)
        with raises(ValueError, match="Eccentricity is >= 1, which is not a"+\
                    " physical value."):
            totest.inspiral_timescale_from_orbital_period(0.5, 0.5, 0.5, 1.5)
        # examples:
        assert totest.inspiral_timescale_from_orbital_period(0.5, 0.5, 0.5,\
               0.5) == approx(1.06110684309e+4, abs=6e-8)
        assert totest.inspiral_timescale_from_orbital_period(1.5, 0.5, 0.5,\
               0.5) == approx(4.45636949266e+3, abs=6e-9)
        assert totest.inspiral_timescale_from_orbital_period(0.5, 1.5, 0.5,\
               0.5) == approx(4.45636949266e+3, abs=6e-9)
        assert totest.inspiral_timescale_from_orbital_period(0.5, 0.5, 1.5,\
               0.5) == approx(1.98647206096e+5, abs=6e-7)
        assert totest.inspiral_timescale_from_orbital_period(0.5, 0.5, 0.5,\
               0.0) == approx(2.96565684095e+4, abs=6e-8)
        assert totest.inspiral_timescale_from_orbital_period(0.5, 0.5, 0.5,\
               1.0e-5) == approx(2.96565682261e+4, abs=6e-8)
        assert totest.inspiral_timescale_from_orbital_period(0.5, 0.5, 0.5,\
               1.0-1.0e-5) == approx(1.90490224255e-12, abs=6e-24)

    def test_spin_stable_mass_transfer(self):
        # missing argument
        with raises(TypeError, match="missing 3 required positional"+\
                    " arguments: 'spin_i', 'star_mass_preMT', and"+\
                    " 'star_mass_postMT'"):
            totest.spin_stable_mass_transfer()
        # bad input
        with raises(TypeError) as error_info:
            totest.spin_stable_mass_transfer(None, 1.0, 1.0)
        assert error_info.value.args[0] == "unsupported operand type(s) for"+\
               " ** or pow(): 'NoneType' and 'int'"
        # examples:
        assert totest.spin_stable_mass_transfer(1.0, None, 1.0) is None
        assert totest.spin_stable_mass_transfer(1.0, 1.0, None) is None
        assert totest.spin_stable_mass_transfer(0.5, 0.9, 1.0) ==\
               approx(0.69217452285, abs=6e-12)
        assert totest.np.isnan(totest.spin_stable_mass_transfer(1.0, 1.0, 0.1))
        assert totest.spin_stable_mass_transfer(1.0, 0.1, 1.0) == 1.0

    def test_check_state_of_star(self, SingleStar):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'star'"):
            totest.check_state_of_star()
        # bad input
        with raises(AttributeError, match="'NoneType' object has no"+\
                    " attribute 'mass'"):
            totest.check_state_of_star(None, None, None)
        # examples:
        assert totest.check_state_of_star(SingleStar) ==\
               "undetermined_evolutionary_state"
        assert totest.check_state_of_star(SingleStar, i=0) ==\
               "undetermined_evolutionary_state"
        for CO,m in zip(['WD', 'NS', 'BH'], [1.0, 1.5, 9.0]):
            SingleStar.mass = m
            SingleStar.state = CO
            assert totest.check_state_of_star(SingleStar, star_CO=True) ==\
                   SingleStar.state

    def test_check_state_of_star_history_array(self, SingleStar):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'star'"):
            totest.check_state_of_star_history_array()
        # bad input
        with raises(TypeError, match="bad operand type for unary -:"+\
                    " 'NoneType'"):
            totest.check_state_of_star_history_array(None, None, None)
        with raises(IndexError, match="list index out of range"):
            totest.check_state_of_star_history_array(SingleStar, N=10)
        # examples:
        assert totest.check_state_of_star_history_array(SingleStar) ==\
               ["undetermined_evolutionary_state"]
        n = 3
        SingleStar.mass_history = n*SingleStar.mass_history
        SingleStar.surface_h1_history = n*SingleStar.surface_h1_history
        SingleStar.center_h1_history = n*SingleStar.center_h1_history
        SingleStar.center_he4_history = n*SingleStar.center_he4_history
        SingleStar.center_c12_history = n*SingleStar.center_c12_history
        SingleStar.log_LH_history = n*SingleStar.log_LH_history
        SingleStar.log_LHe_history = n*SingleStar.log_LHe_history
        SingleStar.log_Lnuc_history = n*SingleStar.log_Lnuc_history
        for i in range(n):
            assert totest.check_state_of_star_history_array(SingleStar, N=i+1)\
                   == (i+1)*["undetermined_evolutionary_state"]
        for CO,m in zip(['WD', 'NS', 'BH'], [1.0, 1.5, 9.0]):
            SingleStar.mass_history[-1] = m
            SingleStar.state = CO
            assert totest.check_state_of_star_history_array(SingleStar,\
                   star_CO=True) == [SingleStar.state]

    def test_get_binary_state_and_event_and_mt_case(self, BinaryStar,\
                                                    monkeypatch):
        def mock_infer_mass_transfer_case(rl_relative_overflow,\
                                          lg_mtransfer_rate, donor_state,\
                                          verbose=False):
            if rl_relative_overflow is not None:
                if rl_relative_overflow>0:
                    return totest.MT_CASE_A
                else:
                    return totest.MT_CASE_UNDETERMINED
            else:
                return totest.MT_CASE_NO_RLO
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'binary'"):
            totest.get_binary_state_and_event_and_mt_case()
        # bad input
        with raises(TypeError, match="argument of type 'NoneType' is not"+\
                    " iterable"):
            totest.get_binary_state_and_event_and_mt_case(BinaryStar)
        # examples: no binary
        assert totest.get_binary_state_and_event_and_mt_case(None) ==\
               [None, None, 'None']
        # examples: not converged
        assert totest.get_binary_state_and_event_and_mt_case(BinaryStar,\
               interpolation_class='not_converged') == [None, None, 'None']
        # examples: initial MT
        assert totest.get_binary_state_and_event_and_mt_case(BinaryStar,\
               interpolation_class='initial_MT') ==\
               ['initial_RLOF', None, 'None']
        # examples: detached
        BinaryStar.star_1.state = "test_state"
        BinaryStar.star_2.state = "test_state"
        assert totest.get_binary_state_and_event_and_mt_case(BinaryStar) ==\
               ['detached', None, 'None']
        BinaryStar.star_1.state_history[0] = "test_state"
        BinaryStar.star_2.state_history[0] = "test_state"
        assert totest.get_binary_state_and_event_and_mt_case(BinaryStar, i=0)\
               == ['detached', None, 'None']
        # examples: mass transfer
        monkeypatch.setattr(totest, "infer_mass_transfer_case",\
                            mock_infer_mass_transfer_case)
        ## donor is star 1
        BinaryStar.rl_relative_overflow_1 = 1.0
        for IC,e in zip([None, 'unstable_MT'], [None, 'oCE1']):
            assert totest.get_binary_state_and_event_and_mt_case(BinaryStar,\
                   interpolation_class=IC) == ['RLO1', e, 'case_A1']
        ## both stars overfill RL leading to double CE initiated by star 1
        BinaryStar.rl_relative_overflow_2 = 1.0
        for IC,e in zip([None, 'unstable_MT'], [None, 'oDoubleCE1']):
            assert totest.get_binary_state_and_event_and_mt_case(BinaryStar,\
                   interpolation_class=IC) == ['contact', e, 'None']
        ## classical CE initiated by star 1
        BinaryStar.star_2.state = "BH"
        assert totest.get_binary_state_and_event_and_mt_case(BinaryStar,\
               interpolation_class='unstable_MT') == ['contact', 'oCE1', 'None']
        BinaryStar.star_2.state = "test_state"
        ## both stars overfill RL leading to double CE initiated by star 2
        BinaryStar.rl_relative_overflow_2 = 2.0
        for IC,e in zip([None, 'unstable_MT'], [None, 'oDoubleCE2']):
            assert totest.get_binary_state_and_event_and_mt_case(BinaryStar,\
                   interpolation_class=IC) == ['contact', e, 'None']
        ## classical CE initiated by star 2
        BinaryStar.star_1.state = "BH"
        assert totest.get_binary_state_and_event_and_mt_case(BinaryStar,\
               interpolation_class='unstable_MT') == ['contact', 'oCE2', 'None']
        BinaryStar.star_1.state = "test_state"
        ## donor is star 2
        BinaryStar.rl_relative_overflow_1 = None
        for IC,e in zip([None, 'unstable_MT'], [None, 'oCE2']):
            assert totest.get_binary_state_and_event_and_mt_case(BinaryStar,\
                   interpolation_class=IC) == ['RLO2', e, 'case_A2']
        # examples: undefined
        BinaryStar.rl_relative_overflow_2 = -1.0
        assert totest.get_binary_state_and_event_and_mt_case(BinaryStar) ==\
               ['undefined', None, 'None']
        # examples: star 2 becomes WD
        BinaryStar.star_2.center_gamma = 10.0
        assert totest.get_binary_state_and_event_and_mt_case(BinaryStar) ==\
               ['undefined', 'CC2', 'None']
        # examples: star 1 becomes WD
        BinaryStar.star_1.center_gamma = 10.0
        assert totest.get_binary_state_and_event_and_mt_case(BinaryStar) ==\
               ['undefined', 'CC1', 'None']
        # examples: no center_gamma_history
        BinaryStar.star_1.center_gamma_history=[]
        BinaryStar.star_2.center_gamma_history=[]
        assert totest.get_binary_state_and_event_and_mt_case(BinaryStar, i=0)\
               == ['detached', None, 'None']

    def test_get_binary_state_and_event_and_mt_case_array(self, BinaryStar,\
                                                          monkeypatch):
        def mock_get_binary_state_and_event_and_mt_case(binary,\
            interpolation_class=None, i=None, verbose=False):
            return ['Unit', 'Test', 'None']
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'binary'"):
            totest.get_binary_state_and_event_and_mt_case_array()
        # bad input
        with raises(TypeError, match="argument of type 'NoneType' is not"+\
                    " iterable"):
            totest.get_binary_state_and_event_and_mt_case_array(BinaryStar)
        with raises(IndexError, match="list index out of range"):
            totest.get_binary_state_and_event_and_mt_case_array(BinaryStar,\
            N=10)
        # examples: no binary
        assert totest.get_binary_state_and_event_and_mt_case_array(None) ==\
               (None, None, 'None')
        # examples: detached
        BinaryStar.star_1.state = "test_state"
        BinaryStar.star_2.state = "test_state"
        assert totest.get_binary_state_and_event_and_mt_case_array(BinaryStar)\
               == ('detached', None, 'None')
        # examples: several mocked entries
        monkeypatch.setattr(totest, "get_binary_state_and_event_and_mt_case",\
                            mock_get_binary_state_and_event_and_mt_case)
        mock_return = mock_get_binary_state_and_event_and_mt_case(BinaryStar)
        for i in range(3):
            assert totest.get_binary_state_and_event_and_mt_case_array(\
                   BinaryStar, N=i+1) == ((i+1)*[mock_return[0]],\
                   (i+1)*[mock_return[1]], (i+1)*[mock_return[2]])

    def test_CO_radius(self):
        # missing argument
        with raises(TypeError, match="missing 2 required positional"+\
                    " arguments: 'M' and 'COtype'"):
            totest.CO_radius()
        # bad input
        with raises(ValueError, match="Compact object mass must be a"+\
                    " positive value"):
            totest.CO_radius(0.0, 'BH')
        with raises(ValueError, match="COtype not in the list of valid"+\
                    " options"):
            totest.CO_radius(1.0, 'TEST')
        # examples:
        for CO,m,r in zip(['WD', 'WD', 'NS', 'NS', 'BH', 'BH'],\
                          [0.5, 1.0, 1.5, 2.0, 7.0, 9.0],\
                          [approx(5.24982189818e-3, abs=6e-15),\
                           approx(4.16678640191e-3, abs=6e-15),\
                           approx(1.79602862151e-5, abs=6e-17),\
                           approx(1.79602862151e-5, abs=6e-17),\
                           totest.Schwarzschild_Radius(7.0),\
                           totest.Schwarzschild_Radius(9.0)]):
            assert totest.CO_radius(m, CO) == r

    def test_He_MS_lifetime(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'mass'"):
            totest.He_MS_lifetime()
        # bad input
        with raises(TypeError, match="'<' not supported between instances"+\
                    " of 'NoneType' and 'float'"):
            totest.He_MS_lifetime(None)
        # examples:
        for m,t in zip([1.0, 3.0, 30.0, 300.0], [1.0e+8,\
                        approx(3.47136114375e+7, abs=6e-5),\
                        approx(6.95883527992e+5, abs=6e-7), 3.0e+5]):
            assert totest.He_MS_lifetime(m) == t

    def test_Schwarzschild_Radius(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'M'"):
            totest.Schwarzschild_Radius()
        # bad input
        with raises(TypeError) as error_info:
            totest.Schwarzschild_Radius(None)
        assert error_info.value.args[0] == "unsupported operand type(s) for"+\
               " *: 'float' and 'NoneType'"
        # examples:
        for m,r in zip([7.0, 70.0], [\
                        approx(2.97147953078e-5, abs=6e-17),\
                        approx(2.97147953078e-4, abs=6e-16)]):
            assert totest.Schwarzschild_Radius(m) == r

    def test_flip_stars(self, BinaryStar):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'binary'"):
            totest.flip_stars()
        # bad input
        with raises(AttributeError, match="'NoneType' object has no"+\
                    " attribute 'star_1'"):
            totest.flip_stars(None)
        # examples: set stars and check swap
        for attr in BinaryStar.star_1.__dict__:
            setattr(BinaryStar.star_1, attr, 1)
        for attr in BinaryStar.star_2.__dict__:
            setattr(BinaryStar.star_2, attr, 2)
        totest.flip_stars(BinaryStar)
        for attr in BinaryStar.star_1.__dict__:
            assert getattr(BinaryStar.star_1, attr) == 2
        for attr in BinaryStar.star_2.__dict__:
            assert getattr(BinaryStar.star_2, attr) == 1
        # examples: binary states
        for s,n in zip(['RLO1', 'RLO2'], ['RLO2', 'RLO1']):
            BinaryStar.state = s
            BinaryStar.state_history[0] = s
            totest.flip_stars(BinaryStar)
            assert BinaryStar.state == n
            assert BinaryStar.state_history[0] == n
        # examples: binary events
        for e,n in zip(['oRLO1', 'oRLO2', 'oCE1', 'oCE2', 'CC1', 'CC2'],\
                       ['oRLO2', 'oRLO1', 'oCE2', 'oCE1', 'CC2', 'CC1']):
            BinaryStar.event = e
            BinaryStar.event_history[0] = e
            totest.flip_stars(BinaryStar)
            assert BinaryStar.event == n
            assert BinaryStar.event_history[0] == n
        # examples: set binary values and check swap
        for v in BinaryStar.__dict__:
            if v[-1]==1:
                setattr(BinaryStar, v, 1)
            elif v[-1]==2:
                setattr(BinaryStar, v, 2)
        totest.flip_stars(BinaryStar)
        for v in BinaryStar.__dict__:
            if v[-1]==1:
                assert getattr(BinaryStar, v) == 2
            elif v[-1]==2:
                assert getattr(BinaryStar, v) == 1

    def test_set_binary_to_failed(self, BinaryStar):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'binary'"):
            totest.set_binary_to_failed()
        # bad input
        with raises(AttributeError, match="'NoneType' object has no"+\
                    " attribute 'state'"):
            totest.set_binary_to_failed(None)
        # examples:
        totest.set_binary_to_failed(BinaryStar)
        assert BinaryStar.state == "ERR"
        assert BinaryStar.event == "FAILED"

    def test_infer_star_state(self):
        # bad input
        with raises(TypeError, match="'<=' not supported between instances"+\
                    " of 'NoneType' and 'float'"):
            totest.infer_star_state(star_CO=True)
        # examples: undetermined
        assert totest.infer_star_state() == totest.STATE_UNDETERMINED
        # examples: compact objects
        for m,CO in zip([0.5*totest.STATE_NS_STARMASS_UPPER_LIMIT,\
                             totest.STATE_NS_STARMASS_UPPER_LIMIT,\
                         2.0*totest.STATE_NS_STARMASS_UPPER_LIMIT],\
                        ["NS", "NS", "BH"]):
            assert totest.infer_star_state(star_mass=m, star_CO=True) == CO
        # examples: loop over all cases
        THNA = totest.THRESHOLD_HE_NAKED_ABUNDANCE
        TCA = totest.THRESHOLD_CENTRAL_ABUNDANCE
        LBT = totest.LOG10_BURNING_THRESHOLD
        for sH1 in [2.0*THNA, THNA, 0.5*THNA]:
            for cH1 in [2.0*TCA, TCA, 0.5*TCA]:
                for cHe4 in [2.0*TCA, TCA, 0.5*TCA]:
                    for cC12 in [2.0*TCA, TCA, 0.5*TCA]:
                        for lgLH in [2.0*LBT, LBT, 0.5*LBT]:
                            for lgLHe in [2.0*LBT, LBT, 0.5*LBT]:
                                if sH1>THNA:
                                    rich = "H-rich"
                                elif sH1<cH1:
                                    rich = "accreted_He"
                                else:
                                    rich = "stripped_He"
                                if cH1<=TCA and cHe4<=TCA and cC12<=TCA:
                                    burn = "Central_C_depletion"
                                elif cH1<=TCA and cHe4<=TCA and cC12>TCA:
                                    burn = "Central_He_depleted"
                                elif cH1>TCA and lgLH>LBT:
                                    burn = "Core_H_burning"
                                elif cH1>TCA and lgLH<=LBT:
                                    burn = "non_burning"
                                elif lgLHe>LBT:
                                    burn = "Core_He_burning"
                                elif lgLH>LBT:
                                    burn = "Shell_H_burning"
                                else:
                                    burn = "non_burning"
                                assert totest.infer_star_state(surface_h1=sH1,\
                                       center_h1=cH1, center_he4=cHe4,\
                                       center_c12=cC12, log_LH=lgLH,\
                                       log_LHe=lgLHe, log_Lnuc=lgLH+lgLHe) ==\
                                       rich+"_"+burn

    def test_infer_mass_transfer_case(self, capsys):
        # missing argument
        with raises(TypeError, match="missing 3 required positional"+\
                    " arguments: 'rl_relative_overflow',"+\
                    " 'lg_mtransfer_rate', and 'donor_state'"):
            totest.infer_mass_transfer_case()
        # examples: no RLO
        assert totest.infer_mass_transfer_case(None, None, None) ==\
               totest.MT_CASE_NO_RLO
        assert totest.infer_mass_transfer_case(None, 1.0, "test") ==\
               totest.MT_CASE_NO_RLO
        assert totest.infer_mass_transfer_case(1.0, None, "test") ==\
               totest.MT_CASE_NO_RLO
        RROT = totest.RL_RELATIVE_OVERFLOW_THRESHOLD
        LMRT = totest.LG_MTRANSFER_RATE_THRESHOLD
        assert totest.infer_mass_transfer_case(min(-1.0, RROT), LMRT, "test")\
               == totest.MT_CASE_NO_RLO
        assert totest.infer_mass_transfer_case(min(-1.0, RROT), LMRT, "test",\
               verbose=True) == totest.MT_CASE_NO_RLO
        assert "checking rl_relative_overflow / lg_mtransfer_rate" in\
               capsys.readouterr().out
        # examples: MT cases
        for ds,c in zip(["test_non_burning", "H-rich_Core_H_burning",\
                         "H-rich_Core_He_burning", "H-rich_Shell_H_burning",\
                         "H-rich_Central_He_depleted",\
                         "H-rich_Central_C_depletion", "H-rich_undetermined",\
                         "stripped_He_Core_He_burning",\
                         "stripped_He_Central_He_depleted",\
                         "stripped_He_Central_C_depletion",\
                         "stripped_He_undetermined", "test_undetermined"],\
                        [totest.MT_CASE_NONBURNING, totest.MT_CASE_A,\
                         totest.MT_CASE_B, totest.MT_CASE_B, totest.MT_CASE_C,\
                         totest.MT_CASE_C, totest.MT_CASE_UNDETERMINED,\
                         totest.MT_CASE_BA, totest.MT_CASE_BB,\
                         totest.MT_CASE_BB, totest.MT_CASE_UNDETERMINED,\
                         totest.MT_CASE_UNDETERMINED]):
            assert totest.infer_mass_transfer_case(2*RROT, 2*LMRT, ds) == c

    def test_cumulative_mass_transfer_numeric(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'MT_cases'"):
            totest.cumulative_mass_transfer_numeric()
        # bad input
        with raises(TypeError, match="object of type 'NoneType' has no len()"):
            totest.cumulative_mass_transfer_numeric(None)
        # examples: no cases
        assert totest.cumulative_mass_transfer_numeric([]) ==\
               [totest.MT_CASE_UNDETERMINED]
        # examples: undetermined
        assert totest.cumulative_mass_transfer_numeric(\
               [totest.MT_CASE_UNDETERMINED]) == [totest.MT_CASE_UNDETERMINED]
        # examples: no RLO
        assert totest.cumulative_mass_transfer_numeric(\
               [totest.MT_CASE_NO_RLO]) == [totest.MT_CASE_NO_RLO]
        # examples: cut out undetermined
        assert totest.cumulative_mass_transfer_numeric([totest.MT_CASE_A,\
               totest.MT_CASE_UNDETERMINED, totest.MT_CASE_A,\
               totest.MT_CASE_B, totest.MT_CASE_A, totest.MT_CASE_A]) ==\
               [totest.MT_CASE_UNDETERMINED, totest.MT_CASE_A,\
                totest.MT_CASE_B, totest.MT_CASE_A]
        # examples: cut out no RLO from strings
        assert totest.cumulative_mass_transfer_numeric(["A", "no_RLO", "A",\
               "B", "A", "A"]) == [totest.MT_CASE_A, totest.MT_CASE_B,\
               totest.MT_CASE_A]

    def test_cumulative_mass_transfer_string(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'cumulative_integers'"):
            totest.cumulative_mass_transfer_string()
        # bad input
        with raises(TypeError, match="object of type 'NoneType' has no len()"):
            totest.cumulative_mass_transfer_string(None)
        with raises(AssertionError): # missing case information
            totest.cumulative_mass_transfer_string([])
        with warns(InappropriateValueWarning, match="Unknown MT case:"):
            # unknown case
            totest.cumulative_mass_transfer_string([-1])
        # examples: undetermined
        assert totest.cumulative_mass_transfer_string(\
               [totest.MT_CASE_UNDETERMINED]) == "?"
        # examples: no RLO
        assert totest.cumulative_mass_transfer_string(\
               [totest.MT_CASE_NO_RLO]) == "no_RLO"
        # examples: cases for both stars
        for c in totest.ALL_RLO_CASES:
            assert totest.cumulative_mass_transfer_string([c, 10+c]) ==\
                   "case_"+totest.MT_CASE_TO_STR[c]+"1/"+\
                   totest.MT_CASE_TO_STR[c]+"2"

    def test_cumulative_mass_transfer_flag(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'MT_cases'"):
            totest.cumulative_mass_transfer_flag()
        # bad input
        with raises(AttributeError, match="'NoneType' object has no"+\
                    " attribute 'copy'"):
            totest.cumulative_mass_transfer_flag(None)
        with warns(EvolutionWarning, match="MT case with unknown donor:"):
            with warns(InappropriateValueWarning):
                totest.cumulative_mass_transfer_flag([30], shift_cases=True)
        # examples:
        assert totest.cumulative_mass_transfer_flag([totest.MT_CASE_A,\
               totest.MT_CASE_B, 10+totest.MT_CASE_A, totest.MT_CASE_A,\
               10+totest.MT_CASE_B, 10+totest.MT_CASE_A]) ==\
               'case_A1/B1/A2/A1/B2/A2'
        assert totest.cumulative_mass_transfer_flag([totest.MT_CASE_A,\
               totest.MT_CASE_B, 10+totest.MT_CASE_A, totest.MT_CASE_A,\
               10+totest.MT_CASE_B, 10+totest.MT_CASE_A],\
               shift_cases=True) == 'case_A1/B1/A2/B1/B2'

    def test_get_i_He_depl(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'history'"):
            totest.get_i_He_depl()
        # bad input
        with raises(AttributeError, match="'NoneType' object has no"+\
                    " attribute 'dtype'"):
            totest.get_i_He_depl(None)
        with raises(TypeError, match="argument of type 'NoneType' is not"+\
                    " iterable"):
            totest.get_i_He_depl(totest.np.array([]))
        # examples: incomplete history to search through
        assert totest.get_i_He_depl(totest.np.array([[]],\
               dtype=([('surface_h1', 'f8')]))) == -1
        # examples: no he depletion
        assert totest.get_i_He_depl(totest.np.array([(1.0, 0.0, 1.0, 0.0, 0.0,\
               1.0, 1.0), (1.0, 0.0, 0.5, 0.5, 0.0, 1.0, 1.0), (1.0, 0.0, 0.2,\
               0.2, 0.0, 1.0, 1.0)], dtype=([('surface_h1', 'f8'),\
               ('center_h1', 'f8'), ('center_he4', 'f8'),\
               ('center_c12', 'f8'), ('log_LH', 'f8'), ('log_LHe', 'f8'),\
               ('log_Lnuc', 'f8')]))) == -1
        # examples: he depletion at 1
        assert totest.get_i_He_depl(totest.np.array([(1.0, 0.0, 1.0, 0.0, 0.0,\
               1.0, 1.0), (1.0, 0.0, 0.0, 0.5, 0.0, 1.0, 1.0), (1.0, 0.0, 0.0,\
               0.0, 0.0, 1.0, 1.0)], dtype=([('surface_h1', 'f8'),\
               ('center_h1', 'f8'), ('center_he4', 'f8'),\
               ('center_c12', 'f8'), ('log_LH', 'f8'), ('log_LHe', 'f8'),\
               ('log_Lnuc', 'f8')]))) == 1

    def test_calculate_Patton20_values_at_He_depl(self, SingleStar):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'star'"):
            totest.calculate_Patton20_values_at_He_depl()
        # bad input
        with raises(AttributeError, match="'NoneType' object has no"+\
                    " attribute 'state_history'"):
            totest.calculate_Patton20_values_at_He_depl(None)
        # examples: no history
        SingleStar.co_core_mass_at_He_depletion = 0.5
        SingleStar.avg_c_in_c_core_at_He_depletion = 0.5
        SingleStar.state_history = None
        totest.calculate_Patton20_values_at_He_depl(SingleStar)
        assert SingleStar.co_core_mass_at_He_depletion is None
        assert SingleStar.avg_c_in_c_core_at_He_depletion is None
        # examples: no He_depleted in history
        SingleStar.co_core_mass_at_He_depletion = 0.5
        SingleStar.avg_c_in_c_core_at_He_depletion = 0.5
        SingleStar.state_history = ["test"]
        totest.calculate_Patton20_values_at_He_depl(SingleStar)
        assert SingleStar.co_core_mass_at_He_depletion is None
        assert SingleStar.avg_c_in_c_core_at_He_depletion is None
        # examples: loop through star types with He depletion
        for s,v in zip(["H-rich_Central_He_depleted",\
                        "stripped_He_Central_He_depleted"], [0.1, 0.2]):
            SingleStar.state_history = ["test", s, s]
            SingleStar.co_core_mass_history = [0.0, v, 1.0]
            SingleStar.avg_c_in_c_core_history = [0.0, v, 1.0]
            totest.calculate_Patton20_values_at_He_depl(SingleStar)
            assert SingleStar.co_core_mass_at_He_depletion == v
            assert SingleStar.avg_c_in_c_core_at_He_depletion == v

    def test_CEE_parameters_from_core_abundance_thresholds(self, monkeypatch,\
                                                           capsys, SingleStar):
        def mock_calculate_core_boundary(donor_mass, donor_star_state,\
            profile, mc1_i=None, core_element_fraction_definition=None,\
            CO_core_in_Hrich_star=False):
            return core_element_fraction_definition
        def mock_calculate_lambda_from_profile(profile, donor_star_state,\
            m1_i=totest.np.nan, radius1=totest.np.nan,\
            common_envelope_option_for_lambda=\
            totest.DEFAULT_CE_OPTION_FOR_LAMBDA,\
            core_element_fraction_definition=0.1, ind_core=None,\
            common_envelope_alpha_thermal=1.0, tolerance=0.001,\
            CO_core_in_Hrich_star=False, verbose=False):
            return core_element_fraction_definition,\
                   core_element_fraction_definition,\
                   core_element_fraction_definition
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'star'"):
            totest.CEE_parameters_from_core_abundance_thresholds()
        # bad input
        with raises(AttributeError, match="'NoneType' object has no"+\
                    " attribute 'mass'"):
            totest.CEE_parameters_from_core_abundance_thresholds(None)
        with raises(TypeError) as error_info:
            totest.CEE_parameters_from_core_abundance_thresholds(SingleStar)
        assert error_info.value.args[0] == "unsupported operand type(s) for"+\
               " ** or pow(): 'float' and 'NoneType'"
        SingleStar.log_R = 0.0
        SingleStar.profile = totest.np.array([(1.0), (1.0), (1.0)],\
                                             dtype=([('mass', 'f8')]))
        with raises(TypeError, match="argument of type 'NoneType' is not"+\
                    " iterable"):
            totest.CEE_parameters_from_core_abundance_thresholds(SingleStar)
        # examples: missing state with profile and verbose
        SingleStar.state = "test_state"
        totest.CEE_parameters_from_core_abundance_thresholds(SingleStar,\
                                                             verbose=True)
        captured_out = capsys.readouterr().out
        for out in ["is not what expected during the "+\
                    "CEE_parameters_from_core_abundance_thresholds.",\
                    "star_state", "m_core_CE_1cent,m_core_CE_10cent,"+\
                    "m_core_CE_30cent,m_core_CE_pure_He_star_10cent",\
                    "r_core_CE_1cent,r_core_CE_10cent,r_core_CE_30cent,"+\
                    "r_core_CE_pure_He_star_10cent", "lambda_CE_1cent,"+\
                    "lambda_CE_10cent,lambda_CE_30cent,"+\
                    "lambda_CE_pure_He_star_10cent"]:
            assert out in captured_out
        for attr in ['m_core_CE_1cent', 'm_core_CE_10cent',\
                     'm_core_CE_30cent', 'm_core_CE_pure_He_star_10cent',\
                     'r_core_CE_1cent', 'r_core_CE_10cent',\
                     'r_core_CE_30cent', 'r_core_CE_pure_He_star_10cent',\
                     'lambda_CE_1cent', 'lambda_CE_10cent',\
                     'lambda_CE_30cent', 'lambda_CE_pure_He_star_10cent']:
            assert totest.np.isnan(getattr(SingleStar, attr))
        monkeypatch.setattr(totest, "calculate_core_boundary",\
                            mock_calculate_core_boundary)
        monkeypatch.setattr(totest, "calculate_lambda_from_profile",\
                            mock_calculate_lambda_from_profile)
        for s in ["test_state", "H-rich_test", "stripped_He_test"]:
            SingleStar.state = s
            totest.CEE_parameters_from_core_abundance_thresholds(SingleStar)
            assert capsys.readouterr().out == ""
            for attr in ['m_core_CE_1cent', 'm_core_CE_10cent',\
                         'm_core_CE_30cent', 'm_core_CE_pure_He_star_10cent',\
                         'r_core_CE_1cent', 'r_core_CE_10cent',\
                         'r_core_CE_30cent', 'r_core_CE_pure_He_star_10cent',\
                         'lambda_CE_1cent', 'lambda_CE_10cent',\
                         'lambda_CE_30cent', 'lambda_CE_pure_He_star_10cent']:
                if s == "test_state":
                    assert totest.np.isnan(getattr(SingleStar, attr))
                elif (("stripped_He" in s) and ("pure_He" not in attr) and\
                      ("m_core" in attr)):
                    assert getattr(SingleStar, attr) == SingleStar.mass
                elif (("stripped_He" in s) and ("pure_He" not in attr) and\
                      ("r_core" in attr)):
                    assert getattr(SingleStar, attr) == 10**SingleStar.log_R
                elif (("stripped_He" in s) and ("pure_He" not in attr) and\
                      ("lambda" in attr)):
                    assert totest.np.isnan(getattr(SingleStar, attr))
                elif "1cent" in attr:
                    assert getattr(SingleStar, attr) == 0.01
                elif "10cent" in attr:
                    assert getattr(SingleStar, attr) == 0.1
                elif "30cent" in attr:
                    assert getattr(SingleStar, attr) == 0.3
                else:
                    assert totest.np.isnan(getattr(SingleStar, attr))
        # examples: no profile
        SingleStar.profile = None
        totest.CEE_parameters_from_core_abundance_thresholds(SingleStar)
        for attr in ['m_core_CE_1cent', 'm_core_CE_10cent',\
                     'm_core_CE_30cent', 'm_core_CE_pure_He_star_10cent',\
                     'r_core_CE_1cent', 'r_core_CE_10cent',\
                     'r_core_CE_30cent', 'r_core_CE_pure_He_star_10cent',\
                     'lambda_CE_1cent', 'lambda_CE_10cent',\
                     'lambda_CE_30cent', 'lambda_CE_pure_He_star_10cent']:
            assert totest.np.isnan(getattr(SingleStar, attr))

    def test_initialize_empty_array(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'arr'"):
            totest.initialize_empty_array()
        # bad input
        with raises(AttributeError, match="'NoneType' object has no"+\
                    " attribute 'copy'"):
            totest.initialize_empty_array(None)
        # examples:
        test_array = totest.np.array([(1.0, "test"), (1.0, "test")],\
                                     dtype=([('mass', 'f8'), ('state', 'U8')]))
        empty_array = totest.initialize_empty_array(test_array)
        for i in range(len(empty_array)):
            assert totest.np.isnan(empty_array['mass'][i]) # nan float
            assert empty_array['state'][i] == 'nan'        # nan str

    def test_calculate_core_boundary(self, star_profile):
        # missing argument
        with raises(TypeError, match="missing 3 required positional"+\
                    " arguments: 'donor_mass', 'donor_star_state', and"+\
                    " 'profile'"):
            totest.calculate_core_boundary()
        # bad input
        with raises(ValueError, match="Not possible to calculate the core"+\
                    " boundary of the donor in CE"):
            totest.calculate_core_boundary(None, None, None)
        with warns(ApproximationWarning, match="Stellar profile columns were"+\
                   " not enough to calculate the core-envelope boundaries"+\
                   " for CE, entire star is now considered an envelope"):
            assert totest.calculate_core_boundary(None, None, None,\
                   core_element_fraction_definition=0.1) == -1
        with raises(AttributeError, match="'NoneType' object has no"+\
                    " attribute 'dtype'"):
            totest.calculate_core_boundary(None, "H-rich_Core_H_burning",\
             None, core_element_fraction_definition=0.1)
        # examples: get He core in H rich star
        assert totest.calculate_core_boundary(None, "H-rich_Core_H_burning",\
               star_profile, core_element_fraction_definition=0.1) == 2
        # examples: get CO core in H rich star
        assert totest.calculate_core_boundary(None, "H-rich_Core_H_burning",\
               star_profile, core_element_fraction_definition=0.3,\
               CO_core_in_Hrich_star=True) == 3
        # examples: get CO core in He star without a match
        assert totest.calculate_core_boundary(None,\
               "stripped_He_Core_He_burning", star_profile,\
               core_element_fraction_definition=0.1) == -1
        # examples: given core mass
        test_mass = totest.np.array([1.0, 0.7, 0.4, 0.1, 0.0])
        assert totest.calculate_core_boundary(test_mass, None, None,\
               mc1_i=0.1) == 4

    def test_period_evol_wind_loss(self):
        # missing argument
        with raises(TypeError, match="missing 4 required positional"+\
                    " arguments: 'M_current', 'M_init', 'Mcomp', and"+\
                    " 'P_init'"):
            totest.period_evol_wind_loss()
        # bad input
        with raises(TypeError) as error_info:
            totest.period_evol_wind_loss(None, None, None, None)
        assert error_info.value.args[0] == "unsupported operand type(s) for"+\
               " +: 'NoneType' and 'NoneType'"
        # examples:
        assert totest.period_evol_wind_loss(0.5, 0.5, 0.5, 0.5) == 0.5
        assert totest.period_evol_wind_loss(1.5, 0.5, 0.5, 0.5) == 0.5/4
        assert totest.period_evol_wind_loss(0.5, 1.5, 0.5, 0.5) == 0.5*4
        assert totest.period_evol_wind_loss(0.5, 0.5, 1.5, 0.5) == 0.5
        assert totest.period_evol_wind_loss(0.5, 0.5, 0.5, 1.5) == 1.5

    def test_separation_evol_wind_loss(self):
        # missing argument
        with raises(TypeError, match="missing 4 required positional"+\
                    " arguments: 'M_current', 'M_init', 'Mcomp', and"+\
                    " 'A_init'"):
            totest.separation_evol_wind_loss()
        # bad input
        with raises(TypeError) as error_info:
            totest.separation_evol_wind_loss(None, None, None, None)
        assert error_info.value.args[0] == "unsupported operand type(s) for"+\
               " +: 'NoneType' and 'NoneType'"
        # examples:
        assert totest.separation_evol_wind_loss(0.5, 0.5, 0.5, 0.5) == 0.5
        assert totest.separation_evol_wind_loss(1.5, 0.5, 0.5, 0.5) == 0.5/2
        assert totest.separation_evol_wind_loss(0.5, 1.5, 0.5, 0.5) == 0.5*2
        assert totest.separation_evol_wind_loss(0.5, 0.5, 1.5, 0.5) == 0.5
        assert totest.separation_evol_wind_loss(0.5, 0.5, 0.5, 1.5) == 1.5

    def test_period_change_stabe_MT(self):
        # missing argument
        with raises(TypeError, match="missing 4 required positional"+\
                    " arguments: 'period_i', 'Mdon_i', 'Mdon_f', and"+\
                    " 'Macc_i'"):
            totest.period_change_stabe_MT()
        # bad input
        with raises(TypeError) as error_info:
            totest.period_change_stabe_MT(None, None, None, None)
        assert error_info.value.args[0] == "unsupported operand type(s) for"+\
               " -: 'NoneType' and 'NoneType'"
        for a,b in zip([-1.0, 2.0, 0.0, 0.0], [0.0, 0.0, -1.0, 2.0]):
            with raises(ValueError) as error_info:
                totest.period_change_stabe_MT(0.5, 0.5, 0.5, 0.5, alpha=a,\
                                              beta=b)
            assert error_info.value.args[0] == "In period_change_stabe_MT,"+\
                   f" mass transfer efficiencies, alpha, beta: {a}, {b} are"+\
                   " not in the [0-1] range."
        with raises(ValueError, match="Donor gains mass from 0.5 to 1.5"):
            totest.period_change_stabe_MT(0.5, 0.5, 1.5, 0.5)
        # examples:
        assert totest.period_change_stabe_MT(0.5, 0.5, 0.5, 0.5) == 0.5
        assert totest.period_change_stabe_MT(1.5, 0.5, 0.5, 0.5) == 1.5
        assert totest.period_change_stabe_MT(0.5, 1.5, 0.5, 0.5) == 0.5
        assert totest.period_change_stabe_MT(0.5, 0.5, 0.5, 1.5) == 0.5
        assert totest.period_change_stabe_MT(0.5, 1.0, 0.5, 1.0) ==\
               approx(1.18518518519, abs=6e-12)
        assert totest.period_change_stabe_MT(0.5, 0.5, 0.5, 0.5, beta=1.0) ==\
               0.5
        assert totest.period_change_stabe_MT(1.5, 0.5, 0.5, 0.5, beta=1.0) ==\
               1.5
        assert totest.period_change_stabe_MT(0.5, 1.5, 0.5, 0.5, beta=1.0) ==\
               approx(0.13385261754, abs=6e-12)
        assert totest.period_change_stabe_MT(0.5, 0.5, 0.5, 1.5, beta=1.0) ==\
               0.5

    def test_linear_interpolation_between_two_cells(self, capsys):
        # missing argument
        with raises(TypeError, match="missing 3 required positional"+\
                    " arguments: 'array_y', 'array_x', and 'x_target'"):
            totest.linear_interpolation_between_two_cells()
        # bad input
        with raises(TypeError, match="'>=' not supported between instances"+\
                    " of 'NoneType' and 'NoneType'"):
            totest.linear_interpolation_between_two_cells(None, None, None)
        # examples: automatic interpolation
        test_array_x = totest.np.array([0.0, 0.2, 1.0])
        test_array_y = totest.np.array([1.0, 0.4, 0.0])
        assert totest.linear_interpolation_between_two_cells(test_array_y,\
               test_array_x, 0.1) == 0.7
        # examples: interpolation over multiple cells with verbose
        assert totest.linear_interpolation_between_two_cells(test_array_y,\
               test_array_x, 0.1, top=2, bot=0, verbose=True) == 0.9
        captured_out = capsys.readouterr().out
        for out in ["linear interpolation",\
                    "x_target, top, bot, len(array_x)",\
                    "x_top, x_bot, y_top, y_bot, y_target"]:
            assert out in captured_out
        # examples: extrapolation
        assert totest.linear_interpolation_between_two_cells(test_array_y,\
               test_array_x, 0.1, top=2) == 0.45
        assert totest.linear_interpolation_between_two_cells(test_array_y,\
               test_array_x, 0.1, bot=1) == 0.45
        with warns(ReplaceValueWarning, match="top=3 is too large, use last"+\
                   " element in array_y"):
            with warns(InterpolationWarning):
                assert totest.linear_interpolation_between_two_cells(\
                       test_array_y, test_array_x, 0.1, top=3) == 0.0
        with warns(ReplaceValueWarning, match="bot=-1 is too small, use"+\
                   " first element"):
            with warns(InterpolationWarning):
                assert totest.linear_interpolation_between_two_cells(\
                       test_array_y, test_array_x, 0.1, top=0) == 1.0
        with warns(InterpolationWarning, match="bot=2 is too large: use y at"+\
                   " top=1"):
            assert totest.linear_interpolation_between_two_cells(\
                   test_array_y, test_array_x, 0.1, top=1, bot=2) == 0.4
        with warns(InterpolationWarning, match="array_x too short, use y at"+\
                   " top=3"):
            assert totest.linear_interpolation_between_two_cells(\
                   totest.np.array([1.0, 0.4, 0.0, -1.0]), test_array_x, 0.1,\
                   top=3) == -1.0

    def test_calculate_lambda_from_profile(self, star_profile, capsys):
        # missing argument
        with raises(TypeError, match="missing 2 required positional"+\
                    " arguments: 'profile' and 'donor_star_state'"):
            totest.calculate_lambda_from_profile()
        # bad input
        with raises(AttributeError, match="'NoneType' object has no"+\
                    " attribute 'dtype'"):
            totest.calculate_lambda_from_profile(None, None)
        with raises(ValueError, match="state test not supported in CEE"):
            totest.calculate_lambda_from_profile(star_profile, "test",\
                                                 ind_core=1)
        with raises(ValueError, match="lambda_CE has a negative value"):
            with warns(EvolutionWarning, match="Ebind_i of the envelope is"+\
                        " positive"):
                totest.calculate_lambda_from_profile(star_profile,\
                                                     "H-rich_test", ind_core=1)
        # examples: get He core in H rich star
        lambda_CE, Mc, Rc = totest.calculate_lambda_from_profile(star_profile,\
                            "test", common_envelope_alpha_thermal=0.1)
        assert lambda_CE == approx(3.16722966582e+33, abs=6e+21)
        assert Mc == approx(0.0, abs=6e-12)
        assert Rc == approx(0.0, abs=6e-12)
        lambda_CE, Mc, Rc = totest.calculate_lambda_from_profile(star_profile,\
                            "test", ind_core=0,\
                            common_envelope_alpha_thermal=0.1)
        assert totest.np.isnan(lambda_CE)
        assert Mc == star_profile['mass'][0]
        assert Rc == star_profile['radius'][0]
        # examples: He core in H rich star
        lambda_CE, Mc, Rc = totest.calculate_lambda_from_profile(star_profile,\
                            "H-rich_test", ind_core=2,\
                            common_envelope_alpha_thermal=0.1, verbose=True)
        assert lambda_CE == approx(3.35757244514e+33, abs=6e+21)
        assert Mc == 0.1
        assert Rc == 0.1
        captured_out = capsys.readouterr().out
        for out in ["lambda_CE_top, lambda_CE_bot",\
                    "m1_i, radius1, len(profile) vs ind_core, mc1_i, rc1_i",\
                    "Ebind_i from profile ", "lambda_CE "]:
            assert out in captured_out
        # examples: CO core in H rich star
        lambda_CE, Mc, Rc = totest.calculate_lambda_from_profile(star_profile,\
                            "H-rich_test", ind_core=2,\
                            common_envelope_alpha_thermal=0.1,\
                            CO_core_in_Hrich_star=True)
        assert lambda_CE == approx(4.25658010995e+32, abs=6e+20)
        assert Mc == approx(1.03333333333, abs=6e-12)
        assert Rc == approx(1.03333333333, abs=6e-12)
        # examples: CO core in He rich star
        assert totest.calculate_lambda_from_profile(star_profile,\
               "H-rich_test", ind_core=2, common_envelope_alpha_thermal=0.1,\
               CO_core_in_Hrich_star=True) ==\
               totest.calculate_lambda_from_profile(star_profile,\
               "stripped_He_test", ind_core=2,\
               common_envelope_alpha_thermal=0.1)

    def test_get_mass_radius_dm_from_profile(self, star_profile):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'profile'"):
            totest.get_mass_radius_dm_from_profile()
        # bad input
        with raises(AttributeError, match="'NoneType' object has no"+\
                    " attribute 'dtype'"):
            totest.get_mass_radius_dm_from_profile(None)
        for f in star_profile.dtype.names:
            with raises(ValueError, match="One or many of the mass and/or"+\
                        " radius needed columns in the profile is not"+\
                        " provided for the CEE"):
                totest.get_mass_radius_dm_from_profile(star_profile[[f]])
        # examples: profile with radius
        with warns(ClassificationWarning, match="Donor mass from the binary"+\
                   " class object and the profile do not agree"):
            with warns(ClassificationWarning, match="Donor radius from the"+\
                       " binary class object and the profile do not agree"):
                d_mass, d_radius, d_dm =\
                 totest.get_mass_radius_dm_from_profile(star_profile[['mass',\
                 'radius']])
                assert totest.np.array_equal(d_mass, star_profile['mass'])
                assert totest.np.array_equal(d_radius, star_profile['radius'])
                assert totest.np.array_equal(d_dm, star_profile['dm'])
        # get total mass and radius
        M = star_profile['mass'].max()
        R = star_profile['radius'].max()
        # examples: profile with log_R
        d_mass, d_radius, d_dm = totest.get_mass_radius_dm_from_profile(\
                                 star_profile[['mass', 'log_R']], m1_i = M,\
                                 radius1 = R)
        assert totest.np.array_equal(d_mass, star_profile['mass'])
        assert totest.np.array_equal(d_radius, star_profile['radius'])
        assert totest.np.array_equal(d_dm, star_profile['dm'])
        # examples: profile with radius and dm (in grams)
        d_mass, d_radius, d_dm = totest.get_mass_radius_dm_from_profile(\
                                 star_profile[['mass', 'dm', 'radius']],\
                                 m1_i = M, radius1 = R)
        assert totest.np.array_equal(d_mass, star_profile['mass'])
        assert totest.np.array_equal(d_radius, star_profile['radius'])
        # convert Msun to grams
        d_dm = d_dm * totest.const.Msun
        assert totest.np.array_equal(d_dm, star_profile['dm'])

    def test_get_internal_energy_from_profile(self, star_profile, monkeypatch):
        def mock_calculate_H2recombination_energy(profile, tolerance=0.001):
            return totest.np.array(len(profile)*[1.0e+99])
        def mock_calculate_recombination_energy(profile, tolerance=0.001):
            return totest.np.array(len(profile)*[1.0e+99])
        # missing argument
        with raises(TypeError, match="missing 2 required positional"+\
                    " arguments: 'common_envelope_option_for_lambda' and"+\
                    " 'profile'"):
            totest.get_internal_energy_from_profile()
        # bad input
        with raises(AttributeError, match="'NoneType' object has no"+\
                    " attribute 'dtype'"):
            totest.get_internal_energy_from_profile(None, None)
        with raises(ValueError, match="unsupported:"+\
                    " common_envelope_option_for_lambda = test"):
            totest.get_internal_energy_from_profile("test",\
            star_profile[['energy']])
        bad_profile = totest.np.array([-1.0, -1.0, -1.0, -1.0],\
                      dtype=([('energy', 'f8')]))
        bad_profile2 = totest.np.empty((2,), dtype=[(f, 'f8') for f in [\
                       'x_mass_fraction_H', 'y_mass_fraction_He',\
                       'neutral_fraction_H', 'neutral_fraction_He',\
                       'avg_charge_He', 'energy']])
        bad_profile2['x_mass_fraction_H'] = totest.np.array([-1.0, -1.0])
        bad_profile2['y_mass_fraction_He'] = totest.np.array([-1.0, -1.0])
        bad_profile2['neutral_fraction_H'] = totest.np.array([0.5, 0.5])
        bad_profile2['neutral_fraction_He'] = totest.np.array([0.5, 0.5])
        bad_profile2['avg_charge_He'] = totest.np.array([0.5, 0.5])
        bad_profile2['energy'] = totest.np.array([0.5, 0.5])
        for CEOFL in ["lambda_from_profile_gravitational_plus_internal",\
                      "lambda_from_profile_gravitational_plus_internal_minus"+\
                      "_recombination"]:
            with raises(ValueError, match="CEE problem calculating internal"+\
                        " energy, giving negative values."):
                totest.get_internal_energy_from_profile(CEOFL, bad_profile)
            with monkeypatch.context() as m:
                m.setattr(totest, "calculate_H2recombination_energy",\
                          mock_calculate_H2recombination_energy)
                m.setattr(totest, "calculate_recombination_energy",\
                          mock_calculate_recombination_energy)
                with raises(ValueError) as error_info:
                    totest.get_internal_energy_from_profile(CEOFL,\
                    star_profile[['energy']])
                assert error_info.value.args[0] == "CEE problem calculating"+\
                       " recombination (and H2 recombination) energy,"+\
                       " remaining internal energy giving negative values."
        # examples: no internal energy
        assert totest.np.array_equal(totest.np.array([0, 0, 0, 0]),\
               totest.get_internal_energy_from_profile(\
               "lambda_from_profile_gravitational", star_profile[['radius']]))
        with warns(ApproximationWarning, match="Profile does not include"+\
                   " internal energy -- proceeding with "+\
                   "'lambda_from_profile_gravitational'"):
            assert totest.np.array_equal(totest.np.array([0, 0, 0, 0]),\
                   totest.get_internal_energy_from_profile("test",\
                   star_profile[['radius']]))
        # examples: with internal energy
        assert totest.np.allclose(totest.get_internal_energy_from_profile(\
               "lambda_from_profile_gravitational_plus_internal",\
               star_profile), totest.np.array([9.99850443e+15, 4.99893173e+15,\
               9.99786347e+14, 9.99786347e+12]))

    def test_calculate_H2recombination_energy(self, star_profile):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'profile'"):
            totest.calculate_H2recombination_energy()
        # bad input
        with raises(AttributeError, match="'NoneType' object has no"+\
                    " attribute 'dtype'"):
            totest.calculate_H2recombination_energy(None)
        with raises(ValueError, match="CEE problem calculating H2"+\
                    " recombination energy, giving negative values"):
            bad_profile = totest.np.array([-1.0, -1.0, -1.0, -1.0],\
                          dtype=([('x_mass_fraction_H', 'f8')]))
            totest.calculate_H2recombination_energy(bad_profile)
        # examples: profile with Hydrogen mass fraction
        assert totest.np.allclose(totest.np.array([1.49557421e+12,\
               1.06826729e+12, 2.13653458e+11, 2.13653458e+09]),\
               totest.calculate_H2recombination_energy(\
               star_profile[['x_mass_fraction_H']]))
        # examples: profile without Hydrogen mass fraction
        with warns(ApproximationWarning, match="Profile does not include"+\
                   " Hydrogen mass fraction calculate H2 recombination"+\
                   " energy -- H2 recombination energy is assumed 0"):
            assert totest.np.array_equal(totest.np.array([0, 0, 0, 0]),\
                   totest.calculate_H2recombination_energy(\
                   star_profile[['radius']]))

    def test_calculate_recombination_energy(self, star_profile):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'profile'"):
            totest.calculate_recombination_energy()
        # bad input
        with raises(AttributeError, match="'NoneType' object has no"+\
                    " attribute 'dtype'"):
            totest.calculate_recombination_energy(None)
        with raises(ValueError, match="CEE problem calculating"+\
                    " recombination energy, giving negative values"):
            bad_profile = totest.np.empty((2,), dtype=[(f, 'f8') for f in [\
                          'x_mass_fraction_H', 'y_mass_fraction_He',\
                          'neutral_fraction_H', 'neutral_fraction_He',\
                          'avg_charge_He']])
            bad_profile['x_mass_fraction_H'] = totest.np.array([-1.0, -1.0])
            bad_profile['y_mass_fraction_He'] = totest.np.array([-1.0, -1.0])
            bad_profile['neutral_fraction_H'] = totest.np.array([0.5, 0.5])
            bad_profile['neutral_fraction_He'] = totest.np.array([0.5, 0.5])
            bad_profile['avg_charge_He'] = totest.np.array([0.5, 0.5])
            totest.calculate_recombination_energy(bad_profile)
        # examples: profile with mass fraction and ionization information
        assert totest.np.allclose(totest.np.array([0.00000000e+00,\
               4.73639308e+12, 7.53791976e+12, 3.29724895e+12]),\
               totest.calculate_recombination_energy(\
               star_profile[['x_mass_fraction_H', 'y_mass_fraction_He', 'neutral_fraction_H', 'neutral_fraction_He', 'avg_charge_He']]))
        # examples: profile without mass fraction and ionization information
        with warns(ApproximationWarning, match="Profile does not include"+\
                   " mass fractions and ionizations of elements to calculate"+\
                   " recombination energy -- recombination energy is assumed"+\
                   " 0"):
            assert totest.np.array_equal(totest.np.array([0, 0, 0, 0]),\
                   totest.calculate_recombination_energy(\
                   star_profile[['radius']]))

    def test_profile_recomb_energy(self):
        # missing argument
        with raises(TypeError, match="missing 5 required positional"+\
                    " arguments: 'x_mass_fraction_H', 'y_mass_fraction_He',"+\
                    " 'frac_HII', 'frac_HeII', and 'frac_HeIII'"):
            totest.profile_recomb_energy()
        # bad input
        with raises(TypeError) as error_info:
            totest.profile_recomb_energy(None, None, None, None, None)
        assert error_info.value.args[0] == "unsupported operand type(s) for"+\
               " *: 'float' and 'NoneType'"
        # examples:
        assert totest.profile_recomb_energy(0.5, 0.5, 0.5, 0.5, 0.5) ==\
               approx(9.49756867371e+12, abs=6e+0)
        assert totest.profile_recomb_energy(0.1, 0.5, 0.5, 0.5, 0.5) ==\
               approx(6.89384394360e+12, abs=6e+0)
        assert totest.profile_recomb_energy(0.5, 0.1, 0.5, 0.5, 0.5) ==\
               approx(4.50323846485e+12, abs=6e+0)
        assert totest.profile_recomb_energy(0.5, 0.5, 0.1, 0.5, 0.5) ==\
               approx(6.89384394360e+12, abs=6e+0)
        assert totest.profile_recomb_energy(0.5, 0.5, 0.5, 0.1, 0.5) ==\
               approx(8.31217893947e+12, abs=6e+0)
        assert totest.profile_recomb_energy(0.5, 0.5, 0.5, 0.5, 0.1) ==\
               approx(5.68862819909e+12, abs=6e+0)

    def test_calculate_binding_energy(self, star_profile, capsys):
        # missing argument
        with raises(TypeError, match="missing 7 required positional"+\
                    " arguments: 'donor_mass', 'donor_radius', 'donor_dm',"+\
                    " 'specific_internal_energy', 'ind_core',"+\
                    " 'factor_internal_energy', and 'verbose'"):
            totest.calculate_binding_energy()
        # bad input
        with raises(TypeError, match="'NoneType' object cannot be"+\
                    " interpreted as an integer"):
            totest.calculate_binding_energy(None, None, None, None, None,\
                                            None, None)
        with raises(ValueError, match="CEE problem calculating"+\
                    " gravitational energy, giving positive values"):
            totest.calculate_binding_energy(-1.0*star_profile['mass'],\
                                            star_profile['radius'],\
                                            star_profile['dm'],\
                                            star_profile['energy'], 1, 0.1,\
                                            False)
        # examples:
        with warns(EvolutionWarning, match="Ebind_i of the envelope is"+\
                   " positive"):
            assert totest.calculate_binding_energy(star_profile['mass'],\
                   -1.0e+60*star_profile['radius'], star_profile['dm'],\
                   star_profile['energy'], 1, 0.5, False) == 4.973e+48
            assert totest.calculate_binding_energy(star_profile['mass'],\
                   star_profile['radius'], star_profile['dm'],\
                   star_profile['energy'], 2, 0.5, False) ==\
                   approx(3.547071313426e+48, abs=6e+36)
        assert totest.calculate_binding_energy(star_profile['mass'],\
               star_profile['radius'], star_profile['dm'],\
               star_profile['energy'], 1, 0.1, True) ==\
               approx(-9.02693714763e+47, abs=6e+35)
        captured_out = capsys.readouterr().out
        for out in ["integration of gravitational energy surface to core"+\
                    " [Grav_energy], integration of internal energy surface"+\
                    " to core [U_i] (0 if not taken into account)",\
                    "Ebind = Grav_energy + factor_internal_energy*U_i  :  "]:
            assert out in captured_out

    def test_calculate_Mejected_for_integrated_binding_energy(self,\
        star_profile):
        # missing argument
        with raises(TypeError, match="missing 4 required positional"+\
                    " arguments: 'profile', 'Ebind_threshold', 'mc1_i', and"+\
                    " 'rc1_i'"):
            totest.calculate_Mejected_for_integrated_binding_energy()
        # bad input
        with raises(AttributeError, match="'NoneType' object has no"+\
                    " attribute 'dtype'"):
            totest.calculate_Mejected_for_integrated_binding_energy(None,\
             None, None, None)
        # examples:
        # get total mass and radius
        M = star_profile['mass'].max()
        R = star_profile['radius'].max()
        assert totest.calculate_Mejected_for_integrated_binding_energy(\
               star_profile, 1.0, 1.0, 1.0, m1_i=M, radius1=R) == 0.0
        with warns(EvolutionWarning, match="partial mass ejected is greater"+\
                   " than the envelope mass"):
            assert totest.calculate_Mejected_for_integrated_binding_energy(\
                   star_profile, 1.0e+50, 1.0, 1.0, m1_i=M, radius1=R) == 0.0

    def test_convert_metallicity_to_string(self):
        # missing argument
        with raises(TypeError, match="missing 1 required positional"+\
                    " argument: 'Z'"):
            totest.convert_metallicity_to_string()
        # bad input
        with raises(ValueError, match="Metallicity None not supported!"+\
                    " Available metallicities in POSYDON v2 are"):
            totest.convert_metallicity_to_string(None)
        # examples:
        for Z,s in zip([2e+00, 1e+00, 4.5e-01, 2e-01, 1e-01, 1e-02, 1e-03,\
                        1e-04], ['2e+00', '1e+00', '4.5e-01', '2e-01',\
                        '1e-01', '1e-02', '1e-03', '1e-04']):
            assert totest.convert_metallicity_to_string(Z) == s

    def test_rotate(self):
        # missing argument
        with raises(TypeError, match="missing 2 required positional"+\
                    " arguments: 'axis' and 'angle'"):
            totest.rotate()
        # bad input
        with raises(TypeError, match="object of type 'NoneType' has no len()"):
            totest.rotate(None, None)
        with raises(ValueError, match="axis should be of dimension 3"):
            totest.rotate(totest.np.array([0]), 0.0)
        with raises(ValueError, match="axis is a point"):
            totest.rotate(totest.np.array([0, 0, 0]), 0.0)
        # examples:
        for x in range(3):
            for y in range(3):
                for z in range(3):
                    if x!=0 or y!=0 or z!=0:
                        # angle=0 -> no rotation
                        assert totest.np.array_equal(totest.rotate(\
                               totest.np.array([x, y, z]), 0.0),\
                               totest.np.array([[1, 0, 0], [0, 1, 0],\
                                                [0, 0, 1]]))
                        # inverse rotation around inverted axis
                        assert totest.np.array_equal(totest.rotate(\
                               totest.np.array([-x, -y, -z]), -1.0),\
                               totest.rotate(totest.np.array([x, y, z]), 1.0))
        # examples: angle=1.0
        s = totest.np.sin(1.0)
        c = totest.np.cos(1.0)
        assert totest.np.array_equal(totest.rotate(totest.np.array([1, 0, 0]),\
               1.0), totest.np.array([[1, 0, 0], [0, c, -s], [0, s, c]]))
        assert totest.np.array_equal(totest.rotate(totest.np.array([0, 1, 0]),\
               1.0), totest.np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]]))
        assert totest.np.array_equal(totest.rotate(totest.np.array([0, 0, 1]),\
               1.0), totest.np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]]))
        assert totest.np.allclose(totest.rotate(totest.np.array([1, 1, 0]),\
               1.0), totest.np.array([[0.77015115, 0.22984885, 0.59500984],\
               [0.22984885, 0.77015115, -0.59500984], [-0.59500984,\
                0.59500984, 0.54030231]]))
        assert totest.np.allclose(totest.rotate(totest.np.array([1, 0, 1]),\
               1.0), totest.np.array([[0.77015115,-0.59500984, 0.22984885],\
               [0.59500984, 0.54030231, -0.59500984], [0.22984885,\
                0.59500984, 0.77015115]]))
        assert totest.np.allclose(totest.rotate(totest.np.array([0, 1, 1]),\
               1.0), totest.np.array([[0.54030231,-0.59500984, 0.59500984],\
               [0.59500984, 0.77015115, 0.22984885], [-0.59500984,\
                0.22984885, 0.77015115]]))


class TestPchipInterpolator2:
    @fixture
    def PchipInterpolator2(self):
        # initialize an instance of the class with defaults
        return totest.PchipInterpolator2([0.0, 1.0], [1.0, 0.0])

    @fixture
    def PchipInterpolator2_True(self):
        # initialize an instance of the class with defaults
        return totest.PchipInterpolator2([0.0, 1.0], [-1.0, 0.0],\
                                         positive=True)

    # test the PchipInterpolator2 class
    def test_init(self, PchipInterpolator2, PchipInterpolator2_True):
        assert isroutine(PchipInterpolator2.__init__)
        # check that the instance is of correct type and all code in the
        # __init__ got executed: the elements are created and initialized
        assert isinstance(PchipInterpolator2, totest.PchipInterpolator2)
        assert isinstance(PchipInterpolator2.interpolator,\
                          totest.PchipInterpolator)
        assert PchipInterpolator2.positive == False
        assert PchipInterpolator2_True.positive

    def test_call(self, PchipInterpolator2, PchipInterpolator2_True):
        assert isroutine(PchipInterpolator2.__call__)
        assert PchipInterpolator2(0.1) == 0.9
        assert PchipInterpolator2_True(0.1) == 0.0
