"""Unit tests of posydon/popsyn/binarypopulation.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

# import the module which will be tested
import posydon.popsyn.binarypopulation as totest

# aliases
np = totest.np
pd = totest.pd
os = totest.os

from inspect import isclass, isroutine

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import approx, fixture, raises

from posydon.binary_evol.binarystar import BinaryStar
from posydon.binary_evol.singlestar import SingleStar
from posydon.binary_evol.simulationproperties import SimulationProperties


# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['BinaryPopulation', 'PopulationManager',
                    'BinaryGenerator',
                    'saved_ini_parameters',
                    'HISTORY_MIN_ITEMSIZE', 'ONELINE_MIN_ITEMSIZE',
                    'STEP_NAMES_LOADING_GRIDS',
                    'default_kwargs',
                    '__authors__', '__credits__',
                    '__builtins__', '__cached__', '__doc__', '__file__',
                    '__loader__', '__name__', '__package__', '__spec__',
                    'np', 'pd', 'os', 'atexit', 'signal', 'traceback',
                    'psutil', 'tqdm',
                    'posydon', 'BinaryStar', 'SimulationProperties',
                    'SingleStar', 'properties_massless_remnant',
                    'generate_independent_samples',
                    'binarypop_kwargs_from_ini', 'simprop_kwargs_from_ini',
                    'get_kick_samples_from_file', 'get_samples_from_file',
                    'get_formation_times',
                    'orbital_period_from_separation',
                    'orbital_separation_from_period',
                    'set_binary_to_failed',
                    'Zsun', 'POSYDONError',
                    'Catch_POSYDON_Warnings', 'Pwarn',
                    ]
        totest_elements = set(dir(totest))
        missing_in_test = set(elements) - totest_elements
        assert len(missing_in_test) == 0, "There are missing objects in "\
                                          +f"{totest.__name__}: "\
                                          +f"{missing_in_test}. Please "\
                                          +"check, whether they have been "\
                                          +"removed on purpose and update "\
                                          +"this unit test."
        new_in_test = totest_elements - set(elements)
        assert len(new_in_test) == 0, "There are new objects in "\
                                      +f"{totest.__name__}: {new_in_test}. "\
                                      +"Please check, whether they have been "\
                                      +"added on purpose and update this "\
                                      +"unit test."


class TestBinaryGenerator:

    @fixture
    def generator(self):
        """Create a BinaryGenerator with seeded RNG."""
        rng = np.random.default_rng(seed=42)
        return totest.BinaryGenerator(RNG=rng, metallicity=1.0,
                                      star_formation='burst',
                                      max_simulation_time=13.8e9)

    @fixture
    def kick_csv(self, tmp_path):
        """CSV with kick columns for file-based sampling."""
        df = pd.DataFrame({
            'm1': [10.0, 20.0],
            'm2': [5.0, 10.0],
            'orbital_period': [10.0, 20.0],
            'eccentricity': [0.1, 0.2],
            's1_natal_kick_velocity': [100.0, 200.0],
            's1_natal_kick_azimuthal_angle': [0.5, 1.0],
            's1_natal_kick_polar_angle': [0.3, 0.6],
            's1_natal_kick_mean_anomaly': [0.1, 0.2],
            's2_natal_kick_velocity': [50.0, 150.0],
            's2_natal_kick_azimuthal_angle': [0.4, 0.8],
            's2_natal_kick_polar_angle': [0.2, 0.5],
            's2_natal_kick_mean_anomaly': [0.05, 0.15],
        })
        path = os.path.join(tmp_path, "kicks.csv")
        df.to_csv(path, index=False)
        return str(path)

    def test_init_default_rng(self):
        gen = totest.BinaryGenerator()
        assert isinstance(gen.RNG, np.random.Generator)
        assert gen.entropy is not None
        assert gen._num_gen == 0

    def test_init_custom_rng(self, generator):
        assert isinstance(generator.RNG, np.random.Generator)
        assert generator._num_gen == 0
        assert generator.star_formation == 'burst'
        assert generator.Z_div_Zsun == 1.0

    def test_init_bad_rng(self):
        with raises(AssertionError):
            totest.BinaryGenerator(RNG="not_a_generator")

    def test_draw_initial_samples_separation(self, generator):
        output = generator.draw_initial_samples(
            orbital_scheme='separation', number_of_binaries=5)
        assert len(output['S1_mass']) == 5
        assert len(output['separation']) == 5
        assert len(output['orbital_period']) == 5
        assert all(output['binary_index'] == np.arange(0, 5))
        assert generator._num_gen == 5

    def test_draw_initial_samples_period(self, generator):
        output = generator.draw_initial_samples(
            orbital_scheme='period', number_of_binaries=3)
        assert len(output['S1_mass']) == 3
        assert all(np.isfinite(output['S1_mass']))

    def test_draw_initial_samples_bad_scheme(self, generator):
        with raises(ValueError, match="Allowed orbital schemes"):
            generator.draw_initial_samples(orbital_scheme='invalid')

    def test_draw_initial_samples_no_number_of_binaries(self, generator):
        """When number_of_binaries not in kwargs, defaults to 1 for kicks."""
        output = generator.draw_initial_samples(orbital_scheme='separation')
        assert len(output['S1_mass']) == 1

    def test_draw_initial_samples_from_file(self, kick_csv):
        """Kick values read from file."""
        rng = np.random.default_rng(seed=42)
        gen = totest.BinaryGenerator(
            RNG=rng, metallicity=1.0,
            sampler=totest.get_samples_from_file,
            star_formation='burst',
            max_simulation_time=13.8e9)
        output = gen.draw_initial_samples(
            orbital_scheme='period',
            number_of_binaries=2,
            read_samples_from_file=kick_csv)
        assert output['S1_natal_kick_velocity'][0] == 100.0
        assert output['S2_natal_kick_velocity'][1] == 150.0

    def test_draw_initial_binary(self, generator):
        binary = generator.draw_initial_binary(
            orbital_scheme='separation', metallicity=1.0)
        assert isinstance(binary, BinaryStar)
        assert isinstance(binary.star_1, SingleStar)
        assert isinstance(binary.star_2, SingleStar)
        assert binary.star_1.mass > 0
        assert binary.event == 'ZAMS'
        assert binary.state == 'detached'

    def test_draw_initial_binary_single_star(self):
        """With binary_fraction_const=0, all draws produce single stars."""
        rng = np.random.default_rng(seed=42)
        gen = totest.BinaryGenerator(
            RNG=rng, metallicity=1.0,
            star_formation='burst',
            max_simulation_time=13.8e9,
            binary_fraction_const=0.0,
            binary_fraction_scheme='const')
        binary = gen.draw_initial_binary(
            orbital_scheme='separation', metallicity=1.0)
        assert isinstance(binary, BinaryStar)
        assert binary.state == 'initially_single_star'
        assert np.isnan(binary.separation)
        assert np.isnan(binary.orbital_period)
        assert np.isnan(binary.eccentricity)

    def test_draw_initial_binary_with_index(self, generator):
        binary = generator.draw_initial_binary(
            orbital_scheme='separation', metallicity=1.0, index=99)
        assert binary.index == 99

    def test_reset_rng(self, generator):
        generator.draw_initial_samples(
            orbital_scheme='separation', number_of_binaries=5)
        assert generator._num_gen == 5
        generator.reset_rng()
        assert generator._num_gen == 0

    def test_get_original_rng(self, generator):
        rng = generator.get_original_rng()
        assert isinstance(rng, np.random.Generator)

    def test_get_binary_by_iter(self, generator):
        binary = generator.get_binary_by_iter(
            n=3, orbital_scheme='separation', metallicity=1.0)
        assert isinstance(binary, BinaryStar)
        assert generator._num_gen == 0

    def test_get_binary_by_iter_zero(self, generator):
        """n=0 skips the warmup sampling."""
        binary = generator.get_binary_by_iter(
            n=0, orbital_scheme='separation', metallicity=1.0)
        assert isinstance(binary, BinaryStar)

    def test_num_gen_increments(self, generator):
        generator.draw_initial_samples(
            orbital_scheme='separation', number_of_binaries=3)
        assert generator._num_gen == 3
        output = generator.draw_initial_samples(
            orbital_scheme='separation', number_of_binaries=1)
        assert output['binary_index'][0] == 3


class TestPopulationManager:

    @fixture
    def manager(self):
        """Create a PopulationManager with minimal kwargs."""
        return totest.PopulationManager(
            RNG=np.random.default_rng(seed=42),
            metallicity=1.0,
            star_formation='burst',
            max_simulation_time=13.8e9)

    @fixture
    def dummy_binary(self):
        """Create a minimal BinaryStar."""
        s1 = SingleStar(mass=10.0, state='H-rich_Core_H_burning',
                        metallicity=1.0)
        s2 = SingleStar(mass=5.0, state='H-rich_Core_H_burning',
                        metallicity=1.0)
        return BinaryStar(star_1=s1, star_2=s2, index=0,
                          state='detached', event='ZAMS',
                          time=0.0, separation=100.0,
                          orbital_period=10.0, eccentricity=0.0)

    @fixture
    def failed_binary(self):
        """Create a binary with FAILED event."""
        s1 = SingleStar(mass=10.0, state='H-rich_Core_H_burning',
                        metallicity=1.0)
        s2 = SingleStar(mass=5.0, state='H-rich_Core_H_burning',
                        metallicity=1.0)
        b = BinaryStar(star_1=s1, star_2=s2, index=1,
                       state='detached', event='ZAMS',
                       time=0.0, separation=100.0,
                       orbital_period=10.0, eccentricity=0.0)
        b.event = 'FAILED'
        return b

    def test_init(self, manager):
        assert manager.binaries == []
        assert manager.indices == []
        assert isinstance(manager.binary_generator, totest.BinaryGenerator)

    def test_init_with_filename(self):
        mgr = totest.PopulationManager(
            file_name='test.h5',
            RNG=np.random.default_rng(seed=42),
            metallicity=1.0,
            star_formation='burst',
            max_simulation_time=13.8e9)
        assert mgr.store_file == 'test.h5'

    def test_init_with_file_sampler(self, tmp_path):
        csv_path = os.path.join(tmp_path, "samples.csv")
        df = pd.DataFrame({
            'm1': [10.0], 'm2': [5.0],
            'orbital_period': [10.0], 'eccentricity': [0.0],
        })
        df.to_csv(csv_path, index=False)
        mgr = totest.PopulationManager(
            read_samples_from_file=str(csv_path),
            RNG=np.random.default_rng(seed=42),
            metallicity=1.0,
            star_formation='burst',
            max_simulation_time=13.8e9)
        assert mgr.binary_generator.sampler is not totest.generate_independent_samples

    def test_append_single(self, manager, dummy_binary):
        manager.append(dummy_binary)
        assert len(manager.binaries) == 1
        assert manager.indices == [0]

    def test_append_list(self, manager, dummy_binary):
        manager.append([dummy_binary])
        assert len(manager.binaries) == 1

    def test_append_invalid(self, manager):
        with raises(ValueError, match="Must be BinaryStar"):
            manager.append("not_a_binary")

    def test_remove_single(self, manager, dummy_binary):
        manager.append(dummy_binary)
        manager.remove(dummy_binary)
        assert len(manager.binaries) == 0

    def test_remove_list(self, manager, dummy_binary):
        manager.append(dummy_binary)
        manager.remove([dummy_binary])
        assert len(manager.binaries) == 0

    def test_remove_invalid(self, manager):
        with raises(ValueError, match="Must be BinaryStar"):
            manager.remove("not_a_binary")

    def test_clear_dfs(self, manager):
        manager.history_dfs = [pd.DataFrame()]
        manager.oneline_dfs = [pd.DataFrame()]
        manager.clear_dfs()
        assert manager.history_dfs == []
        assert manager.oneline_dfs == []

    def test_generate(self, manager):
        binary = manager.generate(orbital_scheme='separation',
                                  metallicity=1.0)
        assert isinstance(binary, BinaryStar)
        assert len(manager.binaries) == 1

    def test_to_df_empty(self, manager):
        assert manager.to_df() is None

    def test_to_df_with_binaries(self, manager):
        manager.generate(orbital_scheme='separation', metallicity=1.0)
        result = manager.to_df()
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0

    def test_to_df_with_selection_accept(self, manager):
        manager.generate(orbital_scheme='separation', metallicity=1.0)
        result = manager.to_df(selection_function=lambda b: True)
        assert isinstance(result, pd.DataFrame)

    def test_to_df_with_selection_reject(self, manager):
        manager.generate(orbital_scheme='separation', metallicity=1.0)
        result = manager.to_df(selection_function=lambda b: False)
        assert result is None

    def test_to_df_with_history_dfs(self, manager):
        dummy_df = pd.DataFrame({'state': ['detached'], 'time': [0.0]},
                                index=[0])
        manager.history_dfs = [dummy_df]
        result = manager.to_df()
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 1

    def test_to_oneline_df_empty(self, manager):
        assert manager.to_oneline_df() is None

    def test_to_oneline_df_with_binaries(self, manager):
        manager.generate(orbital_scheme='separation', metallicity=1.0)
        result = manager.to_oneline_df()
        assert isinstance(result, pd.DataFrame)

    def test_to_oneline_df_with_selection(self, manager):
        manager.generate(orbital_scheme='separation', metallicity=1.0)
        result = manager.to_oneline_df(selection_function=lambda b: True)
        assert isinstance(result, pd.DataFrame)

    def test_to_oneline_df_with_oneline_dfs(self, manager):
        dummy_df = pd.DataFrame({'state_i': ['detached']}, index=[0])
        manager.oneline_dfs = [dummy_df]
        result = manager.to_oneline_df()
        assert isinstance(result, pd.DataFrame)

    def test_find_failed_empty(self, manager):
        assert manager.find_failed() is None

    def test_find_failed_with_binaries(self, manager, dummy_binary,
                                       failed_binary):
        manager.append(dummy_binary)
        manager.append(failed_binary)
        result = manager.find_failed()
        assert len(result) == 1
        assert result[0].event == 'FAILED'

    def test_find_failed_with_dfs_found(self, manager):
        failed_df = pd.DataFrame({'event': ['ZAMS', 'FAILED'],
                                  'time': [0.0, 1.0]}, index=[0, 0])
        manager.history_dfs = [failed_df]
        result = manager.find_failed()
        assert isinstance(result, pd.DataFrame)

    def test_find_failed_with_dfs_none_found(self, manager):
        ok_df = pd.DataFrame({'event': ['ZAMS', 'END'],
                              'time': [0.0, 1.0]}, index=[0, 0])
        manager.history_dfs = [ok_df]
        result = manager.find_failed()
        assert result == []

    def test_breakdown_to_df(self, manager):
        binary = manager.generate(orbital_scheme='separation',
                                  metallicity=1.0)
        assert len(manager.binaries) == 1
        manager.breakdown_to_df(binary)
        assert len(manager.binaries) == 0
        assert len(manager.history_dfs) == 1
        assert len(manager.oneline_dfs) == 1


class TestBinaryPopulation:

    @fixture
    def pop(self):
        """Create a minimal BinaryPopulation."""
        return totest.BinaryPopulation(
            number_of_binaries=3,
            metallicity=1.0,
            star_formation='burst',
            max_simulation_time=13.8e9,
            entropy=12345)

    def test_init_basic(self, pop):
        assert pop.number_of_binaries == 3
        assert pop.metallicity == 1.0
        assert isinstance(pop.population_properties, SimulationProperties)
        assert isinstance(pop.manager, totest.PopulationManager)
        assert isinstance(pop.RNG, np.random.Generator)
        assert pop.comm is None
        assert pop.JOB_ID is None

    def test_init_metallicity_from_list(self):
        pop = totest.BinaryPopulation(
            number_of_binaries=1,
            metallicities=[0.5, 1.0, 2.0],
            metallicity_index=1,
            star_formation='burst',
            max_simulation_time=13.8e9,
            entropy=42)
        assert pop.metallicity == 1.0

    def test_init_mpi_and_jobarray_incompatible(self):
        class FakeComm:
            def Get_rank(self): return 0
            def Get_size(self): return 2
        with raises(ValueError, match="MPI and Job array runs are not compatible"):
            totest.BinaryPopulation(
                number_of_binaries=1,
                comm=FakeComm(),
                JOB_ID=123,
                star_formation='burst',
                max_simulation_time=13.8e9,
                entropy=42)

    def test_init_mpi_no_entropy(self):
        class FakeComm:
            def Get_rank(self): return 0
            def Get_size(self): return 2
        with raises(ValueError, match="requires an entropy value"):
            totest.BinaryPopulation(
                number_of_binaries=1,
                comm=FakeComm(),
                star_formation='burst',
                max_simulation_time=13.8e9)

    def test_init_mpi_with_entropy(self):
        class FakeComm:
            def Get_rank(self): return 0
            def Get_size(self): return 2
        pop = totest.BinaryPopulation(
            number_of_binaries=10,
            comm=FakeComm(),
            star_formation='burst',
            max_simulation_time=13.8e9,
            entropy=42)
        assert isinstance(pop.RNG, np.random.Generator)

    def test_init_job_array_with_entropy(self):
        pop = totest.BinaryPopulation(
            number_of_binaries=10,
            JOB_ID=100,
            RANK=0,
            size=2,
            star_formation='burst',
            max_simulation_time=13.8e9,
            entropy=42)
        assert pop.JOB_ID == 100

    def test_init_job_array_no_entropy(self):
        """JOB_ID without entropy uses JOB_ID as seed."""
        pop = totest.BinaryPopulation(
            number_of_binaries=10,
            JOB_ID=100,
            RANK=0,
            size=2,
            star_formation='burst',
            max_simulation_time=13.8e9)
        assert pop.JOB_ID == 100
        assert isinstance(pop.RNG, np.random.Generator)

    def test_from_ini(self, monkeypatch, tmp_path):
        def mock_binarypop_kwargs(path, verbose=False):
            return {
                'number_of_binaries': 5,
                'metallicity': 1.0,
                'metallicities': [1.0],
                'star_formation': 'burst',
                'max_simulation_time': 13.8e9,
                'entropy': 99,
            }
        def mock_simprop_kwargs(path):
            return {}

        monkeypatch.setattr(totest, 'binarypop_kwargs_from_ini',
                            mock_binarypop_kwargs)
        monkeypatch.setattr(totest, 'simprop_kwargs_from_ini',
                            mock_simprop_kwargs)

        pop = totest.BinaryPopulation.from_ini(str(tmp_path / "fake.ini"))
        assert pop.number_of_binaries == 5
        assert isinstance(pop.population_properties, SimulationProperties)

    def test_close(self, pop):
        pop.close()

    def test_getstate(self, pop):
        state = pop.__getstate__()
        assert isinstance(state, dict)
        assert state['comm'] is None

    def test_getstate_with_steps_loaded(self, pop):
        pop.population_properties.steps_loaded = True
        state = pop.__getstate__()
        assert state['comm'] is None
        assert not pop.population_properties.steps_loaded