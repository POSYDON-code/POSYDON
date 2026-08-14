"""Unit tests of posydon/popsyn/sample_from_file.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

# import the module which will be tested
import posydon.popsyn.sample_from_file as totest

# aliases
os = totest.os
np = totest.np
pd = totest.pd

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import approx, fixture, raises, warns


# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['infer_key', 'get_samples_from_file',
                    'get_kick_samples_from_file', '__authors__',
                    '__builtins__', '__cached__', '__doc__', '__file__',
                    '__loader__', '__name__', '__package__', '__spec__',
                    'os', 'np', 'pd', 'Pwarn',
                    'generate_eccentricities', 'generate_orbital_periods',
                    'generate_orbital_separations', 'generate_primary_masses',
                    'generate_secondary_masses',
                    'PRIMARY_MASS_NAMES', 'SECONDARY_MASS_NAMES',
                    'PERIOD_NAMES', 'SEPARATION_NAMES', 'ECCENTRICITY_NAMES',
                    'PRIMARY_KICK_VELOCITY_NAMES', 'SECONDARY_KICK_VELOCITY_NAMES',
                    'PRIMARY_KICK_AZIMUTHAL_ANGLE_NAMES', 'SECONDARY_KICK_AZIMUTHAL_ANGLE_NAMES',
                    'PRIMARY_KICK_POLAR_ANGLE_NAMES', 'SECONDARY_KICK_POLAR_ANGLE_NAMES',
                    'PRIMARY_KICK_MEAN_ANOMALY_NAMES', 'SECONDARY_KICK_MEAN_ANOMALY_NAMES',
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

class TestFunctions:

    @fixture
    def full_csv(self, tmp_path):
        """CSV with all binary columns."""
        df = pd.DataFrame({
            'm1': [10.0, 20.0, 30.0],
            'm2': [5.0, 10.0, 15.0],
            'orbital_period': [1.0, 10.0, 100.0],
            'orbital_separation': [50.0, 100.0, 200.0],
            'eccentricity': [0.0, 0.1, 0.2],
        })
        path = os.path.join(tmp_path, "full.csv")
        df.to_csv(path, index=False)
        return path

    @fixture
    def minimal_csv(self, tmp_path):
        """CSV with no recognized binary columns."""
        df = pd.DataFrame({
            'col_a': [1.0, 2.0],
            'col_b': [3.0, 4.0],
        })
        path = os.path.join(tmp_path, "minimal.csv")
        df.to_csv(path, index=False)
        return path

    @fixture
    def kick_csv(self, tmp_path):
        """CSV with all kick columns."""
        df = pd.DataFrame({
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
        return path

    @fixture
    def no_kick_csv(self, tmp_path):
        """CSV with no kick columns."""
        df = pd.DataFrame({
            'col_a': [1.0, 2.0],
        })
        path = os.path.join(tmp_path, "no_kicks.csv")
        df.to_csv(path, index=False)
        return path

    # --- infer_key ---

    def test_infer_key(self):
        # exact match
        assert totest.infer_key(
            available_keys=['m1', 'period'],
            allowed_keys=['m1', 'm2']) == 'm1'

        # case-insensitive match
        assert totest.infer_key(
            available_keys=['M1', 'Period'],
            allowed_keys=['m1']) == 'M1'

        # no match
        assert totest.infer_key(
            available_keys=['col_a', 'col_b'],
            allowed_keys=['m1', 'm2']) == ''

        # empty inputs
        assert totest.infer_key(available_keys=[], allowed_keys=['m1']) == ''
        assert totest.infer_key(available_keys=['m1'], allowed_keys=[]) == ''

    # --- get_samples_from_file ---

    def test_get_samples_from_file_missing_kwarg(self):
        with raises(KeyError, match="no 'read_samples_from_file'"):
            totest.get_samples_from_file(orbital_scheme='period')

    def test_get_samples_from_file_not_found(self):
        with raises(FileNotFoundError, match="not found"):
            totest.get_samples_from_file(
                orbital_scheme='period',
                read_samples_from_file='nonexistent.csv')

    def test_get_samples_from_file_bad_scheme(self, full_csv):
        with raises(ValueError, match="Allowed orbital schemes are separation or period."):
            totest.get_samples_from_file(
                orbital_scheme='invalid',
                read_samples_from_file=full_csv)

    def test_get_samples_from_file_period(self, full_csv):
        orb, ecc, m1, m2 = totest.get_samples_from_file(
            orbital_scheme='period',
            read_samples_from_file=full_csv)
        assert len(orb) == 3
        assert len(ecc) == 3
        assert len(m1) == 3
        assert len(m2) == 3
        np.testing.assert_array_equal(orb, [1.0, 10.0, 100.0])
        np.testing.assert_array_equal(m1, [10.0, 20.0, 30.0])
        np.testing.assert_array_equal(ecc, [0.0, 0.1, 0.2])

    def test_get_samples_from_file_separation(self, full_csv):
        orb, ecc, m1, m2 = totest.get_samples_from_file(
            orbital_scheme='separation',
            read_samples_from_file=full_csv)
        np.testing.assert_array_equal(orb, [50.0, 100.0, 200.0])

    def test_get_samples_from_file_missing_columns(self, minimal_csv):
        """File has no recognized columns — triggers random generation."""
        orb, ecc, m1, m2 = totest.get_samples_from_file(
            orbital_scheme='period',
            read_samples_from_file=minimal_csv,
            RNG=np.random.default_rng(seed=42))
        assert len(orb) == 2
        assert len(ecc) == 2
        assert len(m1) == 2
        assert len(m2) == 2

    def test_get_samples_from_file_missing_columns_separation(self, minimal_csv):
        """Separation scheme with no recognized columns."""
        orb, ecc, m1, m2 = totest.get_samples_from_file(
            orbital_scheme='separation',
            read_samples_from_file=minimal_csv,
            RNG=np.random.default_rng(seed=42))
        assert len(orb) == 2

    def test_get_samples_from_file_with_number(self, full_csv):
        """Request more binaries than in file — triggers expansion."""
        orb, ecc, m1, m2 = totest.get_samples_from_file(
            orbital_scheme='period',
            read_samples_from_file=full_csv,
            number_of_binaries=5)
        assert len(orb) == 5
        assert len(ecc) == 5
        assert len(m1) == 5
        assert len(m2) == 5

    def test_get_samples_from_file_with_index(self, full_csv):
        """Request subset with index offset."""
        orb, ecc, m1, m2 = totest.get_samples_from_file(
            orbital_scheme='period',
            read_samples_from_file=full_csv,
            number_of_binaries=2,
            index=1)
        assert len(orb) == 2
        assert orb[0] == 10.0  # second row from original

    # --- get_kick_samples_from_file ---

    def test_get_kick_samples_from_file_missing_kwarg(self):
        with raises(KeyError, match="no 'read_samples_from_file'"):
            totest.get_kick_samples_from_file()

    def test_get_kick_samples_from_file_not_found(self):
        with raises(FileNotFoundError, match="not found"):
            totest.get_kick_samples_from_file(
                read_samples_from_file='nonexistent.csv')

    def test_get_kick_samples_from_file_full(self, kick_csv):
        k1, k2 = totest.get_kick_samples_from_file(
            read_samples_from_file=kick_csv)
        assert k1.shape == (2, 4)
        assert k2.shape == (2, 4)
        assert k1[0, 0] == 100.0  # s1 velocity row 0
        assert k2[1, 0] == 150.0  # s2 velocity row 1

    def test_get_kick_samples_from_file_no_columns(self, no_kick_csv):
        """No kick columns — all set to None arrays."""
        k1, k2 = totest.get_kick_samples_from_file(
            read_samples_from_file=no_kick_csv)
        assert k1.shape == (2, 4)
        assert k2.shape == (2, 4)
        # All values should be None
        assert all(v is None for v in k1.flatten())
        assert all(v is None for v in k2.flatten())

    def test_get_kick_samples_from_file_with_number(self, kick_csv):
        """Request more binaries than in file — triggers expansion."""
        k1, k2 = totest.get_kick_samples_from_file(
            read_samples_from_file=kick_csv,
            number_of_binaries=5)
        assert k1.shape == (5, 4)
        assert k2.shape == (5, 4)

    def test_get_kick_samples_from_file_with_index(self, kick_csv):
        """Request subset with index offset."""
        k1, k2 = totest.get_kick_samples_from_file(
            read_samples_from_file=kick_csv,
            number_of_binaries=1,
            index=1)
        assert k1.shape == (1, 4)
        assert k1[0, 0] == 200.0  # second row velocity
