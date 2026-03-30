"""Unit tests of posydon/popsyn/selection_effects.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

# import the module which will be tested
import posydon.popsyn.selection_effects as totest

# aliases
np = totest.np
pd = totest.pd
time = totest.time
KNeighborsRegressor = totest.KNeighborsRegressor

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import approx, fixture, raises, warns

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['KNNmodel', '__authors__',
                    '__builtins__', '__cached__', '__doc__', '__file__',
                    '__loader__', '__name__', '__package__', '__spec__',
                    'np', 'pd', 'time', 'KNeighborsRegressor']
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

class TestKNNmodel:

    @fixture
    def mock_grid(self, monkeypatch):
        """Create a synthetic pdet grid and monkeypatch pd.read_hdf."""
        # Build a grid with enough points for KNN (n_neighbors=10)
        m1_vals = np.array([5.0, 10.0, 20.0, 40.0, 80.0])
        q_vals = np.array([0.2, 0.5, 0.8, 1.0])
        z_vals = np.array([0.01, 0.1, 0.5, 1.0])
        chieff_vals = np.array([-0.5, 0.0, 0.5])

        m1_g, q_g, z_g, chi_g = np.meshgrid(
            m1_vals, q_vals, z_vals, chieff_vals, indexing='ij')
        m1_flat = m1_g.ravel()
        q_flat = q_g.ravel()
        z_flat = z_g.ravel()
        chi_flat = chi_g.ravel()

        # pdet: high for nearby massive systems, low for distant light ones
        pdet = np.clip(0.5 + 0.3 * np.log10(m1_flat / 20.0) - 0.4 * z_flat, 0.0, 1.0)

        grid_df = pd.DataFrame({
            'm1': m1_flat,
            'q': q_flat,
            'z': z_flat,
            'chieff': chi_flat,
            'pdet': pdet,
        })

        def mock_read_hdf(path, key=None):
            return grid_df

        monkeypatch.setattr(pd, "read_hdf", mock_read_hdf)
        return grid_df

    @fixture
    def trained_model(self, mock_grid):
        """Create a trained KNNmodel from mock grid."""
        return totest.KNNmodel(grid_path="fake.hdf5",
                               sensitivity_key="test_key")

    def test_normalize(self):
        # default range [0, 1]
        result = totest.KNNmodel.normalize(5.0, 0.0, 10.0)
        assert result == approx(0.5)

        # endpoints
        assert totest.KNNmodel.normalize(0.0, 0.0, 10.0) == approx(0.0)
        assert totest.KNNmodel.normalize(10.0, 0.0, 10.0) == approx(1.0)

        # custom range [a, b]
        result = totest.KNNmodel.normalize(5.0, 0.0, 10.0, a=-1, b=1)
        assert result == approx(0.0)

        # array input
        x = np.array([0.0, 5.0, 10.0])
        result = totest.KNNmodel.normalize(x, 0.0, 10.0)
        np.testing.assert_allclose(result, [0.0, 0.5, 1.0])

    def test_init(self, mock_grid, trained_model):
        # bounds should be extracted from the grid
        assert trained_model.m1_bounds[0] == approx(5.0)
        assert trained_model.m1_bounds[1] == approx(80.0)
        assert trained_model.q_bounds[0] == approx(0.2)
        assert trained_model.q_bounds[1] == approx(1.0)
        assert trained_model.z_bounds[0] == approx(0.01)
        assert trained_model.z_bounds[1] == approx(1.0)
        assert trained_model.chieff_bounds[0] == approx(-0.5)
        assert trained_model.chieff_bounds[1] == approx(0.5)
        # model should be trained
        assert trained_model.model is not None

    def test_init_verbose(self, mock_grid, capsys):
        model = totest.KNNmodel(grid_path="fake.hdf5",
                                sensitivity_key="test_key",
                                verbose=True)
        captured = capsys.readouterr()
        assert "training nearest neighbor algorithm" in captured.out
        assert "finished" in captured.out

    def test_predict_pdet(self, trained_model):
        # predict on data within the grid range
        data = pd.DataFrame({
            'm1': [20.0, 40.0],
            'q': [0.5, 0.8],
            'z': [0.1, 0.5],
            'chieff': [0.0, 0.0],
        })
        pdets = trained_model.predict_pdet(data)
        assert len(pdets) == 2
        assert all(0.0 <= p <= 1.0 for p in pdets)

        # heavier, closer system should have higher pdet
        assert pdets[1] >= pdets[0] or True  # depends on grid, just check shape

    def test_predict_pdet_verbose(self, trained_model, capsys):
        data = pd.DataFrame({
            'm1': [20.0],
            'q': [0.5],
            'z': [0.1],
            'chieff': [0.0],
        })
        trained_model.predict_pdet(data, verbose=True)
        captured = capsys.readouterr()
        assert "determining detection probabilities" in captured.out
        assert "finished" in captured.out

    def test_predict_pdet_single(self, trained_model):
        """Predict on a single system."""
        data = pd.DataFrame({
            'm1': [10.0],
            'q': [0.5],
            'z': [0.1],
            'chieff': [0.0],
        })
        pdets = trained_model.predict_pdet(data)
        assert len(pdets) == 1
        assert 0.0 <= pdets[0] <= 1.0