"""Unit tests of posydon/popsyn/synthetic_population.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

# import the module which will be tested
import posydon.popsyn.synthetic_population as totest
# aliases
np = totest.np
pd = totest.pd

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns, approx
from inspect import isroutine, isclass
import warnings
warnings.simplefilter("always")
import os

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['parameter_array', 'DFInterface','History','Oneline',
                    'Population','PopulationIO','PopulationRunner',
                    'Rates','TransientPopulation',
                    '__authors__','__builtins__', '__cached__', '__doc__', 
                    '__file__','__loader__', '__name__', '__package__', '__spec__', 
                    'np', 'pd', 'tqdm', 'os', 'plt', 
                    'Zsun', 'binarypop_kwargs_from_ini',
                    'initial_total_underlying_mass','plot_pop',
                    'convert_metallicity_to_string','Pwarn','cosmology','const',
                    'get_shell_comoving_volume', 'get_comoving_distance_from_redshift',
                    'get_cosmic_time_from_redshift', 'redshift_from_cosmic_time_interpolator',
                    'DEFAULT_SFH_MODEL', 'get_redshift_bin_edges',
                    'get_redshift_bin_centers', 'SFR_per_met_at_z',
                    'BinaryPopulation', 'HISTORY_MIN_ITEMSIZE','ONELINE_MIN_ITEMSIZE'
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

class TestPopulationRunner:
        
    def test_init(self):
        # missing argument
        with raises(TypeError,match="missing 1 required positional argument: 'path_to_ini'"):
            totest.PopulationRunner()
        # bad input
        with raises(ValueError, match="You did not provide a valid path_to_ini!"):
            totest.PopulationRunner("test")

    def test_evolve(self,monkeypatch):
        # mock dependencies
        class DummyPop:
            def __init__(self, **kwargs):
                self.kwargs = kwargs
                self.comm = None
                self.metallicity = kwargs["metallicity"]
            def evolve(self):
                self.evolved = True
            def combine_saved_files(self, *args):
                self.combined = True
        def dummy_kwargs(path):
            return {
                "metallicity": 0.1,
                "temp_directory": "tmp_dir",
                "verbose": False}
        def dummy_kwargs_list(path):
            return {
                "metallicity": [0.1,1.],
                "temp_directory": "tmp_dir",
                "verbose": False}
        def dummy_merge(pop):
            pop.merged = True
        # Mock out functions
        monkeypatch.setattr(totest, "binarypop_kwargs_from_ini", dummy_kwargs)
        monkeypatch.setattr(totest, "BinaryPopulation", DummyPop)
        monkeypatch.setattr(totest, "convert_metallicity_to_string", lambda x: "0.1")
        # overwrite=False, directory doesn't exist
        monkeypatch.setattr(os.path, "exists", lambda path: False)
        run = totest.PopulationRunner("dummy.ini")
        run.merge_parallel_runs = dummy_merge
        run.evolve()
        assert run.binary_populations[0].evolved is True
        assert run.binary_populations[0].merged is True
        # overwrite=False, directory exists
        monkeypatch.setattr(os.path, "exists", lambda path: True)
        monkeypatch.setattr(totest, "binarypop_kwargs_from_ini", dummy_kwargs_list)
        run = totest.PopulationRunner("dummy.ini", verbose=True)
        with raises(FileExistsError, match="tmp_dir"):
            run.evolve(overwrite=False)
        # overwrite=True, directory exists
        removed = {}
        monkeypatch.setattr(os, "removedirs", lambda path: removed.setdefault("called", path))
        run.merge_parallel_runs = dummy_merge
        run.evolve(overwrite=True)
        assert removed["called"] == "0.1_Zsun_tmp_dir"
        assert run.binary_populations[0].evolved is True
        assert run.binary_populations[0].merged is True

    def test_merge_parallel_runs(self, tmp_path, monkeypatch, capsys):
        class DummyPop:
            def __init__(self, metallicity, temp_directory,**kwargs):
                self.metallicity = metallicity
                self.kwargs = {"temp_directory": temp_directory}
                self.combine_args = None
                self.combined = False

            def combine_saved_files(self, out_path, files):
                self.combine_args = (out_path, files)
                self.combined = True
        def dummy_kwargs(path):
            return {
                "metallicity": 0.1,
                "temp_directory": "tmp_dir",
                "verbose": False}

        monkeypatch.setattr(totest, "binarypop_kwargs_from_ini", dummy_kwargs)
        monkeypatch.setattr(totest, "BinaryPopulation", DummyPop)
        monkeypatch.setattr(totest, "convert_metallicity_to_string", 
                            lambda x: str(os.path.join(tmp_path, "0.1")))
        # 1) File exists case: should raise FileExistsError
        pop = DummyPop(metallicity=0.1, temp_directory=str(tmp_path))
        output_file = os.path.join(tmp_path,"0.1_Zsun_population.h5")
        with open(output_file, "w") as f:
            f.write("test")
        run = totest.PopulationRunner("dummy.ini")
        run.verbose = False
        with raises(FileExistsError, match="Files were not merged"):
            run.merge_parallel_runs(pop)

        # 2) Normal merge case
        temp_dir = tmp_path
        file1 = os.path.join(temp_dir,"file1.tmp")
        file2 = os.path.join(temp_dir,"file2.tmp")
        with open(file1, "w") as f:
            f.write("test")
        with open(file2, "w") as f:
            f.write("test")
        pop = DummyPop(metallicity=0.1, temp_directory=str(temp_dir))
        run = totest.PopulationRunner("dummy.ini")
        run.verbose = True
        monkeypatch.setattr(totest, "convert_metallicity_to_string", lambda x: "0.1")
        run.merge_parallel_runs(pop)
        assert pop.combined is True
        out_path, files = pop.combine_args
        assert out_path == "0.1_Zsun_population.h5"
        # Filter out output file if somehow included (defensive)
        filtered_files = [f for f in files if os.path.basename(f) != out_path]
        assert set(os.path.basename(f) for f in filtered_files) == {"file1.tmp", "file2.tmp"}
        captured = capsys.readouterr()
        assert "Merging 3 files..." in captured.out
        assert "Files merged!" in captured.out
        assert f"Removing files in {temp_dir}" in captured.out
        # Since files still exist, directory not removed; remove files manually
        for f in [file1, file2]:
            os.remove(f)
        output_file = os.path.join(temp_dir, "0.1_Zsun_population.h5")
        if os.path.exists(output_file):
            os.remove(output_file)
        assert len(os.listdir(temp_dir)) == 0
        run.verbose = False
        run.merge_parallel_runs(pop)
        assert not os.path.exists(temp_dir)
        monkeypatch.undo()
             
class TestDFInterface:
    
    def test_head_tail_select(self, tmp_path):
        # Setup test HDF5 file
        data = pd.DataFrame({
            "index": np.repeat(np.arange(5), 2),
            "time": np.random.rand(10),
            "value": np.random.rand(10)
        })
        hdf_path = os.path.join(tmp_path,"test.h5")
        data.to_hdf(hdf_path, key="history", format="table", index=False)

        dfi = totest.DFInterface()
        dfi.filename = str(hdf_path)
        dfi.chunksize = 3

        head = dfi.head("history", n=3)
        tail = dfi.tail("history", n=2)
        subset = dfi.select("history", columns=["time"])

        assert len(head) == 3
        assert len(tail) == 2
        assert "time" in subset.columns
        assert subset.shape[1] == 1

    def test_repr_methods(self, tmp_path):
        df = pd.DataFrame({"index": range(10), "x": np.random.rand(10)})
        path = os.path.join(tmp_path, "test_repr.h5")
        df.to_hdf(path, key="history", format="table", index=False)

        dfi = totest.DFInterface()
        dfi.filename = str(path)

        s = dfi.get_repr("history")
        html = dfi.get_html_repr("history")

        assert isinstance(s, str)
        assert "x" in s
        assert isinstance(html, str)
        assert "<table" in html
        
class TestHistory:
    
    def test_init(self, tmp_path):
        with raises(FileNotFoundError, match="does not exist!"):
            totest.History("nonexistent_file.h5")
            
        df = pd.DataFrame({
            "binary_index": np.repeat(np.arange(3), 2),
            "a": np.random.rand(6),
        })
        file_path = os.path.join(tmp_path, "test_history.h5")
        df.set_index("binary_index", inplace=True)
        df.to_hdf(file_path, key="history", format="table")

        hist = totest.History(str(file_path), verbose=True, chunksize=2)

        assert hist.filename == str(file_path)
        assert hist.lengths is not None
        assert hist.columns == ["a"]
        assert hist.number_of_systems == 3
        assert isinstance(hist.indices, np.ndarray)

    def test_getitem_and_len(self, tmp_path):
        df = pd.DataFrame({
            "binary_index": np.repeat(np.arange(3), 2),
            "a": np.random.rand(6),
        })
        file_path = os.path.join(tmp_path, "test_getitem.h5")
        df.set_index("binary_index", inplace=True)
        df.to_hdf(file_path, key="history", format="table")
        
        lengths_df = pd.DataFrame({'lengths': [2, 2, 2]}, index=[0, 1, 2])
        lengths_df.to_hdf(file_path, key="history_lengths")

        hist = totest.History(str(file_path))
        
        # adding history_lengths
        assert hist.lengths.equals(lengths_df)

        # __len__
        assert len(hist) == 6
        
        # __getitem__ with "none" slice indices
        out_startnone = hist[:3]
        out_stopnone = hist[2:]
        assert not out_startnone.empty
        assert all(i in out_startnone.index for i in range(3))
        assert not out_stopnone.empty
        assert all(i>=2 for i in out_stopnone.index)

        # __getitem__ with int
        out = hist[0]
        assert isinstance(out, pd.DataFrame)

        # __getitem__ with list of int
        out = hist[[0, 1]]
        assert isinstance(out, pd.DataFrame)
        out_none = hist[[]]
        assert isinstance(out_none,pd.DataFrame)
        assert out_none.empty

        # __getitem__ with numpy array of int
        out = hist[np.array([0, 1])]
        assert isinstance(out, pd.DataFrame)
        out_none = hist[np.array([], dtype=int)]
        assert isinstance(out_none,pd.DataFrame)
        assert out_none.empty

        # __getitem__ with bool array
        full_data = pd.read_hdf(file_path, key="history")
        mask = full_data["a"] > -1
        empty_mask = np.array([],dtype=bool)
        out = hist[mask.to_numpy()] 
        out_none = hist[empty_mask]
        assert not out.empty
        assert (out["a"] > -1).all()
        assert isinstance(out_none,pd.DataFrame)
        assert out_none.empty

        # __getitem__ with str column
        out = hist["a"]
        assert "a" in out.columns

        # __getitem__ with list of str columns
        out = hist[["a"]]
        assert "a" in out.columns

        # Invalid column
        with raises(ValueError, match="is not a valid column name"):
            hist["bad_column"]

        # Invalid list of column names
        with raises(ValueError, match="Not all columns in"):
            hist[["a", "bad"]]

        # Invalid type
        with raises(ValueError, match="Invalid key type"):
            hist[None]

    def test_slice(self, tmp_path):
        df = pd.DataFrame({
            "index": np.repeat(np.arange(5), 2),
            "val": np.random.rand(10)
        })
        path = os.path.join(tmp_path, "test_slice.h5")
        df.to_hdf(path, key="history", format="table", index=False)
        hist = totest.History(str(path))
        sliced = hist[1:3]
        assert isinstance(sliced, pd.DataFrame)

    def test_head_tail_repr(self, tmp_path):
        df = pd.DataFrame({
            "index": np.repeat(np.arange(5), 2),
            "val": np.random.rand(10)
        })
        path = os.path.join(tmp_path, "test_repr2.h5")
        df.to_hdf(path, key="history", format="table", index=False)
        hist = totest.History(str(path))

        head = hist.head(n=3)
        tail = hist.tail(n=2)
        rep = repr(hist)
        html = hist._repr_html_()

        assert isinstance(head, pd.DataFrame)
        assert isinstance(tail, pd.DataFrame)
        assert isinstance(rep, str)
        assert isinstance(html, str)
        assert "val" in rep
        assert "<table" in html
        
class TestOneline:
    
    @fixture
    def fix(self):
#         return 
        pass

        pass
    def test_head(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_tail(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_select(self):
        # missing argument
        # bad input
        # examples
        pass
        
class TestPopulation:
    
    @fixture
    def fix(self):
#         return 

        pass
    def test_calculate_underlying_mass(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_export_selection(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_formation_channels(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_calculate_formation_channels(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_columns(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_create_transient_population(self):
        # missing argument
        # bad input
        # examples
        pass
        
class TestTransientPopulation:
    
    @fixture
    def fix(self):
#         return 
        pass

    def test_population(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_columns(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_select(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_get_efficiency_over_metallicity(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_calculate_cosmic_weights(self):
        # missing argument
        # bad input
        # examples
        pass
        
class TestRates:
    
    @fixture
    def fix(self):
        # return 
        pass

    def test_weights(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_z_birth(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_z_events(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_select_rate_slice(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_calculate_intrinsic_rate_density(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_calculate_observable_population(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_observable_population(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_observable_population_names(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_intrinsic_rate_density(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_edges_metallicity_bins(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_centers_metallicity_bins(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_edges_redshift_bins(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_centers_redshift_bins(self):
        # missing argument
        # bad input
        # examples    
        pass
        